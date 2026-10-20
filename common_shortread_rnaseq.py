"""Shared helpers for the short-read RNA-seq step scripts (config, running commands, sample selection)."""
import argparse
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path

# Per-stage output folder names. "raw" = reads exactly as delivered, "trimmed" = fastp output.
STAGES = ("raw", "trimmed")


class Ctx:
    def __init__(self, cfg, dry_run):
        self.cfg = cfg
        self.dry = dry_run
        self.work = Path(cfg["work_dir"])
        self.threads = str(cfg.get("threads", 8))

    def path(self, *parts):
        return self.work.joinpath(*parts)


def load_config(path):
    with open(path) as fh:
        cfg = json.load(fh)
    for key in ("work_dir", "samples"):
        if key not in cfg:
            sys.exit(f"config is missing required key: {key}")
    names = [s.get("name") for s in cfg["samples"]]
    if len(set(names)) != len(names) or None in names:
        sys.exit("config: every sample needs a unique 'name'")
    return cfg


def run(ctx, cmd, **kw):
    """Run a command (list form, no shell). Returns the exit code."""
    print("$ " + " ".join(shlex.quote(str(c)) for c in cmd), flush=True)
    if ctx.dry:
        return 0
    try:
        return subprocess.run([str(c) for c in cmd], **kw).returncode
    except FileNotFoundError:
        sys.exit(f"[ERROR] tool not found on PATH: {cmd[0]} (activate the conda environment?)")


def need(ctx, *files):
    """Abort with a clear message if an input is missing (skipped in dry run)."""
    for f in files:
        if not ctx.dry and not Path(f).exists():
            sys.exit(f"[ERROR] missing: {f}")


def raise_open_file_limit(n=10000):
    """Same as `ulimit -n 10000` in the original job scripts (STAR sorts via many temp files)."""
    try:
        import resource
        soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
        target = n if hard == resource.RLIM_INFINITY else min(n, hard)
        if soft == resource.RLIM_INFINITY or soft < target:
            resource.setrlimit(resource.RLIMIT_NOFILE, (target, hard))
    except (ImportError, ValueError, OSError):
        pass


def kept_samples(cfg):
    return [s for s in cfg["samples"] if s.get("keep", True)]


def pick_sample(cfg, name):
    """--sample NAME, or the SLURM array task id as an index into kept samples."""
    kept = kept_samples(cfg)
    if name:
        for s in kept:
            if s["name"] == name:
                return s
        sys.exit(f"unknown (or keep=false) sample: {name}")
    task = os.environ.get("SLURM_ARRAY_TASK_ID")
    if task is None:
        sys.exit("give --sample NAME (or run as a SLURM array task)")
    if int(task) >= len(kept):
        sys.exit(f"array task {task} is out of range ({len(kept)} kept samples)")
    return kept[int(task)]


def stage_dir(ctx, base, stage):
    """'star' for raw reads, 'star_trimmed' for trimmed reads (same for fastqc, salmon)."""
    return ctx.path(base if stage == "raw" else f"{base}_{stage}")


def is_paired(sample):
    return bool(sample.get("fastq_2"))


def reads(ctx, sample, stage="raw"):
    """Input FASTQ list for a sample: [R1, R2] (paired) or [R1] (single-end)."""
    if stage == "trimmed":
        n = sample["name"]
        files = [ctx.path("trimmed", f"{n}_1.fq.gz")]
        if is_paired(sample):
            files.append(ctx.path("trimmed", f"{n}_2.fq.gz"))
        return files
    base = Path(ctx.cfg.get("fastq_dir", "."))
    files = [base / sample["fastq_1"]]
    if is_paired(sample):
        files.append(base / sample["fastq_2"])
    if len(files) == 2 and files[0] == files[1]:
        sys.exit(f"[ERROR] {sample['name']}: fastq_1 and fastq_2 are the same file")
    return files


def read_files_command(files):
    """Decompressor for STAR --readFilesCommand (cat for plain text, as in the originals)."""
    return "zcat" if str(files[0]).endswith(".gz") else "cat"


def task_tmp(ctx, tag):
    """Per-task scratch directory so concurrent array tasks cannot collide."""
    base = Path(ctx.cfg.get("scratch_dir") or ctx.path(".tmp"))
    d = base / f"{tag}_{os.environ.get('SLURM_ARRAY_TASK_ID', 'local')}"
    if not ctx.dry:
        d.mkdir(parents=True, exist_ok=True)
    os.environ["TMPDIR"] = str(d)
    return d


def base_parser(doc):
    """Parser with the options every step script shares."""
    ap = argparse.ArgumentParser(description=doc.strip().split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog="\n".join(doc.strip().split("\n\n")[1:]))
    ap.add_argument("--config", required=True, help="JSON config (see config.example.json)")
    ap.add_argument("--dry-run", action="store_true", help="print commands, do not run")
    return ap


def add_stage_arg(p):
    p.add_argument("--stage", choices=STAGES, default="raw",
                   help="which reads to use: raw (as delivered) or trimmed (fastp output)")


def add_submit(sub, choices):
    """Add the `submit` subcommand (SLURM array, one task per kept sample)."""
    p = sub.add_parser("submit", help="submit a per-sample step as a SLURM array")
    p.add_argument("step_name", choices=sorted(choices))
    p.add_argument("step_args", nargs=argparse.REMAINDER,
                   help="extra args passed to the step, e.g. --stage trimmed")
    return p


def submit_array(ctx, args, script):
    """Write and submit a SLURM array script (one task per kept sample) for a step."""
    sl = ctx.cfg.get("slurm", {})
    n = len(kept_samples(ctx.cfg))
    if n == 0:
        sys.exit("no samples to submit")
    conc = sl.get("array_concurrency", 3)
    extra = " ".join(shlex.quote(a) for a in args.step_args)
    log_dir = ctx.path("slurm_logs")
    text = f"""#!/bin/bash
#SBATCH -p {sl.get('partition', 'PARTITION_NAME')}
#SBATCH -N 1
#SBATCH -n {ctx.threads}
#SBATCH --mem={sl.get('mem', '40G')}
#SBATCH -t {sl.get('time', '48:00:00')}
#SBATCH --export=ALL
#SBATCH -J {args.step_name}
#SBATCH --array=0-{n - 1}%{conc}
#SBATCH -o {log_dir}/{args.step_name}_%A_%a.out
#SBATCH -e {log_dir}/{args.step_name}_%A_%a.err

source {sl.get('conda_sh', 'Path_to_conda.sh')}
conda activate {sl.get('conda_env', 'rnaseq')}
python3 {Path(script).resolve()} --config {Path(args.config).resolve()} {args.step_name} {extra}
"""
    sbatch_file = ctx.path("slurm", f"{Path(script).stem}_{args.step_name}.sbatch")
    print(text)
    if ctx.dry:
        return
    log_dir.mkdir(parents=True, exist_ok=True)
    sbatch_file.parent.mkdir(parents=True, exist_ok=True)
    sbatch_file.write_text(text)
    subprocess.run(["sbatch", str(sbatch_file)])
