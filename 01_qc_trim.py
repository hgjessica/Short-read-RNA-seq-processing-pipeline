#!/usr/bin/env python3
"""Step 1: FastQC and adapter trimming (fastp).

    python 01_qc_trim.py --config config.json fastqc [--stage raw|trimmed] --sample sample_1
    python 01_qc_trim.py --config config.json trim --sample sample_1
    python 01_qc_trim.py --config config.json submit trim
    python 01_qc_trim.py --config config.json submit fastqc --stage trimmed

fastqc runs on the raw reads (default) or on the fastp output (--stage trimmed).
trim runs fastp with default settings (automatic adapter detection for paired-end data).
Output: fastqc/<file>_fastqc.zip|html (fastqc_trimmed/ for --stage trimmed),
        trimmed/<sample>_1.fq.gz [_2.fq.gz], trimmed/<sample>_fastp.json|html
"""
import sys

from common import *  # noqa: F401,F403  (helpers shared by all steps)
from common import submit_array


def step_fastqc(ctx, args):
    s = pick_sample(ctx.cfg, args.sample)
    files = reads(ctx, s, args.stage)
    need(ctx, *files)
    out = stage_dir(ctx, "fastqc", args.stage)
    tmp = task_tmp(ctx, "fastqc")
    if not ctx.dry:
        out.mkdir(parents=True, exist_ok=True)
    # FastQC uses one thread per input file, so more threads than files are wasted
    rc = run(ctx, ["fastqc", "-t", min(len(files), int(ctx.threads)),
                   "-d", tmp, "-o", out, "-f", "fastq", *files])
    if rc:
        sys.exit(f"[ERROR] fastqc failed for {s['name']}")
    print(f"[OK] {out}")


def step_trim(ctx, args):
    s = pick_sample(ctx.cfg, args.sample)
    files = reads(ctx, s, "raw")
    outs = reads(ctx, s, "trimmed")
    json_f = ctx.path("trimmed", f"{s['name']}_fastp.json")
    html_f = ctx.path("trimmed", f"{s['name']}_fastp.html")
    need(ctx, *files)
    if outs[0].exists() and not args.force:
        print(f"[SKIP] {s['name']}: {outs[0]} exists (use --force to redo)")
        return
    if not ctx.dry:
        ctx.path("trimmed").mkdir(parents=True, exist_ok=True)
    cmd = ["fastp", "-i", files[0], "-o", outs[0]]
    if len(files) == 2:
        cmd += ["-I", files[1], "-O", outs[1]]
    # fastp uses at most 16 threads
    cmd += ["--thread", min(int(ctx.threads), 16), "--json", json_f, "--html", html_f]
    cmd += ctx.cfg.get("fastp", {}).get("extra_args", [])
    if run(ctx, cmd):
        sys.exit(f"[ERROR] fastp failed for {s['name']}")
    print(f"[OK] {outs[0]}")


def main():
    ap = base_parser(__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("fastqc", help="FastQC on one sample's reads")
    p.add_argument("--sample", help="sample name (default: SLURM array task id)")
    add_stage_arg(p)
    p.set_defaults(fn=step_fastqc)

    p = sub.add_parser("trim", help="fastp adapter trimming for one sample")
    p.add_argument("--sample", help="sample name (default: SLURM array task id)")
    p.add_argument("--force", action="store_true", help="redo even if trimmed output exists")
    p.set_defaults(fn=step_trim)

    add_submit(sub, ["fastqc", "trim"]).set_defaults(fn=lambda c, a: submit_array(c, a, __file__))

    args = ap.parse_args()
    args.fn(Ctx(load_config(args.config), args.dry_run), args)


if __name__ == "__main__":
    main()
