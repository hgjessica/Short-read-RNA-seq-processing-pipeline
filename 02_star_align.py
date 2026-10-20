#!/usr/bin/env python3
"""Step 2: STAR alignment + BAM index.

    python 02_star_align.py --config config.json align --sample sample_1
    python 02_star_align.py --config config.json align --sample sample_1 --stage trimmed
    python 02_star_align.py --config config.json submit align [--stage trimmed]

Works for paired-end (fastq_1 + fastq_2) and single-end (fastq_1 only) samples. The STAR genome
index must already exist (config "star_index"; built with the same STAR version).
Output: star/<sample>/<sample>.Aligned.sortedByCoord.out.bam (+ .bai, Log.final.out)
        (star_trimmed/ for --stage trimmed)
"""
import sys
from pathlib import Path

from common import *  # noqa: F401,F403  (helpers shared by all steps)
from common import submit_array


def step_align(ctx, args):
    s = pick_sample(ctx.cfg, args.sample)
    name = s["name"]
    files = reads(ctx, s, args.stage)
    index = ctx.cfg["star_index"]
    need(ctx, *files, Path(index) / "SA", Path(index) / "genomeParameters.txt")
    out = stage_dir(ctx, "star", args.stage) / name
    bam = out / f"{name}.Aligned.sortedByCoord.out.bam"
    if bam.exists() and not args.force:
        print(f"[SKIP] {name}: {bam} exists (use --force to redo)")
        return
    raise_open_file_limit()
    task_tmp(ctx, "star")
    if not ctx.dry:
        out.mkdir(parents=True, exist_ok=True)
    st = ctx.cfg.get("star", {})
    cmd = ["STAR", "--genomeDir", index,
           "--runThreadN", ctx.threads,
           "--runMode", "alignReads",
           "--readFilesCommand", read_files_command(files),
           "--readFilesIn", *files,
           "--outFileNamePrefix", f"{out}/{name}.",
           "--outSAMtype", "BAM", "SortedByCoordinate",
           "--outSAMunmapped", "Within",
           "--outSAMattributes", *st.get("sam_attributes", ["Standard"])]
    if st.get("limit_bam_sort_ram"):
        cmd += ["--limitBAMsortRAM", st["limit_bam_sort_ram"]]
    cmd += st.get("extra_args", [])
    if run(ctx, cmd):
        sys.exit(f"[ERROR] STAR failed for {name}")
    need(ctx, bam)
    if run(ctx, ["samtools", "index", "-@", ctx.threads, bam]):
        sys.exit(f"[ERROR] samtools index failed for {name}")
    print(f"[OK] {bam}")


def main():
    ap = base_parser(__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("align", help="STAR alignment + samtools index for one sample")
    p.add_argument("--sample", help="sample name (default: SLURM array task id)")
    add_stage_arg(p)
    p.add_argument("--force", action="store_true", help="redo even if the BAM exists")
    p.set_defaults(fn=step_align)

    add_submit(sub, ["align"]).set_defaults(fn=lambda c, a: submit_array(c, a, __file__))

    args = ap.parse_args()
    args.fn(Ctx(load_config(args.config), args.dry_run), args)


if __name__ == "__main__":
    main()
