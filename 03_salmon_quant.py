#!/usr/bin/env python3
"""Step 3: Salmon quantification (mapping-based mode, selective alignment).

    python 03_salmon_quant.py --config config.json quant --sample sample_1
    python 03_salmon_quant.py --config config.json quant --sample sample_1 --stage trimmed
    python 03_salmon_quant.py --config config.json submit quant [--stage trimmed]

Reads the FASTQ files directly (not the STAR BAM). Library type is auto-detected (-l A), which
also reports strandedness in lib_format_counts.json / the log. Paired-end uses -1/-2, single-end -r.
Output: salmon/<sample>/quant.sf (salmon_trimmed/ for --stage trimmed)
"""
import sys

from common import *  # noqa: F401,F403  (helpers shared by all steps)
from common import submit_array


def step_quant(ctx, args):
    s = pick_sample(ctx.cfg, args.sample)
    name = s["name"]
    files = reads(ctx, s, args.stage)
    index = ctx.cfg["salmon_index"]
    need(ctx, *files, index)
    out = stage_dir(ctx, "salmon", args.stage) / name
    if (out / "quant.sf").exists() and not args.force:
        print(f"[SKIP] {name}: {out / 'quant.sf'} exists (use --force to redo)")
        return
    raise_open_file_limit()
    task_tmp(ctx, "salmon")
    if not ctx.dry:
        out.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["salmon", "quant", "-i", index, "-l", ctx.cfg.get("salmon", {}).get("libtype", "A")]
    cmd += ["-1", files[0], "-2", files[1]] if len(files) == 2 else ["-r", files[0]]
    cmd += ["--validateMappings", "--threads", ctx.threads, "-o", out]
    cmd += ctx.cfg.get("salmon", {}).get("extra_args", [])
    if run(ctx, cmd):
        sys.exit(f"[ERROR] salmon failed for {name}")
    need(ctx, out / "quant.sf")
    print(f"[OK] {out / 'quant.sf'}")


def main():
    ap = base_parser(__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)

    p = sub.add_parser("quant", help="salmon quant for one sample")
    p.add_argument("--sample", help="sample name (default: SLURM array task id)")
    add_stage_arg(p)
    p.add_argument("--force", action="store_true", help="redo even if quant.sf exists")
    p.set_defaults(fn=step_quant)

    add_submit(sub, ["quant"]).set_defaults(fn=lambda c, a: submit_array(c, a, __file__))

    args = ap.parse_args()
    args.fn(Ctx(load_config(args.config), args.dry_run), args)


if __name__ == "__main__":
    main()
