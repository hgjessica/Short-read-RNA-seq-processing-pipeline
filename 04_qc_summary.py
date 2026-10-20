#!/usr/bin/env python3
"""Step 4: collect all FastQC reports into one CSV (one row per FASTQ file).

    python 04_qc_summary.py --config config.json summary [--stage raw|trimmed] [--basic-only]

Reads every *_fastqc.zip in fastqc/ (or fastqc_trimmed/) and writes a table with
Total Sequences, Sequence Length and the PASS/WARN/FAIL status of each FastQC module.
--basic-only writes just the first two columns (read count and length).
Output: fastqc/fastqc_summary.csv (fastqc_trimmed/fastqc_summary.csv for --stage trimmed)
"""
import csv
import sys
import zipfile

from common import *  # noqa: F401,F403  (helpers shared by all steps)


def parse_zip(zpath):
    """Return (total_sequences, sequence_length, {module: status}) from one FastQC zip."""
    base = zpath.name[:-len(".zip")]
    total, length, status = "N/A", "N/A", {}
    with zipfile.ZipFile(zpath) as zf:
        for line in zf.read(f"{base}/fastqc_data.txt").decode().splitlines():
            if line.startswith("Total Sequences"):
                total = line.split("\t")[1].strip()
            elif line.startswith("Sequence length"):
                length = line.split("\t")[1].strip()
        for line in zf.read(f"{base}/summary.txt").decode().strip().splitlines():
            parts = line.split("\t")
            if len(parts) >= 2:
                status[parts[1]] = parts[0]
    return total, length, status


def step_summary(ctx, args):
    d = stage_dir(ctx, "fastqc", args.stage)
    zips = sorted(d.glob("*_fastqc.zip")) if d.is_dir() else []
    if not zips:
        sys.exit(f"[ERROR] no *_fastqc.zip files in {d}")
    rows, modules = {}, []
    for z in zips:
        try:
            total, length, status = parse_zip(z)
        except (zipfile.BadZipFile, KeyError) as e:
            print(f"[WARN] skipping {z.name}: {e}")
            continue
        for m in status:
            if m not in modules:
                modules.append(m)
        rows[z.name[:-len("_fastqc.zip")]] = (total, length, status)
    out = d / "fastqc_summary.csv"
    if ctx.dry:
        print(f"[DRY] would write {len(rows)} rows to {out}")
        return
    with open(out, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Sample", "Total Sequences", "Sequence Length"] + ([] if args.basic_only else modules))
        for name in sorted(rows):
            total, length, status = rows[name]
            extra = [] if args.basic_only else [status.get(m, "N/A") for m in modules]
            w.writerow([name, total, length] + extra)
    print(f"[OK] CSV written to {out}")


def main():
    ap = base_parser(__doc__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("summary", help="merge FastQC zips into one CSV")
    add_stage_arg(p)
    p.add_argument("--basic-only", action="store_true", help="only read count and length")
    p.set_defaults(fn=step_summary)
    args = ap.parse_args()
    args.fn(Ctx(load_config(args.config), args.dry_run), args)


if __name__ == "__main__":
    main()
