# Short-read-RNA-seq-processing-pipeline
Python scripts that run a standard bulk RNA-seq workflow one step at a time, locally or as SLURM job arrays: QC and trimming, STAR alignment, Salmon quantification and a FastQC summary table.


| Script | Covers |
|---|---|
| `01_qc_trim.py` | `fastqc` (raw or trimmed reads) and `trim` (fastp) |
| `02_star_align.py` | `align`: STAR to sorted BAM, then `samtools index` |
| `03_salmon_quant.py` | `quant`: Salmon `quant` straight from the FASTQ files |
| `04_qc_summary.py` | `summary`: all FastQC reports merged into one CSV |

`common.py` holds the shared helpers and must stay in the same folder. All paths, sample names and
tool settings live in one JSON config, so the code contains nothing specific to a project.

```
FASTQ ──fastqc──▶ QC reports ──summary──▶ fastqc_summary.csv
  │
  ├──fastp──▶ trimmed FASTQ ──fastqc──▶ (QC of trimmed reads)
  │
  ├── STAR ──▶ sorted + indexed BAM            (alignments, IGV, splicing)
  └── salmon ─▶ quant.sf                       (transcript TPM / counts)
```

Steps 2 and 3 each take `--stage raw` (default, the reads as delivered) or `--stage trimmed`
(fastp output). Both stages can coexist; trimmed results go to `*_trimmed/` folders.

## 1. Requirements

- Python >= 3.8 (standard library only)
- Tools on `PATH`: `fastqc`, `fastp`, `STAR`, `salmon`, `samtools`
- Optionally SLURM (`sbatch`) for the `submit` command

```bash
conda env create -f environment.yml
conda activate rnaseq
```

You also need two **pre-built indexes** for the same genome and annotation (the scripts do not build them):

- a STAR genome index directory (`star_index`), built with the same STAR version you run
  (`environment.yml` pins 2.7.10a, the version the original runs used)
- a Salmon transcriptome index directory (`salmon_index`)

Example commands for building them once (not run by this pipeline):

```bash
STAR --runMode genomeGenerate --genomeDir Path_to_STAR_index --genomeFastaFiles genome.fa \
     --sjdbGTFfile annotation.gtf --sjdbOverhang 99 --runThreadN 16
salmon index -t transcripts.fa -d decoys.txt -i Path_to_salmon_index   # genome-decoy-aware index recommended
```

## 2. Configure

Copy `config.example.json` to `config.json` and fill in every `Path_to_...` value.

| Key | Meaning |
|---|---|
| `work_dir` | All outputs are written here |
| `fastq_dir` | Folder with the input FASTQ files (sample file names are relative to it) |
| `star_index`, `salmon_index` | Pre-built indexes (see above) |
| `threads` | Threads per job (also the SLURM `-n`) |
| `scratch_dir` | Optional. Per-task temp dir (use a disk with plenty of free space) |
| `slurm` | Partition, memory, time limit, array concurrency, conda setup (used by `submit`) |
| `star` | `sam_attributes` (default `["Standard"]`), optional `limit_bam_sort_ram` (bytes), `extra_args` (list) |
| `salmon` | `libtype` (default `A`, auto-detect), `extra_args` |
| `fastp` | `extra_args` (default: none, i.e. fastp defaults) |
| `samples` | One entry per sample, see below |

Sample entries:

```json
{"name": "KO_mock_1", "fastq_1": "KO_mock_1_R1.fastq.gz", "fastq_2": "KO_mock_1_R2.fastq.gz", "keep": true}
```

- Paired-end: give `fastq_1` and `fastq_2`. Single-end: give only `fastq_1`; every step then switches
  to single-end mode automatically (fastp without `-I/-O`, STAR with one file, Salmon `-r`).
- Plain `.fastq`/`.fq` and gzipped files both work (STAR uses `cat` or `zcat` by file extension).
- `keep: false` skips a sample in every per-sample step and in `submit`.
- `name` is used for all output names; keep it free of spaces.

Keep `config.json` out of version control if it contains real names or paths.

## 3. Run

Add `--dry-run` to print the commands without running anything. Per-sample steps take
`--sample NAME`; under SLURM they use the array task id instead.

```bash
P1="python 01_qc_trim.py --config config.json"
P2="python 02_star_align.py --config config.json"
P3="python 03_salmon_quant.py --config config.json"
P4="python 04_qc_summary.py --config config.json"
```

| Step | Command | Underlying tool call | Output |
|---|---|---|---|
| FastQC (raw) | `$P1 submit fastqc` | `fastqc -t <n_files> -d <tmp> -o fastqc -f fastq R1 [R2]` | `fastqc/*_fastqc.zip, .html` |
| Trim | `$P1 submit trim` | `fastp -i R1 -I R2 -o .. -O .. --thread 16 --json .. --html ..` | `trimmed/<sample>_1.fq.gz`, `_2.fq.gz`, `_fastp.json/.html` |
| FastQC (trimmed) | `$P1 submit fastqc --stage trimmed` | same, on `trimmed/` | `fastqc_trimmed/` |
| Align | `$P2 submit align [--stage trimmed]` | `STAR --runMode alignReads --readFilesCommand cat/zcat --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard`, then `samtools index` | `star/<sample>/<sample>.Aligned.sortedByCoord.out.bam` + `.bai`, `Log.final.out`, `SJ.out.tab` |
| Quantify | `$P3 submit quant [--stage trimmed]` | `salmon quant -i <index> -l A -1 R1 -2 R2 --validateMappings --threads <n> -o ..` | `salmon/<sample>/quant.sf` |
| QC summary | `$P4 summary [--stage trimmed] [--basic-only]` | (reads the FastQC zips) | `fastqc/fastqc_summary.csv` |

Notes:
- `submit STEP [args]` writes `slurm/<script>_<STEP>.sbatch` (one array task per kept sample, at most
  `array_concurrency` at once), prints it, and calls `sbatch`. Run a step directly (without `submit`)
  to test one sample first.
- Steps skip a sample whose final output already exists; pass `--force` to redo it. Trim, align and quant
  need the previous stage to be finished: wait for one array to complete before submitting the next.
- Each task uses its own `TMPDIR` below `scratch_dir`, so concurrent tasks do not collide. The scripts also
  raise the open-file limit to 10000 (STAR needs this for sorting).
- Alignment (STAR) and quantification (Salmon) both start from the FASTQ files and are independent
  of each other; run either, or both.

## 4. Things that commonly go wrong

- **Paired vs single-end.** Single-end data needs `-r` in Salmon and one file in STAR, and
  mixing the two up fails or silently mis-pairs reads. Leave `fastq_2` out of the sample entry for
  single-end libraries. Make sure `fastq_1` and `fastq_2` are really R1 and R2 (the scripts only stop
  if they are the same file).
- **Index and version.** STAR refuses an index built with a different version. The STAR and Salmon
  indexes must come from the same genome build and annotation as each other and as any GTF you use later,
  or gene IDs will not match.
- **Strandedness.** Salmon is run with `-l A` (auto-detect). Check `salmon/<sample>/lib_format_counts.json`
  to see what it found. STAR is run without any strand-specific option, so check the protocol
  before counting genes from the BAMs.
- **Memory vs `limit_bam_sort_ram`.** The original runs used `--limitBAMsortRAM 214012628346` (about 200 GB)
  on a large-memory node but requested a much smaller SLURM `--mem`. Set `limit_bam_sort_ram` (bytes)
  at or below your job memory, or leave it unset so STAR picks a default.
  Large human genomes need roughly 30-40 GB of RAM for STAR itself.
- **Unique-mapping runs.** One of the original single-end runs used stricter settings. To reproduce them,
  set `"star": {"sam_attributes": ["Standard", "NM"], "extra_args": ["--outFilterMismatchNmax", "5",
  "--outFilterMultimapNmax", "1", "--winAnchorMultimapNmax", "1"]}`.
- **Shell-script typos the Python version avoids.** In the originals, a missing trailing backslash before
  `--limitBAMsortRAM` made that flag run as its own (failing) command; one output prefix did not match its
  `mkdir`d folder; and a Salmon sample list used `<id>.fastq` for R1 but `<id>_2.fastq` for R2.
  Here, R1/R2 are named explicitly per sample and inputs are checked before a tool starts.
- **Reads from different stages.** Trimmed files are always named `<sample>_1.fq.gz` / `_2.fq.gz`
  (single-end: `_1` only), no matter how the raw files were named.
- **Full node scratch.** If node-local `/tmp` is full, set `scratch_dir` to a larger disk.
- **Java.** FastQC needs Java. The conda package installs it; on a cluster with modules, load Java first.
- **Re-running.** Because finished samples are skipped, a corrupt or truncated earlier output (for example
  after a killed job) will also be skipped. Delete the sample's output folder or use `--force`.
