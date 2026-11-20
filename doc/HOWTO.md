# How to use segdup-caller (open-source fork)

This guide covers installing and running the fork that replaces Sentieon with GATK, minimap2, samtools and bcftools. For what each gene module does and how to read its results, see [general.md](general.md) and the per-gene documents.

> **Research use only.** The tool is not validated for clinical use. This fork also changes the small-variant caller and the copy-number caller compared with upstream, and has **not been re-validated**.

## 1. What you need

| Item | Requirement |
|---|---|
| Tools on `PATH` | GATK >= 4.2, minimap2 >= 2.24, samtools >= 1.16, bcftools >= 1.16 |
| Python | >= 3.9 (installed by the package: pysam, pyyaml, colorlog, packaging, scipy, pandas, vcflib) |
| Java | Needed by GATK (it must be available to GATK) |
| Reference | GRCh38, hg19 or b37 FASTA, with `.fai` and a GATK `.dict` |
| Input reads | Aligned, coordinate-sorted, **duplicate-marked**, indexed BAM/CRAM from whole-genome sequencing, 20x minimum, 30x ideal |
| Long reads (optional) | Aligned PacBio HiFi or ONT BAM/CRAM |

No Sentieon software, license, `SENTIEON_LICENSE` variable or model bundles are needed.

The reference build is detected automatically from the FASTA. Contig names must match the BAM header.

## 2. Install

```bash
git clone <your fork url> segdup-caller
cd segdup-caller
conda env create -f environment.yml     # python and the package only; gatk, minimap2, samtools, bcftools must already be on PATH
conda activate segdup-caller
segdup-caller --version
```

Without conda: `pip install ./segdup-caller`. In both cases the four tools are not installed for you.

Check the tools are found and recent enough:

```bash
gatk --version; minimap2 --version; samtools --version | head -1; bcftools --version | head -1
```

## 3. Prepare the reference (once)

```bash
samtools faidx GRCh38.fa
samtools dict GRCh38.fa -o GRCh38.dict      # GATK looks for GRCh38.dict next to GRCh38.fa
```

The caller exits at startup, with the exact command to run, if the `.dict` is missing. A CRAM input needs the same reference it was written with.

## 4. Prepare the reads

Reads are **not** aligned by this tool. Align and mark duplicates first, for example:

```bash
bwa-mem2 mem -t 16 -R '@RG\tID:s1\tSM:NA12878\tPL:ILLUMINA' GRCh38.fa r1.fq.gz r2.fq.gz \
  | samtools fixmate -m - - | samtools sort -@8 - \
  | samtools markdup - NA12878.bam
samtools index NA12878.bam
```

Requirements that the caller checks or relies on:

- the index (`.bai` for BAM, `.crai` for CRAM) exists next to the file;
- the header has an `@RG` line, with an `SM` tag unless you pass `--sample_name`; HaplotypeCaller needs read groups;
- duplicates are flagged (0x400). Flagged duplicates are dropped during liftover re-alignment.

## 5. Run it

Minimal run, all supported genes:

```bash
segdup-caller -s NA12878.bam -r GRCh38.fa -o results/
```

Selected genes only:

```bash
segdup-caller -s NA12878.bam -r GRCh38.fa -g SMN1,CYP2D6,GBA -o results/
```

Short plus long reads:

```bash
segdup-caller -s NA12878.sr.bam -l NA12878.hifi.bam -r GRCh38.fa -o results/
```

One merged VCF in addition to the per-gene VCFs:

```bash
segdup-caller -s NA12878.bam -r GRCh38.fa -o results/ --merge_vcf
```

### Options

| Option | Meaning |
|---|---|
| `-s, --short` | Short-read BAM/CRAM (required) |
| `-l, --long` | Long-read BAM/CRAM (optional) |
| `-r, --reference` | Reference FASTA (required) |
| `-o, --outdir` | Output directory (required) |
| `-g, --genes` | Comma-separated genes. Default: all in the config. `segdup-caller -h` lists them |
| `--sample_name` | Sample name. Default: `SM` tag of the short-read BAM |
| `--sr_prefix`, `--lr_prefix` | Prefix for short/long result files (defaults `sr`, `lr`) |
| `-t, --threads` | Total threads, capped at the CPU count |
| `-w, --workers` | Genes processed in parallel (default 4) |
| `--sex` | `male`, `female`, `M`, `F`, `XY`, `XX`; affects X/Y-linked genes |
| `--merge_vcf` | Also write `<sample>.merged.vcf.gz` |
| `--keep_temp` | Keep the temporary directory (useful for debugging) |
| `--config FILE` | Use your own gene YAML instead of the shipped one (replaces it entirely) |
| `--set KEY.PATH=VALUE` | Override one scalar in the config; repeatable |
| `--gatk-java-options "..."` | JVM options for every GATK call. Default `-Xmx4g -XX:ParallelGCThreads=1` |
| `--sr_model`, `--lr_model` | Accepted and **ignored** (old Sentieon options) |

### Threads and workers

`--workers` genes run at once, and at most `--threads` tool processes run concurrently across all workers. Each GATK call uses one pair-HMM thread and by default 4 GB of heap, so memory is roughly `workers x 4 GB`.

```bash
# 16-core machine, 64 GB RAM
segdup-caller -s s.bam -r ref.fa -o out -w 4 -t 16
# small machine
segdup-caller -s s.bam -r ref.fa -o out -w 2 -t 4 --gatk-java-options "-Xmx2g"
```

### Tuning with `--set`

Keys are paths into the gene YAML (`genecaller/data/<build>/genes.yaml`). The `main:` section holds global settings:

```bash
--set main.min_map_qual=40        # MAPQ threshold for confident segments
--set main.cnv_bin_size=500       # depth-caller bin width in bp (default 1000)
--set main.gatk_min_mapq=20       # HaplotypeCaller --minimum-mapping-quality (default 10)
```

Only scalar values can be overridden this way; copy the YAML and use `--config` for larger changes.

## 6. Outputs

In `--outdir`:

- a **YAML summary**: copy-number states, gene-conversion events, star alleles or genotypes where a gene module reports them, run options, resource use and software versions (segdup-caller, gatk, minimap2, samtools, bcftools);
- **per-gene VCFs** (bgzipped, tabix-indexed), including phased VCFs;
- optionally `<sample>.merged.vcf.gz` with an `INFO/GENE` tag on every record;
- a `<sample>_tmp...` directory, removed unless `--keep_temp` is set.

Clinical-significance text in the YAML is informational only.

## 7. How it works, briefly

1. Per gene, reads from the gene and its paralog regions are extracted (samtools) and re-mapped onto the paralog (minimap2) to "lift over" ambiguous reads.
2. Copy numbers come from a **depth-based caller** (`calling/cnv.py`): `samtools bedcov` over 1 kb bins, normalised to the median depth of well-mapped autosomal bins, combined with allele-depth evidence in the gene models. This replaces Sentieon CNVscope. With `--cnv_caller cnvkit` the same bins are GC/edge/repeat-corrected and segmented by [CNVkit](https://cnvkit.readthedocs.io) (`cnvkit.py` >= 0.9.12 on PATH, flat reference, Haar segmentation) and re-anchored to the sample's genome-wide diploid depth; the default stays the built-in caller.
3. Small variants are called per region and per copy number with `gatk HaplotypeCaller --sample-ploidy N`; with known variants available, genotyping is given-alleles.
4. Variants are phased by the built-in k-modes (short reads) or anchor (long reads) phasers.

## 8. Differences from upstream and known limits

- No ML rescoring of variant calls (DNAscope + model). Expect more false positives and negatives in difficult segdup sequence than upstream; filter with GATK annotations if needed.
- The default copy-number caller is a simpler heuristic: no GC or mappability correction (`--cnv_caller cnvkit` adds it; that backend is also not benchmarked against CNVscope or the default). Low-coverage or uneven-coverage samples will be less reliable, and the calls were not benchmarked against CNVscope.
- Long-read small-variant calling uses HaplotypeCaller, which performs poorly on long reads. Treat long-read VCFs with caution.
- Accuracy and concordance figures in `general.md` and the gene documents are upstream's.

## 9. Troubleshooting

| Message or symptom | Cause and fix |
|---|---|
| `no 'gatk' found in the PATH` | Make sure GATK is on `PATH` (load the module or add it to `PATH`). Same for `minimap2`, `samtools`, `bcftools` |
| `the pipeline requires gatk version '4.2' or later` | Upgrade the tool |
| `Reference sequence dictionary does not exist` | Run `samtools dict ref.fa -o ref.dict` |
| `Index file ... .bai does not exist` | `samtools index` the input |
| `SM tag not found` | Add an `SM` to `@RG`, or pass `--sample_name` |
| `No gene configuration is shipped for reference build` | The build was not recognised. Provide `--config` |
| GATK `Contig ... not found` / `contigs differ` | BAM and reference use different contig names (`chr1` vs `1`) or builds |
| `Too few well-mapped bins ... diploid depth` | Very low coverage or a reference whose autosomes are not recognisable by name |
| Out-of-memory / killed | Lower `--workers`, or lower `--gatk-java-options` heap per call |
| Need to see commands | Run with `--keep_temp` and check the log; commands are logged at debug level (`--log_level debug`) |

Each run creates a fresh temporary directory under `--outdir`, so a failed run is rerun from the start.

## 10. Testing the installation

```bash
pip install pytest
pytest tests -q
```

The tests use a small simulated genome and need only samtools and bcftools (minimap2 for one test); GATK is stubbed. They do **not** prove that GATK integration works on real data. Do a real run on a reference sample, such as a GIAB genome at 30x, over a few genes (`-g GBA,SBDS,CYP2D6`), and check that the VCFs are produced and the copy numbers match what you expect.
