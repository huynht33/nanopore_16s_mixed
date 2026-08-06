# Nanopore full-length 16S analysis of a mixed anaerobic consortium

This repository contains the reproducible workflow, reference sequences,
metadata, derived abundance tables, and R Markdown notebook used to estimate
relative abundance in a five-species anaerobic consortium from Oxford Nanopore
full-length 16S amplicon reads.

The five target species, in the fixed plotting order, are:

1. *Bacteroides thetaiotaomicron*
2. *Cutibacterium acnes*
3. *Fusobacterium nucleatum*
4. *Prevotella melaninogenica*
5. *Veillonella parvula*

## Analysis design

The workflow uses three complementary analyses:

- **Custom Emu:** primary targeted abundance analysis against one
  Plasmidsaurus Linear/PCR consensus 16S reference per expected species.
- **Broad Emu:** secondary screening against a published broad Emu database to
  identify unexpected taxa and reads not accounted for by the five targets.
- **minimap2:** transparent read-level cross-check against the five custom
  references. Reads are reported as confident, ambiguous, or unassigned rather
  than forcing ambiguous alignments to one species.

Two custom-Emu percentages are retained. The all-usable-read percentage keeps
unmapped and filtered reads in the denominator and is the QC-aware result. The
renormalized percentage uses only reads assigned to the five expected species
and is the targeted consortium-composition result. The notebook labels these
denominators explicitly.

## Repository contents

```text
config/                         Configurable filtering and assignment thresholds
environments/                   Conda specification and recorded versions
fastq/raw/                      Place raw FASTQ files here (not tracked by Git)
metadata/                       Sample keys, exclusions, and expected ratios
references/clean/               Five cleaned single-record reference FASTAs
references/combined/            Combined FASTA, taxids, checksums, and verification
results/combined/               Derived machine-readable analysis tables
results/combined_rmarkdown/     Tables exported by the notebook
results/figures_rmarkdown/      PNG and PDF figures exported by the notebook
scripts/                        Numbered Bash, Python, and R workflow scripts
nano16s_emu_figures.Rmd         Inspectable RStudio analysis and plotting notebook
```

Raw/filtered FASTQs, Emu databases, SAM/PAF alignments, logs, NanoPlot reports,
and RStudio state are intentionally excluded. They are large, regenerable, or
environment-specific. No script edits an original FASTQ or original reference
FASTA.

## Reproduce the figures in RStudio

The compact derived inputs needed by the notebook are tracked, including the
separately processed sample 16 baseline. No FASTQ or Emu database is needed for
this route.

```bash
git clone https://github.com/huynht33/nanopore_16s_mixed.git
cd nanopore_16s_mixed
```

Open `nano16s_emu_figures.Rmd` in RStudio. Install the required packages once:

```r
install.packages(c(
  "rmarkdown", "knitr", "readr", "dplyr", "tidyr", "ggplot2",
  "stringr", "scales", "purrr"
))
```

Then click **Knit**, run chunks individually, or render from R:

```r
rmarkdown::render("nano16s_emu_figures.Rmd")
```

The notebook writes tables to `results/combined_rmarkdown/` and figures to
`results/figures_rmarkdown/`. Recorded package versions are in
`reports/R_sessionInfo.txt` and `environments/r_package_versions.tsv`.

## Reproduce from FASTQ files

### 1. Create the software environment

Do not install into Conda base.

```bash
conda env create -f environments/emu16s_environment.yml
conda activate nanopore16s
```

The environment pins the principal versions used for the analysis: Emu 3.6.2,
minimap2 2.31, seqkit 2.13.0, samtools 1.24, and Python 3.12. The scripts inspect
installed command help before using version-dependent Emu options.

### 2. Supply a broad Emu database

The broad database is not stored in Git. It must contain at least
`species_taxid.fasta` and `taxonomy.tsv` (or the alternate metadata layout
accepted by the installed Emu version). Copy the example resource config:

```bash
cp config/local_resources.sh.example config/local_resources.sh
```

Edit `BROAD_EMU_DB` in the copied file. The original analysis reused a legacy
broad database from the preceding pure-isolate project. Its upstream release
and creation date were not recorded, so broad-database assignments should be
treated as database-version dependent. Obtain and document a published Emu
database before reproducing that secondary analysis.

### 3. Add FASTQs without changing them

Place any combination of `.fastq`, `.fastq.gz`, `.fq`, or `.fq.gz` files under
`fastq/raw/`. Files are discovered recursively and sample IDs are derived from
filenames after removing the FASTQ extension.

```bash
mkdir -p fastq/raw
rsync -av --progress /path/to/fastqs/ fastq/raw/
```

The analysis originally expected 20 files. To intentionally analyze another
count, review the manifest and set `ALLOW_NON_20_FASTQ=1`.

### 4. Validate references and build the custom database

The five cleaned reference FASTAs and verified NCBI taxids are tracked. Running
step 01 validates the sequences again and resolves each species exactly against
the supplied broad database taxonomy. It stops on missing or ambiguous matches.

```bash
export PROJECT_ROOT="$PWD"
scripts/01_prepare_references.sh
scripts/02_build_custom_emu_database.sh
```

Verified mappings used here are taxids 818, 1747, 851, 28132, and 29466,
respectively. They remain reviewable in
`references/combined/taxonomy_verification.tsv`; the scripts do not silently
hardcode them during a rebuild.

### 5. Inventory FASTQs and review metadata

```bash
scripts/03_inventory_fastq.sh
```

Review the generated `metadata/fastq_manifest.tsv`, edit
`metadata/sample_metadata.tsv`, confirm `metadata/sample_exclusions.tsv`, and
edit expected ratios where known. The four validation mixtures use CFU-adjusted
percentages calculated from input volumes and measured strain concentrations,
not unadjusted volume percentages. When review is complete:

```bash
touch metadata/.metadata_reviewed
```

### 6. Run the controlled workflow

```bash
scripts/run_pipeline.sh --project-root "$PWD" --run-id "$(date +%Y%m%d_%H%M%S)"
```

The runner stops on missing validation flags, absent FASTQs, duplicate sample
IDs, or an existing run ID. It never deletes previous outputs. The configured
defaults are 1,200-1,800 bp, mean read quality Q10 when Chopper is available,
eight threads, and Emu minimum PID 80%. Thresholds live in
`config/project_config.sh`.

Sample 12 was excluded for having only four short reads. Sample 16 had an
atypical raw length distribution and was therefore excluded from the original
cross-sectional run, then analyzed separately after the standard filters
retained 1,977 near-full-length reads. To reproduce a reviewed single-sample
exploration, see:

```bash
scripts/09_explore_single_sample.sh --help
```

## Known-ratio interpretation

Samples 4-7 are the validation mixtures. Their expected values are based on
input volume multiplied by the measured CFU/mL for each strain. CFU-input
fractions are useful biological comparators but are not expected to equal 16S
read fractions exactly because extraction efficiency, 16S copy number, primer
matching, PCR amplification, and sequencing/assignment can differ by organism.

The notebook reports mean absolute error in percentage points, RMSE,
correlations, organism-specific recovery, systematic bias, and Bray-Curtis
similarity. Samples with known absent species are used to estimate empirical
false-positive backgrounds. No universal detection threshold is declared in
advance.

Relative abundance is compositional: a decrease in one organism can make the
others' percentages increase even if their absolute biomass is unchanged.
Longitudinal P1-normalized plots therefore describe change in read composition,
not absolute CFU or cell abundance.

## Reproducibility notes

- All executable workflow Bash scripts use `set -euo pipefail`.
- Output and command logs are created beneath `logs/` and per-sample result
  directories.
- Original FASTQs and original reference FASTAs are never overwritten.
- Custom, broad, and minimap2 outputs occupy separate directories.
- Reference lengths, GC content, N counts, sequence MD5s, and taxid checks are
  recorded under `references/combined/`.
- The broad Emu database is the only major external resource whose exact legacy
  release was unavailable; record its path, date, and checksums for every new
  run.

## Software and citation

The command-line analysis used Emu v3.6.2 with the Nanopore `map-ont` preset and
minimap2 v2.31-r1302. See `environments/software_versions.tsv` for the complete
record. Emu is described in Curry et al., *Nature Methods* (2022),
<https://doi.org/10.1038/s41592-022-01520-4>.
