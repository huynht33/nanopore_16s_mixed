This directory stores the reproducible software specification and recorded
versions for the command-line and R analyses.

Create a separate environment; do not install packages into Conda base:

```bash
conda env create -f environments/emu16s_environment.yml
conda activate nanopore16s
```

`software_versions.tsv` records command-line versions observed on Superdome.
`r_package_versions.tsv` and `../reports/R_sessionInfo.txt` record the local
RStudio figure-analysis session. NanoPlot and Chopper were present and used, but
their exact versions were not captured in the synchronized logs; they are left
unpinned in the environment specification instead of inventing a version.
