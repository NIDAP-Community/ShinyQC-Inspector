# ShinyQC Inspector

ShinyQC Inspector explores expression-PCA structure alongside technical QC metrics.
It provides interactive PCA panels, ranks samples by QC/PCA atypicality, and can
export HTML reports. Batch utilities generate the same report format for multiple
cohorts, such as TCGA.

## Quick start

Use R 4.5.1 (the version recorded in `renv.lock`) and restore the locked project
environment:

```r
renv::restore()
shiny::runApp()
```

The local `renv/library` directory is intentionally ignored. A new checkout
therefore needs outbound package access for `renv::restore()` before the app is
run.

The app accepts a normalized expression matrix with a `Gene` column, sample
metadata, and QC metadata. The selected sample-ID columns must match the count
matrix sample names after the app's underscore-to-hyphen normalization.

## Batch reports

Create a manifest, then render per-cohort reports:

```sh
Rscript scripts/create_tcga_manifest.R /path/to/tcga_root /path/to/tcga_manifest.tsv
Rscript scripts/batch_qc_reports.R /path/to/tcga_manifest.tsv /path/to/report_output
Rscript scripts/summarize_tcga_qc_pc_axes_from_html.R /path/to/report_output /path/to/summary_output
```

`batch_qc_reports.R --help` describes the required manifest columns. The detailed
analysis design and its provisional thresholds are in
[`docs/tcga_qc_bias_scoring_plan.md`](docs/tcga_qc_bias_scoring_plan.md).

## Validation

Run the repository smoke test from the project root:

```sh
Rscript tests/run_smoke_tests.R
```

It generates a temporary synthetic cohort, renders a report, summarizes the
resulting report, and checks that the Shiny app and helper functions load. This
does not replace review with real cohort data.
