# ShinyQC Inspector handoff

## Current scope

This branch adds a cohort-aware QC/PCA reporting workflow for the ShinyQC
Inspector. It is intended to prioritize samples and technical QC axes for
review; it is not yet an automated sample-exclusion policy.

### Included work

- Interactive Shiny tables for QC/PCA outlier ranking and QC metric redundancy.
- HTML report export from the app.
- Batch report generation from a cohort manifest.
- TCGA manifest creation and two cross-cohort summary utilities.
- An `renv.lock` environment declaration and a data-free smoke test.
- The analysis plan in `tcga_qc_bias_scoring_plan.md`.

## Run and validate

```sh
Rscript tests/run_smoke_tests.R
```

The smoke test creates temporary synthetic counts, sample metadata, and QC
metadata; renders an HTML report; summarizes that report; and loads the Shiny
app's helper functions. It does not test a live browser session or replace an
analysis with real data.

The lockfile is committed but the local `renv/library` and package cache are
not. A new checkout must run `renv::restore()` before running the app or smoke
test. During this handoff, code-level checks passed from installed packages,
but a full locked-environment restore was not completed because a Posit Package
Manager binary download stalled after redirecting. Treat the restore as pending
until it completes successfully on a network-enabled machine.

For a real cohort, follow the commands in the README. Raw TCGA inputs are not
included in this repository. At handoff, the expected `/Volumes/SEQC2` source
mount was unavailable, so a current full-cohort rerun has not been performed.

## Decisions still required

- Confirm expression filtering, normalization, and transformation before PCA.
- Review ACC first, then set the QC-correlation grouping threshold and
  cohort-specific score thresholds.
- Keep score output as ranked QC atypicality until independent evidence supports
  any `Pass`, `Review`, or `Fail` classification.
- Define how samples with biological and technical associations on the same PC
  will be reported as mixed or confounded rather than automatically excluded.

## Review focus

Review the score-construction rule, missing-data handling, sample-ID alignment,
and the interpretation of technical versus biological PC associations. The
existing reports may be summarized from HTML, but they are historical outputs;
reproducibility of their source data requires the external TCGA mount.

## Repository hygiene

`.Rproj.user`, `.posit`, `rsconnect`, and analysis outputs are local state and
are deliberately excluded from version control. The project source, scripts,
documentation, and `renv` lockfile are the handoff artifacts.
