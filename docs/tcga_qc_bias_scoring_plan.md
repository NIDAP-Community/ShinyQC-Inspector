# TCGA QC Bias Scoring Plan

## Goal

Develop a cohort-aware QC scoring workflow that identifies samples whose expression PCA position is plausibly driven by technical QC biases, while avoiding over-counting correlated QC metrics.

The workflow should produce:

- Per-cohort HTML reports.
- Per-cohort sample-level QC scores and suggested `Pass`, `Review`, or `Fail` flags.
- Per-cohort QC metric association summaries for PC1 and PC2.
- Per-cohort clinical/pathology context summaries for PC1 and PC2.
- Per-cohort QC redundancy groups that collapse correlated metrics into interpretable technical axes.
- A combined cross-cohort summary of recurring and cohort-specific technical bias patterns.

## Key Assumptions

- Technical bias patterns may differ by TCGA cohort.
- Expression PCs may reflect real clinical, pathological, or biological structure as well as technical effects.
- PCA sign is arbitrary, so scoring should use absolute associations with PC1, PC2, and PCA distance rather than hard-coded quadrant labels.
- Multiple QC metrics may represent the same underlying technical process, so correlated metrics should be grouped before score construction.
- Clinical/pathology variables should contextualize PC interpretation but should not contribute to technical QC scores.
- Final QC flags should be cohort-specific, not learned from one cohort and blindly applied to all cohorts.

## Per-Cohort Workflow

For each TCGA cohort:

1. Read normalized counts, sample metadata, and QC metadata.
2. Align sample IDs across expression, sample metadata, and QC metadata.
3. Run PCA on filtered expression data.
4. Quantify association between each QC metric and PC1/PC2.
5. Quantify association between selected clinical/pathology context variables and PC1/PC2.
6. Identify QC metrics that strongly track with PC1 or PC2.
7. Compute QC metric redundancy groups from QC-QC correlations.
8. Select representative metrics for each redundancy group.
9. Build a cohort-specific QC score using representative metrics.
10. Assign sample-level `Pass`, `Review`, or `Fail` labels.
11. Save the per-cohort HTML report.

## QC Association Analysis

Numeric QC metrics:

- Use Spearman correlation against PC1 and PC2.
- Track `rho`, p-value, and FDR-adjusted p-value.
- Candidate inclusion rule:
  - `abs(rho) >= 0.4`
  - FDR `< 0.05`

Categorical QC metrics:

- Use Kruskal-Wallis testing of PC1 and PC2 by category.
- Report an effect-size style statistic.
- Treat categorical metrics as annotations or stratification variables unless a clear scoring rule is defined.

## Clinical/Pathology Context Analysis

Clinical, pathological, and biological sample metadata should be tested in parallel with QC metrics to help interpret what PC1 and PC2 represent.

Candidate context variables include:

- Demographics such as `gender` and `age_at_initial_pathologic_diagnosis`.
- Stage and disease status such as `ajcc_pathologic_tumor_stage`, `clinical_M`, `pathologic_T`, `pathologic_N`, `tumor_status`, and `residual_tumor`.
- Treatment/outcome context such as `new_tumor_event_dx_indicator` and `treatment_outcome_first_course`.
- Pathology severity fields such as `weiss_score_overall`, `nuclear_grade_III_IV`, `necrosis`, `mitotic_rate`, `mitoses_per_50_hpf`, and `invasion_of_tumor_capsule` where available.
- Cohort-specific biology such as `history_adrenal_hormone_excess` in ACC.
- Collection context such as `tissue_source_site`.

Use the same PC1/PC2 association methods as the QC association analysis:

- Numeric context variables use Spearman correlation.
- Categorical context variables use Kruskal-Wallis effect size.
- Missing metadata placeholders such as `_Not_Available_`, `_Unknown_`, and `_Not_Applicable_` should be treated as missing.

These variables are reported separately from QC metrics. They can support labels such as `likely biological`, `likely technical`, `mixed/confounded`, or `unclear`, but they should not be included in the QC outlier score or the QC redundancy grouping.

## QC Redundancy Analysis

For numeric QC metrics:

1. Compute pairwise Spearman correlations.
2. Use absolute correlation for redundancy:
   - candidate threshold: `abs(rho) >= 0.7` or `0.8`
3. Cluster correlated metrics into redundancy groups.
4. Pick one representative metric per group using:
   - strongest PC1/PC2 association,
   - interpretability,
   - consistency across related cohorts.

This prevents over-counting one technical process just because it appears in several correlated QC fields.

## Sample Scoring

For each cohort:

1. Use representative numeric QC metrics from PC-associated redundancy groups.
2. Compute robust z-scores per metric.
3. Weight each metric by its strongest absolute association with PC1 or PC2.
4. Sum weighted absolute robust z-scores:

```text
QC score = sum(abs(robust_z(metric)) * strongest_abs_PC_association)
```

Candidate sample labels:

```text
Pass   score < 3
Review score 3 to 6
Fail   score > 6
```

These thresholds are provisional and should be reviewed after looking across cohorts.

## Cross-Cohort Summary

After per-cohort reports are generated, create a combined table:

```text
cohort | strongest_QC_axis | representative_metrics | flagged_samples | notes
```

Examples of likely technical axes:

- Read length / insert distance.
- RNA integrity / coverage bias.
- 3-prime / 5-prime bias.
- Base allocation / genomic feature composition.
- Strandedness or sense/antisense balance.
- Sequencing center or flowcell effects.

## Report Additions

Each cohort HTML report should include:

- Top QC outlier scores.
- QC associations with PC1 and PC2.
- Clinical/pathology context associations with PC1 and PC2.
- QC redundancy groups.
- Representative metrics used for scoring.
- Sample-level `Pass`, `Review`, or `Fail` table.
- PCA panels for QC metrics.
- PCA panels for clinical/pathology context variables.

The combined output folder should remain:

```text
/Users/maggiec/GitHub/Maggie/IODC/ShinyQC_Output
```

## Execution Checkpoint

Before running all cohorts again:

1. Review this plan.
2. Confirm scoring thresholds.
3. Confirm correlation threshold for redundancy grouping.
4. Confirm whether categorical metrics should contribute to scores or only annotate plots.
5. Regenerate ACC first.
6. Review ACC output format.
7. Run remaining cohorts after approval.
