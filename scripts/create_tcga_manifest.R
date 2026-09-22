#!/usr/bin/env Rscript

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/create_tcga_manifest.R tcga_root manifest_output.tsv [selected_vars] [context_vars]\n\n",
    "Example:\n",
    "  Rscript scripts/create_tcga_manifest.R \\\n",
    "    /Volumes/SEQC2/IODC/IODC-0001-TCGA \\\n",
    "    /Users/maggiec/GitHub/Maggie/IODC/ShinyQC_Output/tcga_manifest.tsv\n\n",
    "Use selected_vars=ALL_QC, or omit selected_vars, to include every QC column.\n",
    "Use context_vars=ALL_CONTEXT to include every sample metadata column, or omit context_vars to use the default clinical/pathology context set.\n",
    sep = ""
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 4 || any(args %in% c("-h", "--help"))) {
  usage()
  quit(status = ifelse(length(args) %in% 2:4, 0, 1))
}

tcga_root <- args[[1]]
manifest_output <- args[[2]]
selected_vars <- if (length(args) == 3) {
  args[[3]]
} else {
  "ALL_QC"
}
default_context_vars <- c(
  "gender",
  "age_at_initial_pathologic_diagnosis",
  "ajcc_pathologic_tumor_stage",
  "clinical_M",
  "pathologic_T",
  "pathologic_N",
  "tumor_status",
  "new_tumor_event_dx_indicator",
  "treatment_outcome_first_course",
  "residual_tumor",
  "histologic_diagnosis",
  "weiss_score_overall",
  "nuclear_grade_III_IV",
  "necrosis",
  "mitotic_rate",
  "mitoses_per_50_hpf",
  "invasion_of_tumor_capsule",
  "history_adrenal_hormone_excess",
  "tissue_source_site"
)
context_vars <- if (length(args) == 4) {
  args[[4]]
} else {
  paste(default_context_vars, collapse = ",")
}

if (!dir.exists(tcga_root)) {
  stop("TCGA root does not exist: ", tcga_root)
}

cohort_dirs <- list.dirs(tcga_root, full.names = TRUE, recursive = FALSE)
cohort_dirs <- cohort_dirs[grepl("IODC-0001-TCGA-", basename(cohort_dirs))]

find_one <- function(cohort_dir, pattern) {
  matches <- list.files(cohort_dir, pattern = pattern, full.names = TRUE)
  if (length(matches) == 0) {
    return(NA_character_)
  }
  matches[[1]]
}

manifest <- data.frame(
  counts_file = character(),
  sample_meta_file = character(),
  qc_file = character(),
  sample_meta_col = character(),
  qc_col = character(),
  selected_vars = character(),
  context_vars = character(),
  report_name = character()
)

for (cohort_dir in cohort_dirs) {
  cohort <- basename(cohort_dir)
  counts_file <- find_one(cohort_dir, "\\.normalized_counts\\.tsv$")
  sample_meta_file <- find_one(cohort_dir, "\\.meta\\.sample\\.tsv$")
  qc_file <- find_one(cohort_dir, "\\.qc\\.tsv$")

  if (any(is.na(c(counts_file, sample_meta_file, qc_file)))) {
    warning("Skipping ", cohort, ": missing normalized_counts, meta.sample, or qc file")
    next
  }

  manifest <- rbind(
    manifest,
    data.frame(
      counts_file = counts_file,
      sample_meta_file = sample_meta_file,
      qc_file = qc_file,
      sample_meta_col = "sample",
      qc_col = "Sample",
      selected_vars = selected_vars,
      context_vars = context_vars,
      report_name = paste0(cohort, "_QC_report")
    )
  )
}

dir.create(dirname(manifest_output), showWarnings = FALSE, recursive = TRUE)
write.table(manifest, manifest_output, sep = "\t", row.names = FALSE, quote = FALSE)

cat("Wrote", nrow(manifest), "manifest rows to", manifest_output, "\n")
