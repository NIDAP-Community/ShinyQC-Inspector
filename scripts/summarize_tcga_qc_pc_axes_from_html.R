#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(rvest)
})

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/summarize_tcga_qc_pc_axes_from_html.R report_dir output_dir\n\n",
    "Reads local ShinyQC HTML reports and summarizes QC metrics significantly associated with PC1/PC2.\n\n",
    "Outputs:\n",
    "  tcga_qc_pc_axis_all_associations_from_html.tsv\n",
    "  tcga_qc_pc_axis_significant_from_html.tsv\n",
    "  tcga_qc_pc_axis_parameter_summary_from_html.tsv\n",
    "  tcga_qc_pc_axis_parameter_summary_from_html.md\n",
    sep = ""
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || any(args %in% c("-h", "--help"))) {
  usage()
  quit(status = ifelse(length(args) == 2, 0, 1))
}

report_dir <- args[[1]]
output_dir <- args[[2]]

report_files <- list.files(
  report_dir,
  pattern = "^IODC-0001-TCGA-.*_QC_report[.]html$",
  full.names = TRUE
)

if (length(report_files) == 0) {
  stop("No TCGA report HTML files found in: ", report_dir)
}

cohort_from_file <- function(file_path) {
  gsub("_QC_report[.]html$", "", basename(file_path))
}

to_logical <- function(value) {
  value <- tolower(trimws(as.character(value)))
  value %in% c("true", "t", "1", "yes")
}

format_cohort_list <- function(values, limit = 12) {
  values <- sort(unique(values))
  if (length(values) == 0) {
    return("")
  }
  if (length(values) <= limit) {
    return(paste(values, collapse = ", "))
  }
  paste(paste(values[seq_len(limit)], collapse = ", "), "...")
}

safe_max <- function(values) {
  if (length(values) == 0 || all(is.na(values))) {
    return(NA_real_)
  }
  max(values, na.rm = TRUE)
}

safe_best_cohort <- function(cohorts, values) {
  if (length(values) == 0 || all(is.na(values))) {
    return(NA_character_)
  }
  cohorts[which.max(ifelse(is.na(values), -Inf, values))][[1]]
}

markdown_table <- function(data, columns) {
  data <- data[, columns, drop = FALSE]
  header <- paste(names(data), collapse = " | ")
  separator <- paste(rep("---", length(columns)), collapse = " | ")
  body <- apply(data, 1, function(row) paste(row, collapse = " | "))
  paste(c(header, separator, body), collapse = "\n")
}

read_association_table <- function(file_path) {
  html <- read_html(file_path)
  tables <- html_table(html, fill = TRUE)

  if (length(tables) < 2) {
    stop("Expected at least two tables in report: ", file_path)
  }

  association_table <- tables[[2]]
  required_cols <- c(
    "metric",
    "type",
    "n",
    "PC1_assoc",
    "PC1_FDR",
    "PC2_assoc",
    "PC2_FDR",
    "strongest_PC",
    "strongest_abs_assoc",
    "include_in_score"
  )
  missing_cols <- setdiff(required_cols, names(association_table))
  if (length(missing_cols) > 0) {
    stop(
      "Association table in ",
      basename(file_path),
      " is missing columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  association_table %>%
    mutate(
      cohort = cohort_from_file(file_path),
      n = as.integer(n),
      PC1_assoc = as.numeric(PC1_assoc),
      PC1_FDR = as.numeric(PC1_FDR),
      PC2_assoc = as.numeric(PC2_assoc),
      PC2_FDR = as.numeric(PC2_FDR),
      strongest_abs_assoc = as.numeric(strongest_abs_assoc),
      include_in_score = to_logical(include_in_score)
    ) %>%
    select(
      cohort,
      metric,
      type,
      n,
      PC1_assoc,
      PC1_FDR,
      PC2_assoc,
      PC2_FDR,
      strongest_PC,
      strongest_abs_assoc,
      include_in_score
    )
}

all_associations <- bind_rows(lapply(report_files, read_association_table)) %>%
  mutate(
    strongest_PC_FDR = pmin(PC1_FDR, PC2_FDR, na.rm = TRUE),
    magnitude_threshold = ifelse(type == "numeric", 0.4, 0.1),
    significant_PC_axis = strongest_abs_assoc >= magnitude_threshold & strongest_PC_FDR < 0.05,
    significance_rule = ifelse(
      type == "numeric",
      "abs(Spearman rho) >= 0.4 and FDR < 0.05",
      "Kruskal-Wallis effect >= 0.1 and FDR < 0.05"
    )
  )

significant_associations <- all_associations %>%
  filter(significant_PC_axis) %>%
  arrange(metric, cohort)

parameter_summary <- all_associations %>%
  group_by(metric, type) %>%
  summarise(
    n_cohorts_tested = n_distinct(cohort),
    n_cohorts_significant_PC_axis = n_distinct(cohort[significant_PC_axis]),
    n_cohorts_in_score = n_distinct(cohort[include_in_score]),
    best_abs_PC_axis_assoc = safe_max(strongest_abs_assoc),
    best_PC_axis_cohort = safe_best_cohort(cohort, strongest_abs_assoc),
    significant_cohorts = format_cohort_list(cohort[significant_PC_axis]),
    in_score_cohorts = format_cohort_list(cohort[include_in_score]),
    .groups = "drop"
  ) %>%
  mutate(
    recurrence_fraction = round(n_cohorts_significant_PC_axis / n_cohorts_tested, 3),
    best_abs_PC_axis_assoc = round(best_abs_PC_axis_assoc, 3)
  ) %>%
  arrange(desc(n_cohorts_significant_PC_axis), desc(n_cohorts_in_score), metric)

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write_tsv(
  all_associations,
  file.path(output_dir, "tcga_qc_pc_axis_all_associations_from_html.tsv")
)
write_tsv(
  significant_associations,
  file.path(output_dir, "tcga_qc_pc_axis_significant_from_html.tsv")
)
write_tsv(
  parameter_summary,
  file.path(output_dir, "tcga_qc_pc_axis_parameter_summary_from_html.tsv")
)

top_parameters <- parameter_summary %>%
  filter(n_cohorts_significant_PC_axis > 0) %>%
  head(30)

markdown_lines <- c(
  "# TCGA QC Parameters Significantly Associated With PCA Axes",
  "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Scope",
  "",
  "This summary was extracted from the local ShinyQC HTML reports because the `/Volumes/SEQC2` TCGA raw-data mount was not available during this run.",
  "",
  "It summarizes QC metrics that significantly track PC1 or PC2. It does not include the separate PCA-distance outlier-score test; that requires recomputing from the raw normalized counts and QC files.",
  "",
  "## Significance Rule",
  "",
  "- Numeric metrics: absolute Spearman correlation with PC1 or PC2 at least 0.4 and FDR < 0.05.",
  "- Categorical metrics: Kruskal-Wallis effect size with PC1 or PC2 at least 0.1 and FDR < 0.05.",
  "",
  "## Most Recurrent Significant QC Parameters",
  "",
  markdown_table(
    top_parameters,
    c(
      "metric",
      "type",
      "n_cohorts_tested",
      "n_cohorts_significant_PC_axis",
      "n_cohorts_in_score",
      "best_abs_PC_axis_assoc",
      "best_PC_axis_cohort",
      "significant_cohorts"
    )
  ),
  "",
  "## Output Tables",
  "",
  "- `tcga_qc_pc_axis_all_associations_from_html.tsv`: every metric from every cohort report.",
  "- `tcga_qc_pc_axis_significant_from_html.tsv`: significant cohort-metric associations only.",
  "- `tcga_qc_pc_axis_parameter_summary_from_html.tsv`: recurrence summary by QC metric."
)

writeLines(markdown_lines, file.path(output_dir, "tcga_qc_pc_axis_parameter_summary_from_html.md"))

message("Read ", length(report_files), " reports")
message("Wrote outputs to: ", output_dir)
