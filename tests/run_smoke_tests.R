#!/usr/bin/env Rscript

# Minimal, data-free regression test for the app and batch report pipeline.

file_argument <- commandArgs(trailingOnly = FALSE)
script_path <- sub("^--file=", "", file_argument[grep("^--file=", file_argument)][[1]])
repo_root <- normalizePath(file.path(dirname(script_path), ".."))
setwd(repo_root)

run_rscript <- function(script, args) {
  result <- system2(
    file.path(R.home("bin"), "Rscript"),
    c("--vanilla", script, args),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(result, "status")
  if (!is.null(status) && status != 0) {
    stop(paste(c("Command failed:", result), collapse = "\n"))
  }
  invisible(result)
}

test_dir <- tempfile("shinyqc-smoke-")
dir.create(test_dir, recursive = TRUE)
on.exit(unlink(test_dir, recursive = TRUE, force = TRUE), add = TRUE)

set.seed(42)
samples <- paste0("S", sprintf("%02d", seq_len(12)))
counts <- data.frame(
  Gene = paste0("G", seq_len(30)),
  matrix(
    round(rexp(30 * length(samples), rate = 1 / 100) + 10),
    nrow = 30,
    dimnames = list(NULL, samples)
  ),
  check.names = FALSE
)
sample_metadata <- data.frame(
  sample = samples,
  stage = rep(c("I", "II", "III"), each = 4),
  age = seq(45, 78, length.out = length(samples))
)
qc_metadata <- data.frame(
  Sample = samples,
  RIN = round(c(seq(9.5, 7, length.out = 10), 4, 3), 2),
  coverage = round(c(seq(95, 75, length.out = 10), 45, 35), 2),
  constant_metric = 1,
  Encoding = rep(c("A", "B"), length.out = length(samples))
)

write.table(counts, file.path(test_dir, "counts.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(sample_metadata, file.path(test_dir, "sample.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(qc_metadata, file.path(test_dir, "qc.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)

manifest <- data.frame(
  counts_file = file.path(test_dir, "counts.tsv"),
  sample_meta_file = file.path(test_dir, "sample.tsv"),
  qc_file = file.path(test_dir, "qc.tsv"),
  sample_meta_col = "sample",
  qc_col = "Sample",
  selected_vars = "RIN,coverage,constant_metric,Encoding",
  context_vars = "stage,age",
  report_name = "IODC-0001-TCGA-TEST_QC_report"
)
manifest_file <- file.path(test_dir, "manifest.tsv")
report_dir <- file.path(test_dir, "reports")
summary_dir <- file.path(test_dir, "summary")
write.table(manifest, manifest_file, sep = "\t", row.names = FALSE, quote = FALSE)

run_rscript("scripts/batch_qc_reports.R", c(manifest_file, report_dir))
report_file <- file.path(report_dir, "IODC-0001-TCGA-TEST_QC_report.html")
stopifnot(file.exists(report_file), file.info(report_file)$size > 0)

run_rscript("scripts/summarize_tcga_qc_pc_axes_from_html.R", c(report_dir, summary_dir))
summary_file <- file.path(summary_dir, "tcga_qc_pc_axis_parameter_summary_from_html.tsv")
stopifnot(file.exists(summary_file))
summary_table <- read.delim(summary_file, check.names = FALSE)
constant_row <- summary_table[summary_table$metric == "constant_metric", , drop = FALSE]
stopifnot(nrow(constant_row) == 1L, is.na(constant_row$best_abs_PC_axis_assoc))

app_environment <- new.env(parent = globalenv())
app <- source("app.R", local = app_environment)$value
stopifnot(inherits(app, "shiny.appobj"), is.function(app_environment$server))

pca_df <- data.frame(PC1 = seq(-2, 2, length.out = nrow(qc_metadata)), PC2 = seq(2, -2, length.out = nrow(qc_metadata)))
association_table <- app_environment$qc_pc_association_table(
  qc_metadata,
  c("RIN", "constant_metric", "Encoding"),
  pca_df
)
stopifnot(nrow(association_table) == 3L)

cat("ShinyQC smoke tests: OK\n")
