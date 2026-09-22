#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
})

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/summarize_tcga_qc_pca_associations.R manifest.tsv output_dir\n\n",
    "Outputs:\n",
    "  tcga_qc_pca_all_associations.tsv\n",
    "  tcga_qc_pca_significant_associations.tsv\n",
    "  tcga_qc_pca_parameter_summary.tsv\n",
    "  tcga_qc_pca_parameter_summary.md\n",
    sep = ""
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2 || any(args %in% c("-h", "--help"))) {
  usage()
  quit(status = ifelse(length(args) == 2, 0, 1))
}

manifest_file <- args[[1]]
output_dir <- args[[2]]

read_data_file <- function(file_path) {
  if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
    read_delim(file_path, delim = ",", show_col_types = FALSE)
  } else if (grepl("\\.(tsv|txt)$", file_path, ignore.case = TRUE)) {
    read_delim(file_path, delim = "\t", show_col_types = FALSE)
  } else {
    stop("Unsupported file type: ", file_path)
  }
}

robust_z_score <- function(x) {
  x <- as.numeric(x)
  center <- median(x, na.rm = TRUE)
  scale <- mad(x, constant = 1.4826, na.rm = TRUE)

  if (!is.finite(scale) || scale == 0) {
    scale <- sd(x, na.rm = TRUE)
  }

  if (!is.finite(scale) || scale == 0) {
    return(rep(0, length(x)))
  }

  (x - center) / scale
}

kruskal_effect <- function(values, groups) {
  keep <- !is.na(values) & !is.na(groups)
  values <- values[keep]
  groups <- as.factor(groups[keep])
  n <- length(values)
  k <- length(unique(groups))

  if (n < 3 || k < 2 || k >= n) {
    return(c(effect = NA_real_, p_value = NA_real_))
  }

  test <- kruskal.test(values ~ groups)
  effect <- (as.numeric(test$statistic) - k + 1) / (n - k)
  effect <- max(0, min(1, effect))

  c(effect = effect, p_value = test$p.value)
}

spearman_result <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]
  y <- y[keep]

  if (length(x) < 3 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(c(assoc = NA_real_, p_value = NA_real_))
  }

  result <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  c(assoc = as.numeric(result$estimate), p_value = result$p.value)
}

selected_manifest_vars <- function(row, qc_vars) {
  selected_vars <- character(0)
  if ("selected_vars" %in% names(row) && !is.na(row$selected_vars) && row$selected_vars != "") {
    selected_vars <- trimws(strsplit(row$selected_vars, ",")[[1]])
  }
  if (length(selected_vars) == 0 || identical(toupper(selected_vars), "ALL_QC")) {
    selected_vars <- qc_vars
  }
  selected_vars
}

cohort_from_report <- function(row) {
  if ("report_name" %in% names(row) && !is.na(row$report_name) && row$report_name != "") {
    return(gsub("_QC_report$", "", row$report_name))
  }
  tools::file_path_sans_ext(basename(row$counts_file))
}

analyze_cohort <- function(row) {
  cohort <- cohort_from_report(row)
  message("Analyzing ", cohort)

  counts_file <- row$counts_file
  sample_meta_file <- row$sample_meta_file
  qc_file <- row$qc_file
  sample_meta_col <- row$sample_meta_col
  qc_col <- row$qc_col

  nc <- read_data_file(counts_file)
  numeric_columns <- sapply(nc, is.numeric) | names(nc) == "Gene"
  nc <- nc[, numeric_columns]
  names(nc) <- ifelse(names(nc) == "Gene", "Gene", gsub("_", "-", names(nc)))

  if (!"Gene" %in% names(nc)) {
    stop("The 'Gene' column is missing from: ", counts_file)
  }

  nc <- column_to_rownames(nc, "Gene")
  metadata <- read_data_file(sample_meta_file)
  qc <- read_data_file(qc_file)
  qc_vars <- setdiff(names(qc), qc_col)
  selected_vars <- selected_manifest_vars(row, qc_vars)

  samples <- colnames(nc)
  edf_orig <- as.data.frame(lapply(nc, as.numeric), check.names = FALSE)
  edf_filt <- edf_orig[rowMeans(edf_orig) != 0, ] %>% select(all_of(samples))

  met_filt <- metadata %>% filter(.data[[sample_meta_col]] %in% samples)
  sample_df <- merge(qc, met_filt, by.x = qc_col, by.y = sample_meta_col)
  sample_df <- sample_df %>%
    filter(.data[[qc_col]] %in% samples) %>%
    distinct(.data[[qc_col]], .keep_all = TRUE) %>%
    arrange(match(.data[[qc_col]], samples))

  if (nrow(sample_df) == 0) {
    stop("No samples matched after merging metadata and QC files for: ", counts_file)
  }

  selected_vars <- selected_vars[selected_vars %in% names(sample_df)]
  edf_filt <- edf_filt %>% select(all_of(sample_df[[qc_col]]))
  tedf <- t(edf_filt)
  tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf), drop = FALSE]
  tedf <- tedf[, apply(tedf, 2, var) != 0, drop = FALSE]

  pca <- prcomp(tedf, scale. = TRUE)
  pca_df <- dplyr::select(as.data.frame(pca$x), PC1, PC2)
  pca_outlier_score <- sqrt(robust_z_score(pca_df$PC1)^2 + robust_z_score(pca_df$PC2)^2)
  percent_var <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

  rows <- lapply(selected_vars, function(metric) {
    values <- sample_df[[metric]]
    metric_type <- if (is.numeric(values)) "numeric" else "categorical"

    if (metric_type == "numeric") {
      pc1 <- spearman_result(values, pca_df$PC1)
      pc2 <- spearman_result(values, pca_df$PC2)
      outlier <- spearman_result(values, pca_outlier_score)
    } else {
      pc1 <- kruskal_effect(pca_df$PC1, values)
      pc2 <- kruskal_effect(pca_df$PC2, values)
      outlier <- kruskal_effect(pca_outlier_score, values)
      names(pc1) <- c("assoc", "p_value")
      names(pc2) <- c("assoc", "p_value")
      names(outlier) <- c("assoc", "p_value")
    }

    data.frame(
      cohort = cohort,
      metric = metric,
      type = metric_type,
      n = sum(!is.na(values)),
      n_levels = if (metric_type == "categorical") length(unique(values[!is.na(values)])) else NA_integer_,
      PC1_assoc = pc1["assoc"],
      PC1_p = pc1["p_value"],
      PC2_assoc = pc2["assoc"],
      PC2_p = pc2["p_value"],
      PCA_outlier_assoc = outlier["assoc"],
      PCA_outlier_p = outlier["p_value"],
      PC1_percent_var = percent_var[1],
      PC2_percent_var = percent_var[2],
      stringsAsFactors = FALSE
    )
  })

  bind_rows(rows)
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
    return("")
  }
  cohorts[which.max(ifelse(is.na(values), -Inf, values))][1]
}

markdown_table <- function(data, columns) {
  data <- data[, columns, drop = FALSE]
  header <- paste(names(data), collapse = " | ")
  separator <- paste(rep("---", length(columns)), collapse = " | ")
  body <- apply(data, 1, function(row) paste(row, collapse = " | "))
  paste(c(header, separator, body), collapse = "\n")
}

if (!file.exists(manifest_file)) {
  stop("Manifest file does not exist: ", manifest_file)
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
manifest <- read_data_file(manifest_file)

required_cols <- c("counts_file", "sample_meta_file", "qc_file", "sample_meta_col", "qc_col")
missing_cols <- setdiff(required_cols, names(manifest))
if (length(missing_cols) > 0) {
  stop("Manifest is missing required columns: ", paste(missing_cols, collapse = ", "))
}

all_associations <- bind_rows(lapply(seq_len(nrow(manifest)), function(row_index) {
  analyze_cohort(manifest[row_index, , drop = FALSE])
}))

all_associations <- all_associations %>%
  group_by(cohort) %>%
  mutate(
    PC1_FDR = p.adjust(PC1_p, method = "BH"),
    PC2_FDR = p.adjust(PC2_p, method = "BH"),
    PCA_outlier_FDR = p.adjust(PCA_outlier_p, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    strongest_PC = ifelse(abs(PC1_assoc) >= abs(PC2_assoc), "PC1", "PC2"),
    strongest_PC_assoc = ifelse(
      is.na(PC1_assoc) & is.na(PC2_assoc),
      NA_real_,
      pmax(abs(PC1_assoc), abs(PC2_assoc), na.rm = TRUE)
    ),
    strongest_PC_FDR = ifelse(
      is.na(PC1_FDR) & is.na(PC2_FDR),
      NA_real_,
      pmin(PC1_FDR, PC2_FDR, na.rm = TRUE)
    ),
    magnitude_threshold = ifelse(type == "numeric", 0.4, 0.1),
    significant_PC_axis = !is.na(strongest_PC_assoc) &
      !is.na(strongest_PC_FDR) &
      strongest_PC_assoc >= magnitude_threshold &
      strongest_PC_FDR < 0.05,
    significant_PCA_outlier = !is.na(PCA_outlier_assoc) &
      !is.na(PCA_outlier_FDR) &
      abs(PCA_outlier_assoc) >= magnitude_threshold &
      PCA_outlier_FDR < 0.05,
    significant_any = significant_PC_axis | significant_PCA_outlier
  ) %>%
  mutate(
    strongest_PC_assoc = ifelse(is.finite(strongest_PC_assoc), strongest_PC_assoc, NA_real_),
    strongest_PC_FDR = ifelse(is.finite(strongest_PC_FDR), strongest_PC_FDR, NA_real_),
    significant_PC_axis = ifelse(is.na(significant_PC_axis), FALSE, significant_PC_axis),
    significant_PCA_outlier = ifelse(is.na(significant_PCA_outlier), FALSE, significant_PCA_outlier),
    significant_any = ifelse(is.na(significant_any), FALSE, significant_any)
  ) %>%
  mutate(across(
    c(
      PC1_assoc,
      PC1_p,
      PC2_assoc,
      PC2_p,
      PCA_outlier_assoc,
      PCA_outlier_p,
      PC1_percent_var,
      PC2_percent_var,
      PC1_FDR,
      PC2_FDR,
      PCA_outlier_FDR,
      strongest_PC_assoc,
      strongest_PC_FDR
    ),
    ~ signif(.x, 5)
  ))

significant_associations <- all_associations %>%
  filter(significant_any) %>%
  arrange(metric, cohort)

parameter_summary <- all_associations %>%
  group_by(metric, type) %>%
  summarise(
    n_cohorts_tested = n_distinct(cohort),
    n_cohorts_significant_any = n_distinct(cohort[significant_any]),
    n_cohorts_significant_PC_axis = n_distinct(cohort[significant_PC_axis]),
    n_cohorts_significant_PCA_outlier = n_distinct(cohort[significant_PCA_outlier]),
    best_abs_PC_axis_assoc = safe_max(strongest_PC_assoc),
    best_PC_axis_cohort = safe_best_cohort(cohort, strongest_PC_assoc),
    best_abs_PCA_outlier_assoc = safe_max(abs(PCA_outlier_assoc)),
    best_PCA_outlier_cohort = safe_best_cohort(cohort, abs(PCA_outlier_assoc)),
    significant_any_cohorts = format_cohort_list(cohort[significant_any]),
    significant_outlier_cohorts = format_cohort_list(cohort[significant_PCA_outlier]),
    .groups = "drop"
  ) %>%
  mutate(
    recurrence_fraction = n_cohorts_significant_any / n_cohorts_tested,
    best_abs_PC_axis_assoc = round(best_abs_PC_axis_assoc, 3),
    best_abs_PCA_outlier_assoc = round(best_abs_PCA_outlier_assoc, 3),
    recurrence_fraction = round(recurrence_fraction, 3)
  ) %>%
  arrange(desc(n_cohorts_significant_any), desc(n_cohorts_significant_PCA_outlier), metric)

write_tsv(
  all_associations,
  file.path(output_dir, "tcga_qc_pca_all_associations.tsv")
)
write_tsv(
  significant_associations,
  file.path(output_dir, "tcga_qc_pca_significant_associations.tsv")
)
write_tsv(
  parameter_summary,
  file.path(output_dir, "tcga_qc_pca_parameter_summary.tsv")
)

top_any <- parameter_summary %>%
  filter(n_cohorts_significant_any > 0) %>%
  head(20)
top_outlier <- parameter_summary %>%
  filter(n_cohorts_significant_PCA_outlier > 0) %>%
  arrange(desc(n_cohorts_significant_PCA_outlier), desc(n_cohorts_significant_any), metric) %>%
  head(20)

markdown_lines <- c(
  "# TCGA QC Parameters Associated With PCA Structure",
  "",
  paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Definition",
  "",
  "For each TCGA cohort, PCA was recomputed from the normalized count matrix using the same sample matching rules as the ShinyQC reports.",
  "",
  "Numeric QC metrics were tested with Spearman correlation against PC1, PC2, and a PCA outlier score defined as robust distance from the PC1/PC2 center.",
  "",
  "Categorical QC metrics were tested with Kruskal-Wallis effect size against PC1, PC2, and the PCA outlier score.",
  "",
  "A metric was called significant for a cohort if FDR < 0.05 and the association magnitude was at least 0.4 for numeric metrics or 0.1 for categorical metrics.",
  "",
  "## Most Recurrent QC Parameters Associated With PC1/PC2 Or PCA Outlier Score",
  "",
  markdown_table(
    top_any,
    c(
      "metric",
      "type",
      "n_cohorts_tested",
      "n_cohorts_significant_any",
      "n_cohorts_significant_PC_axis",
      "n_cohorts_significant_PCA_outlier",
      "best_abs_PC_axis_assoc",
      "best_abs_PCA_outlier_assoc"
    )
  ),
  "",
  "## Most Recurrent QC Parameters Associated With PCA Outlier Score",
  "",
  markdown_table(
    top_outlier,
    c(
      "metric",
      "type",
      "n_cohorts_tested",
      "n_cohorts_significant_PCA_outlier",
      "best_abs_PCA_outlier_assoc",
      "best_PCA_outlier_cohort",
      "significant_outlier_cohorts"
    )
  ),
  "",
  "## Output Tables",
  "",
  "- `tcga_qc_pca_all_associations.tsv`: every cohort-metric test.",
  "- `tcga_qc_pca_significant_associations.tsv`: cohort-metric tests passing the significance rule.",
  "- `tcga_qc_pca_parameter_summary.tsv`: recurrence summary by QC metric."
)

writeLines(markdown_lines, file.path(output_dir, "tcga_qc_pca_parameter_summary.md"))

message("Wrote outputs to: ", output_dir)
