#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(htmltools)
  library(plotly)
  library(readr)
  library(tibble)
})

usage <- function() {
  cat(
    "Usage:\n",
    "  Rscript scripts/batch_qc_reports.R manifest.tsv output_dir\n\n",
    "Manifest columns:\n",
    "  counts_file       Path to normalized counts TSV/CSV/TXT with Gene column\n",
    "  sample_meta_file  Path to sample metadata TSV/CSV/TXT\n",
    "  qc_file           Path to QC metadata TSV/CSV/TXT\n",
    "  sample_meta_col   Column in sample metadata matching count sample names\n",
    "  qc_col            Column in QC metadata matching count sample names\n",
    "Optional columns:\n",
    "  selected_vars     Comma-separated QC columns to plot/score, or ALL_QC\n",
    "  context_vars      Comma-separated sample metadata columns to annotate PCs, or ALL_CONTEXT\n",
    "  report_name       Output HTML filename without extension\n",
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

safe_name <- function(value) {
  value <- gsub("[^A-Za-z0-9._-]+", "_", value)
  value <- gsub("_+", "_", value)
  trimws(value, whitespace = "_")
}

html_table <- function(data) {
  if (nrow(data) == 0) {
    return(tags$p(class = "empty-table", "No variables available for this section."))
  }

  tags$table(
    class = "score-table",
    tags$thead(tags$tr(lapply(names(data), tags$th))),
    tags$tbody(
      lapply(seq_len(nrow(data)), function(row_index) {
        tags$tr(
          lapply(data[row_index, , drop = FALSE], function(value) {
            tags$td(as.character(value))
          })
        )
      })
    )
  )
}

parse_manifest_vars <- function(row, column, default_vars, all_vars, all_token) {
  vars <- character(0)
  if (column %in% names(row) && !is.na(row[[column]]) && row[[column]] != "") {
    vars <- trimws(strsplit(as.character(row[[column]]), ",")[[1]])
  }
  if (length(vars) == 0) {
    vars <- default_vars
  }
  if (length(vars) == 1 && identical(toupper(vars), all_token)) {
    vars <- all_vars
  }
  unique(vars[nzchar(vars)])
}

missing_metadata_tokens <- c(
  "",
  "NA",
  "N/A",
  "_Not_Available_",
  "_Not_Applicable_",
  "_Unknown_",
  "_Not_Evaluated_",
  "[Not Available]",
  "[Not Applicable]",
  "[Unknown]",
  "[Not Evaluated]"
)

clean_metadata_values <- function(values) {
  if (!is.character(values)) {
    return(values)
  }

  values <- trimws(values)
  values[values %in% missing_metadata_tokens] <- NA_character_
  type.convert(values, as.is = TRUE)
}

clean_context_columns <- function(sample_df, context_vars) {
  for (var in context_vars) {
    sample_df[[var]] <- clean_metadata_values(sample_df[[var]])
  }
  sample_df
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

spearman_effect <- function(values, target) {
  keep <- !is.na(values) & !is.na(target)
  values <- values[keep]
  target <- target[keep]

  if (length(values) < 3 || length(unique(values)) < 2 || length(unique(target)) < 2) {
    return(c(effect = NA_real_, p_value = NA_real_))
  }

  test <- suppressWarnings(cor.test(values, target, method = "spearman", exact = FALSE))
  c(effect = as.numeric(test$estimate), p_value = test$p.value)
}

pc_association_table <- function(sample_df, selected_vars, pca_df, score_numeric = FALSE) {
  if (length(selected_vars) == 0) {
    return(data.frame(
      metric = character(0),
      type = character(0),
      n = integer(0),
      PC1_assoc = numeric(0),
      PC1_FDR = numeric(0),
      PC2_assoc = numeric(0),
      PC2_FDR = numeric(0),
      strongest_PC = character(0),
      strongest_abs_assoc = numeric(0),
      include_in_score = logical(0)
    ))
  }

  rows <- lapply(selected_vars, function(var) {
    values <- sample_df[[var]]

    if (is.numeric(values)) {
      pc1 <- spearman_effect(values, pca_df$PC1)
      pc2 <- spearman_effect(values, pca_df$PC2)

      data.frame(
        metric = var,
        type = "numeric",
        n = sum(!is.na(values)),
        PC1_assoc = pc1["effect"],
        PC1_p = pc1["p_value"],
        PC2_assoc = pc2["effect"],
        PC2_p = pc2["p_value"]
      )
    } else {
      pc1 <- kruskal_effect(pca_df$PC1, values)
      pc2 <- kruskal_effect(pca_df$PC2, values)

      data.frame(
        metric = var,
        type = "categorical",
        n = sum(!is.na(values)),
        PC1_assoc = pc1["effect"],
        PC1_p = pc1["p_value"],
        PC2_assoc = pc2["effect"],
        PC2_p = pc2["p_value"]
      )
    }
  })

  association_table <- bind_rows(rows)
  association_table$PC1_FDR <- p.adjust(association_table$PC1_p, method = "BH")
  association_table$PC2_FDR <- p.adjust(association_table$PC2_p, method = "BH")
  association_table$strongest_PC <- ifelse(
    abs(association_table$PC1_assoc) >= abs(association_table$PC2_assoc),
    "PC1",
    "PC2"
  )
  association_table$strongest_abs_assoc <- ifelse(
    is.na(association_table$PC1_assoc) & is.na(association_table$PC2_assoc),
    NA_real_,
    pmax(abs(association_table$PC1_assoc), abs(association_table$PC2_assoc), na.rm = TRUE)
  )
  min_fdr <- ifelse(
    is.na(association_table$PC1_FDR) & is.na(association_table$PC2_FDR),
    NA_real_,
    pmin(association_table$PC1_FDR, association_table$PC2_FDR, na.rm = TRUE)
  )
  association_table$include_in_score <- (
    score_numeric &
      association_table$type == "numeric" &
      association_table$strongest_abs_assoc >= 0.4 &
      min_fdr < 0.05
  )
  association_table$include_in_score[is.na(association_table$include_in_score)] <- FALSE

  association_table %>%
    arrange(desc(include_in_score), desc(strongest_abs_assoc)) %>%
    mutate(
      PC1_assoc = round(PC1_assoc, 3),
      PC1_FDR = signif(PC1_FDR, 3),
      PC2_assoc = round(PC2_assoc, 3),
      PC2_FDR = signif(PC2_FDR, 3),
      strongest_abs_assoc = round(strongest_abs_assoc, 3)
    ) %>%
    select(
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

qc_pc_association_table <- function(sample_df, selected_vars, pca_df) {
  pc_association_table(sample_df, selected_vars, pca_df, score_numeric = TRUE)
}

qc_redundancy_table <- function(sample_df, selected_vars, association_table, threshold = 0.7) {
  numeric_vars <- selected_vars[sapply(sample_df[, selected_vars, drop = FALSE], is.numeric)]
  usable_vars <- numeric_vars[
    sapply(numeric_vars, function(var) {
      values <- sample_df[[var]]
      sum(!is.na(values)) >= 3 && is.finite(sd(values, na.rm = TRUE)) && sd(values, na.rm = TRUE) > 0
    })
  ]

  if (length(usable_vars) == 0) {
    return(data.frame(
      group = integer(0),
      representative_metric = character(0),
      n_metrics = integer(0),
      metrics = character(0),
      max_pairwise_abs_corr = numeric(0),
      strongest_PC = character(0),
      strongest_abs_assoc = numeric(0),
      include_in_score = logical(0)
    ))
  }

  cor_mat <- suppressWarnings(cor(
    sample_df[, usable_vars, drop = FALSE],
    method = "spearman",
    use = "pairwise.complete.obs"
  ))
  cor_mat[!is.finite(cor_mat)] <- 0
  diag(cor_mat) <- 1

  if (length(usable_vars) == 1) {
    groups <- setNames(1, usable_vars)
  } else {
    groups <- cutree(
      hclust(as.dist(1 - abs(cor_mat)), method = "average"),
      h = 1 - threshold
    )
  }

  numeric_assoc <- association_table %>%
    filter(metric %in% usable_vars) %>%
    mutate(
      assoc_rank = ifelse(is.na(strongest_abs_assoc), -Inf, strongest_abs_assoc),
      include_rank = ifelse(include_in_score, 1, 0)
    )

  bind_rows(lapply(sort(unique(groups)), function(group_id) {
    members <- names(groups)[groups == group_id]
    member_assoc <- numeric_assoc %>%
      filter(metric %in% members) %>%
      arrange(desc(include_rank), desc(assoc_rank), metric)
    representative <- member_assoc$metric[1]
    member_cor <- cor_mat[members, members, drop = FALSE]
    max_pairwise <- if (length(members) > 1) {
      max(abs(member_cor[upper.tri(member_cor)]), na.rm = TRUE)
    } else {
      NA_real_
    }

    data.frame(
      group = group_id,
      representative_metric = representative,
      n_metrics = length(members),
      metrics = paste(sort(members), collapse = ", "),
      max_pairwise_abs_corr = round(max_pairwise, 3),
      strongest_PC = member_assoc$strongest_PC[1],
      strongest_abs_assoc = member_assoc$strongest_abs_assoc[1],
      include_in_score = any(member_assoc$include_in_score),
      stringsAsFactors = FALSE
    )
  })) %>%
    arrange(desc(include_in_score), desc(strongest_abs_assoc), group)
}

make_plot <- function(pca_df, sample_df, sample_col, variable) {
  plot_df <- pca_df
  plot_df$variable <- sample_df[[variable]]
  plot_df$sample <- sample_df[[sample_col]]
  plot_df$name <- variable
  numeric_values <- if (is.numeric(plot_df$variable)) {
    plot_df$variable[is.finite(plot_df$variable)]
  } else {
    numeric(0)
  }

  if (is.numeric(plot_df$variable) && length(numeric_values) > 0 && length(unique(numeric_values)) > 1) {
    p <- ggplot(
      plot_df,
      aes(
        x = PC1,
        y = PC2,
        text = paste("Sample:", sample, "<br>Value:", variable)
      )
    ) +
      theme_bw() +
      geom_point(aes(color = variable), size = 1) +
      scale_color_gradient2(
        low = "#2b83ba",
        mid = "grey",
        high = "#d7191c",
        midpoint = median(plot_df$variable, na.rm = TRUE),
        limits = range(numeric_values, na.rm = TRUE),
        oob = scales::squish
      ) +
      facet_wrap(~name)
  } else {
    plot_df$variable <- as.factor(plot_df$variable)
    p <- ggplot(
      plot_df,
      aes(
        x = PC1,
        y = PC2,
        text = paste("Sample:", sample, "<br>Value:", variable)
      )
    ) +
      theme_bw() +
      geom_point(aes(color = variable), size = 1) +
      facet_wrap(~name)
  }

  ggplotly(p, tooltip = "text")
}

make_report <- function(row, output_dir) {
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
  metadata_vars <- setdiff(names(metadata), sample_meta_col)

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

  edf_filt <- edf_filt %>% select(all_of(sample_df[[qc_col]]))
  tedf <- t(edf_filt)
  tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
  tedf <- tedf[, apply(tedf, 2, var) != 0]

  pca <- prcomp(tedf, scale. = TRUE)
  pca_df <- dplyr::select(as.data.frame(pca$x), PC1, PC2)
  pca_df$sample <- rownames(pca_df)

  selected_vars <- parse_manifest_vars(row, "selected_vars", qc_vars, qc_vars, "ALL_QC")
  selected_vars <- selected_vars[selected_vars %in% names(sample_df)]
  context_vars <- parse_manifest_vars(
    row,
    "context_vars",
    default_context_vars,
    metadata_vars,
    "ALL_CONTEXT"
  )
  context_vars <- context_vars[context_vars %in% names(sample_df)]
  context_vars <- setdiff(context_vars, selected_vars)
  sample_df <- clean_context_columns(sample_df, context_vars)

  association_table <- qc_pc_association_table(sample_df, selected_vars, pca_df)
  context_association_table <- pc_association_table(sample_df, context_vars, pca_df) %>%
    select(-include_in_score)
  redundancy_table <- qc_redundancy_table(sample_df, selected_vars, association_table)
  score_vars <- redundancy_table %>%
    filter(include_in_score) %>%
    pull(representative_metric)
  if (length(score_vars) == 0) {
    score_vars <- selected_vars[sapply(sample_df[, selected_vars, drop = FALSE], is.numeric)]
  }
  pca_score <- sqrt(robust_z_score(pca_df$PC1)^2 + robust_z_score(pca_df$PC2)^2)

  if (length(score_vars) == 0) {
    score_table <- data.frame(
      sample = sample_df[[qc_col]],
      PC1 = round(pca_df$PC1, 3),
      PC2 = round(pca_df$PC2, 3),
      pca_outlier_score = round(pca_score, 3)
    ) %>%
      arrange(desc(pca_outlier_score)) %>%
      head(25)
  } else {
    qc_z <- sapply(score_vars, function(var) robust_z_score(sample_df[[var]]))
    if (length(score_vars) == 1) {
      qc_z <- matrix(qc_z, ncol = 1)
      colnames(qc_z) <- score_vars
    }
    score_weights <- association_table$strongest_abs_assoc[
      match(score_vars, association_table$metric)
    ]
    score_weights[!is.finite(score_weights) | is.na(score_weights)] <- 1
    weighted_qc_z <- sweep(abs(qc_z), 2, score_weights, `*`)
    strongest_index <- max.col(abs(qc_z), ties.method = "first")
    strongest_parameter <- colnames(qc_z)[strongest_index]
    strongest_z <- qc_z[cbind(seq_len(nrow(qc_z)), strongest_index)]

    score_table <- data.frame(
      sample = sample_df[[qc_col]],
      qc_outlier_score = round(rowSums(weighted_qc_z, na.rm = TRUE), 3),
      pca_outlier_score = round(pca_score, 3),
      PC1 = round(pca_df$PC1, 3),
      PC2 = round(pca_df$PC2, 3),
      strongest_parameter = strongest_parameter,
      strongest_z = round(strongest_z, 3)
    ) %>%
      arrange(desc(qc_outlier_score), desc(pca_outlier_score)) %>%
      head(25)
  }

  report_name <- if ("report_name" %in% names(row) && !is.na(row$report_name) && row$report_name != "") {
    safe_name(row$report_name)
  } else {
    safe_name(tools::file_path_sans_ext(basename(counts_file)))
  }
  report_file <- file.path(output_dir, paste0(report_name, ".html"))
  libdir <- paste0(report_name, "_files")

  report_plots <- lapply(selected_vars, function(var) {
    tags$div(
      class = "plot-panel",
      tags$h3(var),
      make_plot(pca_df, sample_df, qc_col, var)
    )
  })
  context_plots <- lapply(context_vars, function(var) {
    tags$div(
      class = "plot-panel",
      tags$h3(var),
      make_plot(pca_df, sample_df, qc_col, var)
    )
  })

  report <- tagList(
    tags$head(
      tags$title("ShinyQC Inspector report"),
      tags$style(HTML("
        body { font-family: Arial, sans-serif; margin: 24px; color: #222; }
        h1 { margin-bottom: 0; }
        .subtitle { color: #555; margin-top: 4px; }
        .plot-grid {
          display: grid;
          grid-template-columns: repeat(5, minmax(0, 1fr));
          gap: 12px;
          align-items: start;
        }
        .plot-panel {
          min-width: 0;
          break-inside: avoid;
        }
        .plot-panel h3 {
          font-size: 16px;
          margin: 0 0 6px 0;
          overflow-wrap: anywhere;
        }
        .plot-panel .html-widget {
          width: 100% !important;
          height: 260px !important;
        }
        @media (max-width: 1600px) {
          .plot-grid { grid-template-columns: repeat(4, minmax(0, 1fr)); }
        }
        @media (max-width: 1400px) {
          .plot-grid { grid-template-columns: repeat(3, minmax(0, 1fr)); }
        }
        @media (max-width: 1100px) {
          .plot-grid { grid-template-columns: repeat(2, minmax(0, 1fr)); }
        }
        @media (max-width: 720px) {
          .plot-grid { grid-template-columns: 1fr; }
        }
        .score-table { border-collapse: collapse; margin: 18px 0 28px 0; }
        .score-table th, .score-table td {
          border: 1px solid #ccc;
          padding: 6px 8px;
          font-size: 13px;
          text-align: left;
        }
        .score-table th { background: #f1f1f1; }
      "))
    ),
    tags$h1("ShinyQC Inspector report"),
    tags$p(class = "subtitle", paste("Saved", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))),
    tags$p(
      paste("Counts:", counts_file),
      tags$br(),
      paste("Sample metadata:", sample_meta_file),
      tags$br(),
      paste("QC metadata:", qc_file),
      tags$br(),
      paste("Sample metadata column:", sample_meta_col),
      tags$br(),
      paste("QC metadata column:", qc_col),
      tags$br(),
      paste("Selected QC parameters:", paste(selected_vars, collapse = ", ")),
      tags$br(),
      paste("Clinical/pathology context variables:", paste(context_vars, collapse = ", "))
    ),
    tags$h2("Top QC Outlier Scores"),
    html_table(score_table),
    tags$h2("QC Associations With PC1/PC2"),
    tags$p(
      "Numeric metrics use Spearman correlation. Categorical metrics use a Kruskal-Wallis effect size. ",
      "Numeric metrics are marked include_in_score when abs(association) >= 0.4 and FDR < 0.05 for PC1 or PC2."
    ),
    html_table(association_table),
    tags$h2("Clinical/Pathology Context With PC1/PC2"),
    tags$p(
      "These sample metadata variables are tested against PC1 and PC2 to help distinguish biological or clinical structure from technical QC bias. ",
      "They are not used in the QC outlier score."
    ),
    html_table(context_association_table),
    tags$h2("QC Redundancy Groups"),
    tags$p(
      "Numeric QC metrics are clustered by absolute Spearman correlation. ",
      "The representative metric is the strongest PC-associated member of each group."
    ),
    html_table(redundancy_table),
    tags$h2("PCA QC Plots"),
    tags$div(class = "plot-grid", report_plots),
    tags$h2("PCA Clinical/Pathology Context Plots"),
    tags$div(class = "plot-grid", context_plots)
  )

  save_html(report, report_file, libdir = libdir)
  report_file
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

for (row_index in seq_len(nrow(manifest))) {
  cat("Rendering report", row_index, "of", nrow(manifest), "...\n")
  report_path <- make_report(manifest[row_index, , drop = FALSE], output_dir)
  cat("Saved:", report_path, "\n")
}
