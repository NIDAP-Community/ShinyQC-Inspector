library(shiny)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(RCurl)
library(plotly)
library(dplyr)
library(edgeR)
library(scales)
library(gridExtra)
library(tibble)
library(grid)
library(gridExtra)
library(readr)

# Define user interface

options(shiny.maxRequestSize = 1024 * 1024 * 1024)

ui <- fluidPage(

  titlePanel("ShinyQC Inspector"),
  # Custom CSS to style the subtitle and additional text

  tags$style(HTML("
    .subtitle {
      font-size: 16px; /* smaller than the title */
      color: #333; /* dark grey color */
      margin-bottom: 20px; /* space below the subtitle */
    }
    .info-text {
      font-size: 14px; /* even smaller text */
      color: #666; /* lighter grey */
      line-height: 1.6; /* increased line height for better readability */
      margin-bottom: 15px; /* space below each paragraph */
    }
    .multicol {
      column-count: 6;
      column-gap: 20px;
    }
  ")),

  tags$h3("This application allows the user to inspect various QC parameters",
          class = "subtitle"),

  tags$h4("Normalized data should contain a \"Gene\" column", class = "info-text"),

  tags$h4("Accepts comma (\".csv\") or tab-delimited (\".tsv\" or \".txt\") files", class = "info-text"),

  #textInput("filename", "Filename", value = "shiny.html"),
  tags$style(type = 'text/css', "
             .multicol {
               column-count: 6;
               column-gap: 20px;
             }
             "),

  fluidRow(
    column(width = 12,
           fileInput("file1", "Choose File for Normalized Gene Expression Data", accept = ".tsv,.csv,.txt"),
           fileInput("file2", "Choose File for Sample Metadata", accept = ".csv,.tsv,.txt"),
           fileInput("file3", "Choose File for QC Metadata", accept = ".csv,.tsv,.txt"),
           actionButton("read_data", "Read Data"),
           div(class = 'multicol',
               checkboxGroupInput("vars", "Select Columns:", choices = NULL)
           ),
           selectInput("col_file2", "Select sample column from Sample Metadata:", choices = NULL),
           selectInput("col_file3", "Select sample column from QC Metadata:", choices = NULL),
           actionButton("submit", "Run Analysis"),
           downloadButton('downloadData', 'Download selected samples'),
           textInput("report_dir", "HTML report output folder", value = file.path(dirname(getwd()), "ShinyQC_Output")),
           actionButton("saveReport", "Save HTML report"),
           textOutput("report_save_status"),
           tags$hr(),
           tags$h4("Top QC Outlier Scores"),
           tableOutput("qc_score_table")
    )
  ),
  uiOutput("dynamic_plots")
)

# Helper function to read data based on file extension
read_data_file <- function(file_path) {
  # Determine the separator based on the file extension
  if (grepl("\\.csv$", file_path)) {
    read_delim(file_path, delim = ",", show_col_types = FALSE)
  } else if (grepl("\\.(tsv|txt)$", file_path)) {  # Modified to accept .txt files as well
    read_delim(file_path, delim = "\t", show_col_types = FALSE)
  } else {
    stop("Unsupported file type")
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

qc_pc_association_table <- function(sample_df, selected_vars, pca_df) {
  rows <- lapply(selected_vars, function(var) {
    values <- sample_df[[var]]

    if (is.numeric(values)) {
      pc1 <- suppressWarnings(cor.test(values, pca_df$PC1, method = "spearman", exact = FALSE))
      pc2 <- suppressWarnings(cor.test(values, pca_df$PC2, method = "spearman", exact = FALSE))

      data.frame(
        metric = var,
        type = "numeric",
        n = sum(!is.na(values)),
        PC1_assoc = as.numeric(pc1$estimate),
        PC1_p = pc1$p.value,
        PC2_assoc = as.numeric(pc2$estimate),
        PC2_p = pc2$p.value
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


# Define server logic

server <- function(input, output, session) {

  # Reactive Variables to store the data
  Normalized_Counts <- reactiveVal()
  Metadata <- reactiveVal()
  Sampinfo <- reactiveVal()
  Column_Names <- reactiveVal()
  Selected_Cols <- reactiveVal()

  observeEvent(input$read_data, {
    req(input$file1, input$file2, input$file3)

    withProgress(message = 'Reading files...', value = 0, {
      progress <- shiny::Progress$new()

      # Read files using the helper function and update progress
      nc <- read_data_file(input$file1$datapath)

      # Remove non-numeric columns except for the 'Gene' column
      numeric_columns <- sapply(nc, is.numeric) | names(nc) == "Gene"
      nc <- nc[, numeric_columns]
      names(nc) <- ifelse(names(nc) == "Gene", "Gene", gsub("_", "-", names(nc)))

      # Ensure 'Gene' column exists and is set as row names
      if ("Gene" %in% names(nc)) {
        nc <- column_to_rownames(nc, "Gene")
        progress$inc(1/3, "Finished reading Normalized Counts")
      } else {
        progress$close()
        stop("The 'Gene' column is missing from the dataset.")
      }

      progress$inc(1/3, "Finished reading Normalized Counts")

      meta_data <- read_data_file(input$file2$datapath)
      updateSelectInput(session, "col_file2", choices = names(meta_data))
      progress$inc(1/3, "Finished reading Meta Data")

      samp_info <- read_data_file(input$file3$datapath)
      updateSelectInput(session, "col_file3", choices = names(samp_info))
      progress$inc(1/3, "Finished reading QC Metadata")

      # Store read data in reactive variables
      # Assuming you have defined these reactive variables somewhere in your app
      Metadata(meta_data)
      Sampinfo(samp_info)
      Normalized_Counts(nc)

      on.exit(progress$close())
    })

    # Extract column names
    column_names <- c(colnames(Metadata()), colnames(Sampinfo()))
    print(column_names)

    # List columns you want to preselect
    selected_cols <- colnames(Sampinfo())

    # Update the UI choices
    # updateCheckboxGroupInput(session,
    #                          "vars",
    #                          choices = column_names)
    updateCheckboxGroupInput(session,
                             "vars",
                             choices = column_names,
                             selected = selected_cols)
    Column_Names(column_names)
    #Selected_Cols(selected_cols)
  })

  observeEvent(input$submit, {
    print("Inside observeEvent")
    edf.orig <- as.data.frame(lapply(Normalized_Counts(), as.numeric), check.names = FALSE)
    idx <- rowMeans(edf.orig) != 0
    edf.filt <- edf.orig[idx,]

    samples <- colnames(Normalized_Counts())
    edf.filt <- edf.filt %>% select(all_of(samples))

    col_file2 <- input$col_file2   # Sample metadata
    col_file3 <- input$col_file3   # QC metadata
    met.filt <-  Metadata() %>% filter(.data[[col_file2]] %in% samples)


	    Sample.df <- merge(Sampinfo(), met.filt,
	                       by.x = col_file3, by.y = col_file2)
	    Sample.df <- Sample.df %>% filter(.data[[col_file3]] %in% samples) %>%
	      distinct(.data[[col_file3]], .keep_all = TRUE) %>%
	      arrange(match(.data[[col_file3]], samples))
    edf.filt <- edf.filt %>% select(all_of(Sample.df[[col_file3]]))
    head(edf.filt)
    tedf <- t(edf.filt)
    tedf <- tedf[, colSums(is.na(tedf)) != nrow(tedf)]
    tedf <- tedf[, apply(tedf, 2, var) != 0]

    print("Before PCA....")

    withProgress(message = 'Running PCA...', value = 0, {
      progress <- shiny::Progress$new()
      progress$set() # Set total progress steps
      pca <- prcomp(tedf, scale.=TRUE)
      progress$inc(0.5, "PCA Completed")
      on.exit(progress$close())
    })

    print("After PCA....")
    pca.df <- dplyr::select(as.data.frame(pca$x), PC1, PC2)
    pca.df$sample <- rownames(pca.df)
    print(head(pca.df))

    qc_score_table <- function() {
      req(input$vars)

      selected_vars <- unique(input$vars)
      selected_vars <- selected_vars[selected_vars %in% names(Sample.df)]
      association_table <- qc_pc_association_table(Sample.df, selected_vars, pca.df)
      redundancy_table <- qc_redundancy_table(Sample.df, selected_vars, association_table)
      score_vars <- redundancy_table %>%
        filter(include_in_score) %>%
        pull(representative_metric)

      if (length(score_vars) == 0) {
        score_vars <- selected_vars
      }
      if (length(score_vars) > 0) {
        score_vars <- score_vars[sapply(Sample.df[, score_vars, drop = FALSE], is.numeric)]
      }

      pca_score <- sqrt(robust_z_score(pca.df$PC1)^2 + robust_z_score(pca.df$PC2)^2)

      if (length(score_vars) == 0) {
        return(data.frame(
          sample = Sample.df[[col_file3]],
          PC1 = round(pca.df$PC1, 3),
          PC2 = round(pca.df$PC2, 3),
          pca_outlier_score = round(pca_score, 3)
        ) %>%
          arrange(desc(pca_outlier_score)) %>%
          head(25))
      }

      qc_z <- sapply(score_vars, function(var) robust_z_score(Sample.df[[var]]))
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

      data.frame(
        sample = Sample.df[[col_file3]],
        qc_outlier_score = round(rowSums(weighted_qc_z, na.rm = TRUE), 3),
        pca_outlier_score = round(pca_score, 3),
        PC1 = round(pca.df$PC1, 3),
        PC2 = round(pca.df$PC2, 3),
        strongest_parameter = strongest_parameter,
        strongest_z = round(strongest_z, 3)
      ) %>%
        arrange(desc(qc_outlier_score), desc(pca_outlier_score)) %>%
        head(25)
    }

    output$qc_score_table <- renderTable({
      qc_score_table()
    })

    # Define a function to plot PCA

    plotPCA <- function(qc){

      req(input$vars) # Ensure variable is available

      pca.df$variable <- Sample.df[[qc]]
      pca.df$sample <- Sample.df[[col_file3]]
      pca.df$name <- qc
      pca.df <- pca.df %>% arrange(variable)

      plotcolors <- c("darkred","cadetblue","coral","deeppink",
                      "darkblue","darkgoldenrod","darkolivegreen3", "dodgerblue",
                      "darkorange", "forestgreen", "firebrick", "orchid",
                      "gold", "mediumturquoise", "saddlebrown", "darkviolet", "lightcoral",
                      "limegreen", "deepskyblue", "tomato", "mediumslateblue", "darkgoldenrod",
                      "mediumseagreen", "lightsalmon", "darkolivegreen", "mediumpurple", "sienna")

      num_groups <- length(unique(pca.df$variable))

      generateColors <- function(n) {
        hues <- seq(0, 1, length.out = n + 1)
        colors <- hsv(h = hues[1:n], s = 0.6, v = 0.9) # Adjust s and v if needed
        return(colors)
      }

      if (num_groups > length(plotcolors)) {
        colnum <- num_groups - length(plotcolors)
        colors <- generateColors(colnum)
        plotcolors <- c(plotcolors, colors)
      }

      perc.var <- (pca$sdev^2/sum(pca$sdev^2))*100
      pc.x.lab <- paste0("PC1 ", round(perc.var[1], 2),"%")
      pc.y.lab <- paste0("PC2 ", round(perc.var[2], 2),"%")


      if(class(pca.df$variable) %in% c("factor","character")){
        p <- ggplot(pca.df, aes(x=PC1, y=PC2,
                                text = paste("Sample:", sample, "<br>Value:", variable))) +
          theme_bw() +
          theme(strip.text = element_text(size = 20, color = "white"),
                strip.background = element_rect(fill = "blue"),
                legend.title=element_blank(),
                legend.position="right",
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank()) +
          xlab(pc.x.lab) + ylab(pc.y.lab) +
          geom_point(aes(color=variable), size=1) +
          scale_colour_manual(values = plotcolors) +
          facet_wrap(~name)

        if(qc == "flowcell"){
          p1 <- ggplotly(p, tooltip = "text", source = "plot1Source") %>%
            layout(dragmode = "lasso")
          p1 <- event_register(p1, 'plotly_selected')
          return(p1)
        } else {
          return(p)
        }
      } else {
        return(ggplotly(ggplot(pca.df, aes(x=PC1, y=PC2,
                                           text = paste("Sample:", sample, "<br>Value:", variable))) + #theme_common +
                          theme_bw() +
                          theme(strip.text = element_text(size = 20,color = "white"),
                                strip.background = element_rect(fill = "blue"),
                                legend.title=element_blank(),
                                legend.position="right",
                                panel.grid.major = element_blank(),
                                panel.grid.minor = element_blank(),
                                panel.background = element_blank()) +
                          xlab(pc.x.lab) + ylab(pc.y.lab) +
                          geom_point(aes(color=variable), size=1) +
                          scale_color_gradient2(low = "#2b83ba",
                                                mid = "grey",
                                                high = "#d7191c",
                                                midpoint = median(pca.df$variable),
                                                limits = c(min(pca.df$variable),
                                                           max(pca.df$variable)),
                                                oob = scales::squish) +
                          facet_wrap(~name), tooltip = "text"))
	      }
	    }

	    output$report_save_status <- renderText("")

	    observeEvent(input$saveReport, {
	      req(input$vars)

	      report_dir <- input$report_dir
	      if (is.null(report_dir) || trimws(report_dir) == "") {
	        report_dir <- file.path(dirname(getwd()), "ShinyQC_Output")
	      }

	      dir.create(report_dir, showWarnings = FALSE, recursive = TRUE)

	      report_name <- paste0(
	        "ShinyQC_report_",
	        format(Sys.time(), "%Y%m%d_%H%M%S"),
	        ".html"
	      )
	      report_path <- file.path(report_dir, report_name)
	      libdir <- paste0(tools::file_path_sans_ext(report_name), "_files")

	      score_table <- qc_score_table()
	      association_table <- qc_pc_association_table(Sample.df, input$vars, pca.df)
	      redundancy_table <- qc_redundancy_table(Sample.df, input$vars, association_table)
	      score_table_html <- htmltools::tags$table(
	        class = "score-table",
	        htmltools::tags$thead(
	          htmltools::tags$tr(lapply(names(score_table), htmltools::tags$th))
	        ),
	        htmltools::tags$tbody(
	          lapply(seq_len(nrow(score_table)), function(row_index) {
	            htmltools::tags$tr(
	              lapply(score_table[row_index, , drop = FALSE], function(value) {
	                htmltools::tags$td(as.character(value))
	              })
	            )
	          })
	        )
	      )
	      association_table_html <- htmltools::tags$table(
	        class = "score-table",
	        htmltools::tags$thead(
	          htmltools::tags$tr(lapply(names(association_table), htmltools::tags$th))
	        ),
	        htmltools::tags$tbody(
	          lapply(seq_len(nrow(association_table)), function(row_index) {
	            htmltools::tags$tr(
	              lapply(association_table[row_index, , drop = FALSE], function(value) {
	                htmltools::tags$td(as.character(value))
	              })
	            )
	          })
	        )
	      )
	      redundancy_table_html <- htmltools::tags$table(
	        class = "score-table",
	        htmltools::tags$thead(
	          htmltools::tags$tr(lapply(names(redundancy_table), htmltools::tags$th))
	        ),
	        htmltools::tags$tbody(
	          lapply(seq_len(nrow(redundancy_table)), function(row_index) {
	            htmltools::tags$tr(
	              lapply(redundancy_table[row_index, , drop = FALSE], function(value) {
	                htmltools::tags$td(as.character(value))
	              })
	            )
	          })
	        )
	      )

	      report_plots <- lapply(input$vars, function(var) {
	        plot_obj <- plotPCA(var)
	        if (inherits(plot_obj, "ggplot")) {
	          plot_obj <- ggplotly(plot_obj, tooltip = "text")
	        }

	        htmltools::tags$div(
	          class = "plot-panel",
	          htmltools::tags$h3(var),
	          plot_obj
	        )
	      })

	      report <- htmltools::tagList(
	        htmltools::tags$head(
	          htmltools::tags$title("ShinyQC Inspector report"),
	          htmltools::tags$style(htmltools::HTML("
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
	        htmltools::tags$h1("ShinyQC Inspector report"),
	        htmltools::tags$p(
	          class = "subtitle",
	          paste("Saved", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
	        ),
	        htmltools::tags$p(
	          paste("Sample metadata column:", col_file2),
	          htmltools::tags$br(),
	          paste("QC metadata column:", col_file3),
	          htmltools::tags$br(),
	          paste("Selected QC parameters:", paste(input$vars, collapse = ", "))
	        ),
	        htmltools::tags$h2("Top QC Outlier Scores"),
	        score_table_html,
	        htmltools::tags$h2("QC Associations With PC1/PC2"),
	        htmltools::tags$p(
	          "Numeric metrics use Spearman correlation. Categorical metrics use a Kruskal-Wallis effect size. ",
	          "Numeric metrics are marked include_in_score when abs(association) >= 0.4 and FDR < 0.05 for PC1 or PC2."
	        ),
	        association_table_html,
	        htmltools::tags$h2("QC Redundancy Groups"),
	        htmltools::tags$p(
	          "Numeric QC metrics are clustered by absolute Spearman correlation. ",
	          "The representative metric is the strongest PC-associated member of each group."
	        ),
	        redundancy_table_html,
	        htmltools::tags$h2("PCA QC Plots"),
	        htmltools::tags$div(class = "plot-grid", report_plots)
	      )

	      htmltools::save_html(report, report_path, libdir = libdir)
	      output$report_save_status <- renderText({
	        paste("Saved report:", report_path)
	      })
	    })

	    # Define a reactive variable to store selected samples
	    selected_data <- reactiveVal(data.frame())

    # Dynamically render plots based on selected variables
    output$dynamic_plots <- renderUI({
      req(input$vars)

      # Create a list of plotly outputs
      plots_output_list <- lapply(input$vars, function(var) {
        plotlyOutput(paste0("plot_", var))
      })

      # Break plots list into chunks and put each chunk in a column
      num_plots <- length(plots_output_list)
      num_cols <- 3
      plots_per_col <- ceiling(num_plots / num_cols)

      columns <- lapply(1:num_cols, function(col_num) {
        start_idx = ((col_num-1) * plots_per_col) + 1
        end_idx = min(col_num * plots_per_col, num_plots)
        column(4, plots_output_list[start_idx:end_idx])
      })

      # Return the organized columns to the UI
      do.call(fluidRow, columns)
    })


    # Create separate renderPlotly functions for each variable
    observe({
      req(input$vars)
      lapply(input$vars, function(var) {
        output_name <- paste0("plot_", var)

        output[[output_name]] <- renderPlotly({
          plotPCA(var)
        })
      })
    })

    # observeEvent(input$saveBtn, {
    #   session$sendCustomMessage(type = 'invokeSaveHTML', message = 'dummy')
    # })

    # Observe selected points in the plot
    observeEvent(event_data("plotly_selected", source = "plot1Source"), {
      selected <- event_data("plotly_selected", source = "plot1Source")
      print("hello")
      print(selected)

      if (!is.null(selected)) {
        selected_x <- round(selected$x, digits = 4)
        selected_y <- round(selected$y, digits = 4)

        pca.df$pc1 <- round(pca.df$PC1, digits = 4)
        pca.df$pc2 <- round(pca.df$PC2, digits = 4)
        # Use the x and y values to match the corresponding sample names from pca.df
        selected_samples <- pca.df$sample[pca.df$pc1 %in% selected_x &
                                            pca.df$pc2 %in% selected_y]

        # Update the reactive variable with the selected samples
        selected_data(selected_samples) # Store selected samples

        # Print the selected sample names
        print(selected_samples)
      }
    })

    #Define download handler for selected samples
    output$downloadData <- downloadHandler(
      filename = function() {
        cat("Filename function triggered\n")  # debug line
        paste("selected_samples_", Sys.Date(), ".txt", sep = "")
      },
      content = function(file) {
        cat("Content function triggered\n")  # debug line
        print(selected_data())
        if (!is.null(selected_data())) {
          write.table(selected_data(), file, sep = "\t", row.names = FALSE, col.names = FALSE)
        }
      }
    )

  })
}

# Run the application
shinyApp(ui = ui, server = server)
