library(shiny)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(RCurl)
library(plotly)
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
           downloadButton('downloadData', 'Download selected samples')
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
    edf.orig <- as.data.frame(lapply(Normalized_Counts(),as.numeric))
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
      arrange(match(.data[[col_file3]], samples))
    edf.filt <- edf.filt %>% select(Sample.df[[col_file3]])
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
