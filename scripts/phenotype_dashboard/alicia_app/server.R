# Install necessary packages if not already installed
# install.packages(c("shiny", "plotly", "cluster", "reshape2", "DBI", "RMySQL", "factoextra", "umap", "dendextend", "ggdendro"))

# Load libraries
library(shiny)
library(plotly)
library(cluster)
library(reshape2)
library(DBI)
library(RMySQL)
library(factoextra)
library(umap)
library(dendextend)
library(ggdendro)

# Define Server Logic
server <- function(input, output, session) {
  
  # Establish a connection to the MySQL database
  con <- dbConnect(
    RMySQL::MySQL(),
    dbname = "IMPCDb",
    host = "localhost",
    port = 3306,
    user = "root",
    password = "KCL2024!"
  )
  
  # Ensure the database connection is closed when the app stops
  onStop(function() {
    dbDisconnect(con)
  })
  
  # Dynamically populate dropdowns for Figure 1
  observe({
    # Populate knockout mouse options
    gene_choices <- dbGetQuery(con, "SELECT DISTINCT gene_symbol FROM Genes ORDER BY gene_symbol ASC;")
    if (nrow(gene_choices) > 0) {
      updateSelectInput(session, "genotype_mouse", choices = gene_choices$gene_symbol)
    }
    
    # Populate mouse strain options
    mouse_strain <- dbGetQuery(con, "SELECT DISTINCT mouse_strain FROM Analyses;")
    if (nrow(mouse_strain) > 0) {
      updateSelectInput(session, "genotype_mouse_strain", choices = c("All", mouse_strain$mouse_strain))
    }
    
    # Populate life stage options
    life_stage <- dbGetQuery(con, "SELECT DISTINCT mouse_life_stage FROM Analyses;")
    if (nrow(life_stage) > 0) {
      updateSelectInput(session, "genotype_life_stage", choices = c("All", life_stage$mouse_life_stage))
    }
  })
  
  # Dynamically populate dropdowns for Figure 2
  observe({
    # Populate phenotype group options
    groups <- dbGetQuery(con, "SELECT DISTINCT group_name FROM Groupings ORDER BY group_name ASC;")
    groups$group_name <- str_to_title(groups$group_name)
    updateSelectInput(session, "phenotype_group", choices = groups$group_name)
  })
  
  observe({
    # Populate phenotypes based on the selected group
    req(input$phenotype_group)
    phenotype_data <- dbGetQuery(con, sprintf("
      SELECT P.parameter_name, COUNT(A.p_value) AS data_count
      FROM Parameters P
      JOIN ParameterGroupings PG ON P.parameter_id = PG.parameter_id
      LEFT JOIN Analyses A ON P.parameter_id = A.parameter_id
      JOIN Groupings GR ON PG.group_id = GR.group_id
      WHERE GR.group_name = '%s'
      GROUP BY P.parameter_name ORDER BY P.parameter_name ASC;",
                                              input$phenotype_group))
    
    available_data <- phenotype_data %>%
      filter(data_count > 0) %>%
      pull(parameter_name)
    
    if (length(available_data) > 0) {
      updateSelectInput(session, "phenotype", choices = available_data)
      output$phenotype_explanation <- renderText("Only phenotypes with data available are shown.")
    } else {
      updateSelectInput(session, "phenotype", choices = c("No data available"))
      output$phenotype_explanation <- renderText("No data is available for the selected parameter group.")
    }
  })
  
  # Dynamically populate dropdowns for Figure 3
  observe({
    # Populate gene options for user-specific clustering
    all_genes <- dbGetQuery(con, "SELECT DISTINCT gene_symbol FROM Genes ORDER BY gene_symbol ASC;")
    updateSelectizeInput(session, "user_genes", choices = all_genes$gene_symbol, server = TRUE)
  })
  
  # Visualization 1: Statistical Scores for Selected Knockout Mouse
  
  # Render the UI container based on the selected plot type
  output$plot_container <- renderUI({
    if (input$genotype_plot_type == "All Phenotypes") {
      div(
        style = "overflow-x: auto; overflow-y: hidden; height: 700px;",
        plotlyOutput("mouse_genotype_plot", width = "5000px", height = "100%")
      )
    } else {
      plotlyOutput("mouse_genotype_plot", width = "100%", height = "700px")
    }
  })
  
  # Render the plot for Figure 1: Statistical Scores for Selected Knockout Mouse
  output$mouse_genotype_plot <- renderPlotly({
    req(input$genotype_mouse)
    
    base_query <- sprintf("
    SELECT A.p_value, A.parameter_id, P.parameter_name 
    FROM Analyses A
    JOIN Parameters P ON A.parameter_id = P.parameter_id
    WHERE A.gene_accession_id IN (
      SELECT gene_accession_id FROM Genes WHERE gene_symbol = '%s'
    )
    AND A.p_value IS NOT NULL
    AND A.p_value > 0
    AND P.parameter_name IS NOT NULL",
                          input$genotype_mouse)
    
    if (input$genotype_mouse_strain != "All") {
      base_query <- paste0(base_query, " AND A.mouse_strain = '", input$genotype_mouse_strain, "'")
    }
    if (input$genotype_life_stage != "All") {
      base_query <- paste0(base_query, " AND A.mouse_life_stage = '", input$genotype_life_stage, "'")
    }
    
    final_query <- paste0(base_query, " ORDER BY A.p_value ASC;")
    
    data <- tryCatch(
      dbGetQuery(con, final_query),
      error = function(e) {
        message("Error fetching data: ", e$message)
        return(NULL)
      }
    )
    
    output$no_data_message <- renderUI({
      if (is.null(data) || nrow(data) == 0) {
        tagList(
          h3("No data available for the selected parameters. Please adjust your filters.",
             style = "color: red; text-align: center;")
        )
      } else {
        NULL
      }
    })
    
    if (!is.null(data) && nrow(data) > 0) {
      data <- data %>%
        filter(!is.na(parameter_name) & tolower(parameter_name) != "na") %>%
        group_by(parameter_id) %>%
        summarize(
          parameter_name = str_to_title(first(parameter_name)),
          p_value = mean(p_value),
          Threshold = ifelse(any(p_value < input$genotype_threshold), "Significant", "Not Significant")
        ) %>%
        arrange(p_value)
      
      if (input$genotype_plot_type == "Top 25 Phenotypes") {
        data <- data[1:min(25, nrow(data)), ]
      }
      
      p <- ggplot(data, aes(
        x = reorder(parameter_name, p_value),
        y = -log2(p_value),
        fill = Threshold,
        text = paste0(
          "Parameter Name: ", parameter_name,
          "<br>P-value: ", signif(p_value, digits = 3),
          "<br>Threshold: ", Threshold
        )
      )) +
        geom_bar(stat = "identity", position = position_dodge(width = 2), show.legend = TRUE, width = 0.6) +
        scale_fill_manual(values = c("Significant" = "palegreen3", "Not Significant" = "indianred3")) +
        labs(
          title = paste(input$genotype_plot_type, "for", input$genotype_mouse),
          x = "Knockout Mouse Phenotype",
          y = "Phenotype Significance (-log2(p-value))"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 0, vjust=1, size = 10, face = "bold"),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 15, face = "bold", hjust = 0.5)
        ) +
        geom_hline(yintercept = -log2(input$genotype_threshold), linetype = "dashed", color = "black")
      
      ggplotly(p, tooltip = "text")
    } else {
      NULL
    }
  })
  }
