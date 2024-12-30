library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
library(cluster)  # For clustering
library(reshape2) # For data manipulation
library(DBI)
library(RMySQL)
library(stringr)
library(uwot)

# Define UI
ui <- fluidPage(
  titlePanel("Statistical Analysis and Visualisation of Knockout Mice"),
  
  sidebarLayout(
    sidebarPanel(
      
      # Conditional inputs for each tab
      conditionalPanel(
        condition = "input.tabs == 'mouse_phenotype_tab'",
        selectInput("selected_mouse", "Select Knockout Mouse:",
                    choices = NULL, selected = NULL),
        sliderInput("mouse_threshold", "Significance Threshold (p-value):",
                    min = 0, max = 1, value = 0.05, step = 0.01),
        selectInput("plot_type", "Select Plot Type:",
                    choices = c("Bar Plot", "Dot Plot"), selected = "Bar Plot"),
        selectInput("mouse_strain", "Select Mouse Strain:",
                    choices = c("All", "129SV", "C57BL")),
        selectInput("life_stage", "Select Mouse Life Stage:",
                    choices = c("All", "E12.5", "E15.5", "E18.5", "E9.5", "Early adult", "Late adult", "Middle aged adult"))
      ),
      
      conditionalPanel(
        condition = "input.tabs == 'phenotype_mice_tab'",
        selectInput("selected_phenotype_group", "Select Parameter Grouping:",
                    choices = NULL, selected = NULL),
        selectInput("selected_phenotype", "Select Phenotype:",
                    choices = NULL, selected = NULL)
        
      ),
      
      conditionalPanel(
        condition = "input.tabs == 'clusters_tab'",
        selectInput("cluster_method", "Clustering Method:",
                    choices = c("Hierarchical", "PCA", "UMAP"),
                    selected = "Hierarchical"),
        numericInput("num_clusters", "Number of Clusters (K-Means):",
                     value = 3, min = 2, max = 50, step = 1),
        selectInput("gene_subset", "Subset of Genes:",
                    choices = c("All genes", 
                                "Genes with significant phenotypes (p<0.05)", 
                                "User-specific genes"),
                    selected = "All genes"),
        textInput("user_genes", 
                  "Enter gene symbols (comma-separated) if 'User-specific genes' is selected:"),
        selectInput("mouse_strain", "Select Mouse Strain:", choices = NULL, selected = "All"),
        selectInput("life_stage", "Select Mouse Life Stage:", choices = NULL, selected = "All")
        )
      ),
      
      conditionalPanel(
        condition = "input.tabs == 'gene_disease_tab'",
        selectInput("disease", "Select Disease:", choices = NULL),
        uiOutput("gene_select_ui")  # Dynamically generated selectInput
      )
    ),
    
    mainPanel(
      tabsetPanel(
        id = "tabs",
        
        tabPanel("Phenotype Scores for Knockout Mouse",
                 value = "mouse_phenotype_tab",
                 uiOutput("no_data_message"),
                 plotlyOutput("mouse_phenotype_plot", height = "900px", width = "150%"),
                 downloadButton("download_mouse_data", "Download Mouse Data"),
                 ),
        
        tabPanel("Knockout Mice for Phenotype",
                 value = "phenotype_mice_tab",
                 plotOutput("phenotype_mouse_plot"),
                 downloadButton("download_phenotype_data", "Download Phenotype Data")),
        
        tabPanel("Gene Clusters",
                 value = "clusters_tab",
                 plotOutput("gene_cluster_plot"),
                 downloadButton("download_cluster_data", "Download Cluster Data")),
        
        tabPanel("Gene-Disease Associations",
                 value = "gene_disease_tab",  
                 plotOutput("gene_disease_plot"),
                 downloadButton("download_gene_disease_data", "Download Gene-Disease Data"))
      )
    )
  )

# Define Server
server <- function(input, output, session) {
  
  # Connect to MySQL database
  con <- dbConnect(
    RMySQL::MySQL(),
    dbname = "IMPCDb",
    host = "localhost",
    port = 3306,
    user = "root",
    password = "mahiat123"
  )
  
  onStop(function() {
    dbDisconnect(con)
  })
  
  # Populate dropdowns
  
  observe({
    gene_choices <- dbGetQuery(con, "SELECT DISTINCT gene_symbol FROM Genes;")
    updateSelectInput(session, "selected_mouse", choices = gene_choices$gene_symbol)
  })
  
  observe({
    groups <- dbGetQuery(con, "SELECT DISTINCT group_id FROM ParameterGroupings;")
    updateSelectInput(session, "selected_phenotype_group", choices = groups$group_id)
  })
  
  observe({
    req(input$selected_phenotype_group)  
    phenotype_choices <- dbGetQuery(con, sprintf("
    SELECT DISTINCT P.parameter_name
    FROM Parameters P
    JOIN ParameterGroupings PG ON P.parameter_id = PG.parameter_id
    WHERE PG.group_id = '%s';", input$selected_phenotype_group))
    updateSelectInput(session, "selected_phenotype", choices = phenotype_choices$parameter_name)
  })
  
  observe({
    gene_choices <- dbGetQuery(con, "SELECT DISTINCT gene_symbol FROM Genes;")
    updateSelectInput(session, "selected_mouse", choices = gene_choices$gene_symbol)
  })
  
  observe({
    mouse_strains <- dbGetQuery(con, "SELECT DISTINCT mouse_strain FROM Analyses")
    updateSelectInput(session, "mouse_strain", choices = c("All", sort(mouse_strains$mouse_strain)))
  })
  
  observe({
    life_stages <- dbGetQuery(con, "SELECT DISTINCT mouse_life_stage FROM Analyses")
    updateSelectInput(session, "life_stage", choices = c("All", sort(life_stages$mouse_life_stage)))
  })
  
  # Visualisation 1: Phenotype Scores for Knockout Mouse
  
  output$mouse_phenotype_plot <- renderPlotly({
    req(input$selected_mouse)  # Ensure a gene symbol is selected
    
    query <- sprintf("
    SELECT A.p_value, P.parameter_name 
    FROM Analyses A
    JOIN Parameters P ON A.parameter_id = P.parameter_id
    WHERE A.gene_accession_id IN (
      SELECT gene_accession_id FROM Genes WHERE gene_symbol = '%s'
    ) 
    AND A.p_value IS NOT NULL 
    AND A.p_value > 0 
    AND P.parameter_name IS NOT NULL",
                     input$selected_mouse)
    
    # Add filters for mouse strain and life stage
    if (input$mouse_strain != "All") {
      query <- paste0(query, " AND A.mouse_strain = '", input$mouse_strain, "'")
    }
    if (input$life_stage != "All") {
      query <- paste0(query, " AND A.mouse_life_stage = '", input$life_stage, "'")
    }
    
    query <- paste0(query, " ORDER BY A.p_value ASC;")
    
    # Fetch the data from the database
    data <- dbGetQuery(con, query)

    if (nrow(data) == 0) {
      output$no_data_message <- renderUI({
        tagList(
          h3("No data available for the selected mouse. Please adjust your filters.", style = "color: red; text-align: center;")
        )
      })
      return(NULL)
    }
    
    # Remove rows where `parameter_name` is NA
    data <- data %>%
      dplyr::filter(!(is.na(parameter_name) | parameter_name == "NA"))
    
    # Aggregate data and calculate thresholds
    data <- data %>%
      group_by(parameter_name) %>%
      summarize(
        p_value = mean(p_value),
        Threshold = ifelse(any(p_value < input$mouse_threshold), "Significant", "Not Significant")
      ) %>%
      ungroup()
    
    # Create the plot based on the user's choice of plot type
    p <- if (input$plot_type == "Bar Plot") {
      ggplot(data, aes(x = reorder(parameter_name, p_value), y = p_value, fill = Threshold, text = paste("p-value:", p_value, "<br>Phenotype:", parameter_name))) +
        geom_bar(stat = "identity", show.legend = FALSE) +
        scale_fill_manual(values = c("Significant" = "palegreen3", "Not Significant" = "indianred3")) +
        labs(
          title = paste("The Phenotype Scores for Knockout Mouse:", input$selected_mouse),
          subtitle = paste("Showing phenotypes with p-value <= ", input$mouse_threshold),
          x = "Knockout Mouse Phenotype", 
          y = "p-value for Phenotype Association"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5),
          axis.text.y = element_text(size = 10)
        )
    } else {
      ggplot(data, aes(x = reorder(parameter_name, p_value), y = p_value, color = Threshold, text = paste("p-value:", p_value, "<br>Phenotype:", parameter_name))) +
        geom_point(size = 4) +
        scale_color_manual(values = c("Significant" = "palegreen3", "Not Significant" = "indianred3")) +
        labs(
          title = paste("The Phenotype Scores for Knockout Mouse:", input$selected_mouse),
          subtitle = paste("Showing phenotypes with p-value <= ", input$mouse_threshold),
          x = "Knockout Mouse Phenotype", 
          y = "p-value for Phenotype Association"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5),
          axis.text.y = element_text(size = 10)
        ) +
        geom_hline(yintercept = input$mouse_threshold, linetype = "dashed", color = "black")
    }
    
    # Convert ggplot to plotly for interactivity (hover functionality)
    ggplotly(p, tooltip = "text")
  })
  
  # Visualisation 2: Scores of All Knockout Mice for a Selected Phenotype
  
  output$phenotype_group_plot <- renderPlot({
    req(input$selected_phenotype, input$selected_phenotype_group)
    
    # Query to fetch p-values for the selected phenotype and group
    query <- sprintf("
    SELECT A.p_value, G.gene_accession_id
    FROM Analyses A
    JOIN Parameters P ON A.parameter_id = P.parameter_id
    JOIN Genes G ON A.gene_accession_id = G.gene_accession_id
    JOIN ParameterGroupings PG ON P.parameter_id = PG.parameter_id
    WHERE P.parameter_name = '%s' AND PG.group_id = '%s'
    ORDER BY A.p_value ASC;", 
                     input$selected_phenotype, input$selected_phenotype_group)
    
    # Fetch data from the database
    data <- dbGetQuery(con, query)
    
    # If no data is available, display a message
    if (nrow(data) == 0) {
      plot.new()
      title("No data available for this selection.")
      return()
    }
    
    # Plotting the p-value for each gene on the x-axis
    ggplot(data, aes(x = gene_accession_id, y = p_value)) +
      geom_point(size = 3, color = "blue") +
      labs(title = paste("Phenotype Data for:", input$selected_phenotype),
           x = "Gene Accessions", y = "p-value") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8),
            axis.text.y = element_text(size = 10)) +
      theme_minimal()
  })
  
  # Visualisation 3  
  analysis_data <- reactive({
    # Base query
    query <- "SELECT gene_accession_id, parameter_id, ROUND(AVG(p_value), 6) AS avg_p_value 
            FROM Analyses 
            WHERE 1=1"
    # Add filters for life stage
    if (input$life_stage != "All") {
      query <- paste0(query, " AND mouse_life_stage = '", input$life_stage, "'")
    }
    # Add filters for strain
    if (input$mouse_strain != "All") {
      query <- paste0(query, " AND mouse_strain = '", input$mouse_strain, "'")
    }
    # Final group by & order
    query <- paste0(query, " GROUP BY gene_accession_id, parameter_id 
                           ORDER BY avg_p_value ASC;")
    # Execute the query
    df <- dbGetQuery(con, query)
    # Return the raw aggregated data
    return(df)
  })
  
  pca_matrix <- reactive({
    df <- analysis_data()
    
    if (nrow(df) == 0) {
      return(NULL)
    }
    
    # Pivot wider: one row per gene, columns are parameter_ids
    wide_df <- df %>%
      tidyr::pivot_wider(
        names_from = parameter_id,
        values_from = avg_p_value
      )
    
    # Convert gene_accession_id into rownames
    # (Assumes 'gene_accession_id' is a column in df)
    wide_df <- as.data.frame(wide_df)
    rownames(wide_df) <- wide_df$gene_accession_id
    
    # Remove the gene_accession_id column now that it’s the rowname
    wide_df <- wide_df[, !names(wide_df) %in% "gene_accession_id"]
    
    return(wide_df)
  })
  
  output$gene_cluster_plot <- renderPlot({
    req(pca_matrix())           # Wait until the data is available
    data_wide <- pca_matrix()   # This is the wide gene-by-parameter matrix
    
    # If no data:
    if (is.null(data_wide) || nrow(data_wide) == 0) {
      plot.new()
      title("No data available after filtering.")
      return()
    }
    
    # --- Filter for gene_subset ---
    if (input$gene_subset == "Genes with significant phenotypes (p<0.05)") {
      # Keep rows (genes) that have p-values < 0.05
      data_wide <- data_wide[rowMeans(data_wide < 0.05, na.rm = TRUE) > 0.2,]
    } else if (input$gene_subset == "User-specific genes") {
      # Filter by user-specified genes
      user_genes <- strsplit(input$user_genes, ",")[[1]]
      data_wide <- data_wide[rownames(data_wide) %in% user_genes, ]
    }
    
    # Run PCA
    pca <- prcomp(data_wide, scale = TRUE)
    
    # --- Clustering ---
    if (input$cluster_method == "Hierarchical") {
      distance_matrix <- dist(pca$x)
      cluster <- hclust(distance_matrix)
      plot(cluster)
    }
    
    else if (input$cluster_method == "PCA") {
      # PCA plot
      pca_data <- data.frame(pca$x[, 1:2])
      ggplot(pca_data, aes(x = PC1, y = PC2)) +
        geom_point() +
        labs(title = "PCA of Gene Expression Data")
    }
    
    else if (input$cluster_method == "UMAP") {
      umap_data <- uwot::umap(pca$x)
      umap_df <- data.frame(umap_data)
      
      ggplot(umap_df, aes(x = V1, y = V2)) +
        geom_point() +
        labs(title = "UMAP of Gene Expression Data")
    }
  })

}

# Run the application
shinyApp(ui = ui, server = server)




