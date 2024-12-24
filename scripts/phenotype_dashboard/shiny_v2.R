#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#


library(shiny)
library(ggplot2)
library(dplyr)
library(cluster)  #for clustering
library(reshape2) #for data manipulation 

install.packages("DBI")
install.packages("RMySQL")
library(DBI)
library(RMySQL)

# Define UI
ui <- fluidPage(
  titlePanel("Mouse Phenotype Visualisation Dashboard"),
  sidebarLayout(
    sidebarPanel(
      # We’ll create conditional panels to show inputs 
      # only relevant for the selected tab.
      
      # 1) Inputs for 'Phenotype Scores for Knockout Mouse'
      conditionalPanel(
        condition = "input.tabs == 'mouse_phenotype_tab'",
        selectInput("selected_mouse", "Select Knockout Mouse:",
                    choices = NULL, selected = NULL),
        sliderInput("mouse_threshold", "Significance Threshold (p-value):",
                    min = 0, max = 1, value = 0.05, step = 0.01,
                    helpText = "Highlight phenotypes below this p-value")
      ),
      
      # 2) Inputs for 'Statistical scores of all knockout mice for a selected phenotype'
      conditionalPanel(
        condition = "input.tabs == 'phenotype_mice_tab'",
        selectizeInput("phenotype", "Select Phenotype:",
          choices = NULL,       # Populated in server
          selected = NULL, multiple = FALSE,
          options = list(maxOptions = 1000),  # or however many we want to show at once
          server = TRUE
        # If you want an optional gene filter, add another select or textInput here
      ),
      
      # 3) Inputs for 'Gene Clusters'
      conditionalPanel(
        condition = "input.tabs == 'clusters_tab'",
        selectInput("cluster_method", "Clustering Method:",
                    choices = c("Hierarchical", "K-Means", "PCA"),
                    selected = "Hierarchical"),
        numericInput("num_clusters", "Number of Clusters (K-Means):",
                     value = 3, min = 2, max = 10, step = 1),
        selectInput("gene_subset", "Subset of Genes:",
                    choices = c("All genes", 
                                "Genes with significant phenotypes (p<0.05)", 
                                "User-specific genes"),
                    selected = "All genes"),
        # If user-specific genes are allowed, you might include a multi-select input:
        textInput("user_genes", 
                  "Enter gene symbols (comma-separated) if 'User-specific genes' is selected:")
      )
    ),
    mainPanel(
      tabsetPanel(
        id = "tabs",   # important for conditionalPanel logic
        # 1) Phenotype Scores for a Single Knockout Mouse
        tabPanel(
          title = "Phenotype Scores for Knockout Mouse",
          value = "mouse_phenotype_tab",
          plotOutput("mouse_phenotype_plot"),
          downloadButton("download_mouse_data", "Download Mouse Data")
        ),
        
        # 2) Scores of All Knockout Mice for a Selected Phenotype
        tabPanel(
          title = "Knockout Mice for Phenotype",
          value = "phenotype_mice_tab",
          plotOutput("phenotype_mouse_plot"),
          downloadButton("download_phenotype_data", "Download Phenotype Data")
        ),
        
        # 3) Gene Clusters
        tabPanel(
          title = "Gene Clusters",
          value = "clusters_tab",
          plotOutput("gene_cluster_plot"),
          downloadButton("download_cluster_data", "Download Cluster Data")
        )
      )
    )
    )
  ),

# Define Server
server <- function(input, output, session) {
  
  # Connect to MySQL database
  # NOTE: It's recommended to store credentials in environment variables
  con <- dbConnect(
    RMySQL::MySQL(),
    dbname   = "IMPCDb",
    host     = "localhost",
    port     = 3306,
    user     = "root",
    password = "password"
  )
  
  onStop(function() {
    dbDisconnect(con)
  })
  
  #------------------------------------------------------
  # Populate dropdowns
  #------------------------------------------------------
  
  # For knockout mice
  observe({
    # Suppose your Genes table has: gene_accession_id or gene_symbol 
    # Let’s assume you want to use gene_accession_id uniquely
    mouse_choices <- dbGetQuery(con, "SELECT DISTINCT gene_accession_id FROM Genes;")
    updateSelectInput(session, "selected_mouse",
                      choices = mouse_choices$gene_accession_id)
  })
  
  # For phenotypes
  observe({
    pheno_choices <- dbGetQuery(con, "SELECT DISTINCT parameter_name FROM Parameters;")
    updateSelectInput(session, "selected_phenotype",
                      choices = pheno_choices$parameter_name)
  })
  
  #------------------------------------------------------
  # 1) Phenotype Scores for a Single Knockout Mouse
  #------------------------------------------------------
  output$mouse_phenotype_plot <- renderPlot({
    req(input$selected_mouse)
    
    # Get data for the selected mouse
    query <- sprintf("
      SELECT A.p_value, P.parameter_name 
      FROM Analyses A
      JOIN Parameters P ON A.parameter_id = P.parameter_id
      JOIN Genes G ON A.gene_id = G.gene_id
      WHERE G.gene_accession_id = '%s'
      ORDER BY A.p_value ASC;",
                     input$selected_mouse
    )
    data <- dbGetQuery(con, query)
    
    if (nrow(data) == 0) {
      plot.new()
      title("No data for the selected mouse.")
      return()
    }
    
    # Optional color-coding:
    # If p_value < input$mouse_threshold => color differently
    data <- data %>%
      mutate(color_flag = ifelse(p_value < input$mouse_threshold, "Significant", "Not Sig"))
    
    ggplot(data, aes(x = reorder(parameter_name, p_value), y = p_value, fill = color_flag)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("Significant" = "red", "Not Sig" = "blue")) +
      labs(
        title = paste("Phenotype Scores for", input$selected_mouse),
        x = "Phenotype",
        y = "p-value",
        fill = "Group"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$download_mouse_data <- downloadHandler(
    filename = function() {
      paste("knockout_mouse_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(input$selected_mouse)
      query <- sprintf("
        SELECT A.p_value, P.parameter_name 
        FROM Analyses A
        JOIN Parameters P ON A.parameter_id = P.parameter_id
        JOIN Genes G ON A.gene_id = G.gene_id
        WHERE G.gene_accession_id = '%s'
        ORDER BY A.p_value ASC;",
                       input$selected_mouse
      )
      data <- dbGetQuery(con, query)
      write.csv(data, file, row.names = FALSE)
    }
  )
  
  #------------------------------------------------------
  # 2) Scores of All Knockout Mice for a Selected Phenotype
  #------------------------------------------------------
  output$phenotype_mouse_plot <- renderPlot({
    req(input$selected_phenotype)
    
    query <- sprintf("
      SELECT A.p_value, G.gene_accession_id 
      FROM Analyses A
      JOIN Parameters P ON A.parameter_id = P.parameter_id
      JOIN Genes G ON A.gene_id = G.gene_id
      WHERE P.parameter_name = '%s'
      ORDER BY A.p_value ASC;",
                     input$selected_phenotype
    )
    data <- dbGetQuery(con, query)
    
    if (nrow(data) == 0) {
      plot.new()
      title("No data available for this phenotype.")
      return()
    }
    
    # We create a scatter plot
    # Color or shape for significance
    data <- data %>%
      mutate(sig_flag = ifelse(p_value < 0.05, "Significant", "Not Sig"))
    
    ggplot(data, aes(x = gene_accession_id, y = p_value, color = sig_flag)) +
      geom_point(size = 3, alpha = 0.7) +
      labs(
        title = paste("Scores of All Knockout Mice for Phenotype:", input$selected_phenotype),
        x = "Knockout Mice (gene_accession_id)",
        y = "p-value",
        color = "Status"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  output$download_phenotype_data <- downloadHandler(
    filename = function() {
      paste("phenotype_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(input$selected_phenotype)
      query <- sprintf("
        SELECT A.p_value, G.gene_accession_id 
        FROM Analyses A
        JOIN Parameters P ON A.parameter_id = P.parameter_id
        JOIN Genes G ON A.gene_id = G.gene_id
        WHERE P.parameter_name = '%s'
        ORDER BY A.p_value ASC;", input$selected_phenotype)
      data <- dbGetQuery(con, query)
      write.csv(data, file, row.names = FALSE)
    }
  )
  
  #------------------------------------------------------
  # 3) Gene Clusters
  #------------------------------------------------------
  output$gene_cluster_plot <- renderPlot({
    
    # Step 1: Query all relevant data
    query <- "
      SELECT G.gene_accession_id, P.parameter_name, A.p_value
      FROM Analyses A
      JOIN Parameters P ON A.parameter_id = P.parameter_id
      JOIN Genes G ON A.gene_id = G.gene_id;
    "
    data <- dbGetQuery(con, query)
    
    if (nrow(data) == 0) {
      plot.new()
      title("No data available for clustering.")
      return()
    }
    
    # Step 2: Subset genes if required
    # "All genes" = no filter
    # "Genes with significant phenotypes (p<0.05)" => keep only rows with p_value < 0.05
    # "User-specific genes" => parse input$user_genes
    if (input$gene_subset == "Genes with significant phenotypes (p<0.05)") {
      data <- data %>% filter(p_value < 0.05)
      # This will cause duplication if the same gene appears multiple times with multiple phenotypes.
      # You might want to define "significant if ANY phenotype < 0.05" => group_by(gene_accession_id)
    } else if (input$gene_subset == "User-specific genes") {
      # parse user_genes as comma separated
      user_list <- unlist(strsplit(input$user_genes, "\\s*,\\s*"))
      data <- data %>% filter(gene_accession_id %in% user_list)
    }
    
    if (nrow(data) == 0) {
      plot.new()
      title("No genes left after the subset selection.")
      return()
    }
    
    # Step 3: Pivot the data (gene by parameters -> p_value)
    cluster_data <- dcast(data, gene_accession_id ~ parameter_name, value.var = "p_value", fill = 0)
    
    # Step 4: Perform clustering based on input$cluster_method
    # We exclude the 'gene_accession_id' column
    mat <- cluster_data[, -1, drop = FALSE]
    
    if (input$cluster_method == "Hierarchical") {
      dist_matrix <- dist(mat)
      hc <- hclust(dist_matrix, method = "ward.D")
      plot(hc, main = "Hierarchical Clustering of Genes",
           xlab = "Genes", ylab = "Distance", sub = "")
      
    } else if (input$cluster_method == "K-Means") {
      set.seed(123)  # for reproducibility
      km <- kmeans(mat, centers = input$num_clusters, nstart = 5)
      
      # A simple 2D plot if the data can be reduced or if we want to do a principal component approach
      # For example, do PCA for visualization:
      pca <- prcomp(mat, scale. = TRUE)
      pca_data <- data.frame(pca$x[, 1:2], cluster = factor(km$cluster),
                             gene = cluster_data$gene_accession_id)
      
      ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster, label = gene)) +
        geom_point(size = 3) +
        geom_text(hjust = 1.2, vjust = 0) +
        labs(title = paste("K-Means Clustering (K=", input$num_clusters, ") + PCA for Genes", sep=""))
      
    } else if (input$cluster_method == "PCA") {
      # Purely do PCA for dimension reduction
      pca <- prcomp(mat, scale. = TRUE)
      pca_data <- data.frame(pca$x[, 1:2],
                             gene = cluster_data$gene_accession_id)
      
      ggplot(pca_data, aes(x = PC1, y = PC2, label = gene)) +
        geom_point(size = 3) +
        geom_text(hjust = 1.2, vjust = 0) +
        labs(title = "PCA of Genes by Phenotype Scores")
    }
  })
  
  output$download_cluster_data <- downloadHandler(
    filename = function() {
      paste("gene_cluster_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      query <- "
        SELECT G.gene_accession_id, P.parameter_name, A.p_value
        FROM Analyses A
        JOIN Parameters P ON A.parameter_id = P.parameter_id
        JOIN Genes G ON A.gene_id = G.gene_id;
      "
      data <- dbGetQuery(con, query)
      write.csv(data, file, row.names = FALSE)
    }
  )
}

#Run Shiny App
shinyApp(ui = ui, server = server)



#parse(file = "/Users/sanjanasrinivasan/Desktop/DCDM_IMPC_Project/scripts/phenotype_dashboard/app.R")






