library(shiny)
library(ggplot2)
library(dplyr)
library(cluster)  # For clustering
library(reshape2) # For data manipulation
library(DBI)
library(RMySQL)

# Define UI
ui <- fluidPage(
  titlePanel("Mouse Phenotype Visualisation Dashboard"),
  
  sidebarLayout(
    sidebarPanel(
      textInput("gene_search", "Search Gene Symbol:", ""),
      
      # Conditional inputs for each tab
      conditionalPanel(
        condition = "input.tabs == 'mouse_phenotype_tab'",
        selectInput("selected_mouse", "Select Knockout Mouse:",
                    choices = NULL, selected = NULL),
        sliderInput("mouse_threshold", "Significance Threshold (p-value):",
                    min = 0, max = 1, value = 0.05, step = 0.01)
      ),
      
      conditionalPanel(
        condition = "input.tabs == 'phenotype_mice_tab'",
        selectInput("selected_phenotype", "Select Phenotype:",
                    choices = NULL, selected = NULL)
      ),
      
      conditionalPanel(
        condition = "input.tabs == 'clusters_tab'",
        selectInput("cluster_method", "Clustering Method:",
                    choices = c("Hierarchical", "K-Means", "PCA"),
                    selected = "Hierarchical"),
        numericInput("num_clusters", "Number of Clusters (K-Means):",
                     value = 3, min = 2, max = 50, step = 1),
        selectInput("gene_subset", "Subset of Genes:",
                    choices = c("All genes", 
                                "Genes with significant phenotypes (p<0.05)", 
                                "User-specific genes"),
                    selected = "All genes"),
        textInput("user_genes", 
                  "Enter gene symbols (comma-separated) if 'User-specific genes' is selected:")
      ),
      
      conditionalPanel(
        condition = "input.tabs == 'gene_disease_tab'",
        selectInput("disease", "Select Disease:", choices = NULL),
        uiOutput("gene_select_ui")  # Dynamically generated selectInput
      ),
    ),
    
    mainPanel(
      tabsetPanel(
        id = "tabs",
        
        tabPanel("Phenotype Scores for Knockout Mouse",
                 value = "mouse_phenotype_tab",
                 plotOutput("mouse_phenotype_plot"),
                 downloadButton("download_mouse_data", "Download Mouse Data")),
        
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
    gene_choices <- dbGetQuery(con, "SELECT DISTINCT gene_id FROM Genes;")
    updateSelectInput(session, "selected_mouse", choices = gene_choices$gene_id)
  })
  
  observe({
    phenotype_choices <- dbGetQuery(con, "SELECT DISTINCT parameter_name FROM Parameters;")
    updateSelectInput(session, "selected_phenotype", choices = phenotype_choices$parameter_name)
  })
  
  observe({
    diseases <- dbGetQuery(con, "SELECT DISTINCT disease_term FROM Diseases;")
    updateSelectInput(session, "disease", choices = diseases$disease_term)
  })
  
  output$gene_select_ui <- renderUI({
    gene_choices <- dbGetQuery(con, "SELECT DISTINCT gene_symbol FROM Genes;")
    req(nrow(gene_choices) > 0)
    selectInput("selected_mouse", "Select Gene:", choices = gene_choices$gene_symbol)
  })
  
  # Search Gene Symbol functionality
  observe({
    gene_search_results <- dbGetQuery(con, paste0(
      "SELECT DISTINCT gene_symbol FROM Genes WHERE gene_symbol LIKE '%", 
      input$gene_search, "%'"))
    updateSelectInput(session, "mouse_id", choices = gene_search_results$gene_symbol)
  })
  
  # Visualisation 1: Phenotype Scores for a Single Knockout Mouse
  output$mouse_phenotype_plot <- renderPlot({
    req(input$selected_mouse)
    
    query <- sprintf("
      SELECT A.p_value, P.parameter_name 
      FROM Analyses A
      JOIN Parameters P ON A.parameter_id = P.parameter_id
      WHERE A.gene_id IN (
        SELECT gene_id FROM Genes WHERE gene_id = '%s'
      ) AND A.p_value <= %f
      ORDER BY A.p_value ASC;",
                     input$selected_mouse, input$mouse_threshold)
    
    data <- dbGetQuery(con, query)
    
    if (nrow(data) == 0) {
      plot.new()
      title("No phenotypes meet the threshold for this knockout mouse.")
      return()
    }
    
    data <- data %>%
      mutate(color_flag = ifelse(p_value < 0.05, "Significant", "Not Significant"))
    
    ggplot(data, aes(x = reorder(parameter_name, p_value), y = p_value, fill = color_flag)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c("Significant" = "red", "Not Significant" = "blue")) +
      labs(title = paste("Phenotype Scores for", input$selected_mouse),
           x = "Phenotypes", y = "p-value") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Visualisation 2: Scores of All Knockout Mice for a Selected Phenotype
  output$phenotype_mouse_plot <- renderPlot({
    req(input$selected_phenotype)
    
    query <- sprintf("
      SELECT A.p_value, G.gene_id 
      FROM Analyses A
      JOIN Parameters P ON A.parameter_id = P.parameter_id
      JOIN Genes G ON A.gene_id = G.gene_id
      WHERE P.parameter_name = '%s'
      ORDER BY A.p_value ASC;",
                     input$selected_phenotype)
    
    data <- dbGetQuery(con, query)
    
    if (nrow(data) == 0) {
      plot.new()
      title("No data available for this phenotype.")
      return()
    }
    
    data <- data %>%
      mutate(sig_flag = ifelse(p_value < 0.05, "Significant", "Not Significant"))
    
    ggplot(data, aes(x = gene_id, y = p_value, color = sig_flag)) +
      geom_point(size = 3) +
      labs(title = paste("Scores of All Knockout Mice for Phenotype:", input$selected_phenotype),
           x = "Knockout Mice", y = "p-value") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Gene-Disease Associations Plot
  
  output$gene_disease_plot <- renderPlot({
    req(input$selected_mouse)  # Ensure that the gene is selected
    
    # Query to fetch the gene-disease associations
    query <- sprintf("
        SELECT Genes.gene_symbol, Diseases.disease_term, PhenodigmScores.phenodigm_score
        FROM PhenodigmScores
        JOIN Genes ON PhenodigmScores.gene_id = Genes.gene_id
        JOIN Diseases ON PhenodigmScores.disease_id = Diseases.disease_id
        WHERE Genes.gene_symbol = '%s'", input$selected_mouse)
    
    # Fetch data from the database
    data <- dbGetQuery(con, query)
    if (nrow(data) == 0) {
      plot.new()
      title("No Gene-Disease association data available for this gene.")
      return()
    }
    
    # Create a bar plot using ggplot2
    ggplot(data, aes(x = reorder(disease_term, phenodigm_score), y = phenodigm_score, fill = disease_term)) +
      geom_bar(stat = "identity") +
      labs(title = paste("Gene-Disease Associations for", input$selected_mouse),
           x = "Disease Term", y = "Phenodigm Score") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  
  # Visualisation 3: Gene Clusters
  output$gene_cluster_plot <- renderPlot({
    query <- "
      SELECT G.gene_id, P.parameter_name, A.p_value
      FROM Analyses A
      JOIN Parameters P ON A.parameter_id = P.parameter_id
      JOIN Genes G ON A.gene_id = G.gene_id;"
    
    data <- dbGetQuery(con, query)
    
    if (nrow(data) == 0) {
      plot.new()
      title("No data available for clustering.")
      return()
    }
    
    if (input$gene_subset == "Genes with significant phenotypes (p<0.05)") {
      data <- data %>% filter(p_value < 0.05)
    } else if (input$gene_subset == "User-specific genes") {
      user_genes <- unlist(strsplit(input$user_genes, "\\s*,\\s*"))
      data <- data %>% filter(gene_id %in% user_genes)
    }
    
    cluster_data <- dcast(data, gene_id ~ parameter_name, value.var = "p_value", fill = 0)
    mat <- cluster_data[, -1]
    
    if (input$cluster_method == "Hierarchical") {
      hc <- hclust(dist(mat), method = "ward.D")
      plot(hc, main = "Hierarchical Clustering of Genes", xlab = "Genes", ylab = "Distance")
    } else if (input$cluster_method == "K-Means") {
      km <- kmeans(mat, centers = input$num_clusters)
      ggplot(data.frame(cluster = km$cluster), aes(x = seq_along(cluster), y = cluster)) +
        geom_point() + labs(title = "K-Means Clustering")
    } else if (input$cluster_method == "PCA") {
      pca <- prcomp(mat, scale. = TRUE)
      ggplot(data.frame(pca$x), aes(x = PC1, y = PC2)) + geom_point() +
        labs(title = "PCA Clustering")
    }
  })
}

shinyApp(ui = ui, server = server)
