library(shiny)
library(ggplot2)
library(dplyr)
library(cluster)  # For clustering
library(reshape2) # For data manipulation
library(DBI)
library(RMySQL)
library(stringr)


# Define UI
ui <- fluidPage(
  titlePanel("Statistical Analysis and Visualization of Knockout Mouses"),
  
  sidebarLayout(
    sidebarPanel(
      
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
        selectInput("selected_phenotype_group", "Select Parameter Grouping:",
                    choices = NULL, selected = NULL),
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
                 plotOutput("mouse_phenotype_plot", height = "700px"),
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
    password = "Llama123@"
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
      req(input$selected_phenotype_group)  # Ensure the group is selected
      phenotype_choices <- dbGetQuery(con, sprintf("
    SELECT DISTINCT P.parameter_name
    FROM Parameters P
    JOIN ParameterGroupings PG ON P.parameter_id = PG.parameter_id
    WHERE PG.group_id = '%s';", input$selected_phenotype_group))
      updateSelectInput(session, "selected_phenotype", choices = phenotype_choices$parameter_name)
    })
    
  
  observe({
    diseases <- dbGetQuery(con, "SELECT DISTINCT disease_term FROM Diseases;")
    updateSelectInput(session, "disease", choices = diseases$disease_term)
  })

  # Visualisation 1: Phenotype Scores for Knockout Mouse
  
  output$mouse_phenotype_plot <- renderPlot({
    req(input$selected_mouse)  # Ensure a gene symbol is selected
    
    query <- sprintf("
    SELECT A.p_value, P.parameter_name 
    FROM Analyses A
    JOIN Parameters P ON A.parameter_id = P.parameter_id
    WHERE A.gene_accession_id IN (
      SELECT gene_accession_id FROM Genes WHERE gene_symbol = '%s'
    ) AND A.p_value IS NOT NULL AND A.p_value > 0
    ORDER BY A.p_value ASC;",
                     input$selected_mouse, input$mouse_threshold)  # Use the selected gene symbol
    
    data <- dbGetQuery(con, query)
    
    if (nrow(data) == 0) {
      plot.new()
      title("No phenotypes meet the threshold for this knockout mouse.")
      return()
    }
    
    data <- data %>%
      mutate(Threshold = ifelse(p_value < input$mouse_threshold, "Significant", "Not Significant"))
    
    # Adjust the plot height
    ggplot(data, aes(x = reorder(parameter_name, p_value), y = p_value, fill = Threshold)) +
      geom_bar(stat = "identity", show.legend = TRUE, width = 0.7) +  # Bar plot with adjusted bar width
      scale_fill_manual(values = c("Significant" = "palegreen3", "Not Significant" = "indianred3")) + 
      labs(
        title = paste("The Phenotype Scores for Knockout Mouse:", input$selected_mouse),
        subtitle = paste("Showing phenotypes with p-value <= ", input$mouse_threshold),
        x = "Knockout Mouse Phenotype", 
        y = "p-value for Phenotype Association"
      ) +
      theme_minimal() +  # Use minimal theme for a cleaner look
      theme(
        axis.text.x = element_text(angle = 90, hjust = 1, face = "bold"),  
        axis.title.x = element_text(size = 12, face = "bold"),  
        axis.title.y = element_text(size = 12, face = "bold"),  
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
        plot.subtitle = element_text(size = 12, hjust = 0.5),  # Subtitle styling
        axis.text.y = element_text(size = 10)  # Adjust size of y-axis labels 
      ) +
      geom_hline(yintercept = input$mouse_threshold, linetype = "dashed", color = "black") +  # Threshold line
      scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) # Wrap long labels to prevent overlap
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
  
  
  # Gene-Disease Associations Plot
  
  output$gene_disease_plot <- renderPlot({
    req(input$selected_mouse)  # Ensure that the gene is selected
    
    # Query to fetch the gene-disease associations
    query <- sprintf("
        SELECT Genes.gene_symbol, Diseases.disease_term, PhenodigmScores.phenodigm_score
        FROM PhenodigmScores
        JOIN Genes ON PhenodigmScores.gene_accession_id = Genes.gene_accession_id
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
      SELECT G.gene_accession_id, P.parameter_name, A.p_value
      FROM Analyses A
      JOIN Parameters P ON A.parameter_id = P.parameter_id
      JOIN Genes G ON A.gene_accession_id = G.gene_accession_id;"
    
    data <- dbGetQuery(con, query)
    
    if (nrow(data) == 0) {
      plot.new()
      title("No data available for clustering.")
      return()
    }
    
    #Filter subset of genes if specific 
    if (input$gene_subset == "Genes with significant phenotypes (p<0.05)") {
      data <- data %>% filter(p_value < 0.05)
    } else if (input$gene_subset == "User-specific genes") {
      user_genes <- unlist(strsplit(input$user_genes, "\\s*,\\s*"))
      data <- data %>% filter(gene_accession_id %in% user_genes)
    }
    
    #Pivot the data for clustering
    cluster_data <- dcast(data, gene_accession_id ~ parameter_name, value.var = "p_value", fill = 0)
    mat <- cluster_data[, -1]
    
    if (input$cluster_method == "Hierarchical") {
      hc <- hclust(dist(mat), method = "ward.D2")
      ggdendro <- as.dendrogram(hc)
      plot(ggdendro, main = "Hierarchical Clustering of Genes", 
           xlab = "Genes", ylab = "Distance", cex = 0.7)
    } else if (input$cluster_method == "K-Means") {
      km <- kmeans(mat, centers = input$num_clusters)
      kmeans_data <- data.frame(
        Gene = cluster_data$gene_accession_id,
        Cluster = factor(km$cluster)
      )
      ggplot(kmeans_data, aes(x = seq_along(Cluster), y = Cluster, color = Cluster)) +
        geom_point(size = 3) +
        scale_color_manual(values = rainbow(input$num_clusters)) +
        labs(
          title = paste("K-Means Clustering with", input$num_clusters, "Clusters"),
          x = "Gene Index",
          y = "Cluster ID",
          color = "Cluster"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10),
          legend.position = "right"
        )
    } else if (input$cluster_method == "PCA") {
      # PCA computation
      pca <- prcomp(mat, scale. = TRUE)
      pca_data <- data.frame(pca$x[, 1:2], gene = cluster_data$gene_accession_id)
      
      # Add K-Means clustering for color coding
      km <- kmeans(mat, centers = input$num_clusters)
      pca_data$cluster <- factor(km$cluster)  # Add cluster information
      
      # Set dynamic title based on gene subset selection
      gene_subset_label <- switch(input$gene_subset,
                                  "All genes" = "All Genes",
                                  "Genes with significant phenotypes (p<0.05)" = "Significant Genes",
                                  "User-specific genes" = "User-Selected Genes")
      plot_title <- paste("PCA Clustering of", gene_subset_label)
      
      ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster, label = gene)) +
        geom_point(size = 3, alpha = 0.8) +
        scale_color_manual(values = rainbow(input$num_clusters)) +  # Color palette for clusters
        labs(
          title = plot_title,
          x = "Principal Component 1",
          y = "Principal Component 2",
          color = "Cluster"
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10),
          legend.position = "right",
          panel.border = element_rect(color = "black", fill = NA, size = 1.5),  # Add bold border
          panel.grid.major = element_line(size = 0.5, linetype = "dotted", color = "gray80"),
          panel.grid.minor = element_blank()  # Hide minor grids for a cleaner look
        )     }
  })
}

shinyApp(ui = ui, server = server)
