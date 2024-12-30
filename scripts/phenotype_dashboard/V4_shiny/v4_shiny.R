library(shiny)
library(ggplot2)
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
                    choices = c("Bar Plot", "Dot Plot"), selected = "Bar Plot")
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
                 plotOutput("mouse_phenotype_plot", height = "900px", width = "150%"),
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
    diseases <- dbGetQuery(con, "SELECT DISTINCT disease_term FROM Diseases;")
    updateSelectInput(session, "disease", choices = diseases$disease_term)
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
  
  output$mouse_phenotype_plot <- renderPlot({
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
  AND P.parameter_name IS NOT NULL 
  ORDER BY A.p_value ASC;",
                     input$selected_mouse)
    
    # Fetch the data from the database
    data <- dbGetQuery(con, query)
    
    # Ensure `data` is a dataframe and exclude NA values from `parameter_name`
    if (!is.data.frame(data) || nrow(data) == 0) {
      plot.new()
      title("No phenotypes meet the threshold for this knockout mouse.")
      return()
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
    
    # Plot based on the selected plot type
    if (input$plot_type == "Bar Plot") {
      # Bar Plot
      ggplot(data, aes(x = reorder(parameter_name, p_value), y = p_value, fill = Threshold)) +
        geom_bar(stat = "identity", position = position_dodge(width = 2), show.legend = TRUE, width = 0.6) +
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
          axis.text.y = element_text(size = 10),
        ) +
        geom_hline(yintercept = input$mouse_threshold, linetype = "dashed", color = "black")
      
    } else if (input$plot_type == "Dot Plot") {
      
      p <- ggplot(data, aes(x = reorder(parameter_name, p_value), y = p_value, color = Threshold, text = parameter_name)) +
        geom_point(size = 4) +  # Adjust size of the points
        scale_color_manual(values = c("Significant" = "palegreen3", "Not Significant" = "indianred3")) + 
        labs(
          title = paste("The Phenotype Scores for Knockout Mouse:", input$selected_mouse),
          subtitle = paste("Showing phenotypes with p-value <= ", input$mouse_threshold),
          x = NULL,  # Remove x-axis label
          y = "p-value for Phenotype Association"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_blank(),  # Remove x-axis labels
          axis.title.x = element_blank(),  # Remove x-axis title
          axis.title.y = element_text(size = 12, face = "bold"),  
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
          plot.subtitle = element_text(size = 12, hjust = 0.5),  
          axis.text.y = element_text(size = 10),
        ) +
        geom_hline(yintercept = input$mouse_threshold, linetype = "dashed", color = "black")
      
      # If you want to check the plot without interactivity, comment out ggplotly and use print(p) directly
      print(p)  # Check the plot without plotly to see if it renders normally
      
      # Make the plot interactive using plotly
      ggplotly(p, tooltip = "text")  # This will allow you to hover and see the parameter name
    }
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
  
  
  #Visualisation 3  
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
      # Keep rows (genes) where ANY parameter’s p-value is < 0.05
      keep_rows <- apply(data_wide, 1, function(x) any(x < 0.05, na.rm = TRUE))
      data_wide <- data_wide[keep_rows, , drop = FALSE]
      
    } else if (input$gene_subset == "User-specific genes") {
      user_genes <- unlist(strsplit(input$user_genes, "\\s*,\\s*"))
      # rownames(data_wide) are the gene_accession_ids
      data_wide <- data_wide[rownames(data_wide) %in% user_genes, , drop = FALSE]
    }
    
    # Check if anything remains
    if (nrow(data_wide) == 0) {
      plot.new()
      title("No genes left after filtering.")
      return()
    }
    
    # Scale the data
    mat_scaled <- scale(data_wide)
    
    
    #Clustering
    if (input$cluster_method == "Hierarchical") {
      
      # Hierarchical clustering
      hc <- hclust(dist(mat_scaled), method = "ward.D2")
      plot(as.dendrogram(hc),
           main = "Hierarchical Clustering of Genes",
           xlab = "Genes", ylab = "Distance", cex = 0.7)
      
    } else if (input$cluster_method == "PCA") {
      
      # Run PCA
      pca <- prcomp(mat_scaled, scale. = FALSE)  # mat_scaled is already scaled
      pca_data <- data.frame(pca$x[, 1:2])
      pca_data$gene <- rownames(mat_scaled)
      
      # K-means clustering for coloring
      km <- kmeans(mat_scaled, centers = input$num_clusters)
      pca_data$cluster <- factor(km$cluster)
      
      # Title
      gene_subset_label <- switch(
        input$gene_subset,
        "All genes" = "All Genes",
        "Genes with significant phenotypes (p<0.05)" = "Significant Genes",
        "User-specific genes" = "User-Selected Genes"
      )
      plot_title <- paste("PCA Clustering of", gene_subset_label)
      
      # Plot
      ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster, label = gene)) +
        geom_point(size = 3, alpha = 0.8) +
        scale_color_manual(values = rainbow(input$num_clusters)) +
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
          panel.border = element_rect(color = "black", fill = NA, size = 1.5),
          panel.grid.major = element_line(size = 0.5, linetype = "dotted", color = "gray80"),
          panel.grid.minor = element_blank()
        )
      
    } else if (input$cluster_method == "UMAP") {
      
      # UMAP
      umap_result <- umap(mat_scaled, n_neighbors = 15, min_dist = 0.1)
      umap_data <- data.frame(
        UMAP1 = umap_result$layout[, 1],
        UMAP2 = umap_result$layout[, 2],
        gene  = rownames(mat_scaled)
      )
      
      # K-means for color
      km <- kmeans(mat_scaled, centers = input$num_clusters)
      umap_data$cluster <- factor(km$cluster)
      
      # Title
      gene_subset_label <- switch(
        input$gene_subset,
        "All genes" = "All Genes",
        "Genes with significant phenotypes (p<0.05)" = "Significant Genes",
        "User-specific genes" = "User-Selected Genes"
      )
      
      ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = cluster, label = gene)) +
        geom_point(size = 3, alpha = 0.8) +
        scale_color_manual(values = rainbow(input$num_clusters)) +
        labs(
          title = paste("UMAP Clustering of", gene_subset_label),
          x = "UMAP 1",
          y = "UMAP 2",
          color = "Cluster"
        ) +
        theme_minimal() +
        theme(
          plot.title   = element_text(size = 16, face = "bold", hjust = 0.5),
          axis.title   = element_text(size = 12, face = "bold"),
          axis.text    = element_text(size = 10),
          panel.border = element_rect(color = "black", fill = NA, size = 1.5)
        )
    }
  })
}

shinyApp(ui = ui, server = server)
