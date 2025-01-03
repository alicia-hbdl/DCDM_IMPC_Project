library(shiny)
library(ggplot2)
library(plotly)
library(dplyr)
library(cluster)  # For clustering
library(reshape2) # For data manipulation
library(DBI)
library(RMySQL)
library(stringr)
library(factoextra)
library(umap)

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
                    choices = c("All", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")),
        selectInput("selected_phenotype", "Select Phenotype:",
                    choices = NULL, selected = NULL),
        # Display a message explaining why some parameter groups are missing
        textOutput("phenotype_explanation"),
      ),
      
      conditionalPanel(
        condition = "input.tabs == 'clusters_tab'",
        selectInput("cluster_method", "Clustering Method:",
                    choices = c("Hierarchical", "PCA", "UMAP"),
                    selected = "Hierarchical"),
        # Only show number of clusters if PCA or UMAP is chosen
        conditionalPanel(
          condition = "input.cluster_method == 'PCA' || input.cluster_method == 'UMAP'",
          numericInput("num_clusters", "Number of Clusters (K-Means):",
                       value = 3, min = 2, max = 50, step = 1)
        ),
        selectInput("gene_subset", "Subset of Genes:",
                    choices = c("All genes", 
                                "Genes with significant phenotypes (p<0.05)", 
                                "User-specific genes"),
                    selected = "All genes"),
        # Only show user_genes selection if "User-specific genes" is chosen
        conditionalPanel(
          condition = "input.gene_subset == 'User-specific genes'",
          selectizeInput("user_genes", 
                         "Select Gene Symbols:", 
                         choices = NULL,   # We'll update via server
                         multiple = TRUE,
                         options = list(placeholder = 'Select at least 3 genes',
                                        maxOptions = 10))  # tune as needed
        ),
        # We rename the strain/life_stage *for cluster tab* to avoid ID conflicts
        selectInput("cluster_mouse_strain", "Select Mouse Strain:", 
                    choices = NULL, 
                    selected = "All"),
        selectInput("cluster_life_stage",   "Select Mouse Life Stage:", 
                    choices = NULL, 
                    selected = "All")
        )
      ),
    
    mainPanel(
      tabsetPanel(
        id = "tabs",
        
        tabPanel("Phenotype Scores for Knockout Mouse",
                 value = "mouse_phenotype_tab",
                 uiOutput("no_data_message"),
                 plotlyOutput("mouse_phenotype_plot", height = "900px", width = "160%"),
                 downloadButton("download_mouse_data", "Download Mouse Data")),
        
        tabPanel("Knockout Gene Scores for Phenotype",
                 value = "phenotype_mice_tab",
                 plotlyOutput("phenotype_mouse_plot", height = "900px", width = "160%"),
                 downloadButton("download_phenotype_data", "Download Phenotype Data")),
        
        tabPanel("Gene Clusters",
                 value = "clusters_tab",
                 plotlyOutput("gene_cluster_plot", height = "2000px"),
                 downloadButton("download_cluster_data", "Download Cluster Data"))
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
  
  ## Populate dropdowns
  
  # Knockout mouse (phenotype tab)
  observe({
    gene_choices <- dbGetQuery(con, "SELECT DISTINCT gene_symbol FROM Genes;")
    updateSelectInput(session, "selected_mouse", choices = gene_choices$gene_symbol)
  })
  
  # Parameter grouping
  observe({
    groups <- dbGetQuery(con, "SELECT DISTINCT group_id FROM ParameterGroupings;")
    updateSelectInput(session, "selected_phenotype_group", choices = groups$group_id)
  })
  
  # Phenotype choices based on group
  observe({
    req(input$selected_phenotype_group)  # Ensure the input is valid
    
    # Query to fetch phenotypes and their data availability
    query <- sprintf("
    SELECT P.parameter_name, COUNT(A.p_value) AS data_count
    FROM Parameters P
    JOIN ParameterGroupings PG ON P.parameter_id = PG.parameter_id
    LEFT JOIN Analyses A ON P.parameter_id = A.parameter_id
    WHERE PG.group_id = '%s'
    GROUP BY P.parameter_name;", input$selected_phenotype_group)
    
    phenotype_data <- dbGetQuery(con, query)
    
    # Filter to include only parameter names with available data
    available_data <- phenotype_data %>%
      filter(data_count > 0) %>%
      pull(parameter_name)
    
    # Update dropdown based on available data
    if (length(available_data) > 0) {
      updateSelectInput(session, "selected_phenotype", choices = available_data)
      explanation <- "Only phenotypes with data available are shown."
    } else {
      updateSelectInput(session, "selected_phenotype", choices = c("No data available"))
      explanation <- "No data is available for the selected parameter group."
    }
    
    # Update the explanation text
    output$phenotype_explanation <- renderText({
      explanation
    })
  })
  
  # Disease
  observe({
    diseases <- dbGetQuery(con, "SELECT DISTINCT disease_term FROM Diseases;")
    updateSelectInput(session, "disease", choices = diseases$disease_term)
  })
  
  # Mouse strains for phenotype tab
  observe({
    mouse_strains <- dbGetQuery(con, "SELECT DISTINCT mouse_strain FROM Analyses")
    updateSelectInput(session, "mouse_strain", choices = c("All", sort(mouse_strains$mouse_strain)))
  })
  
  # Mouse life stages for phenotype tab
  observe({
    life_stages <- dbGetQuery(con, "SELECT DISTINCT mouse_life_stage FROM Analyses")
    updateSelectInput(session, "life_stage", choices = c("All", sort(life_stages$mouse_life_stage)))
  })
  
  # 5) cluster_mouse_strain & cluster_life_stage (used in Gene Clusters tab)
  observe({
    mouse_strains <- dbGetQuery(con, "SELECT DISTINCT mouse_strain FROM Analyses;")
    updateSelectInput(session, "cluster_mouse_strain", 
                      choices = c("All", sort(mouse_strains$mouse_strain)))
  })
  
  observe({
    life_stages <- dbGetQuery(con, "SELECT DISTINCT mouse_life_stage FROM Analyses;")
    updateSelectInput(session, "cluster_life_stage", 
                      choices = c("All", sort(life_stages$mouse_life_stage)))
  })
  
  # Populate user_genes (for user-specific genes) with all gene symbols
  observe({
    # We'll re-use the Genes table
    all_genes <- dbGetQuery(con, "SELECT DISTINCT gene_symbol FROM Genes;")
    # Update the selectizeInput with all possible gene symbols
    updateSelectizeInput(session, "user_genes", choices = all_genes$gene_symbol, server = TRUE)
  })
  
  # Visualisation 1: Statistical Scores for Selected Knockout Mouse
  
  output$mouse_phenotype_plot <- renderPlotly({
    req(input$selected_mouse)  
    
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
    
    
    data <- data %>%
      group_by(parameter_name) %>%
      summarize(
        p_value = mean(p_value),
        Threshold = ifelse(any(p_value < input$mouse_threshold), "Significant", "Not Significant")
      ) %>%
      ungroup()
    
    p <- if (input$plot_type == "Bar Plot") {
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
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7, face = "bold"),  
          axis.title.x = element_text(size = 12, face = "bold"),  
          axis.title.y = element_text(size = 12, face = "bold"),  
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),  
          plot.subtitle = element_text(size = 12, hjust = 0.5),  
          axis.text.y = element_text(size = 10),
        ) +
        geom_hline(yintercept = input$mouse_threshold, linetype = "dashed", color = "black")
    } else 
      
      if (input$plot_type == "Dot Plot") {
      ggplot(data, aes(x = parameter_name, y = p_value, color = Threshold, text = paste("p-value:", p_value, "<br>Phenotype:", parameter_name))) +
        geom_point(size = 2.5) +
        scale_color_manual(values = c("Significant" = "palegreen3", "Not Significant" = "indianred3")) +
        labs(
          title = paste("The Phenotype Scores for Knockout Mouse:", input$selected_mouse),
          subtitle = paste("Showing phenotypes with p-value <= ", input$mouse_threshold),
          x = "Knockout Mouse Phenotype", 
          y = "p-value for Phenotype Association"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5),
          axis.text.y = element_text(size = 10)
        ) +
        geom_hline(yintercept = input$mouse_threshold, linetype = "dashed", color = "black")
    }
  
    ggplotly(p, tooltip = "text")
  })

  
  # Visualisation 2: Statistical Scores of All Knockout Mice for a Selected Phenotype
  
  output$phenotype_mouse_plot <- renderPlotly({
    req(input$selected_phenotype, input$selected_phenotype_group)
    
    # Query to fetch p-values for the selected phenotype and group
    query <- sprintf("
    SELECT A.p_value, G.gene_symbol
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
    
    # Group and summarize data for bar plot
    data <- data %>%
      group_by(gene_symbol) %>%
      summarize(
        p_value = mean(p_value),
        Threshold = ifelse(any(p_value < input$mouse_threshold), "Significant", "Not Significant")
      ) %>%
      ungroup()
    
        p <- ggplot(data, aes(x = reorder(gene_symbol, p_value), y = p_value, color = Threshold, text = paste("p-value:", p_value, "<br>Gene:", gene_symbol))) +
        geom_point(size = 1.5) +
        scale_color_manual(values = c("Significant" = "palegreen3", "Not Significant" = "indianred3")) +
        labs(
          title = paste("Gene Knockout Scores for Selected Phenotype:", input$selected_phenotype),
          subtitle = paste("Showing genes with p-value <= ", input$mouse_threshold),
          x = "Gene Symbol", 
          y = "p-value for Gene Association"
        ) +
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 4, face = "bold"),
          axis.title.x = element_text(size = 12, face = "bold"),
          axis.title.y = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 12, hjust = 0.5),
          axis.text.y = element_text(size = 10)
        ) +
        geom_hline(yintercept = input$mouse_threshold, linetype = "dashed", color = "black")
  
      ggplotly(p, tooltip = "text")
      })
  
  
  #Visualisation 3  
  analysis_data <- reactive({
    # Base query
    query <- "SELECT gene_accession_id, parameter_id, ROUND(AVG(p_value), 6) AS avg_p_value 
            FROM Analyses 
            WHERE 1=1"
    # Add filters for strain
    if (input$cluster_mouse_strain != "All") {
      query <- paste0(query, " AND mouse_strain = '", input$cluster_mouse_strain, "'")
    }
    # Add filters for strain
    if (input$cluster_life_stage != "All") {
      query <- paste0(query, " AND mouse_life_stage = '", input$cluster_life_stage, "'")
    }
    # Final group by & order
    query <- paste0(query, " GROUP BY gene_accession_id, parameter_id 
                           ORDER BY avg_p_value ASC;")
    # Execute the query
    df <- dbGetQuery(con, query)
    
    df
  })
  
  # pca_matrix -> pivot to wide
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
    
    # Convert gene_accession_id into rownames, assuming 'gene_accession_id' is a column in df
    wide_df <- as.data.frame(wide_df)
    rownames(wide_df) <- wide_df$gene_accession_id
    # Remove the gene_accession_id column now that it’s the rowname
    wide_df <- wide_df[, !names(wide_df) %in% "gene_accession_id"]
    
    wide_df
  })
  
  # Render the cluster plot
  output$gene_cluster_plot <- renderPlotly({
    req(pca_matrix())           # Wait until the data is available
    data_wide <- pca_matrix()   # This is the wide gene-by-parameter matrix
    
    # If no data:
    if (is.null(data_wide) || nrow(data_wide) == 0) {
      plot.new()
      title("No data available after filtering.")
      return()
    }
    
    # Additional filters on subset
    if (input$gene_subset == "Genes with significant phenotypes (p<0.05)") {
      # Keep rows (genes) where ANY parameter’s p-value is < 0.05
      keep_rows <- apply(data_wide, 1, function(x) any(x < 0.05, na.rm = TRUE))
      data_wide <- data_wide[keep_rows, , drop = FALSE]
      
    } else if (input$gene_subset == "User-specific genes") {
      # If user selected fewer than 2 genes, show a message
      if (length(input$user_genes) < 3) {
        plot.new()
        title("Please select at least 2 genes for PCA/UMAP.")
        return(NULL)
      }
      data_wide <- data_wide[rownames(data_wide) %in% input$user_genes, , drop = FALSE]
    }
    
    # Check if anything remains
    if (nrow(data_wide) == 0) {
      plot.new()
      title("No genes left after filtering.")
      return(NULL)
    }
    
    # Scale the data
    mat_scaled <- scale(data_wide)
    
    
    #Clustering
    if (input$cluster_method == "Hierarchical") {
      # from earlier code. Make sure you have at least 2 rows to cluster:
      if (nrow(mat_scaled) < 2) {
        plot.new()
        title("Not enough data for hierarchical clustering (need >= 2 rows).")
        return(NULL)
      }
      # Perform hierarchical clustering
      hc <- hclust(dist(mat_scaled), method = "ward.D")
      
      # 3. Convert hclust -> dendrogram -> apply "hang" to shorten tips
      #    "hang" sets how far tips hang below the rest of the dendrogram.
      #    0.1 or 0.2 often works well. Smaller => shorter vertical lines to labels.
      library(dendextend)  # install.packages("dendextend") if needed
      dend <- as.dendrogram(hc)
      dend <- hang.dendrogram(dend, hang = 0.2)  # tweak 0.2 as desired
      
      # 2) Convert to a ggdendro-friendly structure
      library(ggdendro)   # install.packages("ggdendro") if not installed
      dend_data <- dendro_data(dend, type = "rectangle")
      
      str(dend_data$labels)
      
      set.seed(123)  # For reproducibility
      
      # Generate a color palette for the branches
      branch_colors <- sample(colors(), nrow(dend_data$segments), replace = TRUE)
      
      p <- ggplot() +
        # Segments for the branches
        geom_segment(
          data = dend_data$segments,
          aes(x = x, y = y, xend = xend, yend = yend),
          color = branch_colors,  # Custom color for the branches
          size = 0.8  # Thicker branches for better visibility
        ) +
      
        geom_text(
          data = dend_data$labels %>%
            mutate(y = y - max(dend_data$segments$y) * 0.2),  # Shift labels further down
          aes(x = x, y = y, label = label),
          size = 2.5,
          hjust = 0.5,  # Center the text horizontally
          color = "black"
        ) +
        # Flip coords to mimic the usual horizontal R dendrogram look
        coord_flip() +
        # Reverse y so the tree grows "up" (optional but common)
        scale_y_reverse(expand = c(0.2, 0)) +
        # Minimal theme
        theme_minimal(base_size = 10) +
        labs(
          title = "Hierarchical Clustering of Knockout Mouse Genes",
          x = "Genes of Knockout Mouse",
          y = "Cluster Distance"
        ) +
        theme(
          plot.title   = element_text(size = 16, face = "bold", hjust = 0.5), 
          axis.text.y  = element_text(size = 8),  # Axis label styling
          axis.ticks.y = element_blank(),  # Remove y-axis ticks for a cleaner look
          axis.title.y = element_blank(),  # Remove the y-axis title
          panel.grid.major = element_line(color = "grey90", size = 1),  # Light grid lines for better readability
          panel.grid.minor = element_blank()
        )

    } else if (input$cluster_method == "PCA") {
      
      # Run PCA
      pca <- prcomp(mat_scaled, scale. = FALSE)  # mat_scaled is already scaled
      pca_data <- data.frame(pca$x[, 1:2])
      pca_data$gene <- rownames(mat_scaled)
      
      # K-means clustering 
      km <- kmeans(mat_scaled, centers = input$num_clusters)
      pca_data$cluster <- factor(km$cluster)
      
      # Set dynamic graph title
      gene_subset_label <- switch(
        input$gene_subset,
        "All genes" = "All Genes",
        "Genes with significant phenotypes (p<0.05)" = "Significant Genes",
        "User-specific genes" = "User-Selected Genes"
      )
      plot_title <- paste("PCA Clustering of", gene_subset_label)
      
      # Plot
      p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = cluster, text = paste0("Gene: ", gene, "<br>Cluster: ", cluster))) +
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
      ggplotly(p, tooltip = "text")

      
    } else if (input$cluster_method == "UMAP") {
      # UMAP
      # Make sure mat_scaled is a matrix
      mat_scaled <- as.matrix(mat_scaled)
      
      umap_result <- umap(mat_scaled, n_neighbors = 15, min_dist = 0.1)
      umap_data <- data.frame(
        UMAP1 = umap_result$layout[, 1],
        UMAP2 = umap_result$layout[, 2],
        gene  = rownames(mat_scaled)
      )
      
      # K-means for colouring
      km <- kmeans(mat_scaled, centers = input$num_clusters)
      umap_data$cluster <- factor(km$cluster)
      
      # Set dynamic graph title
      gene_subset_label <- switch(
        input$gene_subset,
        "All genes" = "All Genes",
        "Genes with significant phenotypes (p<0.05)" = "Significant Genes",
        "User-specific genes" = "User-Selected Genes"
      )
      
      p <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = cluster, text = paste0("Gene: ", gene, "<br>Cluster: ", cluster))) +
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
      ggplotly(p, tooltip = "text")
    }
  })
}

shinyApp(ui = ui, server = server)
