# Install packages
# install.packages(c("shiny", "tidyverse", "plotly", "cluster", "reshape2", "DBI", "RMySQL", "factoextra", "umap", "dendextend", "ggdendro"))

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

# Define UI
ui <- fluidPage(
  titlePanel("Statistical Analysis and Visualisation of Knockout Mice"),
  
  sidebarLayout(
    sidebarPanel(
      
      # Dropdown menu for figure 1 
      conditionalPanel(
        condition = "input.tabs == 'mouse_genotype_tab'",
        selectInput("genotype_mouse", 
                    "Select Knockout Mouse:", 
                    choices = NULL, selected = NULL),  # Dynamically populated
        sliderInput("genotype_threshold", 
                    "Significance Threshold (p-value):", 
                    min = 0, max = 1, value = 0.05, step = 0.01),
        selectInput("genotype_plot_type", 
                    "Select Plot Type:", 
                    choices = c("All Phenotypes", "Top 25 Phenotypes"), selected = "All Phenotypes"),
        selectInput("genotype_mouse_strain", 
                    "Select Mouse Strain:", 
                    choices = NULL, selected = "All"), # Dynamically populated
        selectInput("genotype_life_stage", 
                    "Select Mouse Life Stage:", 
                    choices = NULL, selected = "All")  # Dynamically populated
      ),
      
      # Dropdown menu for figure 2
      conditionalPanel(
        condition = "input.tabs == 'mouse_phenotype_tab'",
        sliderInput("phenotype_threshold", 
                    "Significance Threshold (p-value):", 
                    min = 0, max = 1, value = 0.05, step = 0.01),
        selectInput("phenotype_plot_type", 
                    "Select Plot Type:", 
                    choices = c("All Genotypes", "Top 25 Genotypes"), selected = "All Genotypes"),
        selectInput("phenotype_group", 
                    "Select Phenotype Group:", 
                    choices = NULL, selected = NULL),  # Dynamically populated
        selectInput("phenotype", 
                    "Select Phenotype:", 
                    choices = NULL, selected = NULL),  # Dynamically populated
        textOutput("phenotype_explanation")
      ),
      
      # Dropdown menu for figure 3
      conditionalPanel(
        condition = "input.tabs == 'clustering_tab'",
        selectInput("cluster_method", "Clustering Method:", choices = c("Hierarchical", "PCA", "UMAP"), selected = "Hierarchical"),
        conditionalPanel(
          condition = "input.cluster_method == 'PCA' || input.cluster_method == 'UMAP'",
          numericInput("num_clusters", "Number of Clusters (K-Means):", value = 3, min = 2, max = 50, step = 1)
        ),
        selectInput("gene_subset", "Subset of Genes:", choices = c("All genes", "Genes with significant phenotypes (p<0.05)", "User-specific genes"), selected = "All genes"),
        conditionalPanel(
          condition = "input.gene_subset == 'User-specific genes'",
          selectizeInput("cluster_genotype", "Select Gene Symbols:", choices = NULL, multiple = TRUE, options = list(placeholder = 'Select at least 3 genes', maxOptions = 10))
        ),
        selectInput("cluster_mouse_strain", "Select Mouse Strain:", choices = NULL, selected = "All"),
        selectInput("cluster_life_stage", "Select Mouse Life Stage:", choices = NULL, selected = "All")
      )
    ),
    
    mainPanel(
      tabsetPanel(
        id = "tabs",
        
        # Tab for figure 1
        tabPanel(
          "Phenotype Scores for Knockout Mouse",
          value = "mouse_genotype_tab",
          uiOutput("no_data_message"),  # Dynamic message output
          uiOutput("plot_container"),  # Placeholder for dynamically rendered plot output
          downloadButton("download_mouse_data", "Download Mouse Data")  # Button for downloading mouse data
        ),
        
        # Tab for figure 2 
        tabPanel(
          "Knockout Gene Scores for Phenotype",
          value = "mouse_phenotype_tab",
          plotlyOutput("phenotype_mouse_plot", height = "700px", width = "100%"),  # Plot for gene scores
          downloadButton("download_phenotype_data", "Download Phenotype Data")  # Button for downloading phenotype data
        ),
        
        # Tab for figure 3
        tabPanel(
          "Gene Clusters",
          value = "clustering_tab",
          plotlyOutput("gene_cluster_plot", height = "2000px"),  # Plot for gene clusters
          downloadButton("download_cluster_data", "Download Cluster Data")  # Button for downloading cluster data
        )
      )
    )
  )
)
