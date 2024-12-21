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
install.packages("RSQLite")
library(DBI)
library(RSQLite)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Mouse Phenotype Visualisation Dashboard"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          selectInput("mouse_id", "Select Knockout Mouse:",
                      choices = NULL,  # Populate dynamically in server
                      selected = NULL),
          selectInput("phenotype", "Select Phenotype:",
                      choices = NULL,  # Populate dynamically in server
                      selected = NULL),
          sliderInput("threshold", "Set Significant Threshold:",
                      min = 0, max = 100, value = 50)
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel( #this neatly organises the 3 visualisations, making the app intuititive to use
            tabPanel("Phenotype Scores for Knockout Mouse",
                     plotOutput("mouse_phenotype_plot"),
                     downloadButton("download_mouse_data", "Download Knockout Mouse Data")),
            tabPanel("Knockout Mice for Phenotype",
                     plotOutput("phenotype_mouse_plot"),
                     downloadButton("download_phenotype_data", "Download Phenotype Data")),
            tabPanel("Gene Clusters",
                     plotOutput("gene_cluster_plot"),
                     downloadButton("download_cluster_data", "Download Cluster Data"))
          )
        )
    )
)

# Define server logic required to draw a histogram
# Server
server <- function(input, output, session) {
  
  # Connect to the database
  con <- dbConnect(SQLite(), "path_to_your_database.sqlite")
  
  # Ensure connection is closed when the app stops
  onStop(function() {
    dbDisconnect(con)
  })
  
  # Populate dropdown with available knockout mice
  observe({
    gene_symbols <- dbGetQuery(con, "SELECT DISTINCT gene_symbol FROM Genes")
    updateSelectInput(session, "mouse_id", choices = gene_symbols$gene_symbol)
  })
  
  # Visualisation test 1: Phenotype Scores for Knockout Mouse
  # Query and plot data for the selected mouse
  output$mouse_phenotype_plot <- renderPlot({
    req(input$mouse_id)  # Ensure a mouse ID is selected
    
    query <- sprintf("
      SELECT Analyses.p_value, Parameters.parameter_name 
      FROM Analyses 
      JOIN Parameters ON Analyses.parameter_id = Parameters.parameter_id 
      WHERE Analyses.gene_id IN (
        SELECT gene_id 
        FROM Genes 
        WHERE gene_symbol = '%s'
      ) AND Analyses.p_value <= %f
      ORDER BY Analyses.p_value ASC;", input$mouse_id, input$threshold)
    
    data <- dbGetQuery(con, query)
    
    ggplot(data, aes(x = reorder(parameter_name, p_value), y = p_value)) +
      geom_bar(stat = "identity", fill = "blue") +
      labs(title = paste("Phenotype Scores for", input$mouse_id),
           x = "Phenotypes", y = "p-value") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  # Populate dropdown with available phenotypes
  observe({
    phenotypes <- dbGetQuery(con, "SELECT DISTINCT parameter_name FROM Parameters")
    updateSelectInput(session, "phenotype", choices = phenotypes$parameter_name)
  })
  
  #Download Test 1
  output$download_mouse_data <- downloadHandler(
    filename = function() {
      paste("knockout_mouse_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(input$mouse_id)  # Ensure a mouse ID is selected
      query <- sprintf("
      SELECT Analyses.p_value, Parameters.parameter_name 
      FROM Analyses 
      JOIN Parameters ON Analyses.parameter_id = Parameters.parameter_id 
      WHERE Analyses.gene_id IN (
        SELECT gene_id 
        FROM Genes 
        WHERE gene_symbol = '%s'
      )
      ORDER BY Analyses.p_value ASC;", input$mouse_id)
      data <- dbGetQuery(con, query)
      write.csv(data, file, row.names = FALSE)
    }
  )
  
  # Visualisation 2: Scores of All Knockout Mice for a selected Phenotype
  output$phenotype_mouse_plot <- renderPlot({
    req(input$phenotype)  # Ensure a phenotype is selected
    
    query <- sprintf("
      SELECT Analyses.p_value, Genes.gene_symbol 
      FROM Analyses 
      JOIN Parameters ON Analyses.parameter_id = Parameters.parameter_id 
      JOIN Genes ON Analyses.gene_id = Genes.gene_id 
      WHERE Parameters.parameter_name = '%s'
      ORDER BY Analyses.p_value ASC;", input$phenotype)
    
    data <- dbGetQuery(con, query)
    
    ggplot(data, aes(x = reorder(gene_symbol, p_value), y = p_value)) +
      geom_boxplot() +
      labs(title = paste("Scores of All Knockout Mice for Phenotype:", input$phenotype),
           x = "Knockout Mice", y = "p-value") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  #Download Test 2
  output$download_phenotype_data <- downloadHandler(
    filename = function() {
      paste("phenotype_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      req(input$phenotype)  # Ensure a phenotype is selected
      query <- sprintf("
      SELECT Analyses.p_value, Genes.gene_symbol 
      FROM Analyses 
      JOIN Parameters ON Analyses.parameter_id = Parameters.parameter_id 
      JOIN Genes ON Analyses.gene_id = Genes.gene_id 
      WHERE Parameters.parameter_name = '%s'
      ORDER BY Analyses.p_value ASC;", input$phenotype)
      data <- dbGetQuery(con, query)
      write.csv(data, file, row.names = FALSE)
    }
  )
  
  #Visualisation Test 3: Gene clusters
  #Execute the query to retrieve the full dataset and perform clustering
  output$gene_cluster_plot <- renderPlot({
    query <- "
      SELECT Genes.gene_symbol, Parameters.parameter_name, Analyses.p_value 
      FROM Analyses 
      JOIN Parameters ON Analyses.parameter_id = Parameters.parameter_id 
      JOIN Genes ON Analyses.gene_id = Genes.gene_id;"
    
    data <- dbGetQuery(con, query)
    
    # Pivot the data for clustering
    cluster_data <- dcast(data, gene_symbol ~ parameter_name, value.var = "p_value", fill = 0)
    
    # Perform clustering
    dist_matrix <- dist(cluster_data[, -1])
    hc <- hclust(dist_matrix, method = "ward.D")
    
    # Plot dendrogram
    plot(hc, main = "Cluster of Genes with Similar Phenotype Scores",
         xlab = "Genes", ylab = "Distance", sub = "")
  })
  
  #Download Test 3
  output$download_cluster_data <- downloadHandler(
    filename = function() {
      paste("gene_cluster_data-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      query <- "
      SELECT Genes.gene_symbol, Parameters.parameter_name, Analyses.p_value 
      FROM Analyses 
      JOIN Parameters ON Analyses.parameter_id = Parameters.parameter_id 
      JOIN Genes ON Analyses.gene_id = Genes.gene_id;"
      data <- dbGetQuery(con, query)
      write.csv(data, file, row.names = FALSE)
    }
  )
  
}
  


# Run the application 
shinyApp(ui = ui, server = server)
