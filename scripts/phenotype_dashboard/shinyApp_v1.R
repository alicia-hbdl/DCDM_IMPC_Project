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

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Phenotype Scores for Knockout Mouse"),

    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
          selectInput("mouse_id", "Select Knockout Mouse:",
                      choices = NULL,  # Populated in server
                      selected = NULL),
          selectInput("phenotype", "Select Phenotype:",
                      choices = NULL,  # Populated in server
                      selected = NULL),
          sliderInput("threshold", "Set Significant Threshold:",
                      min = 0, max = 1, value = 0.05, step = 0.01)
        ),

        # Show a plot of the generated distribution
        mainPanel(
          tabsetPanel( #this neatly organises the 3 visualisations, making the app intuitive to use
            # 1) Phenotype Scores for a Single Knockout Mouse
            tabPanel("Phenotype Scores for Knockout Mouse",
                     plotOutput("mouse_phenotype_plot"),
                     downloadButton("download_mouse_data", "Download Mouse Data")),
            
            # 2) Scores of All Knockout Mice for a Selected Phenotype
            tabPanel("Knockout Mice for Phenotype",
                     plotOutput("phenotype_mouse_plot"),
                     downloadButton("download_phenotype_data", "Download Phenotype Data")),
            
            # 3) Gene Clusters
            tabPanel("Gene Clusters",
                     plotOutput("gene_cluster_plot"),
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
    dbname   = "IMPCDb",
    host     = "localhost",
    port     = 3306,
    user     = "root",
    password = "password"
  )
  
  # Ensure connection is closed when app stops
  onStop(function() {
    dbDisconnect(con)
  })
  
  # Populate knockout mice dropdown
  observe({
    gene_symbols <- dbGetQuery(con, "SELECT DISTINCT gene_symbol FROM Genes;")
    updateSelectInput(session, "mouse_id", choices = gene_symbols$gene_symbol)
  })
  
  # Populate phenotypes dropdown
  observe({
    phenotypes <- dbGetQuery(con, "SELECT DISTINCT parameter_name FROM Parameters;")
    updateSelectInput(session, "phenotype", choices = phenotypes$parameter_name)
  })
  
  # Visualisation 1: Phenotype Scores for a Single Knockout Mouse
  output$mouse_phenotype_plot <- renderPlot({
    req(input$mouse_id)  # Ensure a mouse ID is selected
    
    # Use the threshold from the slider as a float in the SQL query
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
    
    # If no data returned, display a message or skip plotting
    if (nrow(data) == 0) {
      plot.new()
      title("No phenotypes meet the selected threshold for this knockout mouse.")
      return()
    }
    
    ggplot(data, aes(x = reorder(parameter_name, p_value), y = p_value)) +
      geom_bar(stat = "identity", fill = "blue") +
      labs(title = paste("Phenotype Scores for", input$mouse_id),
           x = "Phenotypes", y = "p-value") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  #Download data for mouse phenotype plot
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
  
  # Visualisation 2: Scores of All Knockout Mice for a Selected Phenotype
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
    
    if (nrow(data) == 0) {
      plot.new()
      title("No data available for this phenotype.")
      return()
    }
    
    ggplot(data, aes(x = reorder(gene_symbol, p_value), y = p_value)) +
      geom_boxplot(fill = "blue", alpha = 0.6) +
      labs(title = paste("Scores of All Knockout Mice for Phenotype:", input$phenotype),
           x = "Knockout Mice", y = "p-value") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  })
  
  #Download data for phenotype mouse plot
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
  
  #Visualisation 3: Gene clusters
  output$gene_cluster_plot <- renderPlot({
    query <- "
      SELECT Genes.gene_symbol, Parameters.parameter_name, Analyses.p_value 
      FROM Analyses 
      JOIN Parameters ON Analyses.parameter_id = Parameters.parameter_id 
      JOIN Genes ON Analyses.gene_id = Genes.gene_id;"
    
    data <- dbGetQuery(con, query)
    
    if (nrow(data) == 0) {
      plot.new()
      title("No data available for clustering.")
      return()
    }
    
    # Pivot the data for clustering (reshape to wide format)
    cluster_data <- dcast(data, gene_symbol ~ parameter_name, value.var = "p_value", fill = 0)
    
    # Perform clustering
    # Exclude the gene_symbol column from the distance matrix
    dist_matrix <- dist(cluster_data[, -1])
    hc <- hclust(dist_matrix, method = "ward.D")
    
    # Plot dendrogram
    plot(hc, main = "Cluster of Genes with Similar Phenotype Scores",
         xlab = "Genes", ylab = "Distance", sub = "")
  })
  
  #Download data for gene clusters
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
