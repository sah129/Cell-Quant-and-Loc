library(shiny)
library(shinyjs)
library(shinythemes)
library(shinyFiles)
library(EBImage)
library(ggplot2)
library(DT)
library(zip)

source("src/Main.R")
source("src/shiny_functions.R")
source("src/appLite-functions.R")


ui <- fluidPage(
  shinyjs::useShinyjs(),
  
  titlePanel( h1( "CellQuant: O'Donnell Lab", align = "center") , windowTitle = "O'Donnell Lab"),
  fluidRow(
    column(11,  shinyDirButton('datasetpath', 'Select a directory containing the images for the pipeline', 'Please select a folder', FALSE, style = "width: 95%"), align = "center"),
    
    column(1, actionButton("process_dataset", "Run"))),
  tags$hr(),
  fluidRow(column(10, verbatimTextOutput("stream_main"))),
  uiOutput("results"),
  tags$hr(),
  fluidRow(column(12, h6("Created by Sarah Hawbaker, O'Donnell Lab, University of Pittsburgh, 2020.", align = "center")))
)
server <- function(input, output,session) 
{
  v <- reactiveValues(res = NULL, log = "")
  volumes <- c(def_root = "C:/Users/Sarah's Computer/Documents/GitHub/Cell-Quant-and-Loc/Datasets/", 
               def_dataset = "C:/Users/Sarah's Computer/Documents/GitHub/Cell-Quant-and-Loc/Datasets/Dummy 1B")
  
  shinyDirChoose(input, 'datasetpath',
                 roots = volumes,
                 filetypes=c(',', 'txt'),
                 session = session)

  dpath <- reactive(return( parseDirPath(volumes["def_root"], input$datasetpath)))  
  

  observeEvent(input$process_dataset, {
     progress <- shiny::Progress$new()
     on.exit(progress$close())
     withCallingHandlers(
       {
         shinyjs::html("stream_main", "")
         progress$set(message = "Running pipeline...", value = 0)
         v$res <- pipeline(dpath(),testing=FALSE, gui=TRUE, progress=progress, interactive = FALSE) 
         output$stream_main <- renderText("")
         output$results <- finish_screen()
       },
       message = function(m) {
         shinyjs::html(id = "stream_main", html = m$message, add = TRUE) #save the logfile somewhere
         v$log <- paste0(v$log, m$message)
       })
     
   }, once = TRUE)
  
output$results_log <- renderText(v$log)

 
output$downloadData <- downloadHandler(
    filename = function() 
      {
      paste("resultsDownload", "zip", sep=".")
    },
    content = function(file) 
      {
      zip(file, c("FinalOutput/Images", "FinalOutput/Individual Spreadsheets", "FinalOutput/Aggregated Spreadsheets"))
    },
    contentType = "application/zip"
  )



  
  

}

shinyApp(ui, server)