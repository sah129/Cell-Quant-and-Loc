
library(shiny)
library(dplyr)
library(stringr)
library(gridExtra)
library(shinyFiles)
source('src/functions.R')
source('src/Main.R')

options(shiny.maxRequestSize = 100*1024^2)

ui <- fluidPage(
  
  titlePanel("Automated Cell Quant - O'Donnell Lab", windowTitle = "Automated Cell Quant - O'Donnell Lab"),
  sidebarLayout(
    sidebarPanel(
      h4("Pipeline Options", align = "center"),
      fluidRow(
        column(1, strong("1.")),
        
        column(3, shinyDirButton('input_dir', 'Input Folder', 'Please select a folder')),
        column(6, textOutput("input_dir_text")),  
        
        column(2, actionButton("inputdir_help", "?"))),
      fluidRow(
      column(1, strong("1.")),
      
      column(3, shinyDirButton('output_dir', 'Output Folder', 'Please select a folder')),
      column(6, textOutput("output_dir_text")), 
      column(2, actionButton("outputdir_help", "?"))),
    
      fluidRow(
        column(1, strong("2.")),
        column(3, textInput("cmac_chan", "CMAC", value = "", placeholder = "1")),
        column(3, textInput("gfp_chan", "GFP", value = "", placeholder = "2")),
        column(3, textInput("dic_chan", "DIC", value = "", placeholder = "3")),
        column(2, br(), actionButton("inputchannels_help", "?"))),
      fluidRow(
        column(1, strong("3.")),
        column(9, textInput("cutoff_value", "Cell size cutoff", placeholder="100")),
        column(2, br(), actionButton("cutoffvalue_help", "?"))),
      
  
      
      fluidRow(
        column(1, strong("4.")),
        column(9, radioButtons("algchoose", "Membrane Detection Algorithm", choices = c("GFP", "DIC"), selected = "GFP", inline = TRUE)),
        column(2, br(), actionButton("algchoose_help", "?"))),
      fluidRow(
        column(1, strong("5.")),
        
        column(9, radioButtons("factorchoose", "Factor:", choices = c("1", "2", "4", "8", "16"), selected = "5", inline = TRUE)),
        
        
        column(2, br(), actionButton("factorchoose_help", "?"))),
      fluidRow(actionButton("run", "Run Pipeline"), align = "center"),
      width=4),
    mainPanel(
   
      
      width = 8)
    
    
    
    
  )
  
  )

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  values <- reactiveValues()
  volumes = c(home = getwd())
  

  shinyDirChoose(input, 'input_dir', roots= volumes, session = session)

  dpath <- reactive(return(parseDirPath(volumes, input$input_dir)))

  output$input_dir_text <- renderText({dpath()})
  
  shinyDirChoose(input, 'output_dir', roots= volumes, session = session)
  
  opath <- reactive(return(parseDirPath(volumes, input$output_dir)))
  
  output$output_dir_text <- renderText({opath()})
  
  
  # help buttons
  observeEvent(input$inputdir_help, {
    showModal(modalDialog(
      title = "modal title input",
      "this is the help text for input",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  observeEvent(input$outputdir_help, {
    showModal(modalDialog(
      title = "modal title output",
      "this is the help text for output",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  observeEvent(input$inputchannels_help, {
    showModal(modalDialog(
      title = "modal title channel",
      "this is the help text for channel",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  observeEvent(input$cutoffvalue_help, {
    showModal(modalDialog(
      title = "modal title cutoff",
      "this is the help text for cutoff",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  observeEvent(input$algchoose_help, {
    showModal(modalDialog(
      title = "modal title algchoose",
      "this is the help text for algchoose",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  observeEvent(input$factorchoose_help, {
    showModal(modalDialog(
      title = "modal title factorchoose",
      "this is the help text for factorchoose",
      easyClose = TRUE,
      footer = NULL
    ))
  })
  

  
  
  
  observeEvent(input$run,
     {

       if (is.null(dpath()) || is.null(input$factorchoose) || is.null(input$gfp_chan) || is.null(input$cmac_chan) || is.null(input$dic_chan) || is.null(input$algchoose) || is.null(input$cutoff_value))
       {
         print('null')
         return(NULL)
       }
       
       progress <- shiny::Progress$new()
       
       on.exit(progress$close())
       
       progress$set(message = "Running Test Cases...", value = 0)
       
       values$res <- pipeline(dpath(),
                              gui=TRUE, 
                              progress=progress, 
                              factor = input$factorchoose, 
                              gfp_chan = input$gfp_chan,
                              dic_chan = input$dic_chan,
                              cmac_chan = input$cmac_chan,
                              alg = input$algchoose,
                              cutoff = input$cutoff_value,
                              outpath = opath()) 
       
       showModal(modalDialog(
         title = "Pipeline Complete!",
         "Some stuff",
         easyClose = TRUE,
         footer = NULL
       ))
       
     })

  # output$contents <- renderUI(
  
  #  if(is.null(values$res))
  #       return(NULL)
  #else
  #  return(uiOutput(textOutput("DONE")))
  #  )
  
  
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
