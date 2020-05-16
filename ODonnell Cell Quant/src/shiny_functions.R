

get_mpi_table <- function(res)
{
  renderDT({
    if(is.null(res)) return()
    df <- as.data.frame(res$FCM[,"membrane.a.b.mean"])
   # ids <- as.numeric(names(df))
    names(df)[1] <- "mpi"

   
    datatable(df, 
              class = 'cell-border stripe',
              colnames = c( "Cell ID", "MPI"),
              options = list(
                searching = FALSE,
                lengthChange = FALSE,
                scrollCollapse = FALSE
                
              ),
              
    ) %>% 
      formatStyle(columns=colnames(df), 
                  target = 'row',
                  backgroundColor = "aquamarine4") %>% 
      formatRound(columns = colnames(df), 4)
    
  })
}

break_file_name <- function(s)
{
  sp <- strsplit(s, "_")
  sp <- unlist(sp)
  

  s1 <- paste(sp[1], sp[2], sep = "_")
  
  s2 <- sp[-1]
  s2 <- s2[-1]
  s2 <- paste(s2, collapse = "_")
  
  HTML(paste(s1, s2, sep="<br/>"))
}


add_result_ui <-function()
{
  renderUI({
  fluidRow(
    column(2, h4("Results:"),uiOutput("image_select")),
    column(10,
           tabsetPanel(
             tabPanel("Image Result", 
                      br(),
                      sidebarLayout(
                    sidebarPanel(   
                      checkboxGroupInput("toggle_selections",
                                             label = "",
                                             choices = c("Membranes" = "mem_select",
                                                         "Vacuoles" = "vac_select",
                                                         "Labels" = "label_select",
                                                         "Removed" = "rem_select")),
                      tags$hr(),
                      radioButtons("channel_sel", 
                                   label = "", 
                                   choices = c("CMAC" = "cmac",
                                               "GFP" = "gfp",
                                               "DIC" = "dic")), width = 2),
                   
                      #displayOutput("img_result")),
                    mainPanel(plotOutput("img_result_plot"),align="left", width = 4))),
             tabPanel("Summary", fluidRow(br(),
                                          column(3, uiOutput("sum_text")), 
                                          column(7, plotOutput("sum_hist")))),
             tabPanel("Table", fluidRow(
                                        br(),
                                        column(4, br(), plotOutput("labeled_img")), 
                                        column(3,DTOutput("mpi_table")))),
             tabPanel("Log", verbatimTextOutput("results_log"))),
    ))
  })
}
get_image_plot <-function(res, sel, chan)
{
 

  renderPlot(
    {
     
      if(is.null(res)) 
      {
        return()
      }
      par(mar=rep(0,4))
      if(chan == "cmac")
        plot( res$channels$cmac)
      else if(chan == "gfp")
        plot( res$channels$gfp)
      else if(chan == "dic")
        plot( res$channels$dic)
      else
        plot(res$channels$gfp)

      par(new=TRUE)
      if(!is.null(sel))
      {
       
        if("mem_select" %in% sel)
          sapply(res$mem_pts, function(x){points(x, type = 'l', col="white")})
        if("vac_select" %in% sel)
          sapply(res$vac_pts, function(x){points(x, type = 'l', col="yellow")})
        if("label_select" %in% sel)
          text(x = res$labels$label_pts[,"membrane.0.m.cx"], 
               y = res$labels$label_pts[,"membrane.0.m.cy"], 
               labels = res$labels$labels, 
               col = "red", 
               pos = c(2,3), 
               vfont = c("sans serif", "bold"))
        if("rem_select" %in% sel)
          sapply(res$removed_puts, function(x){points(x, type = 'l', col="red")})
      
      }
      par(new = TRUE)


    })

  
}

get_final_labeled <- function(res)
{
  renderPlot({
  plot(res$channels$gfp)
  sapply(res$mem_pts, function(x){points(x, type = 'l', col="white")})
  sapply(res$vac_pts, function(x){points(x, type = 'l', col="yellow")})
  text(x = res$labels$label_pts[,"membrane.0.m.cx"], 
         y = res$labels$label_pts[,"membrane.0.m.cy"], 
         labels = res$labels$labels, 
         col = "red", 
         pos = c(2,3), 
         vfont = c("sans serif", "bold"))
  })
}
get_sum_text <- function(res)
{
  renderUI({
    HTML(paste0("Number of cell membranes: ", 
           length(res$FCM[,1]), 
           br(),"Average MPI over all cells: ", 
           format(mean(res$FCM[,"membrane.a.b.mean"]), nsmall = 3),
           br(),"Number budding cells removed: "), length(res$removed_puts))
    
    
  })
}
get_hist <- function(res)
{
  renderPlot({
  df <- as.data.frame(res$FCM[,"membrane.a.b.mean"])
  names(df)[1] <- "mpi"
  ggplot(data = df, mapping = aes(x=mpi))+ 
    geom_histogram(binwidth = .05, 
                   fill="aquamarine4",
                   alpha = 0.5,
                   color = "azure4")+
    labs(title = paste("Mean Pixel Intensity: ",res$filename),
         y = "Frequency",
         x = "Mean Pixel Intensity")+
    theme(plot.title = element_text(hjust = 0.5))
  })
}


############## old functions ###############
get_image_plot_old <-function(res, sel, chan, b)
{
  
  
  renderPlot(
    {
      
      if(is.null(res)) 
      {
        return()
      }
      if(is.null(b))
      {
        b$xmin <-0
        b$xmax <- dim(res$channels$cmac)[1]
        b$ymin <-0
        b$ymax <- dim(res$channels$cmac)[2]
      }
      if(chan == "cmac")
        plot( res$channels$cmac[b$xmin:b$xmax,b$ymin:b$ymax,])
      else if(chan == "gfp")
        plot( res$channels$gfp[b$xmin:b$xmax,b$ymin:b$ymax,])
      else if(chan == "dic")
        plot( res$channels$dic[b$xmin:b$xmax,b$ymin:b$ymax])
      else
        plot(res$channels$gfp[b$xmin:b$xmax,b$ymin:b$ymax,])
      
      if(!is.null(sel))
      {
        
        if("mem_select" %in% sel)
          sapply(res$mem_pts, function(x){points(x, type = 'l', col="white")})
        if("vac_select" %in% sel)
          sapply(res$vac_pts, function(x){points(x, type = 'l', col="white")})
        if("label_select" %in% sel)
          text(x = res$labels$label_pts[,"membrane.0.m.cx"], 
               y = res$labels$label_pts[,"membrane.0.m.cy"], 
               labels = res$labels$labels, 
               col = "red", 
               pos = c(2,3), 
               vfont = c("sans serif", "bold"))
        if("rem_select" %in% sel)
          sapply(res$removed_puts, function(x){points(x, type = 'l', col="red")})
        
      }
      
      
    })
  
  
} 