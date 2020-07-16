
source("src/functions.R")


pipeline <- function(datasetpath, testing, gui, progress)
{
  if(testing)
  {
    return(readRDS("Demo/Saved Results/presentation_results.rds"))
  }

  imageset <- read_in_imageset(datasetpath)
  results = list()
  for( row in 1:nrow(imageset))
  {
    if(gui)
      progress$inc(1/nrow(imageset), detail = paste0(imageset[row,"file"], "(", row, "/",nrow(imageset),")" ))
    
    channels <- read_in_channels(imageset[row,], datasetpath)
    img_gray <- convert_to_grayscale(channels)
    membranes <- detect_membranes(img_gray, offset, sigma)
    vacuoles <- find_vacuoles(membranes, img_gray)
    res <- exclude_and_bind(membranes, vacuoles)
    final<-tidy_up(membranes,vacuoles,res)
    
    #writeImage(fillHull(final$membranes),paste0("Masks/",imageset[row, "file"], "_pm_mask.png"), type = "png", quality = 100)
    
    
    tiff(filename = paste0("Output/",imageset[row, "file"], "_final_results.tiff"))
    
    get_display_img(df = final$df,
                    membranes = final$membranes, 
                    col_membranes = 'white', 
                    vacuoles = final$vacuoles, 
                    col_vacuoles ='yellow', 
                    removed = membranes$removed,
                    closed_vacuoles = TRUE, 
                    img = channels$gfp, 
                    showRemoved = TRUE, 
                    showMemLabels = TRUE, 
                    showVacLabels = FALSE)
        dev.off()
    
    write.csv(final$df, paste0("Output/",imageset[row, "file"], '_results.csv'), row.names=FALSE)
    results[[row]] <- list(df = final$df,
                           img_gray = img_gray,
                           channels = channels, 
                           filename = imageset[row, "file"],
                          membranes = final$membranes,
                          vacuoles = final$vacuoles,
                           mem_pts = ocontour(final$membranes),
                           vac_pts = ocontour(final$vacuoles),
                           removed_pts = ocontour(membranes$removed))
    
    
  
   
  }
  message("End of main")
 # if(gui)
    return(results)
  }


pipeline_git1 <- function(datasetpath, testing, gui, progress)
{
  if(testing)
  {
    return(readRDS("Demo/Saved Results/presentation_results.rds"))
  }
  
  imageset <- read_in_imageset_git1_tiffs(datasetpath)
  results = list()
  for( row in 1:nrow(imageset))
  {
    if(gui)
      progress$inc(1/nrow(imageset), detail = paste0(imageset[row,"filename"], "(", row, "/",nrow(imageset),")" ))
    
    channels <- read_in_channels_git1(imageset[row,], datasetpath)
    img_gray <- convert_to_grayscale_git1(channels)
    membranes <- detect_membranes_git1(img_gray, channels)
    vacuoles <- find_vacuoles(membranes, img_gray, channels)
    res <- exclude_and_bind(membranes, vacuoles)
    
   # first_pass <- get_first_pass(membranes,vacuoles,res, renum = TRUE)
    
    final<-tidy_up(membranes,vacuoles,res)
    
    if(nrow(final$df)==0)
    {
      mem_pts = NULL #list()
      vac_pts = NULL
    }
    else
    {
      mem_pts = ocontour(final$membranes)
      vac_pts = ocontour(final$vacuoles)
     
    }
    if(length(membranes$removed) > 0)
    {
      removed_pts = ocontour(membranes$removed)
    }
    else
    {
      removed_pts = NULL
    }
  
    #writeImage(fillHull(final$membranes),paste0("Masks/",imageset[row, "file"], "_pm_mask.png"), type = "png", quality = 100)
    
    
    tiff(filename = paste0("FinalOutput/Images/",imageset[row, "filename"], "_final_results.tiff"))
    
    get_display_img(df = final$df,
                    membranes = final$membranes, 
                    col_membranes = 'white', 
                    vacuoles = final$vacuoles, 
                    col_vacuoles ='yellow', 
                    removed = membranes$removed,
                    closed_vacuoles = TRUE, 
                    img = channels$gfp, 
                    showRemoved = TRUE, 
                    showMemLabels = TRUE, 
                    showVacLabels = FALSE)
    dev.off()
    
    write.csv(final$df, paste0("FinalOutput/Spreadsheets/",imageset[row, "filename"], '_results.csv'), row.names=FALSE)
    results[[row]] <- list(df = final$df,
                           img_gray = img_gray,
                           channels = channels, 
                           filename = imageset[row, "filename"],
                           membranes = final$membranes,
                           vacuoles = final$vacuoles,
                           mem_pts = mem_pts,
                           vac_pts = vac_pts,
                           removed_pts = removed_pts)
    
    
    
    
  }
  message("End of main")
  # if(gui)
  return(results)
}

pipeline_git1_interactive <- function(datasetpath, testing, gui, progress, offset, sigma)
{
  if(testing)
  {
    return(readRDS("Demo/Saved Results/presentation_results.rds"))
  }
  
  imageset <- read_in_imageset_git1_tiffs(datasetpath)
  results = list()
  for( row in 1:nrow(imageset))
  {
    if(gui)
      progress$inc(1/nrow(imageset), detail = paste0(imageset[row,"filename"], "(", row, "/",nrow(imageset),")" ))
    
    channels <- read_in_channels_git1(imageset[row,], datasetpath)
    img_gray <- convert_to_grayscale_git1(channels)
    membranes <- detect_membranes_git1(img_gray, channels, offset, sigma)
    removed <- membranes$removed
    vacuoles <- find_vacuoles(membranes, img_gray, channels)
    res <- exclude_and_bind(membranes, vacuoles)
    
    first_pass <- get_first_pass(membranes,vacuoles,res)
    #final<-tidy_up(membranes,vacuoles,res)
    
    #writeImage(fillHull(final$membranes),paste0("Masks/",imageset[row, "file"], "_pm_mask.png"), type = "png", quality = 100)
    
    
    tiff(filename = paste0("Demo/Output/",imageset[row, "filename"], "_final_results.tiff"))
    
    get_display_img(df = first_pass$df,
                    membranes = first_pass$membranes, 
                    col_membranes = 'white', 
                    vacuoles = first_pass$vacuoles, 
                    col_vacuoles ='yellow', 
                    removed = removed,
                    closed_vacuoles = TRUE, 
                    img = channel(channels$gfp, "asgreen"), 
                    showRemoved = TRUE, 
                    showMemLabels = TRUE, 
                    showVacLabels = FALSE)
    dev.off()
    
   # write.csv(final$df, paste0("Git1Output/",imageset[row, "filename"], '_results.csv'), row.names=FALSE)
    results[[row]] <- list(df = first_pass$df,
                           img_gray = img_gray,
                           channels = channels, 
                           filename = imageset[row, "filename"],
                           membranes = first_pass$membranes,
                           vacuoles = first_pass$vacuoles,
                           mem_pts = ocontour(first_pass$membranes),
                           vac_pts = ocontour(first_pass$vacuoles),
                           removed = removed,
                           removed_pts = ocontour(removed))
    
    
    
    
  }
  message("End of main")
  # if(gui)
  return(results)
}


#res <- pipeline_git1("Datasets/Git1_pngs", testing=FALSE, gui=FALSE, progress=NULL)

#saveRDS(res, "Saved Results/presentation_results.rds")


#res <- readRDS(file = 'results-04-03-2020.rds')



