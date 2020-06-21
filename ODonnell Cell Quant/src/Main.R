
source("src/functions.R")


pipeline <- function(datasetpath, testing, gui, progress)
{
  if(testing)
  {
    return(readRDS("Saved Results/presentation_results.rds"))
  }

  imageset <- read_in_imageset(datasetpath)
  results = list()
  for( row in 1:nrow(imageset))
  {
    if(gui)
      progress$inc(1/nrow(imageset), detail = paste0(imageset[row,"file"], "(", row, "/",nrow(imageset),")" ))
    
    channels <- read_in_channels(imageset[row,], datasetpath)
    img_gray <- convert_to_grayscale(channels)
    membranes <- detect_membranes(img_gray)
    vacuoles <- find_vacuoles(membranes, img_gray)
    res <- exclude_and_bind(membranes, vacuoles)
    final<-tidy_up(membranes,vacuoles,res)
    
    writeImage(fillHull(final$membranes),paste0("Masks/",imageset[row, "file"], "_pm_mask.png"), type = "png", quality = 100)
    
    
   # tiff(filename = paste0("Output/",imageset[row, "file"], "_final_results.tiff"))
    
  #  get_display_img(df = final$df,
   #                 membranes = final$membranes, 
    #                col_membranes = 'white', 
     #               vacuoles = final$vacuoles, 
      #              col_vacuoles ='yellow', 
       #             removed = membranes$removed,
        #            closed_vacuoles = TRUE, 
         #           img = channels$gfp, 
          #          showRemoved = TRUE, 
           #         showMemLabels = TRUE, 
            #        showVacLabels = FALSE)
#    dev.off()
    
 #   write.csv(final$df, paste0("Output/",imageset[row, "file"], '_results.csv'), row.names=FALSE)
  #  results[[row]] <- list(df = final$df,
   #                        img_gray = img_gray,
    #                       channels = channels, 
     #                      filename = imageset[row, "file"],
      #                    membranes = final$membranes,
       #                   vacuoles = final$vacuoles,
        #                   mem_pts = ocontour(final$membranes),
         #                  vac_pts = ocontour(final$vacuoles),
          #                 removed_pts = ocontour(membranes$removed))
    
    
  
   
  }
  message("End of main")
 # if(gui)
    return(results)
  }

#res <- pipeline("Datasets/Diploid/Diploid_Dataset_all/FirstThree", testing=FALSE, gui=FALSE, progress=NULL)

#saveRDS(res, "Saved Results/presentation_results.rds")


#res <- readRDS(file = 'results-04-03-2020.rds')



