setwd("C:/Users/sarah/Desktop/Capstone/Main pipeline")

source("src/functions.R")


pipeline <- function(datasetpath, testing)#, progress)
{
  if(testing)
  {
    return(readRDS("Saved Results/demo_five_results.rds"))
  }

  imageset <- read_in_imageset(datasetpath)
  results = list()
  for( row in 1:nrow(imageset))
  {

   # progress$inc(1/nrow(imageset), detail = paste0(imageset[row,"file"], "(", row, "/",nrow(imageset),")" ))
    
    
    
    channels <- read_in_channels(imageset[row,], datasetpath)
    img <- convert_to_grayscale(channels)
    img <- scale_intensities(img)
    cells <- detect_cells(img)
  
   # if(!(length(table(cells$membranes))) == length(table(cells$cell_bodies)) == length(table(cells$inner)))
    #  print("PROBLEM")
    
    
    res <- find_vacuoles(cells, img)

    
    final <- display_output(res$membranes, channels$gfp, res$vacuoles, labeled = FALSE)
    
    labels <- roi_labels_demo(res$membranes, res$FCM)
  
    writeImage(x = final,
               type = "tiff",
               bits.per.sample = 16,
               compression = "none",
               quality = 100,
               files = paste0("Output/",imageset[row, "file"], "_img.tiff"))

    
    write_to_csv(res$FCM, imageset[row, "file"])
    
  #  results[[row]] <- list(channels = channels, 
   #                        #final = final,
    #                       labels = labels,
     #                      filename = imageset[row, "file"],
                          # membranes = res$membranes,
                           #vacuoles = res$vacuoles,
      #                     mem_pts = ocontour(res$membranes),
         #                  vac_pts = ocontour(res$vacuoles),
       #                   removed_puts = ocontour(cells$removed),
        #                   FCM = res$FCM)
    gc()
  }
  message("End of main")

  #return(results)
  }


#saveRDS(res, "Saved Results/demo_dataset_three.rds")

#res <- readRDS(file = 'results-04-03-2020.rds')


