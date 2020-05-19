




# Function to sort batch of imagesets into an N x 4 array where N = # of files.
# Array has structure [ [filename] [cmac image] [gfp image] [dic image] ]

### add in nd2 functionality as well as extention parameters
read_in_imageset <- function(datapath)
{
  
  dirs <- list.dirs(path = datapath)
  dirs <- dirs[-1]
  
  imagesets <- matrix(NA, length(dirs), 4)
  colnames(imagesets) <- c("file", "cmac", "gfp", "dic")
  
  for(dir in seq_along(dirs))
  {
    imagesets[dir,"file"] = str_remove(dirs[dir], paste0(datapath, "/"))
    
    
    imagesets[dir, "cmac"] <- list.files(
      path = file.path(dirs[dir]),
      pattern='cmac.*?\\.tif',
      ignore.case = TRUE)
    
    imagesets[dir, "gfp"] <- list.files(
      path = file.path(dirs[dir]),
      pattern='gfp.*?\\.tif',
      ignore.case = TRUE)
    
    imagesets[dir, "dic"] <- list.files(
      path = file.path(dirs[dir]),
      pattern='dic.*?\\.tif',
      ignore.case = TRUE)
    
  }
  return(imagesets)
}




# Function to read in channels for image.  Returns list of 3 channels
read_in_channels <- function(imageset, datasetpath)
{
  message("#####################################################")
  message(paste0("Examining image: ", imageset["file"]))
  dic_path=paste0(datasetpath,"/", imageset["file"], "/", imageset["dic"])
  gfp_path=paste0(datasetpath, "/", imageset["file"], "/", imageset["gfp"])
  cmac_path=paste0(datasetpath, "/", imageset["file"], "/", imageset["cmac"])
  
  
  
  
  dic = readImage(file.path(dic_path))
  gfp = readImage(file.path(gfp_path))
  cmac = readImage(file.path(cmac_path))
  
  list(cmac = cmac, gfp = gfp, dic = dic)
  
}

write_to_csv <-function(FCM, filename)
{
  
  mean_pixel_intensity <- FCM[,"membrane.a.b.mean"]
  ids <- as.numeric(names(mean_pixel_intensity))
  
  
  
  mean_pixel_intensity <- as.data.frame( mean_pixel_intensity) 
  mean_pixel_intensity$cell_id = ids
  mean_pixel_intensity = rev(mean_pixel_intensity)
  
  
  write.table(mean_pixel_intensity, 
              file = paste0("Output/",filename, "_quant.csv"),
              sep=",", row.names = FALSE, quote=FALSE)
  
}
