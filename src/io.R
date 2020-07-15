




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

read_in_imageset_git1 <- function(datapath)
{
  
  dirs <- list.dirs(path = datapath)
  dirs <- dirs[-1]
  
  imagesets <- matrix(NA, length(dirs), 3)
  colnames(imagesets) <- c("file", "cmac", "gfp")
  
  for(dir in seq_along(dirs))
  {
    imagesets[dir,"file"] = str_remove(dirs[dir], paste0(datapath, "/"))
    
    
    imagesets[dir, "cmac"] <- list.files(
      path = file.path(dirs[dir]),
      pattern='.*?cmac\\.png',
      ignore.case = TRUE)
    
    imagesets[dir, "gfp"] <- list.files(
      path = file.path(dirs[dir]),
      pattern='.*?gfp\\.png',
      ignore.case = TRUE)
  
    
  }
  return(imagesets)
}

read_in_imageset_git1_tiffs <- function(datapath)
{
  
  files <- list.files(path = datapath)

  imagesets <- matrix(NA, length(files), 2)
  colnames(imagesets) <- c("filepath", "filename")
  
  for(file in seq_along(files))
  {
    imagesets[file,"filename"] = str_remove(files[file], ".tif")
    
    imagesets[file,"filepath"]= files[file]
    
    
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
  
  
  
  
  dic = readImage(file.path(dic_path), as.is = TRUE)
  gfp = readImage(file.path(gfp_path),  as.is = TRUE)
  cmac = readImage(file.path(cmac_path),  as.is = TRUE)
  
  list(cmac = cmac, gfp = gfp, dic = dic)
  
}


read_in_channels_git1 <- function(imageset, datasetpath)
{
  message("#####################################################")
  message(paste0("Examining image: ", imageset["filename"]))
  gfp_path=paste0(datasetpath, "/", imageset["file"], "/", imageset["gfp"])
  cmac_path=paste0(datasetpath, "/", imageset["file"], "/", imageset["cmac"])
  
  
  gfp = readImage(file.path(paste0(datasetpath, "/", imageset["filepath"])))[,,2]
  cmac = readImage(file.path(paste0(datasetpath, "/", imageset["filepath"])))[,,1]
  
  ref_gfp = readImage(file.path(paste0(datasetpath, "/", imageset["filepath"])),  as.is = TRUE)[,,2]
  ref_cmac = readImage(file.path(paste0(datasetpath, "/", imageset["filepath"])),  as.is = TRUE)[,,1]
  
  list(cmac = cmac, gfp = gfp, ref_cmac = ref_cmac, ref_gfp = ref_gfp)
  
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
