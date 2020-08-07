
exclude_and_bind_puncta <- function(cells, fc, gfp)
{
  df <- data.frame(matrix(NA, nrow = length(table(cells)), ncol = 4))
  names(df) <- c('CellID', 'puncta_count', 'pm_center_x', 'pm_center_y')
  
  l = length(table(cells))-1
  empty_cells <- vector("list", l) # cells containing no detected vacuoles
  
  pl = Image(data=0, dim = dim(gfp))
  # Loop through list of PMs
  for(i in seq(1:l))
  {
    cell <- cells == i
    puncta <- cell * gfp
    puncta <- normalize(puncta)
    
    
    pos <- which(puncta > 0)
    pixels <- puncta[pos]
    
    
    puncta <- puncta > quantile(pixels, 0.98)
   # puncta <- puncta > 0.8
    
    puncta <- bwlabel(puncta)
    pl = pl + puncta
    puncta_count <- length(table(puncta)[-1])
    if(puncta_count == 0)
      empty_cells[i] = i
    else
    {
      df[i,] <- list(CellID = i,
                     puncta_count = puncta_count,
                     pm_center_x = fc[i,'cell.0.m.cx'],
                     pm_center_y = fc[i,'cell.0.m.cy'])
      
    }
  }
  list(df=df, empty_cells = unlist(empty_cells), pl = pl)
}

read_in_channels_puncta <- function(imageset, datasetpath)
{
  message("#####################################################")
  message(paste0("Examining image: ", imageset["filename"]))
  gfp_path=paste0(datasetpath, "/", imageset["file"], "/", imageset["gfp"])
  cmac_path=paste0(datasetpath, "/", imageset["file"], "/", imageset["cmac"])
  
  gfp = readImage(file.path(paste0(datasetpath, "/", imageset["filepath"])))[,,2]
  cmac = readImage(file.path(paste0(datasetpath, "/", imageset["filepath"])))[,,1]
  ip = readImage(file.path(paste0(datasetpath, "/", imageset["filepath"])))[,,3]
  
  ref_gfp = readImage(file.path(paste0(datasetpath, "/", imageset["filepath"])),  as.is = TRUE)[,,2]
  ref_cmac = readImage(file.path(paste0(datasetpath, "/", imageset["filepath"])),  as.is = TRUE)[,,1]
  ref_ip = readImage(file.path(paste0(datasetpath, "/", imageset["filepath"])),  as.is = TRUE)[,,3]
  
  list(cmac = cmac, gfp = gfp, ref_cmac = ref_cmac, ref_gfp = ref_gfp, ip = ip, ref_ip = ref_ip)
}

get_boundaries <- function(channels)
{
  n<-normalize(channels$cmac)
  #g<-gblur(n,sigma=5)
  #gt<-thresh(g)
  
  gt<-n > otsu(n)

  
  
  cells <- bwlabel(gt)
  fc<-computeFeatures.shape(cells)
  sel <- which(fc[,"s.area"]<100)
  cells <- rmObjects(cells, sel)
  
  l = length(table(cells))-1
  empty_cells <- vector("list", l) # cells containing no detected vacuoles
  
  
  for(i in seq(1:l))
  {
    cell <- cells == i
    if(all(fillHull(cell) == cell))
    {
      empty_cells[i] = i
    }
  }
  empty_cells <- unlist(empty_cells)
  cells <- rmObjects(cells, empty_cells)
  cells<-fillHull(cells)
  cells<-dilate(cells, kern = makeBrush(5, shape = 'disc'))
  return(cells)
}


read_in_imageset_files <- function(datapath)
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