#functionsfor pipeline 

library("EBImage")
library("genefilter")
library("stringr")
# defaults for image frames
cmac_channel = 1
gfp_channel = 2
dic_channel = 3


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


####I/O######


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


# Converts all channels to grayscale
convert_to_grayscale<- function(img)
{
  colorMode(img$gfp) <- 'Grayscale'
  colorMode(img$cmac) <- 'Grayscale'
  
  gfp_gray <- getFrame(img$gfp, 2, type = 'render')
  cmac_gray <- getFrame(img$cmac, 3, type = 'render')
  
  
  #gfp_gray = channel(img$gfp, 'gray')
  #cmac_gray = channel(img$cmac, 'gray')
  
  # RBioformats has .nd2 capability, add in later
  imgn = combine(cmac_gray, gfp_gray, img$dic)
  return(imgn)
}

# normalizes CMAC and GFP channels
scale_intensities <- function(img)
{
  img[,,cmac_channel] = normalize(img[,,cmac_channel])
  img[,,gfp_channel] = normalize(img[,,gfp_channel])
  
  return(img)
}


######CELLS############


detect_cells <-function(img)
{
  message("########################CELLS########################")
  
  cells_blurred = gblur(img[,,gfp_channel], sigma = 2)
  ct = thresh(cells_blurred, offset = 0.005)
  
  cm = bwlabel(ct)
  FS = computeFeatures.shape(cm)
  sel <- which(FS[,"s.area"] < 500)
  membranes <-rmObjects(cm, sel)
  
  
  cell_bodies <- fillHull(membranes)
  cell_bodies <- bwlabel(cell_bodies)
  

  FCCB <- computeFeatures.moment(cell_bodies)
  res <- remove_edge_membranes(membranes, FCCB, dim(img[,,cmac_channel]))

  cell_bodies <- res$cell_bodies
  
  
  cell_inner <- fillHull(cell_bodies) - cell_bodies
  cell_inner <- bwlabel(cell_inner)

  
  
  membranes <- cell_bodies - cell_inner
  
  FCCI <- computeFeatures.shape(cell_inner)
  message(paste0("Number of cells detected on first pass: ", length(table(cell_bodies))))
  message(paste0("Number of inner detected on first pass: ", length(table(cell_inner))))
  
  list(
    cell_inner = cell_inner,
    removed = res$removed,
    FCCI = FCCI, 
    membranes = membranes, 
    #FCM = FCM, 
    cell_bodies = cell_bodies, 
    FCCB = res$FCCB)
  
}


remove_edge_membranes <-function(cell_bodies, FMM, dim_img)
{
  
 
  
  contours <- ocontour(cell_bodies)
  bound <- list(l = 3,
                r = dim_img[1]-3,
                t = 3,
                b = dim_img[2]-3)
  # find more efficient way to do this
  left <- lapply(contours, function(x){min(x[,1])})
  right <- lapply(contours, function(x){max(x[,1])})
  top <- lapply(contours, function(x){min(x[,2])})
  bottom <- lapply(contours, function(x){max(x[,2])})
  edge_cells <- c(which(left < bound$l),which(right > bound$r),  which(top < bound$t),which(bottom > bound$b))
  edge_cells <- as.numeric(names(edge_cells))
  edge_cells <- unique(edge_cells) #line might be unecessary

  removed_ind = which(!(seq(1, length(table(cell_bodies))) %in% edge_cells))

  removed = rmObjects(cell_bodies, removed_ind)
  
  cell_bodies <- rmObjects(cell_bodies, edge_cells)
  cell_bodies <- bwlabel(cell_bodies)
  FCCB <- computeFeatures.shape(cell_bodies)#(cell_bodies, ref = img[,,gfp_channel], xname = "seg")
  
  
  list(removed = removed, cell_bodies = cell_bodies, FCCB = FCCB)
  
}


#####VACUOLES###########

# Finds vacuoles from CMAC image. This is done by creating a moving window with w,h dimensions based on the biggest minimum cell radius in the image.
#  This allows us to detect vacuoles within the cell without picking up noise in the cytoplasm.  We then mask the vacuoles using the cells,
# in order to exclude vacuoles that reside in cells we have already excluded from our working set. 
find_vacuoles <- function(cell_info, img)
{
  message("######################VACUOLES#######################")
  
  b <- max(cell_info$FCCI[,"s.radius.min"])
 
  vmask = thresh(img[,,cmac_channel],w=b,h=b,offset = sd(img[,,cmac_channel]))
  vmask = opening(vmask, makeBrush(5,shape="disc"))
  vmask = fillHull(vmask)
  vmask = bwlabel(vmask)

  message(paste0("Number of vacuoles detected on first pass: ", format(length(table(vmask)), nsmall = 4)))
  
  vmask <- vmask*(cell_info$cell_inner) #to ensure no intersection
  vmask <- bwlabel(vmask)
  
  message(paste0("Number of vacuoles after segmask: ", length(table(vmask))))
  
  FCM <- computeFeatures(cell_info$membranes, ref = img[,,gfp_channel], xname = "membrane")
  FV <- computeFeatures(vmask, ref= img[,,gfp_channel], xname = "vac")

  res<-list(vacuoles = vmask, 
       membranes = cell_info$membranes,
       FCM = FCM,
         FV = FV,
      # FCCB = FCCB,
       cell_inner = cell_info$cell_inner,
       cell_bodies = cell_info$cell_bodies,
       FCCI = cell_info$FCCI)
  
 # res <- remove_empty_cells(cell_info, vmask, img)
  #message(paste0("Number of vacuoles after multiples removed: ", length(table(res$vacuoles))))
  return(res)

}


detect_cells_old <-function(img)
{
  message("########################CELLS########################")
  
  cells_blurred = gblur(img[,,gfp_channel], sigma = 2)
  ct = thresh(cells_blurred, offset = 0.005)
  
  cm = bwlabel(ct)
  FS = computeFeatures.shape(cm)
  sel <- which(FS[,"s.area"] < 500)
  membranes <-rmObjects(cm, sel)
  
  # membranes <- bwlabel(membranes)
  
  
  cell_bodies <- fillHull(membranes)
  cell_bodies <- bwlabel(cell_bodies)
  
  
  FCCB <- computeFeatures.shape(cell_bodies)
  
  
  
  
  cell_inner <- fillHull(membranes) - membranes
  cell_inner <- bwlabel(cell_inner)
  cell_inner <- fillHull(cell_inner)
  
  cell_inner <- cell_inner*cell_bodies
  cell_inner <- bwlabel(cell_inner)
  
  
  FCI <- computeFeatures.shape(cell_inner)
  sel <- which(FCI[,"s.area"] < 500)
  cell_inner <-rmObjects(cell_inner, sel)
  
  #  cell_inner <- bwlabel(cell_inner)
  
  
  membranes <- cell_bodies - cell_inner
  # membranes <- bwlabel(membranes)
  
  #FCM <- computeFeatures(membranes, ref = img[,,gfp_channel], xname = "membrane")
  
  
  
  FCCI <- computeFeatures.shape(cell_inner)
  message(paste0("Number of cells detected on first pass: ", length(table(cell_bodies))))
  message(paste0("Number of inner detected on first pass: ", length(table(cell_inner))))
  
  list(
    cell_inner = cell_inner,
   # removed = mask$removed,
    FCCI = FCCI, 
    membranes = membranes, 
    #FCM = FCM, 
    cell_bodies = cell_bodies, 
    FCCB = FCCB)
  
}

remove_empty_cells <- function(cell_info, vmask, img)
{
  
  l = length(table(cell_info$cell_bodies))
  l = l-1
  
  empty_cells <- vector("list", l)
  
  for(i in seq(1:l))
  {
    
    cseg = cell_info$cell_bodies == i
    vcount <- table(cseg*vmask)
    vcount <- vcount[-1]
    #message(vcount)
    if( length(vcount) == 0)
    {
      empty_cells[i] = i
    }
  }
  
  empty_cells <- unlist(empty_cells)
  
  
  message(paste0("Number of empty cells: ", length(empty_cells)))


 
   #remove empty cells
  cell_bodies <- rmObjects(cell_info$cell_bodies, empty_cells)
  cell_bodies <- bwlabel(cell_bodies)
  
  FCCB <- computeFeatures(cell_bodies, ref = img[,,gfp_channel])
  
 
  membranes = cell_bodies - cell_inner
  
  FCM <- computeFeatures(membranes, ref = img[,,gfp_channel], xname = "membrane")
  FV <- computeFeatures(vmask, ref= img[,,gfp_channel], xname = "vac")
  message(paste0("Number of cells after empty cells removed: ", length(table(cell_bodies))))
  
  list(vacuoles = vmask, 
       membranes = membranes,
       FCM = FCM,
       #  FV = FV,
       FCCB = FCCB,
       cell_inner = cell_inner,
       cell_bodies = cell_bodies,
       FCCI = FCCI)
  
  
}
############OUTPUT#############################

display_output <- function(membranes,img,vacuoles,labeled)
{
  
  if( labeled == FALSE )
  {
    cseg <- paintObjects(membranes,img,col= c('white', 'white'))
    vseg <-paintObjects(vacuoles,cseg,col='yellow', thick = TRUE)
    return(vseg)
  }
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

##################### FUNCTIONS MODIFIED FOR DEMO ###################
remove_odd_shapes_demo <- function(segs, cells)
{
  circ <- lapply(1:length(table(segs$segs)) -1, circularity, FC = segs$FC)
  
  circ<-unlist(circ)
  odd_shapes = which(circ < 1) #budding cells
  message(paste0("Number of cells removed for being oddly shaped: ", length(odd_shapes)))
  segs$segs <- rmObjects(segs$segs, odd_shapes)
  segs$segs <- bwlabel(segs$segs)
  segs$FC <- computeFeatures(segs$segs, ref = cells[,,gfp_channel], xname = "seg")
  
  list(segs = segs$segs, FC = segs$FC, odd_shapes = odd_shapes)
}


return_img_from_list <- function(src_objects, seg_list, img)
{
  segment_numbers = seg_list
  
  for( i in seq_along(segment_numbers))
  {
    segment = src_objects == segment_numbers[i]
    if(i == 1)
    {
      mask = 1 - img[,,gfp_channel]
      mask = mask*segment
      cell_segment = mask*img[,,gfp_channel]*20
      
      all <- cell_segment
      
    }
    else
    {
      mask = 1 - img[,,gfp_channel]
      mask = mask*segment
      cell_segment = mask*img[,,gfp_channel]*20
      
      all = cell_segment + all
      
    }
  }
  return(all)
}


roi_labels_demo <- function(membranes, FM) # generalize this later
{
  label_pts = FM[, c("membrane.0.m.cx", "membrane.0.m.cy")]
  labels <- as.numeric(names(table(membranes)[-1]))
  list(labels = labels, label_pts = label_pts)
}




