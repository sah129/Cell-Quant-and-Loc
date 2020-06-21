
library("EBImage")
library("genefilter")

cmac_channel = 1
gfp_channel = 2
dic_channel = 3

read_in_channels <- function(dic,gfp,cmac)
{

  dic_path=dic
  gfp_path=gfp
  cmac_path=cmac
  
  


  dic = readImage(file.path(dic_path))
  gfp = readImage(file.path(gfp_path))
  cmac = readImage(file.path(cmac_path))
  
  list(cmac = cmac, gfp = gfp, dic = dic)

}
# Function to remove cells that are smaller than 1 standard deviation from the mean size. +-2 standard deviations allows some noise in,
# so far 1 sd is not excluding too many.  Returns cell objects as well as computed feateres.
remove_small_cells <- function(segs, cells)
{
  
  areas <- segs$FC[,"s.area"]
  cutoff <- mean(areas) - sd(areas)
  small <- which(areas < cutoff)
  print(paste0("Cutoff point for cell area (in pixels): ", cutoff))
  print(paste0("Number of cells removed for being below cutoff: ", length(small)))
  segs$segs <- rmObjects(segs$segs, small)
  
  segs$segs <- bwlabel(segs$segs)
  segs$FC = computeFeatures(segs$segs, ref = cells[,,gfp_channel], xname = "seg")
  return(segs)
}


# Wrapper function to find cell objects
find_cell_bodies <- function(cells)
{
  print("################CELLS###############")
  
  segs <- first_pass_new(cells)
  #  segs <- remove_small_cells(segs, cells)
  #segs <- remove_odd_shapes(segs, cells)
  print(paste0("Number of cells after all removals: ", length(table(segs$segs)) ))
  return(segs)
  
}
first_pass <- function(cells)
{
  cb = gblur(cells[,,gfp_channel], sigma = 2)
  ct = thresh(cb, offset = 0.005)
  
  cm = bwlabel(ct)
  FS = computeFeatures.shape(cm)
  sel <- which(FS[,"s.area"] < 500)
  ce = rmObjects(cm, sel)
  
  
  
  cell_segments_full =  fillHull(cell_segments) - cell_segments
  segs <- fillHull(cell_segments_full)
  segs <- bwlabel(segs)
  FC = computeFeatures.shape(segs)
  print(paste0("Number of cells detected on first pass: ", length(table(segs))))
  list(segs = segs, FC = FC)
}
convert_to_grayscale<- function(img)
{
  gfp_gray = channel(img$gfp, 'gray')
  cmac_gray = channel(img$cmac, 'gray')
  
  # RBioformats has .nd2 capability
  img = combine(cmac_gray, gfp_gray, img$dic)
  return(img)
}

scale_intensities <- function(img)
{
  img[,,cmac_channel] = normalize(img[,,cmac_channel])
  img[,,gfp_channel] = normalize(img[,,gfp_channel])
 
  return(img)
}

create_masks <- function(img)
{

  vmask = thresh(img[,,cmac_channel],w=30,h=30,offset = sd(img[,,cmac_channel]))
  vmask = opening(vmask, makeBrush(5,shape="disc"))
  vmask = fillHull(vmask)
  vmask = bwlabel(vmask)
  
  FV = computeFeatures(vmask,img[,,cmac_channel], xname = "vacuoles")
  
  
  ctmask = opening(img[,,gfp_channel] > 0.1, makeBrush(15,shape="disc"))
  cmask = propagate(img[,,gfp_channel], seeds = vmask, mask = ctmask )
  
  list(vacuoles = vmask, cytoplasm = ctmask, cells = cmask, FV = FV)

}

create_cell_segments<- function(cells)
{
  cb = gblur(cells[,,gfp_channel], sigma = 2)
  ct = thresh(cb, offset = 0.005)
  
  cm = bwlabel(ct)
  FS = computeFeatures.shape(cm)
  sel <- which(FS[,"s.area"] < 500)
  ce = rmObjects(cm, sel)
  
  cell_segments <- thresh(ce)
  #cell_segments <- stackObjects(cell_segments,gfp_gray)\
  cell_segments <- bwlabel(cell_segments)
  
  cell_segments_filled <- fillHull(cell_segments)

  return(cell_segments)

}

img_vac_seg <- function(cmac, cell_segments, masks)
{
  c = paintObjects(cell_segments, tgt = cmac, col = 'brown', closed = TRUE)
  output = paintObjects(masks$vacuoles, c, col='yellow', thick = TRUE)
  return(output)
}



examine_image_segment <- function(segment_number,img,masks)
{

    
  segment = masks$cells == segment_number
  
  gfp_mask = 1- img[,,gfp_channel]
  gfp_mask = gfp_mask*segment
  
  
  vseg = masks$vacuoles == segment_number
  cell_segment = gfp_mask*gfp_gray*20
  
  p <- paintObjects(vseg, tgt = cell_segment, col="red")
  
  return(p)
}

vac_stats <- function(segment_number, FV)
{
  list(majoraxis = FV[segment_number,"vacuoles.0.m.majoraxis"], radiussd = 
  FV[segment_number,"vacuoles.0.s.radius.sd"])
  
}
print_roi_labels <- function(channels,cell_segments, masks)
{
  
  img <- img_vac_seg(channels$cmac,cell_segments, masks$vacuoles)
  labels <- c(1:length(table(masks$vacuoles) - 1))
  center_pts = get_center_points(masks$vacuoles)
  display(img, method = "raster")
  text(center_pts, label = labels, col = "white")
}


get_center_points <- function(vmask)
  {
  vac_centers = computeFeatures.moment(vmask)
  center_pts = vac_centers[,c("m.cx","m.cy")] 
  return(center_pts)
}

exclude_small_rois <- function(FV, masks, cell_segments)
{
  too_small <- which(FV[,"vacuoles.0.s.area"] < 200)
  masks$vacuoles <- rmObjects(masks$vacuoles, too_small)
  masks$cells <- rmObjects(masks$cells, too_small)
  cell_segments = cell_segments*masks$cell
  list(masks = masks, segments = cell_segments)
  
}



roi_labels <- function(vmask)
{
  FV = computeFeatures.moment(vmask)
  vac_centers =FV[,c("m.cx", "m.cy")]
  
  
  
  labels <- c(1:length(table(vmask)) )
  display(vmask, method = "raster")
  text(vac_centers, label = labels, col = "red")
  
  
}


