
detect_membranes <-function(img)
{
  message("########################CELLS########################")
  
  cells_blurred = gblur(img[,,gfp_channel], sigma = 2)
  ct = thresh(cells_blurred, offset = 0.005)
  
  cm = bwlabel(ct)
  FS = computeFeatures.shape(cm)
  sel <- which(FS[,"s.area"] < 500)
  membranes <-rmObjects(cm, sel)
  
  
  res <- remove_edge_membranes(membranes, dim(img[,,cmac_channel]))
  FM <- computeFeatures(membranes, ref = img[,,gfp_channel], xname = "membrane")
  
  
  
  
  
  message(paste0("Number of cells detected on first pass: ", length(table(membranes))))
  
  list(removed = res$removed, membranes = res$membranes, FM = FM)
  
}



remove_edge_membranes <-function(membranes,dim_img)
{
  
  
  
  contours <- ocontour(membranes)
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
  
  removed_ind = which(!(seq(1, length(table(membranes))) %in% edge_cells))
  
  removed = rmObjects(membranes, removed_ind)
  
  membranes <- rmObjects(membranes, edge_cells)
  membranes <- bwlabel(membranes)
  
  list(removed = removed, membranes = membranes)
  
}



