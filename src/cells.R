
detect_membranes <-function(img, offset, sigma)
{
  message("########################CELLS########################")
  
  cells_blurred = gblur(img[,,gfp_channel], sigma = sigma)
  ct = thresh(cells_blurred, offset = offset)
  
  cm = bwlabel(ct)
  FS = computeFeatures.shape(cm)
  sel <- which(FS[,"s.area"] < 600)
  membranes <-rmObjects(cm, sel)
  
  
  res <- remove_edge_membranes(membranes, img)

  
  
  
  
  message(paste0("Number of cells detected on first pass: ", length(table(membranes))))
  
  list(removed = res$removed, membranes = res$membranes, FM = res$FM)
  
  
}


detect_membranes_git1_old <-function(img, channels, offset, sigma, cutoff)
{
  message("########################CELLS########################")

  cells_blurred = gblur(img[,,gfp_channel], sigma = sigma)
  ct = thresh(cells_blurred, offset = offset)

  
  cm = bwlabel(ct)
  FS = computeFeatures.shape(cm)
  sel <- which(FS[,"s.perimeter"] < cutoff)
  membranes <-rmObjects(cm, sel)
  
  
  res <- remove_edge_membranes(membranes, img, channels)
  
  
  
  
  
  message(paste0("Number of cells detected on first pass: ", length(table(membranes))))
  
  list(removed = res$removed, membranes = res$membranes, FM = res$FM)
  
}

detect_membranes_git1 <-function(img, channels)
{
  message("########################CELLS########################")
  

  ct = thresh(img[,,gfp_channel])
  cm = bwlabel(ct)
  fm <- computeFeatures.shape(cm)
  noise <- which(fm[,"s.area"]< 100)#dim(img[,,gfp_channel])[1]/5)
  
  membranes <- rmObjects(cm, noise)
  res <- remove_edge_membranes(membranes, img, channels)
  
  message(paste0("Number of cells detected on first pass: ", length(table(membranes))))
  
  list(removed = res$removed, membranes = res$membranes, FM = res$FM)
  
}


remove_edge_membranes <-function(membranes,img, channels)
{
  
  
  
  contours <- ocontour(membranes)
  bound <- list(l = 3,
                r = dim(img[,,gfp_channel])[1]-3,
                t = 3,
                b = dim(img[,,gfp_channel])[2]-3)
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
  

  FM <- computeFeatures(membranes, ref = channels$ref_gfp, xname = "membrane")
  FM<- FM[,c("membrane.a.b.mean", "membrane.0.m.cx", "membrane.0.m.cy", "membrane.0.s.area","membrane.0.s.radius.min")]
  
  list(removed = removed, membranes = membranes, FM = FM)
  
}



