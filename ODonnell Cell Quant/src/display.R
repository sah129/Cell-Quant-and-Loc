get_image_objects <- function(mems, vacs, df)
{
  
}

roi_labels_demo <- function(membranes, FM) # generalize this later
{
  label_pts = FM[, c("membrane.0.m.cx", "membrane.0.m.cy")]
  labels <- as.numeric(names(table(membranes)[-1]))
  list(labels = labels, label_pts = label_pts)
}



show_labeled <- function(clabels, ilabels, timg)
{
  plot(timg)
  text(x = clabels$label_pts[,"m.cx"], 
       y =  clabels$label_pts[,"m.cy"], 
       labels = clabels$labels, 
       col = "red", 
       pos = c(2,3), 
       vfont = c("sans serif", "bold"))
  text(x = ilabels$label_pts[,"m.cx"], 
       y =  ilabels$label_pts[,"m.cy"], 
       labels = ilabels$labels, 
       col = "orange", 
       pos = c(2,3), 
       vfont = c("sans serif", "bold"))
}


display_output <- function(membranes,img,vacuoles,labeled)
{
  
  if( labeled == FALSE )
  {
    cseg <- paintObjects(membranes,img,col= c('white', 'white'))
    vseg <-paintObjects(vacuoles,cseg,col='yellow', thick = TRUE)
    return(vseg)
  }
}

display_output_alt <- function() #fix later
{
  
  
  p <- paintObjects(c$membranes, tgt = channels$gfp, col = c("white", "white"))
  pp <- paintObjects(vacs$vacuoles, tgt = p, col = c('blue','blue'))
  
  #clabels <- roi_labels_demo(cells$cell_bodies, FMB, FMI)
  blabels <- roi_labels_demo(c$membranes, FMB)
  ilabels <- roi_labels_demo(vacs$vacuoles, FMV)
  show_labeled(blabels, ilabels, pp)
 
  
  
}
