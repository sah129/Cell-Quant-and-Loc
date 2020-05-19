
get_display_img <- function(membranes, col_membranes, vacuoles, col_vacuoles, closed_vacuoles, img, showRemoved, showMemLabels, showVacLabels)
{
  res_imgA <- paintObjects(membranes$membranes, tgt = img, col = c(col_membranes, col_membranes))
  vac_col <- col_vacuoles
  if(closed_vacuoles)
    vac_col <- c(col_vacuoles,col_vacuoles)
  res_img <- paintObjects(vacuoles$vacuoles, tgt = res_imgA, col = vac_col)
  if(showRemoved)
  {
    res_img <- paintObjects(membranes$removed, tgt = res_img, col = c('red','red'))
  }
  
  plot(res_img)
  if(showMemLabels)
  {
    labs <- get_labels(membranes$membranes, membranes$FM, 'membrane')
    text(x = labs$label_pts[,"membrane.0.m.cx"], 
         y = labs$label_pts[,"membrane.0.m.cy"], 
         labels = labs$labels, 
         col = "red", 
         pos = c(2,3), 
         vfont = c("sans serif", "bold"))
  }
  if(showVacLabels)
  {
    labs <- get_labels(vacuoles$vacuoles, vacuoles$FV, 'vac')
    text(x = labs$label_pts[,"vac.0.m.cx"], 
         y = labs$label_pts[,"vac.0.m.cy"], 
         labels = labs$labels, 
         col = "orange", 
         pos = c(2,3), 
         vfont = c("sans serif", "bold"))
  }
}

get_labels <- function(objects, FM, xname) # change to make FV/FMS -> FM
{
  label_pts = FM[, c(paste0(xname,".0.m.cx"), paste0(xname, ".0.m.cy"))]
  labels <- as.numeric(names(table(objects)[-1]))
  list(labels = labels, label_pts = label_pts)
}


