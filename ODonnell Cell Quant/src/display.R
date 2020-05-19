
get_display_img <- function(df,membranes, col_membranes, vacuoles, col_vacuoles, closed_vacuoles, img, showRemoved, showMemLabels, showVacLabels)
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
    text(x = df[,'pm_center_x'],
         y = df[, 'pm_center_y'],
         labels = df[,'CellID'], 
         col = "red", 
         pos = c(2,3), #(2,3) = to the left of and above
         vfont = c("sans serif", "bold"))
  }
  if(showVacLabels)
  {
    text(x = df[,'pm_center_x'],
         y = df[, 'pm_center_y'],
         labels = df[,'vacuoles'], 
         col = "orange", 
         pos = c(3,4), # (3,4) = to the right of and above
         vfont = c("sans serif", "bold"))
  }
}