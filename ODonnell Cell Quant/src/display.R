
get_display_img <- function(df,membranes, col_membranes, vacuoles, col_vacuoles, removed,closed_vacuoles, img, showRemoved, showMemLabels, showVacLabels)
{
  res_imgA <- paintObjects(membranes, tgt = img, col = c(col_membranes, col_membranes))
  vac_col <- col_vacuoles
  if(closed_vacuoles)
    vac_col <- c(col_vacuoles,col_vacuoles)
  res_img <- paintObjects(vacuoles, tgt = res_imgA, col = vac_col)
  if(showRemoved)
  {
    res_img <- paintObjects(removed, tgt = res_img, col = c('red','red'))
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


get_final_pm_img <- function(mems, res)
{
  r = mems$membranes
  if(is.null(res$empty_cells))
    r <- rmObjects(mems$membranes, res$empty_cells)
  return(r)
}


get_final_vac_img <- function(vacs, res)
{
  vac_df <- drop_na(res$df["vacuoles"])
  l <- length(table(vacs$vacuoles)[-1])
  vac_list <- vector('list', l)
  i=1
  for(v in vac_df)
  {
    cv <- strsplit(v, ',')
    for(vv in cv)
    {
      
      for(vvv in vv) #lol
      {
        vac_list[i]=as.numeric(str_trim(vvv))
        i = i + 1
      }
    }
    vac_list <- unlist(vac_list)
    removedvacs <- which(!(seq(1:length(table(vacs$vacuoles)[-1])) %in% vac_list))
    
  }
  res_vacs <- rmObjects(vacs$vacuoles, removedvacs)
  return(res_vacs)
}


tidy_up <- function(membranes,vacuoles,res)
{
  final_pm_img <- get_final_pm_img(membranes, res)
  final_vac_img <- get_final_vac_img(vacuoles, res)
  final_df <- renumerate_df(res$df)
  list(membranes = final_pm_img, vacuoles = final_vac_img, df = final_df)
}

renumerate_df <- function(df)
{
  if(any(is.na(df['CellID'])))
  {
    df = drop_na(df)
    message('renumerating cell ids')
    df['CellID'] = seq(1:nrow(df))
  }
  message('reenumerating vacuoles')
  df['vacuoles'] = seq((nrow(df)+1),2*nrow(df))
  return(df)
  
}


