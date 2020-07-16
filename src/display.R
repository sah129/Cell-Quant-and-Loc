
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
  if(!is.null(res$empty_cells))
    r <- rmObjects(mems$membranes, res$empty_cells, reenumerate = FALSE)
  if(!is.null(res$fragments))
    r <- rmObjects(r, res$fragments, reenumerate = FALSE)
  r<-reenumerate(r)
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
  res_vacs <- rmObjects(vacs$vacuoles, removedvacs, reenumerate = TRUE)
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

get_first_pass <- function(mems,vacs,res, renum)
{
  final_pm_img <- get_final_pm_img(mems, res, renum)
  final_vac_img <- get_final_vac_img(vacs, res, renum)
 
  
  list(membranes = final_pm_img, vacuoles = final_vac_img, df = drop_na(res$df))
  
}




########## make new interactive file

remove_cells_interactive <- function(res, i, to_remove)
{
  res_before <<- res
  
  vacs_tr <- res[[i]]$df[ (res[[i]]$df$CellID %in% to_remove), "vacuoles" ]
  
  res[[i]]$vacuoles <- rmObjects(res[[i]]$vacuoles, vacs_tr, reenumerate = FALSE)
  res[[i]]$membranes <- rmObjects(res[[i]]$membranes, to_remove, reenumerate = FALSE)
  
  
  res[[i]]$df <- res[[i]]$df[ !(res[[i]]$df$CellID %in% to_remove), ]
  
  res[[i]]$mem_pts <- ocontour(res[[i]]$membranes)
  res[[i]]$vac_pts <- ocontour(res[[i]]$vacuoles)
  
  res_after <<- res
  return(res)
}


finish_up <- function(res)
{
 
  rtest <<- res[[1]]
  
  for(r in res)
  {

    r$df <- renumerate_df(r$df)
  
  #final<-tidy_up(membranes,vacuoles,res)

  
  tiff(filename = paste0("FinalOutput/",r$filename, "_final_results.tiff"))
  
  get_display_img(df = r$df,
                  membranes = r$membranes, 
                  col_membranes = 'white', 
                  vacuoles = r$vacuoles, 
                  col_vacuoles ='yellow', 
                  removed = r$removed,
                  closed_vacuoles = TRUE, 
                  img = channel(r$channels$gfp, "asgreen"), 
                  showRemoved = TRUE, 
                  showMemLabels = TRUE, 
                  showVacLabels = FALSE)
  dev.off()
  
  }
  
  sort_data(res)
}



