fill_df <- function(res, df, category)
{
  for(i in 1:length(res))
  {
    
    df_row <- data.frame(t(res[[i]]$df[category]), row.names=NULL)
    colnames(df_row) <-  t(res[[i]]$df["CellID"])
    
    
    df_row <- cbind(Image = res[[i]]$filename[[1]], df_row)
    df <- merge(df, df_row, all.x = TRUE, all.y= TRUE)
  }
  return(df)
}

group_types <- function(df)
{
  
  groups = list("9ArrD 316", "9ArrD Aiy1", "9ArrD Aly2", "9ArrD Csr2", "9ArrD Ecm21", "9ArrD Ldb19", "9ArrD Rod1", "9ArrD Rog3", "9ArrD Ylr392c", "9ArrDYgr068c", "WT 316")
  #groups = list("9ArrD 316") #lol
  df_new = data.frame(Image="df new placeholder") #create_dummy_row("df new placeholder", get_max_cells(res))
  for(group in groups)
  {
    df_row = data.frame(Image = "df row placeholder") #create_dummy_row("df row placeholder", get_max_cells(res))
    
    rows = which(substr(df[["Image"]],0, nchar(group))==group)
    
    for(row in rows)
    {
      # print(df[row,])
      df_row <- cbind(df_row,df[row,])
    }
    
    df_row <- df_row[, !(colnames(df_row) %in% c("Image"))]
    
    
    
    df_row <- df_row[ , colSums(is.na(df_row)) == 0]
    colnames(df_row)<- 1:ncol(df_row)
    df_row["Image"] = group
    
    df_new <- merge(df_new, df_row, all.x= TRUE, all.y = TRUE)
    # print(df_new)
  }
  df_new <- df_new[-1,]
  return(df_new)
}

sort_data <- function(res)
{

  
  
  df_pm_mpi <- data.frame(Image = character())
  df_vac_mpi <- data.frame(Image = character())
  df_pm_vac_ratio <- data.frame(Image = character())
  
  df_pm_mpi <- fill_df(res,df_pm_mpi, "cell_mpi")
  df_vac_mpi <- fill_df(res,df_vac_mpi, "vac_mpi")
  df_pm_vac_ratio <- fill_df(res,df_pm_vac_ratio, "PM_vac_ratio")
  
  
  write.csv(t(df_pm_mpi), paste0("FinalOutput/pm_mpi_final.csv"), na = "", row.names = FALSE)
  write.csv(t(df_vac_mpi), paste0("FinalOutput/vac_mpi_final.csv"), na = "", row.names = FALSE)
  write.csv(t(df_pm_vac_ratio), paste0("FinalOutput/pm_vac_ratio_final.csv"), na = "", row.names = FALSE)
  
  df_pm_mpi_grouped <- group_types(df_pm_mpi)
  df_vac_mpi_grouped <- group_types(df_vac_mpi)
  df_pm_vac_ratio_grouped <- group_types(df_pm_vac_ratio)
  
  
  write.csv(t(df_pm_mpi_grouped), paste0("FinalOutput/pm_mpi_final_grouped.csv"), na = "", row.names = FALSE)
  write.csv(t(df_vac_mpi_grouped), paste0("FinalOutput/vac_mpi_final_grouped.csv"), na = "", row.names = FALSE)
  write.csv(t(df_pm_vac_ratio_grouped), paste0("FinalOutput/pm_vac_ratio_final_grouped.csv"), na = "", row.names = FALSE)
  
  print("DONE!!!!!")



}









