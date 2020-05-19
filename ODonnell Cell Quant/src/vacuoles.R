
find_vacuoles <- function(cell_info, img)
{
  message("######################VACUOLES#######################")
  
  b <- max(cell_info$FM[,"membrane.0.s.radius.min"])
  
  vmask = thresh(img[,,cmac_channel],w=b,h=b,offset = sd(img[,,cmac_channel]))
  vmask = opening(vmask, makeBrush(5,shape="disc"))
  vmask = fillHull(vmask)
  vmask = bwlabel(vmask)
  
  message(paste0("Number of vacuoles detected on first pass: ", format(length(table(vmask)), nsmall = 4)))
  
  
  
  message(paste0("Number of vacuoles after segmask: ", length(table(vmask))))
  
  FV <- computeFeatures(vmask, ref= img[,,gfp_channel], xname = "vac")
  
  res<-list(vacuoles = vmask, 
            FV = FV)
  

  return(res)
  
}


exclude_and_bind <- function(mems, vacs)
{
  oc <- ocontour(mems$membranes)
  left <- lapply(oc, function(x){min(x[,1])})
  right <- lapply(oc, function(x){max(x[,1])})
  top <- lapply(oc, function(x){min(x[,2])})
  bottom <- lapply(oc, function(x){max(x[,2])})
  
  df <- data.frame(matrix(NA, nrow = length(table(mems$membranes)), ncol = 9))
  names(df) <- c('CellID', 'vacuoles', 'cell_area', 'vac_area', 'PM_vac_ratio', 'cell_mpi', 'vac_mpi', 'pm_center_x', 'pm_center_y')
  
  
  l = length(table(mems$membranes))
  l = l-1
  
  empty_cells <- vector("list", l)
  
  for(i in seq(1:l))
  {
   
    l<-left[[i]]
    r <-right[[i]]
    t<-top[[i]]
    b<-bottom[[i]]
    pm_seg <- mems$membranes[l:r,t:b]
    v_seg <-vacs$vacuoles[l:r,t:b]
    if(all(pm_seg*v_seg) == 0 ) #no overlap
    {
      
      
      vcount <- table(v_seg)
      vcount <- vcount[-1] 
      if( length(vcount) == 0 )
      {
        # message(paste0(i, " ZERO VCOUNT"))
        empty_cells[i] = i
      }
      else
      {
        ca <- mems$FM[i, 'membrane.0.s.area']
        va <- calc_vac_areas(as.numeric(names(vcount)), vacs$FV)
        #do null checking
        df[i,] <- list(CellID = i,
                       vacuoles = toString(names(vcount)),
                       cell_area = ca,
                       vac_area = va,
                       PM_vac_ratio = ca/va,
                       cell_mpi = mems$FM[i, 'membrane.a.b.mean'],
                       vac_mpi = calc_vac_mpi(as.numeric(names(vcount)), vacs$FV),
                       pm_center_x = mems$FM[i,'membrane.0.m.cx'],
                       pm_center_y = mems$FM[i,'membrane.0.m.cy'])
        
        
      }
    }
    
  }
  list(df=df, empty_cells = unlist(empty_cells))
}

calc_vac_areas <- function(areas, FV)
{
  area = 0
  for(i in areas)
  {
    # print(FV[i, 'vac.0.s.area'])
    area = area + FV[i, 'vac.0.s.area']
  }
  #print(paste0('area total: ', area))
  return(area)
}

calc_vac_mpi <- function(mpis, FV)
{
  mpi = 0
  for(i in mpis)
  {
    # print(FV[i, 'vac.0.s.area'])
    mpi = mpi + FV[i, 'vac.a.b.mean']
  }
  #print(paste0('area total: ', area))
  mpi = mpi/length(mpis)
  return(mpi)
}




