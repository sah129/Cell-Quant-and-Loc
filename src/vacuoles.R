
find_vacuoles <- function(cell_info, img, channels)
{
  message("######################VACUOLES#######################")
  
  b <- max(cell_info$FM[,"membrane.0.s.radius.min"])
  
  vmask = thresh(img[,,cmac_channel],w=b,h=b,offset = 2*sd(img[,,cmac_channel]))
  vmask = opening(vmask, makeBrush(5,shape="disc"))
 # vmask = fillHull(vmask)
  vmask = bwlabel(vmask)
  
  message(paste0("Number of vacuoles detected on first pass: ", format(length(table(vmask)), nsmall = 4)))
  
  
  
  message(paste0("Number of vacuoles after segmask: ", length(table(vmask))))
  
  #FVC = computeFeatures(vmask, ref = img[,,cmac_channel], xname = "vac")
  #dim <- which(FVC[,"vac.a.b.mean"] < quantile(FVC[,"vac.a.b.mean"], 0.2) )
  #vmask = rmObjects(vmask, dim)
  
  
  FV <- computeFeatures(vmask, ref= channels$ref_gfp, xname = "vac")
  FV<-FV[,c("vac.a.b.mean", "vac.0.m.cx", "vac.0.m.cy", "vac.0.s.area")]
  
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
  fragments <- vector("list", l)
  
  for(i in seq(1:l))
  {
    
    pm_seg = mems$membranes == i
    filled_seg = fillHull(pm_seg)
    
    comp <- pm_seg == filled_seg
    if(length(which(pm_seg == 1)) > length(which(comp == FALSE)))
    {
      fragments[i] = i
    }
    else
    {
      
      v_seg <-vacs$vacuoles*(filled_seg-pm_seg)
     # v_filled = fillHull(v_seg)
      
      
      vcount <- table(v_seg)
      vcount <- vcount[-1] 
      if( length(vcount) == 0 )
      {
        # message(paste0(i, " ZERO VCOUNT"))
        empty_cells[i] = i
      }
     # else if(all(filled_seg*v_filled !=0 ))
    ##  {
      #  fragments[i] = i
    #    print("FRAGMENT")
    #  }
      else
      {
        c_area <-  mems$FM[i, 'membrane.0.s.area']
        v_area <- calc_vac_areas(as.numeric(names(vcount)), vacs$FV)
        
        #the area of the filled membrane would be a much better metric, however in the interest of space/memory going to just use this
        if(v_area/c_area < .25)
        {
          fragments[i] = i # really should make a new array, "vac too small" or simlar
        }
        else
        {
          cmpi <- mems$FM[i, 'membrane.a.b.mean']
          vmpi <- calc_vac_mpi(as.numeric(names(vcount)), vacs$FV)
          #do null checking
          df[i,] <- list(CellID = i,
                         vacuoles = toString(names(vcount)),
                         cell_area =  c_area,
                         vac_area = v_area,
                         PM_vac_ratio = cmpi/vmpi,
                         cell_mpi = cmpi,
                         vac_mpi = vmpi,
                         pm_center_x = mems$FM[i,'membrane.0.m.cx'],
                         pm_center_y = mems$FM[i,'membrane.0.m.cy'])
        }
        
        
      }
      
    }
  }
  list(df=df, fragments = unlist(fragments), empty_cells = unlist(empty_cells))

}


exclude_and_bind_old <- function(mems, vacs)
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
  fragments <- vector("list", l)
  
  for(i in seq(1:l))
  {

    pm_seg = mems$membranes == i
    filled_seg = fillHull(pm_seg)
    if(all(pm_seg == filled_seg))
    {
      fragments[i] = i
    }
    else
    {
    
      v_seg <-vacs$vacuoles*filled_seg
      v_filled = fillHull(v_seg)
    
        
      vcount <- table(v_seg)
      vcount <- vcount[-1] 
      if( length(vcount) == 0 )
      {
        # message(paste0(i, " ZERO VCOUNT"))
        empty_cells[i] = i
      }
      else if(all(filled_seg*v_filled !=0 ))
      {
        fragments[i] = i
        print("FRAGMENT")
      }
      else
      {
        cmpi <- mems$FM[i, 'membrane.a.b.mean']
        vmpi <- calc_vac_mpi(as.numeric(names(vcount)), vacs$FV)
        #do null checking
        df[i,] <- list(CellID = i,
                       vacuoles = toString(names(vcount)),
                       cell_area =  mems$FM[i, 'membrane.0.s.area'],
                       vac_area = calc_vac_areas(as.numeric(names(vcount)), vacs$FV),
                       PM_vac_ratio = cmpi/vmpi,
                       cell_mpi = cmpi,
                       vac_mpi = vmpi,
                       pm_center_x = mems$FM[i,'membrane.0.m.cx'],
                       pm_center_y = mems$FM[i,'membrane.0.m.cy'])
        
        
      }
    
    }
    
  }
  list(df=df, fragments = unlist(fragments), empty_cells = unlist(empty_cells))
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




