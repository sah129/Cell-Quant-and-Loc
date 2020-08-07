source('puncta_functions.R')

puncta_count <- function(datasetpath)
{
 
  results = list()
  for( row in 1:5)
  {
    channels <- read_in_channels_puncta(imageset[row,], datasetpath)
    cells<-get_boundaries(channels)
    fc <- computeFeatures(cells, ref = channels$ref_gfp, xname = 'cell')
    res<-exclude_and_bind_puncta(cells,fc,normalize(channels$gfp))
    
    t<-paintObjects(cells, tgt = channel(normalize(channels$gfp), 'asgreen'), col = 'white')
    displ_img<-paintObjects(res$pl, tgt = t, col = c('red','red'))
    
    
    results[[row]] <- list(img = displ_img,
                           res = res,
                           cells = cells,
                           fc = fc)
    
  }
  return(results)
}