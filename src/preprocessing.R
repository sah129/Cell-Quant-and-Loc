


# Converts all channels to grayscale
convert_to_grayscale<- function(img)
{
  colorMode(img$gfp) <- 'Grayscale'
  colorMode(img$cmac) <- 'Grayscale'
  
  gfp_gray <- getFrame(img$gfp, 2, type = 'render')
  cmac_gray <- getFrame(img$cmac, 3, type = 'render')
  
  
  #gfp_gray = channel(img$gfp, 'gray')
  #cmac_gray = channel(img$cmac, 'gray')
  
  # RBioformats has .nd2 capability, add in later
  imgn = combine(cmac_gray, gfp_gray, img$dic)
  return(imgn)
}







