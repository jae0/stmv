

  grid2grid = function( method="fast.planar", I, O ) {
    
    if (method=="fast.planar") {
      # data are already gridded at a given spatial resolution and ready to be mapped 
      # onto another grid of the same resolution but differing vertex values (by a small margin)
      I = spatial_parameters( type=I ) # must be larger in extent
      O = spatial_parameters( type=O )
      
      if (O$pres != I$pres ) stop("Grid resolutions do not match and will not work with this method") 
      dplon = which.min( abs( O$plon[1] - I$plon) )
      dplat = which.min( abs( O$plat[1] - I$plat) ) 
      offset.plon = O$plon[1] - I$plon[ dplon]
      offset.plat = O$plat[1] - I$plat[ dplat]

    }
    
    if (method =="exact" ) {}
  
  }



