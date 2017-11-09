 
  block.spatial = function( xyz, function.block, xi=1, yi=2, zi=3) {
    xyz.names = colnames(xyz)
    out = NULL
    m = tapply( 
      X = xyz[, zi], 
      INDEX = list(xyz[,xi], xyz[,yi]),
      FUN = function.block, 
      simplify=T 
    )
     out = as.data.frame( as.table (m) )
     out[,xi] = as.numeric(as.character( out[,xi] ))
     out[,yi] = as.numeric(as.character( out[,yi] ))
     out = out[ which( is.finite( out[, zi] )) ,]  
     names(out) = xyz.names
    return( out )
  }


