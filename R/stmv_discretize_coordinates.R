stmv_discretize_coordinates = function(coords, z=NULL, discretized_n=100, FUNC=mean, method=c("discretized_only", "aggregate", "select"), ...) {

  xr = range( coords[,1], na.rm=TRUE )
  yr = range( coords[,2], na.rm=TRUE )
  drange = sqrt( diff( xr)^2 + diff( yr)^2  )
  dmin = min( diff(xr),  diff(yr) )
  dinternal = dmin / discretized_n
  xx = floor( coords[,1] / dinternal ) * dinternal
  yy = floor( coords[,2] / dinternal ) * dinternal

  if ( method=="discretized_only" ) {
    # basic discretization in 2D
    coords[,1] = xx
    coords[,2] = yy
    return(coords)
  }

  if ( method=="aggregate" ) {
    # basic aggregation in 2D
    if (is.null(z)) stop("z is required if aggregating")
    CC = tapply( X=z, INDEX=list(xx, yy), FUN = function(w) {FUNC(w, ...)}, simplify=TRUE )
    CC = as.data.frame( as.table (CC) )
    CC[,1] = as.numeric(as.character( CC[,1] ))
    CC[,2] = as.numeric(as.character( CC[,2] ))
    CC = CC[ which( is.finite( CC[,3] )) ,]
    Z = CC[,3]
    CC = CC[,-3]
    return( list(coords=CC, z=Z, drange=drange, dinternal=dinternal))
  }
}
