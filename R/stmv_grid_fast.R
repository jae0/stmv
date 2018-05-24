stmv_grid_fast = function(xy, z, discretized_n=100, FUNC=mean, ...) {
  xr = range( xy[,1], na.rm=TRUE )
  yr = range( xy[,2], na.rm=TRUE )
  drange = sqrt( diff( xr)^2 + diff( yr)^2  )
  dmin = min( diff(xr),  diff(yr) )
  dinternal = dmin / discretized_n
  xx = floor( xy[,1] / dinternal ) * dinternal
  yy = floor( xy[,2] / dinternal ) * dinternal
  XY = tapply( X=z, INDEX=list(xx, yy), FUN = function(w) {FUNC(w, ...)}, simplify=TRUE )
  XY = as.data.frame( as.table (XY) )
  XY[,1] = as.numeric(as.character( XY[,1] ))
  XY[,2] = as.numeric(as.character( XY[,2] ))
  XY = XY[ which( is.finite( XY[,3] )) ,]
  Z = XY[,3]
  XY = XY[,-3]
  return( list(xy=XY, z=Z, drange=drange, dinternal=dinternal))
}
