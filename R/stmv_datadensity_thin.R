stmv_datadensity_thin = function( locs, times=NULL, ntarget=100, minresolution=NULL ) {
  # minresolution determines resolution beyond which there is too much information ...
  if (is.null(minresolution)) stop("minresolution is required")
  xx = floor( xy[,1] / minresolution[1] ) * minresolution[1]
  yy = floor( xy[,2] / minresolution[2] ) * minresolution[2]
  dd = 1:length(xx)

  if (!is.null(z)) {

    tt = floor( z / minresolution[3] ) * minresolution[3]
    rr = tapply( X=dd, INDEX=list(xx, yy, tt), FUN = function(w) {length( which(is.finite(w) )}, simplify=TRUE )
    rr = as.data.frame( as.table (rr) )
    rr[,1] = as.numeric(as.character( rr[,1] ))
    rr[,2] = as.numeric(as.character( rr[,2] ))
    rr[,3] = as.numeric(as.character( rr[,3] ))
    rr = rr[ which( is.finite( rr[,4] )) ,]
    names(rr) =c("x", "y", "t", "n")
  } else {

    rr = tapply( X=dd, INDEX=list(xx, yy), FUN = function(w) {length( which(is.finite(w) )}, simplify=TRUE )
    rr = as.data.frame( as.table (rr) )
    rr[,1] = as.numeric(as.character( rr[,1] ))
    rr[,2] = as.numeric(as.character( rr[,2] ))
    rr = rr[ which( is.finite( rr[,4] )) ,]
    names(rr) =c("x", "y", "n")
  }

  rrsum = sum(rr$n, na.rm=TRUE)
  ntoremove = rrsum - ntarget
  invcount = 1/rr$n * ntarget / rrsum  # proportion to remove to make each cell equal in weight

  keep = NULL
  for (o in 1:nrow(rr)) {
    if (!is.null(z)) {
      oo = which( xx == rr[o,1] & yy == rr[o,2] & tt == rr[o,2] )
    } else {
      oo = which( xx == rr[o,1] & yy == rr[o,2] )
    }
    noo = length(oo)
    if ( noo > 1) {
      okeep = floor(invcount[o] * noo)
      if (okeep > 1) {
        ss = oo[ .Internal( sample( noo, okeep, replace=FALSE, prob=NULL)) ]
        keep = c(keep, ss )
      }
    } else {
      keep = c(keep, noo )
    }
  }
  return( tokeep )
}
