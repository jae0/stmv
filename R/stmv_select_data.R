
stmv_select_data = function( p, Si, localrange ) {

  # obtain indices of data locations withing a given spatial range, optimally determined via variogram
  # faster to take a block .. but easy enough to take circles ...
  E = stmv_error_codes()
  Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )
  Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
  Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )
  Yi = stmv_attach( p$storage.backend, p$ptr$Yi )  # initial indices of good data
  if ( exists("TIME", p$variables)) Ytime = stmv_attach( p$storage.backend, p$ptr$Ytime )

  U = which(
    (abs( Sloc[Si,1] - Yloc[Yi[],1] ) <= localrange ) &
    (abs( Sloc[Si,2] - Yloc[Yi[],2] ) <= localrange )
  )
  ndata = length(U)

  if (ndata < p$stmv_nmin) {

    Sflag[Si] = E[["insufficient_data"]]

  } else if (ndata > p$stmv_nmax) {

    # try to trim
    if ( exists("TIME", p$variables)) {
      iU = stmv_discretize_coordinates( coo=cbind(Yloc[U,], Ytime[U]), ntarget=p$stmv_nmax, minresolution=p$minresolution, method="thin" )
    } else {
      iU = stmv_discretize_coordinates( coo=Yloc[U,], ntarget=p$stmv_nmax, minresolution=p$minresolution, method="thin" )
    }
    ndata = length(iU)
    # print(ndata)

    ntarget = floor( (p$stmv_nmax+p$stmv_nmin)/2 )

    if (ndata < p$stmv_nmin) {
      # return some data
      uu = setdiff( 1:length(U), iU)
      nuu = length(uu)
      iMore = uu[ .Internal( sample( nuu, {ntarget - ndata}, replace=FALSE, prob=NULL)) ]
      U = U[c(iU, iMore)]
      ndata = ntarget
      iMore = nuu = uu = NULL

    } else if (ndata > p$stmv_nmax) {

      # force via a random subsample
      U = U[iU]
      U = U[ .Internal( sample( length(U), ntarget, replace=FALSE, prob=NULL)) ] # simple random
      ndata = ntarget

    } else {

      U = U[iU]

    }

  } else  if (ndata <= p$stmv_nmax & ndata >= p$stmv_nmin) {
    # all good .. nothing to do
  }
  iU = NULL

  return(Yi[U])

}
