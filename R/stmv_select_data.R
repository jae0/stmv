
stmv_select_data = function( p, Si, localrange ) {

  # obtain indices of data locations withing a given spatial range, optimally determined via variogram
  # faster to take a block .. but easy enough to take circles ...
  E = stmv_error_codes()
  Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )
  Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
  Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )
  Yi = stmv_attach( p$storage.backend, p$ptr$Yi )  # initial indices of good data

  # space only, ignore time as it is modelled separately
  # U are indices of Yi
  U = which(
    (abs( Sloc[Si,1] - Yloc[Yi[],1] ) <= localrange ) &
    (abs( Sloc[Si,2] - Yloc[Yi[],2] ) <= localrange )
  )
  # Yuniq are flags indexed on Yi
  Yuniq = !duplicated(Yloc[Yi[U],])  # dups in space only ... leave time alone
  ndata = length( which(Yuniq) )

  if (ndata < p$stmv_nmin) {
    Sflag[Si] = E[["insufficient_data"]]
    return(NULL)
  } else if (ndata > p$stmv_nmax) {
    # try to trim
    iU = NULL
    if ( exists("TIME", p$variables)) {
      if ( grepl("stmv_variogram_resolve_time", p$stmv_fft_filter)) {
        Ytime = stmv_attach( p$storage.backend, p$ptr$Ytime )
        iU = stmv_discretize_coordinates( coo=cbind(Yloc[Yi[U],], Ytime[Yi[U]]), ntarget=floor(p$stmv_nmax*p$stmv_tmax), minresolution=p$minresolution, method="thin" )
      }
    }
    if (is.null(iU)) iU = stmv_discretize_coordinates( coo=Yloc[Yi[U],], ntarget=p$stmv_nmax, minresolution=p$minresolution, method="thin" )

    ndata = length( which(Yuniq[iU] ) )
    ntarget = p$stmv_nmax
    if (ndata < p$stmv_nmin) {
      # return what data we do have .. calling funcction determines what to do with it if anything
      uu = setdiff( 1:length(U), iU)
      iMore = uu[ .Internal( sample( length(uu), {ntarget - ndata}, replace=FALSE, prob=NULL) ) ]
      iU = unique(sort(c(iU, iMore)))
      iMore =  uu = NULL
      Sflag[Si] = E[["insufficient_data"]]
    } else if (ndata > p$stmv_nmax) {
      # force via a random subsample
      iU = iU[ .Internal( sample( length(iU), ntarget, replace=FALSE, prob=NULL)) ] # simple random
    } else {
      # nothing to do
    }
    ndata = length( which(Yuniq[iU] ) )
    U = U[iU]
    iU = NULL
  } else if (ndata <= p$stmv_nmax & ndata >= p$stmv_nmin) {
    # all good .. nothing to do
  }
  return( list( data_index=Yi[U], unique_spatial_locations=ndata ) )
}
