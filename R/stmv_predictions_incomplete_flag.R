stmv_predictions_incomplete_flag  = function( p ) {

  # statistics locations where estimations need to be redone
  E = stmv_error_codes()
  P = stmv_attach( p$storage.backend, p$ptr$P )
  if (ncol(P) == 1 ) {
    noP = which( !is.finite( P[]) )
  } else {
    noP = which( !is.finite( rowSums( P[])) )
  }
  uP = NULL
  if( length(noP)>0 ) {
    Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
    Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )

    Sloc_nplat = ceiling( diff( p$corners$plat) / p$stmv_distance_statsgrid)
    Sloc_nplon = ceiling( diff( p$corners$plon) / p$stmv_distance_statsgrid)

    Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
    uS = array_map( "2->1", floor( cbind(Sloc[,1]-p$origin[1], Sloc[,2]-p$origin[2])/p$stmv_distance_statsgrid)+1, c(Sloc_nplon, Sloc_nplat) )
    uP = array_map( "2->1", floor( cbind(Ploc[noP,1]-p$origin[1], Ploc[noP,2]-p$origin[2])/p$stmv_distance_statsgrid)+1, c(Sloc_nplon, Sloc_nplat) )
    inrange = which( (uP >= min(uS)) & (uP <= max(uS)) )
    if (length( inrange) > 0) uP = uP[inrange]
    uP = unique(uP)
    Sflag[uP] = E[["todo"]]  # force set to redo
  }

  return(uP)
}
