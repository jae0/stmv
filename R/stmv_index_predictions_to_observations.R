
  stmv_index_predictions_to_observations = function(p) {

    Ploc = stmv_attach( p$storage_backend, p$ptr$Ploc )
    Yloc = stmv_attach( p$storage_backend, p$ptr$Yloc )
    coo = Yloc[]
    coo[,1] = trunc( coo[,1] / p$pres  ) * p$pres
    coo[,2] = trunc( coo[,2] / p$pres  ) * p$pres

    iYP = match(
      array_map( "xy->1", coo[], gridparams=p$gridparams ),
      array_map( "xy->1", Ploc[], gridparams=p$gridparams )
    )
    coo = NULL
    if ( p$dimensionality =="space" ) {
      # nothing to do
    }
    if ( p$dimensionality =="space-time" ) {
      Ytime = stmv_attach( p$storage_backend, p$ptr$Ytime )
      iYP = cbind(iYP, match( trunc(Ytime[]), p$yrs ) )
    }
    if ( p$dimensionality =="space-time-cyclic" ) {
      Ytime = stmv_attach( p$storage_backend, p$ptr$Ytime )
      yrs = trunc(Ytime[])
      dyear = Ytime[] - yrs
      dyear_breaks = c(p$dyears, p$dyears[length(p$dyears)]+ diff(p$dyears)[1] )
      dyear_index = as.numeric( cut( dyear, breaks=dyear_breaks, include.lowest=TRUE, ordered_result=TRUE, right=FALSE ) )
      iYP = cbind(iYP, match( yrs, p$yrs ), dyear_index )
      yrs = dyear = dyear_breaks = dyear_index = NULL
    }
    return (iYP)
  }
