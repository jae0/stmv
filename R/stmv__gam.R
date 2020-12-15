
stmv__gam = function( p=NULL, dat=NULL, pa=NULL, ... ) {
  #\\ this is the core engine of stmv .. localised space-time modelling interpolation and prediction
  #\\ simple GAM with spatial weights (inverse distance squared) and ts harmonics
  #\\ family is gaussian("identity") as we are operating upon the link scale by this point

  sdTotal=sd(dat[[p$stmv_variables$Y]], na.rm=T)
  fit = NULL
  if ( exists("stmv_local_model_distanceweighted", p) ) {
    if (p$stmv_local_model_distanceweighted) {
       fit = try( gam( p$stmv_local_modelformula, data=dat, na.action="na.omit", weights=weights ) )
    }
  }

  if (is.null(fit)) fit = try( gam( p$stmv_local_modelformula, data=dat, na.action="na.omit" ) )

  if (is.null(fit)) return(NULL)

  if ( inherits( fit, "try-error") ) return( NULL )

  ss = summary( fit)
#  if (ss$r.sq < p$stmv_rsquared_threshold ) return(NULL)  # smooth/flat surfaces are ok ..

  out = try( predict(  fit, newdata=pa, type="response", se.fit=TRUE ) )

  if ( inherits(out, "try-error") ) return( NULL )

  pa$mean = as.vector(out$fit)
  pa$sd = as.vector(out$se.fit) # this is correct: se.fit== stdev of the mean fit: eg:  https://stat.ethz.ch/pipermail/r-help/2005-July/075856.html

  stmv_stats = list( sdTotal=sdTotal, rsquared=ss$r.sq, ndata=ss$n )

  # lattice::levelplot( mean ~ plon + plat, data=pa[pa$tiyr==min(pa$tiyr),], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, stmv_stats=stmv_stats ) )
}
