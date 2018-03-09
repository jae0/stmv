
stmv__gam = function( p=NULL, dat=NULL, pa=NULL, variablelist=FALSE, ... ) {
  #\\ this is the core engine of stmv .. localised space-time modelling interpolation and prediction
  #\\ simple GAM with spatial weights (inverse distance squared) and ts harmonics
  #\\ family is gaussian("identity") as we are operating upon the link scale by this point
  if (variablelist)  return( c() )

  sdTotal=sd(dat[,p$variable$Y], na.rm=T)

  if ( exists("stmv_local_model_distanceweighted", p) ) {
    if (p$stmv_local_model_distanceweighted) {
      hmod = try( gam( p$stmv_local_modelformula, data=dat, na.action="na.omit", weights=weights, optimizer=p$stmv_gam_optimizer) )
    } else {
      hmod = try( gam( p$stmv_local_modelformula, data=dat, na.action="na.omit", optimizer=p$stmv_gam_optimizer  ) )
    }
  } else {
      hmod = try( gam( p$stmv_local_modelformula, data=dat, na.action="na.omit", optimizer=c("outer", "bfgs")  ) )
  }


  if ( "try-error" %in% class(hmod) ) return( NULL )

  ss = summary(hmod)
  if (ss$r.sq < p$stmv_rsquared_threshold ) return(NULL)

  out = try( predict( hmod, newdata=pa, type="response", se.fit=T ) )

  if ( "try-error" %in% class( out ) ) return( NULL )

  pa$mean = as.vector(out$fit)
  pa$sd = as.vector(out$se.fit) # this is correct: se.fit== stdev of the mean fit: eg:  https://stat.ethz.ch/pipermail/r-help/2005-July/075856.html

  stmv_stats = list( sdTotal=sdTotal, rsquared=ss$r.sq, ndata=ss$n ) # must be same order as p$statsvars

  # lattice::levelplot( mean ~ plon + plat, data=pa[pa$tiyr==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, stmv_stats=stmv_stats ) )
}

