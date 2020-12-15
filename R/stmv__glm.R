
stmv__glm = function(p=NULL, dat=NULL, pa=NULL, ... ) {
  #\\ this is the core engine of stmv .. localised space-time modelling interpolation and prediction
  #\\ simple GAM with spatial weights (inverse distance squared) and ts harmonics
  #\\ operating upon link scale ... family for the local model is forced to be gaussian("identity")

  sdTotal=sd(dat[[p$stmv_variables$Y]], na.rm=T)

  fit  = NULL

  if ( exists("stmv_local_model_distanceweighted", p) ) {
    if (p$stmv_local_model_distanceweighted) {
      fit = try( glm( p$stmv_local_modelformula, data=dat, weights=Y_wgt  ) )
    } else {
      fit = try( glm( p$stmv_local_modelformula, data=dat ) )
    }
  } else {
      fit = try( glm( p$stmv_local_modelformula, data=dat ) )
  }

  if (is.null(fit)) return(NULL)

  if ( inherits(fit, "try-error") ) return( NULL )

  ss = summary(fit)
  rsq = 1 - (ss$deviance/ss$null.deviance)
  if ( rsq < p$stmv_rsquared_threshold ) return(NULL)

  out = try( predict( fit, newdata=pa, type="response", se.fit=TRUE ) )  # already on link scale

  if ( inherits(out, "try-error") ) return( NULL )

  pa$mean = as.vector(out$fit)
  pa$sd = as.vector(out$se.fit) # this is correct: se.fit== stdev of the mean fit: eg:  https://stat.ethz.ch/pipermail/r-help/2005-July/075856.html

  stmv_stats = list( sdTotal=sdTotal, rsquared=rsq, ndata=nrow(dat) ) # pseudo rsquared for logistic .. for poisson {1- logLik(mod) / logLik(mod_saturated)} might be better

  # lattice::levelplot( mean ~ plon + plat, data=pa[pa$tiyr==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, stmv_stats=stmv_stats ) )
}

