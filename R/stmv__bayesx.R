
stmv__bayesx = function( p=NULL, dat=NULL, pa=NULL, variablelist=FALSE, ... ) {
  #\\ this is the core engine of stmv .. localised space-time modelling interpolation and prediction .. using bayesx

  # EG: see: bayesx.term.options( bs="kr", method="REML" )
  #  logzinc ~  sx( x,y, nu=1.5, bs="kr")  # "kr" is perhaps overly smooth  ..  ie guassian process  .. kriging
  #  logzinc ~  sx( x,y, bs="te")  # more detail .. "te" is preferred
  if (variablelist)  return( c() )

  sdTotal=sd(dat[[ p$stmv_variables$Y ]], na.rm=T)

  if ( !exists( "stmv_local_model_bayesxmethod", p) ) p$stmv_local_model_bayesxmethod="MCMC"  # slightly more smoothing than the REML method

  hmod = try( bayesx( p$stmv_local_modelformula, data=dat, method=p$stmv_local_model_bayesxmethod,
                     family="gaussian" ) )

  if ( inherits(hmods, "try-error") ) return( NULL )

  px = predict(hmod)
  ss = summary(lm( px ~ dat[[ p$stmv_variables$Y ]], na.action="na.omit" ))
  if (ss$r.squared < p$stmv_rsquared_threshold ) return(NULL)

  out = try( predict( hmod, newdata=pa, type="response" ) )

  # plot( hmod, term = "sx(x,y)", image=TRUE)
  # summary(hmod)
  # lattice::levelplot( out ~ plon+plat, data=k, aspect="iso" )

  if ( inherits(out, "try-error")) return( NULL )

  pa$mean = as.vector( out )
  pa$sd = 1 # no option right now to estim posterior prediction errors .. may be possible with type="terms" but would be slow to simulate  and do not know how to do it yet .. please fix this ..
  varSpatial = hmod$smooth.hyp[,"Variance"]
  varObs = hmod$fixed.effects[1,"Std. Error"]
  nu = 1.5
  phi=1/hmod$smooth.hyp[,"Smooth Par."]

  stmv_stats = list( sdTotal=sdTotal, rsquared=ss$r.squared, ndata=nrow(dat),
    sdSpatial=sqrt(varSpatial), sdObs=sqrt(varObs), phi=phi, nu=nu )

  # lattice::levelplot( mean ~ plon + plat, data=pa[pa$tiyr==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, stmv_stats=stmv_stats ) )
}
