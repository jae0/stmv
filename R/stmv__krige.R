
stmv__krige = function( p=NULL, dat=NULL, pa=NULL, nu=NULL, phi=NULL, varObs=NULL, varSpatial=NULL, variablelist=FALSE, ... ) {
  #\\ this is the core engine of stmv .. localised space (no-time) modelling interpolation
  # \ as a 2D gaussian process (basically, simple krigimg or TPS -- time is treated as being independent)
  #\\ note: time is not being modelled and treated independently
  #\\      .. you had better have enough data in each time slice ..  essentially this is kriging
  if (variablelist)  return( c() )

  sdTotal = sd(dat[[ p$stmv_variables$Y ]], na.rm=T)

  dat$mean = NA

  pa = data.table(pa)

  pa$mean = NA
  pa$sd = sqrt(varSpatial)  # leave as this as sd estimation is too expensive

  vns = p$stmv_variables$LOCS

  for ( ti in 1:p$nt ) {

    if ( exists("TIME", p$stmv_variables) ) {
      xi = which( dat[[ p$stmv_variables$TIME ]] == p$prediction_ts[ti] )
      pa_i = which( pa[[ p$stmv_variables$TIME ]] == p$prediction_ts[ti] )
    } else {
      xi = 1:nrow(dat) # all data as p$nt==1
      pa_i = 1:nrow(pa)
    }

    fspmodel = try( Krig( dat[xi, ..vns], dat[xi] [[p$stmv_variables$Y]],
      sigma2=varObs, rho=varSpatial , cov.function="stationary.cov",
      Covariance="Matern", range=phi, smoothness=nu) )
    if (inherits(fspmodel, "try-error") )  next()
    dat$mean[xi] = fspmodel$fitted.values
    ss = try( lm( dat$mean[xi] ~ dat[xi] [[p$stmv_variables$Y ]], na.action=na.omit) )
    if ( inherits(ss, "try-error") ) next()
    rsquared = summary(ss)$r.squared
    if (rsquared < p$stmv_rsquared_threshold ) next()
    pa$mean[pa_i] = predict(fspmodel, x=pa[pa_i, ..vns] )
    # pa$sd[pa_i]   = predictSE(fspmodel, x=pa[pa_i, ..vns] ) # SE estimates are slooow

    if ( 0 ){
      # debugging plots
      ti = 1
      vnt = c( p$stmv_variables$LOCS, p$stmv_variables$Y)
      xi = which( dat[[ p$stmv_variables$TIME ]] == p$prediction_ts[ti] )
      mbas = MBA::mba.surf( dat[xi, ..vnt ], 300, 300, extend=TRUE)$xyz.est
      image(mbas)

      surface(fspmodel)
      fsp.p<- predictSurface(fspmodel, nx=500, ny=500 ) # finer grid
      surface(fsp.p, type="I")

      fsp.p2<- predictSurfaceSE(fspmodel)
      surface(fsp.p2, type="C")
    }

  }

  # plot(pred ~ z , dat)
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
  ss = try( lm( dat$mean ~ dat[,p$stmv_variables$Y], na.action=na.omit) )
  if ( inherits(ss, "try-error") ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$stmv_rsquared_threshold ) return(NULL)

  stmv_stats = list( sdTotal=sdTotal, rsquared=rsquared, ndata=nrow(dat) ) # must be same order as p$statsvars
  return( list( predictions=pa, stmv_stats=stmv_stats ) )
}
