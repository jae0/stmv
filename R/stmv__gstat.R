
stmv__gstat = function( p=NULL, dat=NULL, pa=NULL, nu=NULL, phi=NULL, varObs=NULL, varSpatial=NULL, variablelist=FALSE, ...  ) {
  #\\ this is the core engine of stmv .. localised space (no-time) modelling interpolation
  #\\ note: time is not being modelled and treated independently
  #\\      .. you had better have enough data in each time slice ..  essentially this is kriging
  if (variablelist)  return( c() )

  if (!exists( "stmv_gstat_formula", p)) p$stmv_gstat_formula = formula( paste( p$stmv_variables$Y, "~ 1 "))

  sdTotal = sd(dat[,p$stmv_variables$Y], na.rm=T)

  approx_range = phi*sqrt( 8*nu)

  dat$mean = NA
  pa$mean = NA
  pa$sd = sdTotal  # leave as this as sd estimation is too expensive

  for ( ti in 1:p$nt ) {

    if ( exists("TIME", p$stmv_variables) ) {
      xi = which( dat[ , p$stmv_variables$TIME ] == p$prediction_ts[ti] )
      pa_i = which( pa[, p$stmv_variables$TIME]==p$prediction_ts[ti])
    } else {
      xi = 1:nrow(dat) # all data as p$nt==1
      pa_i = 1:nrow(pa)
    }
    xy = dat[xi, c( p$stmv_variables$LOCS, p$stmv_variables$Y) ]

    vMod0 = vgm(psill=varSpatial, model="Mat", range=phi, nugget=varObs, kappa=nu ) # starting model parameters
    gs = gstat(id = "hmk", formula=p$stmv_gstat_formula, locations=~plon+plat, data=xy[xi,], maxdist=approx_range, nmin=p$stmv_nmin, nmax=p$stmv_nmax, force=TRUE, model=vMod0 )
    # this step adds a lot of time ..
    preds = predict(gs, newdata=xy[xi,] )
    dat$mean[xi] = as.vector( preds[,1] )
    ss = try( lm( dat$mean[xi] ~ dat[xi,p$stmv_variables$Y], na.action=na.omit))
    if ( inherits(ss, "try-error") ) next()
    rsquared = summary(ss)$r.squared
    if (rsquared < p$stmv_rsquared_threshold ) next()
    gsp = predict(gs, newdata=pa[pa_i,] ) # slow for large n
    pa$mean[pa_i] = as.vector(gsp[,1] )
    pa$sd[pa_i]   = as.vector(gsp[,2] )

  }

  # plot(pred ~ z , dat)
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
  ss = try( lm( dat$mean ~ dat[,p$stmv_variables$Y], na.action=na.omit))
  if ( inherits(ss, "try-error") ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$stmv_rsquared_threshold ) return(NULL)

  stmv_stats = list( sdTotal=sdTotal, rsquared=rsquared, ndata=nrow(dat) ) # must be same order as p$statsvars
  return( list( predictions=pa, stmv_stats=stmv_stats ) )
}

