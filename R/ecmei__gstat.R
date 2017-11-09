
ecmei__gstat = function( p, dat, pa, nu, phi, varObs, varSpatial ) {
  #\\ this is the core engine of ecmei .. localised space (no-time) modelling interpolation 
  #\\ note: time is not being modelled and treated independently 
  #\\      .. you had better have enough data in each time slice ..  essentially this is kriging 
  
  if (!exists( "ecmei_gstat_formula", p)) p$ecmei_gstat_formula = formula( paste( p$variables$Y, "~ 1 ")) 

  sdTotal = sd(dat[,p$variable$Y], na.rm=T)

  approx_range = phi*sqrt( 8*nu)

  dat$mean = NA
  pa$mean = NA
  pa$sd = sdTotal  # leave as this as sd estimation is too expensive

  for ( ti in 1:p$nt ) {
    
    if ( exists("TIME", p$variables) ) {
      xi = which( dat[ , p$variables$TIME ] == p$prediction.ts[ti] )
      pa_i = which( pa[, p$variables$TIME]==p$prediction.ts[ti])
    } else {
      xi = 1:nrow(dat) # all data as p$nt==1
      pa_i = 1:nrow(pa)
    }
    xy = dat[xi, c( p$variables$LOCS, p$variables$Y) ]

    vMod0 = vgm(psill=varSpatial, model="Mat", range=phi, nugget=varObs, kappa=nu ) # starting model parameters
    gs = gstat(id = "hmk", formula=p$ecmei_gstat_formula, locations=~plon+plat, data=xy[xi,], maxdist=approx_range, nmin=p$n.min, nmax=p$n.max, force=TRUE, model=vMod0 )
    # this step adds a lot of time .. 
    preds = predict(gs, newdata=xy[xi,] )
    dat$mean[xi] = as.vector( preds[,1] )
    ss = lm( dat$mean[xi] ~ dat[xi,p$variables$Y], na.action=na.omit)
    if ( "try-error" %in% class( ss ) ) next()
    rsquared = summary(ss)$r.squared
    if (rsquared < p$ecmei_rsquared_threshold ) next()
    gsp = predict(gs, newdata=pa[pa_i,] ) # slow for large n
    pa$mean[pa_i] = as.vector(gsp[,1] )
    pa$sd[pa_i]   = as.vector(gsp[,2] )

  }

  # plot(pred ~ z , dat)
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
  ss = lm( dat$mean ~ dat[,p$variables$Y], na.action=na.omit)
  if ( "try-error" %in% class( ss ) ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$ecmei_rsquared_threshold ) return(NULL)

  ecmei_stats = list( sdTotal=sdTotal, rsquared=rsquared, ndata=nrow(dat) ) # must be same order as p$statsvars
  return( list( predictions=pa, ecmei_stats=ecmei_stats ) )  
}

