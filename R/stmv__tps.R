
stmv__tps = function( p=NULL, dat=NULL, pa=NULL, lambda=NULL, variablelist=FALSE, ...  ) {
  #\\ this is the core engine of stmv .. localised space (no-time) modelling interpolation
  # \ as a 2D gaussian process (basically, simple krigimg or TPS -- time is treated as being independent)
  #\\ note: time is not being modelled and treated independently
  #\\      .. you had better have enough data in each time slice ..
  if (variablelist)  return( c() )

  sdTotal = sd(dat[,p$variable$Y], na.rm=T)

  dat$mean = NA
  pa$mean = NA
  pa$sd = NA

  for ( ti in 1:p$nt ) {

    if ( exists("TIME", p$variables) ) {
      xi   = which( dat[ , p$variables$TIME] == p$prediction_ts[ti] )
      pa_i = which( pa[, p$variables$TIME] == p$prediction_ts[ti] )
      if (length(xi) < 5 ) {
        # print( ti)
        next()
      }
    } else {
      xi   = 1:nrow(dat) # all data as p$nt==1
      pa_i = 1:nrow(pa)
    }

    ftpsmodel = try( Tps(x=dat[xi, p$variables$LOCS], Y=dat[xi, p$variables$Y], lambda=lambda ) )
    if (inherits(ftpsmodel, "try-error") )  next()
    dat$mean[xi] = ftpsmodel$fitted.values
    ss = lm( dat$mean[xi] ~ dat[xi,p$variables$Y], na.action=na.omit)
    if ( "try-error" %in% class( ss ) ) next()
    rsquared = summary(ss)$r.squared
    if (rsquared < p$stmv_rsquared_threshold ) next()
    pa$mean[pa_i] = predict(ftpsmodel, x=pa[pa_i, p$variables$LOCS] )

    pa$sd[pa_i]   = predictSE(ftpsmodel, x=pa[pa_i, p$variables$LOCS] ) # SE estimates are slooow
    # pa$sd[pa_i] = sd( dat[xi,p$variable$Y], na.rm=T)   ## fix as NA

    if ( 0 ){
      # debugging plots
      surface(ftpsmodel)
    }
  }

  # plot(pred ~ z , dat)
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
  ss = lm( dat$mean ~ dat[,p$variables$Y], na.action=na.omit)
  if ( "try-error" %in% class( ss ) ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$stmv_rsquared_threshold ) return(NULL)

  stmv_stats = list( sdTotal=sdTotal, rsquared=rsquared, ndata=nrow(dat) ) # must be same order as p$statsvars
  return( list( predictions=pa, stmv_stats=stmv_stats ) )
}
