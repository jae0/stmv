
stmv__tps = function( p=NULL, dat=NULL, pa=NULL, lambda=NULL, variablelist=FALSE, ...  ) {
  #\\ this is the core engine of stmv .. localised space (no-time) modelling interpolation
  # \ as a 2D gaussian process (basically, simple krigimg or TPS -- time is treated as being independent)
  #\\ note: time is not being modelled and treated independently
  #\\      .. you had better have enough data in each time slice ..
  if (variablelist)  return( c() )

  vns = p$stmv_variables$LOCS

  sdTotal = sd(dat[[ p$stmv_variables$Y ]], na.rm=T)

  dat$mean = NA

  pa = data.table(pa)
  pa$mean = NA
  pa$sd = NA

  for ( ti in 1:p$nt ) {

    if ( exists("TIME", p$stmv_variables) ) {
      xi   = which( dat[[ p$stmv_variables$TIME]] == p$prediction_ts[ti] )
      pa_i = which( pa[[ p$stmv_variables$TIME]] == p$prediction_ts[ti] )
      if (length(xi) < 5 ) {
        # print( ti)
        next()
      }
    } else {
      xi   = 1:nrow(dat) # all data as p$nt==1
      pa_i = 1:nrow(pa)
    }

    ftpsmodel = try( Tps(x=dat[xi, ..vns], Y=dat[xi] [[p$stmv_variables$Y]], lambda=lambda ) )
    if (inherits(ftpsmodel, "try-error") )  next()
    dat$mean[xi] = ftpsmodel$fitted.values
    ss = try( lm( dat$mean[xi] ~ dat[xi] [[p$stmv_variables$Y]], na.action=na.omit) )
    if ( inherits(ss, "try-error") ) next()
    rsquared = summary(ss)$r.squared
    if (rsquared < p$stmv_rsquared_threshold ) next()
    pa$mean[pa_i] = predict(ftpsmodel, x=pa[pa_i, ..vns] )

    pa$sd[pa_i]   = predictSE(ftpsmodel, x=pa[pa_i, ..vns] ) # SE estimates are slooow
    # pa$sd[pa_i] = sd( [xi] [[p$stmv_variables$Y]], na.rm=T)   ## fix as NA

    if ( 0 ){
      # debugging plots
      surface(ftpsmodel)
    }
  }

  # plot(pred ~ z , dat)
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
  ss = try( lm( dat$mean ~ dat[[ p$stmv_variables$Y ]], na.action=na.omit) )
  if ( inherits(ss, "try-error") ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$stmv_rsquared_threshold ) return(NULL)

  stmv_stats = list( sdTotal=sdTotal, rsquared=rsquared, ndata=nrow(dat) ) # must be same order as p$statsvars
  return( list( predictions=pa, stmv_stats=stmv_stats ) )
}
