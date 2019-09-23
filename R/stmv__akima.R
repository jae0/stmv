
stmv__akima = function( p=NULL,  dat=NULL, pa=NULL,  variablelist=FALSE, ...  ) {
  #\\ this is the core engine of stmv .. localised space (no-time) modelling interpolation
  #\\ note: time is not being modelled and treated independently
  #\\      .. you had better have enough data in each time slice ..  essentially this is cubic b-splines interpolation

  library(akima)  ### slow

  if (variablelist)  return( c() )

  x_r = range(dat[,p$variables$LOCS[1]])
  x_c = range(dat[,p$variables$LOCS[2]])

  nr = round( diff(x_r)/p$pres +1 )
  nc = round( diff(x_c)/p$pres +1 )

  xo = seq(x_r[1], x_r[2], length.out = nr )
  yo = seq(x_c[1], x_c[2], length.out = nc )

  origin=c(x_r[1], x_c[1])
  res=c(p$pres, p$pres)

  dat$mean = NA
  pa$mean = NA
  pa$sd = NA

  sdTotal = sd(dat[,p$variable$Y], na.rm=T)

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

    X = interp( x=dat[xi, p$variables$LOCS[1] ], y=dat[xi, p$variables$LOCS[2]], z=dat[xi, p$variables$Y],
      xo=xo, yo=yo, linear=TRUE, extrap=FLASE, duplicate="mean" )$z  # cannot extrapolate with linear
    # pa$sd[pa_i] = NA  ## fix as NA

    rY = range( dat[ xi, p$variables$Y ], na.rm=TRUE )
    lb = which( X < rY[1] )
    if (length(lb) > 0) X[lb] = rY[1]
    lb = NULL
    ub = which( X > rY[2] )
    if (length(ub) > 0) X[ub] = rY[2]
    ub = NULL

    X_i = array_map( "xy->2", coords=pa[pa_i, p$variables$LOCS], origin=origin, res=res )
    tokeep = which( X_i[,1] >= 1 & X_i[,2] >= 1  & X_i[,1] <= nr & X_i[,2] <= nc )
    if (length(tokeep) < 1) next()
    X_i = X_i[tokeep,]

    pa$mean[pa_i[tokeep]] = X[X_i]
    pa$sd[pa_i[tokeep]] = sd(dat[ xi, p$variable$Y ], na.rm=TRUE)  # timeslice guess

    X = X_i = NULL

    dat[ xi, p$variable$LOCS ] = round( dat[ xi, p$variable$LOCS ] / p$pres  ) * p$pres
    iYP = match(
      stmv::array_map( "xy->1", dat[ xi, p$variable$LOCS ], gridparams=p$gridparams ),
      stmv::array_map( "xy->1", pa[ pa_i , p$variable$LOCS ], gridparams=p$gridparams )
    )
    dat$mean[xi] = pa$mean[pa_i][iYP]
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
