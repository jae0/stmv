
stmv__constant = function( p=NULL,  dat=NULL, pa=NULL, ...  ) {
  #\\ this is the core engine of stmv .. localised space (no-time) modelling interpolation
  #\\ note: time is not being modelled and treated independently
  #\\      .. you had better have enough data in each time slice ..  essentially this is cubic b-splines interpolation

  sdTotal = sd(dat[,p$stmv_variables$Y], na.rm=T)

  vns = p$stmv_variables$LOCS

  pa = data.table(pa)

  x_r = range(pa[[vns[1] ]])
  x_c = range(pa[[vns[2] ]])

  nr = trunc( diff(x_r)/p$pres +1 )
  nc = trunc( diff(x_c)/p$pres +1 )

  xo = seq(x_r[1], x_r[2], length.out = nr )
  yo = seq(x_c[1], x_c[2], length.out = nc )

  dat$mean = NA
  pa$mean = NA
  pa$sd = NA


  for ( ti in 1:p$nt ) {

    if ( exists("TIME", p$stmv_variables) ) {
      xi   = which( dat[[ p$stmv_variables$TIME ]] == p$prediction_ts[ti] )
      pa_i = which( pa[[ p$stmv_variables$TIME ]] == p$prediction_ts[ti] )
      if (length(xi) < 5 ) {
        # print( ti)
        next()
      }
    } else {
      xi   = 1:nrow(dat) # all data as p$nt==1
      pa_i = 1:nrow(pa)
    }

    pa$mean[pa_i] = mean( dat[xi][[  p$stmv_variables$Y ]], na.rm=TRUE )
    pa$sd[pa_i]   = sd  ( dat[xi][[  p$stmv_variables$Y ]], na.rm=TRUE )  # just a crude guess for each timeslice


    iYP = match(
      array_map( "xy->1", dat[ xi, ..vns ], gridparams=p$gridparams ),
      array_map( "xy->1", pa[ pa_i , ..vns ], gridparams=p$gridparams )
    )
    dat$mean[xi] = pa$mean[pa_i][iYP]

  }

  # plot(pred ~ z , dat)
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
  ss = lm( dat$mean ~ dat[[ p$stmv_variables$Y ]], na.action=na.omit)
  if ( inherits(ss, "try-error") ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$stmv_rsquared_threshold ) return(NULL)

  stmv_stats = list( sdTotal=sdTotal, rsquared=rsquared, ndata=nrow(dat) )
  return( list( predictions=pa, stmv_stats=stmv_stats ) )
}
