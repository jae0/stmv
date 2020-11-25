
stmv__kernel = function( p=NULL, dat=NULL, pa=NULL, phi=NULL, nu=NULL, varObs=NULL, varSpatial=NULL,  variablelist=FALSE, ...  ) {
  #\\ this is the core engine of stmv .. localised space (no-time) modelling interpolation
  #\\ note: time is not being modelled and treated independently
  #\\      .. you had better have enough data in each time slice ..  essentially this is cubic b-splines interpolation

  library(fields)

  if (variablelist)  return( c() )

  vns = p$stmv_variables$LOCS
  pa = data.table(pa)

  sdTotal = sd(dat[[ p$stmv_variables$Y ]], na.rm=T)

  if ( grepl("fast_predictions", p$stmv_fft_filter)) {
    # predict only where required
    x_r = range(pa[[ vns[1] ]] )
    x_c = range(pa[[ vns[2] ]] )
  } else  {
    # predict on full data subset
    x_r = range(dat[[ vns[1] ]] )
    x_c = range(dat[[ vns[2] ]] )
  }


  nr = aegis_floor( diff(x_r)/p$pres + 1L )
  nc = aegis_floor( diff(x_c)/p$pres + 1L )

  dat$mean = NA
  pa$mean = NA
  pa$sd = sqrt(varSpatial)

  origin=c(x_r[1], x_c[1])
  res=c(p$pres, p$pres)

  theta = matern_phi2distance( phi=phi, nu=nu, cor=p$stmv_autocorrelation_localrange )

  wght = setup.image.smooth( nrow=nr, ncol=nc,  dx=p$pres, dy=p$pres, theta=theta, xwidth=p$pres, ywidth=p$pres)

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

    X = as.image(
      dat[xi] [[p$stmv_variables$Y]],
      ind=as.matrix(array_map( "xy->2", coords=dat[xi,..vns], origin=origin, res=res )),
      na.rm=TRUE,
      nx=nr,
      ny=nc
    )

    X = fields::image.smooth( X, dx=dx, dy=dy, wght )

    if ( grepl("fast_predictions", p$stmv_fft_filter)) {
      pa$mean[pa_i] = X
      # pa$sd[pa_i[tokeep]] = sd( dat[xi,p$stmv_variables$Y], na.rm=T)   ## fix
      X = NULL
    }  else  {
      X_i = array_map( "xy->2", coords=pa[pa_i, ..vns], origin=origin, res=res )
      tokeep = which( X_i[,1] >= 1 & X_i[,2] >= 1  & X_i[,1] <= nr & X_i[,2] <= nc )
      if (length(tokeep) < 1) next()
      X_i = X_i[tokeep,]
      pa$mean[pa_i[tokeep]] = X[X_i]
      # pa$sd[pa_i[tokeep]] = sd( dat[xi] [[p$stmv_variables$Y]], na.rm=T)   ## fix
      X = X_i = NULL
    }

    # dat[ xi, ..vns ] = aegis_floor( dat[ xi, ..vns ] / p$pres + 1L ) * p$pres
    iYP = match(
      array_map( "xy->1", dat[ xi, ..vns ], gridparams=p$gridparams ),
      array_map( "xy->1", pa[ pa_i , ..vns ], gridparams=p$gridparams )
    )
    dat$mean[xi] = pa$mean[pa_i][iYP]

  }

  # plot(pred ~ z , dat)
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
  ss = try( lm( dat$mean ~ dat[[ p$stmv_variables$Y ]] , na.action=na.omit) )
  if ( inherits(ss, "try-error") ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$stmv_rsquared_threshold ) return(NULL)

  stmv_stats = list( sdTotal=sdTotal, rsquared=rsquared, ndata=nrow(dat) ) # must be same order as p$statsvars
  return( list( predictions=pa, stmv_stats=stmv_stats ) )
}
