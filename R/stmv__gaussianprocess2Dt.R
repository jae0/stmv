
stmv__gaussianprocess2Dt = function(p=NULL, dat=NULL, pa=NULL, localrange=NULL, ...  ) {
  #\\ this is the core engine of stmv .. localised space (no-time) modelling interpolation
  # \ as a 2D gaussian process (basically, simple krigimg or TPS -- time is treated as being independent)
  #\\ note: time is not being modelled and treated independently
  #\\      .. you had better have enough data in each time slice ..  essentially this is kriging

  if (!exists("fields.cov.function", p)) p$fields.cov.function="stationary.cov"
  if (!exists("fields.cov.args", p) & p$fields.Covariance=="Matern") {
    if (!exists("fields.nu", p)) p$fields.nu=0.5  # note: this is the smoothness or shape parameter (fix at 0.5 if not calculated or given -- exponential)
    p$fields.cov.args=list( Covariance=p$fields.Covariance, smoothness=p$fields.nu ) # this is exponential covariance
  }
  if (!exists("fields.Covariance", p)) p$fields.Covariance="Exponential" # note that "Rad.cov" is TPS

  sdTotal=sd(dat[[ p$stmv_variables$Y ]], na.rm=TRUE)

  pa = data.table(pa)

  dat$mean = NA
  pa$mean = NA
  pa$sd = NA

  phi.grid = p$phi.grid * localrange
  vns = p$stmv_variables$LOCS

  for ( ti in 1:p$nt ) {

    if ( exists("TIME", p$stmv_variables) ) {
      xi = which( dat[[ p$stmv_variables$TIME ]] == p$prediction_ts[ti] )
    } else {
      xi = 1:nrow(dat) # all data as p$nt==1
    }

    xy = dat[xi, ..vns]
    z = dat[xi] [[p$stmv_variables$Y]]

    fsp = try( spatialProcess(xy, z, cov.function=p$fields.cov.function, cov.args=p$fields.cov.args ,
      theta.grid=phi.grid, lambda.grid=p$lambda.grid, ngrid = 10, niter = 15, tol = 0.01,
      Distance = "rdist", nstep.cv = 50 ) )

    if (inherits(fsp, "try-error") )  next()
    if ( fsp$converge != 0 ) next()

    fspmodel <- try( Krig( xy, z, cov.function=p$fields.cov.function, cov.args=p$fields.cov.args,
      theta=fsp$pars["theta"], lambda=fsp$pars["lambda"] ) )
    if (inherits(fspmodel, "try-error") )  next()

    dat$mean[xi] = as.vector( predict(fspmodel, x=dat[xi, ..vns] ) )
    ss = lm( dat$mean[xi] ~ dat[xi][[ p$stmv_variables$Y ]], na.action=na.omit)
    if ( "try-error" %in% class( ss ) ) next()
    rsquared = summary(ss)$r.squared
    if (rsquared < p$stmv_rsquared_threshold ) next()

    if ( exists("TIME", p$stmv_variables) ) {
      pa_i = which( pa[, p$stmv_variables$TIME]==p$prediction_ts[ti])
    } else {
      pa_i = 1:nrow(pa)
    }

    pa$mean[pa_i] = predict(fspmodel, x=pa[pa_, ..vns] )
    pa$sd[pa_i]   = predictSE(fspmodel, x=pa[pa_, ..vns] )


    if ( 0 ){
      # debugging plots
      surface(fspmodel)
      fsp.p<- predictSurface(fspmodel, lambda=fsp$pars["lambda"], nx=200, ny=200, )
      surface(fsp.p, type="I")
      fsp.p2<- predictSurfaceSE(fspmodel)
      surface(fsp.p, type="C")
    }

  }

  # plot(pred ~ z , dat)
  ss = lm( dat$mean ~ dat[[p$stmv_variables$Y]], na.action=na.omit)
  if ( inherits(ss, "try-error") ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$stmv_rsquared_threshold ) return(NULL)

  # TODO:: add some more stats: eg. range estimates, nugget/sill, etc..

  stmv_stats = list( sdTotal=sdTotal, rsquared=rsquared, ndata=nrow(dat) )

  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, stmv_stats=stmv_stats ) )
}

