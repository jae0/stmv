

stmv_interpolate = function( ip=NULL, p, debugging=FALSE, runoption="default", ... ) {
  #\\ core function to interpolate (model and predict) in parallel

  if (0) {
    # for debugging  runs ..
    currentstatus = stmv_statistics_status( p=p )
    p = parallel_run( p=p, runindex=list( locs=sample( currentstatus$todo )) )
    ip = 1:p$nruns
    debugging=TRUE
    runoption="default"
  }

  # ---------------------
  # deal with additional passed parameters
  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast=TRUE))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable

  if (exists( "libs", p)) suppressMessages( RLibrary( p$libs ) )

  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns


  #---------------------
  # data for modelling
  S = stmv_attach( p$storage.backend, p$ptr$S )
  Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )
  E = stmv_error_codes()

  Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
  Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )

  Y = stmv_attach( p$storage.backend, p$ptr$Y )

  if (p$nloccov > 0) Ycov = stmv_attach( p$storage.backend, p$ptr$Ycov )
  if ( exists("TIME", p$variables) ) Ytime = stmv_attach( p$storage.backend, p$ptr$Ytime )

  # misc intermediate calcs to be done outside of parallel loops

  # pre-calculate indices and dim for data to use inside the loop
  dat_names = unique( c(  p$variables$Y, p$variable$LOCS, p$variables$local_all,  "weights") )  # excludes p$variables$TIME
  if (p$stmv_local_modelengine %in% c("fft", "tps") ) {
    if ( exists("TIME", p$variables)) {
      dat_names = c(dat_names, p$variables$TIME)
    }
  }
  # unless it is an explicit covariate and not a seasonal component there is no need for it
  # .. prediction grids create these from a time grid on the fly

  dat_nc = length( dat_names )
  iY = which(dat_names== p$variables$Y)
  ilocs = which( dat_names %in% p$variable$LOCS )
  # iwei = which( dat_names %in% "weights" )

  if (p$nloccov > 0) {
    icov = which( dat_names %in% p$variables$local_cov )
    icov_local = which( p$variables$COV %in% p$variables$local_cov )
  }
  if (exists("TIME", p$variables)) {
    ti_cov = setdiff(p$variables$local_all, c(p$variables$Y, p$variables$LOCS, p$variables$local_cov ) )
    itime_cov = which(dat_names %in% ti_cov)
  }

  local_fn = ifelse (p$stmv_local_modelengine=="userdefined", p$stmv_local_modelengine_userdefined, stmv_interpolation_function( p$stmv_local_modelengine ) )

  nip = length(ip)
  if (nip < 100) {
    nlogs = 3
  } else {
    nlogs = p$nlogs
  }
  logpoints  = ip[ floor( seq( from=10, to=(nip-10), length.out=nlogs ) ) ]

  if (debugging) {
    nsavepoints = 3
    savepoints = sample(logpoints, nsavepoints)
  }


  if (runoption=="default") local_corel = p$stmv_range_correlation
  if (runoption=="boostdata") local_corel = p$stmv_range_correlation_boostdata

## drange = max( min( max(p$stmv_distance_scale )),  min(p$pres, p$stmv_distance_scale ))
      # global estimates
  global_nu = median( S[, match("nu", p$statsvars)], na.rm=TRUE )

  if ( global_nu < 0.01 | global_nu > 5) global_nu = 0.5
  global_phi = median( S[, match("phi", p$statsvars)], na.rm=TRUE )
  global_range = matern_phi2distance( phi=global_phi, nu=global_nu, cor=local_corel )
  distance_limits = range( c(p$pres*3,  p$stmv_distance_scale ) )   # for range estimate
  if ( global_range < distance_limits[1] | global_range > distance_limits[2] ) global_range = distance_limits[2]

  i_ndata = match( "ndata", p$statsvars )
  i_rsquared = match("rsquared", p$statsvars )
  i_nu = match("nu",   p$statsvars)
  i_phi = match("phi",   p$statsvars)
  i_sdSpatial = match("sdSpatial",   p$statsvars)
  i_sdObs = match("sdObs",   p$statsvars)

# main loop over each output location in S (stats output locations)
  for ( iip in ip ) {

    if ( iip %in% logpoints )  currentstatus = stmv_logfile(p=p, flag= paste("Interpolation", runoption) )
    Si = p$runs[ iip, "locs" ]

    # print( paste("index =", iip, ";  Si = ", Si ) )
    if (debugging) print( paste("index =", iip, ";  Si = ", Si ) )
    if ( Sflag[Si] == E[["complete"]] ) next()

    nu    = S[Si, i_nu]
    phi   = S[Si, i_phi]

    ii = NULL
    ii = which(
      {abs( Sloc[Si,1] - Sloc[,1] ) <= distance_limits[2]} &
      {abs( Sloc[Si,2] - Sloc[,2] ) <= distance_limits[2]}
    )

    # nu checks
    if ( !is.finite(nu) ) {
      Sflag[Si] %in% E[["variogram_failure"]]
    } else {
      if ( (nu < 0.01 ) | (nu > 5) ) {
        Sflag[Si] %in% E[["variogram_range_limit"]]
      }
    }
    if ( !is.finite(nu) | (nu < 0.01) | (nu > 5) ) if (length(ii) > 0) nu =  median( S[ii, match("nu", p$statsvars)], na.rm=TRUE )
    if ( !is.finite(nu) | (nu < 0.01) | (nu > 5) ) nu = global_nu
    if ( !is.finite(nu) | (nu < 0.01) | (nu > 5) ) nu = 0.5

    # phi checks
    phi_lim = distance_limits / sqrt(8*nu)  # inla-approximation
    if ( !is.finite(phi) ) {
      Sflag[Si] %in% E[["variogram_failure"]]
    } else {
      if ( (phi < phi_lim[1]) | (phi > phi_lim[2]) ) {
        Sflag[Si] %in% E[["variogram_range_limit"]]
      }
    }
    if ( !is.finite(phi) | (phi < phi_lim[1]) | (phi > phi_lim[2]) )  if (length(ii) > 0)  phi =  median( S[ii, match("phi", p$statsvars)], na.rm=TRUE )
    if ( !is.finite(phi) | (phi < phi_lim[1]) | (phi > phi_lim[2]) )  phi = global_phi
    if ( !is.finite(phi) | (phi < phi_lim[1]) | (phi > phi_lim[2]) )  phi = phi_lim[2]


    # range checks
    localrange =  matern_phi2distance( phi=phi, nu=nu, cor=local_corel )
    if ( !is.finite(localrange) )  {
      Sflag[Si] %in% E[["variogram_failure"]]
    } else {
      if ( (localrange < distance_limits[1] ) | (localrange > distance_limits[2]) )  {
        Sflag[Si] %in% E[["variogram_range_limit"]]
      }
    }
    if ( !is.finite(localrange) | (localrange < distance_limits[1] ) | (localrange > distance_limits[2]) ) if (length(ii) > 0) localrange = median( S[ii, match("range", p$statsvars)], na.rm=TRUE )
    if ( !is.finite(localrange) | (localrange < distance_limits[1] ) | (localrange > distance_limits[2]) ) localrange = global_range
    if ( !is.finite(localrange) | (localrange < distance_limits[1] ) | (localrange > distance_limits[2]) ) localrange = distance_limits[2]

    ii = NULL

    if ( Sflag[Si] != E[["todo"]] ) {
      if (exists("stmv_rangecheck", p)) {
        if (p$stmv_rangecheck=="paranoid") {
          if ( Sflag[Si] %in% c( E[["variogram_range_limit"]], E[["variogram_failure"]]) ) {
            if (debugging) message("Error: stmv_rangecheck paranoid")
            if (runoption == "default" )  next()
          }
        }
      }
    }

    # last check .. ndata can change due to random sampling
    yi = stmv_select_data( p=p, Si=Si, localrange=localrange )
    if (is.null( yi )) next()
    ndata = length(yi)
    S[Si, i_ndata] = ndata

    if (ndata < p$stmv_nmin) {
      yi=NULL;
      next()
    }

    # if here then there is something to do
    # NOTE:: yi are the indices of locally useful data
    # p$stmv_distance_prediction determines the data entering into local model construction

    # prep dependent data
    # reconstruct data for modelling (dat)
    dat = matrix( 1, nrow=ndata, ncol=dat_nc )
    dat[,iY] = Y[yi] # these are residuals if there is a global model
    # add a small error term to prevent some errors when duplicate locations exist; localrange offsets to positive values
    dat[,ilocs] = Yloc[yi,] + localrange * runif(2*ndata, -1e-6, 1e-6)

    if (p$nloccov > 0) dat[,icov] = Ycov[yi, icov_local] # no need for other dim checks as this is user provided
    if (exists("TIME", p$variables)) {
      dat[, itime_cov] = as.matrix(stmv_timecovars( vars=ti_cov, ti=Ytime[yi,] ) )
      dat[, p$variables$TIME] = Ytime[yi,] # not sure if this is needed ?...
    }
    dat = as.data.frame(dat)
    names(dat) = dat_names

    # remember that these are crude mean/discretized estimates
    if (debugging) {
      dev.new()
      # check data and statistical locations
      plot( Sloc[,], pch=20, cex=0.5, col="gray")
      points( Yloc[,], pch=20, cex=0.2, col="green")
      points( Yloc[yi,], pch=20, cex=1, col="yellow" )
      points( Sloc[Si,2] ~ Sloc[Si,1], pch=20, cex=5, col="blue" )
    }

    # construct data (including covariates) for prediction locations (pa)
    pa = try( stmv_predictionarea( p=p, sloc=Sloc[Si,], windowsize.half=p$windowsize.half ) )
    if (is.null(pa)) {
      Sflag[Si] = E[["prediction_area"]]
      if (debugging) message( Si )
      if (debugging) message("Error: prediction grid ... null .. this should not happen")
      pa = NULL
      next()
    }
    if ( inherits(pa, "try-error") ) {
      pa = NULL
      Sflag[Si] = E[["prediction_area"]]
      if (debugging) message("Error: prediction grid ... try-error .. this should not happen")
      next()
    }

    if (debugging) {

      dev.new()
      # check that position indices are working properly
      Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
      Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
      Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )
      plot( Yloc[yi,2] ~ Yloc[yi,1], col="red", pch=".",
        ylim=range(c(Yloc[yi,2], Sloc[Si,2], Ploc[pa$i,2]) ),
        xlim=range(c(Yloc[yi,1], Sloc[Si,1], Ploc[pa$i,1]) ) ) # all data
      points( Yloc[yi,2] ~ Yloc[yi,1], col="green" )  # with covars and no other data issues
      points( Sloc[Si,2] ~ Sloc[Si,1], col="blue" ) # statistical locations
      # statistical output locations
      grids= spatial_grid(p, DS="planar.coords" )
      points( grids$plat[floor( (Sloc[Si,2]-p$origin[2])/p$pres) + 1]
            ~ grids$plon[floor( (Sloc[Si,1]-p$origin[1])/p$pres) + 1] , col="purple", pch=25, cex=5 )
      points( Ploc[pa$i,2] ~ Ploc[ pa$i, 1] , col="black", pch=6, cex=0.7 ) # check on pa$i indexing -- prediction locations
    }

    yi = NULL


    # model and prediction .. outputs are in scale of the link (and not response)
    # the following permits user-defined models (might want to use compiler::cmpfun )

    res =NULL
    res = try(
      local_fn (
        p=p,
        dat=dat,
        pa=pa,
        nu=nu,
        phi=phi,
        varObs = S[Si, i_sdObs]^2,
        varSpatial = S[Si, i_sdSpatial]^2,
        sloc = Sloc[Si,],
        distance = localrange
      )
    )

    dat =  NULL
    pa  =  NULL

    if ( is.null(res)) {
      Sflag[Si] = E[["local_model_error"]]   # modelling / prediction did not complete properly
      if (debugging) message("Error: local model error")
      next()
    }

    if ( inherits(res, "try-error") ) {
      Sflag[Si] =  E[["local_model_error"]]   # modelling / prediction did not complete properly
      res = NULL
      if (debugging) message("Error: local model error")
      next()
    }

    if (!exists("predictions", res)) {
      Sflag[Si] =  E[["prediction_error"]]   # modelling / prediction did not complete properly
      res = NULL
      if (debugging) message("Error: prediction error")
      next()
    }

    if (!exists("mean", res$predictions)) {
      Sflag[Si] =  E[["prediction_error"]]   # modelling / prediction did not complete properly
      res = NULL
      if (debugging) message("Error: prediction error")
      next()
    }

    if (length(which( is.finite(res$predictions$mean ))) < 5) {
      Sflag[Si] =  E[["prediction_error"]]   # modelling / prediction did not complete properly
      res = NULL
      if (debugging) {
        message("Error: prediction error")
        browser()
      }
      next()  # looks to be a faulty solution
    }


    # update to rsquared and "range" in stats
    if ( exists("rsquared", res$stmv_stats ) ) {
      if (is.finite(res$stmv_stats$rsquared)) S[Si, i_rsquared] = res$stmv_stats$rsquared
    }
    res$stmv_stats = NULL # reduce memory usage


    sf = try( stmv_predictions_update(p=p, preds=res$predictions ) )
    res = NULL

    if ( is.null(sf) ) {
      Sflag[Si] = E[["prediction_update_error"]]
      sf = NULL
      if (debugging) message("Error: prediction update error .. is null")
      next()
    }
    if ( inherits(sf, "try-error") ) {
      Sflag[Si] = E[["prediction_update_error"]]
      sf = NULL
      if (debugging) message("Error: prediction update error .. try-error")
      next()
    }
    if ( sf=="error" ) {
      Sflag[Si] = E[["prediction_update_error"]]
      sf = NULL
      if (debugging) message("Error: prediction update error .. general")
      next()
    }


    # ----------------------
    # do last. it is an indicator of completion of all tasks
    # restarts would be broken otherwise
    Sflag[Si] = E[["complete"]]  # mark as complete without issues

  }  # end for loop


  return(NULL)



    if (0) {

      require(MBA)
      require(fields)

      # kriged
      fit = Krig( dat[, c("plon", "plat")], dat$z, Covariance="Matern", theta=phi, smoothness=0.5)
      dev.new()
      op = predict(fit)
      tst = cbind( dat[, c("plon", "plat")], op )
      mba.int <- mba.surf( tst, 300, 300, extend=TRUE)$xyz.est
      image(mba.int, xaxs="r", yaxs="r")

      # raw data + mba
      dev.new()
      tst = cbind(  dat[, c("plon", "plat")], dat$z )
      mba.int <- mba.surf( tst, 300, 300, extend=TRUE)$xyz.est
      image(mba.int, xaxs="r", yaxs="r")

      # default
      dev.new()
      tst = cbind( res$predictions$plon,  res$predictions$plat,  res$predictions$mean )
      mba.int <- mba.surf( tst, p$windowsize.half *2, p$windowsize.half *2, extend=TRUE)$xyz.est
      image(mba.int, xaxs="r", yaxs="r")

      # kernel-based
      tst = as.image( Z=dat$z, x=dat[, c("plon", "plat")], nx=300, ny=300, na.rm=TRUE)
      out = fields::image.smooth( tst, theta=phi/300, xwidth=p$pres, ywidth=p$pres )
      dev.new(); image(out)

      print( str(res) )

      lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
      lattice::levelplot( mean ~ plon + plat, data=res$predictions, col.regions=heat.colors(100), scale=list(draw=TRUE) , aspect="iso" )
      for( i in sort(unique(res$predictions[,p$variables$TIME])))  print(lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==i,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )

      dev.new()
      plot(  dat[,iY] ~ dat$yr, col="red"  )
      points( mean~tiyr, res$predictions, pch=20, col="gray", cex=0.5 )

    }


}
