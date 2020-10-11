

stmv_interpolate_lattice = function( ip=NULL, p, debugging=FALSE, just_testing_variablelist=FALSE,  ... ) {
  #\\ core function to interpolate (model and predict) in parallel

  if (0) {
    # for debugging  runs ..
    currentstatus = stmv_statistics_status( p=p )
    p = parallel_run( p=p, runindex=list( locs=sample( currentstatus$todo )) )
    ip = 1:p$nruns
    debugging=TRUE
  }


  Sloc = stmv_attach( p$storage_backend, p$ptr$Sloc )
  Yloc = stmv_attach( p$storage_backend, p$ptr$Yloc )
  Y = stmv_attach( p$storage_backend, p$ptr$Y )

  if (p$nloccov > 0) Ycov = stmv_attach( p$storage_backend, p$ptr$Ycov )
  if ( exists("TIME", p$stmv_variables) ) Ytime = stmv_attach( p$storage_backend, p$ptr$Ytime )

  # misc intermediate calcs to be done outside of parallel loops

  # pre-calculate indices and dim for data to use inside the loop
  dat_names = unique( c(  p$stmv_variables$Y, p$stmv_variables$LOCS, p$stmv_variables$local_all,  "weights") )  # excludes p$stmv_variables$TIME
  # if (p$stmv_local_modelengine %in% c("fft", "tps", "twostep") ) {
    if ( exists("TIME", p$stmv_variables)) {
      dat_names = c(dat_names, p$stmv_variables$TIME)
    }
  # }
  # unless it is an explicit covariate and not a seasonal component there is no need for it
  # .. prediction grids create these from a time grid on the fly

  dat_nc = length( dat_names )

  iY = which(dat_names== p$stmv_variables$Y)
  ilocs = which( dat_names %in% p$stmv_variables$LOCS )

  # iwei = which( dat_names %in% "weights" )
  i_rsquared = match("rsquared", p$statsvars )
  i_localrange = match("localrange", p$statsvars )
  i_nu = match("nu",   p$statsvars)
  i_phi = match("phi",   p$statsvars)
  i_sdSpatial = match("sdSpatial",   p$statsvars)
  i_sdObs = match("sdObs",   p$statsvars)

  distance_limits = range( p$stmv_distance_scale )    # for range estimate

  # error (km) to add to locations to force solutions that are affected by duplicated locations
  dist_error = 1e-6
  if ( exists( "stmv_lowpass_nu", p) & exists( "stmv_lowpass_phi", p) ) {
    dist_error_target = matern_phi2distance( phi=p$stmv_lowpass_phi, nu=p$stmv_lowpass_nu, cor=p$stmv_autocorrelation_localrange )
    if (is.finite(dist_error_target)) dist_error = dist_error_target
  }

  local_fn = ifelse (p$stmv_local_modelengine=="userdefined", p$stmv_local_modelengine_userdefined, stmv_interpolation_function( p$stmv_local_modelengine ) )

  if (p$nloccov > 0) {
    icov = which( dat_names %in% p$stmv_variables$local_cov )
    icov_local = which( p$stmv_variables$COV %in% p$stmv_variables$local_cov )
  }
  if (exists("TIME", p$stmv_variables)) {
    ti_cov = setdiff(p$stmv_variables$local_all, c(p$stmv_variables$Y, p$stmv_variables$LOCS, p$stmv_variables$local_cov ) )
    itime_cov = which(dat_names %in% ti_cov)
  }



  if (just_testing_variablelist) {

    # purpose: to obtain param names from a model run
    # param names vary depending upon the model form
    # only possible to run a model to see waht gets exported ..
    # this is essentially a minimal form of stmv_interpolate  that stops once a solution with params are found



    message("testing a run of the model to check for output")

    p = parallel_run( p=p, runindex=list( locs=sample( stmv_statistics_status( p=p )$todo )) )
    ip = 1:100

    for ( iip in ip ) {
      Si = p$runs[ iip, "locs" ]
      sloc = Sloc[Si,]
      nu    = S[Si, i_nu]
      phi   = S[Si, i_phi]
      localrange = S[Si, i_localrange]   # attached to p$stmv_autocorrelation_localrange
      varObs = S[Si, i_sdObs]^2,
      varSpatial = S[Si, i_sdSpatial]^2,

      localrange_interpolation = ifelse( !exists("stmv_interpolation_basis_distance", p), p$stmv_distance_statsgrid *1.5, p$stmv_interpolation_basis_distance )  # force a simple solution

      data_subset = stmv_select_data( p=p, Si=Si, localrange=localrange_interpolation )
      if (is.null( data_subset )) next()

      unique_spatial_locations = data_subset$unique_spatial_locations
      ndata = length(data_subset$data_index)
      if (unique_spatial_locations < p$stmv_nmin) next()

      dat = matrix( 1, nrow=ndata, ncol=dat_nc )
      dat[,iY] = Y[data_subset$data_index] # these are residuals if there is a global model
        dat[,ilocs] = Yloc[data_subset$data_index,]

      dat[,ilocs] = dat[,ilocs] + localrange_interpolation * runif(2*ndata, -dist_error, dist_error)

      if (p$nloccov > 0) dat[,icov] = Ycov[data_subset$data_index, icov_local] # no need for other dim checks as this is user provided
      if (exists("TIME", p$stmv_variables)) {
        dat[, itime_cov] = as.matrix(stmv_timecovars( vars=ti_cov, ti=Ytime[data_subset$data_index,] ) )
        itt = which(dat_names==p$stmv_variables$TIME)
        dat[, itt ] = Ytime[data_subset$data_index,]
        # crude check of number of time slices
        n_time_slices = stmv_discretize_coordinates( coo=dat[, itt], ntarget=p$nt, minresolution=p$minresolution[3], method="thin"  )
        if ( length(n_time_slices) < p$stmv_tmin )  next()
      }
      dat = data.table(dat)
      names(dat) = dat_names
      dat_range = range( dat[,..iY], na.rm=TRUE )  # used later

      # construct prediction/output grid area ('pa')
      prediction_area = localrange_interpolation
      if (exists("stmv_distance_prediction_limits", p)) {
        prediction_area = min( max( localrange_interpolation, min(p$stmv_distance_prediction_limits) ), max(p$stmv_distance_prediction_limits), na.rm=TRUE )
      }

      windowsize.half =  aegis_floor( prediction_area / p$pres ) + 1L

      pa = try( stmv_predictionarea_lattice( p=p, sloc=sloc, windowsize.half=windowsize.half ) )

      if ( is.null(pa) ) next()
      if ( inherits(pa, "try-error") ) next()

      res =NULL
      res = try(
        local_fn ( p=p, dat=dat, pa=pa,
          nu=nu, phi=phi, localrange=localrange_interpolation, varObs=varObs, varSpatial=varSpatial, sloc=sloc
        )
      )

      if ( is.null(res))  next()
      if ( inherits(res, "try-error") ) next()
      if (!exists("predictions", res))  next()
      if (!exists("mean", res$predictions)) next()
      if (length(which( is.finite(res$predictions$mean ))) < 1) next()

      return( res )
    }
  }


  # ---------------------

  p = parameters_control(p, list(...), control="add") # add passed args to parameter list, priority to args

  if (exists( "libs", p)) suppressMessages( RLibrary( p$libs ) )

  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns


  E = stmv_error_codes()
  Sflag = stmv_attach( p$storage_backend, p$ptr$Sflag )
  S = stmv_attach( p$storage_backend, p$ptr$S )

  i_ndata = match( "ndata", p$statsvars )

  if (length(ip) < 100) {
    nlogs = length(ip) / 5
  } else {
    nlogs = ifelse( length(ip) > (p$nlogs*5), p$nlogs, length(ip) / 5  )
  }
  logpoints  =  sort( sample( ip, round( max(1, nlogs) ) ) )  # randomize

  # main loop over each output location in S (stats output locations)
  for ( iip in ip ) {

    if ( iip %in% logpoints )  slog = stmv_logfile(p=p, flag= paste("Interpolation", p$runoption) )
    Si = p$runs[ iip, "locs" ]

    print( paste("index =", iip, ";  Si = ", Si ) )
    if ( Sflag[Si] == E[["complete"]] ) next()

    sloc = Sloc[Si,]
    nu    = S[Si, i_nu]
    phi   = S[Si, i_phi]
    localrange = S[Si, i_localrange]   # attached to p$stmv_autocorrelation_localrange
    varObs = S[Si, i_sdObs]^2,
    varSpatial = S[Si, i_sdSpatial]^2,

    # range checks
    if ( any( !is.finite( c(localrange, nu, phi) ) ) )  {
      Sflag[Si] %in% E[["variogram_failure"]]
      s1 = abs( sloc[1] - Sloc[,1] )
      s2 = abs( sloc[2] - Sloc[,2] )

      ii = NULL
      if ( !is.finite( localrange ) ) {
        # obtain estimate of localrange from adjoining areas
        localrange = max( distance_limits ) #initial guess
        ii = which( ( s1 <= localrange ) & ( s2 <= localrange ) )
        if (length( ii ) < 1)  next()
        localrange = median( S[ii, i_localrange ], na.rm=TRUE ) # initial local estimate
        ii = which( ( s1 <= localrange ) & ( s2 <= localrange ) )
        if (length( ii ) < 1)  next()
        localrange = median( S[ii, i_localrange ], na.rm=TRUE )  # refined local estimate
      }
      if ( !is.finite( localrange ) ) next()  # last check

      # the above forces localrange to always be available
      if ( !is.finite(nu) )  {
        if (is.null(ii))  ii = which( ( s1 <= localrange ) & ( s2 <= localrange ) ) # update
        nu =  median( S[ii, i_nu ], na.rm=TRUE )
      }
      if ( !is.finite(nu) ) next()  # last check

      # phi meaningful only with some given nu .. restimate from localrange
      if ( !is.finite(phi) ) phi = matern_distance2phi( dis=localrange, nu=nu, cor=p$stmv_autocorrelation_localrange )
      if ( !is.finite(phi) ) next() # last check

      s1 = NULL
      s2 = NULL
      ii = NULL
      gc()
    }

    if ( Sflag[Si] != E[["todo"]] ) {
      if (exists("stmv_rangecheck", p)) {
        if (p$stmv_rangecheck=="paranoid") {
          if ( Sflag[Si] %in% c( E[["variogram_range_limit"]], E[["variogram_failure"]]) ) {
            if (debugging) message("Error: stmv_rangecheck paranoid")
              next()
          }
        }
      }
    }

    if (!exists("stmv_interpolation_basis", p)) p$stmv_interpolation_basis = "correlation"

    if (p$stmv_interpolation_basis == "correlation")  {
      localrange_interpolation = matern_phi2distance( phi=phi, nu=nu, cor=p$stmv_interpolation_basis_correlation )
    } else if (p$stmv_interpolation_basis == "distance")  {
      localrange_interpolation = p$stmv_interpolation_basis_distance
    }

    data_subset = stmv_select_data( p=p, Si=Si, localrange=localrange_interpolation )
    if (is.null( data_subset )) {
      Sflag[Si] = E[["insufficient_data"]]
      next()
    }
    unique_spatial_locations = data_subset$unique_spatial_locations
    if (unique_spatial_locations < p$stmv_nmin) {
      data_subset = NULL;
      Sflag[Si] = E[["insufficient_data"]]
      next()
    }

    # ndata abovev is for unique locations .. now ndata is for dim of input data
    ndata = length(data_subset$data_index)
    S[Si, i_ndata] = ndata

    # if here then there is something to do
    # NOTE:: data_subset$data_index are the indices of locally useful data

    # prep dependent data
    # reconstruct data for modelling (dat)
    dat = matrix( 1, nrow=ndata, ncol=dat_nc )
    dat[,iY] = Y[data_subset$data_index] # these are residuals if there is a global model
    # add a small error term to prevent some errors when duplicate locations exist; localrange_interpolation offsets to positive values

    dat[,ilocs] = Yloc[data_subset$data_index,]
    dat[,ilocs] = dat[,ilocs] + localrange_interpolation * runif(2*ndata, -dist_error, dist_error)

    if (p$nloccov > 0) dat[,icov] = Ycov[data_subset$data_index, icov_local] # no need for other dim checks as this is user provided
    if (exists("TIME", p$stmv_variables)) {
      dat[, itime_cov] = as.matrix(stmv_timecovars( vars=ti_cov, ti=Ytime[data_subset$data_index,] ) )
      itt = which(dat_names==p$stmv_variables$TIME)
      dat[, itt ] = Ytime[data_subset$data_index,]
      # crude check of number of time slices
      n_time_slices = stmv_discretize_coordinates( coo=dat[, itt], ntarget=p$nt, minresolution=p$minresolution[3], method="thin"  )
      if ( length(n_time_slices) < p$stmv_tmin )  next()
    }
    dat = data.table(dat)
    names(dat) = dat_names
    dat_range = range( dat[,..iY], na.rm=TRUE )  # used later

    # remember that these are crude mean/discretized estimates
    if (debugging) {
      dev.new()
      # check data and statistical locations
      plot( Sloc[,], pch=20, cex=0.5, col="gray")
      points( Yloc[,], pch=20, cex=0.2, col="green")
      points( Yloc[data_subset$data_index,], pch=20, cex=1, col="yellow" )
      points( Sloc[Si,2] ~ Sloc[Si,1], pch=20, cex=5, col="blue" )
    }

    # construct prediction/output grid area ('pa')
    prediction_area = localrange_interpolation
    if (exists("stmv_distance_prediction_limits", p)) {
      prediction_area = min( max( localrange_interpolation, min(p$stmv_distance_prediction_limits) ), max(p$stmv_distance_prediction_limits), na.rm=TRUE )
    }

    windowsize.half =  aegis_floor( prediction_area / p$pres ) + 1L

    pa = try( stmv_predictionarea_lattice( p=p, sloc=sloc, windowsize.half=windowsize.half ) )

    if ( is.null(pa) ) {
      Sflag[Si] = E[["prediction_area"]]
      next()
    }

    if ( inherits(pa, "try-error") ) {
      pa = NULL
      Sflag[Si] = E[["prediction_area"]]
      if (debugging) message("Error: prediction grid ... try-error .. this should not happen.  check this")
      next()
    }

    if (debugging) {

      dev.new()
      # check that position indices are working properly
      Ploc = stmv_attach( p$storage_backend, p$ptr$Ploc )
      Sloc = stmv_attach( p$storage_backend, p$ptr$Sloc )
      Yloc = stmv_attach( p$storage_backend, p$ptr$Yloc )
      plot( Yloc[data_subset$data_index,2] ~ Yloc[data_subset$data_index,1], col="red", pch=".",
        ylim=range(c(Yloc[data_subset$data_index,2], Sloc[Si,2], Ploc[pa$i,2]) ),
        xlim=range(c(Yloc[data_subset$data_index,1], Sloc[Si,1], Ploc[pa$i,1]) ) ) # all data
      points( Yloc[data_subset$data_index,2] ~ Yloc[data_subset$data_index,1], col="green" )  # with covars and no other data issues
      points( Sloc[Si,2] ~ Sloc[Si,1], col="blue" ) # statistical locations
      # statistical output locations
      grids= spatial_grid(p, DS="planar.coords" )

      points( grids$plat[aegis_floor( (Sloc[Si,2]-p$origin[2])/p$pres) + 1]
            ~ grids$plon[aegis_floor( (Sloc[Si,1]-p$origin[1])/p$pres) + 1] , col="purple", pch=25, cex=5 )

      points( Ploc[pa$i,2] ~ Ploc[ pa$i, 1] , col="black", pch=6, cex=0.7 ) # check on pa$i indexing -- prediction locations
    }

    data_subset = NULL

    # model and prediction .. outputs are in scale of the link (and not response)
    # the following permits user-defined models (might want to use compiler::cmpfun )

    res =NULL
    res = try(
      local_fn ( p=p, dat=dat, pa=pa, nu=nu, phi=phi, localrange=localrange_interpolation,
        varObs=varObs, varSpatial=varSpatial, sloc=sloc
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

    if (length(which( is.finite(res$predictions$mean ))) < 1) {
      Sflag[Si] =  E[["prediction_error"]]   # modelling / prediction did not complete properly
      res = NULL
      if (debugging) {
        message("Error: prediction error")
        browser()
      }
      next()  # looks to be a faulty solution
    }

    if (grepl( Sflag[Si], paste0( c("local_model_error", "prediction_error", "prediction_area", "insufficient_data" ), collapse=" " )) ) next()



    tolimit = which( res$predictions$mean < dat_range[1])
    if (length(tolimit) > 0 ) {
      res$predictions$mean[tolimit] = dat_range[1]
      res$predictions$sd[tolimit] = NA
    }

    tolimit = which( res$predictions$mean > dat_range[2])
    if (length(tolimit) > 0 ) {
      res$predictions$mean[tolimit] = dat_range[2]
      res$predictions$sd[tolimit] = NA
    }


    # update in stats .. if there are updates
    if (!is.na( res$stmv_stats_by_time )) {
      # only produced when: if ( grepl("stmv_variogram_resolve_time", p$stmv_fft_filter))
      lss = colMeans( res$stmv_stats_by_time, na.rm=TRUE ) ## only with fft method ...  ..
      names(lss) = p$statsvars
      message( "flag: stmv_variogram_resolve_time used in p$stmv_fft_filter ... time-varying variogram parameters ... need to create as system to store/use this properly! .. incomplete")
    }


    for ( vv in 1:length(p$statsvars) ) {
      vn = p$statsvars[vv]
      vi = match(vn, p$statsvars)
      if ( is.finite(res$stmv_stats[[ vn ]] ) ) S[Si, vi] = res$stmv_stats[[ vn ]]
    }

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
      tst = na.omit(tst)
      mba.int <- mba.surf( tst, windowsize.half *2, windowsize.half *2, extend=TRUE)$xyz.est
      image(mba.int, xaxs="r", yaxs="r")

      # kernel-based
      tst = as.image( Z=dat$z, x=dat[, c("plon", "plat")], nx=300, ny=300, na.rm=TRUE)
      out = fields::image.smooth( tst, theta=phi/300, xwidth=p$pres, ywidth=p$pres )
      dev.new(); image(out)

      print( str(res) )

      lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$stmv_variables$TIME]==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
      lattice::levelplot( mean ~ plon + plat, data=res$predictions, col.regions=heat.colors(100), scale=list(draw=TRUE) , aspect="iso" )
      for( i in sort(unique(res$predictions[,p$stmv_variables$TIME])))  print(lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$stmv_variables$TIME]==i,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )

      dev.new()
      plot(  dat[,iY] ~ dat$yr, col="red"  )
      points( mean~tiyr, res$predictions, pch=20, col="gray", cex=0.5 )

    }


}
