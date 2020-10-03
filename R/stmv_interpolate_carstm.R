
stmv_interpolate_carstm = function( ip=NULL, p, debugging=FALSE, eps = 1e-6, global_sppoly=NULL, ... ) {
  #\\ core function to interpolate (model variogram) in parallel
  # carstm has additional polygon prep and classification steps and does not rely upon an autocorrelation range/scale

  if (0) {
    # for debugging  runs ..
    currentstatus = stmv_statistics_status( p=p )
    p = parallel_run( p=p, runindex=list( locs=sample( currentstatus$todo )) )
    # parallel_run( stmv_carstm, p=p, runindex=list( locs=sample( currentstatus$todo )) )
    p$runmode = "carstm"
    ip = 1:p$nruns
    debugging=TRUE
    eps = 1e-6
    stmv_interpolation_basis_distance = 2 * p$stmv_distance_statsgrid  # fixed distance

    require(INLA)
    p$carstm_modelcall = paste(
          'inla(
            formula = ', p$stmv_variables$Y, ' ~ 1
              + f(aui, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2),
            family = "normal",
            data= dat,
            control.compute=list(dic=TRUE, waic=TRUE, cpo=FALSE, config=FALSE),  # config=TRUE if doing posterior simulations
            control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
            control.predictor=list(compute=FALSE, link=1 ),
            control.fixed=H$fixed,  # priors for fixed effects, generic is ok
            control.inla = list(h=1e-4, tolerance=1e-9, cmin=0), # restart=3), # restart a few times in case posteriors are poorly defined
            verbose=TRUE
          ) '
    )

  }




  p = parameters_control(p, list(...), control="add") # add passed args to parameter list, priority to args

  if (exists( "libs", p)) suppressMessages( RLibrary( p$libs ) )

  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns


  #---------------------
  # data for modelling
  S = stmv_attach( p$storage_backend, p$ptr$S )
  Sloc = stmv_attach( p$storage_backend, p$ptr$Sloc )

  Sflag = stmv_attach( p$storage_backend, p$ptr$Sflag )
  E = stmv_error_codes()

  Yloc = stmv_attach( p$storage_backend, p$ptr$Yloc )
  Y = stmv_attach( p$storage_backend, p$ptr$Y )

  if (p$nloccov > 0) Ycov = stmv_attach( p$storage_backend, p$ptr$Ycov )

  # misc intermediate calcs to be done outside of parallel loops

  # pre-calculate indices and dim for data to use inside the loop
  dat_names = unique( c(  p$stmv_variables$Y, p$stmv_variables$LOCS, p$stmv_variables$local_all ) )  # excludes p$stmv_variables$TIME
  # unless it is an explicit covariate and not a seasonal component there is no need for it
  # .. prediction grids create these from a time grid on the fly
  dat_nc = length( dat_names )

  iY = which(dat_names== p$stmv_variables$Y)
  ilocs = which( dat_names %in% p$stmv_variables$LOCS )
  i_ndata = match( "ndata", p$statsvars )
  # i_sdSpatial = match("sdSpatial",   p$statsvars)
  # i_sdObs = match("sdObs",   p$statsvars)
  # i_rsquared = match("rsquared", p$statsvars )
  # i_localrange = match("localrange", p$statsvars )
  # i_nu = match("nu",   p$statsvars)
  # i_phi = match("phi",   p$statsvars)

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


    # obtain indices of data locations within a given spatial range, optimally determined via variogram
    # find data nearest Sloc[Si,] and with sufficient data
    ndata = 0

    localrange_interpolation = ifelse( !exists("stmv_interpolation_basis_distance", p), p$stmv_distance_statsgrid *1.5, p$stmv_interpolation_basis_distance )

    data_subset = NULL
    data_subset = stmv_select_data( p=p, Si=Si, localrange=localrange_interpolation )
    if (is.null( data_subset )) {
      Sflag[Si] = E[["insufficient_data"]]
      next()
    }

    unique_spatial_locations = data_subset$unique_spatial_locations
    ndata = length(data_subset$data_index)
    if (unique_spatial_locations < p$stmv_nmin) {
      data_subset = NULL;
      Sflag[Si] = E[["insufficient_data"]]
      next()
    }
    # ndata abovev is for unique locations .. now ndata is for dim of input data
    S[Si, i_ndata] = ndata

    # prep data for modelling
    # NOTE:: data_subset$data_index are the indices of locally useful data

    # prep dependent data
    # reconstruct data for modelling (dat)
    dat = matrix( 1, nrow=ndata, ncol=dat_nc )
    dat[,iY] = Y[data_subset$data_index] # these are residuals if there is a global model
    dat[,ilocs] = Yloc[data_subset$data_index,]

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
    # convert distance to discretized increments of row/col indices;
    windowsize.half =  floor( localrange_interpolation / p$pres ) + 1L
    # construct data (including covariates) for prediction locations (pa)
    pa = try( stmv_predictionarea( p=p, sloc=Sloc[Si,], windowsize.half=windowsize.half ) )
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

    # determine spatial polygons for prediction .. pa is for space only at this point .. note pa only has static vars and covars
    sppoly = stmv_predictionarea_polygon( pa=pa, dx=p$pres, dy=p$pres, pa_coord_names=p$stmv_variables$LOCS[1:2], pa_proj4string_planar_km=p$aegis_proj4string_planar_km, global_sppoly=global_sppoly )

    if ( exists("TIME", p$stmv_variables) )  pa = try( stmv_predictiontime( p=p, pa=pa ) ) # add time to pa and time varying covars

    if ( is.null(pa) ) {
      Sflag[Si] = E[["prediction_time"]]
      next()
    }

    if ( inherits(pa, "try-error") ) {
      pa = NULL
      Sflag[Si] = E[["prediction_time"]]
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

      points( grids$plat[floor( (Sloc[Si,2]-p$origin[2])/p$pres) + 1]
            ~ grids$plon[floor( (Sloc[Si,1]-p$origin[1])/p$pres) + 1] , col="purple", pch=25, cex=5 )

      points( Ploc[pa$i,2] ~ Ploc[ pa$i, 1] , col="black", pch=6, cex=0.7 ) # check on pa$i indexing -- prediction locations

      print( paste("index =", iip, ";  Sflag = ", names(E)[match(Sflag[Si], E)]  ) )
    }

    # if here then there is something to do


    # remember that these are crude mean/discretized estimates
    if (debugging) {
        dev.new()
        # check data and statistical locations
        plot( Sloc[,], pch=20, cex=0.5, col="gray")
        points( Yloc[,], pch=20, cex=0.2, col="green")
        points( Yloc[data_subset$data_index,], pch=20, cex=1, col="yellow" )
        points( Sloc[Si,2] ~ Sloc[Si,1], pch=20, cex=5, col="blue" )
    }


    res = stmv__carstm( p=NULL, dat=NULL, pa=NULL, sppoly=NULL  )

    dat =  NULL
    pa = NULL
    sppoly = NULL

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

    if ( exists("stmv_localstats", res) ) {
      lss = res$stmv_localstats


      ------

      if (is.finite(lss[["nu"]])) S[Si, i_nu] = lss[["nu"]]
      if (is.finite(lss[["phi"]])) S[Si, i_phi] = lss[["phi"]]
      if (is.finite(lss[["localrange"]])) S[Si, i_localrange] = lss[["localrange"]]
      if (is.finite(lss[["rsquared"]])) S[Si, i_rsquared] = lss[["rsquared"]]
      if (is.finite(lss[["sdSpatial"]])) S[Si, i_sdSpatial] = lss[["sdSpatial"]]
      if (is.finite(lss[["sdObs"]])) S[Si, i_sdObs] = lss[["sdObs"]]
      if (is.finite(lss[["ndata"]])) S[Si, i_ndata] = lss[["ndata"]]
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

    om  = o[[p$stmv_variogram_method]] # save stats

    statvars_scale = c(
      sdTotal =sqrt( o$varZ),
      sdSpatial = sqrt(om$varSpatial) ,
      sdObs = sqrt(om$varObs),
      phi = om$phi,
      nu = om$nu,
      localrange = matern_phi2distance( phi=om$phi, nu=om$nu, cor=p$stmv_autocorrelation_localrange ),
      ndata=ndata
    )

    S[Si,match( names(statvars_scale), p$statsvars )] = statvars_scale

    S[Si, match( names(statsvars_time), p$statsvars )] = statsvars_time



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




    # ----------------------
    # do last. it is an indicator of completion of all tasks
    # restarts would be broken otherwise
    Sflag[Si] = E[["complete"]]  # mark as complete without issues

  }  # end for loop


  return(NULL)

      if (0) {
        # stats
        # p$statsvars = c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "phi", "nu", "localrange" )
        S = stmv_attach( p$storage_backend, p$ptr$S )
        sbox = list(
          plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$stmv_distance_statsgrid ),
          plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$stmv_distance_statsgrid ) )
        # statistics coordinates
        locations = as.matrix( expand_grid_fast( sbox$plons, sbox$plats ))
        levelplot(S[,match("localrange", p$statsvars )]~ locations[,1]+locations[,2])
        levelplot(S[,match("nu", p$statsvars )]~ locations[,1]+locations[,2])
        levelplot(S[,match("sdTotal", p$statsvars )]~ locations[,1]+locations[,2])
        levelplot(S[,match("sdSpatial", p$statsvars )]~ locations[,1]+locations[,2])

      }



}
