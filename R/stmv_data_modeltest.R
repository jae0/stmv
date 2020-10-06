stmv_data_modeltest = function(p, runmode="default", global_sppoly=NULL ) {

  # purpose: to obtain param names from a model run
  # param names vary depending upon the model form
  # only possible to run a model to see waht gets exported ..
  # this is essentially a minimal form of stmv_interpolate  that stops once a solution with params are found

  Sloc = stmv_attach( p$storage_backend, p$ptr$Sloc )

  Yloc = stmv_attach( p$storage_backend, p$ptr$Yloc )
  Y = stmv_attach( p$storage_backend, p$ptr$Y )

  if (p$nloccov > 0) Ycov = stmv_attach( p$storage_backend, p$ptr$Ycov )
  if ( exists("TIME", p$stmv_variables) ) Ytime = stmv_attach( p$storage_backend, p$ptr$Ytime )


  # pre-calculate indices and dim for data to use inside the loop
  dat_names = unique( c(  p$stmv_variables$Y, p$stmv_variables$LOCS, p$stmv_variables$local_all, "weights" ) )  # excludes p$stmv_variables$TIME
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

  if (runmode=="default") {
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

  message("testing a run of the model to check for output")

  p = parallel_run( p=p, runindex=list( locs=sample( stmv_statistics_status( p=p )$todo )) )
  ip = 1:100

  for ( iip in ip ) {

    Si = p$runs[ iip, "locs" ]

    localrange_interpolation = ifelse( !exists("stmv_interpolation_basis_distance", p), p$stmv_distance_statsgrid *1.5, p$stmv_interpolation_basis_distance )

    data_subset = stmv_select_data( p=p, Si=Si, localrange=localrange_interpolation )
    if (is.null( data_subset )) next()

    unique_spatial_locations = data_subset$unique_spatial_locations
    ndata = length(data_subset$data_index)
    if (unique_spatial_locations < p$stmv_nmin) next()

    dat = matrix( 1, nrow=ndata, ncol=dat_nc )
    dat[,iY] = Y[data_subset$data_index] # these are residuals if there is a global model
       dat[,ilocs] = Yloc[data_subset$data_index,]
    if (runmode=="default") {
      dat[,ilocs] = dat[,ilocs] + localrange_interpolation * runif(2*ndata, -dist_error, dist_error)
    }

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

    prediction_area = localrange_interpolation
    if (exists("stmv_distance_prediction_limits", p)) {
      prediction_area = min( max( localrange_interpolation, min(p$stmv_distance_prediction_limits) ), max(p$stmv_distance_prediction_limits), na.rm=TRUE )
    }

    windowsize.half =  aegis_floor( prediction_area / p$pres ) + 1L
    # construct data (including covariates) for prediction locations (pa)
    pa = try( stmv_predictionarea( p=p, sloc=Sloc[Si,], windowsize.half=windowsize.half ) )
    if ( is.null(pa) ) next()
    if ( inherits(pa, "try-error") ) next()

    if ( "carstm" %in% runmode ) {
      # determine spatial polygons for prediction .. pa is for space only at this point .. note pa only has static vars and covars
      sppoly = stmv_predictionarea_polygon( pa=pa, dx=p$pres, dy=p$pres, pa_coord_names=p$stmv_variables$LOCS[1:2], pa_proj4string_planar_km=p$aegis_proj4string_planar_km, global_sppoly=global_sppoly )
     }

    if ( exists("TIME", p$stmv_variables) )  pa = try( stmv_predictiontime( p=p, pa=pa ) ) # add time to pa and time varying covars

    if ( is.null(pa) )  next()
    if ( inherits(pa, "try-error") )  next()

    res =NULL

browser()

    if ( runmode=="default" ) {
      vartot = var(dat[,iY], na.rm=TRUE)
      res = try(
        local_fn (
          p=p,
          dat=dat,
          pa=pa,
          nu=0.5,
          phi=localrange_interpolation/sqrt(8),
          localrange=localrange_interpolation,
          varObs = vartot/2,
          varSpatial = vartot/2,
          sloc = Sloc[Si,]
        )
      )
    } else if ( runmode=="carstm" ) {
      res = local_fn( p=p, dat=dat, pa=pa, sppoly=sppoly  )
    }

    if ( is.null(res))  next()
    if ( inherits(res, "try-error") ) next()
    if (!exists("predictions", res))  next()
    if (!exists("mean", res$predictions)) next()
    if (length(which( is.finite(res$predictions$mean ))) < 1) next()

    return( res )

  }
}