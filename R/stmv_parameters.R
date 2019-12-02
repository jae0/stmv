
stmv_parameters = function( p=list(), ... ) {

  # ---------------------
  # deal with additional passed parameters

  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast=TRUE))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable

  if (!exists("stmvSaveDir", p)) p$stmvSaveDir = file.path(p$data_root, "modelled", p$stmv_variables$Y, p$spatial_domain )
  if ( !file.exists(p$stmvSaveDir)) dir.create( p$stmvSaveDir, recursive=TRUE, showWarnings=FALSE )

  if (!exists("stmv_current_status", p))  p$stmv_current_status = file.path( p$stmvSaveDir, "stmv_current_status" )

  if ( !exists("stmv_distance_scale", p)) {
    message( "||| stmv_distance_scale must be defined" )
    stop()
  }

  if ( !exists("stmv_distance_statsgrid", p)) {
    message( "||| stmv_distance_statsgrid must be defined" )
    stop()
  }


  if( !exists( "storage_backend", p))  p$storage_backend="bigmemory.ram"

  if( !exists( "stmv_variogram_method", p)) p$stmv_variogram_method = "fft"   # note GP methods are slow when there is too much data
  if( !exists( "stmv_variogram_nbreaks_totry", p)) p$stmv_variogram_nbreaks_totry = c( 16, 32, 64, 47, 39, 21, 13 )  # different numbers of nbreaks can influence variogram stabilty

  if( !exists( "stmv_autocorrelation_localrange", p)) p$stmv_autocorrelation_localrange = 0.1   # auto-correlation at which to compute range distance
  if( !exists( "stmv_autocorrelation_fft_taper", p))  p$stmv_autocorrelation_fft_taper = 0.75   # scale at which to mark tapering

  if (!exists( "stmv_global_family", p)) p$stmv_global_family = gaussian(link = "identity")

  if (!exists( "boundary", p)) p$boundary = FALSE
  if (!exists( "stmv_filter_depth_m", p)) p$stmv_filter_depth_m = FALSE # if !FALSE .. depth is given as m so, choose andy stats locations with elevation > 1 m as being on land

  if (!exists( "stmv_nmin_downsize_factor", p)) p$stmv_nmin_downsize_factor = c(1.0, 0.9, 0.8, 0.7, 0.6)

  if (!exists( "stmv_lowpass_nu", p)) p$stmv_lowpass_nu = 0.5 # this is exponential covar
  if (!exists( "stmv_lowpass_phi", p)) p$stmv_lowpass_phi = stmv::matern_distance2phi( distance=ifelse(exists("pres", p), p$pres, stop("'p$pres' needs to be defined")), nu=p$stmv_lowpass_nu, cor=p$stmv_autocorrelation_localrange ) # FFT based method when operating gloablly

  # used by "fields" GRMF functions
  if ( p$stmv_local_modelengine %in% c("gaussianprocess2Dt")) {
    if (!exists("phi.grid", p) ) p$phi.grid = 10^seq( -6, 6, by=0.5)  # maxdist is aprox magnitude of the phi parameter
    if (!exists("lambda.grid", p) ) p$lambda.grid = 10^seq( -9, 3, by=0.5) # ratio of tau sq to sigma sq
  }

  if (!exists("stmv_gam_optimizer", p)) p$stmv_gam_optimizer=c("outer","bfgs")

  if ( p$stmv_local_modelengine %in% c("gam" )) {
    # p$stmv_gam_optimizer=c("outer","optim")
    # if (!exists("stmv_gam_optimizer", p)) p$stmv_gam_optimizer="perf"
  }

  if ( p$stmv_local_modelengine %in% c("krige" )) {
     # nothing to add yet ..
  }



  if (!exists("stmv_variables", p)) p$stmv_variables=list()

  p = stmv_variablelist(p)

  # require knowledge of size of stats output which varies with a given type of analysis
  p$statsvars = c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "phi", "nu", "localrange" )
  if (exists("TIME", p$stmv_variables) )  p$statsvars = c( p$statsvars, "ar_timerange", "ar_1" )
  if (p$stmv_local_modelengine == "userdefined" ) {
    if (exists("stmv_local_modelengine", p) ) {
      if (exists("stmv_local_modelengine_userdefined", p) ) {
        if (class(p$stmv_local_modelengine_userdefined) == "function" ) {
          oo = NULL
          oo = try( p$stmv_local_modelengine_userdefined(variablelist=TRUE), silent=TRUE )
          if ( ! inherits(oo, "try-error") ) p$statsvars = unique( c( p$statsvars, oo ) )
        }
      }
    }
  }

  # determine storage format
  p$libs = unique( c( p$libs, "sp", "rgdal", "interp", "parallel" ) )
  if (!exists("storage_backend", p)) p$storage_backend = storage_backend
  if (any( grepl ("ff", p$storage_backend)))         p$libs = c( p$libs, "ff", "ffbase" )
  if (any( grepl ("bigmemory", p$storage_backend)))  p$libs = c( p$libs, "bigmemory" )
  if (p$storage_backend=="bigmemory.ram") {
    if ( length( unique(p$clusters)) > 1 ) {
      stop( "||| More than one unique cluster server was specified .. the bigmemory RAM-based method only works within one server." )
    }
  }

  # other libs
  if (!exists("stmv_local_modelengine", p)) p$stmv_local_modelengine = "none"
  if (p$stmv_local_modelengine=="bayesx")  p$libs = c( p$libs, "R2BayesX" )
  if (p$stmv_local_modelengine %in% c("gam", "mgcv") )  p$libs = c( p$libs, "mgcv" )
  if (p$stmv_local_modelengine %in% c("inla") )  p$libs = c( p$libs, "INLA" )
  if (p$stmv_local_modelengine %in% c("fft", "gaussianprocess2Dt") )  p$libs = c( p$libs, "fields", "fftwtools", "fftw" )
  if (p$stmv_local_modelengine %in% c("twostep") )  p$libs = c( p$libs, "mgcv", "fields", "fftwtools", "fftw"  )
  if (p$stmv_local_modelengine %in% c("krige") ) p$libs = c( p$libs, "fields" )
  if (p$stmv_local_modelengine %in% c("gstat") ) p$libs = c( p$libs, "gstat" )

  if (!exists("stmv_global_modelengine", p)) p$stmv_global_modelengine = "none"   # required variable .. default to none
  if (p$stmv_global_modelengine %in% c("gam", "mgcv") ) p$libs = c( p$libs, "mgcv" )
  if (p$stmv_global_modelengine %in% c("bigglm", "biglm") ) p$libs = c( p$libs, "biglm" )

  p$libs = unique( p$libs )
  suppressMessages( RLibrary( p$libs ) )


  if ( !exists("minresolution", p)) {
    if ( exists("TIME", p$stmv_variables)) {
      p$minresolution = c(p$pres, p$pres, p$tres)
    } else {
      p$minresolution = c(p$pres, p$pres )
    }
  }

  p$nloccov = 0
  if (exists("local_cov", p$stmv_variables)) p$nloccov = length(p$stmv_variables$local_cov)

  if ( !exists("stmv_force_complete_method", p)) p$stmv_force_complete_method = "linear"  # moving average

  if ( !exists("stmv_rsquared_threshold", p) ) p$stmv_rsquared_threshold = 0  # essentially ignore ..

  return(p)
}
