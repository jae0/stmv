
stmv_parameters = function( p=list(), ... ) {

  # ---------------------
  # deal with additional passed parameters

  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast=TRUE))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable

  if (!exists("stmvSaveDir", p)) p$stmvSaveDir = file.path(p$data_root, "modelled", p$variables$Y, p$spatial.domain )
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

  if (!exists("stmv_clusters", p)) {
    if (exists( "clusters", p)) {
      p$stmv_clusters = p$clusters
    } else {
      warning( "p$stmv_clusters not specified, using 1 core only" )
      p$stmv_clusters = "localhost"
    }
  }


  if( !exists( "storage.backend", p))  p$storage.backend="bigmemory.ram"

  if( !exists( "stmv_variogram_method", p)) p$stmv_variogram_method="geoR"   # note GP methods are slow when there is too much data
  if( !exists( "stmv_range_correlation", p)) p$stmv_range_correlation = 0.1   # auto-correlation value at which to compute a "range" distance for estimation
  if( !exists( "stmv_range_correlation", p)) p$stmv_range_correlation_fft_taper = 0.5   # auto-correlation value at which to compute a tapered "range" distance for estimation

  if( !exists( "stmv_range_correlation_boostdata", p)) p$stmv_range_correlation_boostdata = 0.05   # auto-correlation value at which to compute a "range" distance for estimation

  if (!exists( "stmv_global_family", p)) p$stmv_global_family = gaussian(link = "identity")

  if (!exists( "boundary", p)) p$boundary = FALSE
  if (!exists( "depth.filter", p)) p$depth.filter = FALSE # if !FALSE .. depth is given as m so, choose andy stats locations with elevation > 1 m as being on land

  if (!exists( "stmv_nmin_downsize_factor", p)) p$stmv_nmin_downsize_factor = c(1.0, 0.9, 0.8, 0.7, 0.6, 0.5)

  if (!exists( "stmv_lowpass_phi", p)) p$stmv_lowpass_phi = p$pres*2 # FFT based method when operating gloablly
  if (!exists( "stmv_lowpass_nu", p)) p$stmv_lowpass_nu = 0.5 # this is exponential covar

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

  if (p$stmv_local_modelengine == "twostep") {
    if (exists("stmv_rsquared_threshold", p) ) {
      if (p$stmv_rsquared_threshold > 0) {
        # message( "Ignoring value of p$stmv_rsquared_threshold as it is meaningless with twostep")
        p$stmv_rsquared_threshold = 0  # override :: this is meaningless when broken apart in space and time ..
      }
    }
  }

  if ( p$stmv_local_modelengine %in% c("krige" )) {
     # nothing to add yet ..
  }


  if (!exists("variables", p)) p$variables=list()

  p = stmv_variablelist(p)

  # require knowledge of size of stats output which varies with a given type of analysis
  p$statsvars = c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "range", "phi", "nu" )
  if (exists("TIME", p$variables) )  p$statsvars = c( p$statsvars, "ar_timerange", "ar_1" )
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
  p$libs = unique( c( p$libs, "sp", "rgdal", "parallel" ) )
  if (!exists("storage.backend", p)) p$storage.backend = storage.backend
  if (any( grepl ("ff", p$storage.backend)))         p$libs = c( p$libs, "ff", "ffbase" )
  if (any( grepl ("bigmemory", p$storage.backend)))  p$libs = c( p$libs, "bigmemory" )
  if (p$storage.backend=="bigmemory.ram") {
    if ( length( unique(p$clusters)) > 1 ) {
      stop( "||| More than one unique cluster server was specified .. the bigmemory RAM-based method only works within one server." )
    }
  }

  # other libs
  if (exists("stmv_local_modelengine", p)) {
    if (p$stmv_local_modelengine=="bayesx")  p$libs = c( p$libs, "R2BayesX" )
    if (p$stmv_local_modelengine %in% c("gam", "mgcv") )  p$libs = c( p$libs, "mgcv" )
    if (p$stmv_local_modelengine %in% c("inla") )  p$libs = c( p$libs, "INLA" )
    if (p$stmv_local_modelengine %in% c("fft", "gaussianprocess2Dt") )  p$libs = c( p$libs, "fields" )
    if (p$stmv_local_modelengine %in% c("twostep") )  p$libs = c( p$libs, "mgcv", "fields" )
    if (p$stmv_local_modelengine %in% c("krige") ) p$libs = c( p$libs, "fields" )
    if (p$stmv_local_modelengine %in% c("gstat") ) p$libs = c( p$libs, "gstat" )
  }

  if (exists("stmv_global_modelengine", p)) {
    if (p$stmv_global_modelengine %in% c("gam", "mgcv") ) p$libs = c( p$libs, "mgcv" )
    if (p$stmv_global_modelengine %in% c("bigglm", "biglm") ) p$libs = c( p$libs, "biglm" )
  }

  p$libs = unique( p$libs )
  suppressMessages( RLibrary( p$libs ) )


  if ( !exists("minresolution", p)) {
    if ( exists("TIME", p$variables)) {
      p$minresolution = c(p$pres, p$pres, p$tres)
    } else {
      p$minresolution = c(p$pres, p$pres )
    }
  }

  p$nloccov = 0
  if (exists("local_cov", p$variables)) p$nloccov = length(p$variables$local_cov)

  if ( !exists("stmv_distance_prediction_fraction", p)) p$stmv_distance_prediction_fraction = 1  # fraction of statsgrid (below)

  if ( !exists("stmv_distance_prediction", p)) {
    # this is a half window km
    p$stmv_distance_prediction = p$stmv_distance_statsgrid * p$stmv_distance_prediction_fraction
  }

  # construct prediction/output grid area ('pa')
  if ( !exists("windowsize.half", p)) p$windowsize.half = floor(p$stmv_distance_prediction/p$pres) # convert distance to discretized increments of row/col indices; stmv_distance_prediction = 0.75* stmv_distance_statsgrid (unless overridden)

  return(p)
}
