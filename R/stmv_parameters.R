
stmv_parameters = function( p=list(), ... ) {

  p = parameters_add(p, list(...)) # add passed args to parameter list, priority to args

  if (! exists("data_root", p) ) stop( "||| data_root is required" )
  if (! exists("stmv_variables", p) ) stop( "||| stmv_variables is required" )
  if (! exists("Y", p$stmv_variables) ) stop( "|||  Y is required in p$stmv_variables" )
  if (! exists("spatial_domain", p) ) stop( "||| spatial_domain is required" )
  if (! exists("stmv_distance_statsgrid", p)) stop( "||| stmv_distance_statsgrid must be defined" )

  if (! exists("stmv_model_label", p) ) stmv_model_label = "default"

  p = parameters_add_without_overwriting( p,
    project_class = "stmv",
    stmv_model_label = "default",
    stmv_global_modelengine = "none",
    stmv_local_modelengine = "none"
  )

  p = parameters_add_without_overwriting( p,
    stmvSaveDir = file.path( p$modeldir, p$stmv_model_label, p$project_class, paste(  p$stmv_global_modelengine, p$stmv_local_modelengine, sep="_"), p$stmv_variables$Y, p$spatial_domain)
  )

  if ( !file.exists(p$stmvSaveDir)) dir.create( p$stmvSaveDir, recursive=TRUE, showWarnings=FALSE )

  p = parameters_add_without_overwriting( p,
    stmv_current_status = file.path( p$stmvSaveDir, "stmv_current_status" ),
    stmv_control_file = file.path( p$stmvSaveDir, "stmv_control" )
  )

  p = parameters_add_without_overwriting( p,
    storage_backend="bigmemory.ram",
    stmv_variogram_method = "fft",   # note GP methods are slow when there is too much data
    stmv_variogram_nbreaks_totry = c( 16, 32, 64, 47, 39, 21, 13 ),  # different numbers of nbreaks can influence variogram stabilty
    stmv_autocorrelation_localrange = 0.1,   # auto-correlation at which to compute range distance
    stmv_autocorrelation_fft_taper = 0.75,   # scale at which to mark tapering
    stmv_global_family = gaussian(link = "identity"),
    boundary = FALSE,
    stmv_filter_depth_m = FALSE, # if !FALSE .. depth is given as m so, choose andy stats locations with elevation > 1 m as being on land
    stmv_nmin_downsize_factor = c(1.0, 0.8, 0.6, 0.4, 0.2, 0.1),
    stmv_lowpass_nu = 0.5, # this is exponential covar
    stmv_lowpass_phi = stmv::matern_distance2phi( distance=ifelse(exists("pres", p), p$pres, stop("'p$pres' needs to be defined")), nu=p$stmv_lowpass_nu, cor=p$stmv_autocorrelation_localrange ) # FFT based method when operating gloablly
  )
 

  if ( p$stmv_local_modelengine %in% c("gaussianprocess2Dt")) {
    p$libs = c( p$libs, "fields", "fftwtools", "fftw" )
    p = parameters_add_without_overwriting( p,
      phi.grid = 10^seq( -6, 6, by=0.5),  # maxdist is aprox magnitude of the phi parameter
      lambda.grid = 10^seq( -9, 3, by=0.5) # ratio of tau sq to sigma sq
    )
  }
  if (p$stmv_local_modelengine %in% c("inla") )  p$libs = c( p$libs, "INLA" )
  if (p$stmv_local_modelengine %in% c("twostep") )  p$libs = c( p$libs, "mgcv", "fields", "fftwtools", "fftw"  )
  if (p$stmv_local_modelengine %in% c("krige") ) p$libs = c( p$libs, "fields" )
  if (p$stmv_local_modelengine %in% c("gstat") ) p$libs = c( p$libs, "gstat" )
  if (p$stmv_global_modelengine %in% c("bigglm", "biglm") ) p$libs = c( p$libs, "biglm" )
  if (p$stmv_local_modelengine %in% c("bayesx") )  p$libs = c( p$libs, "R2BayesX" )
  if (p$stmv_local_modelengine %in% c("gam", "mgcv") )  {
    p$libs = c( p$libs, "mgcv" )
    p = parameters_add_without_overwriting( p, stmv_gam_optimizer=c("outer","bfgs") ) # "perf"
  }
  if (p$stmv_local_modelengine %in% c("carstm") ) {
    p$libs = c( p$libs, "INLA", "sf", "sp", "spdep" )
    p = parameters_add_without_overwriting( p,
      stmv_au_distance_reference = "none", # additional filters upon polygons relative to windowsize: "centroid", "inside_or_touches_boundary", completely_inside_boundary"
      stmv_au_buffer_links = 0, # number of additional neighbours to extend beyond initial solution
      pres = 1  # this governs resolution of lattice predictions
    )
    p = parameters_add_without_overwriting( p,
      control.inla.variations = list(
        list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive", h=0.001, cmin=0),
        list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive", h=0.01, cmin=0),
        list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive", h=0.1,  cmin=0),
        list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive") ,
        list( optimise.strategy="smart", stupid.search=TRUE, h=0.0001 ),
        list( optimise.strategy="smart", stupid.search=TRUE, h=0.001 ),
        list( optimise.strategy="smart", stupid.search=TRUE, h=0.01 ),
        list( optimise.strategy="smart", stupid.search=TRUE, h=0.1 )
      )
    )

  }

  # determine storage format
  p$libs = unique( c( p$libs, "sp", "spdep", "sf", "rgdal", "interp", "parallel", "raster", "INLA" ) )

  if (any( grepl ("ff", p$storage_backend)))         p$libs = c( p$libs, "ff", "ffbase" )
  if (any( grepl ("bigmemory", p$storage_backend)))  p$libs = c( p$libs, "bigmemory" )
  if (p$storage_backend=="bigmemory.ram") {
    if ( length( unique(p$clusters)) > 1 ) {
      stop( "||| More than one unique cluster server was specified .. the bigmemory RAM-based method only works within one server." )
    }
  }

  # tidy libs
  p$libs = unique( p$libs )
  suppressMessages( RLibrary( p$libs ) )


  if (!exists("stmv_variables", p)) p$stmv_variables=list()
  p = stmv_variablelist(p)

  minresolution = c(p$pres, p$pres )
  if ( exists("TIME", p$stmv_variables) ) minresolution = c(p$pres, p$pres, p$tres)

  nloccov = 0
  if ( exists("local_cov", p$stmv_variables) ) nloccov = length(p$stmv_variables$local_cov)

  p = parameters_add_without_overwriting( p,
    minresolution = minresolution,
    nloccov = nloccov,
    stmv_force_complete_method = "linear",  # moving average
    stmv_rsquared_threshold = 0  # essentially ignore ..
  )


  return(p)
}
