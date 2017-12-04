
stm_parameters = function( p=NULL, ... ) {

  # ---------------------
  # deal with additional passed parameters
  if ( is.null(p) ) p=list()
  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast=TRUE))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable

  if (!exists("stm_current_status", p))  p$stm_current_status = file.path( p$stmSaveDir, "stm_current_status" )

  if (!exists("clusters", p)) p$clusters = rep("localhost", detectCores() )  # default if not given
  if( !exists( "storage.backend", p))  p$storage.backend="bigmemory.ram"
  if( !exists( "stm_variogram_method", p)) p$stm_variogram_method="fast"   # note GP methods are slow when there is too much data
  if (!exists( "stm_global_family", p)) p$stm_global_family = gaussian()
  if (!exists( "stm_eps", p)) p$stm_eps = 0.001  # distance units for eps noise to permit mesh gen for boundaries
  if (!exists( "stm_quantile_bounds", p)) p$stm_quantile_bounds = c(0.01, 0.99) # remove these extremes in interpolations
  if (!exists( "eps", p)) p$eps = 1e-6 # floating point precision
  if (!exists( "boundary", p)) p$boundary = FALSE
  if (!exists( "depth.filter", p)) p$depth.filter = FALSE # depth is given as log(depth) so, choose andy stats locations with elevation > 1 m as being on land
  if (!exists( "stm_kernelmethods_use_all_data", p)) p$stm_kernelmethods_use_all_data =TRUE ## speed and RAM usage improvement is minimal (if any) when off, leave on or remove option and fix as on
  if (!exists( "stm_multiplier_stage2", p) ) p$stm_multiplier_stage2 = c( 1.1, 1.25 ) # distance multiplier for stage 2 interpolations

  if ( !exists("sampling", p))  {
    # fractions of distance scale  to try in local block search
    p$sampling = c( 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.5, 1.75, 2 )
  }

  # used by "fields" GRMF functions
  if ( p$stm_local_modelengine %in% c("gaussianprocess2Dt")) {
    if (!exists("phi.grid", p) ) p$phi.grid = 10^seq( -6, 6, by=0.5) * p$stm_distance_scale # maxdist is aprox magnitude of the phi parameter
    if (!exists("lambda.grid", p) ) p$lambda.grid = 10^seq( -9, 3, by=0.5) # ratio of tau sq to sigma sq
  }

  if ( p$stm_local_modelengine %in% c("gam" )) {
    # p$stm_gam_optimizer=c("outer","optim")
    if (!exists("stm_gam_optimizer", p)) p$stm_gam_optimizer=c("outer","bfgs")
    # if (!exists("stm_gam_optimizer", p)) p$stm_gam_optimizer="perf"
  }

  if ( p$stm_local_modelengine %in% c("krige" )) {
     # nothing to add yet ..
  }



  # require knowledge of size of stats output which varies with a given type of analysis
  p$statsvars = c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "range", "phi", "nu" ) 
  if (exists("TIME", p$variables) )  p$statsvars = c( p$statsvars, "ar_timerange", "ar_1" )
  if (p$stm_local_modelengine == "userdefined" ) {
    if (exists("stm_local_modelengine", p) ) {
      if (exists("stm_local_modelengine_userdefined_function", p) ) {
        if (class(p$stm_local_modelengine_userdefined_function) == "function" ) {
          oo = NULL
          oo = try( p$stm_local_modelengine_userdefined_function(variablelist=TRUE), silent=TRUE )
          if ( ! inherits(oo, "try-error") ) p$statsvars = unique( c( p$statsvars, oo ) )
        }
      }
    }
  }

  # determine storage format
  p$libs = unique( c( p$libs, "sp", "rgdal", "parallel", "RandomFields", "geoR" ) )
  if (!exists("storage.backend", p)) p$storage.backend = storage.backend
  if (any( grepl ("ff", p$storage.backend)))         p$libs = c( p$libs, "ff", "ffbase" )
  if (any( grepl ("bigmemory", p$storage.backend)))  p$libs = c( p$libs, "bigmemory" )
  if (p$storage.backend=="bigmemory.ram") {
    if ( length( unique(p$clusters)) > 1 ) {
      stop( "||| More than one unique cluster server was specified .. the bigmemory RAM-based method only works within one server." )
    }
  }

  # other libs
  if (exists("stm_local_modelengine", p)) {
    if (p$stm_local_modelengine=="bayesx")  p$libs = c( p$libs, "R2BayesX" )
    if (p$stm_local_modelengine %in% c("gam", "mgcv") )  p$libs = c( p$libs, "mgcv" )
    if (p$stm_local_modelengine %in% c("inla") )  p$libs = c( p$libs, "INLA" )
    if (p$stm_local_modelengine %in% c("fft", "gaussianprocess2Dt") )  p$libs = c( p$libs, "fields" )
    if (p$stm_local_modelengine %in% c("twostep") )  p$libs = c( p$libs, "mgcv", "fields" )
    if (p$stm_local_modelengine %in% c("krige") ) p$libs = c( p$libs, "fields" )
    if (p$stm_local_modelengine %in% c("gstat") ) p$libs = c( p$libs, "gstat" )
  }

  if (exists("stm_global_modelengine", p)) {
    if (p$stm_global_modelengine %in% c("gam", "mgcv") ) p$libs = c( p$libs, "mgcv" )
    if (p$stm_global_modelengine %in% c("bigglm", "biglm") ) p$libs = c( p$libs, "biglm" )
  }


  p$libs = unique( p$libs )
  suppressMessages( RLibrary( p$libs ) )

  p$variables$all = NULL
  if (exists("stm_local_modelformula", p))  {
    p$variables$local_all = all.vars( p$stm_local_modelformula )
    p$variables$local_cov = intersect( p$variables$local_all, p$variables$COV )
    p$variables$all = unique( c( p$variables$all, p$variables$local_all ) )
  }
  if (exists("stm_global_modelformula", p)) {
    p$variables$global_all = all.vars( p$stm_global_modelformula )
    p$variables$global_cov = intersect( p$variables$global_all, p$variables$COV )
    p$variables$all = unique( c( p$variables$all, p$variables$global_all ) )
  }
  p$stm_current_status = file.path( p$stmSaveDir, "stm_current_status" )





  return(p)
}
