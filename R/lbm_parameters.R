

ecmei_parameters = function( p=NULL  ) {
  # some generic defaults
  if (is.null(p)) stop( "Parameter list is not structured properly" )

  if (!exists("ecmei_current_status", p))  p$ecmei_current_status = file.path( p$savedir, "ecmei_current_status" )

  if (!exists("clusters", p)) p$clusters = rep("localhost", detectCores() )  # default if not given
  if( !exists( "storage.backend", p))  p$storage.backend="bigmemory.ram"
  if( !exists( "ecmei_variogram_method", p)) p$ecmei_variogram_method="fast"   # note GP methods are slow when there is too much data
  if (!exists( "ecmei_global_family", p)) p$ecmei_global_family = gaussian()
  if (!exists( "ecmei_eps", p)) p$ecmei_eps = 0.001  # distance units for eps noise to permit mesh gen for boundaries
  if (!exists( "ecmei_quantile_bounds", p)) p$ecmei_quantile_bounds = c(0.01, 0.99) # remove these extremes in interpolations
  if (!exists( "eps", p)) p$eps = 1e-6 # floating point precision
  if (!exists( "boundary", p)) p$boundary = FALSE
  if (!exists( "depth.filter", p)) p$depth.filter = FALSE # depth is given as log(depth) so, choose andy stats locations with elevation > 1 m as being on land
  if (!exists( "ecmei_kernelmethods_use_all_data", p)) p$ecmei_kernelmethods_use_all_data =TRUE ## speed and RAM usage improvement is minimal (if any) when off, leave on or remove option and fix as on
  if (!exists( "ecmei_multiplier_stage2", p) ) p$ecmei_multiplier_stage2 = c( 1.1, 1.25 ) # distance multiplier for stage 2 interpolations

  # used by "fields" GRMF functions
  if ( p$ecmei_local_modelengine %in% c("gaussianprocess2Dt", "gaussianprocess" )) {
    if (!exists("phi.grid", p) ) p$phi.grid = 10^seq( -6, 6, by=0.5) * p$ecmei_distance_scale # maxdist is aprox magnitude of the phi parameter
    if (!exists("lambda.grid", p) ) p$lambda.grid = 10^seq( -9, 3, by=0.5) # ratio of tau sq to sigma sq
  }

  if ( p$ecmei_local_modelengine %in% c("gam" )) {
    # p$ecmei_gam_optimizer=c("outer","optim")
    if (!exists("ecmei_gam_optimizer", p)) p$ecmei_gam_optimizer=c("outer","bfgs")
    # if (!exists("ecmei_gam_optimizer", p)) p$ecmei_gam_optimizer="perf"
  }

  if ( p$ecmei_local_modelengine %in% c("krige" )) {
     # nothing to add yet ..
  }

  if ( p$ecmei_local_modelengine %in% c("spate")) {
    if (!exists( "ecmei_spate_method", p)) p$ecmei_spate_method = "mcmc_fast"  # "mcmc" is too slow for temperature
    if (!exists( "ecmei_spate_boost_timeseries", p)) p$ecmei_spate_boost_timeseries = TRUE 
    if (!exists( "ecmei_spate_nburnin", p)) p$ecmei_spate_nburnin =1000
    if (!exists( "ecmei_spate_nposteriors", p)) p$ecmei_spate_nposteriors = 4000
    if (!exists( "ecmei_spate_nCovUpdates", p)) p$ecmei_spate_nCovUpdates = 20 # no of times to update cov matrix during simulations
  }

  return(p)
}


