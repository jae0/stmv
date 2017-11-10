

emei_parameters = function( p=NULL  ) {
  # some generic defaults
  if (is.null(p)) stop( "Parameter list is not structured properly" )

  if (!exists("emei_current_status", p))  p$emei_current_status = file.path( p$savedir, "emei_current_status" )

  if (!exists("clusters", p)) p$clusters = rep("localhost", detectCores() )  # default if not given
  if( !exists( "storage.backend", p))  p$storage.backend="bigmemory.ram"
  if( !exists( "emei_variogram_method", p)) p$emei_variogram_method="fast"   # note GP methods are slow when there is too much data
  if (!exists( "emei_global_family", p)) p$emei_global_family = gaussian()
  if (!exists( "emei_eps", p)) p$emei_eps = 0.001  # distance units for eps noise to permit mesh gen for boundaries
  if (!exists( "emei_quantile_bounds", p)) p$emei_quantile_bounds = c(0.01, 0.99) # remove these extremes in interpolations
  if (!exists( "eps", p)) p$eps = 1e-6 # floating point precision
  if (!exists( "boundary", p)) p$boundary = FALSE
  if (!exists( "depth.filter", p)) p$depth.filter = FALSE # depth is given as log(depth) so, choose andy stats locations with elevation > 1 m as being on land
  if (!exists( "emei_kernelmethods_use_all_data", p)) p$emei_kernelmethods_use_all_data =TRUE ## speed and RAM usage improvement is minimal (if any) when off, leave on or remove option and fix as on
  if (!exists( "emei_multiplier_stage2", p) ) p$emei_multiplier_stage2 = c( 1.1, 1.25 ) # distance multiplier for stage 2 interpolations

  # used by "fields" GRMF functions
  if ( p$emei_local_modelengine %in% c("gaussianprocess2Dt", "gaussianprocess" )) {
    if (!exists("phi.grid", p) ) p$phi.grid = 10^seq( -6, 6, by=0.5) * p$emei_distance_scale # maxdist is aprox magnitude of the phi parameter
    if (!exists("lambda.grid", p) ) p$lambda.grid = 10^seq( -9, 3, by=0.5) # ratio of tau sq to sigma sq
  }

  if ( p$emei_local_modelengine %in% c("gam" )) {
    # p$emei_gam_optimizer=c("outer","optim")
    if (!exists("emei_gam_optimizer", p)) p$emei_gam_optimizer=c("outer","bfgs")
    # if (!exists("emei_gam_optimizer", p)) p$emei_gam_optimizer="perf"
  }

  if ( p$emei_local_modelengine %in% c("krige" )) {
     # nothing to add yet ..
  }

  if ( p$emei_local_modelengine %in% c("spate")) {
    if (!exists( "emei_spate_method", p)) p$emei_spate_method = "mcmc_fast"  # "mcmc" is too slow for temperature
    if (!exists( "emei_spate_boost_timeseries", p)) p$emei_spate_boost_timeseries = TRUE 
    if (!exists( "emei_spate_nburnin", p)) p$emei_spate_nburnin =1000
    if (!exists( "emei_spate_nposteriors", p)) p$emei_spate_nposteriors = 4000
    if (!exists( "emei_spate_nCovUpdates", p)) p$emei_spate_nCovUpdates = 20 # no of times to update cov matrix during simulations
  }

  return(p)
}


