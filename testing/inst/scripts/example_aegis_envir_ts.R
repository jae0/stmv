

# exploratory / alternate methods .. examples modelling aegis spatio-temporal variables

# ----------------
# temperature

  year.assessment = lubridate::year(lubridate::now())

  testing.root = file.path( "/home", "jae", "tmp", "stmvdev", "snowcrab" ),

  p = aegis.temperature::temperature_parameters(
      data_root = testing.root,
      # data_root = project.datadirectory( "aegis", "temperature" ),
      spatial_domain = "canada.east", # default
      DATA = 'temperature_db( p=p, DS="stmv_inputs" )',
      additional.data=c("groundfish", "snowcrab", "USSurvey_NEFSC", "lobster"),
      inputdata_spatial_discretization_planar_km = p$pres / 100, # controls resolution of data prior to modelling (km .. ie 100 linear units smaller than the final discretization pres)
      yrs = 1950:year.assessment,
      dimensionality="space-time-cyclic",
      stmv_local_modelengine = "gam" ,
      stmv_local_modelformula = formula( paste(
        't ~ s(yr, k=5, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ',
          '+ s(log(z), k=3, bs="ts") + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") ',
          '+ s(log(z), plon, plat, cos.w, sin.w, yr, k=100, bs="ts")') ),
      stmv_gam_optimizer=c("outer", "bfgs"),
      stmv_global_modelengine = "gam",
      stmv_global_modelformula = formula( t ~ s(z, bs="ts" + s(s.range, bs="ts") + s(dZ, bs="ts") + s(ddZ, bs="ts") + s(log(substrate.grainsize), bs="ts")  ) ), # marginally useful .. consider removing it."none",
      stmv_global_family = gaussian(link="identity"),
      stmv_rsquared_threshold = 0.20, # lower threshold
      stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
      stmv_distance_scale = c(35, 45, 55), # km ... approx guess of 95% AC range

    libs = c("mgcv", "spate"),
   #  stmv_distance_prediction = 4, # this is a half window km
    stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    stmv_local_modelengine = "userdefined",
    stmv_local_modelengine_userdefined = stmvdev::stmv__spate,
    stmv_spate_method = "mcmc_fast",
    stmv_spate_boost_timeseries = TRUE,  # use simple GAM spectral contraint to structure timeseries as spate's fft in time seems to cause overfitting ?
    stmv_spate_nburnin = 1000,
    stmv_spate_nposteriors = 1000,
    stmv_spate_nCovUpdates = 20, # no of times to update cov

      clusters = rep("localhost", ncpus)
    )


  )

  # timings in 2017: stage1 ~35 hrs,  stage2 ~7.5 hrs, stage3 ~15+ hrs

    p = stmv( p=p, runmode=c("initialize", "globalmodel" ) ) # no global_model and force a clean restart

  currentstatus = stmv_statistics_status( p=p )
    parallel_run( stmv_interpolate, p=p,
      runindex=list( locs=currentstatus$todo[sample.int(length( currentstatus$todo ))] ),
      local.n.complete=currentstatus["n.complete"]
    )
    stmv_db( p=p, DS="save_current_state" ) # saved current state (internal format)


    currentstatus = stmv_statistics_status( p=p, reset="incomplete" )
    parallel_run( stmv_interpolate, p=p,
      runindex=list( locs=currentstatus$todo[sample.int(length( currentstatus$todo ))] ),
      local.n.complete=currentstatus["n.complete"]
    )
    stmv_db( p=p, DS="save_current_state" ) # saved current state


   currentstatus = stmv_statistics_status( p=p, reset="incomplete" )
    parallel_run( stmv_interpolate, p=p,
      runindex=list( locs= currentstatus$todo[sample.int(length( currentstatus$todo ))] ),
      local.n.complete=currentstatus["n.complete"],
      stmv_local_modelengine = "tps"
    )
    stmv_db( p=p, DS="save_current_state" )


    stmv_db( p=p, DS="stmv.results" ) # save to disk for use outside stmv*, returning to user scale

    # if (really.finished) stmv_db( p=p, DS="cleanup.all" )




  # 2.  collect predictions from stmv and warp/break into sub-areas defined by
  #     p$spatial_domain_subareas = c( "SSE", "SSE.mpa", "snowcrab" )
  temperature_db( p=p, DS="predictions.redo" ) # 10 min
  temperature_db( p=p, DS="stmv.stats.redo" ) # warp to sub grids

  # 3. extract relevant statistics
  # or parallel_runs: ~ 1 to 2 GB / process .. ~ 4+ hr
  temperature_db( p=p, DS="bottom.statistics.annual.redo" )

  # 4. all time slices in array format
  temperature_db( p=p,  DS="spatial.annual.seasonal.redo" )

  # 5. time slice at prediction time of year
  temperature_db( p=p,  DS="timeslice.redo" )

  # 6. complete statistics and warp/regrid database ... ~ 2 min :: only for  default grid . TODO might as well do for each subregion/subgrid
  temperature_db( p=p, DS="complete.redo")


# 7. maps
  # run only on local cores ... file swapping seem to reduce efficiency
  p = aegis.temperature::temperature_parameters(
    yrs=p$yrs,
    runindex=list(yrs=p$yrs),
    clusters=rep("localhost", detectCores() ) )

  temperature_map( p=p )

  # just redo a couple maps for ResDoc in the  SSE domain
  p$spatial_domain = "SSE"
  #p$bstats = "tmean"
  p = spatial_parameters( p=p )  # default grid and resolution
  p$corners = data.frame(plon=c(150, 1022), plat=c(4600, 5320) )

  temperature_map( p=p, DS='climatology' )
  temperature_map( p=p, DS='annual' )
  temperature_map( p=p, DS='climatology' )
