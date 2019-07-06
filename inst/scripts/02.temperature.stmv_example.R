  # 1. stmv interpolations assuming some seasonal pattern
  # twostep:  ~ 160+ hrs


  if (!exists("year.assessment")) {
    year.assessment=lubridate::year(Sys.Date())
    year.assessment=lubridate::year(Sys.Date()) - 1
  }

  p = aegis.temperature::temperature_parameters( yrs=1950:year.assessment )  # these are default years

    scale_ram_required_main_process = 0.8 # GB twostep / fft
    scale_ram_required_per_process  = 1.25 # twostep / fft /fields vario ..  (mostly 0.5 GB, but up to 5 GB) -- 20 hrs
    scale_ncpus = min( parallel::detectCores(), floor( (ram_local()- scale_ram_required_main_process) / scale_ram_required_per_process ) )

    interpolate_ram_required_main_process = 24 # GB twostep / fft
    interpolate_ram_required_per_process  = 1.25 # 1 GB seems enough for twostep / fft /fields vario .. but make 2 in case
    interpolate_ncpus = min( parallel::detectCores(), floor( (ram_local()- interpolate_ram_required_main_process) / interpolate_ram_required_per_process ) )

    nyrs = year.assessment-1950


    p0 = aegis::spatial_parameters( spatial.domain="temperature_test", internal.crs="+proj=utm +ellps=WGS84 +zone=20 +units=km", dres=1/60/4, pres=0.5, lon0=-64, lon1=-62, lat0=44, lat1=45, psignif=2 )

    DATA = list(
      input = stmv::stmv_test_data( datasource="aegis.spacetime", p=p0),
      output = list( LOCS = spatial_grid(p0) )
    )

    DATA$input = lonlat2planar( DATA$input, p0$internal.crs )
    DATA$input = DATA$input[, c("plon", "plat", "tiyr", "z", "t")]   ## yr, cos.w and sin.w are note required as they are computed internally

    p = aegis.bathymetry::bathymetry_parameters(
      p=p0,
      project.mode="stmv",
      data_root = file.path(work_root, "temperature_test"),
      DATA = DATA,
      spatial.domain = p0$spatial.domain,
      pres_discretization_temperature = 1 / 20, # 1==p$pres; controls resolution of data prior to modelling (km .. ie 20 linear units smaller than the final discretization pres)
      yrs = 1950:year.assessment,
      stmv_dimensionality="space-year-season",
      stmv_global_modelengine = "none",
      stmv_global_modelformula = "none",  # only marginally useful .. consider removing it and use "none",
      stmv_global_family ="none",
      stmv_local_modelengine = "twostep" ,
      stmv_local_modelformula_time = formula( paste(
        't',
        '~ s( yr, k=20, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts")  ',
        '+ s( yr, cos.w, sin.w, k=30, bs="ts") ',
        '+ s( log(z), k=3, bs="ts") + s( plon, k=3, bs="ts") + s( plat, k=3, bs="ts")  ',
        '+ s( log(z), plon, plat, k=30, bs="ts")  '
       ) ),
      stmv_twostep_time = "gam",
      stmv_twostep_space = "fft",  # everything else is too slow ...
      stmv_fft_filter="lowpass_matern_tapered",  #  matern, krige (very slow), lowpass, lowpass_matern
  stmv_fft_taper_fraction = sqrt(0.5),  # in local smoothing convolutions taper to this areal expansion factor sqrt( r=0.5 ) ~ 70% of variance in variogram
  stmv_fft_taper_correlation = 0,  # benchmark from which to taper
  stmv_lowpass_nu = 0.5,
  stmv_lowpass_phi = 0.1,  # note: p$pres = 0.2
  stmv_variogram_method = "fft",
  stmv_variogram_nbreaks = 50,
  stmv_range_correlation=0.1,
  stmv_range_correlation_boostdata = c(0.01, 0.001, 0.0001),
      stmv_local_model_distanceweighted = TRUE,
      stmv_rsquared_threshold = 0, # lower threshold .. not used if twostep method
      stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
      stmv_distance_scale = c( 20, 25, 30, 35, 40, 45, 50 ), # km ... approx guess of 95% AC range
      stmv_distance_prediction_fraction = 4/5, # i.e. 4/5 * 5 = 4 km
      stmv_clusters = list( scale=rep("localhost", scale_ncpus), interpolate=rep("localhost", interpolate_ncpus) ),  # ncpus for each runmode
      stmv_nmin = 16*nyrs,  # ~ 1000 min number of data points req before attempting to model timeseries in a localized space .. control no error in local model
      stmv_nmax = 25*nyrs # no real upper bound.. just speed / RAM limits  .. can go up to 10 GB / core if too large
    )


p$spatial.domain.subareas =NULL

# quick look of data
  dev.new(); surface( as.image( Z=DATA$input$z, x=DATA$input[, c("plon", "plat")], nx=p$nplons, ny=p$nplats, na.rm=TRUE) )


# runmode=c( "globalmodel", "scale", "interpolate", "interpolate_boost", "interpolate_force_complete", "save_completed_data")
runmode=c( "globalmodel", "scale", "interpolate", "save_completed_data")
runmode=c(  "interpolate", "save_completed_data")
runmode=c( "interpolate", "interpolate_boost", "save_completed_data")
# runmode=c( "interpolate", "interpolate_boost", "interpolate_force_complete", "save_completed_data")

stmv( p=p, runmode=runmode )  # This will take from 40-70 hrs, depending upon system


predictions = stmv_db( p=p, DS="stmv.prediction", ret="mean" )
statistics  = stmv_db( p=p, DS="stmv.stats" )
locations   = spatial_grid( p )


# comparison
  dev.new(); surface( as.image( Z=predictions, x=locations, nx=p$nplons, ny=p$nplats, na.rm=TRUE) )

  dev.new(); levelplot( predictions ~ locations[,1] + locations[,2], aspect="iso" )
  dev.new(); levelplot( statistics[,7]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu
  dev.new(); levelplot( statistics[,1]  ~ locations[,1] + locations[,2], aspect="iso" ) #sd total
  dev.new(); levelplot( statistics[,8]  ~ locations[,1] + locations[,2], aspect="iso" ) #range


# finished
