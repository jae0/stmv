
# Bathymetry


# Spatial interpolation using stmv

# FFT is the method of choice for speed and ability to capture the variability .. this shold take just a few minutes

scale_ram_required_main_process = 1 # GB twostep / fft
scale_ram_required_per_process  = 1 # twostep / fft /fields vario ..  (mostly 0.5 GB, but up to 5 GB)
scale_ncpus = min( parallel::detectCores(), floor( (ram_local()- scale_ram_required_main_process) / scale_ram_required_per_process ) )

interpolate_ram_required_main_process = 1 # GB twostep / fft
interpolate_ram_required_per_process  = 1.5 # twostep / fft /fields vario ..
interpolate_ncpus = min( parallel::detectCores(), floor( (ram_local()- interpolate_ram_required_main_process) / interpolate_ram_required_per_process ) )


p0 = aegis::spatial_parameters( spatial.domain="bathymetry_example", internal.crs="+proj=utm +ellps=WGS84 +zone=20 +units=km", dres=1/60/4, pres=0.5, lon0=-64, lon1=-62, lat0=44, lat1=45, psignif=2 )

DATA = list(
  input = stmv::stmv_test_data( datasource="aegis.space", p=p0),
  output = list( LOCS = spatial_grid(p0) )
)

DATA$input = lonlat2planar( DATA$input, p0$internal.crs )
DATA$input = DATA$input[, c("plon", "plat", "z")]

p = aegis.bathymetry::bathymetry_parameters(
  p=p0,
  project.mode="stmv",
  data_root = file.path(work_root, "bathymetry_example"),
  DATA = DATA,
  spatial.domain = p0$spatial.domain,
  pres_discretization_bathymetry = p0$pres,
  stmv_dimensionality="space",
  variables = list(Y="z"),  # required as fft has no formulae
  stmv_global_modelengine = "none",  # too much data to use glm as an entry into link space ... use a direct transformation
  stmv_global_modelformula = "none",  # only marginally useful .. consider removing it and use "none",
  stmv_global_family ="none",
  stmv_local_modelengine="fft",
  # stmv_fft_filter = "matern_tapered", #  matern with taper
  stmv_fft_filter = "lowpass_matern_tapered", #  act as a low pass filter first before matern with taper .. depth has enough data for this. Otherwise, use:
  stmv_fft_taper_method = "modelled",  # vs "empirical"
  # stmv_fft_taper_fraction = 0.5,  # if empirical: in local smoothing convolutions taper to this areal expansion factor sqrt( r=0.5 ) ~ 70% of variance in variogram
  stmv_autocorrelation_fft_taper = 0.5,  # benchmark from which to taper
  stmv_lowpass_nu = 0.1,
  stmv_lowpass_phi = matern_distance2phi( distance=0.1, nu=0.1, cor=0.1 ),  # note: p$pres = 0.2
  stmv_variogram_method = "fft",
  stmv_autocorrelation_localrange = 0.1,
  stmv_autocorrelation_interpolation = c(0.25, 0.1, 0.05, 0.01),
  depth.filter = FALSE,  # need data above sea level to get coastline
  stmv_Y_transform =list(
    transf = function(x) {log10(x + 2500)} ,
    invers = function(x) {10^(x) - 2500}
  ), # data range is from -1667 to 5467 m: make all positive valued
  stmv_rsquared_threshold = 0.01, # lower threshold  .. i.e., ignore
  stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
  stmv_distance_scale = c( 10, 20, 30, 40, 50  ), # km ... approx guesses of 95% AC range
  stmv_distance_prediction_fraction = 0.95, # i.e. 4/5 * 5 = 4 km .. relative to stats grid
  stmv_nmin = 200,  # min number of data points5 req before attempting to model in a localized space
  stmv_nmax = 500, # no real upper bound.. just speed /RAM
  stmv_runmode = list(
    # scale = rep("localhost", scale_ncpus),
    interpolate = list(
        cor_0.5 = rep("localhost", interpolate_ncpus),
        cor_0.1 = rep("localhost", interpolate_ncpus),
        cor_0.05 = rep("localhost", max(1, interpolate_ncpus-1)),
        cor_0.01 = rep("localhost", max(1, interpolate_ncpus-2))
      ),  # ncpus for each runmode
    interpolate_force_complete = rep("localhost", max(1, interpolate_ncpus-2)),
    globalmodel = TRUE,
    save_intermediate_results = FALSE,
    save_completed_data = TRUE # just a dummy variable with the correct name
  )  # ncpus for each runmode
)

p$spatial.domain.subareas =NULL

if (0) {
  # to force serial mode
   p$stmv_runmode = list(
    scale=rep("localhost", scale_ncpus),
    interpolate = list(
        cor_0.5 = rep("localhost", 1),
        cor_0.1 = rep("localhost", 1),
        cor_0.05 = rep("localhost", 1),
        cor_0.01 = rep("localhost", 1)
    ),  # ncpus for each runmode
    interpolate_force_complete = rep("localhost", max(1, interpolate_ncpus-2)),
    globalmodel = FALSE,
    save_intermediate_results = FALSE,
    save_completed_data = TRUE # just a dummy variable with the correct name
   )

}

# quick look of data
  dev.new(); surface( as.image( Z=DATA$input$z, x=DATA$input[, c("plon", "plat")], nx=p$nplons, ny=p$nplats, na.rm=TRUE) )


stmv( p=p  )  # This will take from 40-70 hrs, depending upon system


predictions = stmv_db( p=p, DS="stmv.prediction", ret="mean" )
statistics  = stmv_db( p=p, DS="stmv.stats" )
locations   = spatial_grid( p )


# comparison
  dev.new(); surface( as.image( Z=predictions, x=locations, nx=p$nplons, ny=p$nplats, na.rm=TRUE) )

  dev.new(); levelplot( predictions ~ locations[,1] + locations[,2], aspect="iso" )
  dev.new(); levelplot( statistics[,7]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu
  dev.new(); levelplot( statistics[,1]  ~ locations[,1] + locations[,2], aspect="iso" ) #sd total
  dev.new(); levelplot( statistics[,8]  ~ locations[,1] + locations[,2], aspect="iso" ) #range
