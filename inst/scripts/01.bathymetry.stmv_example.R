
# Bathymetry


# Spatial interpolation using stmv
# Total "superhighres": 2-5 GB/process and 4 GB in parent for fft
# gam method requires more ~ 2X
# boundary def takes too long .. too much data to process -- skip
# "highres": ~ 20 hr with 8, 3.2 Ghz cpus on thoth using fft method jc: 2016 or 2~ 6 hr on hyperion
# "superhighres" fft: looks to be the best in performance/quality; req ~5 GB per process req
# FFT is the method of choice for speed and ability to capture the variability
# krige method is a bit too oversmoothed, especially where rapid changes are occuring


scale_ram_required_main_process = 1 # GB twostep / fft
scale_ram_required_per_process  = 1 # twostep / fft /fields vario ..  (mostly 0.5 GB, but up to 5 GB)
scale_ncpus = min( parallel::detectCores(), floor( (ram_local()- scale_ram_required_main_process) / scale_ram_required_per_process ) )

interpolate_ram_required_main_process = 1 # GB twostep / fft
interpolate_ram_required_per_process  = 1 # twostep / fft /fields vario ..
interpolate_ncpus = min( parallel::detectCores(), floor( (ram_local()- interpolate_ram_required_main_process) / interpolate_ram_required_per_process ) )

p0 = aegis::spatial_parameters( spatial.domain="testing", internal.crs="+proj=utm +ellps=WGS84 +zone=20 +units=km", dres=1/60/4, pres=0.5, lon0=-64, lon1=-62, lat0=44, lat1=45, psignif=2 )

DATA = list(
  input = stmv::stmv_test_data( datasource="aegis.space", p=p0),
  output = list( LOCS = spatial_grid(p0) )
)

DATA$input = lonlat2planar( DATA$input, p0$internal.crs )
DATA$input = DATA$input[, c("plon", "plat", "z")]

p = aegis.bathymetry::bathymetry_parameters(
  p=p0,
  project.mode="stmv",
  data_root = file.path(work_root, "test"),
  DATA = DATA,
  spatial.domain = p0$spatial.domain,
  pres_discretization_bathymetry = p0$pres,
  stmv_dimensionality="space",
  variables = list(Y="z"),  # required as fft has no formulae
  stmv_global_modelengine = "none",  # too much data to use glm as an entry into link space ... use a direct transformation
  stmv_global_modelformula = "none",  # only marginally useful .. consider removing it and use "none",
  stmv_global_family ="none",
  stmv_local_modelengine="fft",
  stmv_fft_filter = "matern_tapered", #  act as a low pass filter first before matern .. depth has enough data for this. Otherwise, use:
  stmv_lowpass_nu = 0.5,
  stmv_lowpass_phi = 0.1,  # note: p$pres = 0.2
  stmv_fft_taper_factor = 5,  # in local smoothing convolutions occur of this correlation scale
  stmv_variogram_method = "fft",
  stmv_variogram_nbreaks = 32,
  depth.filter = FALSE,  # need data above sea level to get coastline
  stmv_Y_transform =list(
    transf = function(x) {log10(x + 2500)} ,
    invers = function(x) {10^(x - 2500)}
  ), # data range is from -1667 to 5467 m: make all positive valued
  stmv_rsquared_threshold = 0.75, # lower threshold
  stmv_distance_statsgrid = 3, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
  stmv_distance_scale = c(5, 10, 20, 30, 40, 50, 60), # km ... approx guess of 95% AC range
  stmv_distance_prediction_fraction = 4/5, # i.e. 4/5 * 5 = 4 km
  stmv_nmin = 750,  # min number of data points req before attempting to model in a localized space
  stmv_nmax = 1000, # no real upper bound.. just speed
  stmv_clusters = list( scale=rep("localhost", scale_ncpus), interpolate=rep("localhost", interpolate_ncpus) )  # ncpus for each runmode
)

p$spatial.domain.subareas =NULL


# runmode=c( "globalmodel", "scale", "interpolate", "interpolate_boost", "interpolate_force_complete", "save_completed_data")
# runmode=c( "interpolate", "interpolate_boost", "save_completed_data")
stmv( p=p, runmode=runmode )  # This will take from 40-70 hrs, depending upon system


predictions = stmv_db( p=p, DS="load_results", ret="predictions" )
statistics  = stmv_db( p=p, DS="load_results", ret="statistics" )
locations   = spatial_grid(p )


dev.new(); levelplot( predictions ~ locations[,1] + locations[,2], aspect="iso" )
dev.new(); levelplot( statistics[,7]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu
dev.new(); levelplot( statistics[,1]  ~ locations[,1] + locations[,2], aspect="iso" ) #sd total
dev.new(); levelplot( statistics[,8]  ~ locations[,1] + locations[,2], aspect="iso" ) #range

u = as.image( Z=DATA$input$z, x=DATA$input[, c("plon", "plat")], nx=p$nplons, ny=p$nplats, na.rm=TRUE)
surface(u)

require(MBA)
mba.int <- mba.surf( cbind( c(u$z), expand.grid(u$x, u$y) ), 100, 100, extend=TRUE)$xyz.est
image(mba.int, xaxs="r", yaxs="r")
