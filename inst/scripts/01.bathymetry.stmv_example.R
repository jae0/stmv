
# Bathymetry


# Spatial interpolation using stmv

# FFT is the method of choice for speed and ability to capture the variability .. this shold take just a few minutes

# 2 mins
scale_ram_required_main_process = 1 # GB twostep / fft ---
scale_ram_required_per_process  = 1 # twostep / fft /fields vario ..  (mostly 0.5 GB, but up to 5 GB)
scale_ncpus = min( parallel::detectCores(), floor( (ram_local()- scale_ram_required_main_process) / scale_ram_required_per_process ) )

# 5 mins
interpolate_ram_required_main_process = 1 # GB twostep / fft
interpolate_ram_required_per_process  = 1.5 # twostep / fft /fields vario ..
interpolate_ncpus = min( parallel::detectCores(), floor( (ram_local()- interpolate_ram_required_main_process) / interpolate_ram_required_per_process ) )


p0 = aegis::spatial_parameters( spatial_domain="bathymetry_example",
  aegis_proj4string_planar_km="+proj=utm +ellps=WGS84 +zone=20 +units=km",
  dres=1/60/4, pres=0.5, lon0=-64, lon1=-62, lat0=44, lat1=45, psignif=2 )
# or:  p0 = stmv_test_data( "aegis.test.parameters")

DATA = list(
  input = stmv::stmv_test_data( datasource="aegis.space", p=p0),
  output = list( LOCS = spatial_grid(p0) )
)

DATA$input = lonlat2planar( DATA$input, p0$aegis_proj4string_planar_km )
DATA$input = DATA$input[, c("plon", "plat", "z")]
DATA$input = DATA$input[ which(is.finite(DATA$input$z)), ]

p = aegis.bathymetry::bathymetry_parameters(
  p=p0,  # start with spatial settings of input data
  project_class="stmv",
  data_root = file.path(work_root, "bathymetry_example"),
  DATA = DATA,
  spatial_domain = p0$spatial_domain,
  inputdata_spatial_discretization_planar_km = p0$pres,  # pres = 0.5
  aegis_dimensionality="space",
  stmv_variables = list(Y="z"),  # required as fft has no formulae
  stmv_global_modelengine = "none",  # too much data to use glm as an entry into link space ... use a direct transformation
  stmv_local_modelengine="fft",
  stmv_fft_filter = "matern tapered lowpass modelled fast_predictions", #  matern with taper, fast predictions are sufficient as data density is high
  stmv_lowpass_nu = 0.5,
  stmv_lowpass_phi = stmv::matern_distance2phi( distance=0.2, nu=0.5, cor=0.1 ),  # note: p$pres = 0.5
  stmv_variogram_method = "fft",
  stmv_autocorrelation_fft_taper = 0.75,  # benchmark from which to taper .. high data density truncate local predictions to capture heterogeneity
  stmv_autocorrelation_localrange = 0.1,
  stmv_autocorrelation_basis_interpolation = c(0.5, 0.25, 0.1, 0.01, 0.001),  # start with smaller localrange estimates and expand
  stmv_filter_depth_m = FALSE,  # need data above sea level to get coastline
  stmv_Y_transform =list(
    transf = function(x) {log10(x + 2500)} ,
    invers = function(x) {10^(x) - 2500}
  ), # data range is from -1667 to 5467 m: make all positive valued
  stmv_rsquared_threshold = 0.01, # lower threshold  .. i.e., ignore ... there is no timeseries model, nor a fixed effect spatial "model"
  stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
  stmv_distance_scale = c( 2.5, 5, 10, 20, 40, 80 ), # km ... approx guesses of 95% AC range
  stmv_distance_basis_interpolation = c( 5, 10, 15, 20, 40, 80  ) , # range of permissible predictions km (i.e 1/2 stats grid to upper limit) .. in this case 5, 10, 20
  stmv_distance_prediction_limits =c( 2.5, 10 ), # range of permissible predictions km (i.e 1/2 stats grid to upper limit based upon data density)
  stmv_nmin = 50,  # min number of data points req before attempting to model in a localized space
  stmv_nmax = 50 * 10, # no real upper bound.. just speed /RAM
  stmv_force_complete_method = "fft",
  stmv_runmode = list(
    scale = rep("localhost", scale_ncpus),
    interpolate = list(
      c1 = rep("localhost", interpolate_ncpus),  # ncpus for each runmode
      c2 = rep("localhost", interpolate_ncpus),  # ncpus for each runmode
      c3 = rep("localhost", max(1, interpolate_ncpus-1)),
      c4 = rep("localhost", max(1, interpolate_ncpus-2)),
      c5 = rep("localhost", max(1, interpolate_ncpus-2))
    ),
    # if a good idea of autocorrelation is missing, forcing via explicit distance limits is an option
    # interpolate_distance_basis = list(
    #   d1 = rep("localhost", interpolate_ncpus),
    #   d2 = rep("localhost", interpolate_ncpus),
    #   d3 = rep("localhost", max(1, interpolate_ncpus-1)),
    #   d4 = rep("localhost", max(1, interpolate_ncpus-1)),
    #   d5 = rep("localhost", max(1, interpolate_ncpus-2)),
    #   d6 = rep("localhost", max(1, interpolate_ncpus-2))
    # ),
    interpolate_force_complete = rep("localhost", max(1, interpolate_ncpus-2)),
    globalmodel = FALSE,
    restart_load = FALSE,
    save_completed_data = TRUE # just a dummy variable with the correct name
  )  # ncpus for each runmode
)

p$spatial_domain_subareas =NULL

if (0) {
  # to force serial mode
   p$stmv_runmode = list(
    scale=rep("localhost", scale_ncpus),
    interpolate = rep("localhost", 1),
    interpolate_force_complete = rep("localhost", max(1, interpolate_ncpus-2)),
    globalmodel = FALSE,
    save_intermediate_results = FALSE,
    save_completed_data = TRUE # just a dummy variable with the correct name
   )

}

# quick look of data
  dev.new(); surface( as.image( Z=DATA$input$z, x=DATA$input[, c("plon", "plat")], nx=p$nplons, ny=p$nplats, na.rm=TRUE) )


stmv( p=p  )  # This will take from a few minutes, depending upon system
# stmv_db( p=p, DS="cleanup.all" )


predictions = stmv_db( p=p, DS="stmv.prediction", ret="mean" )
statistics  = stmv_db( p=p, DS="stmv.stats" )
locations =  spatial_grid( p )

# comparison
dev.new(); surface( as.image( Z=predictions, x=locations, nx=p$nplons, ny=p$nplats, na.rm=TRUE) )

(p$statsvars)
# p$statsvars = c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "phi", "nu", "localrange" )
dev.new(); levelplot( predictions[] ~ locations[,1] + locations[,2], aspect="iso" )
dev.new(); levelplot( statistics[,match("nu", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu
dev.new(); levelplot( statistics[,match("sdTotal", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #sd total
dev.new(); levelplot( statistics[,match("localrange", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #localrange

