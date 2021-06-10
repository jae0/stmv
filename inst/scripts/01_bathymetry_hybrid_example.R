
#########################################
### example 2 -- hybrid stmv/car model:

## :: main difference is:  project_class="hybrid"  (vs "stmv" above)

require(aegis.bathymetry)

scale_ram_required_main_process = 1 # GB twostep / fft ---
scale_ram_required_per_process  = 1 # twostep / fft /fields vario ..  (mostly 0.5 GB, but up to 5 GB)
scale_ncpus = min( parallel::detectCores(), floor( (ram_local()- scale_ram_required_main_process) / scale_ram_required_per_process ) )

# 5 mins
interpolate_ram_required_main_process = 1 # GB twostep / fft
interpolate_ram_required_per_process  = 1.5 # twostep / fft /fields vario ..
interpolate_ncpus = min( parallel::detectCores(), floor( (ram_local()- interpolate_ram_required_main_process) / interpolate_ram_required_per_process ) )


p0 = aegis::spatial_parameters( spatial_domain="bathymetry_example",
  aegis_proj4string_planar_km="+proj=utm +ellps=WGS84 +zone=20 +units=km",
  dres=1/60/4, pres=1, lon0=-64, lon1=-62, lat0=44, lat1=45, psignif=2 )

# or:
p0 = stmv_test_data( "aegis.test.parameters")

input = stmv::stmv_test_data( datasource="aegis.space", p=p0)
input = sf::st_as_sf( input, coords=c("lon","lat"), crs=st_crs(projection_proj4string("lonlat_wgs84")) )
input = sf::st_transform( input, crs=st_crs(p0$aegis_proj4string_planar_km) )
input = as.data.frame( cbind( input$z, st_coordinates(input) ) )
names(input) = c("z", "plon", "plat")
input = input[ which(is.finite(input$z)), ]
output = list( LOCS = spatial_grid(p0) )
DATA = list( input = input,  output = output )

input = output = NULL
gc()

p = bathymetry_parameters(
  p=p0,  # start with spatial settings of input data
  project_class="hybrid",
  stmv_model_label="carstm_statsgrid10km",  # this label is used as directory for storage
  data_root = file.path(tempdir(), "bathymetry_example"),
  DATA = DATA,
  spatial_domain = p0$spatial_domain,
  spatial_domain_subareas =NULL,
  inputdata_spatial_discretization_planar_km = p0$pres/5,  # pres = 1, pres defines underlaying lattice resolution 
  aegis_dimensionality="space",
  stmv_variables = list(Y="z"),  # required as fft has no formulae
  stmv_global_modelengine = "none",  # too much data to use glm as an entry into link space ... use a direct transformation
  stmv_local_modelengine="carstm",
  # stmv_au_distance_reference="completely_inside_boundary",
  stmv_au_distance_reference = "none", 
  stmv_au_buffer_links = 1, # number of additional neighbours to extend beyond initial solution
  stmv_filter_depth_m = FALSE,  # need data above sea level to get coastline
  stmv_Y_transform =list(
    transf = function(x) {log10(x + 2500)} ,
    invers = function(x) {10^(x) - 2500}
  ), # data range is from -1667 to 5467 m: make all positive valued
  stmv_rsquared_threshold = 0.01, # lower threshold  .. i.e., ignore ... there is no timeseries model, nor a fixed effect spatial "model"
  stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
  # stmv_distance_scale = c( 2.5, 5, 10, 20, 40, 60, 80 ), # km ... approx guesses of 95% AC range
  stmv_distance_prediction_limits =c( 5, 10 ), # range of permissible predictions km (i.e 1/2 stats grid to upper limit based upon data density)
  stmv_nmin = 50,  # min number of data points req before attempting to model in a localized space
  stmv_nmax = 600, # no real upper bound.. just speed /RAM
  stmv_force_complete_method = "linear_interp",
  stmv_runmode = list(
    carstm = rep("localhost", 1),
    globalmodel = FALSE,
    restart_load = FALSE,
    save_completed_data = TRUE # just a dummy variable with the correct name
  )
)


if (0) {
  # to force parallel mode
   p$stmv_runmode = list(
    carstm = rep("localhost", scale_ncpus),
    globalmodel = FALSE,
    restart_load = FALSE,
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


statsvars = dimnames(statistics)[[2]]

# statsvars = c("sdTotal", "ndata", "fixed_mean", "fixed_sd", "dic", "dic_p_eff",
#   "waic", "waic_p_eff", "mlik", "Expected_number_of_parameters",
#   "Stdev_of_the_number_of_parameters", "Number_of_equivalent_replicates",
#   "Precision_for_the_Gaussian_observations", "Precision_for_aui",
#   "Phi_for_aui", "Precision_for_the_Gaussian_observations_sd", "Precision_for_aui_sd", "Phi_for_aui_sd"
# )

# statsvars = c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "phi", "nu", "localrange" )
dev.new(); levelplot( predictions[] ~ locations[,1] + locations[,2], aspect="iso" )
dev.new(); levelplot( statistics[,match("Phi_for_space", statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu
dev.new(); levelplot( statistics[,match("sdTotal", statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #sd total
dev.new(); levelplot( statistics[,match("rsquared", statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #localrange

