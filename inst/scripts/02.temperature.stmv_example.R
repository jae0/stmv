# 1. stmv interpolations assuming some seasonal pattern
# twostep:  ~ 160+ hrs


year.assessment = 2018
year.start = 1950


scale_ram_required_main_process = 0.8 # GB twostep / fft
scale_ram_required_per_process  = 1.25 # twostep / fft /fields vario ..  (mostly 0.5 GB, but up to 5 GB) -- 20 hrs
scale_ncpus = min( parallel::detectCores(), floor( (ram_local()- scale_ram_required_main_process) / scale_ram_required_per_process ) )

interpolate_ram_required_main_process = 24 # GB twostep / fft
interpolate_ram_required_per_process  = 1.25 # 1 GB seems enough for twostep / fft /fields vario .. but make 2 in case
interpolate_ncpus = min( parallel::detectCores(), floor( (ram_local()- interpolate_ram_required_main_process) / interpolate_ram_required_per_process ) )


# p = aegis.temperature::temperature_parameters( yrs=yrs )  # these are default years

# prep input data:
p0 = aegis::spatial_parameters( spatial.domain="temperature_test",
  internal.crs="+proj=utm +ellps=WGS84 +zone=20 +units=km",
  dres=1/60/4, pres=0.5, lon0=-64, lon1=-62, lat0=44, lat1=45, psignif=2 )
# or:  p0 = stmv_test_data( "aegis.test.paramaters")

DATA = list(
  input = stmv::stmv_test_data( datasource="aegis.spacetime", p=p0),
  output = list(
    LOCS = spatial_grid(p0),
    COV  = list(
      z = stmv::stmv_test_data( datasource="aegis.bathymetry", p=p0)
    )
  )
)
DATA$input = lonlat2planar( DATA$input, p0$internal.crs )
DATA$input = DATA$input[, c("plon", "plat", "tiyr", "z", "t")]   ## yr, cos.w and sin.w are note required as they are computed internally


# quick look of data
dev.new(); surface( as.image( Z=DATA$input$t, x=DATA$input[, c("plon", "plat")], nx=p$nplons, ny=p$nplats, na.rm=TRUE) )


p = aegis.temperature::temperature_parameters(
  p=p0,  # start with spatial settings of input data
  project.mode="stmv",
  data_root = file.path(work_root, "temperature_test"),
  DATA = DATA,
  spatial.domain = p0$spatial.domain,
  spatial.domain.subareas =NULL,  # prevent subgrid estimation
  pres_discretization_temperature = 1 / 10, # 1==p$pres; controls resolution of data prior to modelling (km .. ie 20 linear units smaller than the final discretization pres)
  yrs = year.start:year.assessment,
  stmv_dimensionality="space-year-season",
  stmv_global_modelengine = "none",
  stmv_local_modelengine = "twostep" ,
  stmv_local_modelformula_time = formula( paste(
    't',
    '~ s( yr, k=30, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts")  ',
    '+ s( yr, cos.w, sin.w, k=20, bs="ts") ',
    '+ s( log(z), k=3, bs="ts") + s( plon, k=3, bs="ts") + s( plat, k=3, bs="ts")  ',
    '+ s( log(z), plon, plat, k=20, bs="ts")  '
    ) ),
  stmv_twostep_time = "gam",
  stmv_twostep_space = "fft",  # everything else is too slow ...
  stmv_fft_filter="matern_tapered",  #  matern, krige (very slow), lowpass, lowpass_matern
  # stmv_fft_taper_fraction = sqrt(0.5),  # when stmv_fft_taper_method = "empirical", in local smoothing convolutions taper to this areal expansion factor sqrt( r=0.5 ) ~ 70% of variance in variogram
  # stmv_lowpass_nu = 0.5,
  # stmv_lowpass_phi = 0.5/5,  # note: p$pres = 0.5
  stmv_variogram_method = "fft",
  stmv_fft_taper_method = "modelled",
  stmv_autocorrelation_fft_taper = 0.5,  # benchmark from which to taper
  stmv_autocorrelation_localrange = 0.1,  # for reporting
  stmv_autocorrelation_interpolation = c( 0.25, 0.1, 0.05, 0.01 ),  # range finding
  stmv_local_model_distanceweighted = TRUE,
  depth.filter = 10, # the depth covariate is input as units of depth (m) so, choose stats locations with elevation > 10m as being on land
  stmv_rsquared_threshold = 0.5, # lower threshold for timeseries model
  stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
  stmv_distance_scale = c( 5, 10, 20, 30, 40, 50 ), # km ... approx guess of 95% AC range
  stmv_distance_prediction_fraction = 0.95, #
  stmv_nmin = 200,  # min number of data points req before attempting to model in a localized space .. control no error in local model
  stmv_nmax = 400, # no real upper bound.. just speed / RAM limits  .. can go up to 10 GB / core if too large
  stmv_tmin = floor( (year.assessment - year.start) * 1 ),
  stmv_force_complete_method = "linear",
  stmv_runmode = list(
    scale = rep("localhost", scale_ncpus),
    interpolate = list(
        cor_0.25 = rep("localhost", interpolate_ncpus),
        cor_0.1 = rep("localhost", interpolate_ncpus),
        cor_0.05 = rep("localhost", max(1, interpolate_ncpus-1)),
        cor_0.01 = rep("localhost", max(1, interpolate_ncpus-2))
      ),  # ncpus for each runmode
    interpolate_force_complete = rep("localhost", max(1, interpolate_ncpus-2)),
    globalmodel = FALSE,
    restart_load = FALSE,  # FALSE means redo all, TRUE means update currently saved  instance
    save_completed_data = TRUE # just a dummy variable with the correct name
  )  # ncpus for each runmode
)

# DATA = NULL

stmv( p=p )  # This will take from xx hrs, depending upon system

# stmv_db( p=p, DS="cleanup.all" )

predictions = stmv_db( p=p, DS="stmv.prediction", ret="mean", yr = year.assessment )
statistics  = stmv_db( p=p, DS="stmv.stats" )
locations   = spatial_grid( p )


# stats
(p$statsvars)
# p$statsvars = c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "phi", "nu", "localrange" )
dev.new(); levelplot( predictions[,1] ~ locations[,1] + locations[,2], aspect="iso" )


dev.new(); levelplot( statistics[,match("nu", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu
dev.new(); levelplot( statistics[,match("sdTotal", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #sd total
dev.new(); levelplot( statistics[,match("localrange", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #localrange


for (i in p$yrs) {
  predictions = stmv_db( p=p, DS="stmv.prediction", ret="mean", yr =i )
  for (j in 1:p$nw) print( levelplot( predictions[,j] ~ locations[,1] + locations[,2], aspect="iso" ) )
}


# comparisons
dev.new(); surface( as.image( Z=rowMeans(predictions), x=locations, nx=p$nplons, ny=p$nplats, na.rm=TRUE) )

# finished
