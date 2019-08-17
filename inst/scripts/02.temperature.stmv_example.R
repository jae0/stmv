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
p0 = aegis::spatial_parameters( spatial.domain="temperature_test", internal.crs="+proj=utm +ellps=WGS84 +zone=20 +units=km", dres=1/60/4, pres=0.5, lon0=-64, lon1=-62, lat0=44, lat1=45, psignif=2 )
# or:  p0 = stmv_test_data( "aegis.test.paramaters")

DATA = list(
  input = stmv::stmv_test_data( datasource="aegis.spacetime", p=p0),
  output = list( LOCS = spatial_grid(p0) )
)
DATA$input = lonlat2planar( DATA$input, p0$internal.crs )
DATA$input = DATA$input[, c("plon", "plat", "tiyr", "z", "t")]   ## yr, cos.w and sin.w are note required as they are computed internally


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
  stmv_global_modelformula = "none",  # only marginally useful .. consider removing it and use "none",
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
  stmv_fft_filter="matern_tapered",  #  matern, krige (very slow), lowpass, lowpass_matern
  # stmv_fft_taper_fraction = sqrt(0.5),  # in local smoothing convolutions taper to this areal expansion factor sqrt( r=0.5 ) ~ 70% of variance in variogram
  # stmv_lowpass_nu = 0.5,
  # stmv_lowpass_phi = 0.1,  # note: p$pres = 0.2
  stmv_variogram_method = "fft",
  stmv_autocorrelation_fft_taper = 0.5,  # benchmark from which to taper
  stmv_autocorrelation_localrange=0.1,
  stmv_autocorrelation_interpolation = c(0.5, 0.1, 0.05, 0.01),
  stmv_local_model_distanceweighted = TRUE,
  stmv_rsquared_threshold = 0, # lower threshold .. not used if twostep method
  stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
  stmv_distance_scale = c( 20, 30, 40, 50 ), # km ... approx guess of 95% AC range
  stmv_distance_prediction_fraction = 0.95, #
  stmv_nmin = 16*(year.assessment - year.start),  # ~ 1000 min number of data points req before attempting to model timeseries in a localized space .. control no error in local model
  stmv_nmax = 25*(year.assessment - year.start), # no real upper bound.. just speed / RAM limits  .. can go up to 10 GB / core if too large
  stmv_runmode = list(
    scale = rep("localhost", scale_ncpus),
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


# quick look of data
dev.new(); surface( as.image( Z=DATA$input$t, x=DATA$input[, c("plon", "plat")], nx=p$nplons, ny=p$nplats, na.rm=TRUE) )

stmv( p=p )  # This will take from 40-70 hrs, depending upon system


predictions = stmv_db( p=p, DS="stmv.prediction", ret="mean" )
statistics  = stmv_db( p=p, DS="stmv.stats" )
locations   = spatial_grid( p )


# comparisons
dev.new(); surface( as.image( Z=rowMeans(predictions), x=locations, nx=p$nplons, ny=p$nplats, na.rm=TRUE) )

# stats
(p$statsvars)
# p$statsvars = c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "phi", "nu", "localrange" )
dev.new(); levelplot( predictions[,1] ~ locations[,1] + locations[,2], aspect="iso" )
dev.new(); levelplot( statistics[,match("nu", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) # nu
dev.new(); levelplot( statistics[,match("sdTot", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #sd total
dev.new(); levelplot( statistics[,match("localrange", p$statsvars)]  ~ locations[,1] + locations[,2], aspect="iso" ) #localrange


# finished
