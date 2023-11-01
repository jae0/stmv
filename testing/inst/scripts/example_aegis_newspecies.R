
# this is an example to demonstrate usage of aegis* and stmv* projects
# to model and infer the abundance and distribution of snow crab in Atlantic Canada (Maritimes Region)
# it is equivalent the the methods in bio.snowcrab/inst/scripts/03.abundance.estimation.r but
# shows a more generalizable and flexible method of inference and estimation

# create alternate, new project data directory for outputs
# NOTE storage of stmv outputs are  in p$stmvSaveDir = stmvSaveDir = file.path( p$modeldir, p$stmv_model_label, p$project_class, paste(  p$stmv_global_modelengine, stmv_local_modelengine, sep="_"), p$stmv_variables$Y, p$spatial_domain)


year.assessment = lubridate::year(lubridate::now()

p = bio.snowcrab::snowcrab_parameters(
  project_class="stmv",
  DATA = 'snowcrab_stmv( p=p, DS="stmv_inputs" )',
  data_root = file.path( "/home", "jae", "bio.data", "stmvdev", "snowcrab" ) # or more simply project.datadirectory("stmvdev","snowcrab")
  yrs=1999:year.assessment,
  stmv_variables=list(Y="snowcrab.large.males_abundance"),
  selection=list(
    type = "abundance",
    biologicals =list(
      sex=0, # male
      mat=1, # do not use maturity status in groundfish data as it is suspect ..
      spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
      len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
      ranged_data="len"
    ),
    survey = list( drop.unreliable.zeros.groundfish.data=TRUE) # esp from 1970 to 1999 measurement of invertebrates was sporatic .. zero-values are dropped as they are unreliable
  ),
  stmv_global_family = gaussian(link=log),
  dimensionality="space-time",
  libs = c("mgcv", "spate"),
  stmv_local_modelengine = "userdefined",
  stmv_local_modelengine_userdefined = stmvdev::stmv__spate,
  stmv_spate_method = "mcmc_fast",
  stmv_spate_boost_timeseries = TRUE,  # use simple GAM spectral contraint to structure timeseries as spate's fft in time seems to cause overfitting ?
  stmv_spate_nburnin = 1000,
  stmv_spate_nposteriors = 1000,
  stmv_spate_nCovUpdates = 20, # no of times to update cov
  stmv_distance_statsgrid = 2, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
  stmv_distance_scale = c(50, 60, 70)
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




# o = snowcrab_stmv(p=p, DS="stmv_inputs" )  # create fields for
snowcrab_stmv( p=p, DS="predictions.redo" ) # warp predictions to other grids
snowcrab_stmv( p=p, DS="stmv.stats.redo" ) # warp stats to other grids
snowcrab_stmv( p=p, DS="complete.redo" )
snowcrab_stmv( p=p, DS="baseline.redo" )
snowcrab_stmv( p=p, DS="map.all" )

global_model = stmv_global_model( p=p, DS="global_model")
summary( global_model )
plot(global_model)






p = bio.snowcrab::snowcrab_parameters( p=p, project_class="stmv",
  stmv_variables=list(Y="snowcrab.large.males_presence_absence"),
  selection=list(
    type = "presence_absence",
    biologicals=list(
      sex=0, # male
      mat=1, # do not use maturity status in groundfish data as it is suspect ..
      spec_bio=bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=2526 ),
      len= c( 95, 200 )/10, #  mm -> cm ; aegis_db in cm
      ranged_data="len"
    ),
    survey=list( drop.unreliable.zeros.groundfish.data=TRUE) # esp from 1970 to 1999 measurement of invertebrates was sporatic .. zero-values are dropped as they are unreliable
  ),
  DATA = 'snowcrab_stmv( p=p, DS="stmv_inputs" )',
  dimensionality="space-time",
  stmv_global_family = binomial(),

  stmv_local_modelengine = "twostep",
  # stmv_local_modelformula = formula( paste(
  #   ' snowcrab.large.males_abundance', '~ s(yr, k=10, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ',
  #   ' + s(cos.w, sin.w, yr, bs="ts", k=20) ',
  #   ' + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=20, bs="ts") ' ) ),
  # stmv_twostep_space = "krige",  # other possibilities: fft", "tps"

  stmv_twostep_space = "gam", #  fft, krige (very slow), lowpass, lowpass_fft
  stmv_local_modelformula_space = formula( paste(
    'snowcrab.large.males_abundance', '~ s(log(z), k=3, bs="ts") + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s( log(z), plon, plat, k=27, bs="ts")  ') ),

  stmv_twostep_time = "gam",
  stmv_local_modelformula_time = formula( paste(
    'snowcrab.large.males_abundance', ' ~ s(yr, k=10, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts")  ',
    '+ s( cos.w, sin.w, yr, k=45, bs="ts")') ),

  stmv_gam_optimizer=c("outer", "bfgs"),
  stmv_distance_statsgrid = 2, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** ),
  stmv_distance_scale = c(50, 60, 70)
)


# o = snowcrab_stmv(p=p, DS="stmv_inputs" )  # create fields for

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


snowcrab_stmv( p=p, DS="predictions.redo" ) # warp predictions to other grids
snowcrab_stmv( p=p, DS="stmv.stats.redo" ) # warp stats to other grids
snowcrab_stmv( p=p, DS="complete.redo" )
snowcrab_stmv( p=p, DS="baseline.redo" )
snowcrab_stmv( p=p, DS="map.all" )

global_model = stmv_global_model( p=p, DS="global_model")
summary( global_model )
plot(global_model)


