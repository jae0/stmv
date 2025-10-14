

# STMV analysis of cod

groundfish_species_code = 10
yrs = 1970:2017

ncpus = parallel::detectCores()
if (0) {
  ram_required_main_process = 6 # GB
  ram_required_per_process  = 1  # about 1.2GB on average ..in 2018, for twostep / fft
  ncpu = min( parallel::detectCores(), trunc( (ram_local()-ram_required_main_process) / ram_required_per_process ) )
}

p = aegis_stmv( DS="parameters",
  project_name = "cod_stmv",
  data_root = project.datadirectory( "carstm", "cod_stmv" ),
  stmv_variables = list(Y="totno"),
  groundfish_species_code=groundfish_species_code,
  yrs = yrs,
  DATA = 'aegis_stmv( p=p, DS="stmv_inputs" )',
  aegis_project_datasources = c("speciescomposition"),
  stmv_global_modelengine ="gam",
  stmv_global_family = poisson(link="log"),
  stmv_global_modelformula = formula( paste(
    ' totno ~ offset(log(data_offset)) ',  # CARSTM does log-transformation internally  but stmv does not
    ' + s( t, k = 3, bs = "ts") + s( tsd, k = 3, bs = "ts") + s( tmax, k = 3, bs = "ts") + s( degreedays, k = 3, bs = "ts") ',
    ' + s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts") + s( log(ddZ), k=3, bs="ts") ',
    # ' + s(pca1, k=3, bs="ts") + s(pca2, k=3, bs="ts") ',
    ' + s( log(substrate.grainsize), k=3, bs="ts")    '
    )),
  stmv_local_modelengine = "twostep",
  stmv_twostep_time = "gam",
  stmv_twostep_space = "fft",
  stmv_fft_filter="matern",  #  matern, krige (very slow), lowpass, lowpass_matern
  stmv_gam_optimizer=c("outer", "bfgs") ,
  stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** ),
  stmv_distance_scale = c( 40, 50, 60 ), #likely must be over 30km, so 50 +/- 20km, should likely match the setting in ~ line 256
  stmv_clusters = list( rep("localhost", ncpus), rep("localhost", ncpus), rep("localhost", ncpus) )  # no of cores used made explicit.. must be same length as "stmv_distance_scale"
)

p$selection=list(
  type="number",
  biologicals=list(
    spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code )
  ),
  survey=list(
    data.source = "groundfish",
    yr = p$yrs,  # time frame for comparison specified above
    months=6:8,  # "summer"
    # dyear = c(150,250)/365, # alternate way of specifying season: summer = which( (x>150) & (x<250) ) , spring = which(  x<149 ), winter = which(  x>251 )
    settype = 1, 
    gear = c("Western IIA trawl", "Yankee #36 otter trawl"),
    polygon_enforce=TRUE,  # make sure mis-classified stations or incorrectly entered positions get filtered out
    ranged_data = c("dyear")  # not used .. just to show how to use range_data
  )
)


  if (0) {  # model testing
      require(mgcv)
      o = aegis_stmv(p=p, DS="stmv_inputs" )  # create fields for
      B = o$input
      range( B$totno )
      # p$stmv_global_modelformula = formula( t ~ s(z, bs="ts") + s(s.range, bs="ts") + s(dZ, bs="ts") + s(ddZ, bs="ts")  ) # marginally useful .. consider removing it and use "none",
      # p$stmv_global_family = gaussian(link="identity")
      p = stmv_variablelist(p=p)  # decompose into covariates, etc
      ii = which( is.finite (rowSums(B[ , c(p$stmv_variables$Y, p$stmv_variables$COV) ])) )
    #  wgts = 1/B$b.sdTotal[ii]
    #  wgts = wgts / mean(wgts)  # makes loglik constant
      global_model = try( gam( formula=p$stmv_global_modelformula, data=B[ii,],
          optimizer= p$stmv_gam_optimizer, family=p$stmv_global_family, na.action=na.omit  ) )
      summary( global_model )
      plot(global_model, all.terms=TRUE, trans=bio.snowcrab::inverse.logit, seWithMean=TRUE, jit=TRUE, rug=TRUE )
  }


#--------------------------------------------------------------
# Run the process
#--------------------------------------------------------------

stmv( p=p, runmode=c("globalmodel", "interpolate" )  ) #  for a clean start

aegis_stmv( p=p, DS="predictions.redo" ) # warp predictions to other grids (if any)
aegis_stmv( p=p, DS="stmv.stats.redo" ) # warp stats to other grids (if any)
aegis_stmv( p=p, DS="complete.redo" )
aegis_stmv( p=p, DS="baseline.redo" )
aegis_stmv( p=p, DS="map.all" )

global_model = stmv_global_model( p=p, DS="global_model")
summary( global_model )

par(mar=c(1,1,1,1)) #change plot margins for Rstudio
plot(global_model)



#
