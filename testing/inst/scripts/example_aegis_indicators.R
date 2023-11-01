

# exploratory / alternate methods .. examples using aegis idnicators

# -----------------------------
# speciescomposition
  year.assessment=lubridate::year(Sys.Date())
  # year.assessment=lubridate::year(Sys.Date()) -1
  p = aegis.speciescomposition::speciescomposition_parameters( yrs=1999:year.assessment )

  # varname="pca1"
  # varname="pca2"

  p$varstomodel = c("pca1", "pca2")
  for ( varname in p$varstomodel) {
    print(varname)

    # ~ 5GB / process for spate
    p = aegis.speciescomposition::speciescomposition_parameters(
      yrs=1999:year.assessment,
      dimensionality="space-time",
      project_class = "stmv",
      stmv_variables=list(Y=varname),
      libs = c("mgcv", "spate", "stmvdev"),
      stmv_local_modelengine = "userdefined",
      stmv_local_modelengine_userdefined = stmvdev::stmv__spate,
      stmv_spate_method = "mcmc_fast",
      stmv_distance_statsgrid = 5, # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
      stmv_distance_scale = c(75, 90, 110),
      stmv_spate_boost_timeseries = TRUE,  # use simple GAM spectral contraint to structure timeseries as spate's fft in time seems to cause overfitting ?
      stmv_gam_optimizer=c("outer", "bfgs")  # .. spate uses local gam-based boost of ts, as in 'twostep'
      stmv_local_modelformula = formula( paste(
        varname, '~ s(yr, k=5, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ',
        ' + s(cos.w, sin.w, yr, bs="ts", k=16) ',
        ' + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=16, bs="ts") ' ) ),
        # similar to GAM model .. this is the seed model
      stmv_local_model_distanceweighted = TRUE,
        clusters=rep("localhost", 4),
      stmv_spate_nburnin = 1000,
      stmv_spate_nposteriors = 1000,
      stmv_spate_nCovUpdates = 20 # no of times to update cov matrix during simulations
    )

  # stmv_local_modelengine = "twostep",
  # stmv_local_modelformula = formula( paste(
  #   ' snowcrab.large.males_abundance', '~ s(yr, k=10, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts") ',
  #   ' + s(cos.w, sin.w, yr, bs="ts", k=20) ',
  #   ' + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=20, bs="ts") ' ) ),
  # # stmv_twostep_space = "krige",  # other possibilities: "fft", "tps"

  # stmv_twostep_space = "gam", #  fft, krige (very slow), lowpass, lowpass_fft
  # stmv_local_modelformula_space = formula( paste(
  #   'snowcrab.large.males_abundance', '~ s(log(z), k=3, bs="ts") + s(plon, k=3, bs="ts") + s(plat, k=3, bs="ts") + s( log(z), plon, plat, k=27, bs="ts")  ') ),

  # stmv_twostep_time = "gam",
  # stmv_local_modelformula_time = formula( paste(
  #   'snowcrab.large.males_abundance', ' ~ s(yr, k=10, bs="ts") + s(cos.w, k=3, bs="ts") + s(sin.w, k=3, bs="ts")  ',
  #   '+ s( cos.w, sin.w, yr, k=45, bs="ts")') ),


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


    aegis_db( p=p, DS="predictions.redo" ) # warp predictions to other grids
    aegis_db( p=p, DS="stmv.stats.redo" ) # warp stats to other grids
    aegis_db( p=p, DS="complete.redo" )
    aegis_db( p=p, DS="baseline.redo" )
    aegis_db_map( p=p )

    if (0) {
      global_model = stmv_global_model( p=p, DS="global_model")
      summary( global_model )
      plot(global_model)
    }
  }


# ----------------

# ca1 :

Family: gaussian
Link function: identity

Formula:
ca1 ~ s(t, k = 3, bs = "ts") + s(tsd.climatology, k = 3, bs = "ts") +
    s(tmean.climatology, k = 3, bs = "ts") + s(log(t.range),
    k = 3, bs = "ts") + s(log(b.range), k = 3, bs = "ts") + s(log(z),
    k = 3, bs = "ts") + s(log(dZ), k = 3, bs = "ts") + s(log(ddZ),
    k = 3, bs = "ts") + s(log(substrate.grainsize), k = 3, bs = "ts")

Parametric coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.025523   0.004243   6.015 1.83e-09 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Approximate significance of smooth terms:
                             edf Ref.df       F  p-value
s(t)                       1.998      2  274.88  < 2e-16 ***
s(tsd.climatology)         1.976      2   77.31  < 2e-16 ***
s(tmean.climatology)       2.000      2 7016.86  < 2e-16 ***
s(log(t.range))            1.396      2   68.72  < 2e-16 ***
s(log(b.range))            1.969      2   21.13 4.59e-10 ***
s(log(z))                  1.983      2  112.45  < 2e-16 ***
s(log(dZ))                 1.888      2   14.47 1.98e-07 ***
s(log(ddZ))                1.856      2   12.15 2.05e-06 ***
s(log(substrate.grainsize)) 1.752      2  274.56  < 2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

R-sq.(adj) =  0.642   Deviance explained = 64.2%
GCV = 0.37765  Scale est. = 0.37733   n = 20956
---


# -------------------
#  varname ="ca2"

Family: gaussian
Link function: identity

Formula:
ca2 ~ s(t, k = 3, bs = "ts") + s(tsd.climatology, k = 3, bs = "ts") +
    s(tmean.climatology, k = 3, bs = "ts") + s(log(t.range),
    k = 3, bs = "ts") + s(log(b.range), k = 3, bs = "ts") + s(log(z),
    k = 3, bs = "ts") + s(log(dZ), k = 3, bs = "ts") + s(log(ddZ),
    k = 3, bs = "ts") + s(log(substrate.grainsize), k = 3, bs = "ts")

Parametric coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) -0.11437    0.00408   -28.1   <2e-16

Approximate significance of smooth terms:
                            edf Ref.df       F p-value
s(t)                       1.58      2  313.49 < 2e-16
s(tsd.climatology)         1.99      2  280.46 < 2e-16
s(tmean.climatology)       2.00      2   81.55 < 2e-16
s(log(t.range))            1.99      2   62.90 < 2e-16
s(log(b.range))            1.50      2   12.85 2.7e-07
s(log(z))                  2.00      2 3337.43 < 2e-16
s(log(dZ))                 1.90      2    5.53  0.0029
s(log(ddZ))                1.98      2   30.29 4.3e-14
s(log(substrate.grainsize)) 1.99      2  140.98 < 2e-16

R-sq.(adj) =  0.748   Deviance explained = 74.8%
GCV = 0.22366  Scale est. = 0.22336   n = 13444



# -------------------
#  varname ="pca1"   ; # July 2017 results follow
Family: gaussian
Link function: identity

Formula:
pca1 ~ s(t, k = 3, bs = "ts") + s(tsd.climatology, k = 3, bs = "ts") +
    s(tmean.climatology, k = 3, bs = "ts") + s(log(t.range),
    k = 3, bs = "ts") + s(log(b.range), k = 3, bs = "ts") + s(log(z),
    k = 3, bs = "ts") + s(log(dZ), k = 3, bs = "ts") + s(log(ddZ),
    k = 3, bs = "ts") + s(log(substrate.grainsize), k = 3, bs = "ts")


Parametric coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept) 0.124972   0.000809     154   <2e-16

Approximate significance of smooth terms:
                             edf Ref.df       F p-value
s(t)                       1.998      2   76.58 < 2e-16
s(tsd.climatology)         1.989      2  141.55 < 2e-16
s(tmean.climatology)       1.995      2 5021.31 < 2e-16
s(log(t.range))            1.976      2  105.57 < 2e-16
s(log(b.range))            1.974      2   40.42 < 2e-16
s(log(z))                  1.969      2   25.18 5.8e-12
s(log(dZ))                 0.764      2    1.49   0.041
s(log(ddZ))                1.361      2   15.25 2.6e-09
s(log(substrate.grainsize)) 1.952      2  296.83 < 2e-16

R-sq.(adj) =  0.739   Deviance explained = 73.9%
GCV = 0.0088125  Scale est. = 0.0088014  n = 13444

# -------------------
#  varname ="pca2"
Family: gaussian
Link function: identity

Formula:
pca2 ~ s(t, k = 3, bs = "ts") + s(tsd.climatology, k = 3, bs = "ts") +
    s(tmean.climatology, k = 3, bs = "ts") + s(log(t.range),
    k = 3, bs = "ts") + s(log(b.range), k = 3, bs = "ts") + s(log(z),
    k = 3, bs = "ts") + s(log(dZ), k = 3, bs = "ts") + s(log(ddZ),
    k = 3, bs = "ts") + s(log(substrate.grainsize), k = 3, bs = "ts")

---

Parametric coefficients:
            Estimate Std. Error t value Pr(>|t|)
(Intercept)  0.02133    0.00081    26.3   <2e-16

Approximate significance of smooth terms:
                                edf Ref.df      F p-value
s(t)                       1.79e+00      2  101.1 < 2e-16
s(tsd.climatology)         1.98e+00      2  194.5 < 2e-16
s(tmean.climatology)       1.93e+00      2   45.2 < 2e-16
s(log(t.range))            1.93e+00      2   20.9 4.1e-10
s(log(b.range))            1.91e+00      2   29.3 4.9e-14
s(log(z))                  1.70e+00      2 2350.7 < 2e-16
s(log(dZ))                 4.76e-08      2    0.0 0.77011
s(log(ddZ))                1.93e+00      2    7.5 0.00041
s(log(substrate.grainsize)) 1.42e+00      2  413.4 < 2e-16

R-sq.(adj) =  0.699   Deviance explained = 69.9%
GCV = 0.0088397  Scale est. = 0.0088295  n = 13444



if (p$stmv_local_modelengine %in% c("gaussianprocess") )  p$libs = c( p$libs  )
if (p$stmv_local_modelengine %in% c("splancs") )  p$libs = c( p$libs, "splancs" )

if (p$stmv_local_modelengine %in% c("LaplacesDemon") )  p$libs = c( p$libs, "LaplacesDemonCpp" )
if (p$stmv_local_modelengine %in% c("stan") ) p$libs = c( p$libs, "cmdstanr" )
# if (p$stmv_local_modelengine %in% c("spate") )  p$libs = c( p$libs, "spate" ) # now copied directly into stmv

if (p$stmv_local_modelengine %in% c("stan", "gaussianprocess")) {
  message( "Precompiling stan program ... ")
  p$stanmodel = gaussian_process_stanmodel()  # precompile outside of the loop to speed the rest
}


habitat = stmv__habitat( p, dat, pa ), # TODO

gaussianprocess = stmv__gaussianprocess( p, dat, pa ),  # TODO
LaplacesDemon = stmv__LaplacesDemon( p, dat, pa ),
stan = stmv__stan( p, dat, pa, stanmodel=p$stanmodel ),  ## todo
splancs = stmv__splancs( p, dat, pa ), # TODO
spate = stmv__spate( p, dat, pa, sloc=Sloc[Si,], distance=stmv_distance_cur, nu=nu, phi=phi, varObs=varObs, varSpatial=varSpatial),
