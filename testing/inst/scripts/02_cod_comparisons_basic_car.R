
# ---------------------------
# Atlantic cod comparison of naive strata-based averages vs CAR/BYM Poisson process models
#

RLibrary( "sp", "spdep" )  # standard CRAN libs -- geostatistical support
RLibrary( "aegis", "bio.taxonomy" )# , "carstm") #,  ) # locally developed code

loadfunctions("carstm")  # require(carstm)   # for fast debugging

p = list(
  id ="Atlantic cod summer standardtow",
  speciesname = "Atlantic_cod",
  groundfish_species_code = 10,   #  10= cod
  yrs = 2010:2017,
  season = "summer",
  areal_units_type = "stratanal_polygons_pre2014",  # "pre2014" for older
  trawlable_units = "towdistance"  # <<<<<<<<<<<<<<<<<<
  # trawlable_units = "standardtow"
  # trawlable_units = "sweptarea"
)


# ------------------------------------------------
# construct basic parameter list
p = aegis.survey::survey_parameters(
  p=p,
  selection=list(
    id = gsub("[[:space:]]", "_", paste(p$speciesname, p$season, p$trawlable_units, sep=".") ),
    trawlable_units = p$trawlable_units,
    biologicals=list(
      spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code )
    ),
    survey=list(
      data.source="groundfish",
      yr = p$yrs,      # time frame for comparison
      months=6:8,
      # dyear = c(150,250)/365, #  summer = which( (x>150) & (x<250) ) , spring = which(  x<149 ), winter = which(  x>251 )
      settype = 1,
      gear = c("Western IIA trawl", "Yankee #36 otter trawl"),
      # strata_toremove=c("Gulf", "Georges_Bank", "Spring", "Deep_Water"),  # strata to remove from standard strata-based analysis
      polygon_enforce=TRUE,
      ranged_data="dyear"
    )
  )
)
# basic spatial parameters
# set up default map projection
p$areal_units_proj4string_planar_km = projection_proj4string("omerc_nova_scotia")   # oblique mercator, centred on Scotian Shelf rotated by 325 degrees




## --------------------------------
# this is to make sure/emphasize that all data enters analysis initially ..
# and then sums only after spatial modelling
# you can of course remove some strata if you have a valid concern that they are anomalous
# based on apriori considerations but then this should be done by specifying sppoly more carefully and not here
p$selection$survey$strata_toremove = NULL
set = survey_db( p=p, DS="filter" , add_groundfish_strata=TRUE  )

# ensure we have some estimate of sweptarea and choose the appropriate one based upon which trawlable units we are using
ft2m = 0.3048
m2km = 1/1000
nmi2mi = 1.1507794
mi2ft = 5280
standardtow_sakm2 = (41 * ft2m * m2km ) * ( 1.75 * nmi2mi * mi2ft * ft2m * m2km )  # surface area sampled by a standard tow in km^2  1.75 nm

set$data_offset = switch( p$selection$trawlable_units,
  standardtow =  rep(standardtow_sakm2, nrow(set)) , # "standard tow"
  towdistance = set$sa_towdistance,  # "sa"=computed from tow distance and standard width, 0.011801==),
  sweptarea = set$sa  # swept area based upon stand tow width and variable lenths based upon start-end locations wherever possible

)

set$data_offset[which(!is.finite(set$data_offset))] = median(set$data_offset, na.rm=TRUE )
nd = which(!is.finite(set$totno))
if (length(nd) > 0) set$totno[nd] = NA
set$yr_factor = factor(set$yr)




# ------------------------------------------------
## using the "standard" polygon definitions  .. see https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
# Here we compute surface area of each polygon via projection to utm or some other appropriate projection.
# This adds some variabilty relative to "statanal"
# NOTE polygon areas of strata are predetermined in GSSTRATUM (warning in sq nautical miles)
# ... this can be variable due to projection being used and so requires reprojection to the area
# of interest to compute surface area of polygons
sppoly = maritimes_groundfish_strata( areal_units_timeperiod="pre2014" )
# prep fields required to help extract results from model fits and compute estimates of biomass given mean size and mean numbers
sppoly$au_sa_km2 = st_area( st_transform(sppoly, st_crs(p$areal_units_proj4string_planar_km) ) )  # km^2  .. planar_crs_km should be in km units
sppoly$strata_to_keep = ifelse( as.character(sppoly$AUID) %in% strata_definitions( c("Gulf", "Georges_Bank", "Spring", "Deep_Water") ), FALSE,  TRUE )
sppoly$AUID = factor(sppoly$AUID, levels=levels(sppoly$AUID) ) # make sure in case not all strata are represented in set
sppoly$strata = as.numeric( sppoly$AUID )
sppoly = sppoly[order(sppoly$strata),]
sppoly$roworder = 1:length(sppoly)
if (0) plot(sppoly)


# ------------------------------------------------
# neighbourhood structure
NB_graph = maritimes_groundfish_strata( areal_units_timeperiod="pre2014", returntype="neighbourhoods" )  # customized neighbourshood structure
scaling_factor = bym_scaling_factor(NB_graph)




## --------------------------------
# construct meanweights matrix
mw = meanweights_by_arealunit( set=set, AUID=as.character( sppoly$AUID ), yrs=p$yrs, fillall=TRUE )
# mw = meanweights_by_arealunit_modelled( p=p, redo=TRUE )  ## needed in carstm_output_compute


mwyr = data.frame(cbind(mw$AUID, meanweights=mw$meanweights), stringsAsFactors=FALSE)
names(mwyr) = c("AUID", p$yrs)
for (i in 2:ncol(mwyr)) mwyr[,i] = as.numeric(mwyr[,i])
mwyr = merge(mwyr, sppoly[, c("AUID", "roworder")] , by="AUID",  )
mwyr = mwyr[order(mwyr$roworder),]



# Model comparisons using "standard tow" assumptions
RES = data.frame(yr=p$selection$survey[["yr"]])


if (0) {
  # alternative access via strata_dataformat .. ie. directly from groundfish_survey_db
  p$selection$biologicals$spec = p$groundfish_species_code
  p$selection$biologicals$spec_bio = NULL
  set = strata_dataformat( selection )   # return values in kg or no per set

  # to access via det .. lots of unmeasured organisms results in total catches not being the same as in catch tables ... :(
  set = aegis.survey::survey_db( p=p, DS="filter", add_groundfish_strata=TRUE )   # return values in kg or no per set
  # equivalent to:
  set = aegis.survey::survey_db( p=p, DS="cat.filter", add_groundfish_strata=TRUE )   # return values in kg or no per set
  dim(set) # [1]
  sum(set$totwgt) # [1]
  sum(set$totno) #  [1]

  # this is slightly different due to use of individual data ... sometimes incomplete sub-sampling
  set = aegis.survey::survey_db( p=p, DS="det.filter", add_groundfish_strata=TRUE )   # return values in kg or no per set
  dim(set) # [1] 1684   48
  sum(set$totwgt) # [1] 4640
  sum(set$totno) #  [1] 5882
}




  RLibrary("cmdstanr")
  # BYM using STAN .. does not impute yet ... abandonned for now until it is a bit faster ... :(
  # debugging parameter values
  p$selection$survey$yr = 2010  # just choose a single year for debugging
  # result are stored here defined by p$carstm_data_root .. change if required

  p$carstm_data_root = project.datadirectory( "carstm", p$selection$id )

  loadfunctions("stmvdev")


  stancode = stan_initialize( model_code=stan_model_code( "bym_scaled_hierarchical" ) )
  stancode$compile()

  compiled_model= stan_model( model_code=stan_model_code( "bym_scaled_hierarchical" ) )

  o = bym_fit.stan(
    set=set,
    NB_graph=nb,
    compiled_model=stancode,
    carstm_data_root = p$carstm_data_root,
    # yrs = p$selection$survey$yr,
    yrs = 2010,
    # all after these must be stan-options ..
    # seed = 1111,
    iter = 3000,
    # thin = 10,
    chains = 4, #parallel::detectCores() - 1,
    cores = 4, # parallel::detectCores() - 1,
    warmup = 1000,
    # algorithm = "HMC",
    # algorithm = "NUTS"  #  "HMC, "HMC", etc,
    # pars = NA, # NA == everything
    # control =  list(adapt_delta=0.975, max_treedepth=16),
    # control =  list(adapt_delta=0.99, max_treedepth=15),
    # control =  list(adapt_delta=0.5),
    verbose =  TRUE
  )


  # debugging parameter values
    fit = bym_fit.stan( yrs=2010, loaddata=TRUE, carstm_data_root = p$carstm_data_root)
    pars = c('B[1]',  'rho', 'sigma', 'phi[1]', 'varepsilon[1]', 'mu[1]', 'lp__')
    print(fit, pars =pars)
                 mean se_mean    sd    2.5%     25%     50%     75%   97.5% n_eff   Rhat
B[1]             3.69    0.04  0.06    3.61    3.67    3.70    3.72    3.77     2  28.22
rho              0.51    0.16  0.23    0.28    0.38    0.43    0.56    0.90     2 124.64
sigma            1.74    0.29  0.41    1.29    1.47    1.64    1.91    2.40     2 109.89
phi[1]           0.03    0.13  0.19   -0.28   -0.06    0.08    0.18    0.23     2  26.18
varepsilon[1]   -0.48    0.20  0.29   -0.72   -0.67   -0.62   -0.40    0.04     2  31.12
mu[1]            0.13    0.23  0.34   -0.27   -0.09    0.00    0.23    0.89     2   6.79
lp__          7014.45   55.93 79.22 6893.37 6955.50 7027.99 7089.58 7102.79     2  22.73


    traceplot(fit, pars = pars)   # visualize results

    stan_ac(fit, pars )  # about 10 iteration autocorrelation

    # extract results from model fits and compute estimates of biomass given mean size and mean numbers

    o = bym_extract.stan( yrs=2010, carstm_data_root=p$carstm_data_root,
      strata=sppoly, meanweights=mw$meanweights,
      meanvars=c( "mu", "varepsilon", "phi"), sdvars=c( "varepsilon", "phi") )

    s = summary(fit)$summary
    plot( s[,"Rhat"],  s[,"n_eff"])
    # s[,"n_eff"] # effective sample size
    # s[,"Rhat"]  # scale reduction when swapping halves of different chains .. index of convergence..

    M = stan_extract( as_draws_df( fit$draws() ) )
  
  # other diagnostics
      stan_diag(fit, "treedepth")
      stan_diag(fit, "divergence")
      stan_diag(fit, "stepsize")
      stan_diag(fit, "sample")

      stan_rhat(fit, pars)
      stan_par(fit, pars[1])
      stan_ess(fit, pars=pars[2:3])  # effective sample size
      stan_mcse(fit, pars)

      hist( M[["sigma"]][[1]], "fd" )
      hist( M[["rho"]][[1]], "fd" )
      hist( M[["B[1]"]][[1]], "fd" )

      # see also: stan_plot, stan_scat, stan_hist, stan_dens, stan_ac

      stan_plot(fit, pars[1:4] )
      stan_hist(fit, pars[1:4] )
      stan_dens(fit, pars[1:4] )
      stan_scat(fit, pars[2:3] )

  # # to do loo ... must also output yhats in 'generated quantities' part of stan code
  #   yhat = exp(apply( M[["yhat"]], 2, mean ))
  #   plot( yhat ~ DATAINPUT$Y  )
  #   cor(yhat , DATAINPUT$Y  ) #
  #
  #   log_lik1 <- extract_log_lik(fit, merge_chains = FALSE)
  #   rel_n_eff <- relative_eff(exp(log_lik1))
  #   fit_loo_1 = loo(log_lik1, r_eff = rel_n_eff, cores = 2)
  #
  #   fit_loo_1

    o = bym_extract.stan( yrs=p$selection$survey$yr, carstm_data_root=p$carstm_data_root,
      strata=sppoly, meanweights=mw$meanweights,
      meanvars=c( "mu", "varepsilon", "phi"), sdvars=c( "varepsilon", "phi") )

    plot( o$aggregate_timeseries[,c("yr","totwgt_mean")] )

    # o$aggregate_sims  (order each by rank and plot sequentially) ---

    # a few maps
    sppoly$mu_exp = exp(o$strata_stat[["2010"]]$mu)
    vn = "mu_exp"
    brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
    spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent", xlim=p$boundingbox$xlim, ylim=p$boundingbox$ylim )


    vn = "mu"
    sppoly$mu = o$strata_stat[["2010"]]$mu
    brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
    spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent", xlim=p$boundingbox$xlim, ylim=p$boundingbox$ylim )


}





#  INLA

  # load one solution and examine contents
  fit = bym_fit.inla( p=p, yrs=2010, loaddata=TRUE, carstm_data_root=p$carstm_data_root)
  summary(fit)
  plot(fit, plot.fixed.effects=TRUE)

  # (fit$mode$mode.status > 0)   # make sure Eigenvalues of Hessian are appropriate (>0)
  ps = inla.posterior.sample(n=10, fit)
  postsims = sapply( ps, posterior_extract_from_inla, scaling_factor=scaling_factor)
  postsims_wgt = exp(postsims) * sppoly$au_sa_km2 * mwyr[,as.character(2010)]
  total_wgt = apply(postsims_wgt[sppoly$strata_to_keep,], 2, sum, na.rm=TRUE)
  c( mean(total_wgt), sd(total_wgt) )
# 43411853 23273763

  # Posterior marginals for the fixed effects and hyperparameters
  fit$summary.fixed    # Posterior marginals for the fixed effects
  fit$summary.hyperpar # Posterior marginals for the hyperparameters

  vn = "mu"
  means = matrix(fit$summary.random$AUID$mean, ncol=2, byrow=FALSE)
  colnames(means) = c("iid", "AUID")
  sppoly$mu = means[, "AUID"]
  brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
  spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent", xlim=p$boundingbox$xlim, ylim=p$boundingbox$ylim )

  plot( fit$marginals.hyperpar$"Phi for strata", type="l")  # posterior distribution of phi
  plot( fit$marginals.hyperpar$"Precision for strata", type="l")
  plot( fit$marginals.hyperpar$"Precision for setno", type="l")


  # effective sample size
  # scale reduction when swapping halves of different chains .. index of convergence..
  if ( any(s[,"n_eff"] > 100) | any( abs(1-s[,"Rhat"]) > 0.01) )  warning( paste( "Check this year for convergence issues: ", yr) )
  s = NULL



# get improved estimates for the hyperparameters
# fit = inla.hyperpar(fit , dz = 0.2 , diff.logdens =20)

# effective sample size
# scale reduction when swapping halves of different chains .. index of convergence..
if ( any(s[,"n_eff"] > 100) | any( abs(1-s[,"Rhat"]) > 0.01) )  warning( paste( "Check this year for convergence issues: ", yr) )
s = NULL


dta$strata_iid = dta$strata

  fit = inla (bym_formula,   family="normal",
    data=dta,
    E=dta$data_offset,
    control.compute=list(dic=TRUE, config=TRUE, cpo=TRUE, waic=TRUE), # config=TRUE causes linear predictors to compute predictions/simulations quickly
    control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
    control.predictor=list(compute=TRUE, link =1 ), # compute=TRUE on each data location
    verbose=TRUE
)


posterior_extract_from_inla = function( o, scaling_factor ) {
  # a function that is required to show inla what to compute from the posterior summary
  # .. in this case the mean posson intensities of each stratum
  # o = ps[[1]]
  hyper_names  = names(o$hyperpar )
  latent_names = rownames(o$latent )
  # latent_names are the rownames that will contain info about the indices ..
  intercept=grep("Intercept", latent_names, fixed=TRUE )
  strata_i = grep("AUID", latent_names, fixed=TRUE )
  strata=grep("strata", latent_names, fixed=TRUE )
  strata_phi=grep("Phi for strata", hyper_names, fixed=TRUE )
  strata_prec=grep("Precision for strata", hyper_names, fixed=TRUE )
  # strata_iid_prec=grep("Precision for strata_iid", hyper_names, fixed=TRUE )
  # setno_prec = 0
  out = matrix( o$latent[strata], ncol=2, byrow=FALSE )
  colnames(out) = c("varepsilon", "phi")
  strata_mu = o$latent[strata_i]
  rho = o$hyperpar[[strata_phi]]
  sigma = 1/sqrt(o$hyperpar[[strata_prec]])
  mu = o$latent[intercept] + sigma * (sqrt(1-rho) * out[,"varepsilon"] + sqrt(rho) * out[,"phi"] / sqrt(scaling_factor)  )
  return( mu )
}


p$carstm_data_root = project.datadirectory( "carstm", p$selection$id )


bym_fit.inla(
    set=set,
    AUID=sppoly$AUID,
    NB_graph=nb,
    carstm_data_root = p$carstm_data_root,
    yrs = p$selection$survey$yr,
    posterior_extract_from_inla=posterior_extract_from_inla,
    formula = bym_formula,
    family="poisson",
    data=dta,
    E=dta$data_offset,
    control.compute=list(dic=TRUE, config=TRUE, cpo=TRUE, waic=TRUE), # config=TRUE causes linear predictors to compute predictions/simulations quickly
    control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
    control.predictor=list(compute=TRUE, link =1 ), # compute=TRUE on each data location
    quantiles=NULL,  # no quants to save storage requirements
    # control.inla = list(h =1e-2), # h=0.01 is default step length for gradient calc of hyper params
    # control.fixed = list(expand.factor.strategy='inla') ,
    # working.directory= itmpdir, keep=FALSE,
    verbose=TRUE
)




# ------------------------------------------------
# Model 2: simple factorial operating directly upon totwgt
# this replicates the stratanal estimates but in an explicit model-based framework
# problems:
# - data distribution is not normal .. perhaps lognormal, poisson or overdispersed something

fit.simple.lm = lm(
  formula = totwgt ~ AUID:yr_factor - 1  ,
  data=set
)
u = coef(fit.simple.lm)
w = matrix( unlist(strsplit(names(u), ":")), ncol=2, byrow=TRUE)
w = gsub("AUID", "", w)
w = gsub("yr_factor", "", w)
sid = match( w[,1], as.character(sppoly$AUID) )
yid = match( w[,2], as.character(levels(set$yr_factor)))
out = matrix( NA, ncol=max(yid), nrow=max(sid) )
out[ cbind(sid, yid)] = u
out = out * sppoly$au_sa_km2 / standardtow_sakm2
RES$lm_gaussian_totwgt = NA
RES$lm_gaussian_totwgt = colSums( out[sppoly$strata_to_keep,] )
# 42933873 34568613  9113429 12298849 27173606  8732553 21470826 14873079  # some minor differences
lines( lm_gaussian_totwgt ~ yr, data=RES, lty=2, lwd=2, col="green")
AIC(fit.simple.lm)  # 19178


vn = "lm_gaussian_totwgt"
sppoly[,vn] = out$"2010"
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent", xlim=p$boundingbox$xlim, ylim=p$boundingbox$ylim )



ps = inla.posterior.sample(n=10000, fit)

  posterior_extract_from_inla = function( o, intercept, strata, strata_phi, strata_prec, scaling_factor ) {
    out = matrix( o$latent[strata], ncol=2, byrow=FALSE )
    colnames(out) = c("varepsilon", "phi")
    rho = o$hyperpar[strata_phi]
    sigma = 1/sqrt(o$hyperpar[strata_prec])
    mu = o$latent[intercept] + sigma * (sqrt(1-rho) * out[,"varepsilon"] + sqrt(rho/ scaling_factor) * out[,"phi"]  )
    return( mu )
  }

  hyper_names  = names(ps[[1]]$hyperpar )
  latent_names = rownames(ps[[1]]$latent )
  # latent_names are the rownames that will contain info about the indices ..
  # optimally the grep search should only be done once but doing so would
  # make it difficult to implement in a simple structure/manner ...
  # the overhead is minimal relative to the speed of modelling and posterior sampling

  postsims = sapply( ps, posterior_extract_from_inla,
    intercept=grep("Intercept", latent_names, fixed=TRUE ),
    strata=grep("strata", latent_names, fixed=TRUE ),
    strata_phi=grep("Phi for strata", hyper_names, fixed=TRUE ),
    strata_prec=grep("Precision for strata", hyper_names, fixed=TRUE ),
    scaling_factor=scaling_factor
  )

  postsims_wgt = exp(postsims) * sppoly$au_sa_km2 * mwyr$meanweight

  total_wgt = apply(postsims_wgt[sppoly$strata_to_keep,], 2, sum, na.rm=TRUE)
  c( mean(total_wgt), sd(total_wgt) )
# 6742048 2690330

  plot( fit$marginals.hyperpar$"Phi for strata", type="l")  # posterior distribution of phi
  plot( fit$marginals.hyperpar$"Precision for strata", type="l")
  plot( fit$marginals.hyperpar$"Precision for setno", type="l")



  # get improved estimates for the hyperparameters
  # fit = inla.hyperpar(fit , dz = 0.2 , diff.logdens =20)

  # effective sample size
  # scale reduction when swapping halves of different chains .. index of convergence..
  if ( any(s[,"n_eff"] > 100) | any( abs(1-s[,"Rhat"]) > 0.01) )  warning( paste( "Check this year for convergence issues: ", yr) )
  s = NULL
