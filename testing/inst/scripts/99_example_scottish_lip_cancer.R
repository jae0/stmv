
# There are many ways of running CAR /ICAR /IAR/ BYM
# The following are tests for speed and consistency to choose the best method
# The data are from CARBayes vignette the Scottish Lip cancer data
# These are simple spatial models (no time)

## summary: all give the same results. diseasemapping is the fastest and simplest
##.. all give essentially the same results (for Total population size)


# require(carstm)
loadfunctions("carstm")



# --------------------
# 0. inla  directly .. inla offers several variations of the bym, the most notable being bym and bym2
  reqiure(INLA)

  DD = scottish_lip_cancer_data()
  DD$ridfac = factor( row.names(DD), levels=row.names(DD) )
  DD$rid = as.numeric( DD$ridfac )
  DD = DD[ order(DD$rid),]

  # create neigbourhood structure
  NB_graph <- spdep::poly2nb(DD, row.names=DD$rid )  # snap=1000

  fit.bym2 = inla(
    formula = observed ~ pcaff + f( rid, model="bym2", graph=NB_graph, scale.model=TRUE, constr=TRUE),
    E=DD$expected,
    data = as.data.frame(DD),
    family="poisson",
    control.compute=list(dic=TRUE, config=TRUE, cpo=TRUE, waic=TRUE), # config=TRUE causes linear predictors to compute predictions/simulations quickly
    control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
    control.predictor=list(compute=TRUE, link =1 ), # compute=TRUE on each data location
    quantiles=NULL,  # no quants to save storage requirements
    # control.inla = list(h =1e-2), # h=0.01 is default step length for gradient calc of hyper params
    # control.fixed = list(expand.factor.strategy='inla') ,
    # working.directory= itmpdir, keep=FALSE,
    verbose=TRUE
  )

  summary(fit.bym2)
  # Time used:
  #  Pre-processing    Running inla Post-processing           Total
  #          1.0141          1.9314          0.3357          3.2812
  #
  # Fixed effects:
  #                mean     sd    mode kld
  # (Intercept) -0.2138 0.1299 -0.2115   0
  # pcaff        0.0361 0.0135  0.0364   0
  #
  # Random effects:
  # Name	  Model
  #  rid   BYM2 model
  #
  # Model hyperparameters:
  #                     mean     sd   mode
  # Precision for rid 4.5794 1.4197 3.9899
  # Phi for rid       0.7868 0.1662 0.9627
  #
  # Expected number of effective parameters(std dev): 30.15(3.453)
  # Number of equivalent replicates : 1.857
  #
  # Deviance Information Criterion (DIC) ...............: 297.66
  # Deviance Information Criterion (DIC, saturated) ....: 89.69
  # Effective number of parameters .....................: 30.38
  #
  # Watanabe-Akaike information criterion (WAIC) ...: 293.70
  # Effective number of parameters .................: 20.17
  #
  # Marginal log-Likelihood:  -135.20
  # CPO and PIT are computed
  #
  # Posterior marginals for linear predictor and fitted values computed



  fit.bym = inla(
    formula = observed ~ pcaff + f( rid, model="bym", graph=NB_graph, scale.model=TRUE, constr=TRUE),
    E=DD$expected,
    data = as.data.frame(DD),
    family="poisson",
    control.compute=list(dic=TRUE, config=TRUE, cpo=TRUE, waic=TRUE), # config=TRUE causes linear predictors to compute predictions/simulations quickly
    control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
    control.predictor=list(compute=TRUE, link =1 ), # compute=TRUE on each data location
    quantiles=NULL,  # no quants to save storage requirements
    # control.inla = list(h =1e-2), # h=0.01 is default step length for gradient calc of hyper params
    # control.fixed = list(expand.factor.strategy='inla') ,
    # working.directory= itmpdir, keep=FALSE,
    verbose=TRUE
  )

  summary(fit.bym)
  # slightly faster ..

  # Time used:
  #  Pre-processing    Running inla Post-processing           Total
  #          0.1404          0.5585          0.1218          0.8207
  #
  # Fixed effects:
  #                mean     sd    mode kld
  # (Intercept) -0.1928 0.1223 -0.1952   0
  # pcaff        0.0337 0.0129  0.0344   0
  #
  # Random effects:
  # Name	  Model
  #  rid   BYM model
  #
  # Model hyperparameters:
  #                                           mean       sd    mode
  # Precision for rid (iid component)     1836.043 1792.582 340.047
  # Precision for rid (spatial component)    4.572    1.541   3.889
  #
  # Expected number of effective parameters(std dev): 28.67(3.401)
  # Number of equivalent replicates : 1.953
  #
  # Deviance Information Criterion (DIC) ...............: 297.16
  # Deviance Information Criterion (DIC, saturated) ....: 89.18
  # Effective number of parameters .....................: 28.89
  #
  # Watanabe-Akaike information criterion (WAIC) ...: 294.26
  # Effective number of parameters .................: 19.92
  #
  # Marginal log-Likelihood:  -136.77
  # CPO and PIT are computed
  #
  # Posterior marginals for linear predictor and fitted values computed





# --------------
    # comparison of bym vs bym2: very similar
    # fixed effects are similar:

    # R> fit.bym2$summary.fixed
#                mean      sd     mode       kld
# (Intercept) -0.21375 0.12991 -0.21150 4.486e-06
# pcaff        0.03612 0.01355  0.03643 1.762e-05

    # R> fit.bym$summary.fixed
#               mean      sd     mode       kld
# (Intercept) -0.19278 0.12229 -0.19521 1.407e-05
# pcaff        0.03366 0.01288  0.03442 1.216e-05

# diseasemapping's fixed effects results: are also wuite similar
    # R> fit$inla$summary.fixed
    #                 mean      sd 0.025quant 0.5quant 0.975quant    mode       kld
    # (Intercept) -0.19745 0.12756  -0.447759 -0.19777    0.05426 -0.1984 6.932e-06
    # pcaff        0.03405 0.01343   0.006986  0.03427    0.05991  0.0347 9.435e-06


# -------------------
    # random effects also similar
    # R> cor( fit.bym2$summary.random$rid$mean , fit.bym$summary.random$rid$mean)
    # [1] 0.9383

    # but still significantly different
    # R> fit.bym2$summary.hyperpar
#                     mean     sd   mode
# Precision for rid 4.5794 1.4197 3.9899 == 0.4673 sd
# Phi for rid       0.7868 0.1662 0.9627 (space dominates)

    # R> fit.bym$summary.hyperpar
#                                           mean       sd    mode
# Precision for rid (iid component)     1836.043 1792.582 340.047  == 0.02334
# Precision for rid (spatial component)    4.572    1.541   3.889 == 0.4677  (spatial dominates)

  # diseasemapping's results: are intermediate between bym and bym2
  # Model hyperparameters:
    #                                                   mean       sd 0.025quant 0.5quant 0.975quant
    # Precision for region.indexS (iid component)     353.46 372.4247    21.6270  244.517   1380.499 == 0.053 sd
    # Precision for region.indexS (spatial component)   1.99   0.6735     0.9927    1.883      3.606 == 0.7089 sd  (space dominates)


    # stan solution to bym2:
    #                 mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat

    # intercept: similar
    # B[1]           -0.21    0.00 0.13  -0.47  -0.30  -0.21  -0.12   0.06  5306    1

    # pcaff: similar
    # B[2]            0.03    0.00 0.01   0.01   0.03   0.03   0.04   0.06  4861    1

    # sdSpatial: similar to bym2
    # sigma           0.53    0.00 0.08   0.38   0.47   0.52   0.58   0.71  3263    1

    # Phi for rid: == rho in inla .. similar
    # rho             0.88    0.00 0.14   0.48   0.82   0.93   0.98   1.00  1534    1
    # logit_rho       3.05    0.05 2.28  -0.10   1.50   2.58   4.14   8.81  1967    1

    # varepsilon[1]   0.42    0.01 0.98  -1.53  -0.24   0.43   1.09   2.32 12000    1
    # phi[1]          1.38    0.00 0.41   0.58   1.11   1.37   1.65   2.19 12000    1
    # lp__          751.47    0.17 8.84 733.32 745.63 751.74 757.64 768.24  2857    1
    #

    # bottom line: mostly consistent ... space dominates

  # choose
  # fit = fit.bym2
  # fit = fit.bym

  plot( fit$marginals.hyperpar$"Phi for strata", type="l")  # posterior distribution of phi
  plot( fit$marginals.hyperpar$"Precision for strata", type="l")
  plot( fit$marginals.hyperpar$"Precision for setno", type="l")



# --------------------
# 1. diseasemapping/inla  .. fastest .. ~ 1 min, laplace approximations via inla
  require(diseasemapping)

  DD = scottish_lip_cancer_data()
  DD$rid = row.names(DD)
  DD$log_offset = log(DD$expected)

  #  keep = which ( is.finite(DD$SA + DD$Total))
  # DD = DD[keep,]

  # create neigbourhood structure
  NB_graph <- spdep::poly2nb(DD, row.names = row.names(DD))  # snap=1000


  # CAR using a  full precision matrix form .. takes some time / MVN
  fit = bym(formula = observed ~ offset(log_offset) + pcaff, adjMat=NB_graph, data=as.data.frame(DD),
            family="poisson", region.id="rid", priorCI = list(sdSpatial = c(0.02, 100), sdIndep = c(0.02, 100)),
            control.compute = list(dic = TRUE, waic = TRUE) )

  #  fit$parameters$summary
  #                 mean      sd 0.025quant 0.5quant 0.975quant     mode       kld meanExp
  # (Intercept) -0.19731 0.12748   -0.44740 -0.19767    0.05432 -0.19841 2.540e-14  0.8248
  # pcaff        0.03399 0.01342    0.00696  0.03421    0.05983  0.03464 4.428e-13  1.0382
  # sdSpatial    0.73371      NA    0.52681  0.72882    1.00370  0.77006        NA      NA
  # sdIndep      0.07493      NA    0.02737  0.06416    0.21256  0.13130        NA      NA

  summary(fit$inla)

    # Random effects:
    # Name	  Model
    #  region.indexS   BYM model
    #
    # Model hyperparameters:
    #                                                   mean       sd 0.025quant 0.5quant 0.975quant
    # Precision for region.indexS (iid component)     353.46 372.4247    21.6270  244.517   1380.499
    # Precision for region.indexS (spatial component)   1.99   0.6735     0.9927    1.883      3.606
    #                                                   mode
    # Precision for region.indexS (iid component)     60.705
    # Precision for region.indexS (spatial component)  1.687
    #
    # Expected number of effective parameters(std dev): 30.12(3.45)
    # Number of equivalent replicates : 1.859
    #
    # Deviance Information Criterion (DIC) ...............: 297.34
    # Deviance Information Criterion (DIC, saturated) ....: 89.36
    # Effective number of parameters .....................: 30.28
    #
    # Watanabe-Akaike information criterion (WAIC) ...: 293.30
    # Effective number of parameters .................: 20.03
    #
    # Marginal log-Likelihood:  -154.65
    # Posterior marginals for linear predictor and fitted values computed

    # As disease mapping uses INLA, the comparison with bym and bym2 is direct .. they are similar but not the same:

    # R> fit$inla$summary.fixed
    #                 mean      sd 0.025quant 0.5quant 0.975quant    mode       kld
    # (Intercept) -0.19745 0.12756  -0.447759 -0.19777    0.05426 -0.1984 6.932e-06
    # pcaff        0.03405 0.01343   0.006986  0.03427    0.05991  0.0347 9.435e-06

    # R> fit$inla$summary.hyperpar
    #                                                   mean       sd 0.025quant 0.5quant 0.975quant
    # Precision for region.indexS (iid component)     353.46 372.4247    21.6270  244.517   1380.499
    # Precision for region.indexS (spatial component)   1.99   0.6735     0.9927    1.883      3.606
    #                                                   mode
    # Precision for region.indexS (iid component)     60.705
    # Precision for region.indexS (spatial component)  1.687


  DD$fitted.mean = fit$data$fitted.mean
  DD$log_random.sd = log( fit$data$random.sd )

  DD$fitted.exp = fit$data$fitted.exp
  DD$fitted.sd = fit$data$fitted.sd


  DD$log_random.sd = log( fit$data$random.sd )
  vn="log_random.sd"
  brks = interval_break(X=DD[[vn]], n=length(p$mypalette), style="quantile")
  spplot( DD, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent", xlim=p$ns.box$xlim, ylim=p$ns.box$ylim )

  plot( log_random.sd ~ log(SA), DD)
  plot( log(fitted.exp) ~ log(Total/SA), DD)
  hist( DD$fitted.exp - DD$Total/DD$SA , "fd" )
  hist( DD$log_random.sd , "fd" )
  hist( log(DD$fitted.sd) , "fd" )


  vn="random.sd"
  brks = interval_break(X=DD[[vn]], n=length(p$mypalette), style="quantile")
  spplot( DD, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent", xlim=p$ns.box$xlim, ylim=p$ns.box$ylim )



  # ----------------------
  # 2. reasonably fast .. ~ 5 min, full mcmc

  require(CARBayes)

  # CAR using a  full precision matrix form .. takes some time / MVN

  DD = scottish_lip_cancer_data()
  DD$rid = row.names(DD)
  DD$log_offset = log(DD$expected)

  NB_graph <- spdep::poly2nb(DD, row.names = row.names(DD))  # snap=1000

  # create neigbourhood structure
  NB_mat = nb2mat( NB_graph, style="B" )

  fit = S.CARbym(formula=observed ~ offset(log_offset) + pcaff, family="poisson", W=NB_mat, data=as.data.frame(DD),
                 burnin=50000, n.sample=200000, thin=10)

  #             Median    2.5%  97.5% n.sample % accept n.effective Geweke.diag
  # (Intercept) -0.2153 -0.4755 0.0555    15000     35.5      1343.7        -0.9
  # pcaff        0.0361  0.0062 0.0637    15000     35.5      1114.9         0.9
  # tau2         0.4822  0.2221 0.9721    15000    100.0      3214.9         0.5
  # sigma2       0.0078  0.0020 0.0745    15000    100.0       819.4         1.5
  #
  # DIC =  301       p.d =  31.17       Percent deviance explained =  59.45

  ## tau2 is spatial component in CARBayes .. sqrt(tau2 ) = sdSpatial  in diseasemapping = 0.6944
  ## sigma2 is nonspatial component in CARBayes .. sqrt(sigma2)=  sdIndep  in diseasemapping = 0.0883
  ## intercept is the same in both


  DD$fitted_mean_CARBayes = apply( fit$samples$fitted , 2, mean )
  DD$fitted_sd_CARBayes = apply( fit$samples$fitted , 2, sd )

  plot( DD$Total ~ DD$SA )
  plot( DD$fitted_mean_CARBayes ~ DD$SA )
  plot( DD$fitted_sd_CARBayes ~ DD$SA )




  # ----------------------
  # 3. brms/STAN .. slowest: ~ 30 min, and needs a lot of RAM : 34 GB
  #  .. but extremely flexible approach: full mcmc, can model covariates as smooths, have timeseries autocorrelation, etc.

  require (brms)

  # define correlation model:
  # a BYM is the same as an Intrinsic Autoregressive Model (IAR), which is called "esicar" by brms
  # see http://mc-stan.org/documentation/case-studies/AIR_Stan.html
  # CAR using a  full precision matrix form .. takes some time / MVN

  DD = scottish_lip_cancer_data()
  DD$rid = row.names(DD)
  DD$log_offset = log(DD$expected)

    #  keep = which ( is.finite(DD$SA + DD$Total))
  # DD = DD[keep,]

  # create neigbourhood structure
  NB_graph <- spdep::poly2nb(DD, row.names = row.names(DD))  # snap=1000

  # NB_graph = poly2nb( DD, row.names=DD[["rid"]], snap=1 )
  NB_mat = nb2mat( NB_graph, style="B" )
  W = nb2mat( NB_graph, style="B" )
  W = Matrix::Matrix(W, sparse = TRUE)
  if (!Matrix::isSymmetric(W, check.attributes = FALSE)) stop2("'W' must be symmetric.")
  not_binary = !(W == 0 | W == 1)
  if (any(not_binary)) {
      message("Converting all non-zero values in 'W' to 1")
      W[not_binary] = 1
  }
  brms_data = list( data=as.data.frame(DD), W=W )

  bym_model = cor_car(NB_mat, type="esicar")

  microbenchmark::microbenchmark( {

    fit = brm( observed ~ (1|rid) + offset(log_offset) + pcaff,
             data=as.data.frame(DD),
             family=poisson(link="log"),
             autocor=bym_model,
             chains=4, cores=4, algorithm="sampling", # "sampling" means do mcmc
             iter=6000, warmup=200, thin=10,  #  can increase these if required .. diagnostics seem ok
             control = list(max_treedepth=10, adapt_delta=0.85),
             save_model=file.path(p$work.directory, "bym.stan") )

  }, times=1 ) # 45 seconds


  # Gradient evaluation took 0.0027 seconds
  # 1000 transitions using 10 leapfrog steps per transition would take 27 seconds.

  summary(fit)

  # Correlation Structures:
  #       Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
  # sdcar     0.74      0.14     0.48     1.03       1548 1.00
  #
  # Group-Level Effects:
  # ~rid (Number of levels: 56)
  #               Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
  # sd(Intercept)     0.15      0.10     0.01     0.39       1766 1.00
  #
  # Population-Level Effects:
  #           Estimate Est.Error l-95% CI u-95% CI Eff.Sample Rhat
  # Intercept    -0.20      0.14    -0.46     0.06       2205 1.00
  # pcaff         0.03      0.01     0.01     0.06       2199 1.00

  # VarCorr(fit)  # for mor eprecise output

  # $rid
  # $rid$sd
  #             Estimate Est.Error     Q2.5  Q97.5
  # Intercept   0.1512    0.1042 0.007478 0.3863



  plot(fit, ask = FALSE)

  #sims:
  fit$fit@sim$samples

  # ## extract random effects standard devations and covariance matrices


  head(ranef(fit))
  ## predict responses based on the fitted model
  head(predict(fit))
  ## plot marginal effects of each predictor
  plot(marginal_effects(fit), ask = FALSE)
  ## investigate model fit

  loo(fit)  # lower is better

  WAIC(fit)
  pp_check(fit)

  sigmas = exp(posterior_samples(fit, "car"))
  ggplot(stack(sigmas), aes(values)) +
    geom_density(aes(fill = ind))




  # --------------------------------------------
  # 4. CAR basic form: using a  full precision matrix form and multivariate normal
  # very slow .. only for demo of code structure
  
  library(cmdstanr)

  # CAR using a  full precision matrix form .. takes some time / MVN

  DD = scottish_lip_cancer_data()
  DD$rid = row.names(DD)

  NB_graph <- spdep::poly2nb(DD, row.names = row.names(DD))  # snap=1000

  W = spdep::nb2mat(NB_graph, style="B") # adjacency matrix
  X = model.matrix(~1+pcaff, DD) #  no covariates, intercept only
  # X = model.matrix(~scale(pcaff)-1, DD)  # made up example covariate

  DS =  list(
    K = nrow(X),         # number of observations
    L = ncol(X),         # number of coefficients
    X = X,               # design matrix
    Y = DD$observed,            # observed number of cases
    log_offset = log(DD$expected), # offset.. expected num. cases
    W = W   # adjacency matrix
  )

  MS = stan_initialize( stan_code=stan_model_code( "car_mvn" ) )
  MS$compile()

  nchains = 4
  niter = 6000 # 63 seconds

  microbenchmark::microbenchmark(  {
    fit  = MS$sample( 
      data=DS, 
      iter_warmup = niter,
      iter_sampling = 1000,
      seed = 123,
      chains = 4,
      parallel_chains = 4,  # The maximum number of MCMC chains to run in parallel.
      max_treedepth = 18,
      adapt_delta = 0.975,
      refresh = 500
    }, times=1 
  )


  pars = c('B[1]', 'B[2]', 'alpha', 'tau_sd', 'lp__')
  print(fit, pars=pars)
  traceplot(fit, pars = pars)   # visualize results

  stan_dens(fit)
 
  mcmc = stan_extract( as_draws_df(fit$draws() ) )

  # intercept is off ..

  #          mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
  # B[1]     1.17    0.01 0.40   0.18   1.00   1.23   1.41   1.84   738 1.01
  # B[2]     0.02    0.00 0.01  -0.01   0.01   0.02   0.03   0.04  2247 1.00
  # alpha    0.94    0.00 0.07   0.76   0.92   0.96   0.98   1.00  1902 1.00
  # tau_sd   0.84    0.00 0.14   0.60   0.74   0.83   0.93   1.16  2304 1.00
  # lp__   814.90    0.14 6.91 800.98 810.38 815.16 819.85 827.47  2307 1.00



  # --------------------------------------------
  # 5. CAR sparse form

   library(cmdstanr)


  # CAR using a  full precision matrix form .. takes some time / MVN

  DD = scottish_lip_cancer_data()
  DD$rid = row.names(DD)
  DD$log_offset = log(DD$expected)

  NB_graph <- spdep::poly2nb(DD, row.names = row.names(DD))  # snap=1000

  W = spdep::nb2mat(NB_graph, style="B") # adjacency matrix
  X = model.matrix(~1+pcaff, DD) #  no covariates, intercept only
  # X = model.matrix(~scale(pcaff)-1, DD)  # made up example covariate

  DS =  list(
    K = nrow(X),         # number of observations
    L = ncol(X),         # number of coefficients
    X = X,               # design matrix
    Y = DD$observed,            # observed number of cases
    log_offset = log(DD$expected), # offset.. expected num. cases
    W = W,   # adjacency matrix
    Wn = sum(W) / 2
  )

  MS = stan_initialize( stan_code=stan_model_code( "car_sparse" ) )
  MS$compile()

  microbenchmark::microbenchmark(  {
    fit  = MS$sample( 
      data=DS, 
      iter_warmup = niter,
      iter_sampling = 1000,
      seed = 123,
      chains = 4,
      parallel_chains = 4,  # The maximum number of MCMC chains to run in parallel.
      max_treedepth = 18,
      adapt_delta = 0.975,
      refresh = 500
    }, times=1 
  )

  pars = c('B[1]', 'B[2]', 'alpha', 'tau_sd', 'lp__')
  print(fit, pars=pars)
  traceplot(fit, pars = pars)   # visualize results

  stan_dens(fit)
 
  mcmc = stan_extract( as_draws_df(fit$draws() ) )

 
  #          mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
  # B[1]    -0.26    0.02 0.38  -1.00  -0.48  -0.27  -0.05   0.51   598    1
  # B[2]     0.03    0.00 0.01   0.01   0.02   0.03   0.04   0.06  4138    1
  # alpha    0.95    0.00 0.05   0.82   0.94   0.97   0.99   1.00  3805    1
  # tau_sd   0.85    0.00 0.13   0.64   0.76   0.84   0.93   1.13  4109    1
  # lp__   779.84    0.12 6.67 765.90 775.59 780.14 784.48 791.96  3153    1


  # --------------------------------------------
  # 6. sparse BYM using Carlin's winbugs code
  library(cmdstanr)

  # CAR using a  full precision matrix form .. takes some time / MVN

  DD = scottish_lip_cancer_data()
  DD$rid = row.names(DD)

  NB_graph <- spdep::poly2nb(DD, row.names = row.names(DD))  # snap=1000

  # X = model.matrix(~1, DD) #  no covariates, intercept only
  # X = model.matrix(~scale(pcaff)-1, DD)  # made up example covariate
  X = model.matrix(~pcaff + 1, DD)  # made up example covariate


  DS  = list(
    K = nrow(X),
    L = ncol(X),
    X = X,
    Y = DD$observed,
    log_offset = log( DD$expected),
    W = W   # adjacency matrix
  )

  DS = c( DS, nb2edge(NB_graph) )

  MS = stan_initialize( stan_code=stan_model_code( "bym_basic" ) )
  MS$compile()

  nchains = 4
  niter = 6000 # 24 seconds

  microbenchmark::microbenchmark(  {
    fit  = MS$sample( 
      data=DS, 
      iter_warmup = niter,
      iter_sampling = 1000,
      seed = 123,
      chains = 4,
      parallel_chains = 4,  # The maximum number of MCMC chains to run in parallel.
      max_treedepth = 18,
      adapt_delta = 0.975,
      refresh = 500
    }, times=1 
  )

  pars = c('B', 'varepsilon_sd', 'phi_sd',   'varepsilon[1]', 'phi[1]', 'lp__') # phi_sd==spatial==carsd; varep_sd == nospatial sd
  print(fit, pars =pars)
  traceplot(fit, pars = pars)   # visualize results
  stan_dens(fit)
  print(fit, pars=pars), probs=c(0.025, 0.5, 0.975), digits=3);


  mcmc = stan_extract( as_draws_df(fit$draws() ) )

 
  #       mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
  # B[1]           -0.29    0.00 0.17  -0.63  -0.41  -0.29  -0.17   0.06  3479    1
  # B[2]            0.04    0.00 0.02   0.01   0.03   0.04   0.05   0.08  3277    1
  # varepsilon_sd   0.48    0.00 0.07   0.37   0.43   0.48   0.52   0.63  6666    1
  # phi_sd          0.68    0.00 0.13   0.46   0.58   0.66   0.76   0.97  3351    1
  # varepsilon[1]   0.82    0.01 0.75  -0.66   0.32   0.83   1.32   2.27 12000    1
  # phi[1]          1.21    0.01 0.46   0.32   0.90   1.21   1.52   2.13  5791    1
  # lp__          756.45    0.15 8.67 738.30 750.77 756.75 762.60 772.26  3511    1
  #



  # --------------------------------------------
  # 7. I(C)AR sparse form using "dot_self" and **variance scaling factor**  .. permits more reasonable prior choices

  # as in INLA (due to D. Simpson .. find reference)
  # https://github.com/stan-dev/example-models/blob/master/knitr/car-iar-poisson/icar_stan.Rmd
  # NOTE:  mungeCARdata4stan.R (abstraction is found below) does the node creation. Originally from:
  # NOTE:  internally uses winbugs data formats and scaling factor

  library(cmdstanr)

  library(INLA)

  # CAR using a  full precision matrix form .. takes some time / MVN

  DD = scottish_lip_cancer_data()
  DD$rid = row.names(DD)

  NB_graph <- spdep::poly2nb(DD, row.names = row.names(DD))  # snap=1000

  X = model.matrix(~1, DD) #  no covariates, intercept only
  # X = model.matrix(~scale(rid)-1, DD)  # made up example covariate

  DS  = list(
    K = nrow(X),
    L = ncol(X),
    X = X,
    Y = DD$observed,
    log_offset = log(DD$expected),
    scaling_factor = bym_scaling_factor(NB_graph) # additional parameters required by the bym_scaled model (topology and variance scale factor)
  )
  DS = c( DS, nb2edge(NB_graph) )
 
  MS = stan_initialize( stan_code=stan_model_code( "bym_scaled" ) )
  MS$compile()

  nchains = 4

  niter = 6000 # 22

  microbenchmark::microbenchmark(  {

    fit  = MS$sample( 
      data=DS, 
      iter_warmup = niter,
      iter_sampling = 1000,
      seed = 123,
      chains = 4,
      parallel_chains = 4,  # The maximum number of MCMC chains to run in parallel.
      max_treedepth = 18,
      adapt_delta = 0.975,
      refresh = 500
    }, times=1 
  )

  pars = c('B', 'sigma', 'rho', 'logit_rho', 'varepsilon_sd', 'phi_sd', 'varepsilon[1]', 'phi[1]', 'lp__')
  print(fit, pars=pars )
  traceplot(fit, pars=pars)   # visualize results
  stan_dens(fit)
  print(fit, pars=pars, probs=c(0.025, 0.5, 0.975), digits=3);
  capture.output(print(fit, digits=3, probs=c(0.025, 0.975)), file=ofile);


  mcmc = stan_extract( as_draws_df(fit$draws() ) )

 
#                 mean se_mean   sd   2.5%    25%    50%    75%  97.5% n_eff Rhat
# B[1]           -0.21    0.00 0.13  -0.47  -0.30  -0.21  -0.12   0.06  5306    1
# B[2]            0.03    0.00 0.01   0.01   0.03   0.03   0.04   0.06  4861    1
# sigma           0.53    0.00 0.08   0.38   0.47   0.52   0.58   0.71  3263    1
# rho             0.88    0.00 0.14   0.48   0.82   0.93   0.98   1.00  1534    1
# logit_rho       3.05    0.05 2.28  -0.10   1.50   2.58   4.14   8.81  1967    1
# varepsilon[1]   0.42    0.01 0.98  -1.53  -0.24   0.43   1.09   2.32 12000    1
# phi[1]          1.38    0.00 0.41   0.58   1.11   1.37   1.65   2.19 12000    1
# lp__          751.47    0.17 8.84 733.32 745.63 751.74 757.64 768.24  2857    1
#

get_posterior_mean(fit, pars)




  # cov neighbors
  cov(mcmc$phi[,6],mcmc$phi[,8]);
  cov(mcmc$phi[,6],mcmc$phi[,3]);
  cov(mcmc$phi[,10],mcmc$phi[,22]);
  # cov non-neighbors
  cov(mcmc$phi[,6],mcmc$phi[,54]);
  cov(mcmc$phi[,8],mcmc$phi[,54]);
  cov(mcmc$phi[,2],mcmc$phi[,55]);
  cov(mcmc$phi[,1],mcmc$phi[,55]);
  cov(mcmc$phi[,2],mcmc$phi[,50]);


  ## TODO .. spatial plots of residuals, etc.
