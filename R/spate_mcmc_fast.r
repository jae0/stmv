
spate_mcmc_fast = function (y, yvar, n, SV=NULL, RWCov=NULL, pars=NULL,   
  dt= 1, NPosteriors=2000, BurnIn=1000, NcovUpdates=20, seed=NULL,  Drift=TRUE, Diffusion=TRUE, nu=1  )  {

 # simplified and modified from spate::spate.mcmc 
 # to do both prediction and parameter estimation and speed up as much as possible
 # removed incidence matrix method and covariate predictors

  if (0) {
    y=w; yvar=sdTotal^2; n=nsq; SV=NULL; Drift=TRUE; Diffusion=TRUE; 
    NPosteriors=2000; BurnIn=2000; NcovUpdates=20;
    RWCov=NULL; dt=1; MultCov=0.5; nu=1; seed=1 
  } 

  if (!is.null(seed)) set.seed(seed)
  
  yvar_max = yvar*10

  pnames = c("rho0", "sigma2", "zeta", "rho1", "gamma", "alpha", "muX", "muY", "tau2")
  logVars = c("rho0", "sigma2", "zeta", "rho1", "gamma", "tau2") 
  diffVars = c("rho1", "gamma", "alpha", "zeta" )
  driftVars = c("muX", "muY")

  todropVars = NULL
  if (!Diffusion) todropVars = c( todropVars, diffVars ) 
  if (!Drift) todropVars = c( todropVars, driftVars ) 

  if (is.null(SV)) {
    SV = c(rho0=0.1, sigma2=yvar/4, zeta=0.1, rho1=0.1, gamma=0.5, alpha=0.5, muX=0, muY=0, tau2=yvar/4)
  }

  if (!is.null(todropVars)) {
    pnames = pnames[-todropVars ]
    if (!is.null(SV)) {
      if (length(SV) != length(pnames) ) {
        # Assume we drop the same sequence:
        SV = SV[ -todropVars]
      }
    }
  } 
  nParms = length(pnames)
  indEst = 1:nParms
  logInd = which( pnames %in% logVars )
  
  MultCov=0.5

  if (is.null(RWCov)) RWCov=diag(rep(0.005, nParms))

  Nmc = BurnIn + NPosteriors + 2
  ibreaks = floor(Nmc/NcovUpdates) # no sims for each update to cov matrix
  BurnIn2 = floor(BurnIn/2)
  
  # some type sensitivity due to direct C-calls (below)
  nt = as.integer( dim(y)[1] )
  n = as.integer(n)
  nn = n*n
  NF = nn  # no. fourier components.. at present fixed at all
  Z = spate::spate.init(n=n, T=nt, NF=NF)
  indCosOnly = 1:Z$ns
  NFc = (NF-Z$ns)/2

  w2 = apply(Z$wave^2, 2, sum)
  d = 2

  #precompute a few spectral constants
  pid2 = pi^(d/2)
  nunu2 = 2^(nu - 1) * nu 
  nu2 = 2 * nu
  nud2 = nu + d/2
  ns = Z$ns

  pars = array(0.0, c(nParms, Nmc))
  dimnames(pars)[[1]] = pnames

  pars[, 1] = pars[, 2] = SV # the first is s throwaway to bootstrap the rest

  indNA = is.na(y)
  nNA = length( which(indNA) )
  if (nNA > 0 ) y[indNA]=0 # start from 0 (as a residual, expected value is 0 .. N(0,tau) )

  ypred = array(0, c(nt, nn, NPosteriors ))
  yhat = wFT = y[] * 0.0
  
  AcRate = 0
  
  
  for (i in 2:Nmc) {

      j= i-1
   
      y[indNA] = yhat[indNA] + RcppZiggurat::zrnorm(nNA) * sqrt(pars["tau2",j]) # 4-5 X faster than rnorm

      pars[,i] = pars[,j]
      pars[logInd,i] = log(pars[logInd,i])

      pars[,i] = FastGP::rcpp_rmvnorm(1, as.matrix(RWCov), pars[,i]) 
      # pars[,i] = rmvnorm(1, mean=pars[,i], sigma=as.matrix(RWCov), method="chol")
      
      pars[logInd,i] = exp(pars[logInd,i])

      al = 0

     # NOTE: logInd forces lognormal, positive values      
      if (  
        pars["rho0",i] > 1e-9 & # to avoid divide by zero
        pars["rho1",i] > 1e-9 & # to avoid divide by zero
        pars["tau2",i] < yvar_max &
        pars["sigma2",i] < yvar_max &
        pars["alpha",i] >= 0 & 
        pars["alpha",i] <= pi/2 & 
        pars["gamma",i] >= 1e-9 &
        pars["gamma",i] <= 1e+9  &
        pars["zeta",i] >= 1e-9 )   {

          # posterior likelihood
          # these are (nearly) the default priors used in spate  .. the log(1) are uniform priors are turned off as they are a constant
          priors =  c(
    #          -log(1), # rho0,  i.e, uniform prior ...  
            -1/2 * log(pars["sigma2",i]), # sigma2,
    #          -log(1),   # zeta 
    #          -log(1),   # rho1
            -log(pars["gamma",i]), # gamma
    #          -log(1), # alpha, 
    #          -log(1),  # Pmux
    #          -log(1),  # Pmuy
            -1/2 * log(pars["tau2",i]) # tau2
          )

          wFT[] = .C("TSreal_fft_stm", n=as.integer(n), T= as.integer(nt), yh=as.double(spate::TSmat.to.vect(y)), inverse=1L, 
            indCos=as.integer(Z$indFFT$indCos), indW=as.integer(Z$indFFT$indW), indWCon=as.integer(Z$indFFT$indWCon), 
            NFc=as.integer(NFc), PACKAGE="stm")$yh

          # spectrum: matern spatial and temporal ar1 
          rho0_inv = 1/pars["rho0",i]
          spec = { nunu2 * rho0_inv^nu2 } / { pid2 * {rho0_inv^2 + w2}^nud2 }
          spec[1:Z$ns] = spec[1:Z$ns]/2
          spec = spec * nn / sum(spec) # normalize 

          G11C = 1 
          G11 = 1
          G12 = -1

          if (Diffusion) {
            Tr = cbind( c(cos(pars["alpha",i]), - pars["gamma",i] * sin(pars["alpha",i])),
                        c(sin(pars["alpha",i]),   pars["gamma",i] * cos(pars["alpha",i]))) / pars["rho1",i]
            Sig = solve(t(Tr) %*% Tr) 
            # Sig = .C("inverse_crossproduct", as.matrix(Tr) )  
            # Sig = solve(t(Tr) %*% Tr) = solve(crossprod(A))
            DiffDamp = - dt  * { colSums(Z$wave * crossprod(Sig, Z$wave) ) + pars["zeta",i] } 
            gg = exp(DiffDamp)
            spec = spec * {1 - gg^2} / {-2} / { DiffDamp/dt } 
            G11C = G11C * gg[indCosOnly]    
            G11 =  G11 * gg[Z$indFFT$indCos]   
            G12 =  G12 * gg[Z$indFFT$indCos]
          }

          if (Advection) {
            Adv = c( dt *  crossprod( c(pars["muX",i], pars["muY",i]), Z$wave ) ) # advection
            G11 = G11 * cos(Adv[Z$indFFT$indCos])   
            G12 = G12 * sin(Adv[Z$indFFT$indCos])
          }

         if (i==2){
            priors_j = priors
            G11C_j = G11C
            G11_j = G11
            G12_j = G12
            spec_j = spec
          }

          ffbs = .C("ffbs_spectral_stm", wFT=as.double(wFT), bw=as.double(FALSE), ll=as.double(TRUE), 
            specCosOnly=as.double(spec[indCosOnly]), G11C=as.double(G11C), specCosSine= as.double(spec[Z$indFFT$indCos]), 
            G11=as.double(G11), G12=as.double(G12), specAll=as.double(spec), tau2=as.double(pars["tau2",i]), 
            T=as.integer(nt), NFc=as.integer(NFc), ns=as.integer(Z$ns), PACKAGE="stm" )

          ffbs_j = .C("ffbs_spectral_stm", wFT=as.double(wFT), bw=as.double(FALSE), ll=as.double(TRUE), 
            specCosOnly=as.double(spec_j[indCosOnly]), G11C=as.double(G11C_j), specCosSine= as.double(spec_j[Z$indFFT$indCos]), 
            G11=as.double(G11_j), G12=as.double(G12_j), specAll=as.double(spec_j), tau2=as.double(pars["tau2",j]), 
            T=as.integer(nt), NFc=as.integer(NFc), ns=as.integer(Z$ns), PACKAGE="stm" )

          al = min(1, exp( sum( c(
            {priors - priors_j },
            {ffbs$ll - ffbs_j$ll }, 
            {sum(log(pars[logInd,i])) - sum(log(pars[logInd,j])) })  , 
            na.rm=TRUE) ) )

      }

      # mcmc jump / update
      if (runif(1)  <= al ) {
        # store for next update
        priors_j = priors
        G11C_j = G11C
        G11_j = G11
        G12_j = G12
        spec_j = spec

        yhat = spate::vect.to.TSmat( .C("TSreal_fft_stm", n=as.integer(n), T=as.integer(nt), yh=as.double(ffbs_j$wFT), inverse=0L, 
                indCos=as.integer(Z$indFFT$indCos), indW=as.integer(Z$indFFT$indW), indWCon=as.integer(Z$indFFT$indWCon),
                NFc=as.integer(NFc), PACKAGE="stm" )$yh, T=nt)

        AcRate = AcRate + 1

      } else {

        pars[,i] = pars[,j] # reset to original (ie. no jump)

      }

      if (i > BurnIn & i  <= {NPosteriors+BurnIn} ) {
        ypred[, , i - BurnIn] = yhat 
      } 

      if ( (i > 100) & (AcRate/i <= 0.01) ) {
          message(paste("Acceptance rate for hyperparameters random walk is less than 1% after", 
            i, "iterations. Because of this, the proposal covariance matrix is divided by 100. \n"))
          RWCov = RWCov/10 
      }

      if (i > BurnIn2 ) {
        if ( (i - BurnIn2)%%ibreaks == 0 ) {
          if ( i-BurnIn2 <= BurnIn ) {
            parsC = pars[,(i -BurnIn2):i]
          } else {
            parsC = pars[,(i -BurnIn):i]
          }
          parsC[logInd, ] = log(parsC[logInd, ])
          RWCovP = MultCov * cov(t(parsC))
          eigenv = eigen(RWCovP)$value
          if (is.double(eigenv) & sum(eigenv <= 0) == 0) {
            RWCov = RWCovP
#              message("Estimated proposal covariance for hyperparameters: \n")
#              print(signif(RWCov, digits=2))
            rm(parsC, eigenv, RWCovP); gc()
          } 
        }
      }
 
   } #end for

  spateMCMC = list(Post=pars[,-c(1,Nmc)], ypred=ypred, RWCov=RWCov, BurnIn=BurnIn, 
    Padding=FALSE, indEst=indEst, nu=nu, DataModel="Normal")
  class(spateMCMC)="spateMCMC"


  return(spateMCMC)


}

