
stm_LaplacesDemon_spatemodel = function(Data=list()) {

  #Data$mon.names = c( "LP", paste0("yhat[",1:Data$N,"]" ) )
  Data$mon.names = c( "LP" )
  
  Data$parm.names = as.parm.names(list(
    rho0=0, sigma2=0, zeta=0, rho1=0, gamma=0, alpha=0, muX=0, muY=0, tau2=0 
  ))

  Data$pos = list(
    rho0 = grep("rho0", Data$parm.names),
    sigma2 = grep("sigma2", Data$parm.names),
    zeta = grep("zeta", Data$parm.names),
    rho1 = grep("rho1", Data$parm.names),
    gamma = grep("gamma", Data$parm.names),
    alpha = grep("alpha", Data$parm.names),
    muX = grep("muX", Data$parm.names),
    muY = grep("muY", Data$parm.names),
    tau2 = grep("tau2", Data$parm.names)
  )

  Data$PGF = function(Data) {
    # Parameter Generating Function: initial values .. get them near the center of mass
    rho0=runif(1,0,100)
    sigma2=rgamma (1, 1, 1)
    zeta=runif(1, 1e-8, 1e3)
    rho1=runif(1,0,100)
    gamma=runif(1, 1e-6, 1e6)
    alpha=runif(1, 0, pi/2)
    muX=runif(1, -0.5, 0.5)
    muY=runif(1, -0.5, 0.5)
    tau2=rgamma (1, 1, 1)
    return( c( rho0, sigma2, zeta, rho1, gamma, alpha, muX, muY, tau2 ))
  }
  Data$PGF  = compiler::cmpfun(Data$PGF)
  
  Data$Model = function(parm, Data) {
    rho0 = parm[Data$pos$rho0]= LaplacesDemonCpp::interval_random(parm[Data$pos$rho0], Data$eps, 100, 100 )
    sigma2 = parm[Data$pos$sigma2]= LaplacesDemonCpp::interval_random(parm[Data$pos$sigma2], Data$eps, Data$yvar, Data$yvar )
    zeta =  parm[Data$pos$zeta] = LaplacesDemonCpp::interval_random(parm[Data$pos$zeta], Data$eps, 1e3, 1e3 )
    rho1 = parm[Data$pos$rho1]= LaplacesDemonCpp::interval_random(parm[Data$pos$rho1], Data$eps, 100, 100 )
    gamma =  parm[Data$pos$gamma] = LaplacesDemonCpp::interval_random(parm[Data$pos$gamma], Data$eps, 1/Data$eps, 1/Data$eps )
    alpha = parm[Data$pos$alpha]= LaplacesDemonCpp::interval_random(parm[Data$pos$alpha], Data$eps, pi/2, pi/2 )
    muX =  parm[Data$pos$muX] = LaplacesDemonCpp::interval_random(parm[Data$pos$muX], -0.5, 0.5, 1 )
    muY =  parm[Data$pos$muY] = LaplacesDemonCpp::interval_random(parm[Data$pos$muY], -0.5, 0.5, 1 )
    tau2 = parm[Data$pos$tau2] = LaplacesDemonCpp::interval_random(parm[Data$pos$tau2], Data$eps, Data$yvar, Data$yvar )


    Tr = cbind( c(cos(alpha), - gamma * sin(alpha)),  c(sin(alpha),   gamma * cos(alpha))) / rho1
    Sig = solve(t(Tr) %*% Tr)

    # diffusion and damping
    DiffDamp = - Data$dt  * { colSums(Data$M$wave * crossprod(Sig, Data$M$wave) ) + zeta } 
    Adv = c( Data$dt *  crossprod( c(muX, muY), Data$M$wave ) ) # advection

      # spectrum: matern spatial and temporal ar1 
    rho0_inv = 1/rho0
    spec = { Data$nunu2 * rho0_inv^Data$nu2 } / { Data$pid2 * {rho0_inv^2 + Data$w2}^Data$nud2 }
    spec[Data$indCosOnly] = spec[Data$indCosOnly]/2
    spec = spec * Data$N / sum(spec) * {1 - exp(2 * DiffDamp)}/{-2}/ ( DiffDamp/Data$dt ) # normalize and include diffusion/advection

    RcppZiggurat::zsetseed( floor(runif(1)*1000) )
    fbs = .C( "ffbs_spectral_stm_oneshot", 
      yh=as.double(spate::TSmat.to.vect(Data$y)), 
      ll=as.double(1.0),  # dummy variable .. make sure it is a float 
      tau2=as.double(tau2),
      indCosOnly=as.integer(Data$indCosOnly),
      indCos=as.integer(Data$indCos),
      indW=as.integer(Data$indW),
      indWCon=as.integer(Data$indWCon),
      spec=as.double(spec),
      DiffDamp=as.double(DiffDamp),
      Adv=as.double(Adv),
      T=as.integer(Data$T), 
      n=as.integer(Data$n), 
      NF=as.integer(Data$NF), 
      NFc=as.integer(Data$NFc), 
      ns=as.integer(Data$ns),
      PACKAGE="stm" )

    LL = fbs$ll ## Log Likelihood

    # these are (nearly) the default priors used in spate
    LP = sum(LL, 
      -log(1), # rho0,  i.e, uniform prior ...  
      -1/2 * log(sigma2), # sigma2,
      -log(1),   # zeta 
      -log(1),   # rho1
      -log(gamma), # gamma
      -log(1), # alpha, 
      -log(1),  # Pmux
      -log(1),  # Pmuy
      -1/2 * log(tau2), # tau2
      parm[c(1, 2, 3, 4, 5, 9)] ) ### Log-Posterior

    Modelout = list(LP=LP, Dev=-2*LL, Monitor=c(LP), yhat=fbs$yh, parm=parm)
    
    rm(fbs, spec, DiffDamp, Adv); gc()
    
    return(Modelout)
  }
  
  Data$Model.ML  = compiler::cmpfun( function(...) (Data$Model(...)$Dev / 2) )  # i.e. - log likelihood
  Data$Model.PML = compiler::cmpfun( function(...) (- Data$Model(...)$LP) ) #i.e., - log posterior 
  Data$Model = compiler::cmpfun(Data$Model) #  byte-compiling for more speed .. use RCPP if you want more speed

  return(Data)

} 
