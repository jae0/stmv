
spateML = function (
  par=NULL, w=NULL, wFT=NULL, x=NULL, n, T, NF=n*n, ns=4, nu=1, dt=1,
  logInd = c("rho0", "sigma2", "zeta", "rho1", "gamma", "tau2"),
  wave=wave.numbers(n)$wave, indCos=wave.numbers(n)$indCos, retvalue="loglik" ) {

  if (0) {
    # usage:

    # copied from spate::loglike
    require(spate)
    n = 20 # must be even
    T = 20
    wv = spate::wave.numbers(n)
    ##Specify hyper-parameters
    par = c(rho0=0.1,sigma2=0.2,zeta=0.5,rho1=0.1,gamma=2,alpha=pi/4,muX=0.2,muY=-0.2,tau2=0.01)
    logInd=c("rho0", "sigma2", "zeta", "rho1", "gamma", "tau2") # logInd=c(1,2,3,4,5,9)
    ##Simulate data
    spateSim = spate.sim(par=par,n=n,T=T,seed=4)
    w = spateSim$w
    ##Initial values for optim. This takes a couple of seconds.
    parI = c(rho0=0.2,sigma2=0.1,zeta=0.25,rho1=0.01,gamma=1,alpha=0.3,muX=0,muY=0,tau2=0.005)
    parI[logInd] = log(parI[logInd]) ##Transform to log-scale
    ##Fourier transform needs to be done only once
    wFT = TSmat.to.vect( real.fft.TS(w, n=n, T=T, inv=TRUE) )
    ##ML estimation using optim, takes a couple of seconds
    ##Load the precomputed object a line below to save time
    spateMLE0 <-
         optim(par=parI,loglike,control=list(trace=TRUE,maxit=1000),wFT=wFT,method="L-BFGS-B",
              lower=c(-10,-10,-10,-10,-10,0,-0.5,-0.5,-10),
              upper=c(10,10,10,10,10,pi/2,0.5,0.5,10),negative=TRUE,
              logScale=TRUE,hessian=TRUE,n=n,T=T)

    # this is using a faster version of spate::loglike
    spateMLE = optim(par=parI, spateML, control=list(trace=TRUE,maxit=1000), wFT=wFT, method="L-BFGS-B",
        lower=c(-10,-10,-10,-10,-10,0,-0.5,-0.5,-10), upper=c(10,10,10,10,10,pi/2,0.5,0.5,10),
        wave=wv$wave, indCos=wv$indCos, hessian=TRUE, n=n, T=T)
    mle = spateMLE$par
    mle[logInd] = exp(mle[logInd])
    sd=sqrt(diag(solve(spateMLE$hessian)))

    MleConfInt = data.frame(array(0,c(4,9)))
    colnames(MleConfInt) = names(par)
    rownames(MleConfInt) = c("True","Estimate","Lower","Upper")
    MleConfInt[1,] = par
    MleConfInt[2,] = mle
    MleConfInt[3,] = spateMLE$par-2*sd
    MleConfInt[4,] = spateMLE$par+2*sd
    MleConfInt[c(3,4),logInd] = exp(MleConfInt[c(3,4),logInd])
    cat("\n")
    round(MleConfInt,digits=4)


     GV <- get.propagator(

                          wave=spateFT$wave, 
                          indCos=spateFT$indCos, 
                          zeta=parV[3], 
                          rho1=parV[4], 
                          gamma=parV[5], 
                          alpha=parV[6], 
                          muX=parV[7], 
                          muY=parV[8], 
                          dt=1, 
                          ns=spateFT$ns)

    lpV <- lp
    ffbs <- ffbs(wh, lp=lpV, G=GV, Sigma=diag(specV), 
      H=IPhi, Omega=diag(rep(parV[9], dim(wh)[2])), 
      lglk=TRUE, BwSp=TRUE)
    mllV <- ffbs$ll
    mll <- ffbs(wh, lp=lp, G=G, Sigma=diag(spec), 
      H=IPhi, Omega=diag(rep(parh[9, i], dim(wh)[2])), 
      lglk=TRUE, BwSp=FALSE)$ll

  }

  # note:: this is spate::loglike broken down into basic parts for speed

  par[logInd] = exp(par[logInd])

  error = FALSE
  if (nu < 0) error=TRUE # print("Error: nu needs to be positive")
  if (par["sigma2"] < 0) error = TRUE # print("Error: sigma2 needs to be positive")
  if (par["zeta"] < 0) error = TRUE # print("Error: zeta needs to be positive")
  if (par["gamma"] < 0) error = TRUE # print("Error: gamma needs to be positive")
  if (par["alpha"] < 0 & par["alpha"] >= (pi/2)) error = TRUE # print("Error: alpha needs to be between 0 and pi/2")
  if (error) return(1e12) # large positive value

  if (!is.null(x)) {
      lp = apply( X=x, MARGIN=c(2, 3), FUN=function(x,beta){ x %*% beta}, beta=par[-c(1:9)])  # llik of predictors ...
  } else {
      lp = 0
  }

  # spate::get.propagator.vec .. copied and modified from
  if ( par["rho1"] == 0) {
    Sig <- cbind(c(0, 0), c(0, 0))
  } else {
    Tr <- cbind( c(cos(par["alpha"]), - par["gamma"] * sin(par["alpha"])),
                 c(sin(par["alpha"]),   par["gamma"] * cos(par["alpha"]))) / par["rho1"]
    Sig <- solve(t(Tr) %*% Tr)
  }

  # diffusion and damping
  DiffDamp0 = - { apply(wave*Sig %*% wave, 2, sum) + rep(par["zeta"],NF) } 
  DiffDamp = dt * DiffDamp0

  Adv <- dt * c(par["muX"], par["muY"]) %*% wave  # advection
  G11C <- exp(DiffDamp[1:ns])
  G11 <-  exp(DiffDamp[indCos]) * cos(Adv[indCos])
  G12 <- -exp(DiffDamp[indCos]) * sin(Adv[indCos])

  # copied from spate::innov.spec
  w2 <- apply(wave^2, 2, sum)
  d <- 2
  spec <- {(2^(nu - 1)) * nu * ((1/par["rho0"])^(2 * nu)) } / {(pi^(d/2)) * ((1/par["rho0"])^2 + w2)^(nu + d/2)}
  spec[1:ns] <- spec[1:ns]/2
  spec <- spec * (n * n)/sum(spec)  # normalize spectrum to unit variance
  spec = spec * (1 - exp(2 * DiffDamp))/(-2)/DiffDamp0


  # copied from spate::ffbs_spectral
  if (retvalue=="loglik") {
    ffbsC <- .C("ffbs_spectral", wFT=as.double(wFT), bw=as.double(sum(FALSE)),
      ll = as.double(sum(TRUE)), as.double(spec[1:ns]), as.double(G11C),
      as.double(spec[indCos]), as.double(G11), as.double(G12),
      as.double(spec), as.double(par["tau2"]), as.integer(T), as.integer((NF-ns)/2), as.integer(ns))
    return( -ffbsC$ll )
  }

  if (retvalue=="sample") {
    ffbsC <- .C("ffbs_spectral", wFT=as.double(wFT), bw=as.double(sum(TRUE)),
      ll = as.double(sum(FALSE)), as.double(spec[1:ns]), as.double(G11C),
      as.double(spec[indCos]), as.double(G11), as.double(G12),
      as.double(spec), as.double(par["tau2"]), as.integer(T), as.integer((NF-ns)/2), as.integer(ns))
    return( vect.to.TSmat(ffbsC$wFT, T=T) )
  }

}

