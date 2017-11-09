
spateSample = function( spate.mle, n, T, nu=1, Nsim=100, w=NULL, NFour=NULL ) {

    if (0) {
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
      wFT = spateMLE$wFT
      ##ML estimation using optim, takes a couple of seconds .. fast version
      spateMLE = optim(par=parI, spateML, control=list(trace=TRUE,maxit=1000), wFT=wFT, method="L-BFGS-B",
          lower=c(-10,-10,-10,-10,-10,0,-0.5,-0.5,-10), upper=c(10,10,10,10,10,pi/2,0.5,0.5,10),
          wave=wv$wave, indCos=wv$indCos, hessian=TRUE, n=n, T=T)
      mle = spateMLE$par
      mle[logInd] = exp(mle[logInd])
      sd=sqrt(diag(solve(spateMLE$hessian)))

      sptime.sample = spateSample( spate.mle=spateMLE, w=w, n=n, T=T, Nsim=10 )

    }

    NF <- n * n
    if (!is.null(NFour)) NF <- NFour

    sPred <- 1:(n^2)

    pnames = c("rho0", "sigma2", "zeta", "rho1", "gamma", "alpha", "muX", "muY", "tau2")
    if (exists("SIGMA", spate.mle )) {
      SIGMA=spate.mle$SIGMA
    } else if (exists("hessian", spate.mle )) {
      SIGMA = solve(spate.mle$hessian) # VAR-COV matrix
    }

    parh = mvnfast::rmvn( Nsim, spate.mle$par, SIGMA )
    colnames(parh) = pnames
    spateFT <- spate::spate.init(n = n, T=T, NF=NF)
    ypred <- array(0, c(T, length(sPred), Nsim))
    for (i in 1:Nsim) {
      sim = spateML( par=parh[i,], wFT=spate.mle$wFT, wave=spateFT$wave, indCos=spateFT$indCos, ns=spateFT$ns,
        n=n, T=T,  nu=nu, retvalue="sample")
      ypred[,,i] <- spate::vect.to.TSmat( spate::real.fft.TS( spate::TSmat.to.vect(sim), n=n, T=T, inv=FALSE, indFFT=spateFT$indFFT) )
    }
    return(ypred)
}


