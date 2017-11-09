


spate_ffbs = function (yshort, params, IPhi, Z, N=dim(yshort)[2], T=dim(yshort)[1], dt=1, nu=1, lglk=FALSE, BwSp=TRUE, filt=FALSE, pred=TRUE ) {
  
  # based on spate::ffbs, all in one index method
  
  if ( params["rho1"] == 0) {
    Sig = cbind(c(0, 0), c(0, 0))
  } else {
    Tr = cbind( c(cos(params["alpha"]), - params["gamma"] * sin(params["alpha"])),
                c(sin(params["alpha"]),   params["gamma"] * cos(params["alpha"]))) / params["rho1"]
    Sig = solve(t(Tr) %*% Tr)
  }

  # diffusion and damping
  DiffDamp0 = - { apply(Z$wave*Sig %*% Z$wave, 2, sum) + rep(params["zeta"],NF) } 
  DiffDamp = dt * DiffDamp0

  Adv = dt * c(params["muX"], params["muY"]) %*% Z$wave  # advection

  G <- matrix(0,ncol=NF,nrow=NF) # "propagator matrix" that accounts for advection and diffusion
  diag(G) <- c(exp(DiffDamp[1:Z$ns]),exp(DiffDamp[(Z$ns+1):NF])*cos(Adv[(Z$ns+1):NF]))
  diag(G[Z$indCos,  Z$indCos+1]) <- -exp(DiffDamp[Z$indCos])*sin(Adv[Z$indCos])
  diag(G[Z$indCos+1,Z$indCos  ]) <- +exp(DiffDamp[Z$indCos])*sin(Adv[Z$indCos])

  # spectrum: matern spatial and temporal ar1 
  spec = { 2^(nu - 1) * nu * {1/params["rho0"]}^{2 * nu} } / { (pi^(d/2)) * {{1/params["rho0"]}^2 + w2}^{nu + d/2} }
  spec[1:Z$ns] = spec[1:Z$ns]/2
  spec = spec * nn / sum(spec) * (1 - exp(2 * DiffDamp))/(-2)/DiffDamp0 # normalize and include diffusion/advection
 
  Sigma = diag(spec)

  Omega = diag(rep(params["tau2"], dim(yshort)[2]))


  tIPhi <- t(IPhi)
  mtt1 <- array(0, c(T, NF))
  mtt  <- array(0, c(T + 1, NF))
  Rtt1 <- array(0, c(T, NF, NF))
  Rtt  <- array(0, c(T + 1, NF, NF))
  simAlpha <- array(0, c(T + 1, NF))
  Innt <- array(0, c(T, N))
  PrecInnt <- array(0, c(T, N, N))
  Rtt[1, , ] <- as.matrix(Sigma)
  
  for (t in 1:T) {
    print(t)
    mtt1[t, ]   <- G %*% mtt[t, ]
    Rtt1[t, , ] <- Sigma + G %*% Rtt[t, , ] %*% t(G)
    Rtt1[t, , ] <- (Rtt1[t, , ] + t(Rtt1[t, , ]))/2
    TM1 <- IPhi %*% Rtt1[t, , ]
    Mti <- solve(Omega + TM1 %*% tIPhi)
    PrecInnt[t, , ] <- Mti
    TM2 <- Rtt1[t, , ] %*% tIPhi %*% Mti
    Innt[t, ] <- yshort[t, ] - IPhi %*% mtt1[t, ]
    mtt[t + 1, ] <- mtt1[t, ] + TM2 %*% (Innt[t, ])
    Rtt[t + 1, , ] <- Rtt1[t, , ] - TM2 %*% TM1
    Rtt[t + 1, , ] <- (Rtt[t + 1, , ] + t(Rtt[t + 1, , ]))/2
  }
  
  if (BwSp) {
    Rtt[T + 1, , ] = (Rtt[T + 1, , ] + t(Rtt[T + 1, , ]))/2
    simAlpha[T + 1, ] <- rmvnorm(1, mean = mtt[T + 1, ], sigma = Rtt[T + 1, , ], method = "chol")
    for (t in T:1) {
      tm <- Rtt[t, , ] %*% t(G) %*% solve(Rtt1[t, , ])
      Rt <- Rtt[t, , ] - tm %*% G %*% Rtt[t, , ]
      Rt <- (Rt + t(Rt))/2
      mt <- mtt[t, ] + tm %*% (simAlpha[t + 1, ] - mtt1[t, ])
      simAlpha[t, ] <- rmvnorm(1, mean = mt, sigma = Rt, method = "chol")
    }
  }

  if (lglk) {
    ll <- 0
    llt = rep(0, T)
    for (t in 1:T) {
      llt[t] = determinant(PrecInnt[t, , ], logarithm = TRUE)$modulus[[1]] - t(Innt[t, ]) %*% PrecInnt[t, , ] %*% Innt[t, ]
      ll <- sum( llt, na.rm=TRUE ) 
    }
    ll <- ll/2 - T * NF * log(2 * pi)/2
  }
  
  ret <- list()
  
  if (BwSp | pred){
    simAlpha=simAlpha[-1, ]
    if(class(simAlpha)=="numeric") simAlpha=t(matrix(simAlpha))
    
    if (BwSp) {
      ret <- c(ret, list(simAlpha = simAlpha))
    }
  
    if ( pred ) {
      alphahC <- array(0, c(T, n * n))
      alphahC[, Z$IndFour] <- simAlpha
#      oo = which(!is.finite(alphahC))
#      if (length(oo) > 0)  alphahC[oo] = 0
      y <- vect.to.TSmat(real.fft.TS(TSmat.to.vect(alphahC), n = n, T = T, inv = FALSE, indFFT = Z$indFFT), T = T)
      y = t(IPhi %*% t(simAlpha))
      ret <- c(ret,  ypred=y  )
    }
  } 

  if (lglk) 
    ret <- c(ret, list(ll = ll))
  
  if (filt) 
    ret <- c(ret, list(mtt = mtt[-1, ]))
  

  return(ret)
}
