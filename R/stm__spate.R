
stm__spate = function( p, dat, pa, sloc, distance, nu, phi, varObs, varSpatial ) {

  # require(spate) #\\ SPDE solution via FFT using the spate library
  # sloc=Sloc[Si,]; distance=stm_distance_cur

  sdTotal=sd(dat[,p$variable$Y], na.rm=T)
  ndata = nrow(dat)
  
  TS = stm_timeseries_smooth(p=p, dat=dat, sloc=sloc, distance=distance )
  if( is.null(TS)) return(NULL)
  if ( nrow( TS$datgridded) < p$n.min ) return(NULL)

  # shorten by 1 row and column to make it even
  nsq = TS$pa_w_n-1
  w = matrix( TS$xM[,1:nsq, 1:nsq ], nrow=p$nt )

  SV = c(rho0 = 0.1, 
        sigma2 = varSpatial, zeta = 0.1, rho1 = 0.1, gamma = 0.5, alpha = 0.5, 
        muX = 0, muY = 0, tau2 = varObs)


  if ( p$stm_spate_method=="mcmc" ) {
      
      ntotalSims = p$stm_spate_nburnin + p$stm_spate_nposteriors

      g = try( spate::spate.mcmc( y=w, n=nsq, Padding=FALSE, trace=TRUE, Nmc=ntotalSims, SV=SV,
        adaptive=TRUE, Separable=FALSE, Drift=TRUE, Diffusion=TRUE,  
        BurnIn=p$stm_spate_nburnin, BurnInCovEst=floor(p$stm_spate_nburnin*0.5), NCovEst=p$stm_spate_nposteriors ), silent=TRUE)

      if ("try-error" %in% class(g)) {
        g = try( spate::spate.mcmc( y=w, n=nsq, Padding=FALSE, trace=FALSE, Nmc=ntotalSims, 
          adaptive=TRUE, Separable=FALSE, Drift=TRUE, Diffusion=TRUE,  
          BurnIn=p$stm_spate_nburnin, BurnInCovEst=floor(p$stm_spate_nburnin*0.75), NCovEst=p$stm_spate_nposteriors ), silent=TRUE) # longer burn-in (1000 is default) and alternate rnd seed
      }

      if ("try-error" %in% class(g)) {
        g = try( spate::spate.mcmc( y=w, n=nsq, Padding=FALSE, trace=FALSE, seed=1, Nmc=ntotalSims,
          adaptive=TRUE, Separable=FALSE, Drift=TRUE, Diffusion=TRUE,  
          BurnIn=p$stm_spate_nburnin, BurnInCovEst=p$stm_spate_nburnin, NCovEst=p$stm_spate_nposteriors ), silent=TRUE) # longer burn-in (1000 is default) and alternate rnd seed
      }

      if ("try-error" %in% class(g)) return(NULL)

      spp <- spate::spate.predict(y=w, tPred=(1:p$nt), 
        spateMCMC=g, Nsim=p$stm_spate_nposteriors, BurnIn=p$stm_spate_nburnin, DataModel="Normal", trace=FALSE ) # nu=nu, defulat is to assume nu =1 
      #  DimRed=TRUE, NFour=101
      rm(w); gc()

      # must be same order as p$statsvars
      pmean_vars = c( "rho_0", "zeta", "rho_1", "gamma", "alpha", "mu_x", "mu_y", "sigma^2", "tau^2" ) #reorder
      psd_vars   = c( "rho_0.sd", "zeta.sd", "rho_1.sd", "gamma.sd", "alpha.sd", "mu_x.sd", "mu_y.sd" ) # drop sd of variance terms

      pmean =  apply(g$Post, 1, mean)[pmean_vars] 
      psd =  apply(g$Post, 1, sd)[pmean_vars] 

      names(pmean) = pmean_vars
      names(psd) = paste( pmean_vars, "sd", sep=".")
      psd = psd[psd_vars]


  } else if ( p$stm_spate_method=="likelihood" ) {

      warning( "This is just to show the method for param estimation via max likelihood .. it does not handle missing data..")

      spateFT <- spate::spate.init(n=nsq, T=p$nt, NF=nsq*nsq )
      wFT <- spate::real.fft.TS(TSmat.to.vect(w), n = nsq, T = p$nt, inv = TRUE, indFFT = spateFT$indFFT)

      wv = spate::wave.numbers(nsq)
      ##Specify hyper-parameters
      par = c(rho0=0.1,sigma2=0.2,zeta=0.5,rho1=0.1,gamma=2,alpha=pi/4,muX=0.2,muY=-0.2,tau2=0.01)
      logInd=c("rho0", "sigma2", "zeta", "rho1", "gamma", "tau2") # logInd=c(1,2,3,4,5,9)

      ##Initial values for optim. This takes a couple of seconds.
      parI = SV
      parI[logInd] = log(parI[logInd]) ##Transform to log-scale
      ##Fourier transform needs to be done only once

      g = try( optim(par=parI, spateML, control=list(trace=TRUE,maxit=1000), wFT=wFT, method="L-BFGS-B",
        lower=c(-10,-10,-10,-10,-10,0,-0.5,-0.5,-10), upper=c(10,10,10,10,10,pi/2,0.5,0.5,10),
        wave=wv$wave, indCos=wv$indCos, hessian=TRUE, n=nsq, T=p$nt),  silent=TRUE )

      if ("try-error" %in% class(g)) return(NULL)
    
      spp = spateSample( spate.mle=g, w=w, n=nsq, T=p$nt, Nsim=p$stm_spate_nposteriors )
      rm(w); gc()

      pmean_vars = c( "rho_0", "zeta", "rho_1", "gamma", "alpha", "mu_x", "mu_y", "sigma^2", "tau^2" ) #reorder
      psd_vars   = c( "rho_0.sd", "zeta.sd", "rho_1.sd", "gamma.sd", "alpha.sd", "mu_x.sd", "mu_y.sd" ) # drop sd of variance terms

      pmean =  g$par
      pmean[logInd] = exp(pmean[logInd])
      
      psd = sqrt(diag(solve(g$hessian)))[pmean_vars]
 
      names(pmean) = pmean_vars
      names(psd) = paste( pmean_vars, "sd", sep=".")
      psd = psd[psd_vars]




  } else if ( p$stm_spate_method=="mcmc_fast" ) {

      g = try( spate_mcmc_fast( y=w, yvar=sdTotal^2, n=nsq, NPosteriors=p$stm_spate_nposteriors, BurnIn=p$stm_spate_nburnin, NcovUpdates=p$stm_spate_nCovUpdates, SV=SV ), silent=TRUE)

      if (is.null(g)) return(NULL)

      if ("try-error" %in% class(g)) {
        g = try( spate_mcmc_fast( y=w, yvar=sdTotal^2, n=nsq, NPosteriors=p$stm_spate_nposteriors, BurnIn=floor(1.5*p$stm_spate_nburnin), NcovUpdates=p$stm_spate_nCovUpdates), silent=TRUE) # longer burn-in (1000 is default) and alternate rnd seed
      }

      if ("try-error" %in% class(g)) {
        g = try( spate_mcmc_fast( y=w, yvar=2*sdTotal^2, n=nsq, seed=1, NPosteriors=p$stm_spate_nposteriors, BurnIn=floor(2*p$stm_spate_nburnin), NcovUpdates=p$stm_spate_nCovUpdates), silent=TRUE) 
      }

      if ("try-error" %in% class(g)) return(NULL)
      rm(w); gc()


      # must be same order as p$statsvars
        pmean_vars = dimnames(g$Post)[[1]]
        psd_vars   = paste( pmean_vars, ".sd", sep="") 

        pmean =  apply(g$Post, 1, mean)
        names(pmean) = pmean_vars

        psd =  apply(g$Post, 1, sd)
        names(psd) = psd_vars

        spp = g$ypred
        
        logInd=c("rho0", "sigma2", "zeta", "rho1", "gamma", "tau2") # logInd=c(1,2,3,4,5,9)
        pmean[logInd] = exp(pmean[logInd])
  
  }

  rm(g); gc()

  if (0){

       zlim=range(dat[, p$variables$Y ])
        plot.new()
        par( mfrow=c(1,3))
        YPREDS = apply(spp, c(1,2), mean)
        dmbas0 = MBA::mba.surf( dat[, c( p$variables$LOCS, p$variables$Y) ], 300, 300, extend=TRUE)$xyz.est
        for (i in 1:p$nt) {
           pause(1)
          # data
          di = which( floor(dat[ , p$variables$TIME ]) == floor(p$prediction.ts[i]) )
          if (length(di) > 5) {
            dmbas = MBA::mba.surf( dat[di, c( p$variables$LOCS, p$variables$Y) ], 300, 300, extend=TRUE)$xyz.est
           } else { 
            dmbas = dmbas0 
          }
          surface(dmbas, zlim=zlim)

          # grid boosted
          xi = which( TS$datgridded[ , p$variables$TIME ] == p$prediction.ts[i] )
          mbas = MBA::mba.surf( na.omit(TS$datgridded[xi, c( p$variables$LOCS, p$variables$Y) ]), 300, 300, extend=TRUE)$xyz.est
          surface(mbas, zlim=zlim)

         # spate   
          oo = expand.grid( 1:nsq, 1:nsq)
          oo = as.data.frame( cbind(oo, YPREDS[i,] ) )
          oom = MBA::mba.surf( oo, 300, 300, extend=TRUE)$xyz.est
          surface(oom, zlim=zlim)
        }  # end for  

  } # end debug


  pa$id = array_map( "3->1", 
    coords = round( cbind( 
      ( (pa[,p$variables$TIME ] - p$prediction.ts[1] ) / p$tres) + 1 ,
      ( TS$windowsize.half + (pa[,p$variables$LOCS[1]] - sloc[1]) / p$pres) + 1, 
      ( TS$windowsize.half + (pa[,p$variables$LOCS[2]] - sloc[2]) / p$pres) + 1)), 
    dims=TS$adims )

  iii = which(pa$id >= 1 & pa$id <= prod(TS$adims) )
  if (length(iii) < 10) return(NULL) 

  # means 
  TS$xM[,1:nsq,1:nsq] = apply(spp, c(1,2), mean)
  pa$mean = NA
  pa$mean[iii] = TS$xM[pa$id[iii] ]

  TS$datgridded$mean = NA
  TS$datgridded$mean = TS$xM[TS$datgridded$id]

  # sd
  TS$xM[,1:nsq,1:nsq] = apply(spp, c(1,2), sd)
  pa$sd = NA
  pa$sd[iii] = TS$xM[pa$id[iii]]

  rm(spp); gc()

  # plot(mean ~ z , TS$datgridded)

  ss = lm( TS$datgridded$mean ~ TS$datgridded[,p$variables$Y], na.action=na.omit )
  if ( "try-error" %in% class( ss ) ) return( NULL )
  rm(TS); gc()

  rsquared = summary(ss)$r.squared
  if (rsquared < p$stm_rsquared_threshold ) return(NULL)

  stm_stats = c(list( sdTotal=sdTotal, rsquared=rsquared, ndata=ndata), pmean, psd)  
  
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scalse=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, stm_stats=stm_stats ) )  

}

