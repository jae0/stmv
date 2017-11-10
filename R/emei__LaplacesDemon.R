
emei__LaplacesDemon = function( p, dat, pa, sloc, distance, nu, phi, varObs, varSpatial ) {

  # require(LaplacesDemonCpp)
  # sloc=Sloc[Si,]; distance=emei_distance_cur
 
  sdTotal = sd(dat[, p$variables$Y], na.rm=T) 
  ndata = nrow(dat)
  
  TS = emei_timeseries_smooth(p=p, dat=dat, sloc=sloc, distance=distance )

  SV = c(rho0 = 0.25, 
        sigma2 = mean(0.01, varSpatial, na.rm=TRUE), zeta = 0.25, rho1 = 0.25, gamma = 1, alpha = 0.1, 
        muX = 0, muY = 0, tau2 = mean(0.01, varObs, na.rm=TRUE) )
  

    # shorten by 1 row and column to make it even
    nsq = TS$pa_w_n-1
    w = matrix( TS$xM[,1:nsq, 1:nsq ], nrow=p$nt )
    nburnin = p$emei_spate_nmcmc - p$emei_spate_nposteriors


    if (0) {
      require(LaplacesDemonCpp)
      y=w; yvar=sdTotal^2; n=nsq; Nmc=p$emei_spate_nmcmc;SV=SV; Drift=TRUE; Diffusion=TRUE; NPosteriors=1000; BurnIn=100; NcovUpdates=20; RWCov=NULL; parh=NULL; indEst=1:9; dt=1; MultCov=0.5; ogInd=c(1, 2, 3, 4, 5, 9); nu=1; seed=1
    } 

   # modified to do both prediction and parameter estimation and speed up as much as possible
   # removed incidence matrix method and covariate predictors
   
    indNN = c(1, 2, 3, 4, 5, 9)[-match(logInd, c(1, 2, 3, 4, 5, 9))]

    # some type sensitivity due to direct C-calls (below)
    T = as.integer( dim(y)[1] )
    n = as.integer(n)
    NF = n*n  # no. fourier components.. at present fixed at all
    if (!is.null(seed)) set.seed(seed)
   
      
        Z = spate::spate.init(n=n, T=T, NF=NF)
 
        d = 2
         indNA = is.na(y)
         nNA = length( which(indNA) )
          if (nNA > 0 ) y[indNA]=0 


        # fast spin up
         g = spate_mcmc_fast( y=y, yvar=Data$yvar, n=n, NPosteriors=500, BurnIn=500, NcovUpdates=10, SV=SV )
         y= apply( g$ypred, c(1,2), median )
         Initial.Values =  apply( g$Post[, 501:1000 ], 1, median )
         rm(g); gc()

          Data = list(
            eps = 1e-6, epsb = 1e-1, d=d, N = n*n*T,  # required for LaplacesDemon
            NFc = (NF-Z$ns)/2, T=as.integer(T), 
            n=as.integer(n), 
            ns=as.integer(Z$ns), NF=as.integer(NF), 
            dt=1, w2=apply(Z$wave^2, 2, sum), M=list(wave=Z$wave),  #put matrix inside a list to hide it from LaplacesDEmon
            pid2 = pi^(d/2), nunu2 = 2^(nu - 1) * nu, 
            nu2 = 2 * nu, nud2 = nu + d/2, yvar= sdTotal^2, indCosOnly = as.integer(1:Z$ns), indCos=as.integer(Z$indFFT$indCos), 
            indW=as.integer(Z$indFFT$indW), indWCon=as.integer(Z$indFFT$indWCon), 
            y=spate::TSmat.to.vect(y),
            indNA=indNA,
            nNA=nNA 
          )
          source( "~/bio/emei/R/emei_LaplacesDemon_spatemodel.R")
          source( "~/bio/emei/R/emei_LaplacesDemon_Specification.R")



          Data = emei_LaplacesDemon_spatemodel(Data)
          Data$yhat=Data$y[]*0


          str(Data$Model( parm=Data$PGF(Data), Data ) ) # test to see if return values are sensible



          Model = Data$Model


     
        str(Data$Model( parm=Data$PGF(Data), Data ) ) # test to see if return values are sensible


        # maximum likelihood solution .. kind of slow
        # f.ml = optim( par=Data$PGF(Data), fn=Data$Model.ML, Data=Data, control=list(maxit=5000, trace=1), method="BFGS"  )
        # names(f.ml$par ) = Data$parm.names

        # penalized maximum likelihood .. better but still a little unstable depending on algorithm
        # f.pml = optim( par=Data$PGF(Data), fn=Data$Model.PML, Data=Data,  control=list(maxit=1000, trace=1), method="BFGS" , hessian=FALSE )
        # names(f.pml$par ) = Data$parm.names
        #print(sqrt( diag( solve(f.pml$hessian) )) ) # assymptotic standard errors

        f <- LaplacesDemon(Data$Model, Data=Data, Initial.Values=SV, Iterations=100, Status=10, Thinning=1, Algorithm="CHARM" )


        plot( f$Summary1[1:155, "Mean"] ~ gstat.pred$log_zinc.pred )


        parm0=Data$PGF(Data)

        if (0) {
          Interval=1.0E-6; Iterations=100;Method="SPG";Samples=1000;CovEst="Hessian"; sir=TRUE;Stop.Tolerance=1.0E-5;CPUs=1;Type="PSOCK"
        }

        f = LaplaceApproximation(Data$Model, Data=Data, parm=parm0, Method="BFGS", Iterations=10, CPUs=4, Stop.Tolerance=1.0E-9 )


        if (plotdata) {

          f = LaplacesDemon(Data$Model, Data=Data, Initial.Values=as.initial.values(f), Fit.LA=f,
            Iterations=10000, Thinning=100, Status=1000, Covar=f$Covar, CPUs=8 )

          parm0 = as.initial.values(f)
          f0 = LaplacesDemon(Data$Model, Data=Data, Initial.Values=parm0, CPUs=8 )
          mu = f$Summary1[,1]
          f0 = LaplacesDemon(Data$Model, Data=Data, Initial.Values=as.initial.values(f), Fit.LA=f,
            Iterations=5000, Thinning=1, Status=1000, Algorithm="IM", Specs=list(mu=mu),
            Covar=f$Covar, CPUs=8 )

          f0 = LaplacesDemon(Data$Model, Data=Data, Initial.Values=as.initial.values(f), Fit.LA=f,
            Iterations=10000, Thinning=100, Status=1000, Covar=f$Covar, CPUs=8 )

          Consort(f0)
          plot(f0, Data=Data)
          m = f0$Summary2[grep( "\\<yhat\\>", rownames( f0$Summary2 ) ),]
          m = f$Summary2[grep( "\\<yhat\\>", rownames( f$Summary2 ) ),]
          # m = f$Summary2[grep( "muSpatial", rownames( f$Summary2 ) ),]
          plot( Data$y ~ m[, "Mean"]  )


        }

        out$LaplacesDemon = list( fit=f, vgm=NA, model=Data$Model, range=NA,
          varSpatial=f$Summary2["sigma2", "Mean"] *out$varZ,
          varObs=f$Summary2["tau2", "Mean"]*out$varZ,
          nu=f$Summary2["nu", "Mean"],
          phi = out$maxdist * ( f$Summary2["phi", "Mean"]  / sqrt(2*f$Summary2["nu", "Mean"] ) )
        )   ## need to check parameterization...

        #out$LaplacesDemon$range = geoR::practicalRange("matern", phi=out$LaplacesDemon$phi, kappa=out$LaplacesDemon$nu)
        out$LaplacesDemon$range = matern_phi2distance(phi=out$LaplacesDemon$phi, nu=out$LaplacesDemon$nu)

       # print( out$LaplacesDemon )


        if (plotdata) {
          plot.new()
          x = seq( 0,  out$LaplacesDemon$range * 1.25, length.out=100 )
          svar =  out$LaplacesDemon$varObs + out$LaplacesDemon$varSpatial * (1-geoR::matern( x, phi=out$LaplacesDemon$phi, kappa=out$LaplacesDemon$nu  ))
          plot( svar~x, type="l", ylim=c(0, max(svar)) )
          abline( h=out$LaplacesDemon$varObs + out$LaplacesDemon$varSpatial )
          abline( h=out$LaplacesDemon$varObs )
          abline( v=out$LaplacesDemon$range, col="red"  )
        }

    ss = summary(Hmodel)
    emei_stats = list( sdTotal=sdTotal, rsquared=ss$r.sq, ndata=ss$n ) # must be same order as p$statsvars

    return( list( predictions=newdata, emei_stats=emei_stats ) )


  testing = FALSE
  if (testing) {
    require(LaplacesDemonCpp)
    require(emei)
    n = 12
    nn=n*n
    T=10
    
    SV = c(rho0 = 0.25, 
        sigma2 = 0.2, zeta = 0.25, rho1 = 0.25, gamma = 1, alpha = 0.1, muX = 0, muY = 0, tau2 = 0.2 )
    
    w = spate::spate.sim(par=SV, n=n, T=T)$xi

    NF=n*n
    Z = spate::spate.init(n=n, T=T, NF=NF)
    ns = Z$ns
    NFc = (NF-Z$ns)/2
    indCosOnly = 1:4
 
    dt = 1
    nu = 1
    w2 = apply(Z$wave^2, 2, sum)
    d = 2

    sdTotal = sd(y)

        Z = spate::spate.init(n=n, T=T, NF=NF)
 
        d = 2

        Data = list(
          eps = 1e-6,
          epsb = 1e-1,
          d=d,
          N = n*n*T,  # required for LaplacesDemon
          NFc = (NF-Z$ns)/2,
          T=as.integer(T), 
          n=as.integer(n), 
          ns=as.integer(Z$ns),
          NF=as.integer(NF), 
          dt=1,
          w2=apply(Z$wave^2, 2, sum),
          M=list(wave=Z$wave),  #put matrix inside a list to hide it from LaplacesDEmon
          pid2 = pi^(d/2),
          nunu2 = 2^(nu - 1) * nu, 
          nu2 = 2 * nu,
          nud2 = nu + d/2,
          yvar= sdTotal^2,
          indCosOnly = as.integer(1:Z$ns),
          indCos=as.integer(Z$indFFT$indCos),
          indW=as.integer(Z$indFFT$indW),
          indWCon=as.integer(Z$indFFT$indWCon),
          y=spate::TSmat.to.vect(w)
        )

        Data$indNA = is.na(Data$y)
        Data$nNA = length( which(Data$indNA) )
        if (Data$nNA > 0 ) Data$y[Data$indNA]=0 
        
        rm(Z)

        Data = emei_LaplacesDemon_spatemodel(Data)

     
        str(Data$Model( parm=Data$PGF(Data), Data ) ) # test to see if return values are sensible

        f = LaplacesDemon(Data$Model, Data=Data, Initial.Values=SV, Iterations=100, Status=10, Thinning=1) 
        parm0 = f$Summary1[Data$parm.names,"Mean"]
        f = LaplaceApproximation(Data$Model, Data=Data, parm=parm0, Method="BFGS", Iterations=100, CPUs=4, Stop.Tolerance=1.0E-9 )
        f <- LaplacesDemon(Data$Model, Data=Data, as.initial.values(f),
             Covar=f$Covar, Iterations=1000, Status=270, Thinning=20,
             Algorithm="AMWG", Specs=list(Periodicity=35))

        Consort(f)
          plot(f, Data=Data)
          m = f$Summary2[grep( "\\<yhat\\>", rownames( f$Summary2 ) ),]
          plot( Data$y ~ m[, "Mean"]  )

        plot( f$Summary1[1:155, "Mean"] ~ gstat.pred$log_zinc.pred )


        parm0=Data$PGF(Data)

        if (0) {
          Interval=1.0E-6; Iterations=100;Method="SPG";Samples=1000;CovEst="Hessian"; sir=TRUE;Stop.Tolerance=1.0E-5;CPUs=1;Type="PSOCK"
        }

        f = LaplaceApproximation(Data$Model, Data=Data, parm=parm0, Method="BFGS", Iterations=10, CPUs=4, Stop.Tolerance=1.0E-9 )


  }

}
