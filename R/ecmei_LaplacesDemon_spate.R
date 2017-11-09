
ecmei_LaplacesDemon_spate = function( Model, Data, Initial.Values, Covar=NULL, Iterations=1000, 
     Status=100, Thinning=10, Algorithm="MWG", Specs=NULL, LogFile="", Fit.LA=NULL,  
     Debug=list(DB.chol=FALSE, DB.eigen=FALSE, DB.MCSE=FALSE, DB.Model=FALSE), ... ) {

     # copied from LaplacesDemon, modified for speed with a focus upon spatial/temporal modelling

     testing = FALSE
     
     if (testing) {

          Covar=NULL; Iterations=1000 ;Status=1 ;Thinning=2; Algorithm="MWG"; Specs=NULL; LogFile=""; Fit.LA=NULL
          Debug=list(DB.chol=FALSE, DB.eigen=FALSE, DB.MCSE=FALSE, DB.Model=FALSE)

          Specs = NULL

          if (0) {
              Algorithm="MWG"
              Specs= NULL

              Algorithm="NUTS"
              Specs = list( A=10, delta=1/10, epsilon=1e-6, Lmax=10)

              Algorithm="CHARM"
              Specs= list( alpha.star=0.5)
          }

          require(LaplacesDemonCpp)
          require(ecmei)
          n = 12
          nn=n*n
          T=20

          SV = c(rho0 = 0.25, 
          sigma2 = 0.2, zeta = 0.25, rho1 = 0.25, gamma = 1, alpha = 0.1, muX = 0, muY = 0, tau2 = 0.2 )

          w  = spate::spate.sim(par=SV, n=n, T=T)$xi
          
          toremove = sample.int( nn, size=floor(0.75*nn) ) 
          w[,toremove] =  NA 

          y = w
          sdTotal = sd(y, na.rm=TRUE)
          
          indNA = is.na(y)
          nNA = length( which(indNA) )
          if (nNA > 0 ) y[indNA]=0 

          NF=n*n
          Z = spate::spate.init(n=n, T=T, NF=NF)
          ns = Z$ns
          NFc = (NF-Z$ns)/2
          indCosOnly = 1:4

          dt = 1
          nu = 1
          w2 = apply(Z$wave^2, 2, sum)
          d = 2


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
  
          Data = ecmei_LaplacesDemon_spatemodel(Data)
          Data$yhat=Data$y[]*0

          str(Data$Model( parm=Data$PGF(Data), Data ) ) # test to see if return values are sensible

          Model = Data$Model

     }

     # cat("\nLaplace's Demon was called on ", date(), "\n", sep="",     #      file=LogFile, append=TRUE)
     time1 <- proc.time()
     LDcall <- match.call()

     Iterations <- round(abs(Iterations))

     if(Iterations < 11) {
          Iterations <- 11
          cat("'Iterations' has been changed to ", Iterations, ".\n", sep="", file=LogFile, append=TRUE)}
     
     Status <- round(abs(Status))
     if({Status < 1} || {Status > Iterations}) {
          Status <- Iterations
          cat("'Status' has been changed to ", Status, ".\n", sep="", file=LogFile, append=TRUE)}
     
     Thinning <- round(abs(Thinning))
     if({Thinning < 1} || {Thinning > Iterations}) {
          Thinning <- 1
          cat("'Thinning' has been changed to ", Thinning, ".\n", sep="", file=LogFile, append=TRUE)}


    # Parameters specific to each algorithm
    LIV <- length(Initial.Values)
    ScaleF <- 2.381204 * 2.381204 / LIV

    LDspecs = ecmei_LaplacesDemon_Specification( Algorithm=Algorithm, LogFile=LogFile, Specs=Specs, LIV=LIV, Covar=Covar, ScaleF=ScaleF )
    


     Mo0 <- Model(Initial.Values, Data)
     if(!is.list(Mo0))
          stop("Model must return a list.", file=LogFile, append=TRUE)
     if(length(Mo0) != 5)
          stop("Model must return five components.", file=LogFile, append=TRUE)
     if(any(names(Mo0) != c("LP","Dev","Monitor","yhat","parm")))
          stop("Name mismatch in returned list of Model function.", file=LogFile, append=TRUE)
     if(length(Mo0[["LP"]]) > 1)
          stop("Multiple joint posteriors exist!", file=LogFile, append=TRUE)
     if(!identical(length(Mo0[["Monitor"]]), length(Data[["mon.names"]])))
          stop("Length of mon.names differs from length of monitors.", file=LogFile, append=TRUE)
     
     
     if(!identical(Model(Mo0[["parm"]], Data)[["LP"]], Mo0[["LP"]])) {
          cat("WARNING: LP differs when initial values are held constant.\n", file=LogFile, append=TRUE)
          cat("     Derivatives may be problematic if used.\n", file=LogFile, append=TRUE)}
     

     #########################  Initial Settings  #########################
     Acceptance <- 0
     if(!is.finite(Mo0[["LP"]])) {
          cat("Generating initial values due to a non-finite posterior.\n", file=LogFile, append=TRUE)
          if(!is.null(Data[["PGF"]]))
               Initial.Values <- GIV(Model, Data, PGF=TRUE)
          else Initial.Values <- GIV(Model, Data)
          Mo0 <- Model(Initial.Values, Data)
     }
     if(is.infinite(Mo0[["LP"]]))
          stop("The posterior is infinite!", file=LogFile, append=TRUE)
     if(is.nan(Mo0[["LP"]]))
          stop("The posterior is not a number!", file=LogFile, append=TRUE)
     if(is.na(Mo0[["Dev"]]))
          stop("The deviance is a missing value!", file=LogFile, append=TRUE)
     if(is.infinite(Mo0[["Dev"]]))
          stop("The deviance is infinite!", file=LogFile, append=TRUE)
     if(is.nan(Mo0[["Dev"]]))
          stop("The deviance is not a number!", file=LogFile, append=TRUE)
     if(any(is.na(Mo0[["Monitor"]])))
          stop("Monitored variable(s) have a missing value!", file=LogFile, append=TRUE)
     if(any(is.infinite(Mo0[["Monitor"]])))
          stop("Monitored variable(s) have an infinite value!", file=LogFile, append=TRUE)
     if(any(is.nan(Mo0[["Monitor"]])))
          stop("Monitored variable(s) include a value that is not a number!", file=LogFile, append=TRUE)
     if(Algorithm == "t-walk") {
          Mo0 <- Model(Initial.Values, Data)
          if(any(Mo0[["parm"]] == LDspecs$Specs[["SIV"]]))
              stop("Initial.Values and SIV not unique after model update.", file=LogFile, append=TRUE)
      }
    

     ######################  Laplace Approximation  #######################
     ### Sample Size of Data

     if(!is.null(Data[["N"]])) if(length(Data[["N"]]) == 1) N <- Data[["N"]]

     if(is.null(N))
          stop("Sample size of Data not found in n, N, y, or Y.", file=LogFile, append=TRUE)
     if({all(Initial.Values == 0)} & {N >= 5*length(Initial.Values)}) {
        if (isnull(Fit.LA)) {
          cat("\nLaplace Approximation will be used on initial values.\n", file=LogFile, append=TRUE)
          Fit.LA <- LaplaceApproximation(Model, Initial.Values, Data, Method="SPG", CovEst="Identity", sir=FALSE)
          Covar <- 2.381204 * 2.381204 / length(Initial.Values) * Fit.LA$Covar
          Initial.Values <- Fit.LA$Summary1[1:length(Initial.Values),1]
          cat("The covariance matrix from Laplace Approximation has been scaled\n", file=LogFile, append=TRUE)
          cat("for Laplace's Demon, and the posterior modes are now the initial\n", file=LogFile, append=TRUE)
          cat("values for Laplace's Demon.\n\n", file=LogFile, append=TRUE)
        } else { 
          Covar = Fit.LA$Covar
          Initial.Values = as.initial.values(Fit.LA)
        }
     }



     #########################  Prepare for MCMC  #########################
     Mo0 <- Model(Initial.Values, Data)
     Dev <- matrix(Mo0[["Dev"]], floor(Iterations/Thinning)+1, 1)
     Mon <- matrix(Mo0[["Monitor"]], floor(Iterations/Thinning)+1, length(Mo0[["Monitor"]]), byrow=TRUE)
     thinned <- matrix(Initial.Values, floor(Iterations/Thinning)+1, length(Initial.Values), byrow=TRUE)
     


     ############################  Begin MCMC  ############################

     cat("Algorithm:", LDspecs$alg, "\n", file=LogFile, append= TRUE)
     cat("\nLaplace's Demon is beginning to update...\n", file=LogFile, append=TRUE)
     options(warn=2)
     on.exit(options(warn=0))

     # Status=Status, Thinning=Thinning,
     # Acceptance=Acceptance, Dev=Dev, 
     # LIV=LIV, Mon=Mon, Mo0=Mo0, ScaleF=ScaleF, thinned=thinned,
     # LogFile=LogFile,
     
     Specs=LDspecs$Specs
     DiagCovar=LDspecs$DiagCovar
     tuning=LDspecs$tuning
     VarCov=LDspecs$VarCov 
     parm.names=Data[["parm.names"]]
     Debug=list(DB.chol=FALSE, DB.eigen=FALSE, DB.MCSE=FALSE, DB.Model=FALSE) 
     
     mcmc.out = LDspecs$mcmcFN() 

      if (0) {
      
        microbenchmark::microbenchmark( {

        }, times=3 )
 

     }

     options(warn=0)


     #########################  MCMC is Finished  #########################
     Acceptance <- mcmc.out$Acceptance
     Dev <- mcmc.out$Dev
     DiagCovar <- mcmc.out$DiagCovar
     Mon <- mcmc.out$Mon
     thinned <- mcmc.out$thinned
     VarCov <- mcmc.out$VarCov
     remove(mcmc.out)
     rownames(DiagCovar) <- NULL
     colnames(DiagCovar) <- Data[["parm.names"]]
     thinned <- matrix(thinned[-1,], nrow(thinned)-1, ncol(thinned))
     Dev <- matrix(Dev[-1,], nrow(Dev)-1, 1)
     Mon <- matrix(Mon[-1,], nrow(Mon)-1, ncol(Mon))
     if(is.matrix(VarCov) & !is.list(VarCov)) {
          colnames(VarCov) <- rownames(VarCov) <- Data[["parm.names"]]
     } else if(is.vector(VarCov) & !is.list(VarCov)) {
          names(VarCov) <- Data[["parm.names"]]
     }
     thinned.rows <- nrow(thinned)
     
     ### Warnings (After Updating)
     if(any(Acceptance == 0))
          cat("\nWARNING: All proposals were rejected.\n", file=LogFile, append=TRUE)
     
     ### Real Values
     thinned[which(!is.finite(thinned))] <- 0
     Dev[which(!is.finite(Dev))] <- 0
     Mon[which(!is.finite(Mon))] <- 0
     
     ### Assess Stationarity
     cat("\nAssessing Stationarity\n", file=LogFile, append=TRUE)
     if(thinned.rows %% 10 == 0) thinned2 <- thinned
     if(thinned.rows %% 10 != 0) thinned2 <- thinned[1:(10*trunc(thinned.rows/10)),]
     HD <- BMK.Diagnostic(thinned2, batches=10)
     Ind <- 1 * (HD > 0.5)
     BurnIn <- thinned.rows
     batch.list <- seq(from=1, to=nrow(thinned2), by=floor(nrow(thinned2)/10))
     for (i in 1:9) {
          if(sum(Ind[,i:9]) == 0) {
               BurnIn <- batch.list[i] - 1
               break
          }
      }
     Stat.at <- BurnIn + 1
     rm(batch.list, HD, Ind, thinned2)
     
     ### Assess Thinning and ESS Size for all parameter samples
     cat("Assessing Thinning and ESS\n", file=LogFile, append=TRUE)
     acf.rows <- trunc(10*log10(thinned.rows))
     acf.temp <- matrix(1, acf.rows, LIV)
     ESS1 <- Rec.Thin <- rep(1, LIV)
     for (j in 1:LIV) {
          temp0 <- acf(thinned[,j], lag.max=acf.rows, plot=FALSE)
          if(length(temp0$acf[-1,1,1]) == acf.rows)
               acf.temp[,j] <- abs(temp0$acf[-1,1,1])
          ESS1[j] <- ESS(thinned[,j])
          Rec.Thin[j] <- which(acf.temp[,j] <= 0.1)[1]*Thinning
     }
     
     Rec.Thin[which(is.na(Rec.Thin))] <- nrow(acf.temp)
     
     ### Assess ESS for all deviance and monitor samples
     ESS2 <- ESS(Dev)
     ESS3 <- ESS(Mon)
     
     ### Assess ESS for stationary samples
     if(Stat.at < thinned.rows) {
          ESS4 <- ESS(thinned[Stat.at:thinned.rows,])
          ESS5 <- ESS(Dev[Stat.at:thinned.rows])
          ESS6 <- ESS(Mon[Stat.at:thinned.rows,])
      }


     ### Posterior Summary Table 1: All Thinned Samples
     cat("Creating Summaries\n", file=LogFile, append=TRUE)
     Num.Mon <- ncol(Mon)
     Summ1 <- matrix(NA, LIV, 7, dimnames=list(Data[["parm.names"]], c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     Summ1[,1] <- colMeans(thinned)
     Summ1[,2] <- sqrt(.colVars(thinned))
     Summ1[,3] <- 0
     Summ1[,4] <- ESS1
     Summ1[,5] <- apply(thinned, 2, quantile, c(0.025), na.rm=TRUE)
     Summ1[,6] <- apply(thinned, 2, quantile, c(0.500), na.rm=TRUE)
     Summ1[,7] <- apply(thinned, 2, quantile, c(0.975), na.rm=TRUE)
     for (i in 1:ncol(thinned)) {
          temp <- try(MCSE(thinned[,i]), silent=TRUE)
          if(!inherits(temp, "try-error")) Summ1[i,3] <- temp
          else Summ1[i,3] <- MCSE(thinned[,i], method="sample.variance")
      }
     Deviance <- rep(NA,7)
     Deviance[1] <- mean(Dev)
     Deviance[2] <- sd(as.vector(Dev))
     temp <- try(MCSE(as.vector(Dev)), silent=TRUE)
     if(inherits(temp, "try-error")) temp <- MCSE(as.vector(Dev), method="sample.variance")
     Deviance[3] <- temp
     Deviance[4] <- ESS2
     Deviance[5] <- as.numeric(quantile(Dev, probs=0.025, na.rm=TRUE))
     Deviance[6] <- as.numeric(quantile(Dev, probs=0.500, na.rm=TRUE))
     Deviance[7] <- as.numeric(quantile(Dev, probs=0.975, na.rm=TRUE))
     Summ1 <- rbind(Summ1, Deviance)

     for (j in 1:Num.Mon) {
          Monitor <- rep(NA,7)
          Monitor[1] <- mean(Mon[,j])
          Monitor[2] <- sd(as.vector(Mon[,j]))
          temp <- try(MCSE(as.vector(Mon[,j])), silent=TRUE)
          if(inherits(temp, "try-error")) temp <- MCSE(Mon[,j], method="sample.variance")
          Monitor[3] <- temp
          Monitor[4] <- ESS3[j]
          Monitor[5] <- as.numeric(quantile(Mon[,j], probs=0.025, na.rm=TRUE))
          Monitor[6] <- as.numeric(quantile(Mon[,j], probs=0.500, na.rm=TRUE))
          Monitor[7] <- as.numeric(quantile(Mon[,j], probs=0.975, na.rm=TRUE))
          Summ1 <- rbind(Summ1, Monitor)
          rownames(Summ1)[nrow(Summ1)] <- Data[["mon.names"]][j]
     }
     
     ### Posterior Summary Table 2: Stationary Samples
     Summ2 <- matrix(NA, LIV, 7, dimnames=list(Data[["parm.names"]], c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     if(Stat.at < thinned.rows) {
          thinned2 <- matrix(thinned[Stat.at:thinned.rows,], thinned.rows-Stat.at+1, ncol(thinned))
          Dev2 <- matrix(Dev[Stat.at:thinned.rows,], thinned.rows-Stat.at+1, ncol(Dev))
          Mon2 <- matrix(Mon[Stat.at:thinned.rows,], thinned.rows-Stat.at+1, ncol(Mon))
          Summ2[,1] <- colMeans(thinned2)
          Summ2[,2] <- sqrt(.colVars(thinned2))
          Summ2[,3] <- 0
          Summ2[,4] <- ESS4
          Summ2[,5] <- apply(thinned2, 2, quantile, c(0.025), na.rm=TRUE)
          Summ2[,6] <- apply(thinned2, 2, quantile, c(0.500), na.rm=TRUE)
          Summ2[,7] <- apply(thinned2, 2, quantile, c(0.975), na.rm=TRUE)
          for (i in 1:ncol(thinned2)) {
               temp <- try(MCSE(thinned2[,i]), silent=TRUE)
               if(!inherits(temp, "try-error")) Summ2[i,3] <- temp
               else Summ2[i,3] <- MCSE(thinned2[,i], method="sample.variance")
          }
          Deviance <- rep(NA,7)
          Deviance[1] <- mean(Dev2)
          Deviance[2] <- sd(as.vector(Dev2))
          temp <- try(MCSE(as.vector(Dev2)), silent=TRUE)
          if(inherits(temp, "try-error"))
               temp <- MCSE(as.vector(Dev2), method="sample.variance")
          Deviance[3] <- temp
          Deviance[4] <- ESS5
          Deviance[5] <- as.numeric(quantile(Dev2, probs=0.025, na.rm=TRUE))
          Deviance[6] <- as.numeric(quantile(Dev2, probs=0.500, na.rm=TRUE))
          Deviance[7] <- as.numeric(quantile(Dev2, probs=0.975, na.rm=TRUE))
          Summ2 <- rbind(Summ2, Deviance)
          for (j in 1:Num.Mon) {
               Monitor <- rep(NA,7)
               Monitor[1] <- mean(Mon2[,j])
               Monitor[2] <- sd(as.vector(Mon2[,j]))
               temp <- try(MCSE(as.vector(Mon[,j])), silent=TRUE)
               if(inherits(temp, "try-error")) temp <- MCSE(as.vector(Mon[,j]), method="sample.variance")
               Monitor[3] <- temp
               Monitor[4] <- ESS6[j]
               Monitor[5] <- as.numeric(quantile(Mon2[,j], probs=0.025, na.rm=TRUE))
               Monitor[6] <- as.numeric(quantile(Mon2[,j], probs=0.500, na.rm=TRUE))
               Monitor[7] <- as.numeric(quantile(Mon2[,j], probs=0.975, na.rm=TRUE))
               Summ2 <- rbind(Summ2, Monitor)
               rownames(Summ2)[nrow(Summ2)] <- Data[["mon.names"]][j]
          }
    }

     ### Column names to samples
     if(identical(ncol(Mon), length(Data[["mon.names"]])))
          colnames(Mon) <- Data[["mon.names"]]
     if(identical(ncol(thinned), length(Data[["parm.names"]]))) {
          colnames(thinned) <- Data[["parm.names"]]}
     
     ### Logarithm of the Marginal Likelihood
     LML <- list(LML=NA, VarCov=NA)
     
     if( LDspecs$LogMarginalLikelihood  & {Stat.at < thinned.rows}) {
          cat("Estimating Log of the Marginal Likelihood\n", file=LogFile, append=TRUE)
          LML <- LML(theta=thinned2, LL=as.vector(Dev2)*(-1/2), method="NSIS")
     }
     time2 <- proc.time()
     
     ### Compile Output
     cat("Creating Output\n", file=LogFile, append=TRUE)
     LaplacesDemon.out <- list(
        Acceptance.Rate=round(Acceptance/Iterations,7),
        Algorithm_short=Algorithm,
        Algorithm=LDspecs$alg,
        Call=LDcall,
        Covar=VarCov,
        CovarDHis=DiagCovar,
        Deviance=as.vector(Dev),
        DIC1=c(mean(as.vector(Dev)), var(as.vector(Dev))/2, mean(as.vector(Dev)) + var(as.vector(Dev))/2), 
        DIC2=if (Stat.at < thinned.rows) {
                c(mean(as.vector(Dev2)), var(as.vector(Dev2))/2, mean(as.vector(Dev2)) + var(as.vector(Dev2))/2)
             } else {
                rep(NA,3)
             }, 
        Initial.Values=Initial.Values, 
        Iterations=Iterations, 
        LML=LML[[1]], 
        Minutes=round(as.vector(time2[3] - time1[3]) / 60,2), 
        Model=Model, 
        Monitor=Mon,
        Parameters=LIV,
        Posterior1=thinned,
        Posterior2=if(Stat.at < thinned.rows) {
            thinned[Stat.at:thinned.rows,]
          } else {
            thinned[thinned.rows,]
          },
        Rec.BurnIn.Thinned=BurnIn,
        Rec.BurnIn.UnThinned=BurnIn*Thinning,
        Rec.Thinning=min(1000, max(Rec.Thin)),
        Specs=LDspecs$Specs,
        Status=Status,
        Summary1=Summ1,
        Summary2=Summ2,
        Thinned.Samples=thinned.rows,
        Thinning=Thinning 
     )
     class(LaplacesDemon.out) <- "demonoid"
     cat("\nLaplace's Demon has finished.\n", file=LogFile, append=TRUE)

     return(LaplacesDemon.out)

}


