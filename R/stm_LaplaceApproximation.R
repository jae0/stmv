
stm_LaplaceApproximation <- function(Model, parm, Data, Interval=1.0E-6,
     Iterations=100, Method="BFGS", Samples=1000, CovEst="Hessian",
     sir=TRUE, Stop.Tolerance=1.0E-5, CPUs=1, Type="PSOCK") {
     

     # template for doing laplace approximation using LaplacesDemon 
     # ... not used in LBM as we need simultaenous predictions as well but good for parameter estimation   
     
     testing = FALSE
     
     if (testing) {
          Interval=1.0E-6; Iterations=100; Method="BFGS"; Samples=1000; CovEst="Hessian"; sir=TRUE; 
          Stop.Tolerance=1.0E-5; CPUs=1; Type="PSOCK"
         
          require(stm)
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
         parm = Initial.Values =  apply( g$Post[, 501:1000 ], 1, median )
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
  
          Data = stm_LaplacesDemon_spatemodel(Data)
          Data$yhat=Data$y[]*0

          str(Data$Model( parm=Data$PGF(Data), Data ) ) # test to see if return values are sensible

          Model = Data$Model

     }

     # copied from LaplacesDemon::LaplaceApproximation ..simplified with fewer options 
     # and permits simultaneous prediction in the style of spate::spate.mcmc

     ##########################  Initial Checks  ##########################
     time1 <- proc.time()
     LA.call <- match.call()
   
     if(is.null(Data[["mon.names"]]))
          stop("In Data, mon.names is NULL.")
     if(is.null(Data[["parm.names"]]))
          stop("In Data, parm.names is NULL.")
  

     if({Interval <= 0} | {Interval > 1}) Interval <- 1.0E-6
     Iterations <- min(max(round(Iterations), 10), 1000000)
  
     if(Stop.Tolerance <= 0) Stop.Tolerance <- 1.0E-5
  

     if(!is.null(Data[["N"]])) if(length(Data[["N"]]) == 1) N <- Data[["N"]] 

     m.old <- Model(parm, Data)
     if(!is.list(m.old)) stop("Model must return a list.")
     if(length(m.old) != 5) stop("Model must return five components.")
     if(any(names(m.old) != c("LP","Dev","Monitor","yhat","parm")))
          stop("Name mismatch in returned list of Model function.")
     # if(length(m.old[["LP"]]) > 1) stop("Multiple joint posteriors exist!")
     if(!identical(length(parm), length(m.old[["parm"]])))
          stop("The number of initial values and parameters differs.")
     if(!is.finite(m.old[["LP"]])) {
          cat("Generating initial values due to a non-finite posterior.\n")
          if(!is.null(Data[["PGF"]]))
               Initial.Values <- GIV(Model, Data, PGF=TRUE)
          else Initial.Values <- GIV(Model, Data, PGF=FALSE)
          m.old <- Model(Initial.Values, Data)
          }
     if(!is.finite(m.old[["LP"]])) stop("The posterior is non-finite.")
     if(!is.finite(m.old[["Dev"]])) stop("The deviance is non-finite.")
     parm <- m.old[["parm"]]
     if(!identical(Model(m.old[["parm"]], Data)[["LP"]], m.old[["LP"]])) {
          cat("WARNING: LP differs when initial values are held constant.\n")
          cat("     Derivatives may be problematic if used.\n")}

     ####################  Begin Laplace Approximation  ###################
     cat("Laplace Approximation using BFGS begins...\n")


     m.new <- m.old
     Dev <- matrix(m.old[["Dev"]],1,1)
     parm.old <- parm
     parm.len <- length(as.vector(parm))
     post <- matrix(m.old[["parm"]], 1, parm.len)
     tol.new <- 1
     keepgoing <- TRUE
     g.old <- g.new <- rep(0, parm.len)
     B <- diag(parm.len) #Approximate Hessian
     options(warn=-1)
     
     ## do a burn-in here ...

     for (iter in 2:Iterations) {

          # --------------- addition ----------------
          Data$y[Data$indNA] = m.old$yhat[Data$indNA] + RcppZiggurat::zrnorm(Data$nNA) * sqrt(m.old$parm["tau2"]) # 4-5 X faster than rnorm
          # -----------------------------------------

          ### Print Status
          if(iter %% round(Iterations / 1) == 0)
               cat("Iteration: ", iter, " of ", Iterations, ",   LP: ", round(m.old[["LP"]],1), "\n")
          
          ### Gradient and Direction p
          g.old <- g.new
          g.new <- -1*partial(Model, m.old[["parm"]], Data, Interval)
          p <- as.vector(tcrossprod(g.new, -B))
          p[which(!is.finite(p))] <- 0
          
          ### Step-size Line Search
          Step.Size <- 0.8
          changed <- TRUE
          while(m.new[["LP"]] <= m.old[["LP"]] & changed == TRUE) {
               Step.Size <- Step.Size * 0.2
               s <- Step.Size*p
               prop <- m.old[["parm"]] + s
               changed <- !identical(m.old[["parm"]], prop)
               m.new <- Model(prop, Data)
               if(any(!is.finite(c(m.new[["LP"]], m.new[["Dev"]],
                    m.new[["Monitor"]]))))
                    m.new <- m.old
          }
          
          ### BFGS Update to Approximate Hessian B
          if(m.new[["LP"]] > m.old[["LP"]]) {
               m.old <- m.new
               g <- g.new - g.old
               CC <- sum(s*g) #Curvature condition
               if(CC > 0) {
                    y <- as.vector(crossprod(B, g))
                    DD <- as.double(1 + crossprod(g, y)/CC)
                    B <- B - (tcrossprod(s, y) + tcrossprod(y, s) - DD * tcrossprod(s, s))/CC}
               if(any(!is.finite(B))) B <- diag(parm.len)
          }

          ### Storage
          post <- rbind(post, m.old[["parm"]])
          Dev <- rbind(Dev, m.old[["Dev"]])
          
          ### Tolerance
          tol.new <- sqrt(sum(s*s))
          if(keepgoing == FALSE) tol.new <- 0
          if(tol.new <= Stop.Tolerance) break
     }


     options(warn=0)
     ### Output
     LA <- list(Dev=Dev, iter=iter, parm.len=parm.len,
          parm.new=m.old[["parm"]], parm.old=parm.old, post=post,
          Step.Size=Step.Size, tol.new=tol.new)


     Dev <- as.vector(LA$Dev)
     
     if(is.null(LA$H)) {
          H <- FALSE
     } else {
          H <- LA$H
     }

     iter <- LA$iter
     parm.len <- LA$parm.len
     parm.new <- LA$parm.new
     parm.old <- LA$parm.old
     post <- LA$post
     Step.Size <- LA$Step.Size
     tol.new <- LA$tol.new
     rm(LA)
     if(iter == 1) stop("LaplaceApproximation stopped at iteration 1.")

     if(tol.new <= Stop.Tolerance) {
          converged <- TRUE
     } else {
          converged <- FALSE
     }
     

     ### Column names to samples
     if(ncol(post) == length(Data[["parm.names"]])) colnames(post) <- Data[["parm.names"]]
     rownames(post) <- 1:nrow(post)

     ########################  Covariance Matirx  #########################
     cat("Estimating the Covariance Matrix\n")
     if(all(H == FALSE)) {
          VarCov <- CovEstim(Model, parm.new, Data, Method=CovEst)
     } else {
          VarCov <- -as.inverse(as.symmetric.matrix(H))
          diag(VarCov) <- abs(diag(VarCov))
     }

     #################  Sampling Importance Resampling  ##################
     if({sir == TRUE} & {converged == TRUE}) {
          cat("Sampling from Posterior with Sampling Importance Resampling\n")
          posterior <- SIR(Model, Data, mu=parm.new, Sigma=VarCov, n=Samples, CPUs=CPUs, Type=Type)
          Mon <- matrix(0, nrow(posterior), length(Data[["mon.names"]]))
          dev <- rep(0, nrow(posterior))
          for (i in 1:nrow(posterior)) {
               mod <- Model(posterior[i,], Data)
               dev[i] <- mod[["Dev"]]
               Mon[i,] <- mod[["Monitor"]]
          }
          colnames(Mon) <- Data[["mon.names"]]
     } else {
          if({sir == TRUE} & {converged == FALSE})
               cat("Posterior samples are not drawn due to Converge=FALSE\n")
          posterior <- NA; Mon <- NA
     }

     #####################  Summary, Point-Estimate  ######################
     cat("Creating Summary from Point-Estimates\n")
     Summ1 <- matrix(NA, parm.len, 4, dimnames=list(Data[["parm.names"]], c("Mode","SD","LB","UB")))
     Summ1[,1] <- parm.new
     Summ1[,2] <- sqrt(diag(VarCov))
     Summ1[,3] <- parm.new - 2*Summ1[,2]
     Summ1[,4] <- parm.new + 2*Summ1[,2]
     
     ###################  Summary, Posterior Samples  ####################
     Summ2 <- NA
     if({sir == TRUE} & {converged == TRUE}) {
          cat("Creating Summary from Posterior Samples\n")
          Summ2 <- matrix(NA, ncol(posterior), 7,
               dimnames=list(Data[["parm.names"]],
                    c("Mode","SD","MCSE","ESS","LB","Median","UB")))
          Summ2[,1] <- colMeans(posterior)
          Summ2[,2] <- sqrt(.colVars(posterior))
          Summ2[,3] <- Summ2[,2] / sqrt(nrow(posterior))
          Summ2[,4] <- rep(nrow(posterior), ncol(posterior))
          Summ2[,5] <- apply(posterior, 2, quantile, c(0.025))
          Summ2[,6] <- apply(posterior, 2, quantile, c(0.500))
          Summ2[,7] <- apply(posterior, 2, quantile, c(0.975))
          Deviance <- rep(0, 7)
          Deviance[1] <- mean(dev)
          Deviance[2] <- sd(dev)
          Deviance[3] <- sd(dev) / sqrt(nrow(posterior))
          Deviance[4] <- nrow(posterior)
          Deviance[5] <- as.numeric(quantile(dev, probs=0.025, na.rm=TRUE))
          Deviance[6] <- as.numeric(quantile(dev, probs=0.500, na.rm=TRUE))
          Deviance[7] <- as.numeric(quantile(dev, probs=0.975, na.rm=TRUE))
          Summ2 <- rbind(Summ2, Deviance)
          for (j in 1:ncol(Mon)) {
               Monitor <- rep(NA,7)
               Monitor[1] <- mean(Mon[,j])
               Monitor[2] <- sd(as.vector(Mon[,j]))
               Monitor[3] <- sd(as.vector(Mon[,j])) / sqrt(nrow(Mon))
               Monitor[4] <- nrow(Mon)
               Monitor[5] <- as.numeric(quantile(Mon[,j], probs=0.025,
                    na.rm=TRUE))
               Monitor[6] <- as.numeric(quantile(Mon[,j], probs=0.500,
                    na.rm=TRUE))
               Monitor[7] <- as.numeric(quantile(Mon[,j], probs=0.975,
                    na.rm=TRUE))
               Summ2 <- rbind(Summ2, Monitor)
               rownames(Summ2)[nrow(Summ2)] <- Data[["mon.names"]][j]
          }
     }

     ###############  Logarithm of the Marginal Likelihood  ###############
     LML <- list(LML=NA, VarCov=VarCov)
     if({sir == TRUE} & {converged == TRUE}) {
          cat("Estimating Log of the Marginal Likelihood\n")
          lml <- LML(theta=posterior, LL=(dev*(-1/2)), method="NSIS")
          LML[[1]] <- lml[[1]]
     } else if ({sir == FALSE} & {converged == TRUE}) {
          cat("Estimating Log of the Marginal Likelihood\n")
          LML <- LML(Model, Data, Modes=parm.new, Covar=VarCov, method="LME")
     }
     colnames(VarCov) <- rownames(VarCov) <- Data[["parm.names"]]
     time2 <- proc.time()

     #############################  Output  ##############################
     LA <- list(Call=LA.call,
          Converged=converged,
          Covar=VarCov,
          Deviance=as.vector(Dev),
          History=post,
          Initial.Values=parm,
          Iterations=iter,
          LML=LML[[1]],
          LP.Final=as.vector(Model(parm.new, Data)[["LP"]]),
          LP.Initial=m.old[["LP"]],
          Minutes=round(as.vector(time2[3] - time1[3]) / 60, 2),
          Monitor=Mon,
          Posterior=posterior,
          Step.Size.Final=Step.Size,
          Step.Size.Initial=1,
          Summary1=Summ1,
          Summary2=Summ2,
          Tolerance.Final=tol.new,
          Tolerance.Stop=Stop.Tolerance)
     class(LA) <- "laplace"
     cat("Laplace Approximation is finished.\n\n")
     return(LA)
}

