

# LaplacesDemon core functions copied directly from the tail end of 
# LaplacesDemonCpp::LaplacesDemon 12 July 2017, JC 

# these are here for reference only at the moment

if (0) {


.mcmcadmg <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov, LogFile)
     {
     n <- Specs[["n"]]
     Periodicity <- Specs[["Periodicity"]]
     Acceptance <- matrix(0, 1, LIV)
     k <- 100
     AR <- matrix(0, k, LIV)
     obs.sum <- matrix(Mo0[["parm"]]*n, LIV, 1)
     obs.scatter <- tcrossprod(Mo0[["parm"]])*n
     s <- svd(VarCov)
     U <- diag(s$u)
     tol <- LIV*max(s$d)*.Machine$double.eps
     problem <- any(s$d <= tol)
     DiagCovar <- matrix(diag(VarCov), floor(Iterations/Thinning)+1, LIV, byrow=TRUE)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Random-Scan Componentwise Estimation
          AccRate <- colMeans(AR[(k-min(k,iter)+1):k,,drop=FALSE])
          if(problem == FALSE) lambda <- U*rnorm(LIV, 0, sqrt(0.01 + s$d*exp(2*s$d*(AccRate - 0.3))))
          else lambda <- rnorm(LIV, 0, sqrt(diag(VarCov)))
          AR <- rbind(AR[-1,], 0)
          for (j in sample.int(LIV)) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- prop[j] + lambda[j]
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) {
                    Mo0 <- Mo1
                    Acceptance[j] <- Acceptance[j] + 1}
                    AR[k,j] <- 1}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]
               DiagCovar[t.iter,] <- diag(VarCov)}
          ### Update Sample and Scatter Sum
          obs.sum <- obs.sum + Mo0[["parm"]]
          obs.scatter <- obs.scatter + tcrossprod(Mo0[["parm"]])
          ### Adaptation
          if(iter %% Periodicity == 0) {
               VarCov <- obs.scatter/{n + iter} -
                    tcrossprod(obs.sum/{n + iter})
               diag(VarCov) <- diag(VarCov) + 1e-05
               s <- svd(VarCov)
               U <- diag(s$u)
               tol <- LIV*max(s$d)*.Machine$double.eps
               problem <- any(s$d <= tol)}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)), Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=VarCov)
     return(out)
     }





.mcmcagg <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning, LogFile)
     {
     Grid <- Specs[["Grid"]]
     dparm <- Specs[["dparm"]]
     smax <- Specs[["smax"]]
     CPUs <- Specs[["CPUs"]]
     Packages <- Specs[["Packages"]]
     Dyn.libs <- Specs[["Dyn.libs"]]
     AGGCP <- function(Model, Data, j, Mo0, Grid, tuning, smax)
          {
          G <- length(Grid[[j]])
          x <- Grid[[j]] * sqrt(2) * tuning[j]
          LP.grid <- rep(0, G)
          prop <- Mo0[["parm"]]
          theta <- prop[j] + x
          for (g in 1:G) {
               prop[j] <- theta[g]
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0
               LP.grid[g] <- Mo1[["LP"]]
               theta[g] <- Mo1[["parm"]][j]}
          if(all(!is.finite(LP.grid))) LP.grid <- rep(0, G)
          LP.grid[which(!is.finite(LP.grid))] <- min(LP.grid[which(is.finite(LP.grid))])
          LP.grid <- exp(LP.grid - logadd(LP.grid))
          LP.grid <- LP.grid / sum(LP.grid)
          s <- spline(theta, LP.grid, n=1000)
          s$y <- interval(s$y, 0, Inf, reflect=FALSE)
          if(length(which(s$y > 0)) == 0)
               prop[j] <- theta[which.max(LP.grid)[1]]
          else prop[j] <- sample(s$x, 1, prob=s$y)
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]])))) Mo1 <- Mo0
          else tuning[j] <- min(max(sqrt(sum(LP.grid * x^2)), 1e-10), smax)
          Mo0 <- Mo1
          return(list(Mo0=Mo0, tuning=tuning))
          }
     AGGCPP <- function(Model, Data, j, Mo0, Grid, tuning, smax, cl)
          {
          G <- length(Grid[[j]])
          x <- Grid[[j]] * sqrt(2) * tuning[j]
          LP.grid <- rep(0, G)
          LIV <- length(Mo0[["parm"]])
          prop <- matrix(Mo0[["parm"]], G, LIV, byrow=TRUE)
          prop[, j] <- prop[, j] + x
          Mo1 <- parLapply(cl, 1:G, function(x) Model(prop[x,], Data))
          LP.grid <- as.vector(unlist(lapply(Mo1, function(x) x[["LP"]])))
          prop <- matrix(as.vector(unlist(lapply(Mo1, function(x) x[["parm"]]))), G, LIV, byrow=TRUE)
          theta <- prop[, j]
          if(all(!is.finite(LP.grid))) LP.grid <- rep(0, G)
          LP.grid[which(!is.finite(LP.grid))] <- min(LP.grid[which(is.finite(LP.grid))])
          LP.grid <- exp(LP.grid - logadd(LP.grid))
          LP.grid <- LP.grid / sum(LP.grid)
          s <- spline(theta, LP.grid, n=1000)
          s$y <- interval(s$y, 0, Inf, reflect=FALSE)
          prop <- Mo0[["parm"]]
          if(length(which(s$y > 0)) == 0)
               prop[j] <- theta[which.max(LP.grid)[1]]
          else prop[j] <- sample(s$x, 1, prob=s$y)
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]])))) Mo1 <- Mo0
          else tuning[j] <- min(max(sqrt(sum(LP.grid * x^2)), 1e-10), smax)
          Mo0 <- Mo1
          return(list(Mo0=Mo0, tuning=tuning))
          }
     Acceptance <- matrix(0, 1, LIV)
     Grid.orig <- Grid
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     if(CPUs == 1) {
          for (iter in 1:Iterations) {
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               for (j in sample.int(LIV)) {
                    if(j %in% dparm) Mo0 <- .mcmcggdp(Model, Data, j, Mo0, Grid)
                    else {
                         agg <- AGGCP(Model, Data, j, Mo0, Grid, tuning, smax)
                         Mo0 <- agg$Mo0
                         tuning[j] <- agg$tuning[j]}}
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter/Thinning) + 1
                    thinned[t.iter, ] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter, ] <- Mo0[["Monitor"]]
                    DiagCovar <- rbind(DiagCovar, tuning)}}
          }
     else {
          detectedCores <- detectCores()
          cat("\n\nCPUs Detected:", detectedCores, "\n", file=LogFile, append=TRUE)
          if(CPUs > detectedCores) {
               cat("\nOnly", detectedCores, "will be used.\n", file=LogFile, append=TRUE)
               CPUs <- detectedCores}
          cat("\nLaplace's Demon is preparing environments for CPUs...", file=LogFile, append=TRUE)
          cat("\n##################################################\n", file=LogFile, append=TRUE)
          cl <- makeCluster(CPUs)
          cat("\n##################################################\n", file=LogFile, append=TRUE)
          on.exit(stopCluster(cl))
          varlist <- unique(c(ls(), ls(envir=.GlobalEnv), ls(envir=parent.env(environment()))))
          clusterExport(cl, varlist=varlist, envir=environment())
          clusterSetRNGStream(cl)
          wd <- getwd()
          clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"), envir=environment())
          for (iter in 1:Iterations) {
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               for (j in sample.int(LIV)) {
                    if(j %in% dparm)
                         Mo0 <- .mcmcggdpp(Model, Data, j, Mo0, Grid, cl)
                    else {
                         agg <- AGGCPP(Model, Data, j, Mo0, Grid, tuning, smax, cl)
                         Mo0 <- agg$Mo0
                         tuning[j] <- agg$tuning[j]}
                    }
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter/Thinning) + 1
                    thinned[t.iter, ] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter, ] <- Mo0[["Monitor"]]
                    DiagCovar <- rbind(DiagCovar, tuning)}}}
     DiagCovar <- DiagCovar[-1,]
     out <- list(Acceptance=Iterations, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }




.mcmcahmc <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     epsilon <- Specs[["epsilon"]]
     L <- Specs[["L"]]
     Periodicity <- Specs[["Periodicity"]]
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     DiagCovar <- matrix(epsilon, nrow(thinned), ncol(thinned), byrow=TRUE)
     gr0 <- partial(Model, post[1,], Data)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- post[iter,]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          prop <- post[iter,]
          momentum0 <- rnorm(LIV)
          kinetic0 <- sum(momentum0^2) / 2
          momentum1 <- momentum0 + (epsilon/2) * gr0
          Mo0.1 <- Mo0
          for (l in 1:L) {
               prop <- prop + epsilon * momentum1
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0.1
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0.1
               if(any(Mo0.1[["parm"]] == Mo1[["parm"]])) {
                    nomove <- which(Mo0.1[["parm"]] == Mo1[["parm"]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- try(Model(prop, Data), silent=TRUE)
                    if(inherits(Mo1, "try-error")) Mo1 <- Mo0.1
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                         Mo1 <- Mo0.1}
               Mo0.1 <- Mo1
               prop <- Mo1[["parm"]]
               gr1 <- partial(Model, prop, Data)
               if(l < L) momentum1 <- momentum1 + epsilon * gr1}
          momentum1 <- momentum1 + (epsilon/2) * gr1
          momentum1 <- -momentum1
          kinetic1 <- sum(momentum1^2) / 2
          ### Accept/Reject
          H0 <- -Mo0[["LP"]] + kinetic0
          H1 <- -Mo1[["LP"]] + kinetic1
          delta <- H1 - H0
          alpha <- min(1, exp(-delta))
          if(!is.finite(alpha)) alpha <- 0
          if(runif(1) < alpha) {
               Mo0 <- Mo1
               post[iter,] <- Mo1[["parm"]]
               kinetic0 <- kinetic1
               gr0 <- gr1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          ### Adaptation
          if({iter > 10} & {iter %% Periodicity == 0}) {
               acceptances <- length(unique(post[(iter-9):iter,1]))
               if(acceptances <= 1) epsilon <- epsilon * 0.8
               else if(acceptances > 7) epsilon <- epsilon * 1.2
               if(iter %% Thinning == 0)
                    DiagCovar[t.iter,] <- epsilon}
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }
.mcmcaies <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     Nc <- Specs[["Nc"]]
     Z <- Specs[["Z"]]
     beta <- Specs[["beta"]]
     CPUs <- Specs[["CPUs"]]
     Packages <- Specs[["Packages"]]
     Dyn.libs <- Specs[["Dyn.libs"]]
     Mo0 <- list(Mo0=Mo0)
     if(is.null(Z)) {
          Z <- matrix(Mo0[[1]][["parm"]], Nc, LIV, byrow=TRUE)
          for (i in 2:Nc) {
               if(!is.null(Data[["PGF"]])) {
                    Z[i,] <- GIV(Model, Data, PGF=TRUE)
                    }
               else Z[i,] <- GIV(Model, Data)
               }
          }
     for (i in 2:Nc) Mo0[[i]] <- Model(Z[i,], Data)
     if(CPUs == 1) {
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[[1]][["parm"]]
                    Dev[t.iter] <- Mo0[[1]][["Dev"]]
                    Mon[t.iter,] <- Mo0[[1]][["Monitor"]]}
               for (i in 1:Nc) {
                    ### Propose new parameter values with stretch move
                    z <- 1 / sqrt(runif(1, 1 / beta, beta))
                    s <- sample(c(1:Nc)[-i], 1)
                    prop <- Mo0[[s]][["parm"]] +
                         z*(Mo0[[i]][["parm"]] - Mo0[[s]][["parm"]])
                    if(i == 1 & iter %% Status == 0) 
                         cat(",   Proposal: Multivariate,   LP:",  round(Mo0[[1]][["LP"]],1), "\n", sep="",  file=LogFile, append=TRUE)
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data), silent=TRUE)
                    if(inherits(Mo1, "try-error")) Mo1 <- Mo0[[i]]
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                         Mo1 <- Mo0[[i]]
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- (LIV-1)*log(z) + Mo1[["LP"]] -
                         Mo0[[i]][["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0[[i]] <- Mo1
                         if(i == 1) {
                              Acceptance <- Acceptance + 1
                              if(iter %% Thinning == 0) {
                                   thinned[t.iter,] <- Mo1[["parm"]]
                                   Dev[t.iter] <- Mo1[["Dev"]]
                                   Mon[t.iter,] <- Mo1[["Monitor"]]}
                              }
                         }
                    }
               }
          }
     else {
          detectedCores <- detectCores()
          cat("\n\nCPUs Detected:", detectedCores, "\n", file=LogFile, append=TRUE)
          if(CPUs > detectedCores) {
               cat("\nOnly", detectedCores, "will be used.\n", file=LogFile, append=TRUE)
               CPUs <- detectedCores}
          cat("\nLaplace's Demon is preparing environments for CPUs...", file=LogFile, append=TRUE)
          cat("\n##################################################\n", file=LogFile, append=TRUE)
          cl <- makeCluster(CPUs)
          cat("\n##################################################\n", file=LogFile, append=TRUE)
          on.exit(stopCluster(cl))
          varlist <- unique(c(ls(), ls(envir=.GlobalEnv), ls(envir=parent.env(environment()))))
          clusterExport(cl, varlist=varlist, envir=environment())
          clusterSetRNGStream(cl)
          wd <- getwd()
          clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"), envir=environment())
          model.wrapper <- function(x, ...)
               {
               if(!is.null(Packages)) {
                    sapply(Packages, function(x) library(x, character.only=TRUE,  quietly=TRUE))}
               if(!is.null(Dyn.libs)) {
                    sapply(Dyn.libs, function(x) dyn.load(paste(wd, x, sep = "/")))
                    on.exit(sapply(Dyn.libs, function(x) dyn.unload(paste(wd, x, sep = "/"))))}
               Model(prop[x,], Data)
               }
          prop <- Z
          batch1 <- 1:(Nc/2)
          batch2 <- batch1 + (Nc/2)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[[1]][["parm"]]
                    Dev[t.iter] <- Mo0[[1]][["Dev"]]
                    Mon[t.iter,] <- Mo0[[1]][["Monitor"]]}
               for (i in 1:Nc) {
                    ### Propose new parameter values with stretch move
                    z <- 1 / sqrt(runif(1, 1 / beta, beta))
                    if(i <= (Nc/2)) s <- sample(batch2, 1)
                    else s <- sample(batch1, 1)
                    prop[i,] <- Mo0[[s]][["parm"]] +
                         z*(Mo0[[i]][["parm"]] - Mo0[[s]][["parm"]])
                    if(i == 1 & iter %% Status == 0) 
                         cat(",   Proposal: Multivariate\n", file=LogFile,  append=TRUE)}
               ### Log-Posterior of the proposed state
               Mo1 <- clusterApply(cl, 1:Nc, model.wrapper, Model, Data, prop)
               for (i in 1:Nc) {
                    if(any(!is.finite(c(Mo1[[i]][["LP"]], Mo1[[i]][["Dev"]], Mo1[[i]][["Monitor"]]))))
                         Mo1[[i]] <- Mo0[[i]]
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- (LIV-1)*log(z) + Mo1[[i]][["LP"]] -
                         Mo0[[i]][["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0[[i]] <- Mo1[[i]]
                         if(i == 1) {
                              Acceptance <- Acceptance + 1
                              if(iter %% Thinning == 0) {
                                   thinned[t.iter,] <- Mo1[[i]][["parm"]]
                                   Dev[t.iter] <- Mo1[[i]][["Dev"]]
                                   Mon[t.iter,] <- Mo1[[i]][["Monitor"]]}
                              }
                         }
                    }
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }
.mcmcam <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning, VarCov, LogFile)
     {
     Adaptive <- Specs[["Adaptive"]]
     Periodicity <- Specs[["Periodicity"]]
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     Iden.Mat <- diag(LIV)
     DiagCovar <- matrix(diag(VarCov), floor(Iterations/Periodicity), LIV, byrow=TRUE)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- post[iter,]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          MVNz <- try(rbind(rnorm(LIV)) %*% chol(VarCov), silent=TRUE)
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               MVNz <- as.vector(MVNz)
               prop <- t(post[iter,] + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               prop <- post[iter,]
               j <- ceiling(runif(1,0,LIV))
               prop[j] <- rnorm(1, post[iter,j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               post[iter,] <- Mo1[["parm"]]
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}}
          ### Shrinkage of Adaptive Proposal Variance
          if({Adaptive < Iterations} & {Acceptance > 5} &
               {Acceptance / iter < 0.05}) {
               VarCov <- VarCov * {1 - {1 / Iterations}}
               tuning <- tuning * {1 - {1 / Iterations}}}
          ### Adapt the Proposal Variance
          if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
               ### Covariance Matrix (Preferred if it works)
               VarCov <- {ScaleF * .Call("covar", post[1:iter,], PACKAGE="LaplacesDemonCpp")} +
                    {ScaleF * 1.0E-5 * Iden.Mat}
               a.iter <- floor(iter / Periodicity)
               DiagCovar[a.iter,] <- diag(VarCov)
               ### Univariate Standard Deviations
               for (j in 1:LIV) {
                    tuning[j] <- sqrt(ScaleF * {var(post[1:iter,j])} +
                         ScaleF * 1.0E-5)}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=VarCov)
     return(out)
     }
.mcmcamm <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning, VarCov, LogFile)
     {
     Adaptive <- Specs[["Adaptive"]]
     Block <- Specs[["B"]]
     n <- Specs[["n"]]
     Periodicity <- Specs[["Periodicity"]]
     w <- Specs[["w"]]
     obs.sum <- matrix(Mo0[["parm"]]*n, LIV, 1)
     obs.scatter <- tcrossprod(Mo0[["parm"]])*n
     if(all(upper.triangle(VarCov) == 0)) prop.R <- NULL
     else prop.R <- ScaleF * chol(VarCov)
     tuning <- sqrt(0.0001 * ScaleF)
     DiagCovar <- matrix(diag(VarCov), floor(Iterations/Periodicity), LIV, byrow=TRUE)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values from a mixture
          if(is.null(prop.R) || runif(1) < w) {
               prop <- rnorm(LIV, Mo0[["parm"]], tuning)
               if(iter %% Status == 0) 
                    cat(",   Proposal: Non-Adaptive Component,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)}
          else {
               prop <- Mo0[["parm"]] +
                    as.vector(rbind(rnorm(LIV)) %*% prop.R)
               if(iter %% Status == 0) 
                    cat(",   Proposal: Adaptive Component,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)}
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}}
          ### Update Sample and Scatter Sum
          obs.sum <- obs.sum + Mo0[["parm"]]
          obs.scatter <- obs.scatter + tcrossprod(Mo0[["parm"]])
          ### Adapt the Proposal Variance
          if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
               VarCov <- obs.scatter/{n + iter} -
                    tcrossprod(obs.sum/{n + iter})
               diag(VarCov) <- diag(VarCov) + 1e-05
               a.iter <- floor(iter / Periodicity)
               DiagCovar[a.iter,] <- diag(VarCov)
               prop.R <- try(ScaleF * chol(VarCov), silent=TRUE)
               if(!is.matrix(prop.R)) prop.R <- NULL}
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=VarCov)
     return(out)
     }
.mcmcamm.b <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning, VarCov, LogFile)
     {
     Adaptive <- Specs[["Adaptive"]]
     Block <- Specs[["B"]]
     n <- Specs[["n"]]
     Periodicity <- Specs[["Periodicity"]]
     w <- Specs[["w"]]
     B <- length(Block)
     if(!identical(length(VarCov), B))
          stop("Number of components in Covar differs from number of blocks.")
     obs.scatter <- obs.sum <- prop.R <- list()
     for (b in 1:B) {
          if(!identical(length(Block[[b]]), length(diag(VarCov[[b]]))))
               stop("Diagonal of Covar[[",b,"]] differs from block length.")
          obs.sum[[b]] <- matrix(Mo0[["parm"]][Block[[b]]]*n, length(Block[[b]]), 1)
          obs.scatter[[b]] <- matrix(tcrossprod(Mo0[["parm"]][Block[[b]]])*n, length(Block[[b]]), length(Block[[b]]))
          if(all(upper.triangle(VarCov[[b]]) == 0)) prop.R[[b]] <- NA
          else prop.R[[b]] <- ScaleF * chol(VarCov[[b]])}
     tuning <- sqrt(0.0001 * ScaleF)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Proceed by Block
          for (b in 1:B) {
               ### Propose new parameter values from a mixture
               prop <- Mo0[["parm"]]
               if(any(is.na(prop.R[[b]])) || runif(1) < w) {
                    prop[Block[[b]]] <- rnorm(length(Block[[b]]), Mo0[["parm"]][Block[[b]]], tuning)
                    if(b == 1 & iter %% Status == 0) 
                         cat(",   Proposal: Non-Adaptive Component,   LP:",  round(Mo0[["LP"]],1), "\n", sep="",  file=LogFile, append=TRUE)}
               else {
                    prop[Block[[b]]] <- Mo0[["parm"]][Block[[b]]] +
                         as.vector(rbind(rnorm(length(Block[[b]]))) %*%
                              prop.R[[b]])
                    if(b == 1 & iter %% Status == 0) 
                         cat(",   Proposal: Adaptive Component,   LP:",  round(Mo0[["LP"]],1), "\n", sep="",  file=LogFile, append=TRUE)}
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               log.u <- log(runif(1))
               log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo1
                    Acceptance <- Acceptance + length(Block[[b]]) / LIV
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- Mo1[["parm"]]
                         Dev[t.iter] <- Mo1[["Dev"]]
                         Mon[t.iter,] <- Mo1[["Monitor"]]}}
               ### Update Sample and Scatter Sum
               obs.sum[[b]] <- obs.sum[[b]] + Mo0[["parm"]][Block[[b]]]
               obs.scatter[[b]] <- obs.scatter[[b]] +
                    tcrossprod(Mo0[["parm"]][Block[[b]]])
               ### Adapt the Proposal Variance
               if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
                    VarCov[[b]] <- obs.scatter[[b]]/{n + iter} -
                         tcrossprod(obs.sum[[b]]/{n + iter})
                    diag(VarCov[[b]]) <- diag(VarCov[[b]]) + 1e-05
                    if(b == 1) DiagCovar <- rbind(DiagCovar, rep(0,LIV))
                    DiagCovar[nrow(DiagCovar),Block[[b]]] <- diag(VarCov[[b]])
                    prop.R[[b]] <- try(ScaleF * chol(VarCov[[b]]), silent=TRUE)
                    if(!is.matrix(prop.R[[b]])) prop.R[[b]] <- NA}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=VarCov)
     return(out)
     }
.mcmcamwg <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning, LogFile)
     {
     Periodicity <- Specs[["Periodicity"]]
     Acceptance <- matrix(0, 1, LIV)
     DiagCovar <- matrix(tuning, floor(Iterations/Periodicity), LIV, byrow=TRUE)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Random-Scan Componentwise Estimation
          propdraw <- rnorm(LIV)*tuning
          for (j in sample.int(LIV)) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- prop[j] + propdraw[j]
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < {Mo1[["LP"]] - Mo0[["LP"]]}
               if(u == TRUE) {
                    Mo0 <- Mo1
                    Acceptance[j] <- Acceptance[j] + 1}}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Adapt the Proposal Variance
          if(iter %% Periodicity == 0) {
               size <- 1 / min(100, sqrt(iter))
               Acceptance.Rate <- Acceptance / iter
               log.tuning <- log(tuning)
               tuning.num <- which(Acceptance.Rate > 0.44)
               log.tuning[tuning.num] <- log.tuning[tuning.num] + size
               log.tuning[-tuning.num] <- log.tuning[-tuning.num] - size
               tuning <- exp(log.tuning)
               a.iter <- floor(iter / Periodicity)
               DiagCovar[a.iter,] <- tuning}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)), Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=tuning)
     return(out)
     }
.mcmccharm <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     alpha.star <- Specs[["alpha.star"]]
     if(is.na(alpha.star)) {
          Acceptance <- matrix(0, 1, LIV)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Random-Scan Componentwise Estimation
               theta <- rnorm(LIV)
               theta <- theta / sqrt(sum(theta*theta))
               lambda <- runif(1)
               for (j in sample.int(LIV)) {
                    ### Propose new parameter values
                    prop <- Mo0[["parm"]]
                    prop[j] <- prop[j] + lambda*theta[j]
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data), silent=TRUE)
                    if(inherits(Mo1, "try-error")) Mo1 <- Mo0
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                         Mo1 <- Mo0
                    ### Accept/Reject
                    u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
                    if(u == TRUE) {
                         Mo0 <- Mo1
                         Acceptance[j] <- Acceptance[j] + 1}}
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               }
          ### Output
          out <- list(Acceptance=mean(as.vector(Acceptance)), Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
          return(out)
          }
     else {
          tau <- rep(1, LIV)
          Acceptance <- matrix(0, 1, LIV)
          DiagCovar <- matrix(tau, nrow(thinned), LIV)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Random-Scan Componentwise Estimation
               theta <- rnorm(LIV)
               theta <- theta / sqrt(sum(theta*theta))
               lambda <- runif(1)
               for (j in sample.int(LIV)) {
                    ### Propose new parameter values
                    prop <- Mo0[["parm"]]
                    prop[j] <- prop[j] + tau[j]*lambda*theta[j]
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data), silent=TRUE)
                    if(inherits(Mo1, "try-error")) Mo1 <- Mo0
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                         Mo1 <- Mo0
                    ### Accept/Reject
                    u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
                    if(u == TRUE) {
                         Mo0 <- Mo1
                         Acceptance[j] <- Acceptance[j] + 1
                         tau[j] <- tau[j] + (tau[j] / (alpha.star *
                         (1 - alpha.star))) * (1 - alpha.star) / iter
                         }
                    else {
                         tau[j] <- abs(tau[j] - (tau[j] / (alpha.star *
                         (1 - alpha.star))) * alpha.star / iter)}}
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]
                    DiagCovar[t.iter,] <- tau}
               }
          ### Output
          out <- list(Acceptance=mean(as.vector(Acceptance)), Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
          return(out)
          }
     }
.mcmcdemc <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     Nc <- Specs[["Nc"]]
     Z <- Specs[["Z"]]
     gamma <- Specs[["gamma"]]
     w <- Specs[["w"]]
     const <- 2.381204 / sqrt(2)
     Mo0 <- list(Mo0=Mo0)
     if(is.null(Z)) {
          cat("\nGenerating Z...\n", file=LogFile, append=TRUE)
          Z <- array(0, dim=c(floor(Iterations/Thinning)+1, LIV, Nc))
          for (t in 1:dim(Z)[1]) {
               for (i in 1:Nc) {
                    if(t == 1 & i == 1) {
                         Z[t,,i] <- Mo0[[1]][["parm"]]
                         }
                    else {
                         if(!is.null(Data[["PGF"]])) {
                              Z[t,,i] <- GIV(Model, Data, PGF=TRUE)}
                         else Z[t,,i] <- GIV(Model, Data)
                         }
                    }
               }
          }
     else Z[1,,1] <- Mo0[[1]][["parm"]]
     for (i in 2:Nc) Mo0[[i]] <- Model(Z[1,,i], Data)
     for (iter in 1:Iterations) {
          ### Thinned Iteration
          t.iter <- floor(iter / Thinning) + 1
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               Z[t.iter,,] <- Z[t.iter-1,,]
               thinned[t.iter,] <- Mo0[[1]][["parm"]]
               Dev[t.iter] <- Mo0[[1]][["Dev"]]
               Mon[t.iter,] <- Mo0[[1]][["Monitor"]]}
          omega <- runif(1)
          for (i in 1:Nc) {
               r <- sample(dim(Z)[1], 2)
               s <- sample(c(1:Nc)[-i], 2)
               if(omega > w) {
                    ### Parallel Direction Move
                    prop <- Mo0[[i]][["parm"]] +
                         gamma*(Z[r[1],,s[1]] - Z[r[2],,s[2]]) +
                         runif(LIV, -0.001, 0.001)^LIV
                         }
               else {
                    ### Snooker Move
                    si <- sample(c(1:Nc)[-i], 1)
                    prop <- Mo0[[i]][["parm"]] + const*
                         ({Mo0[[si]][["parm"]] - Z[r[1],,s[1]]} -
                          {Mo0[[si]][["parm"]] - Z[r[2],,s[2]]})}
               if(i == 1 & iter %% Status == 0) 
                    cat(",   Proposal: Multivariate,   LP:", round(Mo0[[1]][["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0[[i]]
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0[[i]]
               ### Accept/Reject
               log.u <- log(runif(1))
               log.alpha <- Mo1[["LP"]] - Mo0[[i]][["LP"]]
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0[[i]] <- Mo1
                    Z[t.iter,,i] <- Mo1[["parm"]]
                    if(i == 1) {
                         Acceptance <- Acceptance + 1
                         if(iter %% Thinning == 0) {
                              thinned[t.iter,] <- Mo1[["parm"]]
                              Dev[t.iter] <- Mo1[["Dev"]]
                              Mon[t.iter,] <- Mo1[["Monitor"]]}
                         }
                    }
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=thinned, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }
.mcmcdram <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning, VarCov, LogFile)
     {
     Adaptive <- Specs[["Adaptive"]]
     DR <- 1
     Periodicity <- Specs[["Periodicity"]]
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     Iden.Mat <- diag(LIV)
     DiagCovar <- matrix(diag(VarCov), floor(Iterations/Periodicity), LIV, byrow=TRUE)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- post[iter,]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          MVNz <- try(rbind(rnorm(LIV)) %*% chol(VarCov), silent=TRUE)
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               MVNz <- as.vector(MVNz)
               prop <- t(post[iter,] + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               prop <- post[iter,]
               j <- ceiling(runif(1,0,LIV))
               prop[j] <- rnorm(1, post[iter,j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               post[iter,] <- Mo1[["parm"]]
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          ### Delayed Rejection: Second Stage Proposals
          else if(log.u >= log.alpha) {
               MVNz <- try(rbind(rnorm(LIV)) %*%
                    chol(VarCov * 0.5), silent=TRUE)
               if(!inherits(MVNz, "try-error") &
                    ((Acceptance / iter) >= 0.05)) {
                    MVNz <- as.vector(MVNz)
                    prop <- t(post[iter,] + t(MVNz))}
               else {
                    prop <- post[iter,]
                    j <- ceiling(runif(1,0,LIV))
                    prop[j] <- rnorm(1, post[iter,j], tuning[j])}
               ### Log-Posterior of the proposed state
               Mo2 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo2, "try-error")) Mo2 <- Mo0
               else if(any(!is.finite(c(Mo2[["LP"]], Mo2[["Dev"]], Mo2[["Monitor"]]))))
                    Mo2 <- Mo0
               ### Accept/Reject
               log.u <- log(runif(1))
               options(warn=-1)
               log.alpha.comp <- log(1 - exp(Mo1[["LP"]] - Mo2[["LP"]]))
               options(warn=0)
               if(!is.finite(log.alpha.comp)) log.alpha.comp <- 0
               log.alpha <- Mo2[["LP"]] + log.alpha.comp  -
                    {Mo0[["LP"]] + log(1 - exp(Mo1[["LP"]] - Mo0[["LP"]]))}
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo2
                    post[iter,] <- Mo2[["parm"]]
                    Acceptance <- Acceptance + 1
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- Mo1[["parm"]]
                         Dev[t.iter] <- Mo1[["Dev"]]
                         Mon[t.iter,] <- Mo1[["Monitor"]]}
                    }
               }
          ### Shrinkage of Adaptive Proposal Variance
          if({Adaptive < Iterations} & {Acceptance > 5} &
               {Acceptance / iter < 0.05}) {
               VarCov <- VarCov * {1 - {1 / Iterations}}
               tuning <- tuning * {1 - {1 / Iterations}}}
          ### Adapt the Proposal Variance
          if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
               ### Covariance Matrix (Preferred if it works)
               VarCov <- {ScaleF * .Call("covar", post[1:iter,], PACKAGE="LaplacesDemonCpp")} +
                    {ScaleF * 1.0E-5 * Iden.Mat}
               a.iter <- floor(iter / Periodicity)
               DiagCovar[a.iter,] <- diag(VarCov)
               ### Univariate Standard Deviations
               for (j in 1:LIV) {
                    tuning[j] <- sqrt(ScaleF * {var(post[1:iter,j])} +
                         ScaleF * 1.0E-5)}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=VarCov)
     return(out)
     }
.mcmcdrm <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning, VarCov, LogFile)
     {
     DR <- 1
     U <- chol(VarCov)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          MVNz <- try(rbind(rnorm(LIV)) %*% U, silent=TRUE)
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               MVNz <- as.vector(MVNz)
               prop <- t(as.vector(Mo0[["parm"]]) + t(MVNz))}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               prop <- Mo0[["parm"]]
               j <- ceiling(runif(1,0,LIV))
               prop[j] <- rnorm(1, Mo0[["parm"]][j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          ### Delayed Rejection: Second Stage Proposals
          else if(log.u >= log.alpha) {
               MVNz <- try(rbind(rnorm(LIV)) %*%
                    chol(VarCov * 0.5), silent=TRUE)
               if(!inherits(MVNz, "try-error") &
                    ((Acceptance / iter) >= 0.05)) {
                    MVNz <- as.vector(MVNz)
                    prop <- t(as.vector(Mo0[["parm"]]) + t(MVNz))}
               else {
                    prop <- Mo0[["parm"]]
                    j <- ceiling(runif(1,0,LIV))
                    prop[j] <- rnorm(1, Mo0[["parm"]][j], tuning[j])}
               ### Log-Posterior of the proposed state
               Mo2 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo2, "try-error")) Mo2 <- Mo0
               else if(any(!is.finite(c(Mo2[["LP"]], Mo2[["Dev"]], Mo2[["Monitor"]]))))
                    Mo2 <- Mo0
               ### Accept/Reject
               log.u <- log(runif(1))
               options(warn=-1)
               log.alpha.comp <- log(1 - exp(Mo1[["LP"]] - Mo2[["LP"]]))
               options(warn=0)
               if(!is.finite(log.alpha.comp)) log.alpha.comp <- 0
               log.alpha <- Mo2[["LP"]] + log.alpha.comp  -
                    {Mo0[["LP"]] + log(1 - exp(Mo1[["LP"]] - Mo0[["LP"]]))}
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo2
                    Acceptance <- Acceptance + 1
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- Mo1[["parm"]]
                         Dev[t.iter] <- Mo1[["Dev"]]
                         Mon[t.iter,] <- Mo1[["Monitor"]]}
                    }
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=VarCov)
     return(out)
     }
.mcmcess <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov, LogFile)
     {
     Block <- Specs[["B"]]
     if(length(Block) == 0) {
          nu <- rnorm(LIV, 0, diag(VarCov))
          U <- chol(VarCov)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Propose new parameter values
               nu <- as.vector(rbind(rnorm(LIV)) %*% U)
               theta <- theta.max <- runif(1, 0, 2*pi)
               theta.min <- theta - 2*pi
               shrink <- TRUE
               log.u <- log(runif(1))
               ### Rejection Sampling
               while (shrink == TRUE) {
                    prop <- Mo0[["parm"]]*cos(theta) + nu*sin(theta)
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data), silent=TRUE)
                    if(inherits(Mo1, "try-error")) Mo1 <- Mo0
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                         Mo1 <- Mo0
                    ### Accept/Reject
                    log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0 <- Mo1
                         if(iter %% Thinning == 0) {
                              thinned[t.iter,] <- Mo0[["parm"]]
                              Dev[t.iter] <- Mo0[["Dev"]]
                              Mon[t.iter,] <- Mo0[["Monitor"]]}
                         shrink <- FALSE
                         }
                    else {
                         if(theta < 0) theta.min <- theta
                         else theta.max <- theta
                         theta <- runif(1, theta.min, theta.max)}}
               }
          }
     else {
          B <- length(Block)
          if(!identical(length(VarCov), B))
               stop("Number of components in Covar differs from number of blocks.")
          nu <- rep(NA, LIV)
          for (b in 1:B) {
               if(!identical(length(Block[[b]]), length(diag(VarCov[[b]]))))
                    stop("Diagonal of Covar[[",b,"]] differs from block length.")
               nu[Block[[b]]] <- rnorm(length(Block[[b]]), 0, diag(VarCov[[b]]))}
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Proceed by Block
               for (b in 1:B) {
                    ### Propose new parameter values
                    blen <- length(Block[[b]])
                    nu[Block[[b]]] <- as.vector(rbind(rnorm(blen)) %*%
                         chol(VarCov[[b]]))
                    theta <- theta.max <- runif(1, 0, 2*pi)
                    theta.min <- theta - 2*pi
                    shrink <- TRUE
                    log.u <- log(runif(1))
                    while (shrink == TRUE) {
                         prop <- Mo0[["parm"]]
                         prop[Block[[b]]] <- Mo0[["parm"]][Block[[b]]]*cos(theta) +
                              nu[Block[[b]]]*sin(theta)
                         ### Log-Posterior of the proposed state
                         Mo1 <- try(Model(prop, Data), silent=TRUE)
                         if(inherits(Mo1, "try-error")) Mo1 <- Mo0
                         else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],  Mo1[["Monitor"]]))))
                              Mo1 <- Mo0
                         ### Accept/Reject
                         log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                         if(!is.finite(log.alpha)) log.alpha <- 0
                         if(log.u < log.alpha) {
                              Mo0 <- Mo1
                              if(iter %% Thinning == 0) {
                                   thinned[t.iter,] <- Mo0[["parm"]]
                                   Dev[t.iter] <- Mo0[["Dev"]]
                                   Mon[t.iter,] <- Mo0[["Monitor"]]}
                              shrink <- FALSE
                              }
                         else {
                              if(theta < 0) theta.min <- theta
                              else theta.max <- theta
                              theta <- runif(1, theta.min, theta.max)}}
                    }
               }
          }
     ### Output
     out <- list(Acceptance=Iterations, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.Call("covar", thinned, PACKAGE="LaplacesDemonCpp"))
     return(out)
     }
.mcmcgibbs <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning, LogFile)
     {
     FC <- Specs[["FC"]]
     MWG <- Specs[["MWG"]]
     if(is.null(MWG)) {
          Acceptance <- Iterations
          MWGlen <- 0}
     else {
          MWGlen <- length(MWG)
          Acceptance <- matrix(0, 1, LIV)}
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Gibbs Sampling of Full Conditionals
          prop <- try(FC(Mo0[["parm"]], Data), silent=TRUE)
          if(inherits(prop, "try-error")) prop <- Mo0[["parm"]]
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          Mo0 <- Mo1
          ### Metropolis-within-Gibbs
          if(MWGlen > 0) {
               ### Random-Scan Componentwise Estimation
               for (j in sample(MWG)) {
                    ### Propose new parameter values
                    prop <- Mo0[["parm"]]
                    prop[j] <- rnorm(1, prop[j], tuning[j])
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data), silent=TRUE)
                    if(inherits(Mo1, "try-error")) Mo1 <- Mo0
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                         Mo1 <- Mo0
                    ### Accept/Reject
                    u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
                    if(u == TRUE) {
                         Mo0 <- Mo1
                         Acceptance[j] <- Acceptance[j] + 1}}
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               }
          }
     if(MWGlen > 0) Acceptance <- mean(as.vector(Acceptance[,MWG]))
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=tuning)
     return(out)
     }
.mcmcgg <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     Grid <- Specs[["Grid"]]
     dparm <- Specs[["dparm"]]
     CPUs <- Specs[["CPUs"]]
     Packages <- Specs[["Packages"]]
     Dyn.libs <- Specs[["Dyn.libs"]]
     if(CPUs == 1) {
          for (iter in 1:Iterations) {
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               for (j in sample.int(LIV)) {
                    if(j %in% dparm) Mo0 <- .mcmcggdp(Model, Data, j, Mo0, Grid)
                    else Mo0 <- .mcmcggcp(Model, Data, j, Mo0, Grid)
                    }
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter/Thinning) + 1
                    thinned[t.iter, ] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter, ] <- Mo0[["Monitor"]]}}
          }
     else {
          detectedCores <- detectCores()
          cat("\n\nCPUs Detected:", detectedCores, "\n", file=LogFile, append=TRUE)
          if(CPUs > detectedCores) {
               cat("\nOnly", detectedCores, "will be used.\n", file=LogFile, append=TRUE)
               CPUs <- detectedCores}
          cat("\nLaplace's Demon is preparing environments for CPUs...", file=LogFile, append=TRUE)
          cat("\n##################################################\n", file=LogFile, append=TRUE)
          cl <- makeCluster(CPUs)
          cat("\n##################################################\n", file=LogFile, append=TRUE)
          on.exit(stopCluster(cl))
          varlist <- unique(c(ls(), ls(envir=.GlobalEnv), ls(envir=parent.env(environment()))))
          clusterExport(cl, varlist=varlist, envir=environment())
          clusterSetRNGStream(cl)
          wd <- getwd()
          clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"), envir=environment())
          for (iter in 1:Iterations) {
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               for (j in sample.int(LIV)) {
                    if(j %in% dparm)
                         Mo0 <- .mcmcggdpp(Model, Data, j, Mo0, Grid, cl)
                    else Mo0 <- .mcmcggcpp(Model, Data, j, Mo0, Grid, cl)
                    }
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter/Thinning) + 1
                    thinned[t.iter, ] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter, ] <- Mo0[["Monitor"]]}}}
     out <- list(Acceptance=Iterations, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }
### Griddy-Gibbs Continuous Parameter (Non-Parallelized)
.mcmcggcp <- function(Model, Data, j, Mo0, Grid)
     {
     G <- length(Grid[[j]])
     LP.grid <- rep(0, G)
     prop <- Mo0[["parm"]]
     theta <- prop[j] + Grid[[j]]
     for (g in 1:G) {
          prop[j] <- theta[g]
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          LP.grid[g] <- Mo1[["LP"]]
          theta[g] <- Mo1[["parm"]][j]}
     if(all(!is.finite(LP.grid))) LP.grid <- rep(0, G)
     LP.grid[which(!is.finite(LP.grid))] <- min(LP.grid[which(is.finite(LP.grid))])
     LP.grid <- exp(LP.grid - logadd(LP.grid))
     LP.grid <- LP.grid / sum(LP.grid)
     s <- spline(theta, LP.grid, n=1000)
     s$y <- interval(s$y, 0, Inf, reflect=FALSE)
     if(length(which(s$y > 0)) == 0)
          prop[j] <- theta[which.max(LP.grid)[1]]
     else prop[j] <- sample(s$x, 1, prob=s$y)
     Mo1 <- try(Model(prop, Data), silent=TRUE)
     if(inherits(Mo1, "try-error")) Mo1 <- Mo0
     else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]])))) Mo1 <- Mo0
     Mo0 <- Mo1
     return(Mo0)
     }
### Griddy-Gibbs Continuous Parameter (Parallelized)
.mcmcggcpp <- function(Model, Data, j, Mo0, Grid, cl)
     {
     G <- length(Grid[[j]])
     LP.grid <- rep(0, G)
     LIV <- length(Mo0[["parm"]])
     prop <- matrix(Mo0[["parm"]], G, LIV, byrow=TRUE)
     prop[, j] <- prop[, j] + Grid[[j]]
     Mo1 <- parLapply(cl, 1:G, function(x) Model(prop[x,], Data))
     LP.grid <- as.vector(unlist(lapply(Mo1, function(x) x[["LP"]])))
     prop <- matrix(as.vector(unlist(lapply(Mo1, function(x) x[["parm"]]))), G, LIV, byrow=TRUE)
     theta <- prop[, j]
     if(all(!is.finite(LP.grid))) LP.grid <- rep(0, G)
     LP.grid[which(!is.finite(LP.grid))] <- min(LP.grid[which(is.finite(LP.grid))])
     LP.grid <- exp(LP.grid - logadd(LP.grid))
     LP.grid <- LP.grid / sum(LP.grid)
     s <- spline(theta, LP.grid, n=1000)
     s$y <- interval(s$y, 0, Inf, reflect=FALSE)
     prop <- Mo0[["parm"]]
     if(length(which(s$y > 0)) == 0)
          prop[j] <- theta[which.max(LP.grid)[1]]
     else prop[j] <- sample(s$x, 1, prob=s$y)
     Mo1 <- try(Model(prop, Data), silent=TRUE)
     if(inherits(Mo1, "try-error")) Mo1 <- Mo0
     else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]])))) Mo1 <- Mo0
     Mo0 <- Mo1
     return(Mo0)
     }
### Griddy-Gibbs Discrete Parameter (Non-Parallelized)
#where j is which parameter, and Grid are discrete values
.mcmcggdp <- function(Model, Data, j, Mo0, Grid)
     {
     G <- length(Grid[[j]])
     LP.grid <- rep(0, G)
     prop <- Mo0[["parm"]]
     theta <- Grid[[j]]
     for (g in 1:G) {
          prop[j] <- theta[g]
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          LP.grid[g] <- Mo1[["LP"]]
          theta[g] <- Mo1[["parm"]][j]}
     if(all(!is.finite(LP.grid))) LP.grid <- rep(0, G)
     LP.grid[which(!is.finite(LP.grid))] <- min(LP.grid[which(is.finite(LP.grid))])
     LP.grid <- exp(LP.grid - logadd(LP.grid))
     LP.grid <- LP.grid / sum(LP.grid)
     prop[j] <- sample(theta, 1, prob=LP.grid)
     Mo1 <- try(Model(prop, Data), silent=TRUE)
     if(inherits(Mo1, "try-error")) Mo1 <- Mo0
     else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]])))) Mo1 <- Mo0
     Mo0 <- Mo1
     return(Mo0)
     }
### Griddy-Gibbs Discrete Parameter (Parallelized)
.mcmcggdpp <- function(Model, Data, j, Mo0, Grid, cl)
     {
     G <- length(Grid[[j]])
     LP.grid <- rep(0, G)
     LIV <- length(Mo0[["parm"]])
     prop <- matrix(Mo0[["parm"]], G, LIV, byrow=TRUE)
     prop[, j] <- prop[, j] + Grid[[j]]
     Mo1 <- parLapply(cl, 1:G, function(x) Model(prop[x,], Data))
     LP.grid <- as.vector(unlist(lapply(Mo1, function(x) x[["LP"]])))
     prop <- matrix(as.vector(unlist(lapply(Mo1, function(x) x[["parm"]]))), G, LIV, byrow=TRUE)
     theta <- prop[, j]
     prop <- Mo0[["parm"]]
     if(all(!is.finite(LP.grid))) LP.grid <- rep(0, G)
     LP.grid[which(!is.finite(LP.grid))] <- min(LP.grid[which(is.finite(LP.grid))])
     LP.grid <- exp(LP.grid - logadd(LP.grid))
     LP.grid <- LP.grid / sum(LP.grid)
     prop[j] <- sample(theta, 1, prob=LP.grid)
     Mo1 <- try(Model(prop, Data), silent=TRUE)
     if(inherits(Mo1, "try-error")) Mo1 <- Mo0
     else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]])))) Mo1 <- Mo0
     Mo0 <- Mo1
     return(Mo0)
     }
.mcmcharm <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     alpha.star <- Specs[["alpha.star"]]
     Block <- Specs[["B"]]
     if(is.na(alpha.star) & {length(Block) == 0}) {
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Propose new parameter values
               theta <- rnorm(LIV)
               d <- theta / sqrt(sum(theta*theta))
               prop <- Mo0[["parm"]] + runif(1) * d
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               log.u <- log(runif(1))
               log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo1
                    Acceptance <- Acceptance + 1
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- Mo1[["parm"]]
                         Dev[t.iter] <- Mo1[["Dev"]]
                         Mon[t.iter,] <- Mo1[["Monitor"]]}
                    }
               }
          ### Output
          out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
          return(out)
          }
     else if(is.na(alpha.star) & {length(Block) > 0}) {
          B <- length(Block)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Proceed by Block
               for (b in 1:B) {
                    ### Propose new parameter values
                    theta <- rnorm(length(Block[[b]]))
                    d <- theta / sqrt(sum(theta*theta))
                    prop <- Mo0[["parm"]]
                    prop[Block[[b]]] <- prop[Block[[b]]] + runif(1) * d
                    if({b == 1} & {iter %% Status == 0}) 
                         cat(",   Proposal: Multivariate,   LP:",  round(Mo0[["LP"]],1), "\n", sep="",  file=LogFile, append=TRUE)
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data), silent=TRUE)
                    if(inherits(Mo1, "try-error")) Mo1 <- Mo0
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                         Mo1 <- Mo0
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0 <- Mo1
                         Acceptance <- Acceptance + length(Block[[b]]) / LIV
                         if(iter %% Thinning == 0) {
                              thinned[t.iter,] <- Mo1[["parm"]]
                              Dev[t.iter] <- Mo1[["Dev"]]
                              Mon[t.iter,] <- Mo1[["Monitor"]]}
                         }
                    }
               }
          ### Output
          out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
          return(out)
          }
     else if(length(Block) == 0) {
          tau <- 1
          DiagCovar <- matrix(tau, nrow(thinned), LIV)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Propose new parameter values
               theta <- rnorm(LIV)
               d <- theta / sqrt(sum(theta*theta))
               prop <- Mo0[["parm"]] + runif(1,0,tau) * d
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               log.u <- log(runif(1))
               log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo1
                    Acceptance <- Acceptance + 1
                    tau <- tau + (tau / (alpha.star *
                         (1 - alpha.star))) * (1 - alpha.star) / iter
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- Mo1[["parm"]]
                         Dev[t.iter] <- Mo1[["Dev"]]
                         Mon[t.iter,] <- Mo1[["Monitor"]]
                         DiagCovar[t.iter,] <- tau}
                    }
               else {
                    tau <- abs(tau - (tau / (alpha.star *
                         (1 - alpha.star))) * alpha.star / iter)
                    if(iter %% Thinning == 0) DiagCovar[t.iter,] <- tau}
               }
          ### Output
          out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
          return(out)
          }
     else {
          B <- length(Block)
          tau <- rep(1,B)
          DiagCovar <- matrix(1, nrow(thinned), LIV)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Proceed by Block
               for (b in 1:B) {
                    ### Propose new parameter values
                    theta <- rnorm(length(Block[[b]]))
                    d <- theta / sqrt(sum(theta*theta))
                    prop <- Mo0[["parm"]]
                    prop[Block[[b]]] <- prop[Block[[b]]] +
                         runif(1,0,tau[b]) * d
                    if({b == 1} & {iter %% Status == 0}) 
                         cat(",   Proposal: Multivariate,   LP:",  round(Mo0[["LP"]],1), "\n", sep="",  file=LogFile, append=TRUE)
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data), silent=TRUE)
                    if(inherits(Mo1, "try-error")) Mo1 <- Mo0
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                         Mo1 <- Mo0
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0 <- Mo1
                         Acceptance <- Acceptance + length(Block[[b]]) / LIV
                         tau[b] <- tau[b] + (tau[b] / (alpha.star *
                              (1 - alpha.star))) * (1 - alpha.star) / iter
                         if(iter %% Thinning == 0) {
                              thinned[t.iter,] <- Mo1[["parm"]]
                              Dev[t.iter] <- Mo1[["Dev"]]
                              Mon[t.iter,] <- Mo1[["Monitor"]]
                              DiagCovar[t.iter, Block[[b]]] <- tau[b]}
                         }
                    else {
                         tau[b] <- abs(tau[b] - (tau[b] / (alpha.star *
                              (1 - alpha.star))) * alpha.star / iter)
                         if(iter %% Thinning == 0)
                              DiagCovar[t.iter, Block[[b]]] <- tau[b]}
                    }
               }
          ### Output
          out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
          return(out)
          }
     }
.mcmchmc <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov, LogFile)
     {
     epsilon <- Specs[["epsilon"]]
     L <- Specs[["L"]]
     gr0 <- partial(Model, Mo0[["parm"]], Data)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          prop <- Mo0[["parm"]]
          momentum0 <- rnorm(LIV)
          kinetic0 <- sum(momentum0^2) / 2
          momentum1 <- momentum0 + (epsilon/2) * gr0
          Mo0.1 <- Mo0
          for (l in 1:L) {
               prop <- prop + epsilon * momentum1
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0.1
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0.1
               if(any(Mo0.1[["parm"]] == Mo1[["parm"]])) {
                    nomove <- which(Mo0.1[["parm"]] == Mo1[["parm"]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- try(Model(prop, Data), silent=TRUE)
                    if(inherits(Mo1, "try-error")) Mo1 <- Mo0.1
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                         Mo1 <- Mo0.1}
               Mo0.1 <- Mo1
               prop <- Mo1[["parm"]]
               gr1 <- partial(Model, prop, Data)
               if(l < L) momentum1 <- momentum1 + epsilon * gr1}
          momentum1 <- momentum1 + (epsilon/2) * gr1
          momentum1 <- -momentum1
          kinetic1 <- sum(momentum1^2) / 2
          ### Accept/Reject
          H0 <- -Mo0[["LP"]] + kinetic0
          H1 <- -Mo1[["LP"]] + kinetic1
          delta <- H1 - H0
          alpha <- min(1, exp(-delta))
          if(!is.finite(alpha)) alpha <- 0
          if(runif(1) < alpha) {
               Mo0 <- Mo1
               kinetic0 <- kinetic1
               gr0 <- gr1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=matrix(epsilon, 1, LIV), Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }
.mcmchmcda <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     A <- Specs[["A"]]
     delta <- Specs[["delta"]]
     epsilon <- Specs[["epsilon"]]
     Lmax <- Specs[["Lmax"]]
     lambda <- Specs[["lambda"]]
     leapfrog <- function(theta, r, grad, epsilon, Model, Data, Mo0)
          {
          rprime <- r + 0.5 * epsilon * grad
          thetaprime <-  theta + epsilon * rprime
          Mo1 <- try(Model(thetaprime, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          thetaprime <- Mo1[["parm"]]
          gradprime <- partial(Model, thetaprime, Data)
          rprime <- rprime + 0.5 * epsilon * gradprime
          out <- list(thetaprime=thetaprime, rprime=rprime, gradprime=gradprime, Mo1=Mo1)
          return(out)
          }
     find.reasonable.epsilon <- function(theta0, grad0, Mo0, Model, Data, LogFile)
          {
          cat("\nFinding a reasonable initial value for epsilon...", file=LogFile, append=TRUE)
          epsilon <- 0.001
          r0 <- runif(length(theta0))
          ### Figure out which direction to move epsilon
          leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data, Mo0)
          if(!is.finite(leap$Mo1[["LP"]]))
               stop("LP is not finite in find.reasonable.epsilon().", file=LogFile, append=TRUE)
          acceptprob <- exp(leap$Mo1[["LP"]] - Mo0[["LP"]] - 0.5 *
               (as.vector(leap$rprime %*% leap$rprime) -
               as.vector(r0 %*% r0)))
          a <- 2 * (acceptprob > 0.5) - 1
          ### Keep moving epsilon in that direction until acceptprob
          ### crosses 0.5
          while (acceptprob^a > 2^(-a)) {
               epsilon <- epsilon * 2^a
               leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data, Mo0)
               if(!is.finite(leap$Mo1[["LP"]]))
                    stop("LP is not finite in find.reasonable.epsilon().", file=LogFile, append=TRUE)
               acceptprob <- exp(leap$Mo1[["LP"]] - Mo0[["LP"]] - 0.5 *
                    (as.vector(leap$rprime %*% leap$rprime) -
                    as.vector(r0 %*% r0)))
               }
          cat("\nepsilon: ", round(max(epsilon,0.001),5), "\n\n", sep="", file=LogFile, append=TRUE)
          return(epsilon)
          }
     gr0 <- partial(Model, Mo0[["parm"]], Data)
     if(is.null(epsilon))
          epsilon <- find.reasonable.epsilon(Mo0[["parm"]], gr0, Mo0, Model, Data, LogFile)
     DiagCovar[1,] <- epsilon
     L <- max(1, round(lambda / epsilon))
     L <- min(L, Lmax)
     ### Dual-Averaging Parameters
     epsilonbar <- 1
     gamma <- 0.05
     Hbar <- 0
     kappa <- 0.75
     mu <- log(10*epsilon)
     t0 <- 10
     ### Begin HMCDA
     for (iter in 1:Iterations) {
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          prop <- Mo0[["parm"]]
          momentum1 <- momentum0 <- runif(LIV)
          joint <- Mo0[["LP"]] - 0.5 * as.vector(momentum0 %*% momentum0)
          L <- max(1, round(lambda / epsilon))
          L <- min(L, Lmax)
          gr1 <- gr0
          Mo0.1 <- Mo0
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Leapfrog Function
          for (l in 1:L) {
               momentum1 <- momentum1 + 0.5 * epsilon * gr1
               prop <- prop + epsilon * momentum1
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0.1
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0.1
               if(any(Mo0.1[["parm"]] == Mo1[["parm"]])) {
                    nomove <- which(Mo0.1[["parm"]] == Mo1[["parm"]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- try(Model(prop, Data), silent=TRUE)
                    if(inherits(Mo1, "try-error")) Mo1 <- Mo0.1
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                         Mo1 <- Mo0.1}
               Mo0.1 <- Mo1
               prop <- Mo1[["parm"]]
               gr1 <- partial(Model, prop, Data)
               momentum1 <- momentum1 + epsilon * gr1}
          ### Accept/Reject
          alpha <- min(1, exp(prop - 0.5 * as.vector(momentum1 %*% momentum1) - joint))
          if(!is.finite(alpha)) alpha <- 0
          if(runif(1) < alpha) {
               Mo0 <- Mo1
               gr0 <- gr1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          ### Adaptation
          if(iter > 1) {
               eta <- 1 / (iter - 1 + t0)
               Hbar <- (1 - eta) * Hbar + eta * (delta - alpha)
               if(iter <= A) {
                    epsilon <- exp(mu - sqrt(iter-1)/gamma * Hbar)
                    eta <- (iter-1)^-kappa
                    epsilonbar <- exp((1 - eta) * log(epsilonbar) +
                         eta * log(epsilon))
                    DiagCovar <- rbind(DiagCovar, epsilon)}
               else epsilon <- epsilonbar}
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }
.mcmcim <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov, LogFile)
     {
     mu <- Specs[["mu"]]
     VarCov <- as.positive.definite(as.symmetric.matrix(VarCov * 1.1))
     Omega <- as.inverse(VarCov)
     U <- chol(VarCov)
     d <- eigen(VarCov, symmetric=TRUE)$values
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          MVNz <- try(rbind(rnorm(LIV)) %*% U, silent=TRUE)
          if(!inherits(MVNz, "try-error")) {
               if(iter %% Status == 0) 
                   cat(",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               prop <- as.vector(mu) + as.vector(MVNz)}
          else {prop <- as.vector(Mo0[["parm"]])}
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          ### Importance Densities (dmvn)
          ss <- prop - mu
          z <- rowSums({ss %*% Omega} * ss)
          d1 <- sum(-0.5 * (LIV * log(2*pi) + sum(log(d))) - (0.5*z))
          ss <- Mo0[["parm"]] - mu
          z <- rowSums({ss %*% Omega} * ss)
          d0 <- sum(-0.5 * (LIV * log(2*pi) + sum(log(d))) - (0.5*z))
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]] + d1 - d0
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.Call("covar", thinned, PACKAGE="LaplacesDemonCpp"))
     return(out)
     }
.mcmcinca <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning, VarCov, LogFile)
     {
     Adaptive <- Specs[["Adaptive"]]
     Periodicity <- Specs[["Periodicity"]]
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     Iden.Mat <- diag(LIV)
     con <- get("con")
     Chains <- get("Chains")
     DiagCovar <- matrix(0, floor(Iterations/Periodicity), LIV)
     ### Store all posteriors
     INCA_iter <- 1
     INCA_first <- TRUE
     tmpMean <- numeric(LIV)
     tmpCov <- matrix(0, LIV, LIV)
     tmpAlpha <- numeric(Periodicity)
     lambda <- ScaleF
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- post[iter,]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          MVNz <- try(rbind(rnorm(LIV)) %*% chol(VarCov), silent=TRUE)
          if(!inherits(MVNz, "try-error") &
               ((Acceptance / iter) >= 0.05)) {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               prop <- as.vector(post[iter,]) + as.vector(MVNz)}
          else {
               if(iter %% Status == 0) 
                    cat(",   Proposal: Single-Component,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               prop <- post[iter,]
               j <- ceiling(runif(1,0,LIV))
               prop[j] <- rnorm(1, post[iter,j], tuning[j])}
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               post[iter,] <- Mo1[["parm"]]
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}}
          ### Save log.alpha
          if({iter %% Periodicity} == 0)
               tmpAlpha[Periodicity] <- min(1, exp(log.alpha))
          else tmpAlpha[(iter %% Periodicity)] <- min(1, exp(log.alpha))
          ### Shrinkage of Adaptive Proposal Variance
          if({iter < Adaptive} & {Acceptance > 5} & {Acceptance / iter < 0.05}) {
               VarCov <- VarCov * {1 - {1 / Iterations}}
               tuning <- tuning * {1 - {1 / Iterations}}}
          ### Adapt the Proposal Variance
          if({iter >= Adaptive} & {iter %% Periodicity == 0}) {
               select_post <- cbind(post[(iter-Periodicity+1):iter,], tmpAlpha)
               ### Ask for last posteriors to hpc_server
               tmp <- unserialize(con)
               ### Send new posteriors matrix to hpc_server      
               serialize(select_post, con)               
               if(is.matrix(tmp) && INCA_first == FALSE) {
                    for (i in 1:nrow(select_post)) {
                         tmpMean <- tmpMean + 1/(INCA_iter+1) *
                              (select_post[i, 1:LIV]-tmpMean)
                         tmpCov <- (INCA_iter-1)/INCA_iter * tmpCov +
                              1/INCA_iter *
                              tcrossprod(select_post[i, 1:LIV]-tmpMean)
                         INCA_iter <- INCA_iter + 1}
                    for (i in 1:nrow(tmp)) {
                         tmpMean <- tmpMean + 1/(INCA_iter+1) *
                              (tmp[i, 1:LIV]-tmpMean)
                         tmpCov <- (INCA_iter-1)/INCA_iter * tmpCov +
                              1/INCA_iter * 
                              tcrossprod(tmp[i, 1:LIV]-tmpMean)
                         INCA_iter <- INCA_iter + 1}
                    eta <- INCA_iter^-0.6
                    m1 <- median(select_post[, LIV+1])
                    m2 <- median(tmp[, LIV+1])
                    lambda <- exp(log(lambda) + eta * (m1 - 0.234))
                    lambda <- exp(log(lambda) + eta * (m2 - 0.234))}
               if(INCA_first == TRUE) {
                    for (i in 1:iter) {
                         tmpMean <- tmpMean + 1/(INCA_iter+1) *
                              (post[i, ]-tmpMean)
                         tmpCov <- (INCA_iter-1)/INCA_iter * tmpCov +
                              1/INCA_iter * 
                              tcrossprod(post[i, ] - tmpMean)
                         INCA_iter <- INCA_iter + 1}
                    INCA_first <- FALSE}
               VarCov <- lambda * (tmpCov + 1e-9 * Iden.Mat)
               a.iter <- floor(iter / Periodicity)
               DiagCovar[a.iter,] <- diag(VarCov)
               ### Univariate Standard Deviations
               tuning <- sqrt(diag(VarCov))}
     }
     ### Output
     out <- list(Acceptance=Acceptance,   Dev=Dev,   DiagCovar=DiagCovar,   Mon=Mon,   thinned=thinned,   VarCov=VarCov)
     return(out)
     }
.mcmcmala <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov, LogFile)
     {
     A <- Specs[["A"]]
     alpha.star <- Specs[["alpha.star"]]
     delta <- Specs[["delta"]]
     gamma.const <- Specs[["gamma"]]
     epsilon <- Specs[["epsilon"]]
     Gamm <- as.positive.definite(VarCov)
     mu <- Mo0[["parm"]]
     sigma2 <- 1 / (LIV*LIV)
     DiagCovar <- matrix(diag(Gamm), nrow(thinned), LIV)
     Iden <- diag(LIV)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          gr <- partial(Model, Mo0[["parm"]], Data)
          Dx <- {delta/max(delta, abs(gr))}*gr
          gamm <- min(gamma.const/iter, 1)
          Lambda <- Gamm + epsilon[2]*Iden
          U <- try(chol(sigma2*Lambda), silent=TRUE)
          if(inherits(U, "try-error"))
               U <- chol(as.positive.definite(sigma2*Lambda))
          prop <- as.vector((Mo0[["parm"]] +
               {sigma2/2}*as.vector(Lambda %*% Dx)*Dx) +
               rbind(rnorm(LIV)) %*% U)
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]
                    DiagCovar[t.iter,] <- diag(Lambda)}}
          ### Adapt Gamma (first, since it uses mu[t] not [t+1])
          xmu <- Mo0[["parm"]] - mu
          Gamm.prop <- Gamm + gamm*{xmu %*% t(xmu) - Gamm}
          norm.Gamm <- norm(Gamm.prop, type="F")
          if(norm.Gamm <= A) Gamm <- Gamm.prop
          else if(!is.finite(norm.Gamm)) Gamm <- sigma2*Iden
          else Gamm <- {A/norm.Gamm}*Gamm.prop
          ### Adapt mu
          mu.prop <- mu + gamm*(Mo0[["parm"]] - mu)
          norm.mu <- sqrt(sum(mu.prop*mu.prop))
          if(norm.mu <= A) mu <- mu.prop
          else mu <- {A/norm.mu}*mu.prop
          ### Adapt sigma
          sigma2 <- interval(sqrt(sigma2) +
               gamm*(min(exp(log.alpha),1) - alpha.star), epsilon[1], A, reflect=FALSE)^2
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=Lambda)
     return(out)
     }
.mcmcmcmcmc <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning, VarCov, LogFile)
     {
     lambda <- Specs[["lambda"]]
     CPUs <- Specs[["CPUs"]]
     Packages <- Specs[["Packages"]]
     Dyn.libs <- Specs[["Dyn.libs"]]
     detectedCores <- detectCores()
     cat("\n\nCPUs Detected:", detectedCores, "\n", file=LogFile, append=TRUE)
     if(CPUs > detectedCores) {
          cat("\nOnly", detectedCores, "will be used.\n", file=LogFile, append=TRUE)
          CPUs <- detectedCores}
     cat("\nLaplace's Demon is preparing environments for CPUs...", file=LogFile, append=TRUE)
     cat("\n##################################################\n", file=LogFile, append=TRUE)
     cl <- makeCluster(CPUs)
     cat("\n##################################################\n", file=LogFile, append=TRUE)
     on.exit(stopCluster(cl))
     varlist <- unique(c(ls(), ls(envir=.GlobalEnv), ls(envir=parent.env(environment()))))
     clusterExport(cl, varlist=varlist, envir=environment())
     clusterSetRNGStream(cl)
     wd <- getwd()
     clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"), envir=environment())
     if(length(lambda) == 1) Temperature <- 1/(1 + lambda*(c(1:CPUs) - 1))
     else if(length(lambda) == LIV) Temperature <- lambda
     else Temperature <- 1/(1 + lambda[1]*(c(1:CPUs) - 1))
     coolest <- which.max(Temperature)[1]
     temp <- Mo0
     Mo0 <- list()
     for (i in 1:CPUs) Mo0[[i]] <- temp
     prop <- matrix(Mo0[[1]][["parm"]], CPUs, LIV, byrow=TRUE)
     Acceptance.swap <- 0
     U <- chol(VarCov)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[[coolest]][["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[[coolest]][["parm"]]
               Dev[t.iter] <- Mo0[[coolest]][["Dev"]]
               Mon[t.iter,] <- Mo0[[coolest]][["Monitor"]]}
          ### Propose new parameter values
          for (i in 1:CPUs)
               prop[i,] <- Mo0[[i]][["parm"]] + rbind(rnorm(LIV)) %*% U
          ### Log-Posterior of the proposed state
          Mo1 <- parLapply(cl, 1:CPUs, function(x)
               try(Model(prop[x,], Data), silent=TRUE))
          for (i in 1:CPUs) {
               if(inherits(Mo1[[i]], "try-error")) Mo1[[i]] <- Mo0[[i]]
               else if(any(!is.finite(c(Mo1[[i]][["LP"]], Mo1[[i]][["Dev"]], Mo1[[i]][["Monitor"]]))))
                    Mo1[[i]] <- Mo0[[i]]}
          ### Accept/Reject
          for (i in 1:CPUs) {
               log.u <- log(runif(1))
               log.alpha <- (Mo1[[i]][["LP"]] - Mo0[[i]][["LP"]]) /
                    Temperature[i]
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0[[i]] <- Mo1[[i]]
                    if(i == coolest) {
                         Acceptance <- Acceptance + 1
                         if(iter %% Thinning == 0) {
                              thinned[t.iter,] <- Mo1[[i]][["parm"]]
                              Dev[t.iter] <- Mo1[[i]][["Dev"]]
                              Mon[t.iter,] <- Mo1[[i]][["Monitor"]]}}}}
          ### Swap
          swap <- sample.int(CPUs, 2)
          log.u <- log(runif(1))
          log.alpha <- {(Mo0[[swap[1]]][["LP"]] - Mo0[[swap[2]]][["LP"]]) /
               Temperature[swap[2]]} +
               {(Mo0[[swap[2]]][["LP"]] - Mo0[[swap[1]]][["LP"]]) /
               Temperature[swap[1]]}
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Acceptance.swap <- Acceptance.swap + 1
               temp <- Mo0[[swap[2]]]
               Mo0[[swap[2]]] <- Mo0[[swap[1]]]
               Mo0[[swap[1]]] <- temp
               if({swap[1] == coolest} & {iter %% Thinning == 0}) {
                    thinned[t.iter,] <- Mo0[[swap[1]]][["parm"]]
                    Dev[t.iter] <- Mo0[[swap[1]]][["Dev"]]
                    Mon[t.iter,] <- Mo0[[swap[1]]][["Monitor"]]}}
          }
     cat("\nSwap Acceptance Rate:", round(Acceptance.swap / Iterations, 5), "\n")
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=VarCov)
     return(out)
     }
.mcmcmtm <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned, tuning, LogFile)
     {
     K <- Specs[["K"]]
     CPUs <- Specs[["CPUs"]]
     Packages <- Specs[["Packages"]]
     Dyn.libs <- Specs[["Dyn.libs"]]
     if(CPUs > 1) {
          detectedCores <- detectCores()
          cat("\n\nCPUs Detected:", detectedCores, "\n", file=LogFile, append=TRUE)
          if(CPUs > detectedCores) {
               cat("\nOnly", detectedCores, "will be used.\n", file=LogFile, append=TRUE)
               CPUs <- detectedCores}
          cat("\nLaplace's Demon is preparing environments for CPUs...", file=LogFile, append=TRUE)
          cat("\n##################################################\n", file=LogFile, append=TRUE)
          cl <- makeCluster(CPUs)
          cat("\n##################################################\n", file=LogFile, append=TRUE)
          on.exit(stopCluster(cl))
          varlist <- unique(c(ls(), ls(envir=.GlobalEnv), ls(envir=parent.env(environment()))))
          clusterExport(cl, varlist=varlist, envir=environment())
          clusterSetRNGStream(cl)
          wd <- getwd()
          clusterExport(cl, varlist=c("Packages", "Dyn.libs", "wd"), envir=environment())}
     Acceptance <- matrix(0, 1, LIV)
     Mo1 <- list()
     for (k in 1:K) Mo1[[k]] <- Mo0
     LW <- LP <- rep(0, K)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Random-Scan Componentwise Estimation
          for (j in sample.int(LIV)) {
               ### Propose new parameter values
               prop1 <- matrix(Mo0[["parm"]], K, LIV, byrow=TRUE)
               prop1[,j] <- rnorm(K, prop1[,j], tuning[j])
               ### Log-Posterior of the proposed states
               if(CPUs == 1) {
                    ### Non-parallel
                    for (k in 1:K) {
                         Mo1[[k]] <- try(Model(prop1[k,], Data), silent=TRUE)
                         if(inherits(Mo1[[k]], "try-error"))
                              Mo1[[k]] <- Mo0
                         else if(any(!is.finite(c(Mo1[[k]][["LP"]],  Mo1[[k]][["Dev"]], Mo1[[k]][["Monitor"]]))))
                              Mo1[[k]] <- Mo0
                         LP[k] <- LW[k] <- Mo1[[k]][["LP"]]
                         prop1[k,] <- Mo1[[k]][["parm"]]}
                    }
               else {
                    ### Parallel
                    Mo1 <- parLapply(cl, 1:K, function(x)
                         try(Model(prop1[x,], Data), silent=TRUE))
                    for (k in 1:K) {
                         if(inherits(Mo1[[k]], "try-error"))
                              Mo1[[k]] <- Mo0
                         else if(any(!is.finite(c(Mo1[[k]][["LP"]],  Mo1[[k]][["Dev"]], Mo1[[k]][["Monitor"]]))))
                              Mo1[[k]] <- Mo0
                         LP[k] <- LW[k] <- Mo1[[k]][["LP"]]
                         prop1[k,] <- Mo1[[k]][["parm"]]}
                    }
               ### Normalize Weights
               w <- exp(LW - logadd(LW))
               if(all(w == 0)) w <- rep(1/K, K)
               ### Sample a Proposal
               prop5 <- Mo0[["parm"]]
               prop2 <- sample(prop1[,j], size=1, prob=w)
               prop5[j] <- prop2
               ### Create Reference Set
               Mo2 <- try(Model(prop5, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo2 <- Mo0
               prop3 <- c(rnorm(K-1, Mo2[["parm"]][j], tuning[j]), Mo2[["parm"]][j])
               prop4 <- prop1
               prop4[,j] <- prop3
               ### Calculate Acceptance Probability
               numerator <- logadd(LP)
               denom <- rep(0, K)
               if(CPUs == 1) {
                    ### Non-parallel
                    for (k in 1:K) {
                         Mo1[[k]] <- try(Model(prop4[k,], Data), silent=TRUE)
                         if(inherits(Mo1[[k]], "try-error")) Mo1[[k]] <- Mo0
                         else if(any(!is.finite(c(Mo1[[k]][["LP"]],  Mo1[[k]][["Dev"]], Mo1[[k]][["Monitor"]]))))
                              Mo1[[k]] <- Mo0
                         denom[k] <- Mo1[[k]][["LP"]]}
                    }
               else {
                    ### Parallel
                    Mo1 <- parLapply(cl, 1:K, function(x)
                         try(Model(prop4[x,], Data), silent=TRUE))
                    for (k in 1:K) {
                         if(inherits(Mo1[[k]], "try-error")) Mo1[[k]] <- Mo0
                         else if(any(!is.finite(c(Mo1[[k]][["LP"]],  Mo1[[k]][["Dev"]], Mo1[[k]][["Monitor"]]))))
                              Mo1[[k]] <- Mo0
                         denom[k] <- Mo1[[k]][["LP"]]}}
               denom <- logadd(denom)
               ### Accept/Reject
               u <- log(runif(1)) < (numerator - denom)
               if(u == TRUE) {
                    Mo0 <- Mo2
                    Acceptance[j] <- Acceptance[j] + 1}}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)), Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }


.mcmcmwg <- function(Model, Data, Iterations, Status, Thinning, Specs,
                     Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning, LogFile) {
     Acceptance <- matrix(0, 1, LIV)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Random-Scan Componentwise Estimation
          propdraw <- rnorm(LIV)*tuning
          for (j in sample.int(LIV)) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- prop[j] + propdraw[j]
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < {Mo1[["LP"]] - Mo0[["LP"]]}
               if(u == TRUE) {
                    Mo0 <- Mo1
                    Acceptance[j] <- Acceptance[j] + 1
               }
          }
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]
          }
     }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)), Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=tuning)
     return(out)
     }


.mcmcnuts <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     A <- Specs[["A"]]
     delta <- Specs[["delta"]]
     epsilon <- Specs[["epsilon"]]
     Lmax <- Specs[["Lmax"]]
     post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     leapfrog <- function(theta, r, grad, epsilon, Model, Data, Mo0)
          {
          rprime <- r + 0.5 * epsilon * grad
          thetaprime <-  theta + epsilon * rprime
          Mo1 <- try(Model(thetaprime, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          thetaprime <- Mo1[["parm"]]
          gradprime <- partial(Model, thetaprime, Data)
          rprime <- rprime + 0.5 * epsilon * gradprime
          out <- list(thetaprime=thetaprime, rprime=rprime, gradprime=gradprime, Mo1=Mo1)
          return(out)
          }
     stop.criterion <- function(thetaminus, thetaplus, rminus, rplus)
          {
          thetavec <- thetaplus - thetaminus
          criterion <- (thetavec %*% rminus >= 0) &&
               (thetavec %*% rplus >= 0)
          return(criterion)
          }
     build.tree <- function(theta, r, grad, logu, v, j, epsilon, joint0, Mo0)
          {
          if(j == 0) {
               ### Base case: Take a single leapfrog step in direction v
               leap <- leapfrog(theta, r, grad, v*epsilon, Model, Data, Mo0)
               rprime <- leap$rprime
               thetaprime <- leap$thetaprime
               Mo1 <- leap$Mo1
               gradprime <- leap$gradprime
               joint <- Mo1[["LP"]] - 0.5 * as.vector(rprime %*% rprime)
               ### Is the new point in the slice?
               nprime <- logu < joint
               ### Is the simulation wildly inaccurate?
               sprime <- logu - 1000 < joint
               # Set the return values---minus=plus for all things here, # since the "tree" is of depth 0.
               thetaminus <- thetaprime
               thetaplus <- thetaprime
               rminus <- rprime
               rplus <- rprime
               gradminus <- gradprime
               gradplus <- gradprime
               ### Compute the acceptance probability
               alphaprime <- min(1, exp(Mo1[["LP"]] - 0.5 *
                    as.vector(rprime %*% rprime) - joint0))
               nalphaprime <- 1}
          else {
               # Recursion: Implicitly build the height j-1 left and
               # right subtrees
               tree <- build.tree(theta, r, grad, logu, v, j-1, epsilon, joint0)
               thetaminus <- tree$thetaminus
               rminus <- tree$rminus
               gradminus <- tree$gradminus
               thetaplus <- tree$thetaplus
               rplus <- tree$rplus
               gradplus <- tree$gradplus
               thetaprime <- tree$thetaprime
               gradprime <- tree$gradprime
               Mo1 <- tree$Mo1
               nprime <- tree$nprime
               sprime <- tree$sprime
               alphaprime <- tree$alphaprime
               nalphaprime <- tree$nalphaprime
               ### If the first subtree stopping criterion is met, then stop
               if(sprime == 1) {
                    if(v == -1) {
                         tree <- build.tree(thetaminus, rminus, gradminus,  logu, v, j-1, epsilon, joint0)
                         thetaminus <- tree$thetaminus
                         rminus <- tree$rminus
                         gradminus <- tree$gradminus
                         thetaprime2 <- tree$thetaprime
                         gradprime2 <- tree$gradprime
                         Mo12 <- tree$Mo1
                         nprime2 <- tree$nprime
                         sprime2 <- tree$sprime
                         alphaprime2 <- tree$alphaprime
                         nalphaprime2 <- tree$nalphaprime
                         }
                    else {
                         tree <- build.tree(thetaplus, rplus, gradplus,  logu, v, j-1, epsilon, joint0)
                         thetaplus <- tree$thetaplus
                         rplus <- tree$rplus
                         gradplus <- tree$gradplus
                         thetaprime2 <- tree$thetaprime
                         gradprime2 <- tree$gradprime
                         Mo12 <- tree$Mo1
                         nprime2 <- tree$nprime
                         sprime2 <- tree$sprime
                         alphaprime2 <- tree$alphaprime
                         nalphaprime2 <- tree$nalphaprime
                         }
                    ### Choose a subtree to propagate a sample up from
                    temp <- nprime2 / (nprime + nprime2)
                    if(!is.finite(temp)) temp <- 0
                    if(runif(1) < temp) {
                         thetaprime <- thetaprime2
                         gradprime <- gradprime2
                         Mo1 <- Mo12}
                    ### Update the number of valid points
                    nprime <- nprime + nprime2
                    ### Update the stopping criterion
                    sprime <- sprime && sprime2 &&
                         stop.criterion(thetaminus, thetaplus, rminus,  rplus)
                    ### Update acceptance probability statistics
                    alphaprime <- alphaprime + alphaprime2
                    nalphaprime <- nalphaprime + nalphaprime2}}
          out <- list(thetaminus=thetaminus, rminus=rminus, gradminus=gradminus, thetaplus=thetaplus, rplus=rplus, gradplus=gradplus, thetaprime=thetaprime, gradprime=gradprime, Mo1=Mo1, nprime=nprime, sprime=sprime, alphaprime=alphaprime, nalphaprime=nalphaprime)
          return(out)
          }
     find.reasonable.epsilon <- function(theta0, grad0, Mo0, Model, Data, LogFile)
          {
          cat("\nFinding a reasonable initial value for epsilon...", file=LogFile, append=TRUE)
          epsilon <- 0.001
          r0 <- runif(length(theta0))
          ### Figure out which direction to move epsilon
          leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data, Mo0)
          if(!is.finite(leap$Mo1[["LP"]]))
               stop("LP is not finite in find.reasonable.epsilon().", file=LogFile, append=TRUE)
          acceptprob <- exp(leap$Mo1[["LP"]] - Mo0[["LP"]] - 0.5 *
               (as.vector(leap$rprime %*% leap$rprime) -
               as.vector(r0 %*% r0)))
          a <- 2 * (acceptprob > 0.5) - 1
          ### Keep moving epsilon in that direction until acceptprob
          ### crosses 0.5
          while (acceptprob^a > 2^(-a)) {
               epsilon <- epsilon * 2^a
               leap <- leapfrog(theta0, r0, grad0, epsilon, Model, Data, Mo0)
               if(!is.finite(leap$Mo1[["LP"]]))
                    stop("LP is not finite in find.reasonable.epsilon().", file=LogFile, append=TRUE)
               acceptprob <- exp(leap$Mo1[["LP"]] - Mo0[["LP"]] - 0.5 *
                    (as.vector(leap$rprime %*% leap$rprime) -
                    as.vector(r0 %*% r0)))
               }
          cat("\nepsilon: ", round(max(epsilon,0.001),5), "\n\n", sep="", file=LogFile, append=TRUE)
          return(epsilon)
          }
     Count <- 0
     evals <- 0
     grad <- partial(Model, post[1,], Data)
     if(is.null(epsilon))
          epsilon <- find.reasonable.epsilon(post[1,], grad, Mo0, Model, Data, LogFile)
     DiagCovar[1,] <- epsilon
     ### Dual-Averaging Parameters
     epsilonbar <- 1
     gamma <- 0.05
     Hbar <- 0
     kappa <- 0.75
     mu <- log(10*epsilon)
     t0 <- 10
     ### Reset Dev, Mon, and thinned
     if(A < Iterations) {
          Dev <- matrix(Dev[1:(floor((Iterations-A)/Thinning)+1),])
          Mon <- matrix(Mo0[["Monitor"]], floor((Iterations-A)/Thinning)+1, length(Mo0[["Monitor"]]), byrow=TRUE)
          thinned <- matrix(0, floor((Iterations-A)/Thinning)+1, LIV)}
     ### Begin NUTS
     for (iter in 2:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Current Posterior
          if(iter > 1) post[iter,] <- post[iter-1,]
          ### Save Thinned Samples
          if(iter > A) {
               if((iter-A) %% Thinning == 0) {
                    thinned[((iter-A)/Thinning+1),] <- post[iter,]
                    Dev[((iter-A)/Thinning+1)] <- Mo0[["Dev"]]
                    Mon[((iter-A)/Thinning+1),] <- Mo0[["Monitor"]]}}
          else if(A >= Iterations) {
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- post[iter,]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}}
          prop <- post[iter,]
          r0 <- runif(LIV) ### r0 is momenta
          ### Joint log-probability of theta and momenta r
          joint <- Mo0[["LP"]] - 0.5 * as.vector(r0 %*% r0)
          ### Resample u ~ U([0, exp(joint)])
          logu <- joint - rexp(1)
          ### Initialize Tree
          thetaminus <- prop
          thetaplus <- prop
          rminus <- r0
          rplus <- r0
          gradminus <- grad
          gradplus <- grad
          j <- 0 ### Initial height j=0
          n <- 1 ### Initially, the only valid point is the initial point
          s <- 1 ### Loop until s == 0
          while (s == 1) {
               ### Choose a direction: -1=backward, 1=forward.
               v <- 2*(runif(1) < 0.5) - 1
               ### Double the size of the tree.
               if(v == -1) {
                    tree <- build.tree(thetaminus, rminus, gradminus, logu, v, j, epsilon, joint)
                    thetaminus <- tree$thetaminus
                    rminus <- tree$rminus
                    gradminus <- tree$gradminus
                    thetaprime <- tree$thetaprime
                    gradprime <- tree$gradprime
                    Mo1 <- tree$Mo1
                    nprime <- tree$nprime
                    sprime <- tree$sprime
                    alpha <- tree$alphaprime
                    nalpha <- tree$nalphaprime}
               else {
                    tree <- build.tree(thetaplus, rplus, gradplus, logu, v, j, epsilon, joint)
                    thetaplus <- tree$thetaplus
                    rplus <- tree$rplus
                    gradplus <- tree$gradplus
                    thetaprime <- tree$thetaprime
                    gradprime <- tree$gradprime
                    Mo1 <- tree$Mo1
                    nprime <- tree$nprime
                    sprime <- tree$sprime
                    alpha <- tree$alphaprime
                    nalpha <- tree$nalphaprime}
               ### Accept/Reject
               Count <- Count + 1
               if((sprime == 1) && (runif(1) < nprime/n)) {
                    post[iter,] <- thetaprime
                    Mo0 <- Mo1
                    grad <- gradprime
                    Acceptance <- Acceptance + 1
                    if(iter > A) {
                         if((iter-A) %% Thinning == 0) {
                              thinned[((iter-A)/Thinning+1),] <- Mo1[["parm"]]
                              Dev[((iter-A)/Thinning+1)] <- Mo1[["Dev"]]
                              Mon[((iter-A)/Thinning+1),] <- Mo1[["Monitor"]]}}
                    else if(A >= Iterations) {
                         if(iter %% Thinning == 0) {
                              thinned[t.iter,] <- Mo1[["parm"]]
                              Dev[t.iter] <- Mo1[["Dev"]]
                              Mon[t.iter,] <- Mo1[["Monitor"]]}}}
               ### Update number of observed valid points
               n <- n + nprime
               ### Decide if it is time to stop
               s <- sprime &&
                    stop.criterion(thetaminus, thetaplus, rminus, rplus)
               ### Increment depth
               j <- j + 1
               if(j*j >= Lmax) s <- 0}
          ### Adaptation of epsilon
          eta <- 1 / (iter - 1 + t0)
          Hbar <- (1 - eta) * Hbar + eta * (delta - alpha / nalpha)
          if(iter <= A) {
               epsilon <- exp(mu - sqrt(iter-1)/gamma * Hbar)
               eta <- (iter-1)^-kappa
               epsilonbar <- exp((1 - eta) * log(epsilonbar) +
                    eta * log(epsilon))
               DiagCovar <- rbind(DiagCovar, epsilon)}
          else epsilon <- epsilonbar
          }
     Acceptance <- round(Acceptance / Count * Iterations)
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }
.mcmcohss <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov, LogFile)
     {
     A <- Specs[["A"]]
     n <- Specs[["n"]]
     w <- 0.05 # as with Roberts & Rosenthal
     decomp.freq <- max(floor(Iterations / Thinning / 100), 10)
     S.eig <-try(eigen(VarCov), silent=TRUE)
     if(inherits(S.eig, "try-error")) S.eig <- NULL
     tuning <- 1 #Tuning
     edge.scale <- 5 #Tuning
     t.iter <- 1
     DiagCovar <- matrix(0, floor(Iterations/Thinning), LIV, byrow=TRUE)
     if(A > Iterations)
          post <- matrix(Mo0[["parm"]], Iterations, LIV, byrow=TRUE)
     else post <- matrix(Mo0[["parm"]], A, LIV, byrow=TRUE)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Eigenvectors of the Sample Covariance Matrix
          if({iter %% decomp.freq == 0} & {iter > 2} & {iter <= A}) {
               VarCov2 <- try({VarCov*n +  .Call("covar", post[1:(iter-1),], PACKAGE="LaplacesDemonCpp")*(iter-1)}/{n+iter-1}, silent=TRUE)
               if(inherits(VarCov2, "try-error")) VarCov2 <- VarCov
               if(!is.symmetric.matrix(VarCov2))
                    VarCov2 <- as.symmetric.matrix(VarCov2)
               if(!is.positive.definite(VarCov2))
                    VarCov2 <- as.positive.definite(VarCov2)
               S.eig <- try(eigen(VarCov2), silent=TRUE)
               if(inherits(S.eig, "try-error")) S.eig <- eigen(VarCov)}
          ### Hypercube or Eigenvector
          if(runif(1) < w || is.null(S.eig)) {
               vals <- rep(tuning, LIV)
               vecs <- diag(1, nrow=LIV)
               }
          else {
               vals <- S.eig$values
               vecs <- S.eig$vectors}
          ### Slice Interval
          Mo0.1 <- try(Model(Mo0[["parm"]], Data), silent=TRUE)
          if(inherits(Mo0.1, "try-error")) Mo0.1 <- Mo0
          Mo0 <- Mo0.1
          y.slice <- Mo0[["LP"]] - rexp(1)
          L <- -1 * runif(LIV)
          U <- L + 1
          ### Rejection Sampling
          repeat {
               wt <- runif(LIV, min=L, max=U)
               v <- as.numeric(vecs %*% {edge.scale * wt * vals})
               prop <- Mo0[["parm"]] + v
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               if(Mo1[["LP"]] >= y.slice) break
               else if(all(abs(wt) < 1e-100)) {
                    Mo1 <- Mo0
                    break}
               L[wt < 0] <- wt[wt < 0]
               U[wt > 0] <- wt[wt > 0]}
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo1[["parm"]]
               Dev[t.iter] <- Mo1[["Dev"]]
               Mon[t.iter,] <- Mo1[["Monitor"]]
               DiagCovar[floor(iter / Thinning),] <- diag(S.eig$vectors)}
          Mo0 <- Mo1
          if(iter <= A) post[iter,] <- Mo0[["parm"]]
          }
     if(A > 0)
          VarCov <- VarCov2
     ### Output
     out <- list(Acceptance=Iterations, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=VarCov)
     return(out)
     }
.mcmcpcn <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov, LogFile)
     {
     beta <- Specs[["beta"]]
     U <- chol(VarCov)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          prop <- as.vector(sqrt(1 - beta*beta)*Mo0[["parm"]] +
               beta*(rbind(rnorm(LIV)) %*% U))
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=VarCov)
     return(out)
     }
.mcmcram <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov, LogFile)
     {
     alpha.star <- Specs[["alpha.star"]]
     Dist <- Specs[["Dist"]]
     gamma <- Specs[["gamma"]]
     if(!is.symmetric.matrix(VarCov)) {
          cat("\nAsymmetric Covar, correcting now...\n", file=LogFile, append=TRUE)
          VarCov <- as.symmetric.matrix(VarCov)}
     if(!is.positive.definite(VarCov)) {
          cat("\nNon-Positive-Definite Covar, correcting now...\n", file=LogFile, append=TRUE)
          VarCov <- as.positive.definite(VarCov)}
     Iden.Mat <- diag(LIV)
     S.z <- try(t(chol(VarCov)), silent=TRUE)
     if(!inherits(S.z, "try-error")) S <- S.z
     else S <- Iden.Mat
     DiagCovar <- matrix(diag(VarCov), floor(Iterations/Thinning), LIV, byrow=TRUE)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose New Parameter Values
          if(Dist == "t") U <- rt(LIV, df=5)
          else U <- rnorm(LIV)
          prop <- Mo0[["parm"]] + rbind(U) %*% S
          ### Log-Posterior
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          ### Accept/Reject
          log.u <- log(runif(1))
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log.u < log.alpha) {
               Mo0 <- Mo1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}}
          ### Adaptation
          eta <- min(1, LIV*iter^(-gamma))
          VarCov.test <- S %*% {Iden.Mat +
               eta*(min(1, exp(log.alpha)) - alpha.star) *
               U %*% t(U) / sum(U*U)} %*% t(S)
          if(missing(VarCov.test) || !all(is.finite(VarCov.test)) ||
               !is.matrix(VarCov.test)) {VarCov.test <- VarCov}
          if(!is.symmetric.matrix(VarCov.test))
               VarCov.test <- as.symmetric.matrix(VarCov.test)
          if(is.positive.definite(VarCov.test)) {
               S.z <- try(t(chol(VarCov)), silent=TRUE)
               if(!inherits(S.z, "try-error")) {
                    VarCov <- VarCov.test
                    S <- S.z}}
          DiagCovar[floor(iter / Thinning),] <- diag(VarCov)
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=VarCov)
     return(out)
     }
.mcmcrdmh <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     Acceptance <- matrix(0, 1, LIV)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          s <- sample(c(-1,1), LIV, replace=TRUE)
          u1 <- runif(LIV, -1, 1)
          epsilon1 <- u1^s
          ### Random-Scan Componentwise Estimation
          for (j in sample.int(LIV)) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- prop[j]*epsilon1[j]
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               epsilon2 <- log(abs(Mo1[["parm"]][j] / Mo0[["parm"]][j]))
               if(!is.finite(epsilon2)) epsilon2 <- 0
               ### Accept/Reject
               u2 <- log(runif(1)) < (epsilon2 + Mo1[["LP"]] - Mo0[["LP"]])
               if(u2 == TRUE) {
                    Mo0 <- Mo1
                    Acceptance[j] <- Acceptance[j] + 1}}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)), Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }
.mcmcrefractive <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned, LogFile)
     {
     Adaptive <- Specs[["Adaptive"]]
     m <- Specs[["m"]]
     w <- Specs[["w"]]
     r <- Specs[["r"]]
     alpha.star <- 0.65
     if(Adaptive < Iterations) DiagCovar <- matrix(w, nrow(thinned), LIV)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          prop <- Mo0[["parm"]]
          p <- rnorm(LIV)
          a <- 1
          g <- partial(Model, prop, Data)
          for (i in 1:m) {
               if(t(p) %*% g > 0) {
                    u <- g / sqrt(sum(g*g))
                    r1 <- 1
                    r2 <- r
                    }
               else {
                    u <- -g / sqrt(sum(g*g))
                    r1 <- r
                    r2 <- 1}
               cos.theta.1 <- (t(p) %*% u) / sqrt(sum(p*p))
               cos.2.theta.1 <- cos.theta.1 * cos.theta.1
               cos.2.theta.2 <- 1 - (r1^2 / r2^2)*(1 - cos.2.theta.1)
               if(cos.2.theta.2 > 0) cos.theta.2 <- sqrt(cos.2.theta.2)
               else cos.theta.2 <- -sqrt(abs(cos.2.theta.2))
               if(cos.2.theta.2 < 0) p <- as.vector(p - 2*(t(p) %*% u) %*% u)
               else {
                    p <- (r1 / r2)*p -
                         sqrt(sum(p*p))*((r1 / r2)*cos.theta.1 -
                         cos.theta.2)*u
                    a <- (r1 / r2)^(LIV-1)*(cos.theta.1 / cos.theta.2)*a}
               prop <- prop + w*p
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               prop <- Mo1[["parm"]]}
          ### Accept/Reject
          log.alpha <- Mo1[["LP"]] - Mo0[["LP"]] + exp(a)
          if(!is.finite(log.alpha)) log.alpha <- 0
          if(log(runif(1)) < log.alpha) {
               Mo0 <- Mo1
               Acceptance <- Acceptance + 1
               if(Adaptive < Iterations)
                    w <- w + (w / (alpha.star * (1 - alpha.star))) *
                         (1 - alpha.star) / iter
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]
                    if(Adaptive < Iterations) DiagCovar[t.iter,] <- w}
               }
          else if(Adaptive < Iterations) {
               w <- abs(w - (w / (alpha.star * (1 - alpha.star))) *
                    alpha.star / iter)
               if(iter %% Thinning == 0) DiagCovar[t.iter,] <- w}
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }
.mcmcrj <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     bin.n <- Specs[["bin.n"]]
     bin.p <- Specs[["bin.p"]]
     parm.p <- Specs[["parm.p"]]
     selectable <- Specs[["selectable"]]
     selected <- Specs[["selected"]]
     cur.parm <- cur.sel <- selected
     cur.parm[which(selectable == 0)] <- 1
     nonzero.post <- rep(0, LIV)
     p <- parm.p
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose a variable to include/exclude
          v.change <- sample(LIV, 1, prob=selectable)
          prop.sel <- cur.sel
          prop.parm <- cur.parm
          ### Change proposed size, but not above bin.n
          if(sum(cur.sel) < bin.n) {
               prop.sel[v.change] <- 1 - prop.sel[v.change]
               prop.parm[v.change] <- 1 - prop.parm[v.change]}
          else if(prop.sel[v.change] == 1) 
               prop.parm[v.change] <- prop.sel[v.change] <- 0
          ### Priors
          prior.cur <- sum(dbern(cur.sel, p[which(selectable == 1)], log=TRUE), dbinom(sum(cur.sel), bin.n, bin.p, log=TRUE))
          prior.prop <- sum(dbern(prop.sel, p[which(selectable == 1)], log=TRUE), dbinom(sum(prop.sel), bin.n, bin.p, log=TRUE))
          ### Hit-And-Run Proposal Parameters
          theta <- rnorm(LIV)
          theta <- theta / sqrt(sum(theta*theta))
          lambda <- runif(1)
          ### Random-Scan Componentwise Estimation (Within-Model)
          for (j in sample(which(cur.parm == 1))) {
               ### Propose new parameter values
               temp.post <- Mo0[["parm"]]
               temp.post[which(temp.post == 0)] <- nonzero.post[which(temp.post == 0)]
               temp.post[which(cur.parm == 0)] <- 0
               prop <- Mo0[["parm"]] <- temp.post
               prop[j] <- prop[j] + lambda*theta[j]
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject (Within-Model Move)
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) Mo0 <- Mo1
               if(Mo0[["parm"]][j] != 0) nonzero.post[j] <- Mo0[["parm"]][j]
               Acceptance <- Acceptance + (u * (1 / sum(cur.parm)))}
          ### Random-Scan Componentwise Estimation (Between-Models)
          prop <- Mo0[["parm"]]
          prop[v.change] <- prop.sel[v.change]*(prop[v.change] +
               lambda*theta[v.change])
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          ### Accept/Reject (Between-Models Move)
          u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]] + prior.prop -
               prior.cur)
          if(u == TRUE) {
               Mo0 <- Mo1
               cur.sel <- prop.sel
               cur.parm <- prop.parm}
          if(Mo0[["parm"]][v.change] != 0)
               nonzero.post[v.change] <- Mo0[["parm"]][v.change]
          Acceptance <- Acceptance + (u * (1 / sum(prop.parm)))
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }
.mcmcrss <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, thinned, LogFile)
     {
     m <- Specs[["m"]]
     w <- Specs[["w"]]
     reflections <- 0
     Norm <- function(x) return(sqrt(sum(x*x)))
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          prop <- Mo0[["parm"]]
          g <- partial(Model, prop, Data)
          p <- rnorm(LIV)
          reflections <- 0
          ### Take m Steps
          for (i in 1:m) {
               prop <- prop + w*p
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               prop <- Mo1[["parm"]]
               ### Reflect at boundary
               if(Mo0[["LP"]] > Mo1[["LP"]]) {
                    reflections <- reflections + 1
                    p <- p - 2*g*{(t(p) %*% g) / Norm(g)^2}}}
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          prop <- Mo1[["parm"]]
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo1[["parm"]]
               Dev[t.iter] <- Mo1[["Dev"]]
               Mon[t.iter,] <- Mo1[["Monitor"]]}
          Mo0 <- Mo1
          }
     ### Output
     out <- list(Acceptance=Iterations, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }
.mcmcrwm <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning, VarCov, LogFile)
     {
     Block <- Specs[["B"]]
     if(length(Block) == 0) {
          U <- chol(VarCov)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Propose new parameter values
               prop <- as.vector(Mo0[["parm"]] +
                    rbind(rnorm(LIV)) %*% U)
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               log.u <- log(runif(1))
               log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
               if(!is.finite(log.alpha)) log.alpha <- 0
               if(log.u < log.alpha) {
                    Mo0 <- Mo1
                    Acceptance <- Acceptance + 1
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- Mo1[["parm"]]
                         Dev[t.iter] <- Mo1[["Dev"]]
                         Mon[t.iter,] <- Mo1[["Monitor"]]}
                    }
               }
          }
     else {
          B <- length(Block)
          if(!identical(length(VarCov), B))
               stop("Number of components in Covar differs from number of blocks.")
          for (b in 1:B) {
               if(!identical(length(Block[[b]]), length(diag(VarCov[[b]]))))
                    stop("Diagonal of Covar[[",b,"]] differs from block length.")}
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, sep="", file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0[["parm"]]
                    Dev[t.iter] <- Mo0[["Dev"]]
                    Mon[t.iter,] <- Mo0[["Monitor"]]}
               ### Proceed by Block
               for (b in 1:B) {
                    ### Propose new parameter values
                    prop <- Mo0[["parm"]]
                    prop[Block[[b]]] <- Mo0[["parm"]][[Block[[b]]]] +
                         rbind(rnorm(length(Block[[b]]))) %*%
                         chol(VarCov[[b]])
                    if({b == 1} & {iter %% Status == 0})
                         cat(",   Proposal: Multivariate,   LP:",  round(Mo0[["LP"]],1), "\n", sep="",  file=LogFile, append=TRUE)
                    ### Log-Posterior of the proposed state
                    Mo1 <- try(Model(prop, Data), silent=TRUE)
                    if(inherits(Mo1, "try-error")) Mo1 <- Mo0
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                         Mo1 <- Mo0
                    ### Accept/Reject
                    log.u <- log(runif(1))
                    log.alpha <- Mo1[["LP"]] - Mo0[["LP"]]
                    if(!is.finite(log.alpha)) log.alpha <- 0
                    if(log.u < log.alpha) {
                         Mo0 <- Mo1
                         Acceptance <- Acceptance + 1
                         if(iter %% Thinning == 0) {
                              thinned[t.iter,] <- Mo1[["parm"]]
                              Dev[t.iter] <- Mo1[["Dev"]]
                              Mon[t.iter,] <- Mo1[["Monitor"]]}
                         }
                    }
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=VarCov)
     return(out)
     }
.mcmcsamwg <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning, parm.names, LogFile)
     {
     Dyn <- Specs[["Dyn"]]
     Periodicity <- Specs[["Periodicity"]]
     Acceptance <- matrix(0, 1, LIV)
     for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
          Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
     Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
     staticparms <- c(1:LIV)[-as.vector(Dyn)]
     DiagCovar <- matrix(tuning, floor(Iterations/Periodicity), LIV, byrow=TRUE)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Select Order of Parameters
          if(length(staticparms) == 1) staticsample <- staticparms
          else staticsample <- sample(staticparms)
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          else dynsample <- as.vector(apply(Dyn, 1, sample))
          totsample <- c(staticsample, dynsample)
          ### Componentwise Estimation
          for (j in totsample) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) {
                    Mo0 <- Mo1
                    Acceptance[j] <- Acceptance[j] + 1}}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Adapt the Proposal Variance
          if(iter %% Periodicity == 0) {
               size <- 1 / min(100, sqrt(iter))
               Acceptance.Rate <- Acceptance / iter
               log.tuning <- log(tuning)
               tuning.num <- which(Acceptance.Rate > 0.44)
               log.tuning[tuning.num] <- log.tuning[tuning.num] + size
               log.tuning[-tuning.num] <- log.tuning[-tuning.num] - size
               tuning <- exp(log.tuning)
               a.iter <- floor(iter / Periodicity)
               DiagCovar[a.iter,] <- tuning}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)), Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=tuning)
     return(out)
     }
.mcmcsgld <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     epsilon <- Specs[["epsilon"]]
     file <- Specs[["file"]]
     Nr <- Specs[["Nr"]]
     Nc <- Specs[["Nc"]]
     size <- Specs[["size"]]
     Acceptance <- matrix(Iterations, 1, LIV)
     con <- file(file, open="r")
     on.exit(close(con))
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Sample Data
          seek(con, 0)
          skip.rows <- sample.int(Nr - size, size=1)
          Data[["X"]] <- matrix(scan(file=con, sep=",", skip=skip.rows, nlines=size, quiet=TRUE), size, Nc, byrow=TRUE)
          ### Propose new parameter values
          g <- partial(Model, Mo0[["parm"]], Data)
          eta <- rnorm(LIV, 0, epsilon[iter])
          prop <- Mo0[["parm"]] + {epsilon[iter]/2}*g + eta
          ### Log-Posterior of the proposed state
          Mo1 <- try(Model(prop, Data), silent=TRUE)
          if(inherits(Mo1, "try-error")) Mo1 <- Mo0
          else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
               Mo1 <- Mo0
          Mo0 <- Mo1
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }
.mcmcslice <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     m <- Specs[["m"]]
     w <- Specs[["w"]]
     uni.slice <- function(x0, g, j, w=1, m=Inf, lower=-Inf, upper=Inf, gx0=NULL, Data)
          {
          logy <- gx0[["LP"]] - rexp(1)
          ### Find the initial sampling interval
          u <- runif(1,0,w)
          L <- x0[j] - u
          R <- x0[j] + (w-u)
          ### Expand the interval until its ends are outside the slice, ### or until the limit on steps is reached
          intL <- intR <- x0
          intL[j] <- L
          intR[j] <- R
          ### Unlimited number of steps
          if(is.infinite(m)) { 
               repeat {
                    if(L <= lower) break
                    intL[j] <- L
                    gL <- g(intL, Data)
                    if(!is.finite(gL[["LP"]])) {
                         L <- L + w
                         break}
                    if(gL[["LP"]] <= logy) break
                    L <- L - w}
               repeat {
                    if(R >= upper) break
                    intR[j] <- R
                    gR <- g(intR, Data)
                    if(!is.finite(gR[["LP"]])) {
                         R <- R - w
                         break}
                    if(gR[["LP"]] <= logy) break
                    R <- R + w}
               }
          else if(m > 1) {
               ### Limited number of steps
               J <- floor(runif(1,0,m))
               K <- (m-1) - J
               while (J > 0) {
                    if(L <= lower) break
                    intL[j] <- L
                    gL <- g(intL, Data)
                    if(!is.finite(gL[["LP"]])) {
                         L <- L + w
                         break}
                    if(gL[["LP"]] <= logy) break
                    L <- L - w
                    J <- J - 1}
               while (K > 0) {
                    if(R >= upper) break
                    intR[j] <- R
                    gR <- g(intR, Data)
                    if(!is.finite(gR[["LP"]])) {
                         R <- R - w
                         break}
                    R <- R + w
                    K <- K - 1}
               }
          ### Shrink the interval to lower and upper bounds
          if(L < lower) L <- lower
          if(R > upper) R <- upper
          #### Sample from the interval, shrinking it on each rejection
          prop <- x0
          repeat { 
               x1 <- runif(1,L,R)
               prop[j] <- x1
               gx1 <- g(prop, Data)
               if(is.finite(gx1[["LP"]])) {
                    x1 <- gx1[["parm"]][j]
                    if(gx1[["LP"]] >= logy) break
                    if(x1 > x0[j]) R <- x1
                    else L <- x1}
               }
          return(gx1)
          }
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Random-Scan Componentwise Estimation
          for (j in sample.int(LIV)) {
               ### Univariate Slice Sampling
               Mo0 <- Mo1 <- uni.slice(x0=Mo0[["parm"]], g=Model, j=j, #w=1, m=Inf, lower=-Inf, upper=Inf, gx0=Mo0, Data)
                    w=w[j], m=m[j], lower=-Inf, upper=Inf, gx0=Mo0, Data)}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=Iterations, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }
.mcmcsmwg <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning, parm.names, LogFile)
     {
     Dyn <- Specs[["Dyn"]]
     Acceptance <- matrix(0, 1, LIV)
     for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
          Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
     Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
     staticparms <- c(1:LIV)[-as.vector(Dyn)]
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Select Order of Parameters
          if(length(staticparms) == 1) staticsample <- staticparms
          else staticsample <- sample(staticparms)
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          else dynsample <- as.vector(apply(Dyn, 1, sample))
          totsample <- c(staticsample, dynsample)
          ### Componentwise Estimation
          for (j in totsample) {
               ### Propose new parameter values
               prop <- Mo0[["parm"]]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) {
                    Mo0 <- Mo1
                    Acceptance[j] <- Acceptance[j] + 1}}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance)), Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=tuning)
     return(out)
     }
.mcmcthmc <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     epsilon <- Specs[["epsilon"]]
     L <- Specs[["L"]]
     Temperature <- Specs[["Temperature"]]
     gr <- partial(Model, Mo0[["parm"]], Data)
     sqrt.Temp <- sqrt(Temperature)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Propose new parameter values
          prop <- Mo0[["parm"]]
          momentum1 <- momentum0 <- rnorm(LIV)
          kinetic0 <- sum(momentum0^2) / 2
          Mo0.1 <- Mo0
          for (l in 1:L) {
               if(2*(l-1) < L) momentum1 <- momentum1 * sqrt.Temp
               else momentum1 <- momentum1 / sqrt.Temp
               momentum1 <- momentum1 + (epsilon/2) * gr
               prop <- prop + epsilon * momentum1
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0.1
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0.1
               if(any(Mo0.1[["parm"]] == Mo1[["parm"]])) {
                    nomove <- which(Mo0.1[["parm"]] == Mo1[["parm"]])
                    momentum1[nomove] <- -momentum1[nomove]
                    prop[nomove] <- prop[nomove] + momentum1[nomove]
                    Mo1 <- try(Model(prop, Data), silent=TRUE)
                    if(inherits(Mo1, "try-error")) Mo1 <- Mo0.1
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                         Mo1 <- Mo0.1}
               Mo0.1 <- Mo1
               prop <- Mo1[["parm"]]
               gr <- partial(Model, prop, Data)
               momentum1 <- momentum1 + (epsilon/2) * gr
               if(2*l > L) momentum1 <- momentum1 / sqrt.Temp
               else momentum1 <- momentum1 * sqrt.Temp}
          momentum1 <- -momentum1
          kinetic1 <- sum(momentum1^2) / 2
          ### Accept/Reject
          H0 <- -Mo0[["LP"]] + kinetic0
          H1 <- -Mo1[["LP"]] + kinetic1
          delta <- H1 - H0
          alpha <- min(1, exp(-delta))
          if(!is.finite(alpha)) alpha <- 0
          if(runif(1) < alpha) {
               Mo0 <- Mo1
               kinetic0 <- kinetic1
               Acceptance <- Acceptance + 1
               if(iter %% Thinning == 0) {
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]}
               }
          }
     ### Output
     out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=matrix(epsilon, 1, LIV), Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
     return(out)
     }
.mcmctwalk <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, LogFile)
     {
     xp0 <- SIV <- Specs[["SIV"]]
     n1 <- Specs[["n1"]]
     at <- Specs[["at"]]
     aw <- Specs[["aw"]]
     IntProd <- function(x) {return(sum(x*x))}
     DotProd <- function(x, y) {return(sum(x*y))}
     Simh1 <- function(dim, pphi, x, xp, beta)
          {
          phi <- runif(dim) < pphi
          rt <- NULL
          for (i in 1:dim)
               if(phi[i])
                    rt <- append(rt, xp[i] + beta*(xp[i] - x[i]))
               else
                    rt <- append(rt, x[i])
          return(list(rt=rt, nphi=sum(phi)))
          }
     Simfbeta <- function(at)
          {
          if(runif(1) < (at-1)/(2*at))
               return(exp(1/(at + 1)*log(runif(1))))
          else
               return(exp(1/(1 - at)*log(runif(1))))
          }
     Simh2 <- function(dim, pphi, aw, x, xp)
          {
          u <- runif(dim)
          phi <- runif(dim) < pphi
          z <- (aw/(1+aw))*(aw*u^2 + 2*u -1)
          z <- z*phi
          return(list(rt=x + (x - xp)*z, nphi=sum(phi)))
          }
     Simh3 <- function(dim, pphi, x, xp)
          {
          phi <- runif(dim) < pphi
          sigma <- max(phi*abs(xp - x))
          x + sigma*rnorm(dim)*phi
          return(list(rt=x + sigma*rnorm(dim)*phi, nphi=sum(phi), sigma=sigma))
          }
     G3U <- function(nphi, sigma, h, x, xp)
          {
          if(nphi > 0)
               return((nphi/2)*log(2*pi) + nphi*log(sigma) +
                    0.5*IntProd(h - xp)/(sigma^2))
          else
               return(0) 
          }
     Simh4 <- function(dim, pphi, x, xp)
          {
          phi <- runif(dim) < pphi
          sigma <- max(phi*abs(xp - x))/3
          rt <- NULL
          for (i in 1:dim)
               if(phi[i])
                    rt <- append(rt, xp[i] + sigma*rnorm(1))
               else
                    rt <- append(rt, x[i])
          return(list(rt=rt, nphi=sum(phi), sigma=sigma))
          }
     G4U <- function(nphi, sigma, h, x, xp)
          {
          if(nphi > 0)
               return((nphi/2)*log(2*pi) + nphi*log(sigma) +
                    0.5*IntProd((h - x))/(sigma^2))
          else
               return(0)
          }
     OneMove <- function(dim, Model, Data, x, U, xp, Up, at=at, aw=aw, pphi=pphi, F1=0.4918, F2=0.9836, F3=0.9918, Mo0.1, Mo0.2)
          {
          dir <- runif(1) ### Determine which set of points
          ker <- runif(1) ### Choose a kernel
          if(ker < F1) {
               ### Kernel h1: Traverse
               funh <- 1
               if(dir < 0.5) {
                    beta <- Simfbeta(at)
                    tmp <- Simh1(dim, pphi, xp, x, beta)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    y  <- x
                    propU <- U
                    Mo1.2 <- try(Model(yp, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.2, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.2[["LP"]]) &
                              identical(yp, as.vector(Mo1.2[["LP"]])))
                              check2 <- TRUE}
                    if(check1 & check2) {
                         propUp <- Mo1.2[["LP"]] * -1 ### Symmetric Proposal
                         if(nphi == 0)
                              A <- 1 ### Nothing moved
                         else
                              A <- exp((U - propU) + (Up - propUp) +
                                   (nphi-2)*log(beta))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               else {
                    beta <- Simfbeta(at)
                    tmp <- Simh1(dim, pphi, x, xp, beta)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    yp  <- xp
                    propUp <- Up
                    Mo1.1 <- try(Model(y, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.1, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.1[["LP"]]) &
                              identical(y, as.vector(Mo1.1[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2) {
                         propU <- Mo1.1[["LP"]] * -1 ### Symmetric Proposal
                         if(nphi == 0)
                              A <- 1 ### Nothing moved
                         else
                              A <- exp((U - propU) + (Up - propUp) +
                                   (nphi-2)*log(beta))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          else if(ker < F2) {
               ### Kernel h2: Walk
               funh <- 2
               if(dir < 0.5) {
                    ### x as pivot
                    tmp <- Simh2(dim, pphi, aw, xp, x)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    y  <- x
                    propU <- U
                    Mo1.2 <- try(Model(yp, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.2, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.2[["LP"]]) &
                              identical(yp, as.vector(Mo1.2[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(yp, y)) {
                         propUp <- Mo1.2[["LP"]] * -1
                         A <- exp((U - propU) + (Up - propUp))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               else {
                    ### xp as pivot
                    tmp <- Simh2(dim, pphi, aw, x, xp)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    yp  <- xp
                    propUp <- Up
                    Mo1.1 <- try(Model(y, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.1, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.1[["LP"]]) &
                              identical(y, as.vector(Mo1.1[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(yp, y)) {
                         propU <- Mo1.1[["LP"]] * -1
                         A <- exp((U - propU) + (Up - propUp))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          else if(ker < F3) {
               ### Kernel h3: Blow
               funh <- 3
               if(dir < 0.5) {
                    ### x as pivot
                    tmp <- Simh3(dim, pphi, xp, x)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    y  <- x
                    propU <- U
                    Mo1.2 <- try(Model(yp, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.2, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.2[["LP"]]) &
                              identical(yp, as.vector(Mo1.2[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(yp, x)) {
                         propUp <- Mo1.2[["LP"]] * -1
                         W1 <- G3U(nphi, sigma,  yp, xp,  x)
                         W2 <- G3U(nphi, sigma,  xp, yp,  x)
                         A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               else {
                    ### xp as pivot
                    tmp <- Simh3(dim, pphi, x, xp)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    yp  <- xp
                    propUp <- Up
                    Mo1.1 <- try(Model(y, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.1, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.1[["LP"]]) &
                              identical(y, as.vector(Mo1.1[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(y, xp)) {
                         propU <- Mo1.1[["LP"]] * -1
                         W1 <- G3U(nphi, sigma, y, x, xp)
                         W2 <- G3U(nphi, sigma, x, y, xp)
                         A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          else {
               ## Kernel h4: Hop
               funh <- 4
               if(dir < 0.5) {
                    ### x as pivot
                    tmp <- Simh4(dim, pphi, xp, x)
                    yp <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    y  <- x
                    propU <- U
                    Mo1.2 <- try(Model(yp, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.2, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.2[["LP"]]) &
                              identical(yp, as.vector(Mo1.2[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(yp, x)) {
                         propUp <- Mo1.2[["LP"]] * -1
                         W1 <- G4U(nphi, sigma, yp, xp, x)
                         W2 <- G4U(nphi, sigma, xp, yp, x)
                         A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
                    else {
                         propUp <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               else {
                    ### xp as pivot
                    tmp <- Simh4(dim, pphi, x, xp)
                    y <- tmp$rt
                    nphi <- tmp$nphi
                    sigma <- tmp$sigma
                    yp  <- xp
                    propUp <- Up
                    Mo1.1 <- try(Model(y, Data), silent=TRUE)
                    check1 <- check2 <- FALSE
                    if(!inherits(Mo1.1, "try-error")) {
                         check1 <- TRUE
                         if(is.finite(Mo1.1[["LP"]]) &
                              identical(y, as.vector(Mo1.1[["parm"]])))
                              check2 <- TRUE}
                    if(check1 & check2 & !identical(y, xp)) {
                         propU <- Mo1.1[["LP"]] * -1
                         W1 <- G4U(nphi, sigma, y, x, xp)
                         W2 <- G4U(nphi, sigma, x, y, xp)
                         A <- exp((U - propU) + (Up - propUp) + (W1 - W2))}
                    else {
                         propU <- NULL
                         A <- 0  ### Out of support, not accepted
                         }
                    }
               }
          if(check1 & check2 & is.finite(A) & (dir < 0.5))
               Mo0.2 <- Mo1.2
          else if(check1 & check2 & is.finite(A) & (dir >= 0.5))
               Mo0.1 <- Mo1.1
          else if(!is.finite(A)) A <- 0
          return(list(y=y, propU=propU, yp=yp, propUp=propUp, A=A, funh=funh, nphi=nphi, Mo0.1=Mo0.1, Mo0.2=Mo0.2))
          }
     Runtwalk <- function(Iterations, dim, x0, xp0, pphi, at, aw, F1=0.4918, F2=F1+0.4918, F3=F2+0.0082, Model, Data, Status, Thinning, Acceptance, Dev, Mon, Mo0, thinned, LogFile)
          {
          x <- x0 ### Primary vector of initial values
          xp <- xp0 ### Secondary vector of initial values
          Mo0.1 <- try(Model(x, Data), silent=TRUE)
          Mo0.2 <- try(Model(xp, Data), silent=TRUE)
          if(inherits(Mo0.1, "try-error") | inherits(Mo0.2, "try-error"))
               stop("Error in estimating the log-posterior.", file=LogFile, append=TRUE)
          if(any(!is.finite(c(Mo0.1[["LP"]], Mo0.2[["LP"]]))))
               stop("The log-posterior is non-finite.", file=LogFile, append=TRUE)
          if(identical(x, as.vector(Mo0.1[["parm"]])) &
               identical(xp, as.vector(Mo0.2[["parm"]]))) {
               U <- Mo0.1[["LP"]] * -1
               Up <- Mo0.2[["LP"]] * -1}
          else {
               cat("\nInitial values are out of support.", file=LogFile, append=TRUE)
               cat("\n  Initial.Values=", x, file=LogFile, append=TRUE)
               cat("\n SIV=", xp, file=LogFile, append=TRUE)
               stop("Try re-specifying initial values.", file=LogFile, append=TRUE)}
          if(any(abs(x - xp) <= 0))
               stop("\nBoth vectors of initial values are not unique.", file=LogFile, append=TRUE)
          Acceptance <- 0
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, ",   Proposal: Multivariate Subset,   LP:", round(Mo0.1[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo0.1[["parm"]]
                    Dev[t.iter] <- Mo0.1[["Dev"]]
                    Mon[t.iter,] <- Mo0.1[["Monitor"]]}
               ### Assign x and xp
               x <- as.vector(Mo0.1[["parm"]])
               xp <- as.vector(Mo0.2[["parm"]])
               ### Propose New Parameter Values
               move <- OneMove(dim=dim, Model, Data, x, U, xp, Up, at=at, aw=aw, pphi=pphi, F1=F1, F2=F2, F3=F3, Mo0.1=Mo0.1, Mo0.2=Mo0.2)
               ### Accept/Reject
               if(runif(1) < move$A) {
                    Mo0.1 <- move$Mo0.1
                    Mo0.2 <- move$Mo0.2
                    Acceptance <- Acceptance + 1 #move$nphi/dim
                    if(iter %% Thinning == 0) {
                         thinned[t.iter,] <- move$Mo0.1[["parm"]]
                         Dev[t.iter] <- move$Mo0.1[["Dev"]]
                         Mon[t.iter,] <- move$Mo0.1[["Monitor"]]}
                    x <- move$y
                    U <- move$propU
                    xp <- move$yp
                    Up <- move$propUp
                    }
               }
          out <- list(Acceptance=Acceptance, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=.colVars(thinned))
          return(out)
          }
     out <- Runtwalk(Iterations=Iterations, dim=LIV, x0=Mo0[["parm"]], xp0=xp0, pphi=min(LIV, n1)/LIV, at=6, aw=1.5, Model=Model, Data=Data, Status=Status, Thinning=Thinning, Acceptance=Acceptance, Dev=Dev, Mon=Mon, Mo0=Mo0, thinned=thinned, LogFile=LogFile)
     ### Output
     return(out)
     }
.mcmcuess <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, VarCov, LogFile)
     {
     A <- Specs[["A"]]
     Block <- Specs[["B"]]
     m <- Specs[["m"]]
     n <- Specs[["n"]]
     w <- 0.05
     B <- length(Block)
     if(B == 0) {
          decomp.freq <- max(LIV * floor(Iterations / Thinning / 100), 10)
          S.eig <-try(eigen(VarCov), silent=TRUE)
          if(inherits(S.eig, "try-error")) S.eig <- NULL
          obs.sum <- matrix(Mo0[["parm"]]*n, LIV, 1)
          obs.scatter <- tcrossprod(Mo0[["parm"]])*n
          DiagCovar <- matrix(0, floor(Iterations / Thinning), LIV)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               ### Eigenvectors of the Sample Covariance Matrix
               if({iter %% decomp.freq == 0} & {iter > 1} & {iter < A}) {
                    VarCov <- obs.scatter/{n + iter} -
                         tcrossprod(obs.sum/{n + iter})
                    S.eig <- eigen(VarCov)}
               ### Non-Adaptive or Adaptive
               if(runif(1) < w || is.null(S.eig)) {
                    v <- rnorm(LIV)
                    v <- v / sqrt(sum(v*v))
                    }
               else {
                    which.eig <- floor(1 + LIV * runif(1))
                    v <- S.eig$vectors[,which.eig] *
                         sqrt(abs(S.eig$values[which.eig]))}
               ### Slice Interval
               Mo0.1 <- try(Model(Mo0[["parm"]], Data), silent=TRUE)
               if(inherits(Mo0.1, "try-error")) Mo0.1 <- Mo0
               Mo0 <- Mo0.1
               y.slice <- Mo0[["LP"]] - rexp(1)
               L <- -runif(1)
               U <- L + 1
               if(m > 0) {
                    L.y <- try(Model(Mo0[["parm"]] + v*L, Data)[["LP"]], silent=TRUE)
                    if(inherits(L.y, "try-error")) L.y <- Mo0[["LP"]]
                    else if(!is.finite(L.y)) L.y <- Mo0[["LP"]]
                    U.y <- try(Model(Mo0[["parm"]] + v*U, Data)[["LP"]], silent=TRUE)
                    if(inherits(U.y, "try-error")) U.y <- Mo0[["LP"]]
                    else if(!is.finite(U.y)) U.y <- Mo0[["LP"]]
                    step <- 0
                    while({L.y > y.slice || U.y > y.slice} && step < m) {
                         step <- step + 1
                         if(runif(1) < 0.5) {
                              L <- L - 1
                              L.y <- try(Model(Mo0[["parm"]] + v*L, Data)[["LP"]],  silent=TRUE)
                              if(inherits(L.y, "try-error")) L.y <- Mo0[["LP"]]
                              else if(!is.finite(L.y)) L.y <- Mo0[["LP"]]
                              }
                         else {
                              U <- U + 1
                              U.y <- try(Model(Mo0[["parm"]] + v*U, Data)[["LP"]],  silent=TRUE)
                              if(inherits(U.y, "try-error")) U.y <- Mo0[["LP"]]
                              else if(!is.finite(U.y)) U.y <- Mo0[["LP"]]}}}
               ### Rejection Sampling
               repeat {
                    prop.offset <- runif(1, min=L, max=U)
                    prop <- Mo0[["parm"]] + prop.offset * v
                    Mo1 <- try(Model(prop, Data), silent=TRUE)
                    if(inherits(Mo1, "try-error")) Mo1 <- Mo0
                    else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                         Mo1 <- Mo0
                    prop <- Mo1[["parm"]]
                    if(Mo1[["LP"]] >= y.slice) break
                    else if(abs(prop.offset < 1e-100)) {
                         Mo1 <- Mo0
                         break}
                    if(prop.offset < 0) L <- prop.offset
                    else U <- prop.offset}
               ### Save Thinned Samples
               if(iter %% Thinning == 0) {
                    t.iter <- floor(iter / Thinning) + 1
                    thinned[t.iter,] <- Mo1[["parm"]]
                    Dev[t.iter] <- Mo1[["Dev"]]
                    Mon[t.iter,] <- Mo1[["Monitor"]]
                    DiagCovar[t.iter-1,] <- diag(S.eig$vectors)}
               obs.sum <- obs.sum + Mo1[["parm"]]
               obs.scatter <- obs.scatter + tcrossprod(Mo1[["parm"]])
               Mo0 <- Mo1}
          }
     else {
          S.eig <- obs.sum <- obs.scatter <- list()
          decomp.freq <- rep(0, length(B))
          for (b in 1:B) {
               decomp.freq[b] <- max(length(Block[[b]]) *
                    floor(Iterations / Thinning / 100), 10)
               S.eig[[b]] <-try(eigen(VarCov[[b]]), silent=TRUE)
               if(inherits(S.eig[[b]], "try-error")) S.eig[[b]] <- NULL
               obs.sum[[b]] <- matrix(Mo0[["parm"]][Block[[b]]]*n, length(Block[[b]]), 1)
               obs.scatter[[b]] <- tcrossprod(Mo0[["parm"]][Block[[b]]])*n
          }
          DiagCovar <- matrix(0, floor(Iterations / Thinning), LIV)
          for (iter in 1:Iterations) {
               ### Print Status
               if(iter %% Status == 0)
                    cat("Iteration: ", iter, ",   Proposal: Multivariate,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
               ### Proceed by Block
               for (b in 1:B) {
                    ### Eigenvectors of the Sample Covariance Matrix
                    if({iter %% decomp.freq[b] == 0} & {iter > 1} &
                         {iter < A}) {
                         VarCov[[b]] <- obs.scatter[[b]]/{n + iter} -
                              tcrossprod(obs.sum[[b]]/{n + iter})
                         S.eig[[b]] <- eigen(VarCov[[b]])}
                    ### Non-Adaptive or Adaptive
                    if(runif(1) < w || is.null(S.eig[[b]])) {
                         v <- rnorm(length(Block[[b]]))
                         v <- v / sqrt(sum(v*v))
                         }
                    else {
                         which.eig <- floor(1 + length(Block[[b]]) * runif(1))
                         v <- S.eig[[b]]$vectors[,which.eig] *
                              sqrt(abs(S.eig[[b]]$values[which.eig]))}
                    ### Slice Interval
                    Mo0.1 <- try(Model(Mo0[["parm"]], Data), silent=TRUE)
                    if(inherits(Mo0.1, "try-error")) Mo0.1 <- Mo0
                    Mo0 <- Mo0.1
                    y.slice <- Mo0[["LP"]] - rexp(1)
                    L <- -runif(1)
                    U <- L + 1
                    if(m > 0) {
                         prop <- Mo0[["parm"]]
                         prop[Block[[b]]] <- prop[Block[[b]]] + v*L
                         L.y <- try(Model(prop, Data)[["LP"]], silent=TRUE)
                         if(inherits(L.y, "try-error")) L.y <- Mo0[["LP"]]
                         else if(!is.finite(L.y)) L.y <- Mo0[["LP"]]
                         prop <- Mo0[["parm"]]
                         prop[Block[[b]]] <- prop[Block[[b]]] + v*U
                         U.y <- try(Model(prop, Data)[["LP"]], silent=TRUE)
                         if(inherits(U.y, "try-error")) U.y <- Mo0[["LP"]]
                         else if(!is.finite(U.y)) U.y <- Mo0[["LP"]]
                         step <- 0
                         while({L.y > y.slice || U.y > y.slice} && step < m) {
                              step <- step + 1
                              if(runif(1) < 0.5) {
                                   L <- L - 1
                                   prop <- Mo0[["parm"]]
                                   prop[Block[[b]]] <- prop[Block[[b]]] + v*L
                                   L.y <- try(Model(prop, Data)[["LP"]],  silent=TRUE)
                                   if(inherits(L.y, "try-error")) L.y <- Mo0[["LP"]]
                                   else if(!is.finite(L.y)) L.y <- Mo0[["LP"]]
                                   }
                              else {
                                   U <- U + 1
                                   prop <- Mo0[["parm"]]
                                   prop[Block[[b]]] <- prop[Block[[b]]] + v*U
                                   U.y <- try(Model(prop, Data)[["LP"]],  silent=TRUE)
                                   if(inherits(U.y, "try-error")) U.y <- Mo0[["LP"]]
                                   else if(!is.finite(U.y)) U.y <- Mo0[["LP"]]}}}
                    ### Rejection Sampling
                    repeat {
                         prop.offset <- runif(1, min=L, max=U)
                         prop <- Mo0[["parm"]]
                         prop[Block[[b]]] <- prop[Block[[b]]] + prop.offset*v
                         Mo1 <- try(Model(prop, Data), silent=TRUE)
                         if(inherits(Mo1, "try-error")) Mo1 <- Mo0
                         else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]],  Mo1[["Monitor"]]))))
                              Mo1 <- Mo0
                         prop <- Mo1[["parm"]]
                         if(Mo1[["LP"]] >= y.slice) break
                         else if(abs(prop.offset < 1e-100)) {
                              Mo1 <- Mo0
                              break}
                         if(prop.offset < 0) L <- prop.offset
                         else U <- prop.offset}
                    ### Save Thinned Samples
                    if(iter %% Thinning == 0) {
                         t.iter <- floor(iter / Thinning) + 1
                         thinned[t.iter,] <- Mo1[["parm"]]
                         Dev[t.iter] <- Mo1[["Dev"]]
                         Mon[t.iter,] <- Mo1[["Monitor"]]
                         DiagCovar[t.iter-1,Block[[b]]] <- diag(S.eig[[b]]$vectors)}
                    obs.sum[[b]] <- obs.sum[[b]] + Mo1[["parm"]][Block[[b]]]
                    obs.scatter[[b]] <- obs.scatter[[b]] +
                         tcrossprod(Mo1[["parm"]][Block[[b]]])
                    Mo0 <- Mo1}
               }
          }
     ### Output
     out <- list(Acceptance=Iterations, Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=VarCov)
     return(out)
     }
.mcmcusamwg <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning, parm.names, LogFile)
     {
     Dyn <- Specs[["Dyn"]]
     Periodicity <- Specs[["Periodicity"]]
     Fit <- Specs[["Fit"]]
     Begin <- Specs[["Begin"]]
     Acceptance <- matrix(0, 1, LIV)
     for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
          Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
     Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
     Dyn <- matrix(Dyn[-c(1:(Begin-1)),], nrow(Dyn)-Begin+1, ncol(Dyn))
     n.samples <- nrow(Fit$Posterior1)
     mults <- Iterations / n.samples
     samps <- rep(1:n.samples, each=mults)
     if(Iterations != length(samps))
          stop("Iterations not a multiple of posterior samples.", file=LogFile, append=TRUE)
     ivs <- Mo0[["parm"]]
     post <- Fit$Posterior1[samps,]
     post[1,as.vector(Dyn)] <- ivs[as.vector(Dyn)]
     DiagCovar <- matrix(tuning, floor(Iterations/Periodicity), LIV, byrow=TRUE)
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Store Current Posterior
          if(iter > 1) post[iter,as.vector(Dyn)] <- post[iter-1,as.vector(Dyn)]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- post[iter,]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Select Order of Parameters
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          else dynsample <- as.vector(apply(Dyn, 1, sample))
          ### Componentwise Estimation
          for (j in dynsample) {
               ### Propose new parameter values
               prop <- post[iter,]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) {
                    Mo0 <- Mo1
                    post[iter,] <- Mo0[["parm"]]
                    Acceptance[j] <- Acceptance[j] + 1}}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Adapt the Proposal Variance
          if(iter %% Periodicity == 0) {
               size <- 1 / min(100, sqrt(iter))
               Acceptance.Rate <- Acceptance / iter
               log.tuning <- log(tuning)
               tuning.num <- which(Acceptance.Rate > 0.44)
               log.tuning[tuning.num] <- log.tuning[tuning.num] + size
               log.tuning[-tuning.num] <- log.tuning[-tuning.num] - size
               tuning <- exp(log.tuning)
               a.iter <- floor(iter / Periodicity)
               DiagCovar[a.iter,] <- tuning}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance[dynsample])), Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=tuning)
     return(out)
     }
.mcmcusmwg <- function(Model, Data, Iterations, Status, Thinning, Specs, Acceptance, Dev, DiagCovar, LIV, Mon, Mo0, ScaleF, thinned, tuning, parm.names, LogFile)
     {
     Dyn <- Specs[["Dyn"]]
     Fit <- Specs[["Fit"]]
     Begin <- Specs[["Begin"]]
     Acceptance <- matrix(0, 1, LIV)
     for (k in 1:ncol(Dyn)) {for (t in 1:nrow(Dyn)) {
          Dyn[t,k] <- which(parm.names == Dyn[t,k])}}
     Dyn <- matrix(as.numeric(Dyn), nrow(Dyn), ncol(Dyn))
     Dyn <- matrix(Dyn[-c(1:(Begin-1)),], nrow(Dyn)-Begin+1, ncol(Dyn))
     n.samples <- nrow(Fit$Posterior1)
     mults <- Iterations / n.samples
     samps <- rep(1:n.samples, each=mults)
     if(Iterations != length(samps))
          stop("Iterations not a multiple of posterior samples.", file=LogFile, append=TRUE)
     ivs <- Mo0[["parm"]]
     post <- Fit$Posterior1[samps,]
     post[1,as.vector(Dyn)] <- ivs[as.vector(Dyn)]
     for (iter in 1:Iterations) {
          ### Print Status
          if(iter %% Status == 0)
               cat("Iteration: ", iter, ",   Proposal: Componentwise,   LP:", round(Mo0[["LP"]],1), "\n", sep="", file=LogFile, append=TRUE)
          ### Store Current Posterior
          if(iter > 1) post[iter,as.vector(Dyn)] <- post[iter-1,as.vector(Dyn)]
          ### Save Thinned Samples
          if(iter %% Thinning == 0) {
               t.iter <- floor(iter / Thinning) + 1
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          ### Select Order of Parameters
          if(ncol(Dyn) == 1) dynsample <- sample(Dyn)
          else dynsample <- as.vector(apply(Dyn, 1, sample))
          ### Componentwise Estimation
          for (j in dynsample) {
               ### Propose new parameter values
               prop <- post[iter,]
               prop[j] <- rnorm(1, prop[j], tuning[j])
               ### Log-Posterior of the proposed state
               Mo1 <- try(Model(prop, Data), silent=TRUE)
               if(inherits(Mo1, "try-error")) Mo1 <- Mo0
               else if(any(!is.finite(c(Mo1[["LP"]], Mo1[["Dev"]], Mo1[["Monitor"]]))))
                    Mo1 <- Mo0
               ### Accept/Reject
               u <- log(runif(1)) < (Mo1[["LP"]] - Mo0[["LP"]])
               if(u == TRUE) {
                    Mo0 <- Mo1
                    post[iter,] <- Mo0[["parm"]]
                    Acceptance[j] <- Acceptance[j] + 1}}
          if(iter %% Thinning == 0) {
               thinned[t.iter,] <- Mo0[["parm"]]
               Dev[t.iter] <- Mo0[["Dev"]]
               Mon[t.iter,] <- Mo0[["Monitor"]]}
          }
     ### Output
     out <- list(Acceptance=mean(as.vector(Acceptance[dynsample])), Dev=Dev, DiagCovar=DiagCovar, Mon=Mon, thinned=thinned, VarCov=tuning)
     return(out)
     }

#End

}
