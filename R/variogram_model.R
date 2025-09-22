
variogram_model = function( 
  vgm, method="optim", plotdata=FALSE,  
  weight_by_inverse_distance=TRUE, autocorrelation_localrange=0.1,
  trim =0.75 ) {
   
  # stand alone version
  
  if (0) {

    out = variogram_decorrelated( 
      xy = xyz[, c("x", "y")], 
      z=xyz$z, 
      nbreaks=30
    )
    vgm = out$vgm
    trim = 0.75
  } 

  retain = which( out$vgm$distances < max(out$vgm$distances, na.rm=TRUE)*trim ) 
  if (exists("keep", vgm)) retain = unique( c(retain, which( vgm$keep )) )

  out = list(
    vgm = vgm,
    retain=retain,
    vgm_dist_max = max(vgm$distances[retain]),
    vgm_var_max = max(vgm$sv[retain]),
    summary = list()
  )
 
  vxs = out$vgm$distances[retain] / out$vgm_dist_max
  vgs = out$vgm$sv[retain] / out$vgm_var_max

  if (method == "optim" ) {

    #\\ simple nonlinear least squares fit
    vario_function_phi_nu = function(par, vgs, vxs, w){
      # ie. nu and phi are both estimated
      if (par["nu"] < 0.01) par["nu"] = 0.01
      if (par["phi"] <= 0.001 ) par["phi"] = 0.001
      if (par["nu"] > 5 ) par["nu"] = 5
      if (par["phi"] > 5 ) par["phi"] = 5
      vgm = par["total.var"] *( (1 -par["sigma.sq.fraction"])  +  par["sigma.sq.fraction"]*( 1-stmv_matern(distance=vxs, mRange=par["phi"], mSmooth=par["nu"]) ) )
      obj = sum( w * (vgs - vgm)^2, na.rm=TRUE) # vario normal errors, no weights , etc.. just the line
      return(obj)
    }

    if (weight_by_inverse_distance) {
      w = 1/vxs
    } else {
      w = rep(1, length(vxs ) )
    }

    lower = c(0.1, 0, 0.001, 0.01 ) 
    upper = c(1.9, 1, 10,    10)

    fit = try( optim( 
      par=c(total.var=1, sigma.sq.fraction=0.5, phi=1.1, nu=0.9), 
      vgs=vgs, 
      vxs=vxs, 
      w=w, 
      method="L-BFGS-B", 
      lower=lower, 
      upper=upper, 
      fn = vario_function_phi_nu, 
      control=list(factr=1e-9, maxit = 1000L) 
    ))

    if ( !inherits(fit, "try-error")) if ( fit$convergence != 0 ) class(fit) = "try-error"
    if (exists( "par", fit)) {
      if ( fit$par[["phi"]] < lower[3] | fit$par[["phi"]] > upper[3] ) class(fit) = "try-error"    # give up
      if ( fit$par[["nu"]]  < lower[4] | fit$par[["nu"]]  > upper[4] ) class(fit) = "try-error"    # give up
    }
    
    if ( !inherits(fit, "try-error")) {

      out$summary$convergence = fit$convergence 

      if ( out$summary$convergence==0 ) {
          out$fit = fit
          out$summary$nu =  fit$par[["nu"]] 
          out$summary$phi = fit$par[["phi"]] * vgm_dist_max  
          out$summary$localrange = matern_phi2distance( phi=out$summary$phi, nu=out$summary$nu, cor=autocorrelation_localrange )
          out$summary$varSpatial = fit$par[["total.var"]] * fit$par[["sigma.sq.fraction"]]  * vgm_var_max
          out$summary$varObs = fit$par[["total.var"]] * (1- fit$par[["sigma.sq.fraction"]]) * vgm_var_max
          out$summary$phi_ok = ifelse( out$summary$phi < out$summary$vgm_dist_max*0.99 , TRUE, FALSE )
      }
    }
    
    return(out)
  }


  if (method == "julia_turing") {
   require( JuliaConnectoR )

    # simple GP:
    # diffuse priors centered over expected values
    # phi scaled to distmax: so expected to be in interval (0, 1) .. range from 0.01 to 5 
    # nu expected to be from 0.5 to 5 .. centered on 1 
    # vtot scaled to max of 1 
    # varFraction 

    # JuliaDB not installing due to dependency issue 
    # JuliaDB = juliaImport("JuliaDB")
    # u = JuliaDB$table( vxs )  # send vxs to julia
    # u

      # Turing.jl
      # MCMCChains.jl
      # DynamicPPL.jl
      # AdvancedHMC.jl
      # DistributionsAD.jl
      # Bijectors.jl

    matern_fit = juliaEval('
      # using Zygote, ReverseDiff, ForwardDiff, Tracker , ReverseDiff, DynamicHMC, 
      # using Optim, MCMCChains, DynamicPPL,  DynamicHMC, AdvancedHMC, DistributionsAD, Bijectors, MCMCChains, 
      using SpecialFunctions, Distributions, Statistics
      using Zygote, Tracker, DynamicHMC, AdvancedHMC, DistributionsAD, ForwardDiff, ReverseDiff, Turing
      using Random
      Random.seed!(1);
      Turing.setadbackend(:reversediff)
      # Turing.setadbackend(:forwarddiff)
      # Turing.setadbackend(:zygote)
      # Turing.setadbackend(:tracker)
      
      function matern( nu, hs )
        2.0^(1.0 - nu) / gamma(nu)  * hs^nu * besselk(nu, hs)
      end

      @model function variogram_as_matern(vxs, vgs)
        n = length(vgs)
        eps = 1.0e-12 # to avoid zeros
        nu  ~ truncated( Cauchy( 1.0, 0.5), 0.01, 10.0 )
        phi ~ truncated( Cauchy( 0.0, 0.1), eps, 1.0 )
        vtot ~ truncated( Cauchy( 0.0, 0.5), eps, 10.0 )
        vfsp ~ Beta( 2, 2 )
        yvar ~ truncated( Cauchy( 0.0, 0.5), eps, 10.0 )
        for i in 1:n
          vg = 1.0 - vfsp * matern( nu, sqrt(2.0*nu) * (vxs[i] + eps) / phi )
          if ( vg < 0 ) 
            Turing.@addlogprob! -Inf  
            return
          end
          vgs[i] ~ Normal( vtot * vg, sqrt(yvar) )
        end
      end
    ')

    #  Run sampler, collect results.

      # mle_estimate = optimize(variogram_as_matern(vxs, vgs), MLE())
      # mle_estimate = optimize(variogram_as_matern(vxs, vgs), MLE(), NelderMead())
      # mle_estimate = optimize(variogram_as_matern(vxs, vgs), MLE(), SimulatedAnnealing())
      # mle_estimate = optimize(variogram_as_matern(vxs, vgs), MLE(), ParticleSwarm())
      # mle_estimate = optimize(variogram_as_matern(vxs, vgs), MLE(), Newton())
      # mle_estimate = optimize(variogram_as_matern(vxs, vgs), MLE(), AcceleratedGradientDescent())
      # mle_estimate = optimize(variogram_as_matern(vxs, vgs), MLE(), Newton(), Optim.Options(iterations=10_000, allow_f_increases=true))

      # using StatsBase
      # coeftable(mle_estimate)

      # map_estimate = optimize(variogram_as_matern(vxs, vgs), MAP())

      # c1 = sample(variogram_as_matern(vxs, vgs), SMC(), 1000)
      # c2 = sample(variogram_as_matern(vxs, vgs), PG(10), 1000)
      # c3 = sample(variogram_as_matern(vxs, vgs), HMC(0.1, 5), 1000)
      # c4 = sample(variogram_as_matern(vxs, vgs), Gibbs(PG(10, :m), HMC(0.1, 5, :s²)), 1000)
      # c5 = sample(variogram_as_matern(vxs, vgs), HMCDA(0.15, 0.65), 1000)
      # c6 = sample(variogram_as_matern(vxs, vgs), NUTS(0.65), 1000)

      # c0 = sample(variogram_as_matern(vxs, vgs), NUTS(0.65), 1000, init_params = map_estimate.values.array)
      # SMC: number of particles.
      # PG: number of particles, number of iterations.
      # HMC: leapfrog step size, leapfrog step numbers.
      # Gibbs: component sampler 1, component sampler 2, ...
      # HMCDA: total leapfrog length, target accept ratio.
      # NUTS: number of adaptation steps (optional), target accept ratio.
    
      #  # Summarise results
      #   describe(chn)

      #   # Plot and save results
      #   p = plot(chn)
      #   chn
      
    warmup = 1000L
    niters = 4000L

    Turing = juliaImport("Turing") 

    # chain = Turing$sample( matern_fit(vxs, vgs), Turing$PG(10L), niters, n_adapts=warmup, progress=FALSE, verbose=FALSE, drop_warmup=TRUE )
    
    
    chain = Turing$sample( matern_fit(vxs, vgs), Turing$SMC(), niters, n_adapts=warmup, progress=FALSE, verbose=FALSE, drop_warmup=TRUE )
    
    # chain = Turing$sample( matern_fit(vxs, vgs), Turing$NUTS(0.95), niters, n_adapts=warmup, progress=FALSE, verbose=FALSE, drop_warmup=TRUE )
    
    Turing$summarize( chain)

    # graphical backed (GR) not working with R 
    # JuliaDB = juliaImport("c")
    # JuliaDB$plot( sss )
    
    res = JuliaConnectoR::juliaGet(chain)
    o = as.data.frame(res$value$data[(warmup+1):niters, ,] )
    names( o ) = as.character(unlist(res$value$axes[[2]][["val"]]))
    

    if ( out$summary$convergence==0 ) {
        out$samples = o
        out$summary$nu =  mean(o[,"nu"], na.rm=TRUE)
        out$summary$phi = mean(o[,"phi"], na.rm=TRUE) * out$vgm_dist_max  
        out$summary$localrange = matern_phi2distance( phi=out$summary$phi, nu=out$summary$nu, cor=autocorrelation_localrange )
        out$summary$varSpatial = mean(o[,"vtot"], na.rm=TRUE) * mean(o[,"vfsp"], na.rm=TRUE) * out$vgm_var_max
        out$summary$varObs = mean(o[,"vtot"], na.rm=TRUE) * (1- mean(o[,"vfsp"], na.rm=TRUE) ) * out$vgm_var_max
        out$summary$phi_ok = ifelse( out$summary$phi < out$vgm_dist_max*0.99 , TRUE, FALSE )
    }
    return(out)
  
  }
 

  # to plot:
  # add solution to gstat solution
    xlim= c(0, out$vgm_dist_max*1.1 )
    ylim= c(0, out$vgm_var_max*1.1)
    plot( out$vgm$distances, out$vgm$sv,  xlim=xlim, ylim=ylim, type="n" )
    points( out$vgm$distances, out$vgm$sv, col="black", pch=19  )
    ds = seq( 0, out$vgm_dist_max  , length.out=100 )

    ac = out$summary$varObs + out$summary$varSpatial*(1 - stmv_matern( ds, out$summary$phi, out$summary$nu ) )
    lines( ds, ac, col="red", lwd=4)
    abline( h=0, lwd=1, col="lightgrey" )
    abline( v=0 ,lwd=1, col="lightgrey" )
    abline( h=out$summary$varObs, lty="dashed", col="grey" )
    abline( h=out$summary$varObs + out$summary$varSpatial, lty="dashed", col="grey" )
    abline( v=out$summary$localrange, lty="dashed", col="grey")
 


  if (0) {
 
    install.packages("JuliaConnectoR")
    # juliaEval('using Pkg; Pkg.add(PackageSpec(name = "Flux", version = "0.12"))')

    loadfunctions( c( "aegis", "stmv" ))

    xyz = stmv_test_data( datasource="meuse" )  
    xyz$log_zinc = log( xyz$zinc )
    xyz$z = residuals( lm( log_zinc ~ 1, xyz ) )

    # first a gstat solution:
    library(sp)
    library(gstat)
    xyz2 = xyz 
    coordinates(xyz2) = ~x+y
    vgm1 <- variogram( z~1, xyz2)
    # optimize the value of kappa in a Matern model, using ugly <<- side effect:
    f = function(x) attr(m.fit <<- fit.variogram(vgm1, vgm(,"Mat",nugget=NA,kappa=x)),"SSErr")
    optimize(f, c(0.1, 5))
    plot(vgm1, m.fit)
    # best fit from the (0.3, 0.4, 0.5. ... , 5) sequence:
    (m <- fit.variogram(vgm1, vgm("Mat"), fit.kappa = TRUE))
    attr(m, "SSErr")
  
  }

}
