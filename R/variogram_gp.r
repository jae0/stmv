 
 variogram_gp = function( xy, z, nx=NULL, ny=NULL,   
    plotdata=FALSE, eps=1e-9, autocorrelation_localrange=0.1, 
    internal_scale=NULL, fit_model=TRUE,
    weight_by_inverse_distance=TRUE 
  ) {
 
  # using covariance via Turing

  if (0) {

    loadfunctions( c( "aegis", "stmv" ))

    xyz = stmv_test_data( datasource="meuse" )  
    xyz$log_zinc = log( xyz$zinc )
    xyz$z = residuals( lm( log_zinc ~ 1, xyz ) )
 
      xy=xyz[, c("x", "y")]
      z=xyz$z

    internal_scale = NULL

    autocorrelation_localrange=0.1

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
    # # best fit from the (0.3, 0.4, 0.5. ... , 5) sequence:
    # (m <- fit.variogram(vgm1, vgm("Mat"), fit.kappa = TRUE))
    # attr(m, "SSErr")

    # size_diagonal = sqrt( diff(range(xyz$x))^2 + diff(range(xyz$y))^2) # crude scale of problem (domain size)
    out = variogram_decorrelated( 
      xy=xyz[, c("x", "y")], 
      z=xyz$z 
    )
    # add solution to gstat solution
    xlim= c(0, max(out$vgm$distances) * 1.1 )
    ylim= c(0, max(out$vgm$sv)*1.1)
    # plot( out$vgm$distance, out$vgm$sv,  xlim=xlim, ylim=ylim, type="n" )
    points( out$vgm$distance, out$vgm$sv, col="red", pch=19  )
    
    if (0) {
      plotdata=TRUE
      autocorrelation_localrange=0.1 
      weight_by_inverse_distance=TRUE
      internal_scale=NULL
      fit_model=TRUE
      nx = NULL
    }

  }

  zvar = var( z, na.rm=TRUE )
  if (!is.finite(zvar)) zvar = 0 
  if (zvar == 0) return(NULL)

  zmean = mean( z, na.rm=TRUE)
  zsd = sqrt( zvar )


  out = list(
    Ndata = length(z),
    zmean = zmean,
    zsd = zsd,
    zvar = zvar,  # this is the scaling factor for semivariance .. diving by sd, below reduces numerical floating point issues
    range_crude = sqrt( diff(range(xy[,1]))^2 + diff(range(xy[,2]))^2) / 4  # initial scaling distance beyond which domain edge effects become important
  )
 
 


 # A crude GUESS AT PHI:
  out$internal_scale  = ifelse( is.null(internal_scale), matern_distance2phi( out$range_crude, nu=0.5 ), 1)  # the presumed scaling distance to make calcs use smaller numbers

  out$vgm_dist_max = out$internal_scale
  out$vgm_var_max = zvar
  out$summary = list()
  
  xy = xy /  out$vgm_dist_max
  z = (z - out$zmean) / out$zsd
  vxs = as.matrix( dist( xy, diag=TRUE, upper=TRUE))  # distance matrix between knots
 

 
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
      using SpecialFunctions, Distributions, Statistics, Optim, PDMats, LinearAlgebra
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

      @model function matern_fit(vxs, z)
        n = size(vxs, 1)
        eps = 1.0e-12 # to avoid zeros
        nu  ~ truncated( Cauchy( 1.0, 0.5), 0.01, 10.0 )
        phi ~ truncated( Cauchy( 0.0, 0.1), eps, 1.0 )
        vtot ~ truncated( Cauchy( 0.0, 0.5), eps, 10.0 )
        vfsp ~ Beta( 2, 2 )
        yvar ~ truncated( Cauchy( 0.0, 0.5), eps, 10.0 )
        mu ~ filldist( Cauchy(0.0, 0.5), n) 
        # Sigma = spatial covariance matrix
        Sigma = tzeros(n,n) 
        @.  Sigma =  1.0 - vfsp * matern( nu, sqrt(2.0*nu) * (vxs / phi ) )   
        diag(Sigma) = (1.0 - vfsp) + eps
        @. Sigma =  vtot * Sigma 
        Sigma = Matrix( Hermitian( Sigma ) ) 
        if (!isposdef(Sigma)) 
          Turing.@addlogprob! -Inf  
          return
        end
        # z ~ MvNormal( mu, Sigma )
     
      end
    ')

    test = matern_fit(vxs, z)

    #  Run sampler, collect results.


      # mle_estimate = optimize(matern_fit(vxs, z), MLE())
      # mle_estimate = optimize(matern_fit(vxs, z), MLE(), NelderMead())
      # mle_estimate = optimize(matern_fit(vxs, z), MLE(), SimulatedAnnealing())
      # mle_estimate = optimize(matern_fit(vxs, z), MLE(), ParticleSwarm())
      # mle_estimate = optimize(matern_fit(vxs, z), MLE(), Newton())
      # mle_estimate = optimize(matern_fit(vxs, z), MLE(), AcceleratedGradientDescent())
      # mle_estimate = optimize(matern_fit(vxs, z), MLE(), Newton(), Optim.Options(iterations=10_000, allow_f_increases=true))

      # using StatsBase
      # coeftable(mle_estimate)

      # map_estimate = optimize(matern_fit(vxs, z), MAP())

      # c1 = sample(matern_fit(vxs, z), SMC(), 1000)
      # c2 = sample(matern_fit(vxs, z), PG(10), 1000)
      # c3 = sample(matern_fit(vxs, z), HMC(0.1, 5), 1000)
      # c4 = sample(matern_fit(vxs, z), Gibbs(PG(10, :m), HMC(0.1, 5, :s²)), 1000)
      # c5 = sample(matern_fit(vxs, z), HMCDA(0.15, 0.65), 1000)
      # c6 = sample(matern_fit(vxs, z), NUTS(0.65), 1000)

      # c0 = sample(matern_fit(vxs, z), NUTS(0.65), 1000, init_params = map_estimate.values.array)
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

    # chain = Turing$sample( matern_fit(vxs, z), Turing$PG(10L), niters, n_adapts=warmup, progress=FALSE, verbose=FALSE, drop_warmup=TRUE )
    
    
    chain = Turing$sample( matern_fit(vxs, z), Turing$SMC(), niters, n_adapts=warmup, progress=FALSE, verbose=FALSE, drop_warmup=TRUE )
    
    # chain = Turing$sample( matern_fit(vxs, z), Turing$NUTS(0.95), niters, n_adapts=warmup, progress=FALSE, verbose=FALSE, drop_warmup=TRUE )
    
    Turing$summarize( chain)

    # graphical backed (GR) not working with R 
    # JuliaDB = juliaImport("c")
    # JuliaDB$plot( sss )
    
    res = JuliaConnectoR::juliaGet(chain)
    o = as.data.frame(res$value$data[(warmup+1):niters, ,] )
    names( o ) = as.character(unlist(res$value$axes[[2]][["val"]]))
 
    out$samples = o
    i = grep("mu\\[", names(o))
    mu = o[,i]
    out$mu = apply(mu, 2, mean )
    out$mu_sd = apply(mu, 2, sd )
    
    # plot(z ~ out$mu)
    
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


    if( plotdata ) {
      dev.new()
      # xlim= c(0, max(out$vgm$distance) *1.1 )
      #ylim= c(0, vgm_var_max*1.1)
      # plot( out$vgm$distance, out$vgm$sv,  xlim=xlim, ylim=ylim, type="n" )
      # plot( out$vgm$distance, out$vgm$sv, col="darkslategray"  )
      ds = seq( 0, out$range_crude  , length.out=100 )
      ac = out$summary$varObs + out$summary$varSpatial*(1 - stmv_matern( ds, out$summary$phi, out$summary$nu ) )
      lines( ds, ac, type="l", col="darkgreen", lwd=3 )
      abline( h=0, lwd=1, col="lightgrey" )
      abline( v=0 ,lwd=1, col="lightgrey" )
      abline( h=out$summary$varObs, lty="dashed", col="grey" )
      abline( h=out$summary$varObs + out$summary$varSpatial, lty="dashed", col="grey" )
      abline( v=out$summary$localrange, lty="dashed", col="grey")
    }

}


