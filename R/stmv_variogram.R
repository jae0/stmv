
stmv_variogram = function( XYZ=NULL, xy=NULL, z=NULL, ti=NULL,
  plotdata=FALSE, methods=c("geoR"), discretized_n=64, nbreaks = c( 16, 32, 64, 13, 17, 21, 33, 49 ),
  distance_cutoff=NA, family=gaussian(link="identity"),
  range_correlation=0.1, stmv_autocorrelation_fft_taper=0,
  stanmodel=NULL, modus_operandi="easygoing", discretize_data=FALSE ) {

  #\\ estimate empirical variograms (actually correlation functions)
  #\\ and then model them using a number of different approaches .. using a Matern basis
  #\\ returns empirical variogram and parameter estimates, and the models themselves
  #\\ expect xy = c(plon, plat), z=variable

  #\\ NOTE:: the default parameterization is Wikipedia's paramterization:
  #\\ == Rasmussen, Carl Edward (2006) Gaussian Processes for Machine Learning
  #\\ == RandomFields and geostatsp parameterization ss
  #\\  sigma^2 * (2^{1-nu} / Gamma(nu) ) * (sqrt(2*nu) * ||x|| / phi)^{nu} * K_{nu}( sqrt(2*nu) * ||x|| / phi )
  #\\   where K_{nu} is the Bessel function with smooth nu and phi is known as the range parameter
  # -------------------------


      if ( 0 ) {
        # debugging / comparison of results

        # project.library("aegis", "aegis", "stmv"  )
        loadfunctions( c("aegis", "aegis", "stmv" ))
        XYZ = NULL
        xyz = stmv_test_data( datasource="swiss" )
        xy = xyz[, c("x", "y")]
        mz = log( xyz$rain )
        mm = lm( mz ~ 1 )
        z = residuals( mm)

        plotdata=TRUE
        nbreaks = 15
        distance_cutoff=NA
        family=gaussian(link="identity")
        range_correlation=0.1

        discretize_data=FALSE
        discretized_n=64
        nbreaks = 32
        range_correlation=0.1
        stmv_autocorrelation_fft_taper=0

        microbenchmark::microbenchmark( {gr = stmv_variogram( xy, z, methods="fft", plotdata=FALSE )}, times= 10 )  # 63 milli sec
        gr = stmv_variogram( xy, z, methods="fft", plotdata=TRUE ) # nls


        microbenchmark::microbenchmark( {gr = stmv_variogram( xy, z, methods="optim", plotdata=FALSE )}, times= 10 )  # 63 milli sec
        gr = stmv_variogram( xy, z, methods="optim", plotdata=TRUE ) # nls

        # $optim$vgm_var_max
        # [1] 0.7927

        # $optim$vgm_dist_max
        # [1] 171907

        # $optim$autocorrelation_function
        # [1] "matern"

        # $optim$nu
        # [1] 3

        # $optim$phi
        # [1] 26532

        # $optim$varSpatial
        # [1] 0.4464

        # $optim$varObs
        # [1] 0.1385

        # $optim$range
        # [1] 69508

        # $optim$phi_ok
        # [1] TRUE

        # $optim$objfn
        # [1] 0.1079



        microbenchmark::microbenchmark( {gr = stmv_variogram( xy, z, methods="fast.recursive", plotdata=FALSE )}, times= 10 )  # 138.3 milli sec  .. unstable
        # tests
        gr = stmv_variogram( xy, z, methods="fast.recursive", plotdata=TRUE ) # nls
# $ Ndata              : int 100
# $ varZ               : num 0.578
# $ range_crude        : num 88029
# $ stmv_internal_scale: num 29385
# $ distance_cutoff    : num 132043
# ..$ range     : num 67830
# ..$ nu        : num 1
# ..$ phi       : num 23990
# ..$ varSpatial: num 0.588
# ..$ varObs    : num 0
# ..$ phi_ok  : logi TRUE

      microbenchmark::microbenchmark( {gr = stmv_variogram( xy, z, methods="multipass", plotdata=FALSE )}, times= 10 )  # 31 milli sec
            gr = stmv_variogram( xy, z, methods="multipass", plotdata=TRUE ) # nls
        # $fast$vgm_var_max: 0.7987
        # $fast$vgm_dist_max: 171425
        # $fast$autocorrelation_function: "matern"
        # $fast$nu: 0.5
        # $fast$phi: 29060
        # $fast$varSpatial: 0.5956
        # $fast$varObs: 0
        # $fast$range: 87056
        # $fast$phi_ok: TRUE
        # $fast$objfn: 0.1278

        gr = stmv_variogram( xy, z, methods="fields", plotdata=TRUE ) # ml via profile likelihood
        microbenchmark::microbenchmark( {gr = stmv_variogram( xy, z, methods="fields", plotdata=FALSE )}, times= 10 )  # 1.38 sec
        # fields$range [1] 80510
        # $fields$nu [1] 0.5
        # $fields$phi 26875
        # $fields$varSpatial [1] 0.294
        # $fields$varObs [1] 1.333e-05

        gr = stmv_variogram( xy, z, methods="geoR", plotdata=TRUE ) # ml
        microbenchmark::microbenchmark( {gr = stmv_variogram( xy, z, methods="geoR", plotdata=FALSE )}, times= 10 )  # 45 MILLI Sec   # geoR seems to crash a node ..
       # $geoR$range
        # [1] 60278
        # $geoR$varSpatial
        # [1] 0.572
        # $geoR$varObs
        # [1] 0
        # $geoR$nu
        # [1] 1.127
        # $geoR$phi
        # [1] 21528
        # $geoR$phi_ok
        # [1] TRUE

        # plot( gr$geoR$vgm )
        # lines( gr$geoR$fit, lwd=2, col="slateblue" )
        # xRange = c( 0, max(gr$geoR$range*2.1 ) )
        # yRange = c( 0, max(gr$geoR$vgm$v )*1.05 )
        # plot ( gr$geoR$vgm$v ~ gr$geoR$vgm$u, pch=20, xlim=xRange, ylim=yRange, ylab="Semivariance", xlab="Distance" )
        #   abline( h=0,  col="gray", lwd=2 )
        #   abline( h= (gr$geoR$varSpatial + gr$geoR$varObs), lty="dashed", col="slategray"  )
        #   abline( h=  gr$geoR$varObs , lty="dashed", col="slategray")
        #   abline( v=gr$geoR$range, lty="dotted", col="slateblue" )
        #   abline( v=0,  col="gray", lwd=2 )
        #   x = seq( 0, 2*gr$geoR$range, length.out=100 )
        #   acor = geoR::matern( x, phi=gr$geoR$phi, kappa=gr$geoR$kappa  )
        #   acov = gr$geoR$varObs +  gr$geoR$varSpatial * (1- acor)
        #   lines( acov ~ x , col="blue", lwd=2 )

        microbenchmark::microbenchmark( {gr = stmv_variogram( xy, z, methods="geoR.ML", plotdata=FALSE )}, times= 10 )  # 119 mili sec
        gr = stmv_variogram( xy, z, methods="geoR.ML", plotdata=TRUE ) # ml
        #  unstable

        microbenchmark::microbenchmark( {gr = stmv_variogram( xy, z, methods="gstat", plotdata=FALSE )}, times= 10 )  # 86 milli sec
        gr = stmv_variogram( xy, z, methods="gstat", plotdata=TRUE ) # ml
        # $gstat$range
        # [1] 65461
        # $gstat$nu
        # [1] 1.1
        # $gstat$phi
        # [1] 23333
        # $gstat$varSpatial
        # [1] 0.5776
        # $gstat$varObs
        # [1] 0
        # $gstat$phi_ok
        # [1] TRUE

        microbenchmark::microbenchmark( {gr = stmv_variogram( xy, z, methods="RandomFields", plotdata=FALSE )}, times= 10 ) #  1.349 sec
        gr = stmv_variogram( xy, z, methods="RandomFields", plotdata=TRUE ) # ml via rf
        # $RandomFields$range
        # [1] 65343

        # $RandomFields$varSpatial
        # [1] 0.6561

        # $RandomFields$varObs
        # [1] 2.862e-05

        # $RandomFields$phi
        # [1] 22851

        # $RandomFields$nu
        # [1] 0.8727

        # $RandomFields$error
        # [1] NA

        # $RandomFields$phi_ok
        # [1] TRUE

        gr = stmv_variogram( xy, z, methods="CompRandFld", plotdata=TRUE ) # unstable results
        microbenchmark::microbenchmark( {gr = stmv_variogram( xy, z, methods="CompRandFld", plotdata=FALSE )}, times= 10 )  # 2.34 sec! .. slow!
        # Ndata
        # [1] 100

        # $varZ
        # [1] 0.5782

        # $range_crude
        # [1] 88029

        # $stmv_internal_scale
        # [1] 29385

        # $distance_cutoff
        # [1] 132043

        # $CompRandFld
        # $CompRandFld$varObs
        # [1] 9.182e-12

        # $CompRandFld$varSpatial
        # [1] 0.5211

        # $CompRandFld$nu
        # [1] 0.4918

        # $CompRandFld$phi
        # [1] 44705

        # $CompRandFld$range
        # [1] 134102

        # $CompRandFld$phi_ok
        # [1] FALSE


        gr = stmv_variogram( xy, z, methods="spBayes", plotdata=TRUE ) # mcmc
        microbenchmark::microbenchmark( {gr = stmv_variogram( xy, z, methods="spBayes", plotdata=FALSE )}, times= 10 )  # v. slow! ... 29.3 sec
        # $spBayes$range: 66377
        # $spBayes$varSpatial: 0.4683
        # $spBayes$varObs: 0.1073
        # $spBayes$phi: 24778
        # $spBayes$nu: 2.062

        # out = gsp
        # nd = nrow(out$spBayes$recover$p.theta.samples)
        # rr = rep(NA, nd )
        # for (i in 1:nd) rr[i] = geoR::practicalRange("matern", phi=1/out$spBayes$recover$p.theta.samples[i,3], kappa=out$spBayes$recover$p.theta.samples[i,4] )
        # #  range = matern_phi2distance(phi=phi, nu=nu, cor=range_correlation)

        # hist(rr)  # range estimate

        # hist( out$spBayes$recover$p.theta.samples[,1] ) #"sigma.sq"
        # hist( out$spBayes$recover$p.theta.samples[,2] ) # "tau.sq"
        # hist( out$spBayes$recover$p.theta.samples[,3] ) # 1/phi
        # hist( out$spBayes$recover$p.theta.samples[,4] ) # nu


       gr = stmv_variogram( xy, z, methods="bayesx", plotdata=TRUE )
        microbenchmark::microbenchmark( {gr = stmv_variogram( xy, z, methods="bayesx", plotdata=FALSE )}, times= 10 )  # 90 milli sec
        # $bayesx$range
        # [1] 103783

        # $bayesx$model
        # Call:
        # bayesx(formula = z ~ sx(plon, plat, bs = "kr"), data = xys, family = family$family,
        #     method = "REML")
        # Summary:
        # N = 100  df = 9.725  AIC = 0.547  BIC = 25.882
        # logLik = 9.451  method = REML  family = gaussian

        # $bayesx$varSpatial
        # [1] 1.166

        # $bayesx$varObs
        # [1] 0.3066

        # $bayesx$nu
        # [1] 0.5

        # $bayesx$phi
        # [1] 34644

        # $bayesx$phi_ok
        # [1] TRUE



        gr = stmv_variogram( xy, z, methods="inla", plotdata=TRUE )
        microbenchmark::microbenchmark( {gr = stmv_variogram( xy, z, methods="inla", plotdata=FALSE )}, times= 10 )  # 4.4 sec
        # $inla$range.inla.practical
        # [1] 52481
        # $inla$varSpatial
        # [1] 0.6499
        # $inla$varObs
        # [1] 0.0001482
        # $inla$phi
        # [1] 25407
        # $inla$nu
        # [1] 1
        # $inla$error
        # [1] NA
        # $inla$range
        # [1] 71835
        # $inla$phi_ok
        # [1] TRUE


        # Laplace Approximation:  too slow to use

        #              Mean      SD     MCSE  ESS         LB     Median         UB
        # tausq      4.585e+00 0.00000 0.000000 1000  4.585e+00  4.585e+00  4.585e+00
        # sigmasq   -8.461e-02 0.00000 0.000000 1000 -8.461e-02 -8.461e-02 -8.461e-02
        # phi       -4.978e-01 0.00000 0.000000 1000 -4.978e-01 -4.978e-01 -4.978e-01
        # nu         6.286e-03 0.00000 0.000000 1000  6.286e-03  6.286e-03  6.286e-03
        # Deviance   3.743e+02 0.03871 0.001224 1000  3.743e+02  3.743e+02  3.744e+02
        # LP        -1.954e+02 0.01973 0.000624 1000 -1.954e+02 -1.954e+02 -1.953e+02

        # Laplaces Demon:
        #             Mean        SD      MCSE       ESS         LB     Median         UB
        # tausq      1.187e+00  0.5468  0.16983   8.251    0.35260    1.13508    2.2354
        # sigmasq    8.280e-01  0.5387  0.14286  21.539    0.09564    0.72554    2.4146
        # phi        1.761e+00  0.7113  0.20396  14.262    0.46596    1.77338    3.2806
        # nu         1.415e+00  0.8262  0.24249  10.989    0.21269    1.35122    3.0988
        # Deviance   3.204e+02 36.0896 12.04187   2.287  251.50445  335.42134  366.2237
        # LP        -1.667e+02 19.0692  6.36485   2.197 -192.71295 -174.44162 -130.6786

        # stan:
        #                      mean se_mean    sd    2.5%     25%     50%    75%  97.5% n_eff Rhat
        # sigma_sq             5.73    1.26 13.47    0.08    1.03    2.44   5.65  32.17   115 1.02
        # tau_sq               4.47    1.04  7.66    0.14    0.95    2.08   4.87  23.02    54 1.05
        # phi                  1.32    0.05  0.79    0.09    0.71    1.25   1.80   3.04   275 1.01
        # out$stmv_internal_scale = 30485 ;
        # range = 40240
    } #end if

  # ----- start -----
  

  if (is.null(XYZ) ) XYZ = cbind( xy, z) 

  XYZ = as.data.frame( XYZ )
  names(XYZ) =  c("plon", "plat", "z" ) # arbitrary
  rownames( XYZ) = 1:nrow(XYZ)  # RF seems to require rownames ...
  
  if (discretize_data) XYZ = stmv_discretize_coordinates(coo=XYZ[,c(1,2)], z=XYZ[,3], discretized_n=discretized_n, method="aggregate", FUNC=mean, na.rm=TRUE)

  
  out = list(
    Ndata = nrow(XYZ), 
    stmv_internal_mean = mean( as.vector(XYZ[,3]), na.rm=TRUE ), 
    varZ = var( as.vector(XYZ[,3]), na.rm=TRUE ),  # this is the scaling factor for semivariance .. diving by sd, below reduces numerical floating point issues
    range_crude = sqrt( diff(range(XYZ[,1]))^2 + diff(range(XYZ[,2]))^2) / 4  #initial scaling distance
  )
  
  out$stmv_internal_sd = sqrt( out$varZ )
  out$stmv_internal_scale = matern_distance2phi( out$range_crude, nu=0.5 )  # the presumed scaling distance to make calcs use smaller numbers
  # if max dist not given, make a sensible choice using exponential variogram as a first estimate
  out$distance_cutoff = ifelse( is.na(distance_cutoff), out$range_crude * 1.5, distance_cutoff)

  XYZ[,1] = XYZ[,1] / out$stmv_internal_scale
  XYZ[,2] = XYZ[,2] / out$stmv_internal_scale
  XYZ[,3] = (XYZ[,3] - out$stmv_internal_mean) / out$stmv_internal_sd


  # ------------------------


  if ("fft" %in% methods) {
    
    # empirical variogram by fftw2d
    finished = FALSE
    for ( i in 1:length(nbreaks)) {
      VGM = NULL
      VGM = stmv_variogram_fft( xyz=XYZ, nx=discretized_n, ny=discretized_n, nbreaks=nbreaks[i]  )  
      if (is.null(VGM)) next()
      uu = which( (VGM$vgm$distances < 0.9*max(VGM$vgm$distances) ) & is.finite(VGM$vgm$sv) )
      fit = try( stmv_variogram_optimization( vx=VGM$vgm$distances[uu], vg=VGM$vgm$sv[uu], plotvgm=plotdata,
        stmv_internal_scale=out$stmv_internal_scale, cor=range_correlation  ))
      if ( !inherits(fit, "try-error") ) {
        out$fft = fit$summary
        if (exists("phi", out$fft)) {
          if (is.finite(out$fft$phi)) {
             out$fft$phi  = out$fft$phi  * out$stmv_internal_scale
             out$fft$localrange = out$fft$localrange * out$stmv_internal_scale
             if ( out$fft$phi  < out$distance_cutoff*0.99 ) {
              finished = TRUE
             } 
          }
        }
      }
      if (finished) break()
    }

    if (exists("VGM")) out$vgm = VGM

    if (!finished) return(out)

    if (exists("summary", fit)) {
      if ( !inherits(fit, "try-error") ) {
        
        if (exists("phi", out$fft)) {
          if (is.finite(out$fft$phi)) {
            out$fft$phi_ok = ifelse( out$fft$phi < out$distance_cutoff*0.99, TRUE, FALSE )  # using distance as an upper limit
          }
        }
      }
    }

    if (plotdata) {
      xlim= c(0, fit$summary$vgm_dist_max*1.1)
      ylim= c(0, fit$summary$vgm_var_max*1.1)
      plot( fit$summary$vx, fit$summary$vg, col="green", xlim=xlim, ylim=ylim )
      ds = seq( 0, fit$summary$vgm_dist_max, length.out=100 )
      ac = fit$summary$varObs + fit$summary$varSpatial*(1 - stmv_matern( ds, fit$summary$phi, fit$summary$nu ) )
      lines( ds, ac, col="orange" )
      abline( h=0, lwd=1, col="lightgrey" )
      abline( v=0 ,lwd=1, col="lightgrey" )
      abline( h=fit$summary$varObs, lty="dashed", col="grey" )
      abline( h=fit$summary$varObs + fit$summary$varSpatial, lty="dashed", col="grey" )
      localrange = matern_phi2distance( phi=fit$summary$phi, nu=fit$summary$nu, cor=range_correlation )
      abline( v=localrange, lty="dashed", col="grey")
    }

    return(out)
  }




  # ------------------------


  if ("optim" %in% methods) {
    # spatial discretization

    require( RandomFields )
    VGM = RFvariogram( data=RFspatialPointsDataFrame( coords=XYZ[,c(1,2)], data=XYZ[,3], RFparams=list(vdim=1, n=1) ) )
    # remove the (0,0) point -- force intercept

    todrop = which( !is.finite(VGM@empirical )) # occasionally NaN's are created!
    todrop = unique( c(1, todrop) )
    vg = VGM@empirical[-todrop]
    vx = VGM@centers[-todrop]
    out$vgm = VGM
    fit = try( stmv_variogram_optimization( vx=vx, vg=vg, nu=0.5, plotvgm=plotdata, stmv_internal_scale=out$stmv_internal_scale, cor=range_correlation ))
    if ( !inherits(fit, "try-error") ) {
      out$optim = fit$summary
      if (exists("phi", out$optim)) {
        if (is.finite(out$optim$phi)) {
          out$optim$phi_ok = ifelse( out$optim$phi < out$distance_cutoff*0.99, TRUE, FALSE )
        }
      }
    }

      if (plotdata) {
        xlim= c(0, fit$summary$vgm_dist_max*1.1)
        ylim= c(0, fit$summary$vgm_var_max*1.1)
        plot( fit$summary$vx, fit$summary$vg, col="green", xlim=xlim, ylim=ylim )
        ds = seq( 0, fit$summary$vgm_dist_max, length.out=100 )
        ac = fit$summary$varObs + fit$summary$varSpatial*(1 - stmv_matern( ds, fit$summary$phi, fit$summary$nu ) )
        lines( ds, ac, col="orange" )
        abline( h=0, lwd=1, col="lightgrey" )
        abline( v=0 ,lwd=1, col="lightgrey" )
        abline( h=fit$summary$varObs, lty="dashed", col="grey" )
        abline( h=fit$summary$varObs + fit$summary$varSpatial, lty="dashed", col="grey" )
        localrange = matern_phi2distance( phi=fit$summary$phi, nu=fit$summary$nu, cor=range_correlation )
        abline( v=localrange, lty="dashed", col="grey")
      }

    return(out)
  }



  # -------------------------
  # ------------------------


  if ( "fast.recursive" %in% methods)  {
    # gives a fast stable empirical variogram using nl least squares and a coarse-grained discretization
    # essentially the same as gstat method but coarse grained discretization

    # nu=0.5, discretized_n=100, nbreaks = 13

      # nu=0.5 defaults to exponential and uses gstat
      #\\ NOTE:: the default parameterization is Wikipedia's paramterization:
      #\\ == Rasmussen, Carl Edward (2006) Gaussian Processes for Machine Learning
      #\\  sigma^2 * (2^{1-nu} / Gamma(nu) ) * (sqrt(2*nu) * ||x|| / phi)^{nu} * K_{nu}( sqrt(2*nu) * ||x|| / phi)
      #\\   where K_{nu} is the Bessel function with smooth nu and phi is known as the range parameter
      #\\ As usage sometimes is for high density data, aggregation to a coarse resolution of 'discretized_n' units along
      #\\ the smaller dimension  before computation.
      # -------------------------

      # spatial discretization only
      XYZ = stmv_discretize_coordinates(coo=xy, z=z, discretized_n=discretized_n, method="aggregate", FUNC=mean, na.rm=TRUE)
      names(XYZ) =  c("plon", "plat", "z" ) # arbitrary

      maxdist = out$range_crude   # begin with this (diagonal)

      # gives a fast stable empirical variogram using nl least squares

      XYZ$plon = XYZ$plon / out$stmv_internal_scale  # keeps things smaller in value to avoid floating point issues
      XYZ$plat = XYZ$plat / out$stmv_internal_scale  # keeps things smaller in value to avoid floating point issues

      # empirical variogram
      vEm = gstat::variogram( z ~ 1, locations=~plon+plat, data=XYZ, cutoff=maxdist/out$stmv_internal_scale, width=maxdist/out$stmv_internal_scale/nbreaks, cressie=FALSE )

      vEm$dist0 = vEm$dist * out$stmv_internal_scale

      fit = stmv_variogram_optimization( vx=vEm$dist0, vg=vEm$gamma, plotvgm=FALSE, stmv_internal_scale=out$stmv_internal_scale, cor=range_correlation ) # nu=0.5 == exponential variogram

      if ( fit$summary$phi_ok ) {

        out$fast.recursive = fit$summary
        out$fast.recursive$fit=fit
        out$fast.recursive$vgm=vEm

      } else {
        cnt = 0
        localrange = maxdist

        while ( cnt < 5  ) {
          maxdist = maxdist * 1.25
          if ( localrange > maxdist ) {
            # message ( "Autocorrelation range greater than data range .. retrying a last time at max dist with more data")
            return(NULL)
          }
          # message( "Range longer than distance cutoff ... retrying with a larger distance cutoff")
          cnt = cnt + 1

          vEm = gstat::variogram( z~1, locations=~plon+plat, data=XYZ, cutoff=maxdist/out$stmv_internal_scale, width=maxdist/(out$stmv_internal_scale*nbreaks), cressie=FALSE )  # empirical variogram
          vEm$dist0 = vEm$dist * out$stmv_internal_scale
          vMod0 = gstat::vgm(psill=2*0.75*runif(1), model="Mat", range=2*1*runif(1), nugget=2*0.25*runif(1), kappa=0.5 ) # starting model parameters, 2*X as the expected value of runif(1) = 0.5
          vFitgs =  try( gstat::fit.variogram( vEm, vMod0, fit.kappa =TRUE, fit.sills=TRUE, fit.ranges=TRUE ) )
            # gstat's kappa is the Bessel function's "nu" smoothness parameter
            # gstat::"range" == range parameter == phi
          if (inherits(vFitgs, "try-error") )  return(NULL)

          phi = matern_phi2phi( mRange=vFitgs$range[2], mSmooth=vFitgs$kappa[2], parameterization_input="gstat", parameterization_output="stmv" ) * out$stmv_internal_scale
          nu=vFitgs$kappa[2]

          out$fast.recursive = list( fit=vFitgs, vgm=vEm, nu=vFitgs$kappa[2], phi=phi,
            varSpatial=vFitgs$psill[2], varObs=vFitgs$psill[1]  )
          out$fast.recursive$phi_ok = FALSE
          if ( out$fast.recursive$phi < max(out$fast.recursive$vgm$dist0)*0.99) out$fast.recursive$phi_ok = TRUE
          localrange = matern_phi2distance( phi=phi, nu=nu, cor=range_correlation )
          if ( localrange < maxdist ) break()
        }
      }


      if (plotdata) {
        xub = max(out$distance_cutoff) *1.25
        plot.new()
        plot(vEm, model=vFitgs, add=T)
        plot.new()
        plot( gamma ~ dist0, data=out$fast.recursive$vgm, xlim=c(0,xub),
             ylim=c(0,max(out$fast.recursive$vgm$gamma)*1.1), col="blue", pch=20 )
        abline( h=out$fast.recursive$varSpatial + out$fast.recursive$varObs )
        abline( h=out$fast.recursive$varObs )
        abline( v=localrange )
        abline (v=0)
        x = seq( 0, xub, length.out=100 )
        acor = stmv_matern( x, mRange=out$fast.recursive$phi, mSmooth=out$fast.recursive$nu  )
        acov = out$fast.recursive$varObs + out$fast.recursive$varSpatial * (1- acor)
        lines( acov~x , col="red" )

        if (0) {
          # looks at the predictions
          gs <- gstat(id = "z", formula = z~1, locations=~plon+plat, data=xy, maxdist=distance_cutoff, nmin=10, force=TRUE, model=vFitgs )
          # variogram of residuals
          data(meuse.grid)
          meuse.grid$plon = meuse.grid$x
          meuse.grid$plat = meuse.grid$y

          preds <- predict(gs, newdata=meuse.grid )
          spplot(preds)

        }
      }

      return(out)

      if( 0) {
        xlim= c(0, fit$summary$vgm_dist_max*1.1)
        ylim= c(0, fit$summary$vgm_var_max*1.1)
        plot( fit$summary$vx, fit$summary$vg, col="green", xlim=xlim, ylim=ylim )
        ds = seq( 0, fit$summary$vgm_dist_max, length.out=100 )
        ac = fit$summary$varObs + fit$summary$varSpatial*(1 - stmv_matern( ds, fit$summary$phi, fit$summary$nu ) )
        lines( ds, ac, col="orange" )
        abline( h=0, lwd=1, col="lightgrey" )
        abline( v=0 ,lwd=1, col="lightgrey" )
        abline( h=fit$summary$varObs, lty="dashed", col="grey" )
        abline( h=fit$summary$varObs + fit$summary$varSpatial, lty="dashed", col="grey" )
        abline( v=localrange, lty="dashed", col="grey")
      }
  }


  # ------------------------
  # ------------------------

  if ("CompRandFld" %in% methods) {
    require( CompRandFld )
    stop( "CompRandFld seems unstable .. force stop")
    fit = FitComposite( data=z, coordx=as.matrix(xy/out$stmv_internal_scale), corrmodel="matern",
                       maxdist=out$distance_cutoff/out$stmv_internal_scale, optimizer="BFGS" )
    out$CompRandFld = list(
      varObs=fit$param[["nugget"]],
      varSpatial = fit$param[["sill"]],
      nu = fit$param[["smooth"]]
    )
    out$CompRandFld$phi =  matern_phi2phi( mRange=fit$param[["scale"]], mSmooth=out$CompRandFld$nu, parameterization_input="CompRandFld", parameterization_output="stmv" ) * out$stmv_internal_scale
    out$CompRandFld$phi_ok = ifelse( out$CompRandFld$phi < out$distance_cutoff*0.99, TRUE, FALSE )
    localrange = matern_phi2distance(phi=out$CompRandFld$phi, nu=out$CompRandFld$nu, cor=range_correlation)
    if( plotdata) {
      VGM = EVariogram(data=z, coordx=as.matrix(xy/out$stmv_internal_scale),
                        maxdist=out$distance_cutoff/out$stmv_internal_scale, numbins=nbreaks)
      vg = VGM$variograms
      vx = VGM$centers
      xlim= c(0, max(VGM$centers)*1.1) * out$stmv_internal_scale
      ylim= c(0, max(VGM$variograms)*1.1)
      plot( VGM$centers* out$stmv_internal_scale, VGM$variograms, col="green", xlim=xlim, ylim=ylim )
      ds = seq( 0, max(VGM$centers), length.out=100 ) * out$stmv_internal_scale
      ac = out$CompRandFld$varObs + out$CompRandFld$varSpatial*(1 - stmv_matern( ds, out$CompRandFld$phi, out$CompRandFld$nu ) )
      lines( ds, ac, col="orange" )
      abline( h=0, lwd=1, col="lightgrey" )
      abline( v=0 ,lwd=1, col="lightgrey" )
      abline( h=out$CompRandFld$varObs, lty="dashed", col="grey" )
      abline( h=out$CompRandFld$varObs + out$CompRandFld$varSpatial, lty="dashed", col="grey" )
      abline( v=localrange, lty="dashed", col="grey")
    }
    return(out)
  }


  # ------------------------
  # ------------------------


  if ("fields" %in% methods) {

    require(fields)

    XYZ = stmv_discretize_coordinates(coo=xy, z=z, discretized_n=discretized_n, method="aggregate", FUNC=mean, na.rm=TRUE)
    names(XYZ) =  c("plon", "plat", "z" ) # arbitrary
    xy = XYZ[,c("plon", "plat")]
    z = XYZ$z
    XYZ = NULL
    nu = 0.5 # 0.5 == exponential .. ie. fixed
    res =NULL
    ## NOTE: process variance (rho); range (theta); nugget (sigma**2)
    ## lambda= sigma**2/ rho and rho. Thinking about h as the spatial signal and e as the noise lambda can be interpreted
    ##  as the noise to signal variance ratio in this spatial context
    # MLESpatialProcess is a ML method for Gaussian spatial process
    fsp = MLESpatialProcess(xy/out$stmv_internal_scale, z, cov.function = "stationary.cov",
      # abstol=1e-3,
      # lambda.start=0.5, theta.start=1.0, theta.range=c(0.2, 4.0), gridN=25,
      # optim.args=list(method = "BFGS", control = list(fnscale=-1, parscale=c(0.5, 0.5), ndeps=c(0.05,0.05))),
      cov.args = list(Covariance = "Matern", smoothness = nu)
    )
    if( is.finite(sum(fsp$summary))) res = fsp$summary
    if (is.null(res)) return(NULL)

    # warning ("vgram is really slow ...")

    vg = vgram( xy/out$stmv_internal_scale, z, N=nbreaks, dmax=out$distance_cutoff/out$stmv_internal_scale )
    vgm = Matern( d=vg$centers, range=res[["theta"]], smoothness=nu )
    cvg = data.frame( cbind( x=vg$centers*out$stmv_internal_scale, cvgm= (res[["sigmaMLE"]]^2 + res[["rhoMLE"]] * (1-vgm)) ))
    out$fields = list( fit=fsp, vgm=cvg, nu=nu, phi=NA,
      varSpatial=res[["rhoMLE"]], varObs=res[["sigmaMLE"]]^2  )  # fields::"range" == range parameter == phi

    out$fields$phi = matern_phi2phi( mRange=res[["theta"]], mSmooth=out$fields$nu,
      parameterization_input="fields", parameterization_output="stmv") * out$stmv_internal_scale
    localrange = matern_phi2distance(phi=out$fields$phi, nu=out$fields$nu, cor=range_correlation)
    out$fields$phi_ok = ifelse( out$fields$phi < out$distance_cutoff*0.99, TRUE, FALSE )

    if( plotdata ){
#      xub = max(out$distance_cutoff, max(out$geoR$vgm$u, out$distance_cutoff)) *1.25
      plot.new()

      vEm = try( variog( coords=xy/out$stmv_internal_scale, data=z, uvec=nbreaks, max.dist=out$distance_cutoff/out$stmv_internal_scale ) )
      vEm$u0 = vEm$u * out$stmv_internal_scale
      plot( vEm$v ~ vEm$u0, pch=20 ,
            xlim=c(0, max( c(cvg$x, vEm$u0) )), ylim=c(0, max( c(out$fields$varSpatial + out$fields$varObs, out$varZ, cvg$cvgm, vEm$v) ) ) )

      points( cvg, type="b", ylim=range( c(0, cvg$cvgm) ) )
      abline( h=out$fields$varSpatial + out$fields$varObs)
      abline( h=out$fields$varObs )
      abline( v=localrange )


    #   plot.new()
    #   lambda.MLE<- out$fields$varObs / out$fields$varSpatial  # ratio nugget / sill variance
    #   fsp2<- Krig( xy, z, Covariance="Matern", theta=fsp$pars["theta"], smoothness=nu, lambda= lambda.MLE)
    #   surface(fsp2)

    #   plot.new()
    #   fsp.p<- predictSurface(fsp2, lambda= lambda.MLE, nx=200, ny=200, )
    #   surface(fsp.p, type="I")

    #   plot.new()
    #   fsp.p2<- predictSurfaceSE(fsp2)
    #   surface(fsp.p2, type="C")
    }

    return(out)

  }


  # ------------------------
  # ------------------------


  if ("gstat" %in% methods){
    require(gstat)
    require(sp)

    xy = as.data.frame(xy/out$stmv_internal_scale )
    names(xy) =  c("plon", "plat" ) # arbitrary

    co = out$distance_cutoff / out$stmv_internal_scale
    vEm = try( variogram( z~1, locations=~plon+plat, data=xy, cutoff=co, width=co/nbreaks, cressie=TRUE ) ) # empirical variogram
    if (inherits(vEm, "try-error") ) return(NULL)
    vEm$dist0 = vEm$dist * out$stmv_internal_scale
    vMod0 = vgm(psill=2/3*out$varZ, model="Mat", range=1, nugget=out$varZ/3, kappa=1/2 ) # starting model parameters
    vFitgs =  try( fit.variogram( vEm, vMod0, fit.kappa =TRUE, fit.sills=TRUE, fit.ranges=TRUE ) ) ## gstat's kappa is the Bessel function's "nu" smoothness parameter
    if (inherits(vFitgs, "try-error") )  return(NULL)

    scale = matern_phi2phi( mRange=vFitgs$range[2], mSmooth=vFitgs$kappa[2], parameterization_input="gstat", parameterization_output="stmv" ) * out$stmv_internal_scale

    out$gstat = list( fit=vFitgs, vgm=vEm, nu=vFitgs$kappa[2], phi=scale,
        varSpatial=vFitgs$psill[2], varObs=vFitgs$psill[1]  )  # gstat::"range" == range parameter == phi
    localrange = matern_phi2distance( phi=out$gstat$phi, nu=out$gstat$nu, cor=range_correlation  )
    out$gstat$phi_ok = ifelse( out$gstat$phi < out$distance_cutoff*0.99, TRUE, FALSE )

    if (plotdata) {
      xub = max(out$distance_cutoff) *1.25
      plot.new()
      plot(vEm, model=vFitgs, add=T)
      plot.new()
      plot( gamma ~ dist0, data=out$gstat$vgm, xlim=c(0,xub),
           ylim=c(0,max(out$gstat$vgm$gamma)*1.1), col="blue", pch=20 )
      abline( h=out$gstat$varSpatial + out$gstat$varObs )
      abline( h=out$gstat$varObs )
      abline( v=localrange )
      abline (v=0)
      x = seq( 0, xub, length.out=100 )
      acor = stmv_matern( x, mRange=out$gstat$phi, mSmooth=out$gstat$nu  )
      acov = out$gstat$varObs + out$gstat$varSpatial * (1- acor)
      lines( acov~x , col="red" )

      if (0) {
        # looks at the predictions
        gs <- gstat(id = "z", formula = z~1, locations=~plon+plat, data=xy, maxdist=distance_cutoff, nmin=10, force=TRUE, model=vFitgs )
        # variogram of residuals
        data(meuse.grid)
        meuse.grid$plon = meuse.grid$x
        meuse.grid$plat = meuse.grid$y

        preds <- predict(gs, newdata=meuse.grid )
        spplot(preds)

      }
    }
    return(out)
  }


  # # -------------------------



  if ("geoR" %in% methods) {
    gc() ## crashes often ...
    require( geoR )
    vEm = try( variog( coords=xy/out$stmv_internal_scale, data=z, uvec=nbreaks, max.dist=out$distance_cutoff/out$stmv_internal_scale ) )
    if  (inherits(vEm, "try-error") )  return(NULL)
    vEm$u0 = vEm$u * out$stmv_internal_scale
    gc()

    vMod = try( variofit( vEm, nugget=0.5*out$varZ, kappa=0.5, cov.model="matern",
      ini.cov.pars=c(0.5*out$varZ, 1 ),  limits = pars.limits( phi=c(0.1, 3), kappa=c(0.1, 5), sigmasq=c(0, out$varZ*1.25) ),
      fix.kappa=FALSE, fix.nugget=FALSE, max.dist=out$distance_cutoff/out$stmv_internal_scale, weights="cressie" ) )
      # kappa is the smoothness parameter , also called "nu" by others incl. RF
    gc()

    retry=FALSE
    if  (inherits(vMod, "try-error") ) retry=TRUE
    if (vMod$kappa < 0.1 | vMod$kappa > 5) retry =TRUE

    if (retry) {
      vMod = try( variofit( vEm, nugget=0.5*out$varZ, kappa=0.5, cov.model="matern",
        ini.cov.pars=c(0.5*out$varZ, 1 ),  limits = pars.limits( phi=c(0.1, 3), kappa=c(0.1, 5), sigmasq=c(0, out$varZ*1.25) ),
        fix.kappa=TRUE, fix.nugget=FALSE, max.dist=out$distance_cutoff/out$stmv_internal_scale, weights="cressie" ) )
        # kappa is the smoothness parameter , also called "nu" by others incl. RF
      gc()
    }

    if  (inherits(vMod, "try-error") )  return(NULL)
    scale = matern_phi2phi( mRange=vMod$cov.pars[2], mSmooth=vMod$kappa,
      parameterization_input="geoR", parameterization_output="stmv" ) * out$stmv_internal_scale

    out$geoR = list( fit=vMod, vgm=vEm, model=vMod,
            varSpatial= vMod$cov.pars[1], varObs=vMod$nugget,
            nu=vMod$kappa,  phi=scale )
    out$geoR$phi_ok = ifelse( out$geoR$phi < out$distance_cutoff*0.99, TRUE, FALSE )

    if (plotdata) {
      # not rescaled ...
      xub = max(out$distance_cutoff, out$geoR$vgm$u ) *1.25
      # plot.new()
      # plot(vEm)
      # lines(vMod)
      # plot.new()
      # plot( out$geoR$vgm )
      plot.new()
      plot( out$geoR$vgm$v ~ out$geoR$vgm$u0, pch=20 ,
           xlim=c(0,xub), ylim=c(0, max(out$geoR$varSpatial + out$geoR$varObs, out$varZ)) )
      abline( h=out$geoR$varSpatial + out$geoR$varObs)
      abline( h=out$geoR$varObs )

      x = seq( 0, xub, length.out=100 )
      acor = geoR::matern( x, phi=vMod$cov.pars[2]* out$stmv_internal_scale, kappa=vMod$kappa  )
      acov = out$geoR$varObs +  out$geoR$varSpatial * (1-acor)  ## geoR is 1/2 of gstat and RandomFields gamma's
      lines( acov ~ x , col="orange" )

      # check to see if the same as above .. yes!
      acor2 = stmv_matern( x, mRange=out$geoR$phi, mSmooth=out$geoR$nu  )
      acov2 = out$geoR$varObs +  out$geoR$varSpatial * (1-acor2)  ## geoR is 1/2 of gstat and RandomFields gamma's
      lines( acov2 ~ x , col="red" )
      localrange = matern_phi2distance( phi=out$geoR$phi, nu=out$geoR$nu, cor=range_correlation  )
      abline( v=localrange, col="orange" )

    }
    gc()

    return(out)

  }


  # -------------------------
  # ------------------------


  if ("RandomFields" %in% methods) {
    require( RandomFields ) ## max likilihood
    # rownames( xy) = 1:nrow(xy)  # seems to require rownames ...
   # RandomFields:  Cov(h) = v * Orig(A*h/s) ; s=scale, h=dist, A=aniso, v=variance, Original model (data scale)

   # where nu > 0 and K_nu is the modified Bessel function of second kind and distance r >= 0 between two pointsd
   # The Matern covariance model is given by: C(h) = v * phi(A*h/s).
   #  Cov(r) = 2^{1- nu} Gamma(nu)^{-1} (sqrt{2nu} r)^nu K_nu(sqrt{2nu} r)
   # "phi" = scale / sqrt{2nu} ... confirmed (jc, sep 2017)

   # RFoptions(
   #   allowdistanceZero=TRUE,
    #  modus_operandi="precise", #‘"careless"’,‘"sloppy"’, ‘"easygoing"’, ‘"normal"’, ‘"precise"’,        ‘"pedantic"’, ‘"neurotic"’
   #   bin_dist_factor=out$distance_cutoff/2,
      #bins=nbreaks,
      #critical=TRUE,
   #   approx_zero=0.05, #  Value below which a correlation is considered to be essentially zero.
   #   spConform=TRUE # FALSE is faster
    #)

    model = RMmatern( nu=NA, var=NA, scale=NA) + RMnugget(var=NA)

    o = RFfit(model, x=xy/out$stmv_internal_scale, data=z, allowdistanceZero=TRUE,  modus_operandi=modus_operandi )
    oo=summary(o)

    scale = matern_phi2phi( mRange=oo$param["value", "matern.s"], mSmooth=oo$param["value", "matern.nu"], parameterization_input="RandomFields", parameterization_output="stmv" ) * out$stmv_internal_scale

    out$RandomFields = list ( fit=o, vgm=o[2], model=oo,
              varSpatial=oo$param["value", "matern.var"],
              varObs=oo$param["value", "nugget.var"],
              phi=scale,
              nu=oo$param["value", "matern.nu"], # RF::nu == geoR:: kappa (bessel smoothness param)
              error=NA )
    out$RandomFields$phi_ok = ifelse( out$RandomFields$phi < out$distance_cutoff*0.99, TRUE, FALSE )

    if (plotdata) {
      xub = max(out$distance_cutoff, out$RandomFields$vgm@centers) *1.25
      plot.new()
      py = as.vector(out$RandomFields$vgm@empirical)
      px = out$RandomFields$vgm@centers * out$stmv_internal_scale
      plot(  py ~ px, pch=20, ylim=c(0, max(py)*1.1 )  )
      abline( h=out$RandomFields$varSpatial + out$RandomFields$varObs  )
      abline( h=out$RandomFields$varObs )
      localrange = matern_phi2distance( phi=out$RandomFields$phi, nu=out$RandomFields$nu, cor=range_correlation  )
      abline( v=localrange )

      x = seq( 0, xub, length.out=100 )
      acor = stmv_matern( x, mRange=out$RandomFields$phi, mSmooth=out$RandomFields$nu  )
      acov = out$RandomFields$varObs  + out$RandomFields$varSpatial*(1- acor)
      lines( acov~x , col="red" )
    }
    return(out)

  }


  # -------------------------
  # ------------------------



  if ("geoR.ML" %in% methods) {

    require( geoR )
    vEm = try( variog( coords=xy/out$stmv_internal_scale, data=z, uvec=nbreaks, max.dist=out$distance_cutoff/out$stmv_internal_scale ) )
    if  (inherits(vEm, "try-error") )  return(NULL)
    v0 = try( variofit( vEm, nugget=0.5*out$varZ, kappa=0.5, cov.model="matern",
      ini.cov.pars=c(0.5*out$varZ, 1) ,
      fix.kappa=FALSE, fix.nugget=FALSE, max.dist=out$distance_cutoff/out$stmv_internal_scale, weights="cressie" ) )
      # kappa is the smoothness parameter , also called "nu" by others incl. RF
    if  (inherits(v0, "try-error") )  return(NULL)

    # maximum likelihood method does not work well with Matern
     vMod = try( likfit( coords=as.matrix(xy/out$stmv_internal_scale), data=z, cov.model="matern", ini.cov.pars=v0$cov.pars,
      fix.kappa=FALSE, fix.nugget=FALSE, lik.method = "REML" ) )
# try to add this to make it go faster:  parscale =  c(range=0.1, shape=1, boxcox=1, nugget=out$varZ/100 )
#  and then add to likfit call: control=list(parscale=parscale )
# where they control resolution of optimization

    if  (inherits(vMod, "try-error") ) return (NULL)
    scale = matern_phi2phi( mRange=vMod$cov.pars[2], mSmooth=vMod$kappa, parameterization_input="geoR", parameterization_output="stmv" )*out$stmv_internal_scale
    out$geoR.ML = list( fit=vMod, vgm=vEm, model=vMod,
            varSpatial= vMod$cov.pars[1], varObs=vMod$nugget,
            nu=vMod$kappa,  phi=scale )
    out$geoR.ML$phi_ok = ifelse( out$geoR.ML$phi < out$distance_cutoff*0.99, TRUE, FALSE )

    if (plotdata) {
      xub = max(out$distance_cutoff, out$geoR.ML$vgm$u) *1.25
      # not rescaled ...
      plot.new()
      plot(vEm)
      lines(vMod)
      plot.new()
      plot( out$geoR.ML$vgm )
      plot.new()
      plot( out$geoR.ML$vgm$v ~ I(out$geoR.ML$vgm$u*out$stmv_internal_scale), pch=20 ,
           xlim=c(0,xub), ylim=c(0, max(out$geoR.ML$varSpatial + out$geoR.ML$varObs, out$varZ)) )
      abline( h=out$geoR.ML$varSpatial + out$geoR.ML$varObs)
      abline( h=out$geoR.ML$varObs )
      # abline( v=vMod$practicalRange*out$stmv_internal_scale, col="green" )
      localrange = matern_phi2distance( phi=out$geoR.ML$phi, nu=out$geoR.ML$nu, cor=range_correlation  )
      abline( v=localrange, col="orange" )

      x = seq( 0, xub, length.out=100 )
      acor = geoR::matern( x, phi=vMod$cov.pars[2]*out$stmv_internal_scale, kappa=vMod$kappa  )
      acov = out$geoR.ML$varObs +  out$geoR.ML$varSpatial * (1-acor)  ## geoR.ML is 1/2 of gstat and RandomFields gamma's
      lines( acov ~ x , col="orange" )

      acor2 = stmv_matern( x, mRange=out$geoR.ML$phi, mSmooth=out$geoR.ML$nu  )
      acov2 = out$geoR.ML$varObs +  out$geoR.ML$varSpatial * (1-acor2)  ## geoR.ML is 1/2 of gstat and RandomFields gamma's
      lines( acov2 ~ x , col="red" )

    }

    return(out)

  }


  # -------------------------
  # ------------------------



  if ("spBayes" %in% methods) {
    # note spBayes::phi = 1/ gstat::phi
    require(spBayes)
    library(MBA)
    phibounds = c(1/100, 5 ) ## approximate
    nubounds = c(0.01, 4 )# Finley et al 2007 suggest limiting this to (0,2)
    # Finley, Banerjee Carlin suggest that nu > 2 are indistinguishable .. identifiability problems cause slow solutions
    n.samples = 5000
    starting = list( phi=median(phibounds), sigma.sq=0.51, tau.sq=0.51, nu=1.1  ) # generic start
    tuning   = list( phi=starting$phi/12, sigma.sq=starting$sigma.sq/12, tau.sq=starting$tau.sq/12, nu=starting$nu/12 ) # MH variance to get acceptance rante bet 30-40%
    priors   = list(
      beta.flat = TRUE,
      phi.unif  = phibounds,
      sigma.sq.ig = c(5, 0.5),  # inverse -gamma (shape, scale):: scale identifies centre; shape higher = more centered .. assuming tau ~ sigma
      tau.sq.ig = c(5, 0.5),  # inverse gamma (shape, scale) :: invGamma( 3,1) -> modal peaking < 1, center near 1, long tailed
      nu.unif = nubounds
    )

    model = spLM( z ~ 1, coords=as.matrix(xy)/out$stmv_internal_scale, starting=starting, tuning=tuning, priors=priors, cov.model="matern",
      n.samples=n.samples, verbose=TRUE )

    burn.in <- 0.2*n.samples

    ##recover beta and spatial random effects
    m.1 <- spRecover(model, start=burn.in )

    u = apply(m.1$p.theta.recover.samples, 2, mean)
    scale = matern_phi2phi( mRange=u[["phi"]], mSmooth=u[["nu"]],
      parameterization_input="spBayes", parameterization_output="stmv" ) * out$stmv_internal_scale

    out$spBayes = list( model=model, recover=m.1,
      varSpatial=u[["sigma.sq"]], varObs=u[["tau.sq"]],
      phi=scale, nu=u["nu"] )  # output using geoR nomenclature
    out$spBayes$phi_ok = ifelse( out$spBayes$phi < out$distance_cutoff*0.99, TRUE, FALSE )

    if (plotdata) {
      plot.new()
      x = seq( 0,  max(out$distance_cutoff ) *1.1, length.out=100 )
      acor = stmv_matern( x, mRange=scale, mSmooth=u[["nu"]] )
      acov = u[["tau.sq"]] +  u[["sigma.sq"]] * (1- acor)  ## geoR is 1/2 of gstat and RandomFields gamma's
      plot( acov ~ x , col="orange", type="l", lwd=2, ylim=c(0,max(acov)*1.1) )
      abline( h=u[["tau.sq"]] + u[["sigma.sq"]]  )
      abline( h=u[["tau.sq"]] )
      abline( h=0 )
      abline( v=0 )
      localrange=matern_phi2distance( phi=scale, nu=u[["nu"]], cor=range_correlation  )
      abline( v=localrange )

      if(0) {
        round(summary(m.1$p.theta.recover.samples)$quantiles,2)
        round(summary(m.1$p.beta.recover.samples)$quantiles,2)
        m.1.w.summary <- summary(mcmc(t(m.1$p.w.recover.samples)))$quantiles[,c(3,1,5)]

        plot(z, m.1.w.summary[,1], xlab="Observed w", ylab="Fitted w",
            xlim=range(z), ylim=range(m.1.w.summary), main="Spatial random effects")
        arrows(z, m.1.w.summary[,1], z, m.1.w.summary[,2], length=0.02, angle=90)
        arrows(z, m.1.w.summary[,1], z, m.1.w.summary[,3], length=0.02, angle=90)
        lines(range(z), range(z))

        plot.new()
        obs.surf <-   mba.surf(cbind(xy, z), no.X=100, no.Y=100, extend=T)$xyz.est
        image(obs.surf, xaxs = "r", yaxs = "r", main="Observed response")
        points(xy)
        contour(obs.surf, add=T)
      }
    }

    return(out)

  }


  # -------------------------
  # ------------------------


  if ("inla" %in% methods){
    require(INLA)
    require(lattice)

    inla.setOption(scale.model.default = TRUE)  # better numerical performance of IGMRF models and less dependnence upon hyperpriors

    xys = data.frame( xy / out$stmv_internal_scale)
    locs0  = as.matrix( xys )
    xys$b0 = 1  # intercept for inla

    M0.domain = inla.nonconvex.hull( locs0 )
    MESH = inla.mesh.2d (
      loc=locs0, # locations of data points
      boundary = M0.domain,
       max.n = c(400, 20)
    )

    alpha = 2 # -> alpha-1 == nu (inla fixes it at 1,2,or 3)
    SPDE = inla.spde2.matern( MESH, alpha=alpha )
    spatial.field <- inla.spde.make.index('spatial.field', n.spde=SPDE$n.spde )

    # projection matrix A to translate from mesh nodes to data nodes
    A = inla.spde.make.A( mesh=MESH, loc=locs0 )

    # data stack for occurence (PA)
    Z = inla.stack(
        tag="data",
        data=list( z=z ) ,
        A=list(A, 1 ),
        effects=list( spatial.field=spatial.field, xys )  # b0 is the intercept
    )

    RES <- inla(  z ~ 0 + b0 + f( spatial.field, model=SPDE ), family=family,
        data=inla.stack.data(Z),
        # control.compute=list(dic=TRUE),
        control.results=list(return.marginals.random=TRUE ),
        # control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
        # control.fixed = list(expand.factor.strategy='inla') ,
        control.predictor=list(A=inla.stack.A(Z), compute=TRUE, link=1 ) ,
        # control.inla = list( h=1e-4, tolerance=1e-10),
        # control.inla=list(strategy="laplace", npoints=21, stencil=7 ) ,
        verbose = FALSE
    )

    oo = inla.spde2.result(RES, "spatial.field", SPDE, do.transf=TRUE)

    inames = c( "mode", "mean", "sd", "quant0.025", "quant0.25", "quant0.5",  "quant0.75", "quant0.975", "low", "high" )

    # Range parameter .. ie, sqrt(8)/exp(oo$summary.log.kappa$mean)
    im = oo$marginals.range.nominal[[1]]
    iRange = c( mode=inla.mmarginal( im ), inla.zmarginal( im, silent=TRUE ), as.data.frame(inla.hpdmarginal( 0.95, im )) )

    # "Spatial variance/error ('partial sill variance')"
    im = oo$marginals.variance.nominal[[1]]
    iVar =  c( mode=inla.mmarginal( im ), inla.zmarginal( im, silent=TRUE ), as.data.frame(inla.hpdmarginal( 0.95, im )) )

    im = oo$marginals.kappa[[1]]
#    im[,"x"]  = sqrt(8*)/im[,"x"]    # better to convert here
    iKappa =  c( mode=inla.mmarginal( im ), inla.zmarginal( im, silent=TRUE ), as.data.frame(inla.hpdmarginal( 0.95, im ) ) )

    # tau
    im = oo$marginals.tau[[1]]
    iTau =  c( mode=inla.mmarginal( im ), inla.zmarginal( im, silent=TRUE ), as.data.frame(inla.hpdmarginal( 0.95, im ) ) )

    ## Non-spatial ("observation") error ('nugget variance')
    iprec = grep ( "Precision.*observ.*", names(RES$marginals.hyperpar), ignore.case=TRUE )
    im = inla.tmarginal( function(x) {1/x}, RES$marginals.hyperpar[[ iprec ]] )
    iNugget =  c( mode=inla.mmarginal( im ), inla.zmarginal( im, silent=TRUE ), as.data.frame(inla.hpdmarginal( 0.95, im ) ) )

    inla.summary = as.matrix( rbind( iKappa, iTau, iRange, iVar, iNugget ) )
    rownames( inla.summary) = c( "kappa", "tau", "range", "spatial error", "observation error" )
    colnames( inla.summary) = inames

    scale = matern_phi2phi( mRange=inla.summary[["kappa","mean"]], mSmooth=alpha-1,
      parameterization_input="inla", parameterization_output="stmv" ) * out$stmv_internal_scale

    warning( "Need to check this parameterization for scale.. does not seem right")

    out$inla = list(summary=inla.summary,
      mesh=MESH, res=RES, # range.inla.practical=inla.summary[["range","mean"]]*out$stmv_internal_scale , # about 85.7% cor
      varSpatial=inla.summary[["spatial error","mean"]] ,
      varObs=inla.summary[["observation error","mean"]] ,
      phi=scale, nu=alpha-1, error=NA )

    # out$inla$range = matern_phi2distance( phi=out$inla$phi, nu=out$inla$nu, cor=range_correlation  )
    out$inla$phi_ok = ifelse( out$inla$phi < out$distance_cutoff*0.99, TRUE, FALSE )

    if (plotdata) {
      x = seq( 0,  max(out$distance_cutoff) *1.25, length.out=100 )
      svar =  out$inla$varObs + out$inla$varSpatial * (1- stmv_matern( x, mRange=out$inla$phi, mSmooth=out$inla$nu  ))
      plot( svar~x, type="l" )
      abline( h=out$inla$varObs + out$inla$varSpatial )
      abline( h=out$inla$varObs )
      # abline( v=out$inla$range.inla.practical, col="green"  ) # more like 87%
      localrange = matern_phi2distance( phi=out$inla$phi, nu=out$inla$nu, cor=range_correlation  )
      abline( v=localrange, col="red"  )
    }

    return(out)

  }


  # -------------------------
  # ------------------------


  if ("geostatsp_poisson" %in% methods){
    # poisson!  https://www.jstatsoft.org/article/view/v063i12; page 10
    require(geostatsp)

    n_min = 100
    rx = range(xy[,1])
    ry = range(xy[,2])
    origin = c( rx[1], ry[1] )
    resx = diff(rx) / n_min
    resy = diff(ry) / n_min
    res = c(min(resx, resy), min(resx, resy))
    xy_blocked = array_map( "xy->2", xy, res=res, origin=origin ) * res[1]
    xy_blocked = xy_blocked / out$stmv_internal_scale

    m = tapply( X=z, INDEX=list(xy_blocked[,1], xy_blocked[,2]),
        FUN = function(w) {mean(w, na.rm=TRUE)},
        simplify=TRUE )

    xyz = as.data.frame( as.table (m) )
    xyz[,1] = as.numeric(as.character( xyz[,1] ))
    xyz[,2] = as.numeric(as.character( xyz[,2] ))
    xyz = xyz[ which( is.finite( xyz[,3] )) ,]
    xyz = rasterFromXYZ( xyz )

    fit = geostatsp::lgcp(
      formula = ~ 1,
      data=xyz,
      priorCI = list(range=c(0.01,0.8)*n_min, sd=c(0.01, 2.0)*sqrt(out$varZ) ),  # priors are specified as 95% bounds
      buffer=n_min*0.05,
      grid=xyz
    )

    out$geostatsp = list( mesh=MESH, res=fit,
      # range.inla.practical=inla.summary[["range","mean"]]*out$stmv_internal_scale , # actaully it is p=0.14
      varSpatial=inla.summary[["spatial error","mean"]] ,
      varObs=inla.summary[["observation error","mean"]] ,
      phi=scale, nu=alpha-1, error=NA )

    out$geostatsp$phi_ok = ifelse( out$geostatsp$phi < out$distance_cutoff*0.99, TRUE, FALSE )

    if (plotdata) {
      x = seq( 0,  max(out$distance_cutoff) *1.25, length.out=100 )
      svar =  out$geostatsp$varObs + out$geostatsp$varSpatial * (1- stmv_matern( x, mRange=out$geostatsp$phi, mSmooth=out$geostatsp$nu  ))
      plot( svar~x, type="l" )
      abline( h=out$geostatsp$varObs + out$geostatsp$varSpatial )
      abline( h=out$geostatsp$varObs )
      # abline( v=out$geostatsp$range.inla.practical, col="green"  ) # more like 87%
      localrange = matern_phi2distance( phi=out$geostatsp$phi, nu=out$geostatsp$nu, cor=range_correlation  )
      abline( v=localrange, col="red"  )
    }

    return(out)
  }


  # -------------------------
  # ------------------------


  if ("bayesx" %in% methods){
    library("R2BayesX")
    # by default, bayesx fixes nu=1.5  , see: bayesx.term.options( bs="kr", method="REML" )
    # phi = max(distance) / const, such that Corr(distance=const) = 0.001;
    nu = 0.5
    xys = as.data.frame(xy/out$stmv_internal_scale)
    names(xys) =  c("plon", "plat" ) # arbitrary

    fm <- bayesx( z ~ sx(plon, plat, bs="kr" ), family=family$family, method="REML", data =xys )
    # fm <- bayesx( z ~ sx(plon, plat,  bs="kr" ), family=family$family, method="HMCMC", data =xy )
    # ?bayesx.control
    # warning( "BayesX documentation is not clear if rho is the scale parameter. This seems ok, but should do some more testing." )

    logout = fm$bayesx.run[[1]]
    logout = logout[ grepl( "Parameter rho:", logout) ]
    res = regexpr("Parameter rho: ", logout)
    logout = substring( logout, res[1], res[1]+35 )
    rho = as.numeric( unlist( strsplit( logout, ":" ) )[2] )

    scale = matern_phi2phi( mRange=rho, mSmooth=nu,
      parameterization_input="bayesx", parameterization_output="stmv" ) * out$stmv_internal_scale

    out$bayesx = list( fit=fitted(fm), model=fm,
        varSpatial = fm$smooth.hyp[,"Variance"] ,
        varObs = fm$fixed.effects[1,"Std. Error"] ,
        nu =nu, phi=scale )

    out$bayesx$phi_ok = ifelse( out$bayesx$phi < out$distance_cutoff*0.99, TRUE, FALSE )

    if(0){
      plot( fm, term = "sx(plon,plat)", image=TRUE, contour=TRUE )
      # summary(fm)
      # lattice::levelplot( out ~ plon+plat, data=k, aspect="iso" )
      x = seq( 0, out$distance_cutoff, length.out=100 )
      acor = stmv_matern( x, mRange=out$bayesx$phi, mSmooth=out$bayesx$nu  )
      acov = out$bayesx$varObs  + out$bayesx$varSpatial*(1- acor)
      plot( acov~x , col="red", type="l", xlim=c(0, out$distance_cutoff ),
        ylim=c(0,(out$bayesx$varSpatial + out$bayesx$varObs) *1.25) )
      abline( h=out$bayesx$varSpatial + out$bayesx$varObs  )
      abline( h=out$bayesx$varObs )
      localrange = matern_phi2distance( phi=out$bayesx$phi, nu=out$bayesx$nu, cor=range_correlation  )
      abline( v=localrange )
    }

    return(out)
  }


  # -------------------------
  # ------------------------


  if ("jags.exponential" %in% methods){
    require(rjags)
    require(jagsUI)
    # assume nu = 1 (due to identifiability issues)

    print( "Slow ... 7.5 min for meuse test data" )

    jagsmodel = paste0("
    model{
      for(i in 1:N){
        y[i] ~ dnorm(mu[i], prec)
        mu[i] = beta0 + errorSpatial[i]
        muSpatial[i] = 0
      }
      prec = 1.0/ (tausq + sigmasq )
      invCOVAR = inverse(COVAR)
      errorSpatial ~ dmnorm( muSpatial, invCOVAR)
      for(i in 1:N) {
        COVAR[i,i] = sigmasq
        for(j in 1:(i-1)) {
          COVAR[i,j] = sigmasq * exp(-( DIST[i,j]/phi ))
          COVAR[j,i] = COVAR[i,j]
        }
      }
      tausq = 1/tausq.inv
      tausq.inv ~ dgamma(0.1,0.1)
      sigmasq = 1/sigmasq.inv
      sigmasq.inv ~ dgamma(2,1)
      phi ~ dgamma(1,0.1)
      beta0 ~ dnorm(0,0.0001)
    } ")

    fn = tempfile()
    cat( jagsmodel, file=fn )

    distances = abs( as.matrix(  dist( xy/ out$stmv_internal_scale ) ) )

    Data = list( N=length(z), DIST=distances, y=z )
    fit = jagsUI::jags(data=Data,
       parameters.to.save=c("phi", "sigmasq", "tausq"),
       model.file=fn,
       n.iter=1000,
       n.chains=3,
       n.burnin=100,
       n.thin=5,
       parallel=TRUE,
       DIC=FALSE)

   if (0) {
     summary(fit)
     plot(fit)
     gelman.plot(fit$samples)
    # geweke.plot(fit$samples)
    #update(fit, n.iter=2000, n.thin=20 )
      acf( fit$sims.list$phi)
      acf( fit$sims.list$sigmasq)
      acf( fit$sims.list$tausq)

    }
    print (summary(fit))

    out$jags = list(
      fit = fit, model=jagsmodel,
      phi =  fit$summary["phi", "mean"] * out$stmv_internal_scale,
      sigmasq = fit$summary["sigmasq", "mean"] ,
      tausq = fit$summary["tausq", "mean"]
    )
    localrange = matern_phi2distance( phi=out$jags$phi, nu=0.5, cor=range_correlation  )
    out$jags$phi_ok = ifelse( out$jags$phi < out$distance_cutoff*0.99, TRUE, FALSE )

    return(out)
  }


  # -------------------------

  if ("TMB" %in% methods){

  }


  # --------------------

  if ( grepl("stan", methods) ) {

    library(rstan)
    rstan_options(auto_write = TRUE)
    options(mc.cores = parallel::detectCores())

    inp = list()
    inp$varY = var(z)
    inp$Y = z / inp$varY
    inp$N = length( inp$Y)
    inp$Np = 0
    inp$X=rep(0,inp$N)
    inp$eps = 1e-6
    inp$COVFN=1
    inp$sigma_sq0 = 0.8
    inp$tau_sq0 = 0.2
    inp$dist=as.matrix( dist( as.data.frame( xy/out$stmv_internal_scale )) )
#    inp$vgm_dist_max = max(abs(inp$dist))
#    inp$dist = inp$dist / inp$vgm_dist_max

    if (is.null(stanmodel)) stanmodel=rstan::stan_model( model_code= "

      functions{

        matrix matern_covariance( int N, matrix dist, real phi, real sigma_sq, real tau_sq, real eps, int COVFN) {
          matrix[N,N] S;
          real dist_phi;
          real sqrt3;
          real sqrt5;
          sqrt3=sqrt(3.0);
          sqrt5=sqrt(5.0);

          if (COVFN==1) { // exponential == Matern nu=1/2 , (p=0; nu=p+1/2)
            for(i in 1:(N-1)){
              for(j in (i+1):N){
                dist_phi = fabs(dist[i,j])/phi;
                S[i,j] = sigma_sq * exp(- dist_phi );
            }}


          } else if (COVFN==2) { // Matern nu= 3/2 covariance
            for(i in 1:(N-1)){
              for(j in (i+1):N){
               dist_phi = fabs(dist[i,j])/phi;
               S[i,j] = sigma_sq * (1 + sqrt3 * dist_phi) * exp(-sqrt3 * dist_phi);
            }}


          } else if (COVFN==3) { // Matern nu=5/2 covariance
            for(i in 1:(N-1)){
              for(j in (i+1):N){
                dist_phi = fabs(dist[i,j])/phi;
                S[i,j] = sigma_sq * (1 + sqrt5 *dist_phi + 5* pow(dist_phi,2)/3) * exp(-sqrt5 *dist_phi);
            }}

          } else if (COVFN==4) { // Matern as nu->Inf become Gaussian (aka squared exponential cov)
            for(i in 1:(N-1)){
              for(j in (i+1):N){
                dist_phi = fabs(dist[i,j])/phi;
                S[i,j] = sigma_sq * exp( -pow(dist_phi,2)/2 ) ;
            }}
          }


          for(i in 1:(N-1)){
          for(j in (i+1):N){
            S[j,i] = S[i,j];  // fill upper triangle
          }}

          // create diagonal: nugget(nonspatial) + spatial variance +  eps ensures positive definiteness
          for(i in 1:N) S[i,i] = sigma_sq + tau_sq + eps;
          return(S)   ;
        }

      }

    // ---------------
      data {
        int<lower=1> N;
        // int<lower=0> Np; //number of data locations to predict
        vector[N] Y; // observations
        // vector[N] X; // covariates (predictors)
        matrix[N, N] dist; //distances between points
        real eps;
        // real sigma_sq0;
        // real tau_sq0;
        int<lower=1,upper=4> COVFN;  // Choice of Matern covariance function:
        // 1:nu=1/2 (exponential); 2:nu=3/2; 3:nu=5/2; 4:nu->Inf (gaussian)
      }

     // ---------------
      transformed data{
        vector[N] zeros;
        zeros = rep_vector(0, N);
      }

     // ---------------
      parameters{
        real<lower=0, upper=10> sigma_sq;
        real<lower=0, upper=10> tau_sq;
        real<lower=0, upper=10> phi;
        vector [N] spatialError;
      }

     // ---------------
      transformed parameters{
      }

     // ---------------
      model{
        matrix[N,N] S; // Covariance
        S = matern_covariance( N, dist, phi, sigma_sq, tau_sq, eps, COVFN );
        Y ~ multi_normal_cholesky( spatialError, cholesky_decompose(S) ) ;
        tau_sq ~ cauchy( 0, 0.5) ;
        sigma_sq ~ cauchy( 0, 0.5) ;
        phi ~ normal( 1, 1 ) ;
      }

     // ---------------
      generated quantities{
      //  vector[N] y_pred;
      //  for(i in 1:N)  y_pred[i] = (beta+plogitobs[i]);
      //  for(i in 1:Np) y_pred[N+i] = (beta+plogitpreds[i]);
      }
    "

  )

    # f = optimizing(stanmodel, data=inp, hessian=FALSE )
    # optimizing(stanmodel, data=inp, hessian=FALSE, tol_rel_obj=1e6, algorithm="LBFGS"  )
    # optimizing(stanmodel, data=inp, hessian=FALSE, tol_rel_obj=1e6, algorithm="BFGS"  )

    # f = vb(stanmodel, data=inp)

    f = sampling(stanmodel, data=inp, iter=1000, chains=3)
      # warmup = 200,          # number of warmup iterations per chain
      # control =. list(adapt_delta = 0.9),
      # # refresh = 500,          # show progress every 'refresh' iterations
      # iter = 1000,            # total number of iterations per chain
      # chains = 5,             # number of Markov chains
      # cores = 5              # number of cores (using 2 just for the vignette)
    # f = optimizing(stanmodel, data=inp, hessian=FALSE )

    plot(f)
    print(f)
    traceplot(f)

    # extract samples
    # m2 = as.array(f)
    pred=rstan::extract(f)
    phii = mean(pred$phi) * out$stmv_internal_scale
    rnge = matern_phi2distance( phi=phii, nu=0.5, cor=range_correlation  )
    out$stan$phi_ok = ifelse( out$stan$phi < out$distance_cutoff*0.99, TRUE, FALSE )

    # prob=apply(pred,2,function(x) I(length(x[x>0.10])/length(x) > 0.8)*1)
    return( "Method not yet finished .. speed is the issue (though faster than JAGS) and ML/VB methods are unstable")
  }


  # -------------------------
  # -------------------------

  if ("LaplacesDemon" %in% methods){

    require(LaplacesDemonCpp)
    require( FastGP)

    Data = list(
      eps = 1e-6,
      N = length(z),  # required for LaplacesDemon
      DIST=as.matrix(dist( xy, diag=TRUE, upper=TRUE)), # distance matrix between knots
      Y=z/var(z),
      varZ=var(z)
    )
    Data$DIST = Data$DIST / out$stmv_internal_scale
    Data$mon.names = c( "LP", paste0("yhat[",1:Data$N,"]" ) )
    Data$parm.names = as.parm.names(list(mu=rep(0, Data$N), tausq=0, sigmasq=0, phi=0, nu=0 ))
    Data$pos = list(
      mu = grep("mu", Data$parm.names),
      tausq = grep("tausq", Data$parm.names),
      sigmasq = grep("sigmasq", Data$parm.names),
      phi = grep("phi", Data$parm.names),
      nu = grep("nu", Data$parm.names)
    )
    Data$PGF = compiler::cmpfun(
      function(Data) {
        #initial values .. get them near the center of mass
        mu = rnorm(Data$N, 0, sqrt(Data$varZ))
        va = rcauchy (2, 0, 0.5) # variances
        phi = rnorm (1, 1, 1) # scale
        nu = rnorm(1, 1, 1)
        return( c( mu, va, phi, nu ))
    } )
    Data$Model = compiler::cmpfun(
      function(parm, Data) {
        mu = parm[Data$pos$mu]
        tausq = parm[Data$pos$tausq] = LaplacesDemonCpp::interval(parm[Data$pos$tausq], 0, 10 )
        sigmasq = parm[Data$pos$sigmasq]= LaplacesDemonCpp::interval(parm[Data$pos$sigmasq], 0, 10 )
        phi = parm[Data$pos$phi]= LaplacesDemonCpp::interval(parm[Data$pos$phi], Data$eps, 10 )
        nu = parm[Data$pos$nu] = LaplacesDemonCpp::interval(parm[Data$pos$nu], 0.01, 10 )

        S = sigmasq * (1-stmv_matern( Data$DIST, mRange=phi, mSmooth=nu ) )  # spatial covariance
        diag(S) = diag(S) + tausq + Data$eps
        if ( !is.positive.definite(S)) {
        #   cat("correlation matrix is not positive definite, adding a bit of noise ...\n")
           S = as.positive.definite(S)
        }

        tausq.prior = dcauchy(tausq, 0, 0.5, log=TRUE) # 0-1.55 range
        sigmasq.prior = dcauchy(sigmasq, 0, 0.5, log=TRUE)
        phi.prior = dnorm(phi, 1, 1, log=TRUE)
        nu.prior = dnorm(nu, 1, 1, log=TRUE)
        LL = FastGP::rcpp_log_dmvnorm(S, mu, Data$Y, istoep=FALSE)

#        LL = dmvn( Data$Y, mu, S, log=TRUE) ## Log Likelihood
        LP = sum(LL, sigmasq.prior, tausq.prior, phi.prior, nu.prior) ### Log-Posterior
        yhat = c( FastGP::rcpp_rmvnorm(1,S,mu) ) # rmvn( 1, mu, S )
        Modelout = list(LP=LP, Dev=-2*LL, Monitor=c(LP, yhat), yhat=yhat, parm=parm)
        return(Modelout)
    } )

    parm0=Data$PGF(Data)

    f = LaplaceApproximation(Data$Model, Data=Data, parm=parm0 )


    if (plotdata) {

      f = LaplacesDemon(Data$Model, Data=Data, Initial.Values=as.initial.values(f), Fit.LA=f,
        Iterations=1000, Thinning=100, Status=1000, Covar=f$Covar, CPUs=8 )

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
      plot( Data$Y ~ m[, "Mean"]  )


    }

    out$LaplacesDemon = list( fit=f, vgm=NA, model=Data$Model,
      varSpatial=f$Summary2["sigmasq", "Mean"] ,
      varObs=f$Summary2["tausq", "Mean"],
      nu=f$Summary2["nu", "Mean"],
      phi =f$Summary2["phi", "Mean"] *  out$stmv_internal_scale
    )   ## need to check parameterization...

    out$LaplacesDemon$phi_ok = ifelse( out$LaplacesDemon$phi < out$distance_cutoff*0.99, TRUE, FALSE )

   # print( out$LaplacesDemon )


    if (plotdata) {
      plot.new()
      x = seq( 0, out$LaplacesDemon$range * 1.2, length.out=100 )
      svar =  out$LaplacesDemon$varObs + out$LaplacesDemon$varSpatial * (1-stmv_matern( x, out$LaplacesDemon$phi, out$LaplacesDemon$nu  ))
      plot( svar~x, type="l", ylim=c(0, max(svar)) )
      abline( h=out$LaplacesDemon$varObs + out$LaplacesDemon$varSpatial )
      abline( h=out$LaplacesDemon$varObs )
      localrange = matern_phi2distance( phi=out$LaplacesDemon$phi, nu=out$LaplacesDemon$nu , cor=range_correlation)
      abline( v=localrange, col="red"  )
    }

    return(out)

  }



      if (0) {
        # this is leftover code ... tucked here in case something is useful for RandomFields-based approaches

        require( RandomFields ) ## max likilihood
        z = c(v)

       # RandomFields:  Cov(h) = v * Orig(A*h/s) ; s=scale, h=dist, A=aniso, v=variance, Original model (data scale)
       # where nu > 0 and K_nu is the modified Bessel function of second kind and distance r >= 0 between two pointsd
       # The Matern covariance model is given by: C(h) = v * phi(A*h/s).
       #  Cov(r) = 2^{1- nu} Gamma(nu)^{-1} (sqrt{2nu} r)^nu K_nu(sqrt{2nu} r)
        # Aniso_matern = Aniso=matrix(NA, nc=3, nr=2)    #Aniso_exp = Aniso=matrix(NA, nc=2, nr=2)
        if (RFmodel=="exponential") {
          model = ~ 1@RMfixed(beta=NA) + RMexp( var=NA, scale=NA ) + RMnugget(var=NA )

        }
        if (RFmodel=="matern") {
          model = ~ 1@RMfixed(beta=NA) + RMmatern( var=NA, scale=NA, nu=NA ) + RMnugget(var=NA)
        }




        microbenchmark::microbenchmark(
        {
          pmod =   glm(z~1, family=poisson(link=log), offset=rep( prod(reso), length(z) ) )
          dta = residuals( pmod )

          fit = RFfit(model,
            data = dta,
            x=xy[,1],
            y=xy[,2],
            # use_spam=TRUE,
            # matrix_methods=c(0,1,2), #  0:Choleskey decomposition, 1:SVD,  2:spam (sparse matrix algorithm),  3:OR, 2  4:LU
            # printlevel=5,
            spConform=FALSE, # FALSE is faster
            allowdistanceZero=TRUE,
            # bin_dist_factor=1/3, # empirical variogram is calculated up the distance ‘bin_dist_factor’ times (maximum distance
            # bins=13,
            max_neighbours=2000,
            # practicalrange=2,
            # use_naturalscaling=TRUE,
            modus_operandi="easygoing"  #‘"careless"’,‘"sloppy"’, ‘"easygoing"’, ‘"normal"’, ‘"precise"’, ‘"pedantic"’, ‘"neurotic"’
          )
        }, times=1
        )

        fitsummary = summary(fit)

        vgm=fit[2]
        varObs=fitsummary$param["value", "nugget.var"]

        if (RFmodel=="exponential") {
          varSpatial=fitsummary$param["value", "exp.var"]
          phi=fitsummary$param["value", "exp.s"]
          localrange = geoR::practicalRange("exp", phi=phi, kappa=nu, correlation=0.05 )
        }

        if (RFmodel=="matern") {
          varSpatial=fitsummary$param["value", "matern.var"]
          phi=(fitsummary$param["value", "matern.s"] )/(sqrt(fitsummary$param["value", "matern.nu"]*2) )  # RF parameterizes as scale/sqrt(nu*2)

##*** phi_Wikipedia   =  phi_geostatsp / 2
##*** phi_RandomFields = phi_geostatsp / 2
##*** phi_geoR      = phi_geostatsp / sqrt(8*nu)
##*** scale_inla = alpha =  sqrt(8*nu) / phi_geostatsp
##*** ["range for space"]_inla = phi_geostatsp * length(grid)

          if (!exists("nu")) nu = 1/2
          if ( any( grepl("matern.nu", colnames( fitsummary$param))) ) nu=fitsummary$param["value", "matern.nu"]
          localrange = geoR::practicalRange("matern", phi=phi, kappa=nu, correlation=0.05 )
        }

#      R>         varObs=fitsummary$param["value", "nugget.var"]
# R> varObs
# [1] 0.1924
# R>     varSpatial=fitsummary$param["value", "matern.var"]
# R>           phi=(fitsummary$param["value", "matern.s"] )/(sqrt(fitsummary$param["value", "matern.nu"]*2) )  # RF parameterizes as scale/sqrt(nu*2)
# R>
# R>           if (!exists("nu")) nu = 1/2
# R>           if ( any( grepl("matern.nu", colnames( fitsummary$param))) ) nu=fitsummary$param["value", "matern.nu"]
# R>           range = geoR::practicalRange("matern", phi=phi, kappa=nu, correlation=0.05 )
# R> range
# [1] 45581
# R> phi
# [1] 12117
# R> varSpatial
# [1] 815.6
# R> nu
# [1] 0.8637


        plotdata=FALSE
        if (plotdata) {
          plot.new()
          py = as.vector(vgm$variogram)
          px = vgm@centers
          plot(  py ~ px, pch=20 )
          abline( h=varSpatial + varObs  )
          abline( h=varObs )
          abline( v=localrange )

          x = seq( 0, max(px ), length.out=100 )
          acor = geoR::matern( x, phi=phi, kappa=nu  )
          acov = varObs  + varSpatial*(1- acor)
          lines( acov~x , col="red" )
        }


      }



}
