
stmv_singlepass_fft = function( ip=NULL, p, debugging=FALSE, runoption="default", eps=1e-9, ... ) {
  #\\ core function to interpolate (model and predict) in parallel
  #\\ unfold all main subroutines to max speed and reduce memory usage

  if (0) {
    # for debugging  runs ..
    eps = 1e-9
    currentstatus = stmv_statistics_status( p=p )
    p = parallel_run( p=p, runindex=list( locs=sample( currentstatus$todo )) )
    ip = 1:p$nruns
    debugging=TRUE
    runoption="default"
  }

  if (0) {
    # testing and debugging

    loadfunctions( c("aegis", "stmv"))
    RLibrary(c ("fields", "MBA", "geoR") )

    if (0) {
      XYZ = stmv_test_data( datasource="swiss" )
      mz = log( XYZ$rain )
      mm = lm( mz ~ 1 )
      XYZ$z = residuals( mm)
      XYZ=XYZ[c("x","y","z")]
    }

    if (0) {
      XYZ = stmv_test_data( datasource="meuse" )
      XYZ$z = log(XYZ$elev)
      XYZ=XYZ[, c("x","y","z") ]
    }

  }

  # ---------------------
  # deal with additional passed parameters
  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast=TRUE))
  if (length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable
  if (exists( "libs", p)) suppressMessages( RLibrary( p$libs ) )
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns


  #---------------------
  # data for modelling
  E = stmv_error_codes()
  S = stmv_attach( p$storage.backend, p$ptr$S )
  Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )
  Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
  Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )
  Y = stmv_attach( p$storage.backend, p$ptr$Y )
  Yi = stmv_attach( p$storage.backend, p$ptr$Yi )  # initial indices of good data

  if (p$nloccov > 0) Ycov = stmv_attach( p$storage.backend, p$ptr$Ycov )
  if ( exists("TIME", p$variables) ) Ytime = stmv_attach( p$storage.backend, p$ptr$Ytime )


  # pre-calculate indices and dim for data to use inside the loop
  dat_names = unique( c(  p$variables$Y, p$variable$LOCS, p$variables$local_all,  "weights") )  # excludes p$variables$TIME
    if ( exists("TIME", p$variables)) {
      dat_names = c(dat_names, p$variables$TIME)
    }

  # unless it is an explicit covariate and not a seasonal component there is no need for it
  # .. prediction grids create these from a time grid on the fly
  dat_nc = length( dat_names )
  iY = which(dat_names== p$variables$Y)
  ilocs = which( dat_names %in% p$variable$LOCS )
  # iwei = which( dat_names %in% "weights" )

  if (p$nloccov > 0) {
    icov = which( dat_names %in% p$variables$local_cov )
    icov_local = which( p$variables$COV %in% p$variables$local_cov )
  }
  if (exists("TIME", p$variables)) {
    ti_cov = setdiff(p$variables$local_all, c(p$variables$Y, p$variables$LOCS, p$variables$local_cov ) )
    itime_cov = which(dat_names %in% ti_cov)
  }

  nip = length(ip)
  if (nip < 100) {
    nlogs = 3
  } else {
    nlogs = p$nlogs
  }
  logpoints  = ip[ floor( seq( from=10, to=(nip-10), length.out=nlogs ) ) ]

  ndata_index = match( "ndata", p$statsvars )
  rsquared_index = match("rsquared", p$statsvars )

  stmv_nmins = sort( unique( c(1, p$stmv_nmin_downsize_factor) * p$stmv_nmin ) , decreasing=FALSE )
  stmv_nmax = p$stmv_nmax

# main loop over each output location in S (stats output locations)
  for ( iip in ip ) {

    if ( iip %in% logpoints )  currentstatus = stmv_logfile(p=p, flag= paste("Interpolation", runoption) )
    Si = p$runs[ iip, "locs" ]

    # print( paste("index =", iip, ";  Si = ", Si ) )
    # if (debugging) print( paste("index =", iip, ";  Si = ", Si ) )
    if ( Sflag[Si] == E[["complete"]] ) next()


    # obtain indices of data locations withing a given spatial range, optimally determined via variogram
    # find data nearest Sloc[Si,] and with sufficient data
    ndata = 0

    for ( ntarget in stmv_ntarget ) {
      for ( stmv_distance_cur in stmv_distances )  {
        print(stmv_distance_cur)
        # obtain indices of data locations withing a given spatial range, optimally determined via variogram
        # faster to take a block .. but easy enough to take circles ...
        U = which(
          (abs( Sloc[Si,1] - Yloc[Yi[],1] ) <= stmv_distance_cur ) &
          (abs( Sloc[Si,2] - Yloc[Yi[],2] ) <= stmv_distance_cur )
        )
        ndata = length(U)

        if (ndata < p$stmv_nmin) {

          Sflag[Si] = E[["insufficient_data"]]

        } else if (ndata > p$stmv_nmax) {

          # try to trim
          if ( exists("TIME", p$variables)) {
            iU = stmv_discretize_coordinates( coo=cbind(Yloc[U,], Ytime[U]), ntarget=p$stmv_nmax, minresolution=p$minresolution, method="thin" )
          } else {
            iU = stmv_discretize_coordinates( coo=Yloc[U,], ntarget=p$stmv_nmax, minresolution=p$minresolution, method="thin" )
          }
          ndata = length(iU)
          # print(ndata)

          ntarget = floor( (p$stmv_nmax+p$stmv_nmin)/2 )

          if (ndata < p$stmv_nmin) {
            # return some data
            uu = setdiff( 1:length(U), iU)
            nuu = length(uu)
            iMore = uu[ .Internal( sample( nuu, {ntarget - ndata}, replace=FALSE, prob=NULL)) ]
            U = U[c(iU, iMore)]
            ndata = ntarget
            iMore = nuu = uu = NULL

          } else if (ndata > p$stmv_nmax) {

            # force via a random subsample
            U = U[iU]
            U = U[ .Internal( sample( length(U), ntarget, replace=FALSE, prob=NULL)) ] # simple random
            ndata = ntarget

          } else {

            U = U[iU]

          }

        } else  if (ndata <= p$stmv_nmax & ndata >= p$stmv_nmin) {
          # all good .. nothing to do
        }
        iU = NULL

        yi = Yi[U]

        if ( is.null( yi ) ) next()
        ndata = length(yi)
        if ( ndata >= ntarget ) break()  # innermost loop
      }
      if ( ndata >= ntarget ) break() # middle loop
    }

    if (ndata < p$stmv_nmin ) {
      Sflag[Si] = E[["insufficient_data"]]
      if (debugging) print( paste("index =", iip, ";  insufficient data"  ) )
      next()   #not enough data
    }


    YiU = Yi[U]
    U = NULL
    d1 = d2 = NULL

    # NOTE: this range is a crude estimate that averages across years (if any) ...

    o = NULL
    gc()

    xy = Yloc[yi,]
    z = Y[yi,]

    distance_cutoff=stmv_distance_cur
    discretized_n = stmv_distance_cur / p$pres
    nbreaks=p$stmv_variogram_nbreaks
    range_correlation=p$stmv_range_correlation

    out = list()

    out$Ndata = length(yi)
    out$varZ = var( Y[yi,], na.rm=TRUE )  # this is the scaling factor for semivariance .. diving by sd, below reduces numerical floating point issues

    # reduce scale of distances to reduce large number effects
    out$range_crude = sqrt( diff(range(xy[,1]))^2 + diff(range(xy[,2]))^2) / 4  #initial scaling distance
    out$stmv_internal_scale = matern_distance2phi( out$range_crude, nu=0.5 )  # the presumed scaling distance to make calcs use smaller numbers

    # if max dist not given, make a sensible choice using exponential variogram as a first estimate
    if ( is.na(distance_cutoff) ) {
      out$distance_cutoff = out$range_crude * 1.5
    } else {
      out$distance_cutoff = distance_cutoff
    }

    zmin = min( Y[yi,], na.rm=TRUE )

    xy = xy + out$stmv_internal_scale * runif(2*out$Ndata, -1e-9, 1e-9) # add a small error term to prevent some errors in GRMF methods


    # spatial discretization
    XYZ = stmv_discretize_coordinates(coo=xy, z=Y[yi,], discretized_n=discretized_n, method="aggregate", FUNC=mean, na.rm=TRUE)
    names(XYZ) =  c("plon", "plat", "z" ) # arbitrary

    VGM = stmv_variogram_fft( xyz=XYZ, nx=discretized_n, ny=discretized_n, nbreaks=nbreaks )  # empirical variogram by fftw2d
    uu = which( (VGM$vgm$distances < 0.9*max(VGM$vgm$distances) ) & is.finite(VGM$vgm$sv) )

    fit = try( stmv_variogram_optimization( vx=VGM$vgm$distances[uu], vg=VGM$vgm$sv[uu], plotvgm=plotdata,
      stmv_internal_scale=out$stmv_internal_scale, cor=range_correlation  ))

    if ( !inherits(fit, "try-error") ) {
      out$fft = fit$summary
      if (exists("phi", out$fft)) {
        if (is.finite(out$fft$phi)) {
          out$fft$phi_ok = ifelse( out$fft$phi < out$distance_cutoff*0.99, TRUE, FALSE )
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

    yi = NULL

    gc()

    if ( is.null(o)) {
      Sflag[Si] = E[["variogram_failure"]]
      if (debugging) print( paste("index =", iip, ";  o is null"  ) )
      next()
    }
    if ( inherits(o, "try-error")) {
      Sflag[Si] = E[["variogram_failure"]]
      if (debugging) print( paste("index =", iip, ";  o has try error"  ) )
      next()
    }
    if ( !exists(p$stmv_variogram_method, o)) {
      Sflag[Si] =  E[["variogram_failure"]]
      if (debugging) print( paste("index =", iip, ";  o does not have a solution"  ) )
      next()
    }
    om  = o[[p$stmv_variogram_method]] # save stats

    statvars_scale = c(
      sdTotal =sqrt( o$varZ),
      sdSpatial = sqrt(om$varSpatial) ,
      sdObs = sqrt(om$varObs),
      phi = om$phi,
      nu = om$nu,
      localrange = matern_phi2distance( phi=om$phi, nu=om$nu, cor=p$stmv_range_correlation ),
      ndata=ndata
    )

    S[Si,match( names(statvars_scale), p$statsvars )] = statvars_scale


    # temporal
    if (p$stmv_dimensionality =="space") {
       # nothing to do
    }

    if (p$stmv_dimensionality =="space-year") {

      if (0) {

        # annual ts, seasonally centered and spatially
        ar_timerange = NA
        ar_1 = NA

        pac = res$predictions[ pac_i, ]
        pac$dyr = pac[, p$variables$TIME] - trunc(pac[, p$variables$TIME] )
        piid = which( zapsmall( pac$dyr - p$dyear_centre) == 0 )
        pac = pac[ piid, c(p$variables$TIME, "mean")]
        pac = pac[ order(pac[,p$variables$TIME]),]
        if (length(piid) > 5 ) {
          ts.stat = NULL
          ts.stat = try( stmv_timeseries( pac$mean, method="fft" ) )
          if (!is.null(ts.stat) && !inherits(ts.stat, "try-error") ) {
            ar_timerange = ts.stat$quantilePeriod
            if (all( is.finite(pac$mean))) {
              afin = which (is.finite(pac$mean) )
              if (length(afin) > 5 && var( pac$mean, na.rm=TRUE) > eps ) {
                ar1 = NULL
                ar1 = try( ar( pac$mean, order.max=1 ) )
                if (!inherits(ar1, "try-error")) {
                  if ( length(ar1$ar) == 1 ) {
                    ar_1 = ar1$ar
                  }
                }
              }
            }
            if ( !is.finite(out[["ar_1"]]) ) {
              ar1 = try( cor( pac$mean[1:(length(piid) - 1)], pac$mean[2:(length(piid))], use="pairwise.complete.obs" ) )
              if (!inherits(ar1, "try-error")) ar_1 = ar1
            }
          }

          ### Do the logistic model here ! -- if not already done ..
          if (!exists("ts_K", out)) {
            # model as a logistic with ts_r, ts_K, etc .. as stats outputs

          }
        }
        pac=piid=NULL
        pac_i=NULL
        statsvars_time =c(
          ar_timerange= ar_timerange,
          ar_1 = ar_1
        )
         # save stats
        S[Si, match( names(statsvars_time), p$statsvars )] = statsvars_time
      }
    }

    if (p$stmv_dimensionality=="space-year-season")  {

      if (0) {

        # annual ts, seasonally centered and spatially
        ar_timerange = NA
        ar_1 = NA

        pac = res$predictions[ pac_i, ]
        pac$dyr = pac[, p$variables$TIME] - trunc(pac[, p$variables$TIME] )
        piid = which( zapsmall( pac$dyr - p$dyear_centre) == 0 )
        pac = pac[ piid, c(p$variables$TIME, "mean")]
        pac = pac[ order(pac[,p$variables$TIME]),]
        if (length(piid) > 5 ) {
          ts.stat = NULL
          ts.stat = try( stmv_timeseries( pac$mean, method="fft" ) )
          if (!is.null(ts.stat) && !inherits(ts.stat, "try-error") ) {
            ar_timerange = ts.stat$quantilePeriod
            if (all( is.finite(pac$mean))) {
              afin = which (is.finite(pac$mean) )
              if (length(afin) > 5 && var( pac$mean, na.rm=TRUE) > eps ) {
                ar1 = NULL
                ar1 = try( ar( pac$mean, order.max=1 ) )
                if (!inherits(ar1, "try-error")) {
                  if ( length(ar1$ar) == 1 ) {
                    ar_1 = ar1$ar
                  }
                }
              }
            }
            if ( !is.finite(out[["ar_1"]]) ) {
              ar1 = try( cor( pac$mean[1:(length(piid) - 1)], pac$mean[2:(length(piid))], use="pairwise.complete.obs" ) )
              if (!inherits(ar1, "try-error")) ar_1 = ar1
            }
          }

          ### Do the logistic model here ! -- if not already done ..
          if (!exists("ts_K", out)) {
            # model as a logistic with ts_r, ts_K, etc .. as stats outputs

          }

        }
        pac=piid=NULL
        pac_i=NULL

        statsvars_time =c(
          ar_timerange= ar_timerange,
          ar_1 = ar_1
        )
        # save stats
        S[Si, match( names(statsvars_time), p$statsvars )] = statsvars_time
      }
    }

    Sflag[Si] = E[["complete"]]  # mark as complete without issues

    if (debugging) {
      print( paste("index =", iip, ";  Sflag = ", names(E)[match(Sflag[Si], E)]  ) )
    }


    if ( is.null(o)) {
      Sflag[Si] = E[["variogram_failure"]]
      # if (debugging) print( paste("index =", iip, ";  o is null"  ) )
      next()
    }
    if ( inherits(o, "try-error")) {
      Sflag[Si] = E[["variogram_failure"]]
      # if (debugging) print( paste("index =", iip, ";  o has try error"  ) )
      next()
    }
    if ( !exists(p$stmv_variogram_method, o)) {
      Sflag[Si] =  E[["variogram_failure"]]
      # if (debugging) print( paste("index =", iip, ";  o does not have a solution"  ) )
      next()
    }
    om  = o[[p$stmv_variogram_method]] # save stats

    statvars_scale = c(
      sdTotal =sqrt( o$varZ),
      sdSpatial = sqrt(om$varSpatial) ,
      sdObs = sqrt(om$varObs),
      phi = om$phi,
      nu = om$nu,
      localrange = matern_phi2distance( phi=om$phi, nu=om$nu, cor=p$stmv_range_correlation ),
      ndata=ndata
    )
    S[Si,match( names(statvars_scale), p$statsvars )] = statvars_scale




    distance_limits = range( c(p$pres*3,  p$stmv_distance_scale ) )   # for range estimate





    if (!is.finite(vg$phi)) vg$phi = median(distance_limits)
    if ( vg$phi < distance_limits[1] )  vg$phi = distance_limits[1]
    if ( vg$phi > distance_limits[2] )  vg$phi = distance_limits[2]

    # as range is now set, the following becomes fixed
    if (is.null(ii)) {
      localrange = matern_phi2distance( phi=vg$phi, nu=vg$nu, cor=p$stmv_range_correlation )
      ii = which(
        {abs( Sloc[Si,1] - Sloc[,1] ) <= localrange} &
        {abs( Sloc[Si,2] - Sloc[,2] ) <= localrange}
      )
    }

    # obtain indices of data locations withing a given spatial range, optimally determined via variogram
    # faster to take a block .. but easy enough to take circles ... trim off corners ..
    nu_redo = FALSE
    if (!is.finite(vg$nu)) {
      nu_redo = TRUE
    } else {
      if (vg$nu < 0.25 | vg$nu > 4 )  {
        nu_redo = TRUE
      }
    }

    if (nu_redo) {
      vg$flag = "variogram_failure"
      nu_median =  median( S[ii, match("nu", p$statsvars)], na.rm=TRUE )
      if (!is.finite(nu_median)) {
        vg$nu = 0.5
      } else {
        if (nu_median < 0.25 | nu_median > 4 )  {
          vg$nu = 0.5
        } else {
          vg$nu = nu_median
        }
      }
    }


    phi_redo = FALSE
    dl = distance_limits/sqrt(8*vg$nu)
    if (!is.finite(vg$phi)) {
      phi_redo = TRUE
    } else {
      if (vg$phi < dl[1] | vg$phi > dl[2] )  {
        phi_redo = TRUE
      }
    }

    if (phi_redo) {
      vg$flag = "variogram_failure"
      phi_median =  median( S[ii, match("phi", p$statsvars)], na.rm=TRUE )
      if (is.finite(phi_median)) {
        vg$phi = phi_median
      } else {
        vg$phi = vg$range/sqrt(8*vg$nu)  # approx rule rule in inla for nu=1 (alpha=2)
      }
    }


    if (!is.finite(vg$varSpatial)) {
      vg$flag = "variogram_failure"
      varSP_median = median( S[ii, match("sdSpatial", p$statsvars)], na.rm=TRUE )^2
      if (is.finite(varSP_median)) vg$varSpatial = varSP_median
    }


    if (!is.finite(vg$varObs)) {
      vg$flag = "variogram_failure"
      varObs_median = median( S[ii, match("sdObs", p$statsvars)], na.rm=TRUE )^2
      if (is.finite(varObs_median)) vg$varObs = varObs_median
    }


    if (!is.finite(vg$varTotal)) {
      vg$flag = "variogram_failure"
      varTot_median = median( S[ii, match("sdTotal", p$statsvars)], na.rm=TRUE )^2
      if (is.finite(varTot_median)) vg$varTotal = varTot_median
    }


    if (exists("stmv_rangecheck", p)) {
      if (p$stmv_rangecheck=="paranoid") {
        if ( vg$flag %in% c("variogram_range_limit", "variogram_failure") ) {
          Sflag[Si] = E[[vg$flag]]
          if (runoption == "default" ) next()
        }
      }
    }

    localrange = vg$range
    ndata = vg$range

    useglobal = FALSE
    if (!is.finite( localrange ) ) useglobal =TRUE
    if (!is.finite( ndata ) ) useglobal =TRUE
    if (vg$flag =="variogram_failure") useglobal =TRUE
    if (useglobal) {
      vg = list(
        nu=0.5,
        phi=xx
      )
    }

    if (runoption=="boostdata") localrange = matern_phi2distance( phi=vg$phi, nu=vg$nu, cor=p$stmv_range_correlation_boostdata )

    U = stmv_select_data( p=p, Si=Si, localrange=localrange )

    if ( Sflag[Si] != E[["todo"]] ) {
      if (exists("stmv_rangecheck", p)) {
        if (p$stmv_rangecheck=="paranoid") {
          if ( Sflag[Si] %in% c( E[["variogram_range_limit"]], E[["variogram_failure"]]) ) {
            U = NULL
            # if (debugging) message("Error: stmv_rangecheck paranoid")
            if (runoption == "default" ) next()
          }
        }
      }
    }

    if (is.null( U )) next()

    # last check
    ndata = length(U)
    if (ndata < p$stmv_nmin) next()

    # if here then there is something to do
    # NOTE:: U are the indices of locally useful data



    # prep dependent data
    # reconstruct data for modelling (dat)
    dat = matrix( 1, nrow=ndata, ncol=dat_nc )
    dat[,iY] = Y[U] # these are residuals if there is a global model
    # add a small error term to prevent some errors when duplicate locations exist; localrange offsets to positive values
    dat[,ilocs] = Yloc[U,] + localrange * runif(2*ndata, -1e-6, 1e-6)

    if (p$nloccov > 0) dat[,icov] = Ycov[U, icov_local] # no need for other dim checks as this is user provided
    if (exists("TIME", p$variables)) dat[, itime_cov] = as.matrix(stmv_timecovars( vars=ti_cov, ti=Ytime[U,] ) )
    dat = as.data.frame(dat)
    names(dat) = dat_names

    # not sure if this is needed ?...
    if (p$stmv_local_modelengine %in% c("fft", "tps") ) {
      if ( exists("TIME", p$variables)) {
        dat[, p$variables$TIME] = Ytime[U,]
      }
    }

    # remember that these are crude mean/discretized estimates
    # if (debugging) {
    #   dev.new()
    #   # check data and statistical locations
    #   plot( Sloc[,], pch=20, cex=0.5, col="gray")
    #   points( Yloc[,], pch=20, cex=0.2, col="green")
    #   points( Yloc[U,], pch=20, cex=1, col="yellow" )
    #   points( Sloc[Si,2] ~ Sloc[Si,1], pch=20, cex=5, col="blue" )
    # }


    stmv_distance_prediction = localrange * p$stmv_distance_prediction_fraction # this is a half window km
    # construct prediction/output grid area ('pa')
    windowsize.half = floor(stmv_distance_prediction/p$pres) # convert distance to discretized increments of row/col indices; stmv_distance_prediction = 0.75* stmv_distance_statsgrid (unless overridden)


    # construct data (including covariates) for prediction locations (pa)
    pa = try( stmv_predictionarea( p=p, sloc=Sloc[Si,], windowsize.half=windowsize.half ) )
    if (is.null(pa)) {
      Sflag[Si] = E[["prediction_area"]]
      # if (debugging) message( Si )
      # if (debugging) message("Error: prediction grid ... null .. this should not happen")
      pa = NULL
      next()
    }
    if ( inherits(pa, "try-error") ) {
      pa = NULL
      Sflag[Si] = E[["prediction_area"]]
      # if (debugging) message("Error: prediction grid ... try-error .. this should not happen")
      next()
    }

    # if (debugging) {

    #   dev.new()
    #   # check that position indices are working properly
    #   Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
    #   Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
    #   Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )
    #   plot( Yloc[U,2] ~ Yloc[U,1], col="red", pch=".",
    #     ylim=range(c(Yloc[U,2], Sloc[Si,2], Ploc[pa$i,2]) ),
    #     xlim=range(c(Yloc[U,1], Sloc[Si,1], Ploc[pa$i,1]) ) ) # all data
    #   points( Yloc[U,2] ~ Yloc[U,1], col="green" )  # with covars and no other data issues
    #   points( Sloc[Si,2] ~ Sloc[Si,1], col="blue" ) # statistical locations
    #   # statistical output locations
    #   grids= spatial_grid(p, DS="planar.coords" )
    #   points( grids$plat[floor( (Sloc[Si,2]-p$origin[2])/p$pres) + 1]
    #         ~ grids$plon[floor( (Sloc[Si,1]-p$origin[1])/p$pres) + 1] , col="purple", pch=25, cex=5 )
    #   points( grids$plat[pa$iplat] ~ grids$plon[ pa$iplon] , col="cyan", pch=20, cex=0.01 ) # check on Proc iplat indexing
    #   points( Ploc[pa$i,2] ~ Ploc[ pa$i, 1] , col="black", pch=20, cex=0.7 ) # check on pa$i indexing -- prediction locations
    # }



    # model and prediction .. outputs are in scale of the link (and not response)
    # the following permits user-defined models (might want to use compiler::cmpfun )


  nu=vg$nu
  phi=vg$phi
  varObs=vg$varObs
  varSpatial=vg$varSpatial
  sloc=Sloc[Si,]
  distance=localrange


  #\\ twostep modelling time first as a simple ts and then spatial or spatio-temporal interpolation
  #\\ nu is the bessel smooth param

  # step 1 -- timeseries modelling
  # use all available data in 'dat' to get a time trend .. and assume it applies to the prediction area of interest 'pa'
     # some methods require a uniform (temporal with associated covariates) prediction grid based upon all dat locations

  px = dat # only the static parts .. time has to be a uniform grid so reconstruct below

  ids = array_map( "xy->1", px[, c("plon", "plat")], gridparams=p$gridparams ) # 100X faster than paste / merge
  todrop = which(duplicated( ids) )
  if (length(todrop>0)) px = px[-todrop,]
  ids = todrop=NULL

  # static vars .. don't need to look up
  tokeep = c(p$variables$LOCS )
  if (exists("weights", dat) ) tokeep = c(tokeep, "weights")
  if (p$nloccov > 0) {
    for (ci in 1:p$nloccov) {
      vn = p$variables$local_cov[ci]
      pu = stmv_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
      nts = ncol(pu)
      if ( nts==1 ) tokeep = c(tokeep, vn )
    }
  }
  px = px[ , tokeep ]
  px_n = nrow(px)
  nts = vn = NULL

  # add temporal grid
  if ( exists("TIME", p$variables) ) {
    px = cbind( px[ rep.int(1:px_n, p$nt), ],
                    rep.int(p$prediction.ts, rep(px_n, p$nt )) )
    names(px)[ ncol(px) ] = p$variables$TIME
    px = cbind( px, stmv_timecovars ( vars=p$variables$local_all, ti=px[,p$variables$TIME]  ) )
  }

  if (p$nloccov > 0) {
    # add time-varying covars .. not necessary except when covars are modelled locally
    for (ci in 1:p$nloccov) {
      vn = p$variables$local_cov[ci]
      pu = stmv_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
      nts = ncol(pu)
      if ( nts== 1) {
        # static vars are retained in the previous step
      } else if ( nts == p$ny )  {
        px$iy = px$yr - p$yrs[1] + 1 #yr index
        px[,vn] = pu[ cbind(px$i, px$iy) ]
       } else if ( nts == p$nt) {
        px$it = p$nw*(px$tiyr - p$yrs[1] - p$tres/2) + 1 #ts index
        px[,vn] = pu[ cbind(px$i, px$it) ]
      }
    } # end for loop
    nts = vn = NULL
  } # end if
  rownames(px) = NULL

  # print( "starting gam-timeseries mod/pred")
  ts_preds = NULL

  p$stmv_local_modelformula = p$stmv_local_modelformula_time

  if (p$stmv_twostep_time == "inla" ) ts_preds = stmv__inla_ts( p, dat, px )
  if (p$stmv_twostep_time == "glm" ) ts_preds = stmv__glm( p, dat, px )
  if (p$stmv_twostep_time == "gam" ) ts_preds = stmv__gam( p, dat, px )
  if (p$stmv_twostep_time == "bayesx" ) ts_preds = stmv__bayesx( p, dat, px )

  if (is.null( ts_preds)) {
    # TS modelling failure
    next()
  }

  # range checks
  rY = range( dat[,p$variables$Y], na.rm=TRUE)
  toosmall = which( ts_preds$predictions$mean < rY[1] )
  toolarge = which( ts_preds$predictions$mean > rY[2] )
  if (length(toosmall) > 0) ts_preds$predictions$mean[toosmall] =  rY[1]
  if (length(toolarge) > 0) ts_preds$predictions$mean[toolarge] =  rY[2]

  ts_preds_rsquared = ts_preds$stmv_stats$rsquared  # store for now until return call

  pxts = ts_preds$predictions
  rownames(pxts) = NULL
  ts_preds = NULL

  names(pxts)[which(names(pxts)=="mean")] = p$variables$Y
  names(pxts)[which(names(pxts)=="sd")] = paste(p$variables$Y, "sd", sep=".")

  # if(0){
  #     # debugging plots
  #     for (ti in 1:p$nt){
  #       xi = which( pxts[ , p$variables$TIME ] == p$prediction.ts[ti] )
  #       mbas = MBA::mba.surf( pxts[xi, c( p$variables$LOCS, p$variables$Y) ], 300, 300, extend=TRUE)$xyz.est
  #       image(mbas)
  #     }
  # }


  # step 2 :: spatial modelling .. essentially a time-space separable solution

  dat=pxts


  #\\ this is the core engine of stmv .. localised space (no-time) modelling interpolation
  #\\ note: time is not being modelled and treated independently
  #\\      .. you had better have enough data in each time slice
  #\\ first a low-pass filter as defined by p$stmv_lowpass_nu, p$stmv_lowpass_phi,
  #\\ then a simple covariance filter determined by nu,phi ;; fft (no lowpass)
  #\\ based upon fields::image.smooth and setup.image.smooth
  #\\ now that the local area of interest is stationary, we use local convolutions of highly autocorrelated (short-range)
  #\\ determined by p


  params = list(...)

  sdTotal = sd(dat[,p$variable$Y], na.rm=T)

  # nr .. x/plon
  # nc .. y/plat
  # pa_r = range(pa[,p$variables$LOCS[1]])
  # pa_c = range(pa[,p$variables$LOCS[2]])

  x_r = range(dat[,p$variables$LOCS[1]])
  x_c = range(dat[,p$variables$LOCS[2]])

  nr = floor( diff(x_r)/p$pres ) + 1
  nc = floor( diff(x_c)/p$pres ) + 1

  # final output grid
  x_locs = expand_grid_fast(
    seq( x_r[1], x_r[2], length.out=nr ),
    seq( x_c[1], x_c[2], length.out=nc )
  )
  attr( x_locs , "out.attrs") = NULL
  names( x_locs ) = p$variables$LOCS

  dat$mean = NA
  pa$mean = NA
  pa$sd = sdTotal  # this is ignored with fft

  dx = dy = p$pres

  nr2 = 2 * nr
  nc2 = 2 * nc


  rr = diff( x_r )  # system length scale
  rc = diff( x_c )

  # no of elements
  dr = rr/(nr-1)
  dc = rc/(nc-1)

  # approx sa associate with each datum
  sa = rr * rc
  d_sa = sa/nrow(dat) # sa associated with each datum
  d_length = sqrt( d_sa/pi )  # sa = pi*l^2  # characteristic length scale


  theta.Taper = vgm$distances[ find_intersection( vgm$ac, threshold=p$stmv_fft_taper_correlation ) ]
    theta.Taper = theta.Taper * p$stmv_fft_taper_fraction # fraction of the distance to 0 correlation; sqrt(0.5) = ~ 70% of the variability (associated with correlation = 0.5)


  # constainer for spatial filters
  grid.list = list((1:nr2) * dx, (1:nc2) * dy)
  # dgrid = as.matrix(expand.grid(grid.list))  # a bit slower
  dgrid = expand_grid_fast(  grid.list[[1]],  grid.list[[2]] )
  dimnames(dgrid) = list(NULL, names(grid.list))
  attr(dgrid, "grid.list") = grid.list
  grid.list = NULL

  center = matrix(c((dx * nr), (dy * nc)), nrow = 1, ncol = 2)

  mC = matrix(0, nrow = nr2, ncol = nc2)
  mC[nr, nc] = 1



  if (!exists("stmv_fft_filter",p) ) p$stmv_fft_filter="lowpass" # default in case of no specification

  if ( p$stmv_fft_filter == "lowpass") {
    theta = p$stmv_lowpass_phi
    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=theta, smoothness=p$stmv_lowpass_nu )
    sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2 * nc2)
    sp.covar.kernel = fftwtools::fftw2d(sp.covar) / fftwtools::fftw2d(mC)
  }

  if (p$stmv_fft_filter %in% c("matern") ) {
    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=phi, smoothness=nu )
    sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2 * nc2)
    sp.covar.kernel = fftwtools::fftw2d(sp.covar) / fftwtools::fftw2d(mC)
  }

  if (p$stmv_fft_filter == "lowpass_matern") {
    # both ..
    sp.covar.lowpass = stationary.cov( dgrid, center, Covariance="Matern", theta=p$stmv_lowpass_phi, smoothness=p$stmv_lowpass_nu )
    sp.covar.lowpass = as.surface(dgrid, c(sp.covar.lowpass))$z / (nr2 * nc2)
    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=phi, smoothness=nu )
    sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2 * nc2)
    sp.covar.kernel = (fftwtools::fftw2d(sp.covar.lowpass) / fftwtools::fftw2d(mC) ) * ( fftwtools::fftw2d(sp.covar)/ fftwtools::fftw2d(mC) )
    sp.covar = sp.covar.lowpass = NULL
  }

  if (p$stmv_fft_filter == "matern_tapered") {

    sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=phi, smoothness=nu,
      Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)
    sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2 * nc2)
    sp.covar.kernel = fftwtools::fftw2d(sp.covar) / fftwtools::fftw2d(mC)
    sp.covar = theta.Taper = NULL
  }

  if (p$stmv_fft_filter == "lowpass_matern_tapered") {
    sp.covar.lowpass = stationary.cov( dgrid, center, Covariance="Matern", theta=p$stmv_lowpass_phi, smoothness=p$stmv_lowpass_nu )
    sp.covar.lowpass = as.surface(dgrid, c(sp.covar.lowpass))$z / (nr2 * nc2)

    sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=phi, smoothness=nu,
      Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)
    sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2 * nc2)
    sp.covar.kernel = (fftwtools::fftw2d(sp.covar.lowpass) / fftwtools::fftw2d(mC) ) * ( fftwtools::fftw2d(sp.covar)/ fftwtools::fftw2d(mC) )
    sp.covar = sp.covar.lowpass = theta.Taper = NULL
  }


  if (p$stmv_fft_filter == "normal_kernel") {
      theta = matern_phi2distance( phi=phi, nu=nu, cor=p$stmv_range_correlation )
      xi = seq(-(nr - 1), nr, 1) * dx / theta
      yi = seq(-(nc - 1), nc, 1) * dy / theta
      dd = ((matrix(xi, nr2, nc2)^2 + matrix(yi, nr2, nc2, byrow = TRUE)^2))  # squared distances
      # double.exp: An R function that takes as its argument the _squared_
      # distance between two points divided by the bandwidth. The
      # default is exp( -abs(x)) yielding a normal kernel
      kk = double.exp(dd)
      mK = matrix(kk, nrow = nr2, ncol = nc2)
      sp.covar.kernel = fftwtools::fftw2d(mK) / fftwtools::fftw2d(mC)
      sp.covar.kernel = sp.covar.kernel / (nr2 * nc2)  # kernal weights
  }

  sp.covar = sp.covar.lowpass = dgrid = center =  NULL

  origin=c(x_r[1], x_c[1])
  res=c(p$pres, p$pres)

  zz = matrix(1:(nr*nc), nrow = nr, ncol = nc)
  mY = matrix(0, nrow = nr2, ncol = nc2)
  mN = matrix(0, nrow = nr2, ncol = nc2)

  for ( ti in 1:p$nt ) {

    if ( exists("TIME", p$variables) ) {
      xi   = which( dat[ , p$variables$TIME] == p$prediction.ts[ti] )
      pa_i = which( pa[, p$variables$TIME] == p$prediction.ts[ti] )
      if (length(xi) < 5 ) {
        # print( ti)
        next()
      }
    } else {
      xi   = 1:nrow(dat) # all data as p$nt==1
      pa_i = 1:nrow(pa)
    }


    # u = as.image(
    #   dat[xi, p$variables$Y],
    #   ind=as.matrix(array_map( "xy->2", coords=dat[xi, p$variables$LOCS], origin=origin, res=res )),
    #   na.rm=TRUE,
    #   nx=nr,
    #   ny=nc
    # )



    coo = as.matrix(array_map( "xy->2", coords=dat[xi, p$variables$LOCS], origin=origin, res=resolution ))
    yy = tapply( X=z, INDEX=list(coo[,1], coo[,2]), FUN = function(w) {mean(w, na.rm=TRUE)}, simplify=TRUE )
    nn = tapply( X=z, INDEX=list(coo[,1], coo[,2]), FUN = function(w) {length(w)}, simplify=TRUE )
    mY[as.numeric(dimnames(yy)[[1]]),as.numeric(dimnames(yy)[[2]])] = yy
    mN[as.numeric(dimnames(nn)[[1]]),as.numeric(dimnames(nn)[[2]])] = nn
    yy = nn = NULL
    coo = NULL



    # mY[1:nr,1:nc] = u$z
    mY[!is.finite(mY)] = 0

    #  Nadaraya/Watson normalization for missing values s
    # mN[1:nr,1:nc] = u$weights
    mN[!is.finite(mN)] = 0

    u =NULL

    fY = Re( fftwtools::fftw2d( sp.covar.kernel * fftwtools::fftw2d(mY), inverse = TRUE))[1:nr, 1:nc]  #real amplitudes
    fN = Re( fftwtools::fftw2d( sp.covar.kernel * fftwtools::fftw2d(mN), inverse = TRUE))[1:nr, 1:nc]
    Z = ifelse(fN > eps, (fY / fN), NA)
    fY = fN = NULL

    # low pass filter based upon a global nu,phi .. remove high freq variation
    Z_i = array_map( "xy->2", coords=pa[pa_i,p$variables$LOCS], origin=origin, res=res )

    # bounds check: make sure predictions exist
    Z_i_test = which( Z_i[,1]<1 | Z_i[,2]<1  | Z_i[,1] > nr | Z_i[,2] > nc )

    if (length(Z_i_test) > 0) {
      tokeep = zz[ Z_i[-Z_i_test,] ]
      if ( length(tokeep) > 0) pa$mean[pa_i[tokeep]] = Z[tokeep]
      tokeep = NULL
    } else {
      pa$mean[pa_i] = Z[Z_i]
    }

    # pa$sd[pa_i] = NA  ## fix as NA
    Z = Z_i = Z_i_test = NULL
  }

  stmv_stats = list( sdTotal=sdTotal, rsquared=NA, ndata=nrow(dat) ) # must be same order as p$statsvars

  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  out = ( list( predictions=pa, stmv_stats=stmv_stats ) )

  # override the spatial rsquared with the timeseries model as that is more meaningful (based upon the raw data, rather than a "boosted" series)
  out$stmv_stats$rsquared = ts_preds_rsquared

  if (0) {
    lattice::levelplot( mean ~ plon + plat, data=out$predictions[out$predictions[,p$variables$TIME]==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
    lattice::levelplot( mean ~ plon + plat, data=out$predictions, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
    for( i in sort(unique(out$predictions[,p$variables$TIME])))  print(lattice::levelplot( mean ~ plon + plat, data=out$predictions[out$predictions[,p$variables$TIME]==i,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
  }



    dat =  NULL
    pa  =  NULL

    if ( is.null(res)) {
      Sflag[Si] = E[["local_model_error"]]   # modelling / prediction did not complete properly
      # if (debugging) message("Error: local model error")
      next()
    }

    if ( inherits(res, "try-error") ) {
      Sflag[Si] =  E[["local_model_error"]]   # modelling / prediction did not complete properly
      res = NULL
      # if (debugging) message("Error: local model error")
      next()
    }

    if (!exists("predictions", res)) {
      Sflag[Si] =  E[["prediction_error"]]   # modelling / prediction did not complete properly
      res = NULL
      # if (debugging) message("Error: prediction error")
      next()
    }

    if (!exists("mean", res$predictions)) {
      Sflag[Si] =  E[["prediction_error"]]   # modelling / prediction did not complete properly
      res = NULL
      # if (debugging) message("Error: prediction error")
      next()
    }

    if (length(which( is.finite(res$predictions$mean ))) < 5) {
      Sflag[Si] =  E[["prediction_error"]]   # modelling / prediction did not complete properly
      res = NULL
      # if (debugging) {
      #   message("Error: prediction error")
      #   browser()
      # }
      next()  # looks to be a faulty solution
    }


    if (runoption=="default") {
      # update to rsquared and "ndata" in stats
      S[Si, rsquared_index] = res$stmv_stats$rsquared
      S[Si, ndata_index] = ndata
    }

    res$stmv_stats = NULL # reduce memory usage


    sf = try( stmv_predictions_update(p=p, preds=res$predictions ) )
    res = NULL

    if ( is.null(sf) ) {
      Sflag[Si] = E[["prediction_update_error"]]
      sf = NULL
      # if (debugging) message("Error: prediction update error .. is null")
      next()
    }
    if ( inherits(sf, "try-error") ) {
      Sflag[Si] = E[["prediction_update_error"]]
      sf = NULL
      # if (debugging) message("Error: prediction update error .. try-error")
      next()
    }
    if ( sf=="error" ) {
      Sflag[Si] = E[["prediction_update_error"]]
      sf = NULL
      # if (debugging) message("Error: prediction update error .. general")
      next()
    }


    # ----------------------
    # do last. it is an indicator of completion of all tasks
    # restarts would be broken otherwise
    Sflag[Si] = E[["complete"]]  # mark as complete without issues

  }  # end for loop


  return(NULL)




}
