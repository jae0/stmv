

stmv_interpolate = function( ip=NULL, p, debugging=FALSE, ... ) {
  #\\ core function to interpolate (model and predict) in parallel

  if (0) {
    # for debugging  runs ..
    currentstatus = stmv_db( p=p, DS="statistics.status" )
    p = parallel_run( p=p, runindex=list( locs=sample( currentstatus$todo )) )
    ip = 1:p$nruns
    debugging=TRUE
  }

  # ---------------------
  # deal with additional passed parameters
  if ( is.null(p) ) p=list()
  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast=TRUE))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable

  if (exists( "libs", p)) suppressMessages( RLibrary( p$libs ) )

  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns

  #---------------------
  # data for modelling
  S = stmv_attach( p$storage.backend, p$ptr$S )
  Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )
  E = stmv_error_codes()

  Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
  Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
  Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )

  Y = stmv_attach( p$storage.backend, p$ptr$Y )

  P = stmv_attach( p$storage.backend, p$ptr$P )
  Pn = stmv_attach( p$storage.backend, p$ptr$Pn )
  Psd = stmv_attach( p$storage.backend, p$ptr$Psd )

  if (p$nloccov > 0) Ycov = stmv_attach( p$storage.backend, p$ptr$Ycov )
  if ( exists("TIME", p$variables) ) Ytime = stmv_attach( p$storage.backend, p$ptr$Ytime )

  # misc intermediate calcs to be done outside of parallel loops

  # pre-calculate indices and dim for data to use inside the loop
  dat_names = unique( c(  p$variables$Y, p$variable$LOCS, p$variables$local_all,  "weights") )  # excludes p$variables$TIME
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

  local_fn = switch( p$stmv_local_modelengine,
    bayesx = stmv__bayesx,
    gaussianprocess2Dt = stmv__gaussianprocess2Dt,
    gam = stmv__gam,
    glm = stmv__glm,
    gstat = stmv__gstat,
    krige = stmv__krige,
    fft = stmv__fft,
    tps = stmv__tps,
    twostep = stmv__twostep,
    userdefined = p$stmv_local_modelengine_userdefined
  )

  nip = length(ip)
  if (nip < 100) {
    nlogs = 3
  } else {
    nlogs = p$nlogs
  }
  logpoints  = ip[ floor( seq( from=10, to=(nip-10), length.out=nlogs ) ) ]

  if (debugging) {
    nsavepoints = 3
    savepoints = sample(logpoints, nsavepoints)
  }

# main loop over each output location in S (stats output locations)
  for ( iip in ip ) {

    if ( iip %in% logpoints )  currentstatus = stmv_logfile(p=p)
    Si = p$runs[ iip, "locs" ]
    if ( Sflag[Si] != E[["todo"]] ) next()  # previously attempted .. skip
    print( paste("index =", iip, ";  Si = ", Si ) )

    # obtain indices of data locations withing a given spatial range, optimally determined via variogram

    W = try( stmv_subset_distance( Si=Si, p=p ) )

    if ( is.null(W) ) {
      Sflag[Si] = E[["insufficient_data"]]
      W = WA = NULL
      next()
    }
    if ( inherits(W, "try-error") ) {
      Sflag[Si] = E[["insufficient_data"]]
      W = WA = NULL
      next()
    }

    Sflag[Si] = W[["flag"]]  # update flags

    if ( Sflag[Si] == E[["insufficient_data"]] ) {
      W = WA = NULL
      next()
    }
    if ( Sflag[Si] != E[["todo"]] ) {
      if (exists("stmv_rangecheck", p)) {
        if (p$stmv_rangecheck=="paranoid") {
          W = WA = NULL
          next()
        }
      }
    }

    # if here then W exists and there is something to do
    # NOTE:: W[["U"]] are the indices of locally useful data
    # p$stmv_distance_prediction determine the data entering into local model construction

    # prep dependent data
    # reconstruct data for modelling (dat)
    dat = matrix( 1, nrow=W[["ndata"]], ncol=dat_nc )
    dat[,iY] = Y[W[["U"]]] # these are residuals if there is a global model
    # add a small error term to prevent some errors when duplicate locations exist; stmv_distance_cur offsets to positive values
    dat[,ilocs] = Yloc[W[["U"]],] + W[["stmv_distance_cur"]] * runif(2*W[["ndata"]], -1e-6, 1e-6)

    if (p$nloccov > 0) dat[,icov] = Ycov[W[["U"]], icov_local] # no need for other dim checks as this is user provided
    if (exists("TIME", p$variables)) dat[, itime_cov] = as.matrix(stmv_timecovars( vars=ti_cov, ti=Ytime[W[["U"]], ] ) )
    dat = as.data.frame(dat)
    names(dat) = dat_names


    # remember that these are crude mean/discretized estimates
    nu = phi = varSpatial = varObs = NULL
    if ( !is.na(W[["ores"]])) {
      if (exists("nu", W[["ores"]])) if (is.finite(W[["ores"]][["nu"]])) nu = W[["ores"]][["nu"]]
      if (exists("phi", W[["ores"]])) if (is.finite(W[["ores"]][["phi"]])) if (W[["ores"]][["phi"]] > (p$pres/2)) phi = W[["ores"]][["phi"]]
      if (exists("varSpatial", W[["ores"]])) if (is.finite(W[["ores"]][["varSpatial"]])) varSpatial = W[["ores"]][["varSpatial"]]
      if (exists("varObs", W[["ores"]]))  if (is.finite(W[["ores"]][["varObs"]])) if (W[["ores"]][["varObs"]] > p$eps) varObs = W[["ores"]][["varObs"]]
    }

    if (is.null(nu)) nu = p$stmv_lowpass_nu
    if (is.null(nu)) nu = 0.5
    if (!is.finite(nu)) nu = 0.5

    if (is.null(phi)) phi = W[["stmv_distance_cur"]]/sqrt(8*nu) # crude estimate of phi based upon current scaling  distance approximates the range at 90% autocorrelation(e.g., see Lindgren et al. 2011)
    if (!is.finite(phi)) phi = W[["stmv_distance_cur"]]/sqrt(8*nu)

    if (is.null(varSpatial)) varSpatial =0.5 * var( dat[, p$variables$Y], na.rm=TRUE)
    if (is.null(varObs)) varObs = varSpatial
    if ( !is.na(W[["ores"]]) ) W[["ores"]][["vgm"]] = NULL

    if (debugging) {
      print( "starting interpolation" )
      # check data and statistical locations
      plot( Sloc[,], pch=20, cex=0.5, col="gray")
      points( Yloc[,], pch=20, cex=0.2, col="green")
      points( Yloc[W[["U"]],], pch=20, cex=1, col="yellow" )
      points( Sloc[Si,2] ~ Sloc[Si,1], pch=20, cex=5, col="blue" )
    }

    # construct data (including covariates) for prediction locations (pa)
    pa = try( stmv_predictionarea( p=p, sloc=Sloc[Si,], windowsize.half=p$windowsize.half ) )
    if (is.null(pa)) {
      Sflag[Si] = E[["prediction_area"]]
      message( Si )
      message("Error with issue with prediction grid ... null .. this should not happen")
      pa = W = NULL
      next()
    }
    if ( inherits(pa, "try-error") ) {
      pa = W = NULL
      Sflag[Si] = E[["prediction_area"]]
      next()
    }

    if (debugging) {
      # check that position indices are working properly
      Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
      Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )
      plot( Yloc[W[["U"]],2] ~ Yloc[W[["U"]],1], col="red", pch=".",
        ylim=range(c(Yloc[W[["U"]],2], Sloc[Si,2], Ploc[pa$i,2]) ),
        xlim=range(c(Yloc[W[["U"]],1], Sloc[Si,1], Ploc[pa$i,1]) ) ) # all data
      points( Yloc[W[["U"]],2] ~ Yloc[W[["U"]],1], col="green" )  # with covars and no other data issues
      points( Sloc[Si,2] ~ Sloc[Si,1], col="blue" ) # statistical locations
      # statistical output locations
      grids= aegis::spatial_grid(p, DS="planar.coords" )
      points( grids$plat[round( (Sloc[Si,2]-p$origin[2])/p$pres) + 1]
            ~ grids$plon[round( (Sloc[Si,1]-p$origin[1])/p$pres) + 1] , col="purple", pch=25, cex=5 )
      points( grids$plat[pa$iplat] ~ grids$plon[ pa$iplon] , col="cyan", pch=20, cex=0.01 ) # check on Proc iplat indexing
      points( Ploc[pa$i,2] ~ Ploc[ pa$i, 1] , col="black", pch=20, cex=0.7 ) # check on pa$i indexing -- prediction locations
    }


    if ( exists("force_complete_solution", p)) {
      if (p$force_complete_solution) {
        # augment data with prior estimates and predictions
        dist_fc = floor(W[["stmv_distance_cur"]]/p$pres)
        pa_fc = try( stmv_predictionarea( p=p, sloc=Sloc[Si,], windowsize.half=dist_fc ) )
        if (is.null(pa_fc)) {
          Sflag[Si] = E[["prediction_area"]]
          message( Si )
          message("Error with issue with prediction grid ... null .. this should not happen")
          W = pa_fc = dist_fc = NULL
          next()
        }
        if ( inherits(pa_fc, "try-error") ) {
          W = pa_fc = dist_fc = NULL
          Sflag[Si] = E[["prediction_area"]]
          next()
        }
        if (nrow(pa_fc) > 1) {
          augmented_data = P[pa_fc$i]
          good = which( is.finite(augmented_data))
          ngood = length(good)
          if ( ngood > 1) {
            pa_fc = pa_fc[good,]
            pa_fc[, p$variable$Y] = augmented_data[good] # copy
            if ( (ngood + W[["ndata"]] ) > p$n.max ) {
              nmore = p$n.max - W[["ndata"]]
              keep = .Internal( sample( ngood, nmore, replace=FALSE, prob=NULL) ) # thin
              pa_fc = pa_fc[keep,]
            }
            pa_fc$weights = 1
            dat = rbind(dat, pa_fc[, dat_names])
          }
        }
        pa_fc = keep = augmented_data = nmore = ngood = NULL
      }
    }


    # model and prediction .. outputs are in scale of the link (and not response)
    # the following permits user-defined models (might want to use compiler::cmpfun )
    res =NULL
    res = try(
      local_fn(
        p=p,
        dat=dat,
        pa=pa,
        nu=nu,
        phi=phi,
        varObs=varObs,
        varSpatial=varSpatial,
        sloc=Sloc[Si,],
        distance=W[["stmv_distance_cur"]]
      )
    )

  #
  # if (interp.method == "multilevel.b.splines") {
  #   library(MBA)
  #   out = mba.surf(data, no.X=nr, no.Y=nc, extend=TRUE)
  #   if (0) {
  #     image(out, xaxs = "r", yaxs = "r", main="Observed response")
  #     locs= cbind(data$x, data$y)
  #     points(locs)
  #     contour(out, add=T)
  #   }
  #   return(out$xyz.est)
  # }


    if (debugging) {
      print (str(W) )
      print( str(res) )
      if (0) {
        lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
        lattice::levelplot( mean ~ plon + plat, data=res$predictions, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
        for( i in sort(unique(res$predictions[,p$variables$TIME])))  print(lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==i,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
      }
    }

    if ( is.null(res)) {
      Sflag[Si] = E[["local_model_error"]]   # modelling / prediction did not complete properly
      dat = pa = NULL
      next()
    }

    if ( inherits(res, "try-error") ) {
      Sflag[Si] =  E[["local_model_error"]]   # modelling / prediction did not complete properly
      dat = pa = res = NULL
      next()
    }

    if (!exists("predictions", res)) {
      Sflag[Si] =  E[["prediction_error"]]   # modelling / prediction did not complete properly
      dat = pa = res = NULL
      next()
    }

    if (!exists("mean", res$predictions)) {
      Sflag[Si] =  E[["prediction_error"]]   # modelling / prediction did not complete properly
      dat = pa = res = NULL
      next()
    }

    if (length(which( is.finite(res$predictions$mean ))) < 5) {
      Sflag[Si] =  E[["prediction_error"]]   # modelling / prediction did not complete properly
      dat = pa = res = NULL
      next()  # looks to be a faulty solution
    }

    if (exists( "stmv_quantile_bounds", p)) {
      tq = quantile( Y[W[["U"]]], probs=p$stmv_quantile_bounds, na.rm=TRUE  )
      toolow  = which( res$predictions$mean < tq[1] )
      toohigh = which( res$predictions$mean > tq[2] )
      if (length( toolow) > 0)  res$predictions$mean[ toolow] = tq[1]
      if (length( toohigh) > 0) res$predictions$mean[ toohigh] = tq[2]
      toolow = toohigh = tq= NULL
    }


    # extract stats and compute a few more things
    sf = try( stmv_statistics_update( p=p, res=res, W=W, Si=Si ) )
    if ( is.null(sf) ) {
      Sflag[Si] = E[["statistics_update_error"]]
      res = pa = sf = NULL
      next()
    }
    if ( inherits(sf, "try-error") ) {
      Sflag[Si] = E[["statistics_update_error"]]
      res = pa = sf = NULL
      next()
    }
    if ( sf=="error" ) {
      Sflag[Si] =  E[["statistics_update_error"]]
      res = pa = sf = NULL
      next()
    }

    res$stmv_stats = NULL # reduce memory usage
    sf = try( stmv_predictions_update(p=p, preds=res$predictions ) )
    if ( is.null(sf) ) {
      Sflag[Si] = E[["prediction_update_error"]]
      res = pa = sf = NULL
      next()
    }
    if ( inherits(sf, "try-error") ) {
      Sflag[Si] = E[["prediction_update_error"]]
      res = pa = sf = NULL
      next()
    }
    if ( sf=="error" ) {
      Sflag[Si] = E[["prediction_update_error"]]
      res = pa = sf = NULL
      next()
    }

    res = pa = NULL


    # ----------------------
    # do last. it is an indicator of completion of all tasks
    # restarts would be broken otherwise
    Sflag[Si] = E[["complete"]]  # mark as complete without issues

    if (0) {
      # only for debugging data bigmemory structures
      if ( iip %in% savepoints ) {
        sP = P[]; save( sP, file=p$saved_state_fn$P, compress=TRUE ); sP=NULL
        sPn = Pn[]; save( sPn, file=p$saved_state_fn$Pn, compress=TRUE ); sPn=NULL
        sPsd = Psd[]; save( sPsd, file=p$saved_state_fn$Psd, compress=TRUE ); sPsd=NULL
        sS = S[]; save( sS, file=p$saved_state_fn$stats, compress=TRUE ); sS=NULL
        sSflag = Sflag[]; save( sSflag, file=p$saved_state_fn$sflag, compress=TRUE ); sSflag=NULL
        currentstatus = stmv_logfile(p=p)
      }
    }

  }  # end for loop

  return(NULL)

}
