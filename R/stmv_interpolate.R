

stmv_interpolate = function( ip=NULL, p, debugging=FALSE, runoption="default", ... ) {
  #\\ core function to interpolate (model and predict) in parallel

  if (0) {
    # for debugging  runs ..
    currentstatus = stmv_statistics_status( p=p )
    p = parallel_run( p=p, runindex=list( locs=sample( currentstatus$todo )) )
    ip = 1:p$nruns
    debugging=TRUE
    runoption="default"
  }

  # ---------------------
  # deal with additional passed parameters
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
  Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )

  Y = stmv_attach( p$storage.backend, p$ptr$Y )
  Yi = stmv_attach( p$storage.backend, p$ptr$Yi )  # initial indices of good data

  if (p$nloccov > 0) Ycov = stmv_attach( p$storage.backend, p$ptr$Ycov )
  if ( exists("TIME", p$variables) ) Ytime = stmv_attach( p$storage.backend, p$ptr$Ytime )

  # misc intermediate calcs to be done outside of parallel loops

  # pre-calculate indices and dim for data to use inside the loop
  dat_names = unique( c(  p$variables$Y, p$variable$LOCS, p$variables$local_all,  "weights") )  # excludes p$variables$TIME
  if (p$stmv_local_modelengine %in% c("fft", "tps") ) {
    if ( exists("TIME", p$variables)) {
      dat_names = c(dat_names, p$variables$TIME)
    }
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

  local_fn = ifelse (p$stmv_local_modelengine=="userdefined", p$stmv_local_modelengine_userdefined, stmv_interpolation_function( p$stmv_local_modelengine ) )

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

  # global estimates
  vg_global = list(
    range = max( min( median( S[, match( "range", p$statsvars )], na.rm=TRUE ), max(p$stmv_distance_scale )),  min(p$pres, p$stmv_distance_scale )),
    nu = median( S[, match("nu", p$statsvars)], na.rm=TRUE ),
    varObs = median( S[, match("varObs", p$statsvars)], na.rm=TRUE ),
    varSpatial = median( S[, match("varSpatial", p$statsvars)], na.rm=TRUE ),
    phi = median( S[, match("phi", p$statsvars)], na.rm=TRUE )
  )

  if( vg_global$nu < 0.25 | vg_global$nu > 4) vg_global$nu = 0.5

  range_index = match( "range", p$statsvars )
  ndata_index = match( "ndata", p$statsvars )
  rsquared_index = match("rsquared", p$statsvars )


# main loop over each output location in S (stats output locations)
  for ( iip in ip ) {

    if ( iip %in% logpoints )  currentstatus = stmv_logfile(p=p, flag= paste("Interpolation", runoption) )
    Si = p$runs[ iip, "locs" ]

    # print( paste("index =", iip, ";  Si = ", Si ) )
    if (debugging) print( paste("index =", iip, ";  Si = ", Si ) )
    if ( Sflag[Si] == E[["complete"]] ) next()

    vg = stmv_scale_filter( p=p, Si=Si )
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
    if (useglobal) vg = vg_global

    if (runoption=="boostdata") localrange = matern_phi2distance( phi=vg$phi, nu=vg$nu, cor=p$stmv_range_correlation_boostdata )

    U = stmv_select_data( p=p, Si=Si, localrange=localrange )

    if ( Sflag[Si] != E[["todo"]] ) {
      if (exists("stmv_rangecheck", p)) {
        if (p$stmv_rangecheck=="paranoid") {
          if ( Sflag[Si] %in% c( E[["variogram_range_limit"]], E[["variogram_failure"]]) ) {
            U = NULL
            if (debugging) message("Error: stmv_rangecheck paranoid")
            if (runoption == "default" ) next()
          }
        }
      }
    }

    if(is.null( U )) next()

    # last check
    ndata = length(U)
    if (ndata < p$stmv_nmin) next()

    # if here then there is something to do
    # NOTE:: U are the indices of locally useful data
    # p$stmv_distance_prediction determines the data entering into local model construction

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
    if (debugging) {
      dev.new()
      # check data and statistical locations
      plot( Sloc[,], pch=20, cex=0.5, col="gray")
      points( Yloc[,], pch=20, cex=0.2, col="green")
      points( Yloc[U,], pch=20, cex=1, col="yellow" )
      points( Sloc[Si,2] ~ Sloc[Si,1], pch=20, cex=5, col="blue" )
    }

    # construct data (including covariates) for prediction locations (pa)
    pa = try( stmv_predictionarea( p=p, sloc=Sloc[Si,], windowsize.half=p$windowsize.half ) )
    if (is.null(pa)) {
      Sflag[Si] = E[["prediction_area"]]
      if (debugging) message( Si )
      if (debugging) message("Error: prediction grid ... null .. this should not happen")
      pa = NULL
      next()
    }
    if ( inherits(pa, "try-error") ) {
      pa = NULL
      Sflag[Si] = E[["prediction_area"]]
      if (debugging) message("Error: prediction grid ... try-error .. this should not happen")
      next()
    }

    if (debugging) {

      dev.new()
      # check that position indices are working properly
      Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
      Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
      Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )
      plot( Yloc[U,2] ~ Yloc[U,1], col="red", pch=".",
        ylim=range(c(Yloc[U,2], Sloc[Si,2], Ploc[pa$i,2]) ),
        xlim=range(c(Yloc[U,1], Sloc[Si,1], Ploc[pa$i,1]) ) ) # all data
      points( Yloc[U,2] ~ Yloc[U,1], col="green" )  # with covars and no other data issues
      points( Sloc[Si,2] ~ Sloc[Si,1], col="blue" ) # statistical locations
      # statistical output locations
      grids= spatial_grid(p, DS="planar.coords" )
      points( grids$plat[floor( (Sloc[Si,2]-p$origin[2])/p$pres) + 1]
            ~ grids$plon[floor( (Sloc[Si,1]-p$origin[1])/p$pres) + 1] , col="purple", pch=25, cex=5 )
      points( grids$plat[pa$iplat] ~ grids$plon[ pa$iplon] , col="cyan", pch=20, cex=0.01 ) # check on Proc iplat indexing
      points( Ploc[pa$i,2] ~ Ploc[ pa$i, 1] , col="black", pch=20, cex=0.7 ) # check on pa$i indexing -- prediction locations
    }



    # model and prediction .. outputs are in scale of the link (and not response)
    # the following permits user-defined models (might want to use compiler::cmpfun )

    res =NULL
    res = try(
      local_fn (
        p=p,
        dat=dat,
        pa=pa,
        nu=vg$nu,
        phi=vg$phi,
        varObs=vg$varObs,
        varSpatial=vg$varSpatial,
        sloc=Sloc[Si,],
        distance=localrange
      )
    )

    dat =  NULL
    pa  =  NULL

    if ( is.null(res)) {
      Sflag[Si] = E[["local_model_error"]]   # modelling / prediction did not complete properly
      if (debugging) message("Error: local model error")
      next()
    }

    if ( inherits(res, "try-error") ) {
      Sflag[Si] =  E[["local_model_error"]]   # modelling / prediction did not complete properly
      res = NULL
      if (debugging) message("Error: local model error")
      next()
    }

    if (!exists("predictions", res)) {
      Sflag[Si] =  E[["prediction_error"]]   # modelling / prediction did not complete properly
      res = NULL
      if (debugging) message("Error: prediction error")
      next()
    }

    if (!exists("mean", res$predictions)) {
      Sflag[Si] =  E[["prediction_error"]]   # modelling / prediction did not complete properly
      res = NULL
      if (debugging) message("Error: prediction error")
      next()
    }

    if (length(which( is.finite(res$predictions$mean ))) < 5) {
      Sflag[Si] =  E[["prediction_error"]]   # modelling / prediction did not complete properly
      res = NULL
      if (debugging) {
        message("Error: prediction error")
        browser()
      }
      next()  # looks to be a faulty solution
    }


    if (runoption=="default") {
      # update to rsquared and "range" in stats
      S[Si, rsquared_index] = res$stmv_stats$rsquared
      S[Si, ndata_index] = ndata
    }

    res$stmv_stats = NULL # reduce memory usage


    sf = try( stmv_predictions_update(p=p, preds=res$predictions ) )
    res = NULL

    if ( is.null(sf) ) {
      Sflag[Si] = E[["prediction_update_error"]]
      sf = NULL
      if (debugging) message("Error: prediction update error .. is null")
      next()
    }
    if ( inherits(sf, "try-error") ) {
      Sflag[Si] = E[["prediction_update_error"]]
      sf = NULL
      if (debugging) message("Error: prediction update error .. try-error")
      next()
    }
    if ( sf=="error" ) {
      Sflag[Si] = E[["prediction_update_error"]]
      sf = NULL
      if (debugging) message("Error: prediction update error .. general")
      next()
    }


    # ----------------------
    # do last. it is an indicator of completion of all tasks
    # restarts would be broken otherwise
    Sflag[Si] = E[["complete"]]  # mark as complete without issues

  }  # end for loop


  return(NULL)





    if (0) {

      require(MBA)
      require(fields)

      # kriged
      fit = Krig( dat[, c("plon", "plat")], dat$z, Covariance="Matern", theta=vg$phi, smoothness=0.5)
      x11()
      op = predict(fit)
      tst = cbind( dat[, c("plon", "plat")], op )
      mba.int <- mba.surf( tst, 300, 300, extend=TRUE)$xyz.est
      image(mba.int, xaxs="r", yaxs="r")

      # raw data + mba
      x11()
      tst = cbind(  dat[, c("plon", "plat")], dat$z )
      mba.int <- mba.surf( tst, 300, 300, extend=TRUE)$xyz.est
      image(mba.int, xaxs="r", yaxs="r")

      # mba - pa
      x11()
      tst = cbind( pa$plon, pa$plat, pa$mean )
      mba.int <- mba.surf( tst, 300, 300, extend=TRUE)$xyz.est
      image(mba.int, xaxs="r", yaxs="r")

      # default
      x11()
      tst = cbind( res$predictions$plon,  res$predictions$plat,  res$predictions$mean )
      mba.int <- mba.surf( tst, 300, 300, extend=TRUE)$xyz.est
      image(mba.int, xaxs="r", yaxs="r")

      # kernel-based
      tst = as.image( Z=dat$z, x=dat[, c("plon", "plat")], nx=300, ny=300, na.rm=TRUE)
      out = fields::image.smooth( tst, theta=vg$phi/300, xwidth=p$pres, ywidth=p$pres )
      image(out)

      print( str(res) )

      lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
      lattice::levelplot( mean ~ plon + plat, data=res$predictions, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
      for( i in sort(unique(res$predictions[,p$variables$TIME])))  print(lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==i,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )

      dev.new()
      plot(  dat[,iY] ~ dat$yr, col="red"  )
      points( mean~tiyr, res$predictions, pch=20, col="gray", cex=0.5 )

    }


}
