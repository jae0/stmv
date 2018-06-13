

stmv_interpolate = function( ip=NULL, p, debugging=FALSE, stime=Sys.time(), ... ) {
  #\\ core function to interpolate (model and predict) in parallel

  if (0) {
    # for debugging  runs ..
    currentstatus = stmv_db( p=p, DS="statistics.status" )
    p = parallel_run( p=p, runindex=list( locs=sample( currentstatus$todo )) )
    ip = 1:p$nruns
    debugging=TRUE
    stime=Sys.time()
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

  Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
  Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
  Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )

  Y = stmv_attach( p$storage.backend, p$ptr$Y )

  P = stmv_attach( p$storage.backend, p$ptr$P )
  Pn = stmv_attach( p$storage.backend, p$ptr$Pn )
  Psd = stmv_attach( p$storage.backend, p$ptr$Psd )


  if (p$nloccov > 0) {
    Ycov = stmv_attach( p$storage.backend, p$ptr$Ycov )
  }

  if ( exists("TIME", p$variables) ) {
    Ytime = stmv_attach( p$storage.backend, p$ptr$Ytime )
  }

  Yi = stmv_attach( p$storage.backend, p$ptr$Yi )

  # misc intermediate calcs to be done outside of parallel loops
  upsampling = c(1,0, 1.25, 1.5, 1.75, 2.0, 2.5) * p$stmv_distance_scale
  upsampling = upsampling[ which(upsampling <= p$stmv_distance_max )]

  # pre-calculate indices and dim for data to use inside the loop
  dat_names = unique( c(  p$variables$Y, p$variable$LOCS, p$variables$local_all,  "weights") )
  dat_nc = length( dat_names )
  iY = which(dat_names== p$variables$Y)
  ilocs = which( dat_names %in% p$variable$LOCS )
  # iwei = which( dat_names %in% "weights" )

  if (p$nloccov > 0) icov = which( dat_names %in% p$variables$local_cov )
  if (exists("TIME", p$variables)) {
    ti_cov = setdiff(p$variables$local_all, c(p$variables$Y, p$variables$LOCS, p$variables$local_cov ) )
    itime_cov = which(dat_names %in% ti_cov)
  }

  local_fn = switch( p$stmv_local_modelengine,
    bayesx = stmv__bayesx,
    gaussianprocess2Dt = stmv__gaussianprocess2Dt,
    inla = stmv__inla,
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

    if ( iip %in% logpoints )  currentstatus = stmv_logfile(p=p, stime=stime)

    if (0) {
      # only for debugging
      if ( iip %in% savepoints ) {
        sP = P[]; save( sP, file=p$saved_state_fn$P, compress=TRUE ); sP=NULL
        sPn = Pn[]; save( sPn, file=p$saved_state_fn$Pn, compress=TRUE ); sPn=NULL
        sPsd = Psd[]; save( sPsd, file=p$saved_state_fn$Psd, compress=TRUE ); sPsd=NULL
        sS = S[]; save( sS, file=p$saved_state_fn$stats, compress=TRUE ); sS=NULL
        sSflag = Sflag[]; save( sSflag, file=p$saved_state_fn$sflag, compress=TRUE ); sSflag=NULL
        currentstatus = stmv_logfile(p=p, stime=stime)
      }
    }

    Si = p$runs[ iip, "locs" ]

    if ( Sflag[Si] != 0L ) next()  # previously attempted .. skip
      # 0=to do
      # 1=complete
      # 2=oustide bounds(if any) -
      # 3=shallow(if z is a covariate) -
      # 4=predictionarea not ok,
      # 5=skipped due to insufficient data,
      # 6=skipped .. fast variogram did not work
      # 7=variogram estimated range not ok
      # 8=problem with prediction and/or modelling
      # 9=attempting ... if encountered then it was some general problem  or was interrrupted

    # find data nearest S[Si,] and with sufficient data
    dlon = abs( Sloc[Si,1] - Yloc[Yi[],1] )
    dlat = abs( Sloc[Si,2] - Yloc[Yi[],2] )

    ndata = 0
    for ( stmv_distance_cur in upsampling )  {
      U = which( {dlon  <= stmv_distance_cur} & {dlat <= stmv_distance_cur} )  # faster to take a block
      ndata = length(U)
      if ( ndata >= p$n.min ) break()
    }
    if (ndata < p$n.min) {
      Sflag[Si] = 5L   # skipped .. not enough data
      next()
    }
    
    print( paste(iip, Si, ndata ) )

    if (0) {
      plot( Sloc[,], pch=20, cex=0.5, col="gray")
      points( Yloc[,], pch=20, cex=0.2, col="green")
      points( Yloc[U,], pch=20, cex=1, col="yellow" )
      points( Sloc[Si,2] ~ Sloc[Si,1], pch=20, cex=5, col="blue" )
    }

    # crude (mean) variogram across all time slices
    o = NULL
    o = try( stmv_variogram( xy=Yloc[U,], z=Y[U], methods=p$stmv_variogram_method,
      distance_cutoff=stmv_distance_cur, nbreaks=13 ) )
    if ( is.null(o)) Sflag[Si] = 6L   # fast variogram did not work
    if ( inherits(o, "try-error")) Sflag[Si] = 6L   # fast variogram did not work
    if (exists("stmv_rangecheck", p)) if (p$stmv_rangecheck=="paranoid") if ( Sflag[Si] == 6L ) next()

    ores = NULL
    if ( exists(p$stmv_variogram_method, o)) {
      ores = o[[p$stmv_variogram_method]] # store current best estimate of variogram characteristics
      if ( !exists("range_ok", ores) ){
        Sflag[Si] = 7L
      } else {
        if ( ores[["range_ok"]] ) {
          stmv_distance_cur = ores[["range"]]
          vario_U  = which( {dlon  <= ores[["range"]] } & {dlat <= ores[["range"]]} )
          vario_ndata =length(vario_U)
          if (vario_ndata < p$n.min) {
            # insufficient data at estimated range
            # ..could  stop analysis but this is not necessary .. range identifies the distance
            # at which AC is no differnt from background noise and so addition of more distant data
            # does not alter interpretation.
            # NOTE: this range is a crude estimate that averages across years (if any) ...
            if (exists("stmv_rangecheck", p)) {
              if (p$stmv_rangecheck=="paranoid") {
                Sflag[Si] = 5L
                next()
              }
            }
          }
          U  = vario_U
          ndata = vario_ndata
          vario_U = vario_ndata = NULL
        }
      }
    }

    if (ndata > p$n.max) {
      # U = U[ .Internal( sample( vario_ndata, p$n.max, replace=FALSE, prob=NULL)) ] # simple random
      if (exists("TIME", p$variables)) {
        iU = stmv_discretize_coordinates( coo=Yloc[U,], times=Ytime[YiU, ], ntarget=p$n.max,
          minresolution=p$downsampling_multiplier*c(p$pres, p$pres, p$tres), method="thin" )
      } else {
        iU = stmv_discretize_coordinates( coo=Yloc[U,], ntarget=p$n.max,
          minresolution=p$downsampling_multiplier*c(p$pres, p$pres ), method="thin" )
      }
      U = U[iU]
      ndata = length(U)
    }

    iU=dlon=dlat=o=NULL

    YiU = Yi[U] # YiU and p$stmv_distance_prediction determine the data entering into local model construction
    pa = stmv_predictionarea( p=p, sloc=Sloc[Si,], windowsize.half=p$windowsize.half )
    if (is.null(pa)) {
      Sflag[Si] = 4L
      message( Si )
      message("Error with issue with prediction grid ... null .. this should not happen")
      next()
    }

    if (0) {
      # check that position indices are working properly
      Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
      Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )
      plot( Yloc[U,2]~ Yloc[U,1], col="red", pch=".",
        ylim=range(c(Yloc[U,2], Sloc[Si,2], Ploc[pa$i,2]) ),
        xlim=range(c(Yloc[U,1], Sloc[Si,1], Ploc[pa$i,1]) ) ) # all data
      points( Yloc[YiU,2] ~ Yloc[YiU,1], col="green" )  # with covars and no other data issues
      points( Sloc[Si,2] ~ Sloc[Si,1], col="blue" ) # statistical locations
      # statistical output locations
      grids= aegis::spatial_grid(p, DS="planar.coords" )
      points( grids$plat[round( (Sloc[Si,2]-p$origin[2])/p$pres) + 1]
            ~ grids$plon[round( (Sloc[Si,1]-p$origin[1])/p$pres) + 1] , col="purple", pch=25, cex=5 )

      points( grids$plat[pa$iplat] ~ grids$plon[ pa$iplon] , col="cyan", pch=20, cex=0.01 ) # check on Proc iplat indexing
      points( Ploc[pa$i,2] ~ Ploc[ pa$i, 1] , col="black", pch=20, cex=0.7 ) # check on pa$i indexing -- prediction locations
    }


    # prep dependent data
    # reconstruct data for modelling (dat) and data for prediction purposes (pa)
    dat = matrix( 1, nrow=ndata, ncol=dat_nc )
    dat[,iY] = Y[YiU] # these are residuals if there is a global model
    dat[,ilocs] = Yloc[YiU,] + stmv_distance_cur * runif(2*ndata, -1e-6, 1e-6) # add a small error term to prevent some errors when duplicate locations exist; stmv_distance_cur offsets to positive values


    if (p$nloccov > 0) dat[,icov] = Ycov[YiU,] # no need for other dim checks as this is user provided
    if (exists("TIME", p$variables)) dat[, itime_cov] = as.matrix(stmv_timecovars( vars=ti_cov, ti=Ytime[YiU, ] ) )
    dat = as.data.frame(dat)
    names(dat) = dat_names

    # remember that these are crude mean/discretized estimates
    nu = phi = varSpatial = varObs = NULL
    if (!is.null(ores)) {
      if (exists("nu", ores)) if (is.finite(ores$nu)) nu = ores$nu
      if (exists("phi", ores)) if (is.finite(ores$phi)) if (ores$phi > (p$pres/2)) phi = ores$phi
      if (exists("varSpatial", ores)) if (is.finite(ores$varSpatial)) varSpatial = ores$varSpatial
      if (exists("varObs", ores))  if (is.finite(ores$varObs)) if (ores$varObs > p$eps) varObs = ores$varObs
    }

    if (is.null(nu)) nu = p$stmv_lowpass_nu
    if (is.null(nu)) nu = 0.5
    if (!is.finite(nu)) nu = 0.5

    if (is.null(phi)) phi = stmv_distance_cur/sqrt(8*nu) # crude estimate of phi based upon current scaling  distance approximates the range at 90% autocorrelation(e.g., see Lindgren et al. 2011)
    if (!is.finite(phi)) phi = stmv_distance_cur/sqrt(8*nu)

    if (is.null(varSpatial)) varSpatial =0.5 * var( dat[, p$variables$Y], na.rm=TRUE)

    if (is.null(varObs)) varObs = varSpatial


    # print( "starting interpolation" )


    ores$vgm = NULL # can be large

    # model and prediction .. outputs are in scale of the link (and not response)
    # the following permits user-defined models (might want to use compiler::cmpfun )

    res =NULL
    res = try( local_fn( p=p, dat=dat, pa=pa, nu=nu, phi=phi, varObs=varObs, varSpatial=varSpatial, sloc=Sloc[Si,], distance=stmv_distance_cur ) )

    if (debugging) print( str(res))

    if (0) {
      lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
      lattice::levelplot( mean ~ plon + plat, data=res$predictions, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
      for( i in sort(unique(res$predictions[,p$variables$TIME])))  print(lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==i,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
    }

    if ( is.null(res)) {
      Sflag[Si] = 8L   # modelling / prediction did not complete properly
      dat = pa = res = NULL
      next()
    }

    if ( inherits(res, "try-error") ) {
      Sflag[Si] = 8L   # modelling / prediction did not complete properly
      dat = pa = res = NULL
      next()
    }

    if (!exists("predictions", res)) {
      Sflag[Si] = 8L   # modelling / prediction did not complete properly
      dat = pa = res = NULL
      next()
    }

    if (!exists("mean", res$predictions)) {
      Sflag[Si] = 8L   # modelling / prediction did not complete properly
      dat = pa = res = NULL
      next()
    }

    if (length(which( is.finite(res$predictions$mean ))) < 5) {
      Sflag[Si] = 8L   # modelling / prediction did not complete properly
      dat = pa = res = NULL
      next()  # looks to be a faulty solution
    }

    if (exists( "stmv_quantile_bounds", p)) {
      tq = quantile( Y[YiU], probs=p$stmv_quantile_bounds, na.rm=TRUE  )
      toolow  = which( res$predictions$mean < tq[1] )
      toohigh = which( res$predictions$mean > tq[2] )
      if (length( toolow) > 0)  res$predictions$mean[ toolow] = tq[1]
      if (length( toohigh) > 0) res$predictions$mean[ toohigh] = tq[2]
      toolow = toohigh = tq= NULL
    }

    # stats collator
    if (!exists("stmv_stats",  res) ) res$stmv_stats = list()

    if (!exists("sdSpatial", res$stmv_stats)) {
      # some methods can generate spatial stats simultaneously ..
      # it is faster to keep them all together instead of repeating here
      # field and RandomFields gaussian processes seem most promising ...
      # default to fields for speed:
      res$stmv_stats["sdSpatial"] = NA
      res$stmv_stats["sdObs"] = NA
      res$stmv_stats["range"] = NA
      res$stmv_stats["phi"] = NA
      res$stmv_stats["nu"] = NA
      if ( !is.null(ores)) {
        if ( exists("varSpatial", ores) ) res$stmv_stats["sdSpatial"] = sqrt( ores[["varSpatial"]] )
        if ( exists("varObs", ores) ) res$stmv_stats["sdObs"] = sqrt(ores[["varObs"]])
        if ( exists("range", ores) ) res$stmv_stats["range"] = ores[["range"]]
        if ( exists("phi", ores) ) res$stmv_stats["phi"] = ores[["phi"]]
        if ( exists("nu", ores) ) res$stmv_stats["nu"] = ores[["nu"]]
      }
    }

    if ( exists("TIME", p$variables) ){
      # annual ts, seasonally centered and spatially
      # pa_i = which( Sloc[Si,1]==Ploc[,1] & Sloc[Si,2]==Ploc[,2] )
      pac_i = which( res$predictions$plon==Sloc[Si,1] & res$predictions$plat==Sloc[Si,2] )
      # plot( mean~tiyr, res$predictions[pac_i,])
      # plot( mean~tiyr, res$predictions, pch="." )
      res$stmv_stats["ar_timerange"] = NA
      res$stmv_stats["ar_1"] = NA

      if (length(pac_i) > 5) {
        pac = res$predictions[ pac_i, ]
        pac$dyr = pac[, p$variables$TIME] - trunc(pac[, p$variables$TIME] )
        piid = which( zapsmall( pac$dyr - p$dyear_centre) == 0 )
        pac = pac[ piid, c(p$variables$TIME, "mean")]
        pac = pac[ order(pac[,p$variables$TIME]),]
        if (length(piid) > 5 ) {
          ts.stat = NULL
          ts.stat = try( stmv_timeseries( pac$mean, method="fft" ) )
          if (!is.null(ts.stat) && !inherits(ts.stat, "try-error") ) {
            res$stmv_stats["ar_timerange"] = ts.stat$quantilePeriod
            if (all( is.finite(pac$mean))) {
              afin = which (is.finite(pac$mean) )
              if (length(afin) > 5 && var( pac$mean, na.rm=TRUE) > p$eps ) {
                ar1 = NULL
                ar1 = try( ar( pac$mean, order.max=1 ) )
                if (!inherits(ar1, "try-error")) {
                  if ( length(ar1$ar) == 1 ) {
                    res$stmv_stats["ar_1"] = ar1$ar
                  }
                }
              }
            }
            if ( !is.finite(res$stmv_stats[["ar_1"]]) ) {
              ar1 = try( cor( pac$mean[1:(length(piid) - 1)], pac$mean[2:(length(piid))], use="pairwise.complete.obs" ) )
              if (!inherits(ar1, "try-error")) res$stmv_stats["ar_1"] = ar1
            }
          }

          ### Do the logistic model here ! -- if not already done ..
          if (!exists("ts_K", res$stmv_stats)) {
            # model as a logistic with ts_r, ts_K, etc .. as stats outputs

          }

        }
        pac=piid=NULL
      }
      pac_i=NULL
    }

    # update SD estimates of predictions with those from other locations via the
    # incremental  method ("online algorithm") of mean estimation after Knuth ;
    # see https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
    # update means: inverse-variance weighting
    # see https://en.wikipedia.org/wiki/Inverse-variance_weighting

    npred = nrow(res$predictions)

    if ( ! exists("TIME", p$variables) ) {

      u = which( is.finite( P[res$predictions$i] ) )  # these have data already .. update
      if ( length( u ) > 1 ) {
        ui = res$predictions$i[u]  # locations of P to modify
        Pn[ui] = Pn[ui] + 1 # update counts
        stdev_update =  Psd[ui] + ( res$predictions$sd[u] -  Psd[ui] ) / Pn[ui]
        means_update = ( P[ui] / Psd[ui]^2 + res$predictions$mean[u] / res$predictions$sd[u]^2 ) / ( Psd[ui]^(-2) + res$predictions$sd[u]^(-2) )
        mm = which(is.finite( means_update + stdev_update ))
        if( length(mm)> 0) {
          iumm = ui[mm]
          Psd[iumm] = stdev_update[mm]
          P  [iumm] = means_update[mm]
        }
        stdev_update = NULL
        means_update = NULL
        ui=mm=iumm=NULL
      }

      # first time # no data yet
      v = setdiff(1:npred, u)
      if ( length(v) > 0 ) {
        vi = res$predictions$i[v]
        Pn [vi] = 1
        P  [vi] = res$predictions$mean[v]
        Psd[vi] = res$predictions$sd[v]
      }
      vi = NULL
    }


    if ( exists("TIME", p$variables) ) {
      u = which( is.finite( P[res$predictions$i,1] ) )  # these have data already .. update
      u_n = length( u )
      if ( u_n > 1 ) {  # ignore if only one point .. mostly because it can cause issues with matrix form ..
        # locations of P to modify
        ui = sort(unique(res$predictions$i[u]))
        nc = ncol(P)
        if (p$storage.backend == "ff" ) {
          add.ff(Pn, 1, ui, 1:nc ) # same as Pn[ui,] = Pn[ui]+1 but 2X faster
        } else {
          Pn[ui,] = Pn[ui,] + 1
        }
        stdev_update =  Psd[ui,] + ( res$predictions$sd[u] -  Psd[ui,] ) / Pn[ui,]
        means_update = ( P[ui,] / Psd[ui,]^2 + res$predictions$mean[u] / res$predictions$sd[u]^2 ) /
          ( Psd[ui,]^(-2) + res$predictions$sd[u]^(-2) )

        updates = means_update + stdev_update
        if (!is.matrix(updates)) {
          Sflag[Si] = 8L
          u = u_n = ui = nc =stdev_update = means_update = NULL
          message( Si )
          message( "update of predictions were problematic ... this should not happen, proabbaly due to NA's" )
          next()
        }

        mm = which( is.finite( rowSums(updates)))  # created when preds go outside quantile bounds .. this removes all data from a given location rather than the space-time .. severe but likely due to a poor prediction and so remove all (it is also faster this way as few manipulations)
        if( length(mm)> 0) {
          iumm = ui[mm]
          Psd[iumm,] = stdev_update[mm,]
          P  [iumm,] = means_update[mm,]
          iumm = NULL
        }
        stdev_update = means_update = updates = ui = mm=NULL
      }

      # do this as a second pass in case NA's were introduced by the update .. unlikely , but just in case
      v = which( !is.finite( P[res$predictions$i,1] ) )  # these have data already .. update
      nv = length(v)          # no data yet
      if ( nv > 0 ) {
        vi = sort(unique(res$predictions$i[v]))
        Pn [vi,] = 1
        P  [vi,] = res$predictions$mean[v]
        Psd[vi,] = res$predictions$sd[v]
        vi = NULL
      }
    }

    # save stats
    for ( k in 1: length(p$statsvars) ) {
      if (exists( p$statsvars[k], res$stmv_stats )) {
        S[Si,k] = res$stmv_stats[[ p$statsvars[k] ]]
      }
    }

    if (debugging) {
        v = res$predictions
        if ( exists("TIME", p$variables) ){
          v = v[which( v[,p$variables$TIME]==2000.55),]
        }
        require(lattice)
        print(
          levelplot( mean ~ plon+plat, v, aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=TRUE) )
      )
    }

    if (0) {
      if ("time slice at 2012.05") {
        lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
      }

      if ("all TIME time slices from latest predictions") {
        for( i in sort(unique(res$predictions[,p$variables$TIME])))  {
          print(lattice::levelplot( mean ~ plon + plat, data=res$predictions[res$predictions[,p$variables$TIME]==i,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
        }
      }

      if ("all nt time slices in stored predictions P") {
        for (i in 1:p$nt) {
          print( lattice::levelplot( P[pa$i,i] ~ Ploc[pa$i,1] + Ploc[ pa$i, 2], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
        }
      }

      if ("no time slices in P") {
          print( lattice::levelplot( P[pa$i] ~ Ploc[pa$i,1] + Ploc[ pa$i, 2], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
      }

    }


    # ----------------------
    # do last. it is an indicator of completion of all tasks
    # restarts would be broken otherwise
    Sflag[Si] = 1L  # mark as complete without issues
    res = pa = NULL

  }  # end for loop

  return(NULL)

}
