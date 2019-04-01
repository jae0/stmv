

stmv_scale = function( ip=NULL, p, debugging=FALSE, ... ) {
  #\\ core function to interpolate (model variogram) in parallel

  if (0) {
    # for debugging  runs ..
    currentstatus = stmv_db( p=p, DS="statistics.status" )
    p = parallel_run( p=p, runindex=list( locs=sample( currentstatus$todo )) )
    ip = 1:p$nruns
    debugging=TRUE
  }

  # ---------------------
  # deal with additional passed parameters
  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast=TRUE))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable

  if (exists( "libs", p)) suppressMessages( RLibrary( p$libs ) )

  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns

  if (!exists("distance_scale_current", p)) p$distance_scale_current = p$stmv_distance_scale[1]

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
  # unless it is an explicit covariate and not a seasonal component there is no need for it
  # .. prediction grids create these from a time grid on the fly
  dat_nc = length( dat_names )
  iY = which(dat_names== p$variables$Y)
  ilocs = which( dat_names %in% p$variable$LOCS )

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

  distance_to_upsample = p$distance_scale_current * p$stmv_distance_upsampling_fraction

# main loop over each output location in S (stats output locations)
  for ( iip in ip ) { 

    if ( iip %in% logpoints )  currentstatus = stmv_logfile(p=p)
    Si = p$runs[ iip, "locs" ]
    if ( Sflag[Si] != E[["todo"]] ) next()  # previously attempted .. skip
    if (debugging) print( paste("index =", iip, ";  Si = ", Si ) )

    # obtain indices of data locations withing a given spatial range, optimally determined via variogram
    
    # find data nearest Sloc[Si,] and with sufficient data

    dlon = abs( Sloc[Si,1] - Yloc[Yi[],1] )
    dlat = abs( Sloc[Si,2] - Yloc[Yi[],2] )
    ndata = 0
    for ( stmv_distance_cur in distance_to_upsample )  {
      U = which( {dlon  <= stmv_distance_cur} & {dlat <= stmv_distance_cur} )  # faster to take a block
      ndata = length(U)
      if ( ndata >= p$n.min ) break()
    }

    out = NULL
    out = list(
      flag=E[["todo"]], 
      ndata=ndata, 
      stmv_distance_cur=stmv_distance_cur, 
      ores=NA 
    )  # basic threshold ... now tweak it

    if (out$ndata < p$n.min) {
      # retain crude estimate and run with it
      out$flag = E[["insufficient_data"]]
      return( out )   #not enough data
    }

    # NOTE: this range is a crude estimate that averages across years (if any) ...
    o = NULL
    o = try( stmv_variogram( 
      xy=Yloc[Yi[U],], 
      z=Y[Yi[U],], 
      methods=p$stmv_variogram_method, 
      distance_cutoff=stmv_distance_cur, 
      nbreaks=15, 
      range_correlation=p$stmv_range_correlation 
    ) )

    if ( is.null(o)) out$flag = E[["variogram_failure"]]
    if ( inherits(o, "try-error")) out$flag = E[["variogram_failure"]]
    if ( !exists(p$stmv_variogram_method, o)) out$flag =  E[["variogram_failure"]]
    if ( out$flag ==  E[["variogram_failure"]] ) return(out)

    if ( !exists("range_ok", o[[p$stmv_variogram_method]]) ) {
      out$flag =  E[["variogram_range_limit"]]
      return( out )     #
    }

    if ( !o[[p$stmv_variogram_method]][["range_ok"]] ) {
      # retain crude estimate and run with it
      out$flag =  E[["variogram_range_limit"]]
      return(out)
    }

    vario_stmv_distance_cur = o[[p$stmv_variogram_method]][["range"]]
    vario_U = which( {dlon  <= vario_stmv_distance_cur } & {dlat <= vario_stmv_distance_cur} )  # dlon dlat indexed on Yi
    vario_ndata =length(vario_U)


    if (vario_ndata <= p$n.max & vario_ndata >= p$n.min) {
      # all good .. update and return
      out$U  = Yi[vario_U]
      out$ndata = vario_ndata
      out$ores = o[[p$stmv_variogram_method]]
      out$stmv_distance_cur = vario_stmv_distance_cur
      return( out)
    }

    # last try, we are here because (vario_ndata > p$n.max)
  
    W=out

    if ( is.null(W) ) {
      Sflag[Si] = E[["insufficient_data"]]
      W = WA = NULL
      if (debugging) message("Error: insufficient data .. null")
      next()
    }
    if ( inherits(W, "try-error") ) {
      Sflag[Si] = E[["insufficient_data"]]
      W = WA = NULL
      if (debugging) message("Error: insufficient data .. try-error")
      next()
    }

    Sflag[Si] = W[["flag"]]  # update flags
    if ( Sflag[Si] == E[["insufficient_data"]] ) {
      W = WA = NULL
      if (debugging) message("Error: insufficient data")
      next()
    }

    if ( Sflag[Si] != E[["todo"]] ) {
      if (exists("stmv_rangecheck", p)) {
        if (p$stmv_rangecheck=="paranoid") {
          W = WA = NULL
          if (debugging) message("Error: stmv_rangecheck paranoid")
          next()
        }
      }
    }

    # if here then W exists and there is something to do
    # NOTE:: W[["U"]] are the indices of locally useful data
    # p$stmv_distance_prediction determines the data entering into local model construction

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

    if (debugging) {
      # check that position indices are working properly
      Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
      Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )
      plot( Yloc[W[["U"]],2] ~ Yloc[W[["U"]],1], col="red", pch=".",
        ylim=range(c(Yloc[W[["U"]],2], Sloc[Si,2]) ),
        xlim=range(c(Yloc[W[["U"]],1], Sloc[Si,1]) ) ) # all data
      points( Yloc[W[["U"]],2] ~ Yloc[W[["U"]],1], col="green" )  # with covars and no other data issues
      points( Sloc[Si,2] ~ Sloc[Si,1], col="blue" ) # statistical locations
      # statistical output locations
      grids= aegis::spatial_grid(p, DS="planar.coords" )
      points( grids$plat[round( (Sloc[Si,2]-p$origin[2])/p$pres) + 1]
            ~ grids$plon[round( (Sloc[Si,1]-p$origin[1])/p$pres) + 1] , col="purple", pch=25, cex=5 )

    }


    if ( exists("force_complete_solution", p)) {
      if (p$force_complete_solution) {
        # augment data with prior estimates and predictions
        dist_fc = floor(W[["stmv_distance_cur"]]/p$pres)
        pa_fc = try( stmv_predictionarea( p=p, sloc=Sloc[Si,], windowsize.half=dist_fc ) )
        if (is.null(pa_fc)) {
          Sflag[Si] = E[["prediction_area"]]
          if (debugging) message( Si )
          if (debugging) message("Error with issue with prediction grid ... null .. this should not happen")
          W = pa_fc = dist_fc = NULL
          next()
        }
        if ( inherits(pa_fc, "try-error") ) {
          W = pa_fc = dist_fc = NULL
          Sflag[Si] = E[["prediction_area"]]
          if (debugging) message("Error: prediction grid ... try-error .. this should not happen")
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


       nu=nu,
        phi=phi,
        varObs=varObs,
        varSpatial=varSpatial,
        sloc=Sloc[Si,],
        distance=W[["stmv_distance_cur"]]


    if (debugging) {
      print (str(W) )
    }


    # extract stats and compute a few more things
  # compute and extract misc summary statistics from results
  S = stmv_attach( p$storage.backend, p$ptr$S )
  Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )

  sflag = "good"

  out = list()

  if (exists("stmv_stats",  res)) out = res$stmv_stats

  if (!exists("sdSpatial", out)) {
    # some methods can generate spatial stats simultaneously ..
    # it is faster to keep them all together instead of repeating here
    # field and RandomFields gaussian processes seem most promising ...
    # default to fields for speed:
    out["sdSpatial"] = NA
    out["sdObs"] = NA
    out["range"] = NA
    out["phi"] = NA
    out["nu"] = NA
    if (!is.null(W)) {
      if ( !is.na(W[["ores"]])) {
        if ( exists("varSpatial", W[["ores"]]) ) out["sdSpatial"] = sqrt( W[["ores"]][["varSpatial"]] )
        if ( exists("varObs", W[["ores"]]) )     out["sdObs"] = sqrt(W[["ores"]][["varObs"]])
        if ( exists("range", W[["ores"]]) )      out["range"] = W[["ores"]][["range"]]
        if ( exists("phi", W[["ores"]]) )        out["phi"] = W[["ores"]][["phi"]]
        if ( exists("nu", W[["ores"]]) )         out["nu"] = W[["ores"]][["nu"]]
      }
    }

  }

  if ( exists("TIME", p$variables) ){
    # annual ts, seasonally centered and spatially
    # pa_i = which( Sloc[Si,1]==Ploc[,1] & Sloc[Si,2]==Ploc[,2] )
    pac_i = which(
      (abs( res$predictions$plon - Sloc[Si,1] ) < p$pres) &
      (abs( res$predictions$plat - Sloc[Si,2] ) < p$pres)
    )
    # plot( mean~tiyr, res$predictions[pac_i,])
    # plot( mean~tiyr, res$predictions, pch="." )
    out["ar_timerange"] = NA
    out["ar_1"] = NA

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
          out["ar_timerange"] = ts.stat$quantilePeriod
          if (all( is.finite(pac$mean))) {
            afin = which (is.finite(pac$mean) )
            if (length(afin) > 5 && var( pac$mean, na.rm=TRUE) > p$eps ) {
              ar1 = NULL
              ar1 = try( ar( pac$mean, order.max=1 ) )
              if (!inherits(ar1, "try-error")) {
                if ( length(ar1$ar) == 1 ) {
                  out["ar_1"] = ar1$ar
                }
              }
            }
          }
          if ( !is.finite(out[["ar_1"]]) ) {
            ar1 = try( cor( pac$mean[1:(length(piid) - 1)], pac$mean[2:(length(piid))], use="pairwise.complete.obs" ) )
            if (!inherits(ar1, "try-error")) out["ar_1"] = ar1
          }
        }

        ### Do the logistic model here ! -- if not already done ..
        if (!exists("ts_K", out)) {
          # model as a logistic with ts_r, ts_K, etc .. as stats outputs

        }

      }
      pac=piid=NULL
    }
    pac_i=NULL
  }

  # save stats
  for ( k in 1: length(p$statsvars) ) {
    if (exists( p$statsvars[k], out )) {
      S[Si,k] = out[[ p$statsvars[k] ]]
    }
  }

    # ----------------------
    # do last. it is an indicator of completion of all tasks
    # restarts would be broken otherwise
    Sflag[Si] = E[["complete"]]  # mark as complete without issues

  }  # end for loop

  return(NULL)

}
