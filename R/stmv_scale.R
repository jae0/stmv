
stmv_scale = function( ip=NULL, p, stmv_localrange, debugging=FALSE, eps=1e-6, ... ) {
  #\\ core function to interpolate (model variogram) in parallel

  if (0) {
    # for debugging  runs ..
    currentstatus = stmv_statistics_status( p=p )
    p = parallel_run( p=p, runindex=list( locs=sample( currentstatus$todo )) )
    # parallel_run( stmv_scale, p=p, runindex=list( locs=sample( currentstatus$todo )) )
    ip = 1:p$nruns
    debugging=TRUE
    eps = 1e-6
  }

  p = parameters_add( p, list(...) ) # add passed args to parameter list, priority to args

  if (exists( "libs", p)) suppressMessages( RLibrary( p$libs ) )

  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns


  #---------------------
  # data for modelling
  S = stmv_attach( p$storage_backend, p$ptr$S )
  Sflag = stmv_attach( p$storage_backend, p$ptr$Sflag )
  E = stmv_error_codes()

  Yloc = stmv_attach( p$storage_backend, p$ptr$Yloc )
  Y = stmv_attach( p$storage_backend, p$ptr$Y )

  if (p$nloccov > 0) Ycov = stmv_attach( p$storage_backend, p$ptr$Ycov )

  # misc intermediate calcs to be done outside of parallel loops

  # pre-calculate indices and dim for data to use inside the loop
  dat_names = unique( c(  p$stmv_variables$Y, p$stmv_variables$LOCS, p$stmv_variables$local_all,  "weights") )  # excludes p$stmv_variables$TIME
  # unless it is an explicit covariate and not a seasonal component there is no need for it
  # .. prediction grids create these from a time grid on the fly
  dat_nc = length( dat_names )
  iY = which(dat_names== p$stmv_variables$Y)
  ilocs = which( dat_names %in% p$stmv_variables$LOCS )

  if (length(ip) < 100) {
    nlogs = length(ip) / 5
  } else {
    nlogs = ifelse( length(ip) > (p$nlogs*5), p$nlogs, length(ip) / 5  )
  }
  logpoints  =  sort( sample( ip, round( max(1, nlogs) ) ) )  # randomize
  savepoints =  sort( sample( ip,  3 ) )  

  stmv_ntarget = round( sort( unique( c(1, p$stmv_nmin_downsize_factor) * p$stmv_nmax ), decreasing=TRUE ) ) # largest first
  stmv_ntarget = stmv_ntarget[ stmv_ntarget >=  p$stmv_nmin ]
  stmv_nmax = p$stmv_nmax

# main loop over each output location in S (stats output locations)
  for ( iip in ip ) {
 
    stmv_control_check(p=p)
  
    if ( iip %in% logpoints )  slog = stmv_logfile(p=p, flag="Scale determination")
    if ( iip %in% savepoints )  {
      stmv_db(p=p, DS="save_current_state", runmode=p$current_runmode, datasubset=c("P", "Pn", "Psd", "statistics") )  
    }

    Si = p$runs[ iip, "locs" ]
    if ( Sflag[Si] != E[["todo"]] ) next()  # previously attempted .. skip
    if (debugging) print( paste("index =", iip, ";  Si = ", Si ) )
    # print( paste( iip, names(E)[match(Sflag[Si], unlist(E) ) ] ) )

    # obtain indices of data locations withing a given spatial range, optimally determined via variogram
    # find data nearest Sloc[Si,] and with sufficient data
    ndata = 0
    unique_spatial_locations = 0

    for ( ntarget in stmv_ntarget ) {
        data_subset = NULL
        data_subset = stmv_select_data( p=p, Si=Si, localrange=stmv_localrange )
        if (is.null( data_subset )) next()
        unique_spatial_locations = data_subset$unique_spatial_locations
        ndata = length(data_subset$data_index)
        if ( unique_spatial_locations >= ntarget ) break()   
    }

    # check again
    if ( unique_spatial_locations < ntarget ) next()   

    # if (p$stmv_variogram_method=="inla_nonseparable") {
    #   # TODO
    #   return()
    # }

    # generic ... crude separable approximations
    # spatial first

    o = NULL
    gc()

      
    o = try( stmv_variogram(
      xy=Yloc[data_subset$data_index,],
      z=Y[data_subset$data_index,],
      methods=p$stmv_variogram_method,
      distance_cutoff=stmv_localrange,  # initial guess of effective range
      # discretized_n = trunc(stmv_localrange / p$pres),
      nbreaks=p$stmv_variogram_nbreaks_totry # different number of breaks actually has an influence upon the stability of variograms
    ) )
 
    data_subset  = NULL
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

    if ( o$stmv_internal_variance_total < eps ) return(NULL)
    if ( !om$phi_ok ) return(NULL)

    statvars_scale = c(
      sdTotal = sqrt( o$stmv_internal_variance_total),
      rsquared=NA,
      sdSpatial = sqrt(om$varSpatial) ,
      sdObs = sqrt(om$varObs),
      phi = om$phi,
      nu = om$nu,
      localrange = matern_phi2distance( phi=om$phi, nu=om$nu, cor=p$stmv_autocorrelation_localrange ),
      ndata=ndata
    )


    for ( vv in 1:length(statvars_scale) ) {
      vn = names(statvars_scale)[vv]
      vi = match(vn, p$statsvars)
      if ( is.finite(vi)) {
        if ( is.finite( statvars_scale[[ vn ]] ) ) {
          S[Si, vi] = statvars_scale[[ vn ]]
        }
      }
    }

    # temporal
    if (p$dimensionality =="space") {
       # nothing to do
    }

    if (p$dimensionality =="space-time") {

      if (0) {
        
        dyear_centre = p$dyears[ trunc(p$nw/2) ] + p$tres/2  
        
        # annual ts, seasonally centered and spatially
        ar_timerange = NA
        ar_1 = NA

        pac = res$predictions[ pac_i, ]
        pac$dyr = pac[, p$stmv_variables$TIME] - trunc(pac[, p$stmv_variables$TIME] )
        piid = which( zapsmall( pac$dyr - dyear_centre) == 0 )
        pac = pac[ piid, c(p$stmv_variables$TIME, "mean")]
        pac = pac[ order(pac[,p$stmv_variables$TIME]),]
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

    if (p$dimensionality=="space-time-cyclic")  {

      if (0) {

        dyear_centre = p$dyears[ trunc(p$nw/2) ] + p$tres/2  

        # annual ts, seasonally centered and spatially
        ar_timerange = NA
        ar_1 = NA

        pac = res$predictions[ pac_i, ]
        pac$dyr = pac[, p$stmv_variables$TIME] - trunc(pac[, p$stmv_variables$TIME] )
        piid = which( zapsmall( pac$dyr - dyear_centre) == 0 )
        pac = pac[ piid, c(p$stmv_variables$TIME, "mean")]
        pac = pac[ order(pac[,p$stmv_variables$TIME]),]
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

  }  # end for loop


  return(NULL)

      if (0) {
        # stats
        # p$statsvars = c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "phi", "nu", "localrange" )
        S = stmv_attach( p$storage_backend, p$ptr$S )
        sbox = list(
          plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$stmv_distance_statsgrid ),
          plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$stmv_distance_statsgrid ) )
        # statistics coordinates
        locations = as.matrix( expand_grid_fast( sbox$plons, sbox$plats ))
        levelplot(S[,match("localrange", p$statsvars )]~ locations[,1]+locations[,2])
        levelplot(S[,match("nu", p$statsvars )]~ locations[,1]+locations[,2])
        levelplot(S[,match("sdTotal", p$statsvars )]~ locations[,1]+locations[,2])
        levelplot(S[,match("sdSpatial", p$statsvars )]~ locations[,1]+locations[,2])

      }



}
