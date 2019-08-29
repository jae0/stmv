
stmv_scale = function( ip=NULL, p, debugging=FALSE, eps = 1e-6, ... ) {
  #\\ core function to interpolate (model variogram) in parallel

  if (0) {
    # for debugging  runs ..
    currentstatus = stmv_statistics_status( p=p )
    p = parallel_run( p=p, runindex=list( locs=sample( currentstatus$todo )) )
    # parallel_run( stmv_scale, p=p, runindex=list( locs=sample( currentstatus$todo )) )
    p$runmode = "scale"
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


  #---------------------
  # data for modelling
  S = stmv_attach( p$storage.backend, p$ptr$S )
  Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )
  E = stmv_error_codes()

  Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )
  Y = stmv_attach( p$storage.backend, p$ptr$Y )

  if (p$nloccov > 0) Ycov = stmv_attach( p$storage.backend, p$ptr$Ycov )

  # misc intermediate calcs to be done outside of parallel loops

  # pre-calculate indices and dim for data to use inside the loop
  dat_names = unique( c(  p$variables$Y, p$variable$LOCS, p$variables$local_all,  "weights") )  # excludes p$variables$TIME
  # unless it is an explicit covariate and not a seasonal component there is no need for it
  # .. prediction grids create these from a time grid on the fly
  dat_nc = length( dat_names )
  iY = which(dat_names== p$variables$Y)
  ilocs = which( dat_names %in% p$variable$LOCS )


  nlogs = round( max(1, ifelse( length(ip) > p$nlogs*5, round(p$nlogs/5), length(ip) / p$nlogs ) ))
  logpoints  = ip[ round( seq( from=2, to=length(ip), length.out=nlogs ) ) ]
  if (length(logpoints) > 3) logpoints =  logpoints[ -c(1, length(logpoints)) ]  # drop first and last ones


  stmv_ntarget = sort( unique( c(1, p$stmv_nmin_downsize_factor) * p$stmv_nmax ), decreasing=TRUE )  # largest first
  stmv_ntarget = stmv_ntarget[ stmv_ntarget >=  p$stmv_nmin ]
  stmv_nmax = p$stmv_nmax

  stmv_distances = sort( unique( p$stmv_distance_scale ), decreasing=FALSE ) # smallest first

# main loop over each output location in S (stats output locations)
  for ( iip in ip ) {

    if ( iip %in% logpoints )  slog = stmv_logfile(p=p, flag= paste("Scale determination", p$runmode) )
    Si = p$runs[ iip, "locs" ]
    if ( Sflag[Si] != E[["todo"]] ) next()  # previously attempted .. skip
    if (debugging) print( paste("index =", iip, ";  Si = ", Si ) )
    # print( paste( iip, names(E)[match(Sflag[Si], unlist(E) ) ] ) )

    # obtain indices of data locations withing a given spatial range, optimally determined via variogram
    # find data nearest Sloc[Si,] and with sufficient data
    ndata = 0

    for ( ntarget in stmv_ntarget ) {
      for ( stmv_distance_cur in stmv_distances )  {
        # print(stmv_distance_cur)
        data_subset = NULL
        data_subset = stmv_select_data( p=p, Si=Si, localrange=stmv_distance_cur )

        if (is.null( data_subset )) next()
        unique_spatial_locations = data_subset$unique_spatial_locations
        ndata = length(data_subset$data_index)
        if ( unique_spatial_locations >= ntarget ) break()  # innermost loop
      }
      if ( unique_spatial_locations >= ntarget ) break() # middle loop
    }

    if (unique_spatial_locations < p$stmv_nmin ) {
      Sflag[Si] = E[["insufficient_data"]]
      if (debugging) print( paste("index =", iip, ";  insufficient data"  ) )
      next()   #not enough data
    }

    # if (p$stmv_variogram_method=="inla_nonseparable") {
    #   # TODO
    #   return()
    # }

    # generic ... crude separable approximations
    # spatial first

    o = NULL
    gc()

    for ( ss in 1:length(p$stmv_variogram_nbreaks_totry) ) {  # different number of breaks actually has an influence upon the stability of variograms
      o = try( stmv_variogram(
        xy=Yloc[data_subset$data_index,],
        z=Y[data_subset$data_index,],
        methods=p$stmv_variogram_method,
        distance_cutoff=stmv_distance_cur,
        discretized_n = round(stmv_distance_cur / p$pres),
        nbreaks=p$stmv_variogram_nbreaks_totry[ss]
      ) )
      if ( !is.null(o)) {
        if ( !inherits(o, "try-error")) {
          if ( exists(p$stmv_variogram_method, o)) break()
        }
      }
    }
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

    statvars_scale = c(
      sdTotal =sqrt( o$varZ),
      sdSpatial = sqrt(om$varSpatial) ,
      sdObs = sqrt(om$varObs),
      phi = om$phi,
      nu = om$nu,
      localrange = matern_phi2distance( phi=om$phi, nu=om$nu, cor=p$stmv_autocorrelation_localrange ),
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

  }  # end for loop


  return(NULL)

      if (0) {
        # stats
        # p$statsvars = c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "phi", "nu", "localrange" )
        S = stmv_attach( p$storage.backend, p$ptr$S )
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
