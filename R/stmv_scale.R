

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

    Sflag[Si] = E[["todo"]]

    if (ndata < p$n.min) {
      Sflag[Si] = E[["insufficient_data"]]
      next()   #not enough data
    }

    # NOTE: this range is a crude estimate that averages across years (if any) ...
    o = NULL

    if (p$stmv_variogram_method =="inla_space") {

      o = try( stmv_variogram(
        xy=Yloc[Yi[U],],
        z=Y[Yi[U],],
        methods="inla",
        distance_cutoff=stmv_distance_cur,
        nbreaks=15,
        range_correlation=p$stmv_range_correlation
      ) )

    } else if (p$stmv_variogram_method =="inla_space_year") {


    } else {

      o = try( stmv_variogram(
        xy=Yloc[Yi[U],],
        z=Y[Yi[U],],
        methods=p$stmv_variogram_method,
        distance_cutoff=stmv_distance_cur,
        nbreaks=15,
        range_correlation=p$stmv_range_correlation
      ) )

    }


    if ( is.null(o)) {
      Sflag[Si] = E[["variogram_failure"]]
      next()
    }

    if ( inherits(o, "try-error")) {
      Sflag[Si] = E[["variogram_failure"]]
      next()
    }

    if ( !exists(p$stmv_variogram_method, o)) {
      Sflag[Si] =  E[["variogram_failure"]]
      next()
    }

    if ( !exists("range_ok", o[[p$stmv_variogram_method]]) ) {
      Sflag[Si] =  E[["variogram_range_limit"]]
      next()     #
    }

    if ( !o[[p$stmv_variogram_method]][["range_ok"]] ) {
      # retain crude estimate and run with it
      Sflag[Si] =  E[["variogram_range_limit"]]
      next()
    }

    # save stats
    statvars_scale = c("sdTotal", "sdSpatial" ,"sdObs", "range", "phi", "nu" )
    v = match( p$statsvars, statvars_scale )
    u = which( is.finite(v) )
    for ( ll in 1:length(u) ) {
      k = u[ll]
      if (exists( statvars_scale[ll], out )) {
        S[Si,k] = out[[ statvars_scale[ll] ]]
      }
    }

    # ----------------------
    # do last. it is an indicator of completion of all tasks
    # restarts would be broken otherwise
    Sflag[Si] = E[["complete"]]  # mark as complete without issues

  }  # end for loop

  return(NULL)

}
