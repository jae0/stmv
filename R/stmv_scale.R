

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
  YY1 = Yloc[Yi[],1]
  YY2 = Yloc[Yi[],2]

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
    dlon = abs( Sloc[Si,1] - YY1 )
    dlat = abs( Sloc[Si,2] - YY2 )
    ndata = 0
    for ( stmv_distance_cur in distance_to_upsample )  {
      U = which( {dlon  <= stmv_distance_cur} & {dlat <= stmv_distance_cur} )  # faster to take a block
      ndata = length(U)
      if ( ndata >= p$n.min ) break()
    }
    dlon = NULL
    dlat = NULL

    Sflag[Si] = E[["todo"]]

    if (ndata < p$n.min) {
      Sflag[Si] = E[["insufficient_data"]]
      if (debugging) print( paste("index =", iip, ";  insufficient data"  ) )
      next()   #not enough data
    }

    # NOTE: this range is a crude estimate that averages across years (if any) ...
    o = NULL

    if (p$stmv_variogram_method =="inla") {

      if (p$stmv_variogram_method =="space") {

        o = try( stmv_variogram(
          xy=Yloc[Yi[U],],
          z=Y[Yi[U],],
          methods="inla",
          distance_cutoff=stmv_distance_cur,
          nbreaks=15,
          range_correlation=p$stmv_range_correlation #,  plotdata=TRUE
        ) )

      } else if (p$stmv_variogram_method =="space-year") {

      } else if (p$sp$stmv_dimensionality=="space-year-season")  {

      }

    } else {

      o = try( stmv_variogram(
        xy=Yloc[Yi[U],],
        z=Y[Yi[U],],
        methods=p$stmv_variogram_method,
        distance_cutoff=stmv_distance_cur,
        nbreaks=15,
        range_correlation=p$stmv_range_correlation # ,  plotdata=TRUE
      ) )

    }

    U = NULL

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

    # save stats
    om  = o[[p$stmv_variogram_method]]

    statvars_scale = c(
      sdTotal =sqrt( o$varZ),
      sdSpatial = sqrt(om$varSpatial) ,
      sdObs = sqrt(om$varObs),
      range = om$range,
      phi = om$phi,
      nu = om$nu,
      ndata=ndata
    )

    S[Si,match( names(statvars_scale), p$statsvars )] = statvars_scale

    if (debugging) {
      print( paste("index =", iip, ";  Sflag = ", names(E)[match(Sflag[Si], E)]  ) )
    }

  }  # end for loop


  return(NULL)

}
