

stmv_interpolate_predictions = function( ip=NULL, p, localrange_interpolation=NULL, debugging=FALSE, ... ) {
  #\\ simple brute force linear interpolaion of finalized prediction to fill missing data

  if (0) {
    # for debugging  runs ..
    currentstatus = stmv_statistics_status( p=p )
    p = parallel_run( p=p, runindex=list( locs=sample( currentstatus$todo )) )
    ip = 1:p$nruns
    debugging=TRUE
  }


  p = parameters_add(p, list(...)) # add passed args to parameter list, priority to args

  if (exists( "libs", p)) suppressMessages( RLibrary( p$libs ) )

  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns


  #---------------------
  # data for modelling

  P = stmv_attach( p$storage_backend, p$ptr$P )
  Psd = stmv_attach( p$storage_backend, p$ptr$Psd )
  Ploc = stmv_attach( p$storage_backend, p$ptr$Ploc )

  S = stmv_attach( p$storage_backend, p$ptr$S )
  Sflag = stmv_attach( p$storage_backend, p$ptr$Sflag )
  Sloc = stmv_attach( p$storage_backend, p$ptr$Sloc )

  E = stmv_error_codes()

  if (length(ip) < 100) {
    nlogs = length(ip) / 5
  } else {
    nlogs = ifelse( length(ip) > (p$nlogs*5), p$nlogs, length(ip) / 5  )
  }
  logpoints  =  sort( sample( ip, round( max(1, nlogs) ) ) )  # randomize

# main loop over each output location in S (stats output locations)
  for ( iip in ip ) {
    
    stmv_control_check(p=p)
  
    if ( iip %in% logpoints )  slog = stmv_logfile(p=p, flag= paste("Interpolation", p$runoption) )
    Si = p$runs[ iip, "locs" ]

    print( paste("index =", iip, ";  Si = ", Si ) )
    if ( Sflag[Si] == E[["complete"]] ) next()

    # construct prediction/output grid area ('pa')
    # convert distance to discretized increments of row/col indices;
    if (exists("stmv_distance_prediction_limits", p)) {
      localrange_interpolation = min( max( localrange_interpolation, min(p$stmv_distance_prediction_limits) ), max(p$stmv_distance_prediction_limits), na.rm=TRUE )
    }

    windowsize.half = aegis_floor( localrange_interpolation / p$pres ) + 1L
    # construct data (including static covariates) for prediction locations (pa)
    pa = try( stmv_predictionarea_lattice( p=p, sloc=Sloc[Si,], windowsize.half=windowsize.half ) )
    if ( is.null(pa) ) {
      Sflag[Si] = E[["prediction_area"]]
      next()
    }

    if ( inherits(pa, "try-error") ) {
      pa = NULL
      Sflag[Si] = E[["prediction_area"]]
      if (debugging) message("Error: prediction grid ... try-error .. this should not happen.  check this")
      next()
    }

    if ( exists("TIME", p$stmv_variables) )  pa = try( stmv_predictiontime( p=p, pa=pa ) ) # add time to pa and time varying covars
    if ( is.null(pa) ) {
      Sflag[Si] = E[["prediction_time"]]
      next()
    }

    if ( inherits(pa, "try-error") ) {
      pa = NULL
      Sflag[Si] = E[["prediction_time"]]
      if (debugging) message("Error: prediction grid ... try-error .. this should not happen.  check this")
      next()
    }


    #direct interpolation upon P ..

    for ( ti in 1:p$nt ) {

      if ( exists("TIME", p$stmv_variables) ) {
        pa_i = which( pa[, p$stmv_variables$TIME] == p$prediction_ts[ti] )
        if (length(xi) < 5 ) next()
      } else {
        pa_i = 1:nrow(pa)
      }

      ee = pa$i[pa_i]

      loc = Ploc[ ee, ]
      dat = P[ ee ]

      subdomain_withdata = which(is.finite( dat ))
      subdomain_withoutdata = which(!is.finite( dat ))

      if (length(subdomain_withoutdata) == 0) next()
      if (length(subdomain_withdata) < 5 ) next()

      X = try ( interp::interp(
        x=loc[ subdomain_withdata, 1 ],
        y=loc[ subdomain_withdata, 2 ],
        z= dat[ subdomain_withdata  ],
        xo=loc[ subdomain_withoutdata, 1 ],
        yo=loc[ subdomain_withoutdata, 2 ],
        input="points",
        output="points",
        method="linear",
        extrap=TRUE
      )$z, silent =TRUE )
      if ( inherits(X, "try-error") ) next()
      P[ ee[subdomain_withoutdata], ti ] = X

      dat = Psd[ pa_i ]
      subdomain_withdata = which(is.finite( dat ))
      subdomain_withoutdata = which(!is.finite( dat ))
      X = try( interp::interp(
        x=loc[ subdomain_withdata, 1 ],
        y=loc[ subdomain_withdata, 2 ],
        z= dat[ subdomain_withdata  ],
        xo=loc[ subdomain_withoutdata, 1 ],
        yo=loc[ subdomain_withoutdata, 2 ],
        input="points",
        output="points",
        method="linear",
        extrap=TRUE
      )$z, silent =TRUE )
      if ( inherits(X, "try-error") ) next()
      Psd[ ee[subdomain_withoutdata], ti ] = X

      X = NULL
      pa  =  NULL
      ee = NULL
    }

    # ----------------------
    # do last. it is an indicator of completion of all tasks
    # restarts would be broken otherwise
    Sflag[Si] = E[["complete"]]  # mark as complete without issues

  }  # end for loop


  return(NULL)

}
