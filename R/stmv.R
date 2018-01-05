

stmv = function( p, runmode, DATA=NULL, storage.backend="bigmemory.ram",  debug_plot_variable_index=1 ) {

  if (0) {
    DATA=NULL
    storage.backend="bigmemory.ram"
    debug_plot_variable_index=1
  }
  
  #\\ localized modelling of space and time data to predict/interpolate upon a grid
  #\\ speed ratings: bigmemory.ram (1), ff (2), bigmemory.filebacked (3)
  # -----------------------------------------------------
  if (!exists("stmvSaveDir", p)) p$stmvSaveDir = file.path(p$data_root, "modelled", p$variables$Y, p$spatial.domain )

  if ( !file.exists(p$stmvSaveDir)) dir.create( p$stmvSaveDir, recursive=TRUE, showWarnings=FALSE )

  if ( "initialize" %in% runmode ) {

    message( "||| Initializing data files ... " )
    p$time.start = Sys.time()

    message( " ")
    message( "||| In case something should go wrong, intermediary outputs will be placed at:" )
    message( "|||",  p$stmvSaveDir  )
    message( " ")

    p = stmv_parameters(p=p) # fill in parameters with defaults where possible
    p = stmv_db( p=p, DS="filenames" )

    p$ptr = list() # location for data pointers

    # set up the data and problem using data objects
    tmpfiles = unlist( p$cache)
    for (tf in tmpfiles) if (file.exists( tf)) file.remove(tf)
    
    # permit passing a function rather than data directly .. less RAM usage in parent call
    if (is.null(DATA) ) {
      if (exists("DATA", p)) {
        if (class(p$DATA)=="character") {
          assign("DATA", eval(parse(text=p$DATA) ) )
        } else {
          DATA = p$DATA
        }
      }
    } else {
      if (class(DATA)=="character") {
        assign("DATA", eval(parse(text=DATA) ) )
      } else {
        # nothing to do as we assume DATA is a list
      }
    }

    if (is.null(DATA)) stop( "||| something went wrong with DATA inputs ... ")
    if (class(DATA) != "list") stop( "||| DATA should a list ... ")

    testvars = c(p$variables$Y, p$variables$COV, p$variables$TIME, p$variables$LOC)
    withdata = which(is.finite( (rowSums(DATA$input[, testvars] )) ) )
    if (length(withdata) < 1) stop( "Missing data or insufficient data")
    DATA$input = DATA$input[withdata, ]
    rm(withdata)

    # number of time slices
    if (!exists("nt", p)) {
      p$nt = 1  # default to 1 == no time
      if (exists( "ny", p)) p$nt = p$nt * p$ny  # annual time slices
      if (exists( "nw", p)) p$nt = p$nt * p$nw  # sub-annual time slices
    }

    # prediction times for space.annual methods, treat time as independent timeslices
    if ( !exists("prediction.ts", p)) p$prediction.ts = 1


    message(" ")
    message( "||| Initializing temporary storage of data and outputs files... ")
    message( "||| These are large files (4 to 6 X 5GB), it will take a minute ... ")
    stmv_db( p=p, DS="cleanup" )

    p$nloccov = 0
    if (exists("local_cov", p$variables)) p$nloccov = length(p$variables$local_cov)

    # construct prediction/output grid area ('pa')
    p$windowsize.half = floor(p$stmv_distance_prediction/p$pres) # convert distance to discretized increments of row/col indices

    if (exists("stmv_Y_transform", p)) {
      DATA$input[, p$variables$Y ] = p$stmv_Y_transform[[1]] (DATA$input[, p$variables$Y ] )
    }


    if (  exists("stmv_global_modelengine", p) ) {
      if ( p$stmv_global_modelengine !="none" ) {
      # to add global covariate model (hierarchical)
      # .. simplistic this way but faster ~ kriging with external drift
        if ("globalmodel" %in% runmode) {
          stmv_db( p=p, DS="global_model.redo", B=DATA$input )
        }
      }
    }

    # NOTE:: must not sink the following memory allocation into a deeper funcion as
    # NOTE:: bigmemory RAM seems to lose the pointers if they are not made simultaneously ?

    # init output data objects
    # statistics storage matrix ( aggregation window, coords ) .. no inputs required
    sbox = list(
      plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$stmv_distance_statsgrid ),
      plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$stmv_distance_statsgrid ) )

    # statistics coordinates
    Sloc = as.matrix( expand.grid( sbox$plons, sbox$plats ))
      if (p$storage.backend == "bigmemory.ram" ) {
        tmp_Sloc = big.matrix(nrow=nrow(Sloc), ncol=ncol(Sloc), type="double"  )
        tmp_Sloc[] = Sloc
        p$ptr$Sloc  = bigmemory::describe( tmp_Sloc  )
      }
      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$ptr$Sloc  = p$cache$Sloc
        bigmemory::as.big.matrix( Sloc, type="double", backingfile=basename(p$bm$Sloc), descriptorfile=basename(p$cache$Sloc), backingpath=p$stmvSaveDir )
      }
      if (p$storage.backend == "ff" ) {
        p$ptr$Sloc = ff( Sloc, dim=dim(Sloc), file=p$cache$Sloc, overwrite=TRUE )
      }
    rm( sbox )

    S = matrix( NaN, nrow=nrow(Sloc), ncol=length( p$statsvars ) ) # NA forces into logical
      if (p$storage.backend == "bigmemory.ram" ) {
        tmp_S = big.matrix(nrow=nrow(Sloc), ncol=length( p$statsvars ), type="double"  )
        tmp_S[] = S
        p$ptr$S  = bigmemory::describe( tmp_S )
      }
      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$ptr$S  = p$cache$S
        bigmemory::as.big.matrix( S, type="double", backingfile=basename(p$bm$S), descriptorfile=basename(p$cache$S), backingpath=p$stmvSaveDir )
      }
      if (p$storage.backend == "ff" ) {
        p$ptr$S = ff( S, dim=dim(S), file=p$cache$S, overwrite=TRUE )
      }


    Sflag = matrix( 0L, nrow=nrow(Sloc), ncol=1 )  # 0L is the todo flag
    # 0=to do
    # 1=complete
    # 2=oustide bounds(if any)
    # 3=shallow(if z is a covariate)
    # 4=range not ok, 
    # 5=skipped due to insufficient data, 
    # 6=skipped .. fast variogram did not work
    # 7=variogram estimated range not ok
    # 8=problem with prediction and/or modelling
    # 9=attempting ... if encountered then it was some general problem  or was interrrupted 
      if (p$storage.backend == "bigmemory.ram" ) {
        tmp_Sflag = big.matrix(nrow=nrow(Sloc), ncol=1, type="double" )
        tmp_Sflag[] = 0L # TODO flag
s        p$ptr$Sflag  = bigmemory::describe( tmp_Sflag )
      }
      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$ptr$Sflag  = p$cache$Sflag
        bigmemory::as.big.matrix( Sflag, type="double", backingfile=basename(p$bm$Sflag), descriptorfile=basename(p$cache$Sflag), backingpath=p$stmvSaveDir )
      }
      if (p$storage.backend == "ff" ) {
        p$ptr$Sflag = ff( Sflag, dim=dim(Sflag), file=p$cache$Sflag, overwrite=TRUE )
      }
    rm(S, Sflag, Sloc)

    # data to be worked upon .. either the raw data or covariate-residuals
    Ydata = as.matrix(DATA$input[, p$variables$Y ])
    if (exists("stmv_global_modelengine", p)) {
      if (p$stmv_global_modelengine !="none" ) {
        # at present only those that have a predict and residuals methods ... 
        covmodel = stmv_db( p=p, DS="global_model")
        Ypreds = predict(covmodel, type="link", se.fit=FALSE )  ## TODO .. keep track of the SE
        Ydata  = residuals(covmodel, type="deviance") # ie. link scale .. this is the default but make it explicit
        covmodel =NULL; gc()
      }
    }
    Ypreds = NULL
    Ydata = as.matrix( Ydata )
      if (p$storage.backend == "bigmemory.ram" ) {
        tmp_Y = big.matrix( nrow=nrow(Ydata), ncol=1, type="double"  )
        tmp_Y[] = Ydata
        p$ptr$Y  = bigmemory::describe( tmp_Y )
      }
      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$ptr$Y  = p$cache$Y
        bigmemory::as.big.matrix( Ydata, type="double", backingfile=basename(p$bm$Y), descriptorfile=basename(p$cache$Y), backingpath=p$stmvSaveDir )
      }
      if (p$storage.backend == "ff" ) {
        p$ptr$Y = ff( Ydata, dim=dim(Ydata), file=p$cache$Y, overwrite=TRUE )
      }
    rm(Ydata)

    Y = stmv_attach( p$storage.backend, p$ptr$Y )

   # data coordinates
    Yloc = as.matrix( DATA$input[, p$variables$LOCS ])
      if (p$storage.backend == "bigmemory.ram" ) {
        tmp_Yloc = big.matrix( nrow=nrow(Yloc), ncol=ncol(Yloc), type="double" )
        tmp_Yloc[] = Yloc
        p$ptr$Yloc = bigmemory::describe( tmp_Yloc )
      }
      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$ptr$Yloc  = p$cache$Yloc
        bigmemory::as.big.matrix( Yloc, type="double", backingfile=basename(p$bm$Yloc), descriptorfile=basename(p$cache$Yloc), backingpath=p$stmvSaveDir )
      }
      if (p$storage.backend == "ff" ) {
        p$ptr$Yloc = ff( Yloc, dim=dim(Yloc), file=p$cache$Yloc, overwrite=TRUE )
      }
    rm(Yloc)

    # independent variables/ covariate
    if (exists("COV", p$variables)) {
      Ycov = as.matrix(  DATA$input[ , p$variables$COV ] )
        if (p$storage.backend == "bigmemory.ram" ) {
          tmp_Ycov = big.matrix( nrow=nrow(Ycov), ncol=ncol(Ycov), type="double")
          tmp_Ycov[] = Ycov
          p$ptr$Ycov  = bigmemory::describe( tmp_Ycov )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Ycov  = p$cache$Ycov
          bigmemory::as.big.matrix( Ycov, type="double", backingfile=basename(p$bm$Ycov), descriptorfile=basename(p$cache$Ycov), backingpath=p$stmvSaveDir )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Ycov = ff( Ycov, dim=dim(Ycov), file=p$cache$Ycov, overwrite=TRUE )
        }
      rm(Ycov)
    }

    # data times
    if ( exists("TIME", p$variables) ) {
      Ytime = as.matrix(  DATA$input[, p$variables$TIME ] )
        if (p$storage.backend == "bigmemory.ram" ) {
          tmp_Ytime = big.matrix( nrow=nrow(Ytime), ncol=ncol(Ytime), type="double"  )
          tmp_Ytime[] = Ytime
          p$ptr$Ytime  = bigmemory::describe( tmp_Ytime )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Ytime  = p$cache$Ytime
          bigmemory::as.big.matrix( Ytime, type="double", backingfile=basename(p$bm$Ytime), descriptorfile=basename(p$cache$Ytime), backingpath=p$stmvSaveDir )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Ytime = ff( Ytime, dim=dim(Ytime), file=p$cache$Ytime, overwrite=TRUE )
        }
      Ytime =NULL; gc()
    }

    if (exists("COV", p$variables)) {
      # this needs to be done as Prediction covars need to be structured as lists
      if (!exists("Pcov", p$ptr) ) p$ptr$Pcov = list()
      tmp_Pcov = list()
      for ( covname in p$variables$COV ) {
        Pcovdata = as.matrix( DATA$output$COV[[covname]] )
        attr( Pcovdata, "dimnames" ) = NULL
        if (p$storage.backend == "bigmemory.ram" ) {
          tmp_Pcov[[covname]] = big.matrix( nrow=nrow(Pcovdata), ncol=ncol(Pcovdata), type="double"  )
          tmp_Pcov[[covname]][] = Pcovdata
          p$ptr$Pcov[[covname]]  = bigmemory::describe( tmp_Pcov[[covname]] )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Pcov[[covname]]  = p$cache$Pcov[[covname]]
          bigmemory::as.big.matrix( Pcovdata, type="double", backingfile=basename(p$bm$Pcov[[covname]]), descriptorfile=basename(p$cache$Pcov[[covname]]), backingpath=p$stmvSaveDir )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Pcov[[covname]] = ff( Pcovdata, dim=dim(Pcovdata), file=p$cache$Pcov[[covname]], overwrite=TRUE )
        }
        Pcovdata = NULL; gc()
      }
    }

    # predictions and associated stats
    P = matrix( NaN, nrow=nrow(DATA$output$LOCS), ncol=p$nt )
      # predictions
      if (p$storage.backend == "bigmemory.ram" ) {
        tmp_P = big.matrix( nrow=nrow(P), ncol=ncol(P), type="double" )
        tmp_P[] = P
        p$ptr$P  = bigmemory::describe( tmp_P )
      }
      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$ptr$P  = p$cache$P
        bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$P), descriptorfile=basename(p$cache$P), backingpath=p$stmvSaveDir )
      }
      if (p$storage.backend == "ff" ) {
        p$ptr$P = ff( P, dim=dim(P), file=p$cache$P, overwrite=TRUE )
      }

      # count of prediction estimates
      if (p$storage.backend == "bigmemory.ram" ) {
        tmp_Pn = big.matrix( nrow=nrow(P), ncol=ncol(P), type="double" )
        tmp_Pn[] = P
        p$ptr$Pn = bigmemory::describe( tmp_Pn )
      }
      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$ptr$Pn  = p$cache$Pn
        bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$Pn), descriptorfile=basename(p$cache$Pn), backingpath=p$stmvSaveDir )
      }
      if (p$storage.backend == "ff" ) {
        p$ptr$Pn = ff( P, dim=dim(P), file=p$cache$Pn, overwrite=TRUE )
      }

      # sd of prediction estimates
      if (p$storage.backend == "bigmemory.ram" ) {
        tmp_Psd = big.matrix( nrow=nrow(P), ncol=ncol(P), type="double" )
        tmp_Psd[] = P
        p$ptr$Psd =bigmemory::describe( tmp_Psd )
      }
      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$ptr$Psd  = p$cache$Psd
        bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$Psd), descriptorfile=basename(p$cache$Psd), backingpath=p$stmvSaveDir )
      }
      if (p$storage.backend == "ff" ) {
        p$ptr$Psd = ff( P, dim=dim(P), file=p$cache$Psd, overwrite=TRUE )
      }

    # prediction coordinates
    Ploc = as.matrix( DATA$output$LOCS )
      attr( Ploc, "dimnames" ) = NULL
      if (p$storage.backend == "bigmemory.ram" ) {
        tmp_Ploc = big.matrix( nrow=nrow(Ploc), ncol=ncol(Ploc), type="double" )
        tmp_Ploc[] = Ploc
        p$ptr$Ploc  = bigmemory::describe( tmp_Ploc )
      }
      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$ptr$Ploc  = p$cache$Ploc
        bigmemory::as.big.matrix( Ploc, type="double", backingfile=basename(p$bm$Ploc), descriptorfile=basename(p$cache$Ploc), backingpath=p$stmvSaveDir )
      }
      if (p$storage.backend == "ff" ) {
        p$ptr$Ploc = ff( Ploc, dim=dim(Ploc), file=p$cache$Ploc, overwrite=TRUE )
      }
    Ploc = DATA = NULL; gc()

    if (exists("stmv_global_modelengine", p) ) {
      if (p$stmv_global_modelengine !="none" ) {
        # create prediction suface with covariate-based additive offsets
        if (p$storage.backend == "bigmemory.ram" ) {
          tmp_P0= big.matrix( nrow=nrow(P), ncol=ncol(P) , type="double" )
          tmp_P0[] = P
          p$ptr$P0 = bigmemory::describe(tmp_P0 )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$P0  = p$cache$P0
          bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$P0), descriptorfile=basename(p$cache$P0), backingpath=p$stmvSaveDir )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$P0 = ff( P, dim=dim(P), file=p$cache$P0, overwrite=TRUE )
        }

        if (p$storage.backend == "bigmemory.ram" ) {
          tmp_P0sd= big.matrix( nrow=nrow(P), ncol=ncol(P) , type="double" )
          tmp_P0sd[] = P
          p$ptr$P0sd = bigmemory::describe(tmp_P0sd )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$P0sd  = p$cache$P0sd
          bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$P0sd), descriptorfile=basename(p$cache$P0sd), backingpath=p$stmvSaveDir )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$P0sd = ff( P, dim=dim(P), file=p$cache$P0sd, overwrite=TRUE )
        }
        P=NULL; gc()

        # test to see if all covars are static as this can speed up the initial predictions
        message(" ")
        message( "||| Predicting global effect of covariates at each prediction location ... ")
        message( "||| depending upon the size of the prediction grid and number of cpus (~1hr?).. ")

        p$timec_covariates_0 =  Sys.time()
        nc_cov =NULL
        for (i in p$variables$COV ) {
          pu = stmv_attach( p$storage.backend, p$ptr$Pcov[[i]] )
          nc_cov = c( nc_cov,  ncol(pu) )
        }
        p$all.covars.static = ifelse( any(nc_cov > 1),  FALSE, TRUE )
        pc = p # copy
        if (!pc$all.covars.static) if (exists("clusters.covars", pc) ) pc$clusters = pc$clusters.covars
        # takes about 28 GB per run .. adjust cluster number temporarily
        suppressMessages( stmv_db( p=pc, DS="global.prediction.surface" ) )
        p$time_covariates = round(difftime( Sys.time(), p$timec_covariates_0 , units="hours"), 3)
        message( paste( "||| Time taken to predict covariate surface (hours):", p$time_covariates ) )
      }
    }
    P = NULL; gc() # yes, repeat in case covs are not modelled

    stmv_db( p=p, DS="statistics.Sflag" )
    Y = stmv_attach( p$storage.backend, p$ptr$Y )
    Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )
    Yi = 1:length(Y) # index with useable data
    bad = which( !is.finite( Y[]))
    if (length(bad)> 0 ) Yi[bad] = NA

    # data locations
    bad = which( !is.finite( rowSums(Yloc[])))
    if (length(bad)> 0 ) Yi[bad] = NA

  # data locations
    if (exists("COV", p$variables)) {
      Ycov = stmv_attach( p$storage.backend, p$ptr$Ycov )
      if (length(p$variables$COV)==1) {
        bad = which( !is.finite( Ycov[] ))
      } else {
        bad = which( !is.finite( rowSums(Ycov[])))
      }
      if (length(bad)> 0 ) Yi[bad] = NA
      Yi = na.omit(Yi)
    }

    # data locations
    if (exists("TIME", p$variables)) {
      Ytime = stmv_attach( p$storage.backend, p$ptr$Ytime )
      bad = which( !is.finite( Ytime[] ))
      if (length(bad)> 0 ) Yi[bad] = NA
      Yi = na.omit(Yi)
    }
    bad = NULL

    Yi = as.matrix(Yi)
      if (p$storage.backend == "bigmemory.ram" ) {
        tmp_Yi = big.matrix( nrow=nrow(Yi), ncol=ncol(Yi), type="double" )
        tmp_Yi[] = Yi
        p$ptr$Yi  = bigmemory::describe( tmp_Yi )
      }
      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$ptr$Yi  = p$cache$Yi
        bigmemory::as.big.matrix( Yi, type="double", backingfile=basename(p$bm$Yi), descriptorfile=basename(p$cache$Yi), backingpath=p$stmvSaveDir )
      }
      if (p$storage.backend == "ff" ) {
        p$ptr$Yi = ff( Yi, dim=dim(Yi), file=p$cache$Yi, overwrite=TRUE )
      }
    rm(Yi)

    if ( !exists("stmv_distance_scale", p)) {
      Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )
      p$stmv_distance_scale = min( diff(range( Yloc[,1]) ), diff(range( Yloc[,2]) ) ) / 10
      message( paste( "||| Crude distance scale:", p$stmv_distance_scale, "" ) )
    }
    if ( !exists("stmv_distance_min", p)) p$stmv_distance_min = mean( c(p$stmv_distance_prediction, p$stmv_distance_scale /20 ) )
    if ( !exists("stmv_distance_max", p)) p$stmv_distance_max = mean( c(p$stmv_distance_prediction*10, p$stmv_distance_scale * 2 ) )

    message("||| Finished. ")
    message("||| Once analyses begin, you can view maps from an external R session: ")
    message("|||   p = stmv_db( p=p, DS='load.parameters' )" ) 
    message("|||   stmv(p=p, runmode='debug_pred_static_map', debug_plot_variable_index=1) ")
    message("|||   stmv(p=p, runmode='debug_pred_static_log_map', debug_plot_variable_index=1)")
    message("|||   stmv(p=p, runmode='debug_pred_dynamic_map', debug_plot_variable_index=1)")
    message("|||   stmv(p=p, runmode='debug_stats_map', debug_plot_variable_index=1)")
    message("||| Monitor the status of modelling by looking at the output of the following file:")
    message("||| in linux, you can issue the following command:" )
    message("|||   watch -n 60 cat ",  p$stmv_current_status  )

    p$initialized = TRUE
    p <<- p  # push to parent in case a manual restart is needed
    
    stmv_db( p=p, DS="save.parameters" )  # save in case a restart is required .. mostly for the pointers to data objects
    gc()

  } else {

    if (!exists("time.start", p) ) p$time.start = Sys.time()
    message( "||| Continuing from an interrupted start ..." )
    if (!p$initialized) {
      message( "||| Loading parameters from a saved configuration:", file.path( p$stmvSaveDir, 'p.rdata' ) )
      p = stmv_db( p=p, DS="load.parameters" )
      p <<- p  # push to parent in case a manual restart is needed
    }
  }  # end of intialization of data structures

  if (!exists("initialized", p)) stop( "||| stmv was not initialized properly" )
  if (!p$initialized) stop( "||| stmv was not initialized properly" )

  # -----------------------------------------------------
  if ( "debug" %in% runmode ) {
    if (!p$initialized) {
      message( "||| Loading parameters from a saved configuration:", file.path( p$stmvSaveDir, 'p.rdata' ) )
      p = stmv_db( p=p, DS="load.parameters" )
      p <<- p  # push to parent in case a manual restart is needed
    }
    currentstatus = stmv_db( p=p, DS="statistics.status" )
    p = parallel_run( p=p, runindex=list( locs=sample( currentstatus$todo )) ) # random order helps use all cpus
    print( c( unlist( currentstatus[ c("n.total", "n.shallow", "n.todo", "n.skipped", "n.outside", "n.complete" ) ] ) ) )
    message( "||| Entering browser mode ...")
    browser()
    stmv_interpolate (p=p )
  }
    
  # -----------------------------------------------------

  if ( "debug_pred_static_map" %in% runmode) {
      Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
      P = stmv_attach( p$storage.backend, p$ptr$P )
      lattice::levelplot( (P[,debug_plot_variable_index])~Ploc[,1]+Ploc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso")
  }

  # -----------------------------------------------------
  if ( "debug_pred_static_log_map" %in% runmode) {
      Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
      P = stmv_attach( p$storage.backend, p$ptr$P )
      lattice::levelplot( log(P[,debug_plot_variable_index])~Ploc[,1]+Ploc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso")
  }

  # -----------------------------------------------------
  if ( "debug_pred_dynamic_map" %in% runmode) {
      Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
      P = stmv_attach( p$storage.backend, p$ptr$P )
      for (i in 1:p$nt) {
        print( lattice::levelplot( P[,i] ~ Ploc[,1] + Ploc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
      }
  }

  # -----------------------------------------------------
  if ( "debug_stats_map" %in% runmode) {
      Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
      S = stmv_attach( p$storage.backend, p$ptr$S )
      lattice::levelplot(S[,debug_plot_variable_index]~Sloc[,1]+Sloc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso")
  }



  # -----------------------------------------------------
  if ( "stage1" %in% runmode ) {
    # this is the basic run
    timei1 =  Sys.time()
    currentstatus = stmv_db( p=p, DS="statistics.status" )
    # p <<- p  # push to parent in case a manual restart is possible
    suppressMessages( parallel_run( stmv_interpolate, p=p, runindex=list( locs=sample( currentstatus$todo )) )) # random order helps use all cpus )
    p$time_default = round( difftime( Sys.time(), timei1, units="hours" ), 3 )
    message(" ")
    message( paste( "||| Time taken to complete stage 1 interpolations (hours):", p$time_default, "" ) )
    currentstatus = stmv_db( p=p, DS="statistics.status" )
    print( c( unlist( currentstatus[ c("n.total", "n.shallow", "n.todo", "n.skipped", "n.outside", "n.complete" ) ] ) ) )
    gc()
  }


  # -----------------------------------------------------
  if ( "stage2" %in% runmode ) {
    timei2 =  Sys.time()
    message(" ")
    message( "||| Starting stage 2: more permisssive distance settings:")
    message( "|||  i.e., (stmv_distance_max and stmv_distance_scale) X {", paste0(p$stmv_multiplier_stage2, collapse=" "), "}" )

    for ( mult in p$stmv_multiplier_stage2 ) {
      currentstatus = stmv_db(p=p, DS="statistics.status.reset" )
      if (length(currentstatus$todo) > 0) {
        suppressMessages( parallel_run( stmv_interpolate, p=p, 
          stmv_distance_max=p$stmv_distance_max*mult, 
          stmv_distance_scale=p$stmv_distance_scale*mult,
          runindex=list( locs=sample( currentstatus$todo ) ) # random order helps use all cpus
        ))
      }
    }
    p$time_stage2 = round( difftime( Sys.time(), timei2, units="hours" ), 3)
    message(" ")
    message( paste( "||| Time taken to complete stage 2 interpolations (hours):", p$time_stage2, "" ) )
    currentstatus = stmv_db( p=p, DS="statistics.status" )
    print( c( unlist( currentstatus[ c("n.total", "n.shallow", "n.todo", "n.skipped", "n.outside", "n.complete" ) ] ) ) )
    gc()
  }

  # -----------------------------------------------------
  if ( "stage3" %in% runmode ) {
    timei3 =  Sys.time()
    message(" ")
    message( "||| Starting stage 3: simple TPS-based failsafe method to interpolate all the remaining locations " )
    toredo = stmv_db( p=p, DS="flag.incomplete.predictions" )
    if ( !is.null(toredo) && length(toredo) > 0) {
      Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )
      Sflag[toredo]=0L
      p$stmv_local_modelengine = "tps"
      # ?? why is this here ? p = aegis_parameters( p=p, DS="bathymetry" )
      parallel_run( stmv_interpolate, p=p, runindex=list( locs=sample( toredo )) ) # random order helps use all cpus
    }
    p$time_stage3 = round( difftime( Sys.time(), timei3, units="hours" ), 3)
    message( paste( "||| Time taken to complete stage3 interpolations (hours):", p$time_stage3, "" ) )
    currentstatus = stmv_db( p=p, DS="statistics.status" )
    print( c( unlist( currentstatus[ c("n.total", "n.shallow", "n.todo", "n.skipped", "n.outside", "n.complete" ) ] ) ) )
    p <<- p  # push to parent in case a manual restart is needed
    gc()
  }

  # save again, in case some timings/etc needed in a restart
  # p <<- p  # push to parent in case a manual restart is possible

  # -----------------------------------------------------
  if ( "save" %in% runmode ) {
    message(" ")
    message( "||| Saving predictions to disk .. " )
    stmv_db( p=p, DS="stmv.prediction.redo" ) # save to disk for use outside stmv*, returning to user scale
    message( "||| Saving statistics to disk .. " )
    stmv_db( p=p, DS="stats.to.prediction.grid.redo") # save to disk for use outside stmv*
    message ("||| Finished! ")
    stmv_db( p=p, DS="save.parameters" )  # save in case a restart is required .. mostly for the pointers to data objects

    if ( p$storage.backend !="bigmemory.ram" ) {
      resp = readline( "||| Delete temporary files? Type to confirm <YES>:  ")
      if (resp=="YES") {
        stmv_db( p=p, DS="cleanup" )
      } else {
        message(" ")
        message( "||| Leaving temporary files alone in case you need to examine them or restart a process. ")
        message( "||| You can delete them by running: stmv_db( p=p, DS='cleanup' ), once you are done. ")
      }
    }
    p <<- p  # push to parent in case a manual restart is needed
  }

  p$time_total = round( difftime( Sys.time(), p$time.start, units="hours" ),3)
  message( paste( "||| Time taken for full analysis (hours):", p$time_total ) )
  message( paste( "||| Your parameter 'p' has been updated in case you need to re-run something" ) )
  p <<- p
}
