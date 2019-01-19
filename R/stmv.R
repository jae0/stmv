

stmv = function( p, runmode="interpolate", DATA=NULL,
  use_saved_state=NULL, save_completed_data=TRUE, force_complete_solution=TRUE, nlogs=200,
  debug_plot_variable_index=1, debug_data_source="saved.state", debug_plot_log=FALSE ) {

  if (0) {
    nlogs = 25
    force_complete_solution=FALSE
    use_saved_state=NULL # or "disk"
    DATA=NULL
    storage.backend="bigmemory.ram"
    save_completed_data=TRUE  # export out of stmv system for use outside (e.g., by aegis)
    debug_plot_variable_index=1
    runmode=c("interpolate", "globalmodel")
  }

  #\\ localized modelling of space and time data to predict/interpolate upon a grid
  #\\ speed ratings: bigmemory.ram (1), ff (2), bigmemory.filebacked (3)

  # -----------------------------------------------------

  if (!exists("time_start", p) ) p$time_start = Sys.time()
  if (!exists("stmvSaveDir", p)) p$stmvSaveDir = file.path(p$data_root, "modelled", p$variables$Y, p$spatial.domain )
  if ( !file.exists(p$stmvSaveDir)) dir.create( p$stmvSaveDir, recursive=TRUE, showWarnings=FALSE )

  message( "||| Initializing data files ... " )
  message( "||| In case something should go wrong, intermediary outputs will be placed at:" )
  message( "|||",  p$stmvSaveDir  )
  message( " ")

  p$nlogs = nlogs
  p = stmv_parameters(p=p) # fill in parameters with defaults where required
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
  withdata=NULL

  # number of time slices
  if (!exists("nt", p)) {
    p$nt = 1  # default to 1 == no time
    if (exists( "ny", p)) p$nt = p$nt * p$ny  # annual time slices
    if (exists( "nw", p)) p$nt = p$nt * p$nw  # sub-annual time slices
  }

  # prediction times for space.annual methods, treat time as independent timeslices
  if ( !exists("prediction.ts", p)) p$prediction.ts = 1


  message(" ")
  message( "||| Initializing temporary storage of data and output files... ")
  message( "||| These are large files (4 to 6 X 5GB), it will take a minute ... ")
  message( "||| Try to turn off swap/paging such that only RAM is used. ")

  stmv_db( p=p, DS="cleanup" )


  # construct prediction/output grid area ('pa')
  p$windowsize.half = floor(p$stmv_distance_prediction/p$pres) # convert distance to discretized increments of row/col indices; stmv_distance_prediction = 0.75* stmv_distance_statsgrid (unless overridden)

  if (exists("stmv_Y_transform", p)) {
    DATA$input[, p$variables$Y ] = p$stmv_Y_transform$transf(DATA$input[, p$variables$Y ] )
  }

  if ( any(grepl("globalmodel", runmode) ) ) {
    if ( exists("stmv_global_modelengine", p) ) {
      if ( p$stmv_global_modelengine !="none" ) {
        if ( p$stmv_global_modelformula !="none" ) {
        # to add global covariate model (hierarchical)
        # .. simplistic this way but faster ~ kriging with external drift
          stmv_db( p=p, DS="global_model.redo", B=DATA$input )
          if (any(grepl("globalmodel.only", runmode))) {
            o = stmv_db( p=p, DS="global_model" )
            return( o )
          }
        }
      }
    }
  }

  # NOTE:: must not sink the following memory allocation steps into a deeper function as
  # NOTE:: bigmemory RAM seems to lose the pointers if they are not made simultaneously ?

  # init output data objects
  # statistics storage matrix ( aggregation window, coords ) .. no inputs required
  sbox = Sloc = NULL
    sbox = list(
      plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$stmv_distance_statsgrid ),
      plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$stmv_distance_statsgrid ) )
    # statistics coordinates
    Sloc = as.matrix( expand.grid( sbox$plons, sbox$plats ))
    nSloc = nrow(Sloc)
      if (p$storage.backend == "bigmemory.ram" ) {
        tmp_Sloc = big.matrix(nrow=nSloc, ncol=ncol(Sloc), type="double"  )
        tmp_Sloc[] = Sloc[]
        p$ptr$Sloc  = bigmemory::describe( tmp_Sloc  )
      }
      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$ptr$Sloc  = p$cache$Sloc
        bigmemory::as.big.matrix( Sloc, type="double", backingfile=basename(p$bm$Sloc), descriptorfile=basename(p$cache$Sloc), backingpath=p$stmvSaveDir )
      }
      if (p$storage.backend == "ff" ) {
        p$ptr$Sloc = ff( Sloc, dim=dim(Sloc), file=p$cache$Sloc, overwrite=TRUE )
      }
  sbox = Sloc = NULL


  sS = NULL
  if (!is.null(use_saved_state)) {
    if (use_saved_state=="ram") {
      # nothing needs to be done as pointers are already set up and pointed to the data
    }
    if (use_saved_state=="disk") {
      if (file.exists(p$saved_state_fn$stats)) load( p$saved_state_fn$stats )
      if (is.vector(sS)) sS = as.matrix(sS, nrow=nSloc, ncol=1)
    }
  } else {
    sS = matrix( NaN, nrow=nSloc, ncol=length( p$statsvars ) ) # NA forces into logical
  }
  if (!is.null(sS)) {
    if (p$storage.backend == "bigmemory.ram" ) {
      tmp_S = big.matrix(nrow=nSloc, ncol=length( p$statsvars ), type="double"  )
      tmp_S[] = sS[]
      p$ptr$S  = bigmemory::describe( tmp_S )
    }
    if (p$storage.backend == "bigmemory.filebacked" ) {
      p$ptr$S  = p$cache$S
      bigmemory::as.big.matrix( sS, type="double", backingfile=basename(p$bm$S), descriptorfile=basename(p$cache$S), backingpath=p$stmvSaveDir )
    }
    if (p$storage.backend == "ff" ) {
      p$ptr$S = ff( sS, dim=dim(sS), file=p$cache$S, overwrite=TRUE )
    }
    sS = NULL
  }


  sSflag = NULL
  if (!is.null(use_saved_state)) {
    if (use_saved_state=="ram") {
      # nothing needs to be done as pointers are already set up and pointed to the data
    }
    if (use_saved_state=="disk") {
      if (file.exists(p$saved_state_fn$sflag)) load( p$saved_state_fn$sflag )
      if (is.vector(sSflag)) sSflag = as.matrix(sSflag, nrow=nSloc, ncol=1)
    }
  } else {
    sSflag = matrix( stmv_error_codes()[["todo"]], nrow=nSloc, ncol=1 )
  }
  if (!is.null(sSflag)) {
    if (p$storage.backend == "bigmemory.ram" ) {
      tmp_Sflag = big.matrix(nrow=length(sSflag), ncol=1, type="double" )
      tmp_Sflag[] = sSflag[]
      p$ptr$Sflag  = bigmemory::describe( tmp_Sflag )
    }
    if (p$storage.backend == "bigmemory.filebacked" ) {
      p$ptr$Sflag  = p$cache$Sflag
      bigmemory::as.big.matrix( sSflag, type="double", backingfile=basename(p$bm$Sflag), descriptorfile=basename(p$cache$Sflag), backingpath=p$stmvSaveDir )
    }
    if (p$storage.backend == "ff" ) {
      p$ptr$Sflag = ff( sSflag, dim=dim(sSflag), file=p$cache$Sflag, overwrite=TRUE )
    }
    sSflag = NULL
  }




  # data to be worked upon .. either the raw data or covariate-residuals
  Ydata = as.matrix(DATA$input[, p$variables$Y ])
  if (exists("stmv_global_modelengine", p)) {
    if (p$stmv_global_modelengine !="none" ) {
      if (p$stmv_global_modelengine == "userdefined") {
        covmodel = stmv_db( p=p, DS="global_model")
        # TODO MUST find a generic form as below
        # # Ypreds = predict(covmodel, type="link", se.fit=FALSE )  ## TODO .. keep track of the SE
        if (!exists("predict", p$stmv_global_model)) {
          message( "p$stmv_global_model$predict =' " )
          message( "   predict( global_model, newdata=pa, type='link', se.fit=TRUE )' " )
          message( " where 'global_model', newdata=pa' are required " )
          stop()
        }
        preds = try( eval(parse(text=pp$stmv_global_model$predict )) )
        preds = p$stmv_global_model$predict( covmodel )  # needs to be tested. .. JC
        Ydata  = preds - Ydata # ie. i`nternal (link) scale
        Yq = quantile( Ydata, probs=p$stmv_quantile_bounds )
        Ydata[ Ydata < Yq[1] ] = Yq[1]
        Ydata[ Ydata > Yq[2] ] = Yq[2]
        covmodel =NULL

      } else {
        # at present only those that have a predict and residuals methods ...
        covmodel = stmv_db( p=p, DS="global_model")
        # Ypreds = predict(covmodel, type="link", se.fit=FALSE )  ## TODO .. keep track of the SE

        Ydata  = residuals(covmodel, type="working") # ie. internal (link) scale
        Yq = quantile( Ydata, probs=p$stmv_quantile_bounds )
        Ydata[ Ydata < Yq[1] ] = Yq[1]
        Ydata[ Ydata > Yq[2] ] = Yq[2]
        covmodel =NULL

      }
    }
  }
  # Ypreds = NULL
  Ydata = as.matrix( Ydata )
    if (p$storage.backend == "bigmemory.ram" ) {
      tmp_Y = big.matrix( nrow=nrow(Ydata), ncol=1, type="double"  )
      tmp_Y[] = Ydata[]
      p$ptr$Y  = bigmemory::describe( tmp_Y )
    }
    if (p$storage.backend == "bigmemory.filebacked" ) {
      p$ptr$Y  = p$cache$Y
      bigmemory::as.big.matrix( Ydata, type="double", backingfile=basename(p$bm$Y), descriptorfile=basename(p$cache$Y), backingpath=p$stmvSaveDir )
    }
    if (p$storage.backend == "ff" ) {
      p$ptr$Y = ff( Ydata, dim=dim(Ydata), file=p$cache$Y, overwrite=TRUE )
    }
  Ydata = NULL

  Y = stmv_attach( p$storage.backend, p$ptr$Y )

 # data coordinates
  Yloc = as.matrix( DATA$input[, p$variables$LOCS ])
    if (p$storage.backend == "bigmemory.ram" ) {
      tmp_Yloc = big.matrix( nrow=nrow(Yloc), ncol=ncol(Yloc), type="double" )
      tmp_Yloc[] = Yloc[]
      p$ptr$Yloc = bigmemory::describe( tmp_Yloc )
    }
    if (p$storage.backend == "bigmemory.filebacked" ) {
      p$ptr$Yloc  = p$cache$Yloc
      bigmemory::as.big.matrix( Yloc, type="double", backingfile=basename(p$bm$Yloc), descriptorfile=basename(p$cache$Yloc), backingpath=p$stmvSaveDir )
    }
    if (p$storage.backend == "ff" ) {
      p$ptr$Yloc = ff( Yloc, dim=dim(Yloc), file=p$cache$Yloc, overwrite=TRUE )
    }
  Yloc = NULL

    # independent variables/ covariate
    if (exists("COV", p$variables)) {
      if (length(p$variables$COV) > 0) {
        Ycov = as.matrix(  DATA$input[ , p$variables$COV ] )
        if (p$storage.backend == "bigmemory.ram" ) {
          tmp_Ycov = big.matrix( nrow=nrow(Ycov), ncol=ncol(Ycov), type="double")
          tmp_Ycov[] = Ycov[]
          p$ptr$Ycov  = bigmemory::describe( tmp_Ycov )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Ycov  = p$cache$Ycov
          bigmemory::as.big.matrix( Ycov, type="double", backingfile=basename(p$bm$Ycov), descriptorfile=basename(p$cache$Ycov), backingpath=p$stmvSaveDir )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Ycov = ff( Ycov, dim=dim(Ycov), file=p$cache$Ycov, overwrite=TRUE )
        }
        Ycov= NULL
      }
    }

    # data times
    if ( exists("TIME", p$variables) ) {
      Ytime = as.matrix(  DATA$input[, p$variables$TIME ] )
        if (p$storage.backend == "bigmemory.ram" ) {
          tmp_Ytime = big.matrix( nrow=nrow(Ytime), ncol=ncol(Ytime), type="double"  )
          tmp_Ytime[] = Ytime[]
          p$ptr$Ytime  = bigmemory::describe( tmp_Ytime )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Ytime  = p$cache$Ytime
          bigmemory::as.big.matrix( Ytime, type="double", backingfile=basename(p$bm$Ytime), descriptorfile=basename(p$cache$Ytime), backingpath=p$stmvSaveDir )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Ytime = ff( Ytime, dim=dim(Ytime), file=p$cache$Ytime, overwrite=TRUE )
        }
      Ytime =NULL;

    }


    nPlocs = nrow(DATA$output$LOCS)

    if (exists("COV", p$variables)) {
      if (length(p$variables$COV) > 0) {
        # this needs to be done as Prediction covars need to be structured as lists
        p$ptr$Pcov = list()
        tmp_Pcov = list()
        for ( covname in p$variables$COV ) {
          Pcovdata = as.matrix( DATA$output$COV[[covname]] )
          nPcovlocs = nrow(Pcovdata)
          if (nPcovlocs != nPlocs) {
            message( "Inconsistency between number of prediction locations and prediction covariates: input data needs to be checked:")
            message( "Usually this due to bathymetry and temperature being out of sync")

            print(str(DATA))
            stop()
          }
          attr( Pcovdata, "dimnames" ) = NULL
          if (p$storage.backend == "bigmemory.ram" ) {
            tmp_Pcov[[covname]] = big.matrix( nrow=nPcovlocs, ncol=ncol(Pcovdata), type="double"  )
            tmp_Pcov[[covname]][] = Pcovdata[]
            p$ptr$Pcov[[covname]]  = bigmemory::describe( tmp_Pcov[[covname]] )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$Pcov[[covname]]  = p$cache$Pcov[[covname]]
            bigmemory::as.big.matrix( Pcovdata, type="double", backingfile=basename(p$bm$Pcov[[covname]]), descriptorfile=basename(p$cache$Pcov[[covname]]), backingpath=p$stmvSaveDir )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Pcov[[covname]] = ff( Pcovdata, dim=dim(Pcovdata), file=p$cache$Pcov[[covname]], overwrite=TRUE )
          }
          Pcovdata = NULL
        }
      }
    }



  # predictions and associated stats
  sP = NULL
  if (!is.null(use_saved_state)) {
    if (use_saved_state=="ram") {
      # nothing needs to be done as pointers are already set up and pointed to the data
    }
    if (use_saved_state=="disk") {
      if (file.exists(p$saved_state_fn$P)) load( p$saved_state_fn$P )
      if (is.vector(sP)) sP=as.matrix(sP, nrow=nPlocs, ncol=1)
    }
  } else {
     sP = matrix( NaN, nrow=nPlocs, ncol=p$nt )
  }
  if (!is.null(sP)) {
      if (p$storage.backend == "bigmemory.ram" ) {
        tmp_P = big.matrix( nrow=nrow(sP), ncol=ncol(sP), type="double" )
        tmp_P[] = sP[]
        p$ptr$P  = bigmemory::describe( tmp_P )
      }
      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$ptr$P  = p$cache$P
        bigmemory::as.big.matrix( sP, type="double", backingfile=basename(p$bm$P), descriptorfile=basename(p$cache$P), backingpath=p$stmvSaveDir )
      }
      if (p$storage.backend == "ff" ) {
        p$ptr$P = ff( sP, dim=dim(sP), file=p$cache$P, overwrite=TRUE )
      }
    sP = NULL
  }


  # count of prediction estimates
  sPn = NULL
  if (!is.null(use_saved_state)) {
    if (use_saved_state=="ram") {
      # nothing needs to be done as pointers are already set up and pointed to the data
    }
    if (use_saved_state=="disk") {
      if (file.exists(p$saved_state_fn$Pn)) load( p$saved_state_fn$Pn )
      if (is.vector(sPn)) sPn = as.matrix(sPn, nrow=nPlocs, ncol=1)
    }
  } else {
    sPn = matrix( NaN, nrow=nPlocs, ncol=p$nt )
  }
  if (!is.null(sPn)) {
   if (p$storage.backend == "bigmemory.ram" ) {
      tmp_Pn = big.matrix( nrow=nrow(sPn), ncol=ncol(sPn), type="double" )
      tmp_Pn[] = sPn[]
      p$ptr$Pn = bigmemory::describe( tmp_Pn )
    }
    if (p$storage.backend == "bigmemory.filebacked" ) {
      p$ptr$Pn  = p$cache$Pn
      bigmemory::as.big.matrix( sPn, type="double", backingfile=basename(p$bm$Pn), descriptorfile=basename(p$cache$Pn), backingpath=p$stmvSaveDir )
    }
    if (p$storage.backend == "ff" ) {
      p$ptr$Pn = ff( sPn, dim=dim(sPn), file=p$cache$Pn, overwrite=TRUE )
    }
    sPn = NULL
  }



  # sd of prediction estimates
  sPsd = NULL
  if (!is.null(use_saved_state)) {
    if (use_saved_state=="ram") {
      # nothing needs to be done as pointers are already set up and pointed to the data
    }
    if (use_saved_state=="disk") {
      if (file.exists(p$saved_state_fn$Psd)) load( p$saved_state_fn$Psd )
      if (is.vector(sPsd)) sPsd = as.matrix(sPsd, nrow=nPlocs, ncol=1)
    }
  } else {
    sPsd = matrix( NaN, nrow=nPlocs, ncol=p$nt )
  }
  if (!is.null(sPsd)) {
    if (p$storage.backend == "bigmemory.ram" ) {
      tmp_Psd = big.matrix( nrow=nrow(sPsd), ncol=ncol(sPsd), type="double" )
      tmp_Psd[] = sPsd[]
      p$ptr$Psd =bigmemory::describe( tmp_Psd )
    }
    if (p$storage.backend == "bigmemory.filebacked" ) {
      p$ptr$Psd  = p$cache$Psd
      bigmemory::as.big.matrix( sPsd, type="double", backingfile=basename(p$bm$Psd), descriptorfile=basename(p$cache$Psd), backingpath=p$stmvSaveDir )
    }
    if (p$storage.backend == "ff" ) {
      p$ptr$Psd = ff( sPsd, dim=dim(sPsd), file=p$cache$Psd, overwrite=TRUE )
    }
    sPsd = NULL
  }

    # prediction coordinates
    Ploc = as.matrix( DATA$output$LOCS )
      attr( Ploc, "dimnames" ) = NULL
      if (p$storage.backend == "bigmemory.ram" ) {
        tmp_Ploc = big.matrix( nrow=nrow(Ploc), ncol=ncol(Ploc), type="double" )
        tmp_Ploc[] = Ploc[]
        p$ptr$Ploc  = bigmemory::describe( tmp_Ploc )
      }
      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$ptr$Ploc  = p$cache$Ploc
        bigmemory::as.big.matrix( Ploc, type="double", backingfile=basename(p$bm$Ploc), descriptorfile=basename(p$cache$Ploc), backingpath=p$stmvSaveDir )
      }
      if (p$storage.backend == "ff" ) {
        p$ptr$Ploc = ff( Ploc, dim=dim(Ploc), file=p$cache$Ploc, overwrite=TRUE )
      }
    Ploc = NULL;

    DATA = NULL;

    if (exists("stmv_global_modelengine", p) ) {
      if (p$stmv_global_modelengine !="none" ) {
        # create prediction suface .. additive offsets

        sP0 = NULL
        if (!is.null(use_saved_state)) {
          if (use_saved_state=="ram") {
            # nothing needs to be done as pointers are already set up and pointed to the data
          }
          if (use_saved_state=="disk") {
            if (file.exists(p$saved_state_fn$P0)) load( p$saved_state_fn$P0 )
            if (is.vector(sP0)) sP0 = as.matrix(sP0, nrow=nPlocs, ncol=1)
          }
        } else {
          sP0 = matrix( NaN, nrow=nPlocs, ncol=p$nt )
        }
        if (!is.null(sP0)) {
          if (p$storage.backend == "bigmemory.ram" ) {
            tmp_P0= big.matrix( nrow=nrow(sP0), ncol=ncol(sP0) , type="double" )
            tmp_P0[] = sP0[]
            p$ptr$P0 = bigmemory::describe(tmp_P0 )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$P0  = p$cache$P0
            bigmemory::as.big.matrix( sP0, type="double", backingfile=basename(p$bm$P0), descriptorfile=basename(p$cache$P0), backingpath=p$stmvSaveDir )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$P0 = ff( sP0, dim=dim(sP0), file=p$cache$P0, overwrite=TRUE )
          }
          sP0 = NULL
        }


        sP0sd = NULL
        if (!is.null(use_saved_state)) {
          if (use_saved_state=="ram") {
            # nothing needs to be done as pointers are already set up and pointed to the data
          }
          if (use_saved_state=="disk") {
            if (file.exists(p$saved_state_fn$P0sd)) load( p$saved_state_fn$P0sd )
            if (is.vector(sP0sd)) sP0sd = as.matrix(sP0sd, nrow=nPlocs, ncol=1)
          }
        } else {
          sP0sd = matrix( NaN, nrow=nPlocs, ncol=p$nt )
        }
        if (!is.null(sP0sd)) {
          if (p$storage.backend == "bigmemory.ram" ) {
            tmp_P0sd= big.matrix( nrow=nrow(sP0sd), ncol=ncol(sP0sd) , type="double" )
            tmp_P0sd[] = sP0sd[]
            p$ptr$P0sd = bigmemory::describe(tmp_P0sd )
          }
          if (p$storage.backend == "bigmemory.filebacked" ) {
            p$ptr$P0sd  = p$cache$P0sd
            bigmemory::as.big.matrix( sP0sd, type="double", backingfile=basename(p$bm$P0sd), descriptorfile=basename(p$cache$P0sd), backingpath=p$stmvSaveDir )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$P0sd = ff( sP0sd, dim=dim(sP0sd), file=p$cache$P0sd, overwrite=TRUE )
          }
          sP0sd=NULL;
        }


        # test to see if all covars are static as this can speed up the initial predictions
        message(" ")
        message( "||| Predicting global effect of covariates at each prediction location ... ")
        message( "||| depending upon the size of the prediction grid and number of cpus (~1hr?).. ")

        p$time_covariates_0 =  Sys.time()
        if (exists("COV", p$variables)) {
          if (length(p$variables$COV) > 0) {
            nc_cov =NULL
            for (i in p$variables$COV ) {
              pu = stmv_attach( p$storage.backend, p$ptr$Pcov[[i]] )
              nc_cov = c( nc_cov,  ncol(pu) )
            }
            p$all.covars.static = ifelse( any(nc_cov > 1),  FALSE, TRUE )
          } else {
            p$all.covars.static = TRUE # degenerate case where the model is an intercept-only model (to remove mean effects)
          }
        }

        if (!is.null(use_saved_state)) {
            # nothing needs to be done as pointers are already set up and pointed to the data
        } else {
          pc = p # temp copy
          if (exists("all.covars.static", pc)) if (!pc$all.covars.static) if (exists("clusters.covars", pc) ) pc$clusters = pc$clusters.covars
          # takes up to about 28 GB per run (in temperature).. adjust cluster number temporarily
          # pc = parallel_run( p=pc, runindex=list( it=1:p$nt ) )
          # stmv_predict_globalmodel(p=pc)

          global_model = stmv_db( p=p, DS="global_model")
          if (is.null(global_model)) stop("Global model not found.")

          YYY = predict( global_model, type="link", se.fit=TRUE )  # determine bounds from data
          Yq = quantile( YYY$fit, probs=p$stmv_quantile_bounds )
          YYY = NULL

          parallel_run( FUNC=stmv_predict_globalmodel, p=pc, global_model=global_model, Yq=Yq, runindex=list( pnt=1:p$nt ) )

          p$time_covariates = round(difftime( Sys.time(), p$time_covariates_0 , units="hours"), 3)
          message( paste( "||| Time taken to predict covariate surface (hours):", p$time_covariates ) )
        }
      }
    }

    if (is.null(use_saved_state)) {
      stmv_db( p=p, DS="statistics.Sflag" )
    }

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
      if (length(p$variables$COV) > 0) {
        Ycov = stmv_attach( p$storage.backend, p$ptr$Ycov )
        if (length(p$variables$COV)==1) {
          bad = which( !is.finite( Ycov[] ))
        } else {
          bad = which( !is.finite( rowSums(Ycov[])))
        }
        if (length(bad)> 0 ) Yi[bad] = NA
        Yi = na.omit(Yi)
      }
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
        tmp_Yi[] = Yi[]
        p$ptr$Yi  = bigmemory::describe( tmp_Yi )
      }
      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$ptr$Yi  = p$cache$Yi
        bigmemory::as.big.matrix( Yi, type="double", backingfile=basename(p$bm$Yi), descriptorfile=basename(p$cache$Yi), backingpath=p$stmvSaveDir )
      }
      if (p$storage.backend == "ff" ) {
        p$ptr$Yi = ff( Yi, dim=dim(Yi), file=p$cache$Yi, overwrite=TRUE )
      }
    Yi=NULL

    if ( !exists("stmv_distance_scale", p)) {
      Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )
      p$stmv_distance_scale = min( diff(range( Yloc[,1]) ), diff(range( Yloc[,2]) ) ) / 10
      message( paste( "||| Crude distance scale:", p$stmv_distance_scale, "" ) )
    }
    if ( !exists("stmv_distance_min", p)) p$stmv_distance_min = mean( c(p$stmv_distance_prediction, p$stmv_distance_scale /20 ) )
    if ( !exists("stmv_distance_max", p)) p$stmv_distance_max = mean( c(p$stmv_distance_prediction*10, p$stmv_distance_scale * 2 ) )

    p$stmv_distance_scale0 = p$stmv_distance_scale  # store copy as this can get modified below


    if (p$stmv_local_modelengine == "twostep") {
      if (exists("stmv_rsquared_threshold", p) ) {
        if (p$stmv_rsquared_threshold > 0) {
          message( "Ignoring value of p$stmv_rsquared_threshold as it is meaningless with twostep")
          p$stmv_rsquared_threshold = 0  # override :: this is meaningless when broken apart in space and time ..
        }
      }
    }

  p$upsampling = c(1.0, 1.25, 1.5, 1.75, 2.0) * p$stmv_distance_scale
  p$upsampling = p$upsampling[ which(p$upsampling <= p$stmv_distance_max )]

  if ( exists("TIME", p$variables)) {
    p$minresolution = p$downsampling_multiplier*c(p$pres, p$pres, p$tres)
  } else {
    p$minresolution = p$downsampling_multiplier*c(p$pres, p$pres )
  }


    message("||| Finished preparing data structures ... ")
    # message("||| Once analyses begin, you can view maps from an external R session (e.g. for temperature): ")
    # message("||| p = stmv_db( p=list(data_root=project.datadirectory('aegis', 'temperature'), variables=list(Y='t'), spatial.domain='canada.east' ), DS='load.parameters' )" )
    # message("||| see stmv(p=p, runmode='debug_predictions_map', debug_plot_variable_index=1) # for static maps")
    # message("||| see stmv(p=p, runmode='debug_predictions_map', debug_plot_variable_index=1:p$nt, debug_plot_log=TRUE) # for timeseries  of log(Y)")
    # message("||| see stmv(p=p, runmode='debug_statistics_map', debug_plot_variable_index=1:length(p$statsvars))  ")
    # message("||| print( p$statsvars) # will get you your stats variables " )
    message("||| Monitor the status of modelling by looking at the output of the following file:")
    message("||| in linux, you can issue the following command:" )
    message("||| watch -n 60 cat ",  p$stmv_current_status  )


    # end of intialization of data structures
    # -----------------------------------------------------

  p <<- p  # copy to parent (calling) environment

  if ( "initialize_only" %in% runmode ) return(p)

  file.create( p$stmv_current_status, showWarnings=FALSE )
  currentstatus = stmv_logfile(p = p)  # init log file

  if ( any(grepl("debug", runmode)) ) {

    message( " " )
    message( "||| Debugging from man stmv call." )
    message( "||| To load from the saved state try: stmv_db( p=p, DS='load_saved_state' ) " )
    message( "||| To reset stats: currentstatus = stmv_db( p=p, DS='statistics.status.reset' ) " )

    # -----------------------------------------------------
    if ( "debug" %in% runmode ) {
      currentstatus = stmv_db( p=p, DS="statistics.status" )
      pdeb = parallel_run( p=p, runindex=list( locs=sample( currentstatus$todo )) ) # reconstruct reauired params
      print( c( unlist( currentstatus[ c("n.total", "n.too_shallow", "n.todo", "n.skipped", "n.outside_bounds", "n.complete" ) ] ) ))
      message( "||| Entering browser mode ...")
      p <<- pdeb
      browser()
      debug(stmv_interpolate)
      pdeb$time_start_interpolation_debug = Sys.time()
      stmv_interpolate (p=pdeb, debugging=TRUE )
    }

    # -----------------------------------------------------

    if ( "debug_predictions_map" %in% runmode) {
      if ( is.null(debug_plot_variable_index)) debug_plot_variable_index=1
      if (debug_data_source=="saved.state" ) {
        load( p$saved_state_fn$P )
        if (exists("stmv_global_modelengine", p)) {
          if (p$stmv_global_modelengine !="none" ) {
            load( p$saved_state_fn$P0 )
            sP[] = sP[] + sP0[]
          }
        }
      } else {
        P = stmv_attach( p$storage.backend, p$ptr$P )
        sP = P[]
        if (exists("stmv_global_modelengine", p)) {
          if (p$stmv_global_modelengine !="none" ) {
            P0 = stmv_attach( p$storage.backend, p$ptr$P0 )
            sP = sP[] + P0[]
          }
        }
      }
      sP = p$stmv_global_family$linkinv( sP[] )
      if (exists("stmv_Y_transform", p)) {
        sP = p$stmv_Y_transform$invers (sP[])
      }
      if ( debug_plot_log ) sP = log(sP)
      Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
      for (i in debug_plot_variable_index ) {
        print(
          lattice::levelplot( sP[,i] ~ Ploc[,1]+Ploc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso")
        )
      }
    }

    # -----------------------------------------------------

    if ( "debug_statistics_map" %in% runmode) {
      if (debug_data_source=="saved.state" ) {
        load( p$saved_state_fn$stats )
      } else {
        sS = stmv_attach( p$storage.backend, p$ptr$S )[]
      }
      Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
      if ( debug_plot_log ) sS = log(sS)
      if ( is.null(debug_plot_variable_index)) debug_plot_variable_index=1
      for (i in debug_plot_variable_index ) {
        print(
          lattice::levelplot( sS[,i] ~ Sloc[,1]+Sloc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso")
        )
      }

    }

  }

  # -----------------------------------------------------

  if ( "reset_incomplete_locations" %in% runmode ) {
    # this resets errors flags and areas without viable predictions
    toredo = stmv_db( p=p, DS="flag.incomplete.predictions" )
    if ( !is.null(toredo) && length(toredo) > 0) {
      Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )
      Sflag[toredo] = stmv_error_codes()[["todo"]]
    }
  }

  stmv_db( p=p, DS="save.parameters" )  # save in case a restart is required .. mostly for the pointers to data


  # -----------------------------------------------------

  if ("interpolate" %in% runmode ) {
    p$clusters0 = p$clusters
    sm = sort( unique( c(1, p$sampling[ p$sampling > 1] ) ) )
    print ( "Sampling at the following distance mulitpliers:" )
    print (sm)
    for ( smult in sm) {
      print( paste("Entering interpolation at distance multiplier:", smult ) )
      p$stmv_distance_scale = p$stmv_distance_scale0 * smult
      # p$clusters = p$clusters[-1] # as ram reqeuirements increase drop cpus
      toredo = stmv_db( p=p, DS="flag.incomplete.predictions" )
      if ( !is.null(toredo) && length(toredo) > 0) {
        Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )
        Sflag[toredo] = stmv_error_codes()[["todo"]]
      }
      currentstatus = stmv_db( p=p, DS="statistics.status" )
      p$time_start_interpolation = Sys.time()
      parallel_run( stmv_interpolate, p=p, runindex=list( locs=sample( currentstatus$todo )))
    }
  }

  # --------------------
  if (!is.null(use_saved_state)) {
    if (use_saved_state=="disk") {
      # a penultimate save of data as an internal format, just in case the save step or force complete goes funny
      sP = stmv_attach( p$storage.backend, p$ptr$P )[]
      sPn = stmv_attach( p$storage.backend, p$ptr$Pn )[]
      sPsd = stmv_attach( p$storage.backend, p$ptr$Psd )[]
      sS = stmv_attach( p$storage.backend, p$ptr$S )[]
      sSflag = stmv_attach( p$storage.backend, p$ptr$Sflag )[]
      if (exists("stmv_global_modelengine", p)) {
        if (p$stmv_global_modelengine !="none" ) {
          sP0 = stmv_attach( p$storage.backend, p$ptr$P0 )[]
          sP0sd = stmv_attach( p$storage.backend, p$ptr$P0sd )[]
        }
      }
      save( sP, file=p$saved_state_fn$P, compress=TRUE ); sP = NULL
      save( sPn, file=p$saved_state_fn$Pn, compress=TRUE ); sPn = NULL
      save( sPsd, file=p$saved_state_fn$Psd, compress=TRUE ); sPsd=NULL
      save( sS, file=p$saved_state_fn$stats, compress=TRUE ); sS = NULL
      save( sSflag, file=p$saved_state_fn$sflag, compress=TRUE ); sSflag = NULL

      if (exists("stmv_global_modelengine", p)) {
        if (p$stmv_global_modelengine !="none" ) {
          save( sP0,   file=p$saved_state_fn$P0,   compress=TRUE ); sP0 = NULL
          save( sP0sd, file=p$saved_state_fn$P0sd, compress=TRUE ); sP0sd = NULL
        }
      }
    }
  }

  # --------------------

  if (force_complete_solution) {
    print( "Entering -force complete solution- interpolation stage" )
    # finalize all interpolations where there are missing data/predictions using
    # interpolation based on data and augmented by previous predictions
    # NOTE:: no covariates are used
    p$force_complete_solution = TRUE
    p$clusters = p$clusters0
    p$stmv_distance_scale = p$stmv_distance_scale0
    toredo = stmv_db( p=p, DS="flag.incomplete.predictions" )
    if ( !is.null(toredo) && length(toredo) > 0) {
      Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )
      Sflag[toredo] = stmv_error_codes()[["todo"]]
    }
    currentstatus = stmv_db( p=p, DS="statistics.status" )
    p$stmv_local_modelengine = "fft"
    if (!exists("stmv_fft_filter", p)) p$stmv_fft_filter="spatial.process"  #  fft==spatial.process, krige (very slow), lowpass, lowpass_spatial.process
    if (p$stmv_fft_filter=="lowpass" ) {
      if (!exists("stmv_lowpass_phi", p))  p$stmv_lowpass_phi = p$pres / 5 # FFT-baed methods cov range parameter .. not required for "spatial.process" ..
      if (!exists("stmv_lowpass_nu", p))  p$stmv_lowpass_nu = 0.5  #exponential
    }
    p$time_start_interpolation_force_complete = Sys.time()
    parallel_run( stmv_interpolate, p=p, runindex=list( locs=sample( currentstatus$todo )))
  }


  # -----------------------------------------------------

  if (save_completed_data) stmv_db( p=p, DS="stmv.results" ) # save to disk for use outside stmv*, returning to user scale

  # -----------------------------------------------------

  if ( "cleanup" %in% runmode ) {

    if ( p$storage.backend !="bigmemory.ram" ) {
      resp = readline( "||| Delete temporary files? Type to confirm <YES>:  ")
      if (resp=="YES") {
        stmv_db( p=p, DS="cleanup.all" )
      } else {
        message(" ")
        message( "||| Leaving temporary files alone in case you need to examine them or restart a process. ")
        message( "||| You can delete them by running: stmv_db( p=p, DS='cleanup.all' ), once you are done. ")
      }
    }
  }

  p$time_total = round( difftime( Sys.time(), p$time_start, units="hours" ), 3)
  message( paste( "||| Total time taken (hours):", p$time_total ) )
  message( paste( "||| NOTE:: parameter 'p' has been updated in case you need to re-run something" ) )

}
