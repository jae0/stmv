

stmv = function( p, runmode=NULL, DATA=NULL, nlogs=200, niter=1,
  debug_plot_variable_index=1, debug_data_source="saved.state", debug_plot_log=FALSE, robustify_quantiles=c(0.005, 0.995), ... ) {

  if (0) {
    runmode = NULL
    nlogs = 25
    niter = 1
    DATA=NULL
    debug_plot_variable_index=1
    robustify_quantiles=c(0.005, 0.995)
    # runmode=c("interpolate", "globalmodel")
    # runmode=c("interpolate")
  }

  #\\ localized modelling of space and time data to predict/interpolate upon a grid
  #\\ speed ratings: bigmemory.ram (1), ff (2), bigmemory.filebacked (3)

  # -----------------------------------------------------

  if (!exists("time_start", p) ) p$time_start = Sys.time()

  # if runmode is passed, it overrides everything ..
  if (is.null(runmode)) {
    s_runmode = names(p$stmv_runmode)
    o = rep(TRUE, length(p$stmv_runmode))
    TF = c( sapply(p$stmv_runmode, is.logical ))
    o[TF] = p$stmv_runmode[TF]
    s_runmode = s_runmode[ unlist(o) ]
    runmode = intersect(s_runmode, c( "globalmodel", "scale", "interpolate", "interpolate_force_complete", "save_completed_data", "save_intermediate_results", "restart_load") )
  }
  message( "Runmodes: ", paste(runmode, ", ", sep=" ")  )

  p$nlogs = nlogs

  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add )
  p = stmv_parameters(p=p ) # fill in parameters with defaults where required
  p = stmv_db( p=p, DS="filenames" )
  p$ptr = list() # location for data pointers

  # set up the data and problem using data objects
  tmpfiles = unlist( p$cache)
  for (tf in tmpfiles) if (file.exists( tf)) file.remove(tf)

  message( "||| Initializing data files ... " )
  message( "||| In case something should go wrong, intermediary outputs will be placed at:" )
  message( "||| ",  p$stmvSaveDir  )

  message( "||| Monitor status by looking at the output of the following file:")
  message( "||| in linux, you can issue the following command: \n" )
  message( "||| watch -n 60 cat ",  p$stmv_current_status  )



  # permit passing a function rather than data directly .. less RAM usage in parent call
  if (is.null(DATA) ) {
    if (exists("DATA", p)) {
      if (class(p$DATA)=="character") {
        assign("DATA", eval(parse(text=p$DATA) ) )
      } else {
        DATA = p$DATA
        p$DATA = NULL  # to reduce RAM requirements
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

  message( "\n")
  message( "||| Initializing temporary storage of data and output files... ")
  message( "||| These are large files (4 to 6 X 5GB), it will take a minute ... ")
  message( "||| Try to turn off swap/paging such that only RAM is used. ")

  stmv_db( p=p, DS="cleanup" )

  # if (!extrapolate_predictions) {
  #   # force output covars to be in the same range as input data .. no extrapolation
  #   for ( vn in p$variables$COV ){
  #     vnrange = range( DATA$input[,vn], na.rm=TRUE )
  #     toolow = which( DATA$output$COV[[vn]] < vnrange[1] )
  #     if (length(toolow) > 1) DATA$output$COV[[vn]][toolow] = vnrange[1]
  #     toohigh = which( DATA$output$COV[[vn]] > vnrange[1] )
  #     if (length(toohigh) > 1) DATA$output$COV[[vn]][toohigh] = vnrange[2]
  #   }
  # }


  # Ypreds = NULL
  Ydata = as.matrix( DATA$input[, p$variables$Y ] )
  if (exists("stmv_Y_transform", p)) Ydata = p$stmv_Y_transform$transf( Ydata )
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

  if ( any(grepl("globalmodel", runmode) ) ) {
    if ( exists("stmv_global_modelengine", p) ) {
      if ( p$stmv_global_modelengine !="none" ) {
        if ( p$stmv_global_modelformula !="none" ) {
          stmv_db( p=p, DS="global_model.redo", B=DATA$input )
        }
        if ( any(grepl("globalmodel.only", runmode)))  return( stmv_db( p=p, DS="global_model" ) )
        # model complete .. now predict to get residuals
        if ( p$stmv_global_modelengine == "userdefined") {
          global_model = stmv_db( p=p, DS="global_model")
          # TODO MUST find a generic form as below
          # # Ypreds = predict(global_model, type="link", se.fit=FALSE )  ## TODO .. keep track of the SE
          if (!exists("predict", p$stmv_global_model)) {
            message( "||| p$stmv_global_model$predict =' " )
            message( "|||   predict( global_model, newdata=pa, type='link', se.fit=TRUE )' " )
            message( "||| where 'global_model', newdata=pa' are required " )
            stop()
          }
          preds = try( eval(parse(text=pp$stmv_global_model$predict )) )
          preds = p$stmv_global_model$predict( global_model )  # needs to be tested. .. JC
          Ydata  = Ydata - preds # ie. internalR (link) scale
        } else {
          # at present only those that have a predict and residuals methods ...
          global_model = stmv_db( p=p, DS="global_model")
          Ydata  = residuals(global_model, type="working") # ie. internal (link) scale
        }
        global_model =NULL
      }
    }
  }

  Yq = quantile( Ydata, probs=robustify_quantiles )  # extreme data can make convergence slow and erratic 99.9%CI
  Ydata[ Ydata < Yq[1] ] = Yq[1]
  Ydata[ Ydata > Yq[2] ] = Yq[2]

  Y = stmv_attach( p$storage.backend, p$ptr$Y )
  Y[] = Ydata[]  # update Ydata ...
  Ydata = NULL

  gc()





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


  ##############


  # init output data objects
  # statistics storage matrix ( aggregation window, coords ) .. no inputs required
  sbox = Sloc = NULL
  sbox = list(
    plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$stmv_distance_statsgrid ),
    plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$stmv_distance_statsgrid ) )
  # statistics coordinates
  Sloc = as.matrix( expand_grid_fast( sbox$plons, sbox$plats ))
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
  Sloc = NULL
  sbox = NULL


  sS = matrix( NaN, nrow=nSloc, ncol=length( p$statsvars ) ) # NA forces into logical
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


  sSflag = matrix( stmv_error_codes()[["todo"]], nrow=nSloc, ncol=1 )
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


  nPloc = nrow(DATA$output$LOCS)

  if (exists("COV", p$variables)) {
    if (length(p$variables$COV) > 0) {
      # this needs to be done as Prediction covars need to be structured as lists
      p$ptr$Pcov = list()
      tmp_Pcov = list()
      nc_cov =NULL
      for ( covname in p$variables$COV ) {
        Pcovdata = as.matrix( DATA$output$COV[[covname]] )
        nc_cov = c( nc_cov,  ncol(Pcovdata) )  # test no. cols
        nPcovloc = nrow(Pcovdata)
        if (nPcovloc != nPloc) {
          message( "||| Inconsistency between number of prediction locations and prediction covariates: input data needs to be checked:")
          message( "||| Usually this due to bathymetry and temperature being out of sync")

          print(str(DATA))
          stop()
        }
        attr( Pcovdata, "dimnames" ) = NULL
        if (p$storage.backend == "bigmemory.ram" ) {
          tmp_Pcov[[covname]] = big.matrix( nrow=nPcovloc, ncol=ncol(Pcovdata), type="double"  )
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
      p$all.covars.static = ifelse( any(nc_cov > 1),  FALSE, TRUE )
      nc_cov = NULL
    } else {
      p$all.covars.static = TRUE # degenerate case where the model is an intercept-only model (to remove mean effects)
    }
  }


  # predictions and associated stats
  sP = matrix( NaN, nrow=nPloc, ncol=p$nt )
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


  # count of prediction estimates
  sPn = matrix( NaN, nrow=nPloc, ncol=p$nt )
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

  # sd of prediction estimates
  sPsd = matrix( NaN, nrow=nPloc, ncol=p$nt )
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



  ################


  if (exists("stmv_global_modelengine", p) ) {
    if (p$stmv_global_modelengine !="none" ) {
      sP0 = matrix( NaN, nrow=nPloc, ncol=p$nt )
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

      sP0sd = matrix( NaN, nrow=nPloc, ncol=p$nt )
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
    }
  }


  # --------------------------------
  # completed data structures .. save params and clear up RAM

  DATA = NULL;
  if (!is.character(p$DATA)) p$DATA = NULL  # where data was sent, remove it to free RAM
  file.create( p$stmv_current_status, showWarnings=FALSE )
  p <<- p  # force copy to parent (calling) environment (to remove "p$DATA" )
  stmv_db( p=p, DS="save.parameters" )  # save in case a restart is required .. mostly for the pointers to data
  gc()


  # --------------------------------
  # analysis begins
  # check if resetimation is required

  niterations = 1  # default .. no global model nor no covariates
  if ( exists("COV", p$variables )) {
    if (length(p$variables$COV) > 0) {
      niterations = niter
    }
  }

  stmv_statistics_status( p=p, reset="features" )  # flags/filter stats locations base dupon prediction covariates. .. speed up and reduce storage

  devold = Inf
  eps = 1e-9

  iYP = stmv_index_predictions_to_observations(p)
  iYP_nomatch = which(!is.finite(iYP))

  # NOTE:: must not sink the following memory allocation steps into a deeper function as
  # NOTE:: bigmemory RAM seems to lose the pointers if they are not made simultaneously ?

  # data to be worked upon .. either the raw data or covariate-residuals
  for ( nn in 1:niterations ) {
    message( "Iteration ", nn, " of ", niterations, "\n" )

    if (exists("stmv_global_modelengine", p) ) {
      if (p$stmv_global_modelengine !="none" ) {
        # create prediction suface .. additive offsets
        # test to see if all covars are static as this can speed up the initial predictions
        message ( "\n", "||| Entering Predicting global effect of covariates at each prediction location: ", format(Sys.time()) , "\n" )
        message( "||| This can take a while (usually a few minutes but hours if complex) ... ")
        p$time_covariates_0 =  Sys.time()
        global_model = stmv_db( p=p, DS="global_model")
        if (is.null(global_model)) stop("Global model not found.")
        dev = global_model$deviance
        message("Model Deviance: ", dev)
        if (nn > 1) {
          if (abs(dev - devold)/(0.1 + abs(dev)) < eps ) break() # globalmodel_converged
          inputdata = P[ iYP ]  # P are the spatial/spatio-temporal random effects (on link scale)
          inputdata[ iYP_nomatch ] = 0  # E[RaneFF] = 0
          # return to user scale (that of Y)
          if ( exists( "stmv_global_family", p)) {
            if (p$stmv_global_family != "none") {
              if (exists("linkinv", p$stmv_global_family)) inputdata = p$stmv_global_family$linkinv( inputdata )
            }
          }
          global_model$data[, p$variables$Y ] = global_model$data[, p$variables$Y ] - inputdata  # update to current estimate of fixed effects
          inputdata = global_model$data
          global_model = NULL
          gc()
          global_model = stmv_db( p=p, DS="global_model.redo", B=inputdata,  savedata=FALSE )  # do not overwrite initial model
          inputdata = NULL
        }
        devold = dev
        global_model$data = NULL  # reduce file size/RAM
        gc()
        message("Creating fixed effects predictons")
        parallel_run( FUNC=stmv_predict_globalmodel, p=p, runindex=list( pnt=1:p$nt ), Yq=Yq, global_model=global_model )
        Ydata  = residuals(global_model, type="working") # ie. internal (link) scale

        Y = stmv_attach( p$storage.backend, p$ptr$Y )
        Y[] = Ydata[]  # update Ydata ...
        global_model = NULL
        p$time_covariates = round(difftime( Sys.time(), p$time_covariates_0 , units="hours"), 3)
        message( paste( "||| Time taken to predict covariate surface (hours):", p$time_covariates ) )
      }
    }


    # -----------------------------------------------------
    if ( "scale" %in% runmode ) {
      message ( "\n", "||| Entering spatial scale (variogram) determination: ", format(Sys.time()) , "\n" )
      if ( "restart_load" %in% runmode ) {
        stmv_db(p=p, DS="load_saved_state", runmode="scale", datasubset="statistics" )
        currentstatus = stmv_statistics_status( p=p)
      }
      p$time_start_runmode = Sys.time()
      p$runmode = paste( "scale_", p$stmv_variogram_nbreaks, sep="" )
      p$clusters = p$stmv_runmode[["scale"]] # as ram reqeuirements increase drop cpus
      message( "\n||| Entering <", p$runmode, "> stage: ", format(Sys.time()) , "\n" )
      currentstatus = stmv_statistics_status( p=p, reset=c("insufficient_data", "variogram_failure", "variogram_range_limit", "unknown" ) )
      parallel_run( stmv_scale, p=p, runindex=list( locs=sample( currentstatus$todo )) )
      stmv_db(p=p, DS="save_current_state", runmode="scale", datasubset="statistics") # temp save to disk
      message( "||| Time used for scale estimation: ", format(difftime(  Sys.time(), p$time_start_runmode )), "\n"  )
    }

    # -----------------------------------------------------
    if ("interpolate" %in% runmode ) {
      stmv_db(p=p, DS="load_saved_state", runmode="scale", datasubset="statistics" )
      if ( "restart_load" %in% runmode ) {
        stmv_db(p=p, DS="load_saved_state", runmode="interpolate", datasubset="predictions" )
        currentstatus = stmv_statistics_status( p=p)
      }
      p$time_start_runmode = Sys.time()
      p0 = p
      ncomplete = -1
      for ( j in 1:length(p$stmv_autocorrelation_interpolation) ) {
        p = p0 #reset
        p$local_interpolation_correlation = p$stmv_autocorrelation_interpolation[j]
        p$runmode = paste("interpolate_", p$local_interpolation_correlation, sep="")
        message( "\n||| Entering <", p$runmode, "> stage: ", format(Sys.time()) , "\n" )
        p$clusters = p$stmv_runmode[["interpolate"]][[j]] # as ram reqeuirements increase drop cpus
        currentstatus = stmv_statistics_status( p=p, reset="incomplete" )
        if ( ncomplete == currentstatus$n.complete ) break()
        ncomplete = currentstatus$n.complete
        if ( length(currentstatus$todo) == 0 ) break()
        if ( length(currentstatus$todo) < length(p$clusters)) p$clusters = p$clusters[1] # drop to serial mode .. otherwise negative indexing occurs
        parallel_run( stmv_interpolate, p=p, runindex=list( locs=sample( currentstatus$todo ))  )
        stmv_db(p=p, DS="save_current_state", runmode="interpolate", datasubset="predictions")
      }
      message( paste( "Time used for <interpolate", j, ">: ", format(difftime(  Sys.time(), p$time_start_runmode )), "\n" ) )
      p = p0
    }
    # stmv_db(p=p, DS="load_saved_state", runmode="interpolate" )
    # stmv_db(p=p, DS="save_current_state", runmode="interpolate")
    # P = stmv_attach( p$storage.backend, p$ptr$P )
    # Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
    # lattice::levelplot( P[] ~ Ploc[,1] + Ploc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso" )

  }


  # --------------------


  if ("interpolate_hybrid_boost" %in% runmode) {
    message ( "not complete .. placeholder")
    message( "\n||| Entering <interpolate_hybrid_boost> stage: ", format(Sys.time()),  "\n" )
    # finalize all interpolations where there are missing data/predictions using
    # interpolation based on data
    # NOTE:: no covariates are used
    if ( "restart_load" %in% runmode ) {
      stmv_db(p=p, DS="load_saved_state", runmode="interpolate_hybrid_boost", datasubset="predictions" )
      currentstatus = stmv_statistics_status( p=p)
    }
    p$time_start_runmode = Sys.time()
    p0 = p
    ncomplete = -1
    for ( j in 1:length(p$stmv_autocorrelation_interpolation) ) {
      p = p0 #reset
      p$local_interpolation_correlation = p$stmv_autocorrelation_interpolation[j]
      p$runmode = paste("interpolate_hybrid_boost_", p$local_interpolation_correlation, sep="")
      message( "\n||| Entering <", p$runmode, "> stage: ", format(Sys.time()) , "\n" )
      p$clusters = p$stmv_runmode[["interpolate_hybrid_boost"]] # as ram reqeuirements increase drop cpus
      p$stmv_local_modelengine = "kernel"  # override -- no covariates, basic moving window average (weighted by inverse variance)
      currentstatus = stmv_statistics_status( p=p, reset="incomplete" )
      if ( ncomplete == currentstatus$n.complete ) break()
      ncomplete = currentstatus$n.complete
      if ( length(currentstatus$todo) == 0 ) break()
      if ( length(currentstatus$todo) < length(p$clusters)) p$clusters = p$clusters[1] # drop to serial mode
      parallel_run( stmv_interpolate, p=p, runindex=list( locs=sample( currentstatus$todo ))  )
      stmv_db(p=p, DS="save_current_state", runmode="interpolate_hybrid_boost", datasubset="predictions")
      message( paste( "Time used for <interpolate_hybrid_boost", j, ">: ", format(difftime(  Sys.time(), p$time_start_runmode )), "\n" ) )
    }
    message( paste( "Time used for <interpolate_hybrid_boost", j, ">: ", format(difftime(  Sys.time(), p$time_start_runmode )), "\n" ) )
    p = p0
    # stmv_db(p=p, DS="load_saved_state", runmode="interpolate_hybrid_boost" )
    # stmv_db(p=p, DS="save_current_state", runmode="interpolate_hybrid_boost")
  }


  # --------------------


  if ("interpolate_force_complete" %in% runmode) {
    message( "\n||| Entering <interpolate force complete> stage: ", format(Sys.time()),  "\n" )
    currentstatus = stmv_statistics_status( p=p, reset="incomplete" )
    p$stmv_force_complete_method = "kernel"
    if ( length(currentstatus$todo) > 0 ) stmv_interpolate_force_complete( p=p )
    stmv_db(p=p, DS="save_current_state", runmode="interpolate_force_complete", datasubset="predictions")  # not needed as this is a terminal step .. but to be consistent
    # stmv_db(p=p, DS="load_saved_state", runmode="interpolate_force_complete", datasubset="predictions")  # not needed as this is a terminal step .. but to be consistent
  }

  # -----------------------------------------------------

  if ("save_completed_data" %in% runmode) {
    stmv_db( p=p, DS="stmv.results" ) # save to disk for use outside stmv*, returning to user scale
  }


  message( "||| Temporary files exist in case you need to examine them or restart a process. ")
  message( "||| You can delete them by running: stmv_db( p=p, DS='cleanup.all' ), once you are done. ")
  message( "||| Or by adding 'cleanup.all' to the runmodes"  )


  if ("cleanup.all" %in% runmode) {
    if ( p$storage.backend !="bigmemory.ram" ) {
      resp = readline( "||| Delete temporary files? Type to confirm <YES>:  ")
      if (resp=="YES")  stmv_db( p=p, DS="cleanup.all" )
    }
  }

  p$time_total = round( difftime( Sys.time(), p$time_start, units="hours" ), 3)
  message( "||| Total time taken (hours):", p$time_total  )
  message( "||| NOTE:: parameter 'p' has been updated in case you need to re-run something" )

}
