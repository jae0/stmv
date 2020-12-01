

stmv = function( p, runmode=NULL, DATA=NULL, nlogs=100, niter=1,
  debug_plot_variable_index=1, debug_data_source="saved.state", debug_plot_log=FALSE, robustify_quantiles=c(0.0005, 0.9995), ... ) {

  if (0) {
    runmode = NULL
    nlogs = 10
    niter = 1
    DATA=NULL
    debug_plot_variable_index=1
    robustify_quantiles=c(0.0005, 0.9995)
    global_sppoly=NULL
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
    runmode = intersect(s_runmode, c( "globalmodel", "scale", "interpolate_fast_predictions", "interpolate", "interpolate_distance_basis", "interpolate_predictions", "save_completed_data", "restart_load", "carstm" ) )
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


  testvars = c(p$stmv_variables$Y, p$stmv_variables$COV, p$stmv_variables$TIME, p$stmv_variables$LOCS)
  if (length(testvars)==1 ) {
    withdata = which( is.finite( DATA$input[, testvars]  ) )
  } else {
    withdata = which( is.finite( rowSums(DATA$input[, testvars] ) ) )
  }
  if (length(withdata) < 1) stop( "Missing data or insufficient data")
  DATA$input = DATA$input[withdata, ]
  withdata=NULL

  # number of time slices
  if (!exists("nt", p)) {
    p$nt = 1  # default to 1 == no time
    if (exists( "ny", p)) p$nt = p$nt * p$ny  # annual time slices
    if (exists( "nw", p)) p$nt = p$nt * p$nw  # sub-annual time slices
  }

  if ( !exists("stmv_tmin", p))  p$stmv_tmin = max(1, aegis_floor( p$nt / 5) )  # min no of time slices in twostep modelling

  # prediction times for space.annual methods, treat time as independent timeslices
  if ( !exists("prediction_ts", p)) p$prediction_ts = 1

  global_model_do = FALSE
  if (p$stmv_global_modelengine !="none" ) global_model_do = TRUE

  message( "\n")
  message( "||| Initializing temporary storage of data and output files... ")
  message( "||| These are large files (4 to 6 X 5GB), it will take a minute ... ")
  message( "||| Try to turn off swap/paging such that only RAM is used. ")

  stmv_db( p=p, DS="cleanup" )

  global_sppoly=NULL
  if (exists("sppoly", DATA)) global_sppoly = sppoly

  Ydata = as.matrix( DATA$input[, p$stmv_variables$Y ] )
  Ydata_datarange = range( Ydata, na.rm=TRUE )
  if (!is.null(robustify_quantiles)) Ydata_datarange = quantile( Ydata, probs=robustify_quantiles, na.rm=TRUE )
  Yq_link = p$stmv_global_family$linkfun( Ydata_datarange ) # convert raw data range to link range

  if (exists("stmv_Y_transform", p)) Ydata = p$stmv_Y_transform$transf( Ydata )
  if (p$storage_backend == "bigmemory.ram" ) {
    tmp_Y = big.matrix( nrow=nrow(Ydata), ncol=1, type="double"  )
    tmp_Y[] = Ydata[]
    p$ptr$Y  = bigmemory::describe( tmp_Y )
  }
  if (p$storage_backend == "bigmemory.filebacked" ) {
    p$ptr$Y  = p$cache$Y
    bigmemory::as.big.matrix( Ydata, type="double", backingfile=basename(p$bm$Y), descriptorfile=basename(p$cache$Y), backingpath=p$stmvSaveDir )
  }
  if (p$storage_backend == "ff" ) {
    p$ptr$Y = ff( Ydata, dim=dim(Ydata), file=p$cache$Y, overwrite=TRUE )
  }

  # if there is a global model update (overwrite) Ydata with residuals
  if ( any(grepl("globalmodel", runmode) ) ) {
    if ( global_model_do ) {
      stmv_global_model( p=p, DS="global_model.redo", B=DATA$input )
      if ( any(grepl("globalmodel.only", runmode)))  return( stmv_global_model( p=p, DS="global_model" ) )
      # model complete .. now predict to get residuals
      if ( p$stmv_global_modelengine != "userdefined" ) Ydata=NULL  # only needed with userdefined models
      Y = stmv_attach( p$storage_backend, p$ptr$Y )
      Y[] = stmv_global_model(p=p, DS="residuals", Yq_link=Yq_link, Ydata=Ydata)[]
      gc()
    }
  }
  Ydata = NULL


  # data coordinates
  Yloc = as.matrix( DATA$input[, p$stmv_variables$LOCS ])
  if (p$storage_backend == "bigmemory.ram" ) {
    tmp_Yloc = big.matrix( nrow=nrow(Yloc), ncol=ncol(Yloc), type="double" )
    tmp_Yloc[] = Yloc[]
    p$ptr$Yloc = bigmemory::describe( tmp_Yloc )
  }
  if (p$storage_backend == "bigmemory.filebacked" ) {
    p$ptr$Yloc  = p$cache$Yloc
    bigmemory::as.big.matrix( Yloc, type="double", backingfile=basename(p$bm$Yloc), descriptorfile=basename(p$cache$Yloc), backingpath=p$stmvSaveDir )
  }
  if (p$storage_backend == "ff" ) {
    p$ptr$Yloc = ff( Yloc, dim=dim(Yloc), file=p$cache$Yloc, overwrite=TRUE )
  }
  Yloc = NULL

  # independent stmv_variables/ covariate
  if (exists("COV", p$stmv_variables)) {
    if (length(p$stmv_variables$COV) > 0) {
      Ycov = as.matrix(  DATA$input[ , p$stmv_variables$COV ] )
      if (p$storage_backend == "bigmemory.ram" ) {
        tmp_Ycov = big.matrix( nrow=nrow(Ycov), ncol=ncol(Ycov), type="double")
        tmp_Ycov[] = Ycov[]
        p$ptr$Ycov  = bigmemory::describe( tmp_Ycov )
      }
      if (p$storage_backend == "bigmemory.filebacked" ) {
        p$ptr$Ycov  = p$cache$Ycov
        bigmemory::as.big.matrix( Ycov, type="double", backingfile=basename(p$bm$Ycov), descriptorfile=basename(p$cache$Ycov), backingpath=p$stmvSaveDir )
      }
      if (p$storage_backend == "ff" ) {
        p$ptr$Ycov = ff( Ycov, dim=dim(Ycov), file=p$cache$Ycov, overwrite=TRUE )
      }
      Ycov= NULL
    }
  }

  # data times
  if ( exists("TIME", p$stmv_variables) ) {
    Ytime = as.matrix(  DATA$input[, p$stmv_variables$TIME ] )
    if (p$storage_backend == "bigmemory.ram" ) {
      tmp_Ytime = big.matrix( nrow=nrow(Ytime), ncol=ncol(Ytime), type="double"  )
      tmp_Ytime[] = Ytime[]
      p$ptr$Ytime  = bigmemory::describe( tmp_Ytime )
    }
    if (p$storage_backend == "bigmemory.filebacked" ) {
      p$ptr$Ytime  = p$cache$Ytime
      bigmemory::as.big.matrix( Ytime, type="double", backingfile=basename(p$bm$Ytime), descriptorfile=basename(p$cache$Ytime), backingpath=p$stmvSaveDir )
    }
    if (p$storage_backend == "ff" ) {
      p$ptr$Ytime = ff( Ytime, dim=dim(Ytime), file=p$cache$Ytime, overwrite=TRUE )
    }
    Ytime =NULL;
  }


  Y = stmv_attach( p$storage_backend, p$ptr$Y )
  Yloc = stmv_attach( p$storage_backend, p$ptr$Yloc )
  Yi = 1:length(Y) # index with useable data
  bad = which( !is.finite( Y[]))
  if (length(bad)> 0 ) Yi[bad] = NA

  # data locations
  bad = which( !is.finite( rowSums(Yloc[])))
  if (length(bad)> 0 ) Yi[bad] = NA

  # data locations
  if (exists("COV", p$stmv_variables)) {
    if (length(p$stmv_variables$COV) > 0) {
      Ycov = stmv_attach( p$storage_backend, p$ptr$Ycov )
      if (length(p$stmv_variables$COV)==1) {
        bad = which( !is.finite( Ycov[] ))
      } else {
        bad = which( !is.finite( rowSums(Ycov[])))
      }
      if (length(bad)> 0 ) Yi[bad] = NA
      Yi = na.omit(Yi)
    }
  }

  # data locations
  if (exists("TIME", p$stmv_variables)) {
    Ytime = stmv_attach( p$storage_backend, p$ptr$Ytime )
    bad = which( !is.finite( Ytime[] ))
    if (length(bad)> 0 ) Yi[bad] = NA
    Yi = na.omit(Yi)
  }
  bad = NULL

  Yi = as.matrix(Yi)
  if (p$storage_backend == "bigmemory.ram" ) {
    tmp_Yi = big.matrix( nrow=nrow(Yi), ncol=ncol(Yi), type="double" )
    tmp_Yi[] = Yi[]
    p$ptr$Yi  = bigmemory::describe( tmp_Yi )
  }
  if (p$storage_backend == "bigmemory.filebacked" ) {
    p$ptr$Yi  = p$cache$Yi
    bigmemory::as.big.matrix( Yi, type="double", backingfile=basename(p$bm$Yi), descriptorfile=basename(p$cache$Yi), backingpath=p$stmvSaveDir )
  }
  if (p$storage_backend == "ff" ) {
    p$ptr$Yi = ff( Yi, dim=dim(Yi), file=p$cache$Yi, overwrite=TRUE )
  }
  Yi=NULL


  #############


  nPloc = nrow(DATA$output$LOCS)

  if (exists("COV", p$stmv_variables)) {
    if (length(p$stmv_variables$COV) > 0) {
      # this needs to be done as Prediction covars need to be structured as lists
      p$ptr$Pcov = list()
      tmp_Pcov = list()
      nc_cov =NULL
      for ( covname in p$stmv_variables$COV ) {
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
        if (p$storage_backend == "bigmemory.ram" ) {
          tmp_Pcov[[covname]] = big.matrix( nrow=nPcovloc, ncol=ncol(Pcovdata), type="double"  )
          tmp_Pcov[[covname]][] = Pcovdata[]
          p$ptr$Pcov[[covname]]  = bigmemory::describe( tmp_Pcov[[covname]] )
        }
        if (p$storage_backend == "bigmemory.filebacked" ) {
          p$ptr$Pcov[[covname]]  = p$cache$Pcov[[covname]]
          bigmemory::as.big.matrix( Pcovdata, type="double", backingfile=basename(p$bm$Pcov[[covname]]), descriptorfile=basename(p$cache$Pcov[[covname]]), backingpath=p$stmvSaveDir )
        }
        if (p$storage_backend == "ff" ) {
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
  if (p$storage_backend == "bigmemory.ram" ) {
    tmp_P = big.matrix( nrow=nrow(sP), ncol=ncol(sP), type="double" )
    tmp_P[] = sP[]
    p$ptr$P  = bigmemory::describe( tmp_P )
  }
  if (p$storage_backend == "bigmemory.filebacked" ) {
    p$ptr$P  = p$cache$P
    bigmemory::as.big.matrix( sP, type="double", backingfile=basename(p$bm$P), descriptorfile=basename(p$cache$P), backingpath=p$stmvSaveDir )
  }
  if (p$storage_backend == "ff" ) {
    p$ptr$P = ff( sP, dim=dim(sP), file=p$cache$P, overwrite=TRUE )
  }
  sP = NULL


  # count of prediction estimates
  sPn = matrix( NaN, nrow=nPloc, ncol=p$nt )
  if (p$storage_backend == "bigmemory.ram" ) {
    tmp_Pn = big.matrix( nrow=nrow(sPn), ncol=ncol(sPn), type="double" )
    tmp_Pn[] = sPn[]
    p$ptr$Pn = bigmemory::describe( tmp_Pn )
  }
  if (p$storage_backend == "bigmemory.filebacked" ) {
    p$ptr$Pn  = p$cache$Pn
    bigmemory::as.big.matrix( sPn, type="double", backingfile=basename(p$bm$Pn), descriptorfile=basename(p$cache$Pn), backingpath=p$stmvSaveDir )
  }
  if (p$storage_backend == "ff" ) {
    p$ptr$Pn = ff( sPn, dim=dim(sPn), file=p$cache$Pn, overwrite=TRUE )
  }
  sPn = NULL

  # sd of prediction estimates
  sPsd = matrix( NaN, nrow=nPloc, ncol=p$nt )
  if (p$storage_backend == "bigmemory.ram" ) {
    tmp_Psd = big.matrix( nrow=nrow(sPsd), ncol=ncol(sPsd), type="double" )
    tmp_Psd[] = sPsd[]
    p$ptr$Psd =bigmemory::describe( tmp_Psd )
  }
  if (p$storage_backend == "bigmemory.filebacked" ) {
    p$ptr$Psd  = p$cache$Psd
    bigmemory::as.big.matrix( sPsd, type="double", backingfile=basename(p$bm$Psd), descriptorfile=basename(p$cache$Psd), backingpath=p$stmvSaveDir )
  }
  if (p$storage_backend == "ff" ) {
    p$ptr$Psd = ff( sPsd, dim=dim(sPsd), file=p$cache$Psd, overwrite=TRUE )
  }
  sPsd = NULL


  # prediction coordinates
  Ploc = as.matrix( DATA$output$LOCS )
  attr( Ploc, "dimnames" ) = NULL
  if (p$storage_backend == "bigmemory.ram" ) {
    tmp_Ploc = big.matrix( nrow=nrow(Ploc), ncol=ncol(Ploc), type="double" )
    tmp_Ploc[] = Ploc[]
    p$ptr$Ploc  = bigmemory::describe( tmp_Ploc )
  }
  if (p$storage_backend == "bigmemory.filebacked" ) {
    p$ptr$Ploc  = p$cache$Ploc
    bigmemory::as.big.matrix( Ploc, type="double", backingfile=basename(p$bm$Ploc), descriptorfile=basename(p$cache$Ploc), backingpath=p$stmvSaveDir )
  }
  if (p$storage_backend == "ff" ) {
    p$ptr$Ploc = ff( Ploc, dim=dim(Ploc), file=p$cache$Ploc, overwrite=TRUE )
  }
  Ploc = NULL;



  ################


    if (global_model_do ) {
      sP0 = matrix( NaN, nrow=nPloc, ncol=p$nt )
      if (!is.null(sP0)) {
        if (p$storage_backend == "bigmemory.ram" ) {
          tmp_P0= big.matrix( nrow=nrow(sP0), ncol=ncol(sP0) , type="double" )
          tmp_P0[] = sP0[]
          p$ptr$P0 = bigmemory::describe(tmp_P0 )
        }
        if (p$storage_backend == "bigmemory.filebacked" ) {
          p$ptr$P0  = p$cache$P0
          bigmemory::as.big.matrix( sP0, type="double", backingfile=basename(p$bm$P0), descriptorfile=basename(p$cache$P0), backingpath=p$stmvSaveDir )
        }
        if (p$storage_backend == "ff" ) {
          p$ptr$P0 = ff( sP0, dim=dim(sP0), file=p$cache$P0, overwrite=TRUE )
        }
        sP0 = NULL
      }

      sP0sd = matrix( NaN, nrow=nPloc, ncol=p$nt )
      if (!is.null(sP0sd)) {
        if (p$storage_backend == "bigmemory.ram" ) {
          tmp_P0sd= big.matrix( nrow=nrow(sP0sd), ncol=ncol(sP0sd) , type="double" )
          tmp_P0sd[] = sP0sd[]
          p$ptr$P0sd = bigmemory::describe(tmp_P0sd )
        }
        if (p$storage_backend == "bigmemory.filebacked" ) {
          p$ptr$P0sd  = p$cache$P0sd
          bigmemory::as.big.matrix( sP0sd, type="double", backingfile=basename(p$bm$P0sd), descriptorfile=basename(p$cache$P0sd), backingpath=p$stmvSaveDir )
        }
        if (p$storage_backend == "ff" ) {
          p$ptr$P0sd = ff( sP0sd, dim=dim(sP0sd), file=p$cache$P0sd, overwrite=TRUE )
        }
        sP0sd=NULL;
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
  if ( exists("COV", p$stmv_variables )) {
    if (length(p$stmv_variables$COV) > 0) {
      niterations = niter
    }
  }

  ##############


  sbox = Sloc = NULL
  if ("carstm" %in% runmode ) {
    if (!is.null(global_sppoly)) {
      # Slocs == Plocs == centroids of global_sppoly
      nSloc = nPloc
      Sloc = as.matrix( coordinates( global_sppoly) )  # centroids , nominal
      nSloc = nrow(Sloc)
      if (nSloc != nPloc) stop( "Global polygons provided seem to have mismatched stats locations.")
    }
  }

  if (is.null(Sloc)) {
    # default:  statistics storage matrix ( aggregation window, coords ) .. no inputs required
    sbox = list(
      plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$stmv_distance_statsgrid ),
      plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$stmv_distance_statsgrid ) )
    # statistics coordinates
    Sloc = as.matrix( expand_grid_fast( sbox$plons, sbox$plats ))
    nSloc = nrow(Sloc)
  }


  if (p$storage_backend == "bigmemory.ram" ) {
    tmp_Sloc = big.matrix(nrow=nSloc, ncol=ncol(Sloc), type="double"  )
    tmp_Sloc[] = Sloc[]
    p$ptr$Sloc  = bigmemory::describe( tmp_Sloc  )
  }
  if (p$storage_backend == "bigmemory.filebacked" ) {
    p$ptr$Sloc  = p$cache$Sloc
    bigmemory::as.big.matrix( Sloc, type="double", backingfile=basename(p$bm$Sloc), descriptorfile=basename(p$cache$Sloc), backingpath=p$stmvSaveDir )
  }
  if (p$storage_backend == "ff" ) {
    p$ptr$Sloc = ff( Sloc, dim=dim(Sloc), file=p$cache$Sloc, overwrite=TRUE )
  }
  Sloc = NULL
  sbox = NULL


  sSflag = matrix( stmv_error_codes()[["todo"]], nrow=nSloc, ncol=1 )
  if (p$storage_backend == "bigmemory.ram" ) {
    tmp_Sflag = big.matrix(nrow=length(sSflag), ncol=1, type="double" )
    tmp_Sflag[] = sSflag[]
    p$ptr$Sflag  = bigmemory::describe( tmp_Sflag )
  }
  if (p$storage_backend == "bigmemory.filebacked" ) {
    p$ptr$Sflag  = p$cache$Sflag
    bigmemory::as.big.matrix( sSflag, type="double", backingfile=basename(p$bm$Sflag), descriptorfile=basename(p$cache$Sflag), backingpath=p$stmvSaveDir )
  }
  if (p$storage_backend == "ff" ) {
    p$ptr$Sflag = ff( sSflag, dim=dim(sSflag), file=p$cache$Sflag, overwrite=TRUE )
  }
  sSflag = NULL

  # determine stats to retain / expect
  res = NULL
  if ( any(grepl("carstm", names(p$stmv_runmode))) ) {
    res = stmv_interpolate_polygons( p=p, global_sppoly=global_sppoly, stmv_au_buffer_links=p$stmv_au_buffer_links, stmv_au_distance_reference=p$stmv_au_distance_reference, just_testing_variablelist=TRUE  )
  } else {
    res = stmv_interpolate_lattice( p=p, just_testing_variablelist=TRUE  )
  }

  # str( res )
#  print( str(res$stmv_stats) )

  # require knowledge of size of stats output which varies with a given type of analysis
  p$statsvars = c("sdTotal", "rsquared", "ndata")
  if (!is.null(res))  p$statsvars = names(res$stmv_stats )
  if (! "carstm" %in% runmode ) p$statsvars = c(p$statsvars,  "sdSpatial", "sdObs", "phi", "nu", "localrange" )
  if (exists("TIME", p$stmv_variables) )  p$statsvars = c( p$statsvars, "ar_timerange", "ar_1" )

  p$statsvars = unique( p$statsvars )

  res = NULL


  sS = matrix( NaN, nrow=nSloc, ncol=length( p$statsvars ) ) # NA forces into logical
  if (p$storage_backend == "bigmemory.ram" ) {
    tmp_S = big.matrix(nrow=nSloc, ncol=length( p$statsvars ), type="double"  )
    tmp_S[] = sS[]
    p$ptr$S  = bigmemory::describe( tmp_S )
  }
  if (p$storage_backend == "bigmemory.filebacked" ) {
    p$ptr$S  = p$cache$S
    bigmemory::as.big.matrix( sS, type="double", backingfile=basename(p$bm$S), descriptorfile=basename(p$cache$S), backingpath=p$stmvSaveDir )
  }
  if (p$storage_backend == "ff" ) {
    p$ptr$S = ff( sS, dim=dim(sS), file=p$cache$S, overwrite=TRUE )
  }
  sS = NULL

  # debug(stmv_statistics_status)

  stmv_statistics_status( p=p, reset="features", verbose=FALSE ) # flags/filter stats locations base dupon prediction covariates. .. speed up and reduce storage

  devold = Inf
  eps = 1e-9


  # NOTE:: must not sink the following memory allocation steps into a deeper function as
  # NOTE:: bigmemory RAM seems to lose the pointers if they are not made simultaneously ?

  # data to be worked upon .. either the raw data or covariate-residuals
  nn = 1

  for ( nn in 1:niterations ) {

    #
    message( "Iteration ", nn, " of ", niterations, "\n" )

    if ( global_model_do ) {
      # create prediction suface .. additive offsets
      # test to see if all covars are static as this can speed up the initial predictions
      message ( "\n", "||| Entering Predicting global effect of covariates at each prediction location: ", format(Sys.time()) , "\n" )
      message( "||| This can take a while (usually a few minutes but hours if complex) ... ")
      p$time_covariates_0 =  Sys.time()
      global_model = stmv_global_model( p=p, DS="global_model")
      if (is.null(global_model)) stop("Global model not found.")
      dev = global_model$deviance
      message("Model Deviance: ", dev)
      if (nn > 1) {
        if (abs(dev - devold)/(0.1 + abs(dev)) < eps ) break() # globalmodel_converged
        iYP = stmv_index_predictions_to_observations(p)  #### this only works for lattice mode .. areal units will require an alternate ... TO DO
        inputdata = P[ iYP ]  # P are the spatial/spatio-temporal random effects (on link scale)
        iYP = NULL
        iYP_nomatch = which(!is.finite(iYP))
        inputdata[ iYP_nomatch ] = 0  # E[RaneFF] = 0
        iYP_nomatch = NULL

        if (exists("linkinv", p$stmv_global_family )) {
          # return to user scale (that of Y)
          inputdata = p$stmv_global_family$linkinv( inputdata )
        }
        global_model$model[, p$stmv_variables$Y ] = global_model$model[, p$stmv_variables$Y ] - inputdata  # update to current estimate of fixed effects
        inputdata = global_model$model
        global_model = NULL
        gc()
        global_model = stmv_global_model( p=p, DS="global_model.redo", B=inputdata,  savedata=FALSE )  # do not overwrite initial model
        inputdata = NULL
      }

      devold = dev
      gc()
      message("Creating fixed effects predictons")
      parallel_run( FUNC=stmv_predict_globalmodel, p=p, runindex=list( pnt=1:p$nt ), Yq_link=Yq_link, global_model=global_model )
      stmv_db(p=p, DS="save_current_state", runmode="meanprocess", datasubset="P0")
      stmv_db(p=p, DS="save_current_state", runmode="meanprocess", datasubset="P0sd")

      Ydata  = residuals(global_model, type="working") # ie. internal (link) scale

      # Yq_link:: could operate upon quantiles of residuals but in poor models this can hyper inflate errors and slow down the whole estimation process
      # truncating using data range as a crude approximation of overall residual and prediction scale
      lb = which( Ydata < Yq_link[1])
      ub = which( Ydata > Yq_link[2])
      if (length(lb) > 0) Ydata[lb] = Yq_link[1]
      if (length(ub) > 0) Ydata[ub] = Yq_link[2]

      Y = stmv_attach( p$storage_backend, p$ptr$Y )
      Y[] = Ydata[]  # update "Ydata" ... ( residuals)
      Ydata = NULL
      global_model = NULL
      p$time_covariates = round(difftime( Sys.time(), p$time_covariates_0 , units="hours"), 3)
      message( paste( "||| Time taken to predict covariate surface (hours):", p$time_covariates ) )
    }



    # -----------------------------------------------------
    if ("carstm" %in% runmode ) {
      if (!exists("stmv_au_distance_reference", p)) p$stmv_au_distance_reference="none"
      if (!exists("stmv_au_buffer_links", p)) p$stmv_au_buffer_links=0
      # this models, and predicts in same step (via inla)
      message ( "\n", "||| Entering carstm areal unit modelling: ", format(Sys.time())  )
      if ( "restart_load" %in% runmode ) {
        invisible( stmv_db(p=p, DS="load_saved_state", runmode="carstm", datasubset="statistics" ) )
        invisible( stmv_db(p=p, DS="load_saved_state", runmode="carstm", datasubset="P" ) )
        invisible( stmv_db(p=p, DS="load_saved_state", runmode="carstm", datasubset="Psd" ) )
        invisible( stmv_db(p=p, DS="load_saved_state", runmode="carstm", datasubset="Pn" ) )
        if ( global_model_do ) invisible( stmv_db(p=p, DS="load_saved_state", runmode="meanprocess",  datasubset="P0" ) )
        if ( global_model_do ) invisible( stmv_db(p=p, DS="load_saved_state", runmode="meanprocess",  datasubset="P0sd" ) )
        stmv_statistics_status( p=p, reset=c( "all"), verbose=FALSE  ) # required to start as scale determination uses Sflags too
        stmv_statistics_status( p=p, reset=c( "complete" ), verbose=FALSE  )
        stmv_statistics_status( p=p, reset=c( "incomplete" ), verbose=FALSE  )
        stmv_statistics_status( p=p, reset=c( "features" ), verbose=FALSE  ) # required to start as scale determination uses Sflags too
      }

      p$time_start_runmode = Sys.time()
      p0 = p
      for ( j in 1:length(p$stmv_interpolation_basis_distance_choices) ) {
        p = p0 #reset
        p$stmv_interpolation_basis_distance = p$stmv_interpolation_basis_distance_choices[j]
        # p$runmode = paste("carstm_distance_basis_", p$stmv_interpolation_basis_correlation, sep="")
        # ni = p$stmv_runmode[["interpolate"]]
        # jcpu = ifelse( j > ni, ni, j )
        # p$clusters = p$stmv_runmode[["carstm"]][[jcpu]] # as ram reqeuirements increase drop cpus
        p$runmode = "carstm"
        p$clusters = p$stmv_runmode[["carstm"]] # as ram reqeuirements increase drop cpus
        currentstatus = stmv_statistics_status( p=p, reset="flags", reset_flags=c("insufficient_data",  "unknown" ) )
        parallel_run( stmv_interpolate_polygons, p=p, runmode="carstm", global_sppoly=global_sppoly, stmv_au_buffer_links=p$stmv_au_buffer_links,stmv_au_distance_reference=p$stmv_au_distance_reference, runindex=list( locs=sample( currentstatus$todo )) )
        invisible( stmv_db(p=p, DS="save_current_state", runmode="carstm", datasubset="statistics") )# temp save to disk
        invisible( stmv_db(p=p, DS="save_current_state", runmode="carstm", datasubset="P" ) )
        invisible( stmv_db(p=p, DS="save_current_state", runmode="carstm", datasubset="Psd" ) )
        invisible( stmv_db(p=p, DS="save_current_state", runmode="carstm", datasubset="Pn" ) )
        stmv_statistics_status( p=p, verbose=FALSE ) # quick update before logging
        slog = stmv_logfile(p=p, flag= paste("Carstm phase", p$runmode, "completed ...") ) # final update before continuing
      }

      message( "||| Time used for carstm estimation: ", format(difftime(  Sys.time(), p$time_start_runmode )), "\n"  )
    }



    # -----------------------------------------------------
    if ( "scale" %in% runmode ) {
      message ( "\n", "||| Entering spatial scale (variogram) determination: ", format(Sys.time())  )
      if ( "restart_load" %in% runmode ) {
        invisible( stmv_db(p=p, DS="load_saved_state", runmode="scale", datasubset="statistics" ) )
      }
      p$time_start_runmode = Sys.time()
      p$runmode = "scale"
      p$clusters = p$stmv_runmode[["scale"]] # as ram reqeuirements increase drop cpus
      currentstatus = stmv_statistics_status( p=p, reset="flags", reset_flags=c("insufficient_data", "variogram_failure", "variogram_range_limit", "unknown" ) )
      parallel_run( stmv_scale, p=p, runindex=list( locs=sample( currentstatus$todo )) )
      stmv_db(p=p, DS="save_current_state", runmode="scale", datasubset="statistics") # temp save to disk
      stmv_statistics_status( p=p, verbose=FALSE ) # quick update before logging
      slog = stmv_logfile(p=p, flag= paste("Interpolation phase", p$runmode, "completed ...") ) # final update before continuing
      message( "||| Time used for scale estimation: ", format(difftime(  Sys.time(), p$time_start_runmode )), "\n"  )
    }



    # -----------------------------------------------------
    if ("interpolate" %in% runmode ) {

      invisible( stmv_db(p=p, DS="load_saved_state", runmode="scale", datasubset="statistics" ))
      if ( "restart_load" %in% runmode ) {
        invisible( stmv_db(p=p, DS="load_saved_state", runmode=p$restart_load, datasubset="P" ) )
        invisible( stmv_db(p=p, DS="load_saved_state", runmode=p$restart_load, datasubset="Psd" ) )
        invisible( stmv_db(p=p, DS="load_saved_state", runmode=p$restart_load, datasubset="Pn" ) )
        if ( global_model_do ) invisible( stmv_db(p=p, DS="load_saved_state", runmode="meanprocess",  datasubset="P0" ) )
        if ( global_model_do ) invisible( stmv_db(p=p, DS="load_saved_state", runmode="meanprocess",  datasubset="P0sd" ) )
        stmv_statistics_status( p=p, reset=c( "all"), verbose=FALSE  ) # required to start as scale determination uses Sflags too
        stmv_statistics_status( p=p, reset=c( "complete" ), verbose=FALSE  )
        stmv_statistics_status( p=p, reset=c( "incomplete" ), verbose=FALSE  )
        stmv_statistics_status( p=p, reset=c( "features" ), verbose=FALSE  ) # required to start as scale determination uses Sflags too
      }

      p$time_start_runmode = Sys.time()
      p$stmv_interpolation_basis = "correlation"
      p0 = p
      for ( j in 1:length(p$stmv_autocorrelation_basis_interpolation) ) {
        p = p0 #reset
        p$stmv_interpolation_basis_correlation = p$stmv_autocorrelation_basis_interpolation[j]
        p$runmode = paste("interpolate_correlation_basis_", p$stmv_interpolation_basis_correlation, sep="")
        ni = length( p$stmv_runmode[["interpolate"]] )
        jcpu = ifelse( j > ni, ni, j )
        p$clusters = p$stmv_runmode[["interpolate"]][[jcpu]] # as ram reqeuirements increase drop cpus

        if (exists("stmv_fft_filter", p)) {
          if (grepl("fast_and_exhaustive_predictions", p$stmv_fft_filter)) {
            # do a quick cycle first
            p$stmv_fft_filter = paste( p$stmv_fft_filter, "fast_predictions" )
            message( "\n||| Entering < Fast ", p$runmode, " > : ", format(Sys.time()) )
            currentstatus = stmv_statistics_status( p=p, reset=c( "incomplete" ) ) # flags/filter stats locations base dupon prediction covariates. .. speed up and reduce storage
            if ( currentstatus$n.todo == 0 ) break()
            if ( currentstatus$n.todo < (2*length(p$clusters)) ) p$clusters = p$clusters[1] # drop to serial mode
            invisible( parallel_run( stmv_interpolate_lattice, p=p, runindex=list( locs=sample( currentstatus$todo )) ) )
            invisible( stmv_db(p=p, DS="save_current_state", runmode=p$runmode, datasubset="P" ) )
            invisible( stmv_db(p=p, DS="save_current_state", runmode=p$runmode, datasubset="Psd" ) )
            invisible( stmv_db(p=p, DS="save_current_state", runmode=p$runmode, datasubset="Pn" ) )
            stmv_statistics_status( p=p, verbose=FALSE ) # quick update before logging
            slog = stmv_logfile(p=p, flag= paste("Fast Interpolation correlation basis phase", p$runmode, "completed ...") ) # final update before continuing
            p$stmv_fft_filter = gsub( "fast_predictions", "", p$stmv_fft_filter )  # remove fast ... now moving to slow/exhaustive
          }
        }
        message( "\n||| Entering < Exhaustive ", p$runmode, " > : ", format(Sys.time()) )
        currentstatus = stmv_statistics_status( p=p, reset=c( "incomplete" ) ) # flags/filter stats locations base dupon prediction covariates. .. speed up and reduce storage
        if ( currentstatus$n.todo == 0 ) break()
        if ( currentstatus$n.todo < (2*length(p$clusters)) ) p$clusters = p$clusters[1] # drop to serial mode
        invisible( parallel_run( stmv_interpolate_lattice, p=p, runindex=list( locs=sample( currentstatus$todo )) ) )
        invisible( stmv_db(p=p, DS="save_current_state", runmode=p$runmode, datasubset="P" ) )
        invisible( stmv_db(p=p, DS="save_current_state", runmode=p$runmode, datasubset="Psd" ) )
        invisible( stmv_db(p=p, DS="save_current_state", runmode=p$runmode, datasubset="Pn" ) )
        stmv_statistics_status( p=p, verbose=FALSE ) # quick update before logging
        slog = stmv_logfile(p=p, flag= paste("Exhaustive Interpolation correlation basis phase", p$runmode, "completed ...") ) # final update before continuing

      }
      message( paste( "Time used for <interpolations", ">: ", format(difftime(  Sys.time(), p$time_start_runmode )), "\n" ) )

      invisible( stmv_db(p=p, DS="save_current_state", runmode="interpolate", datasubset="P" ) )
      invisible( stmv_db(p=p, DS="save_current_state", runmode="interpolate", datasubset="Psd" ) )
      invisible( stmv_db(p=p, DS="save_current_state", runmode="interpolate", datasubset="Pn" ) )

      p = p0

    }

    if(0) {
      stmv_db(p=p, DS="load_saved_state", runmode="interpolate" )
      stmv_db(p=p, DS="save_current_state", runmode="interpolate")
      P = stmv_attach( p$storage_backend, p$ptr$P )
      Ploc = stmv_attach( p$storage_backend, p$ptr$Ploc )
      if (length(dim(P)) > 1 ) {
        for (i in 1:p$nt) print(lattice::levelplot( P[,i] ~ Ploc[,1] + Ploc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso" ))
      } else {
        lattice::levelplot( P[] ~ Ploc[,1] + Ploc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso" )
      }
    }



    # -----------------------------------------------------
    if ("interpolate_distance_basis" %in% runmode ) {

      invisible( stmv_db(p=p, DS="load_saved_state", runmode="scale", datasubset="statistics" ))
      if ( "restart_load" %in% runmode ) {
        invisible( stmv_db(p=p, DS="load_saved_state", runmode=p$restart_load, datasubset="P" ) )
        invisible( stmv_db(p=p, DS="load_saved_state", runmode=p$restart_load, datasubset="Psd" ) )
        invisible( stmv_db(p=p, DS="load_saved_state", runmode=p$restart_load, datasubset="Pn" ) )
        if ( global_model_do ) invisible( stmv_db(p=p, DS="load_saved_state", runmode="meanprocess",  datasubset="P0" ) )
        if ( global_model_do ) invisible( stmv_db(p=p, DS="load_saved_state", runmode="meanprocess",  datasubset="P0sd" ) )
        stmv_statistics_status( p=p, reset=c( "all"), verbose=FALSE  ) # required to start as scale determination uses Sflags too
        stmv_statistics_status( p=p, reset=c( "complete" ), verbose=FALSE  )
        stmv_statistics_status( p=p, reset=c( "incomplete" ), verbose=FALSE  )
        stmv_statistics_status( p=p, reset=c( "features" ), verbose=FALSE  ) # required to start as scale determination uses Sflags too
      }

      p$time_start_runmode = Sys.time()
      p$stmv_interpolation_basis = "distance"
      p0 = p
      for ( j in 1:length(p$stmv_distance_basis_interpolation) ) {
        p = p0 #reset
        p$stmv_interpolation_basis_distance = p$stmv_distance_basis_interpolation[j]
        p$runmode = paste("interpolate_distance_basis_", p$stmv_interpolation_basis_distance, sep="")
        ni = length( p$stmv_runmode[["interpolate"]] )
        jcpu = ifelse( j > ni, ni, j )
        p$clusters = p$stmv_runmode[["interpolate_distance_basis"]][[jcpu]] # as ram reqeuirements increase drop cpus
        message( "\n||| Entering <", p$runmode, " > : ", format(Sys.time()) )
        currentstatus = stmv_statistics_status( p=p, reset=c( "incomplete" ) ) # flags/filter stats locations base dupon prediction covariates. .. speed up and reduce storage
        if ( currentstatus$n.todo == 0 ) break()
        if ( currentstatus$n.todo < (2*length(p$clusters)) ) p$clusters = p$clusters[1] # drop to serial mode
        invisible( parallel_run( stmv_interpolate_lattice, p=p, runindex=list( locs=sample( currentstatus$todo ) )  ) )
        invisible( stmv_db(p=p, DS="save_current_state", runmode=p$runmode, datasubset="P" ) )
        invisible( stmv_db(p=p, DS="save_current_state", runmode=p$runmode, datasubset="Psd" ) )
        invisible( stmv_db(p=p, DS="save_current_state", runmode=p$runmode, datasubset="Pn" ) )
        stmv_statistics_status( p=p, verbose=FALSE ) # quick update before logging
        slog = stmv_logfile(p=p, flag= paste("Interpolation distance basis phase", p$runmode, "completed ...") ) # final update before continuing
      }
      message( paste( "Time used for <interpolations", ">: ", format(difftime(  Sys.time(), p$time_start_runmode )), "\n" ) )
      stmv_db(p=p, DS="save_current_state", runmode="interpolate_distance_basis" )
      p = p0
    }

    if(0) {
      stmv_db(p=p, DS="load_saved_state", runmode="interpolate_distance_basis" )
      stmv_db(p=p, DS="save_current_state", runmode="interpolate_distance_basis" )
      P = stmv_attach( p$storage_backend, p$ptr$P )
      Ploc = stmv_attach( p$storage_backend, p$ptr$Ploc )
      if (length(dim(P)) > 1 ) {
        for (i in 1:p$nt) print(lattice::levelplot( P[,i] ~ Ploc[,1] + Ploc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso" ))
      } else {
        lattice::levelplot( P[] ~ Ploc[,1] + Ploc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso" )
      }
    }


  }  #end loop

  # --------------------


  if ("interpolate_hybrid_boost" %in% runmode) {
    message ( "not complete .. placeholder")
    if (0) {
        message( "\n||| Entering <interpolate_hybrid_boost> stage: ", format(Sys.time()),  "\n" )
        # finalize all interpolations where there are missing data/predictions using
        # interpolation based on data
        # NOTE:: no covariates are used
        # invisible( stmv_db(p=p, DS="load_saved_state", runmode="interpolation", datasubset="P" ) )
        # invisible( stmv_db(p=p, DS="load_saved_state", runmode="interpolation", datasubset="Psd" ) )
        # invisible( stmv_db(p=p, DS="load_saved_state", runmode="interpolation", datasubset="Pn" ) )
        # if ( global_model_do ) invisible( stmv_db(p=p, DS="load_saved_state", runmode="meanprocess",  datasubset="P0" ) )
        # if ( global_model_do ) invisible( stmv_db(p=p, DS="load_saved_state", runmode="meanprocess",  datasubset="P0sd" ) )
        if ( "restart_load" %in% runmode ) {
          invisible( stmv_db(p=p, DS="load_saved_state", runmode=p$restart_load, datasubset="P" ) )
          invisible( stmv_db(p=p, DS="load_saved_state", runmode=p$restart_load, datasubset="Psd" ) )
          invisible( stmv_db(p=p, DS="load_saved_state", runmode=p$restart_load, datasubset="Pn" ) )
          if ( global_model_do ) invisible( stmv_db(p=p, DS="load_saved_state", runmode="meanprocess",  datasubset="P0" ) )
          if ( global_model_do ) invisible( stmv_db(p=p, DS="load_saved_state", runmode="meanprocess",  datasubset="P0sd" ) )
          stmv_statistics_status( p=p, reset=c( "all"), verbose=FALSE  ) # required to start as scale determination uses Sflags too
          stmv_statistics_status( p=p, reset=c( "complete" ), verbose=FALSE  )
          stmv_statistics_status( p=p, reset=c( "incomplete" ), verbose=FALSE  )
          stmv_statistics_status( p=p, reset=c( "features" ), verbose=FALSE  ) # required to start as scale determination uses Sflags too
        }
        p$time_start_runmode = Sys.time()
        p0 = p
        for ( j in 1:length(p$stmv_autocorrelation_basis_interpolation) ) {
          p = p0 #reset
          p$local_interpolation_correlation = p$stmv_autocorrelation_basis_interpolation[j]
          p$runmode = paste("interpolate_hybrid_boost_", p$local_interpolation_correlation, sep="")
          ni = length( p$stmv_runmode[["interpolate"]] )
          jcpu = ifelse( j > ni, ni, j )
          p$clusters = p$stmv_runmode[["interpolate_hybrid_boost"]][[jcpu]] # as ram reqeuirements increase drop cpus
          message( "\n||| Entering <", p$runmode, "> stage: ", format(Sys.time()) , "\n" )
          p$stmv_local_modelengine = "kernel"  # override -- no covariates, basic moving window average (weighted by inverse variance)
          currentstatus = stmv_statistics_status( p=p, reset=c( "incomplete" ) ) # flags/filter stats locations base dupon prediction covariates. .. speed up and reduce storage
          if ( currentstatus$n.todo == 0 ) break()
          if ( currentstatus$n.todo < (2*length(p$clusters)) ) p$clusters = p$clusters[1] # drop to serial mode
          invisible( parallel_run( stmv_interpolate_lattice, p=p, runindex=list( locs=sample( currentstatus$todo ))  ) )
          invisible( stmv_db(p=p, DS="save_current_state", runmode=p$runmode, datasubset="P" ) )
          invisible( stmv_db(p=p, DS="save_current_state", runmode=p$runmode, datasubset="Psd" ) )
          invisible( stmv_db(p=p, DS="save_current_state", runmode=p$runmode, datasubset="Pn" ) )
          stmv_statistics_status( p=p, verbose=FALSE ) # quick update before logging
          slog = stmv_logfile(p=p, flag= paste("Interpolation phase", p$runmode, "completed ...") ) # final update
        }
        message( paste( "Time used for <interpolate_hybrid_boost", j, ">: ", format(difftime(  Sys.time(), p$time_start_runmode )), "\n" ) )
        stmv_db(p=p, DS="save_current_state", runmode="interpolate_hybrid_boost")  # final save for this runmode
        # stmv_db(p=p, DS="load_saved_state", runmode="interpolate_hybrid_boost" )
        p = p0
    }
  }


  # --------------------


  if ("interpolate_predictions" %in% runmode) {

      message( "\n||| Entering <interpolate force complete> stage: ", format(Sys.time()),  "\n" )
      invisible( stmv_db(p=p, DS="load_saved_state", runmode="scale", datasubset="statistics" ))
      if ( "restart_load" %in% runmode ) {
        invisible( stmv_db(p=p, DS="load_saved_state", runmode="interpolate_predictions", datasubset="P" ) )
        invisible( stmv_db(p=p, DS="load_saved_state", runmode="interpolate_predictions", datasubset="Psd" ) )
        invisible( stmv_db(p=p, DS="load_saved_state", runmode="interpolate_predictions", datasubset="Pn" ) )
        if ( global_model_do ) invisible( stmv_db(p=p, DS="load_saved_state", runmode="meanprocess",  datasubset="P0" ) )
        if ( global_model_do ) invisible( stmv_db(p=p, DS="load_saved_state", runmode="meanprocess",  datasubset="P0sd" ) )
        stmv_statistics_status( p=p, reset=c( "all"), verbose=FALSE  ) # required to start as scale determination uses Sflags too
        stmv_statistics_status( p=p, reset=c( "complete" ), verbose=FALSE  )
        stmv_statistics_status( p=p, reset=c( "incomplete" ), verbose=FALSE  )
        stmv_statistics_status( p=p, reset=c( "features" ), verbose=FALSE  ) # required to start as scale determination uses Sflags too
      }
      p$time_start_runmode = Sys.time()
      p0 = p
      for ( j in 1:length(p$stmv_distance_scale) ) {
        p = p0 #reset
        p$stmv_interpolation_basis_distance = p$stmv_distance_scale[j]
        p$runmode = "interpolate_predictions"
        ni = length( p$stmv_runmode[["interpolate"]] )
        jcpu = ifelse( j > ni, ni, j )
        p$clusters = p$stmv_runmode[["interpolate"]][[jcpu]] # as ram reqeuirements increase drop cpus
        stmv_statistics_status( p=p, reset=c( "complete" ), verbose=FALSE )
        currentstatus = stmv_statistics_status( p=p, reset=c( "incomplete" ) ) # flags/filter stats locations base dupon prediction covariates. .. speed up and reduce storage
        if ( currentstatus$n.todo == 0 ) break()
        if ( currentstatus$n.todo < (2*length(p$clusters)) ) p$clusters = p$clusters[1] # drop to serial mode
        invisible( parallel_run( stmv_interpolate_predictions, p=p, runindex=list( locs=sample( currentstatus$todo )) ) )
        stmv_statistics_status( p=p, verbose=FALSE ) # quick update before logging
        slog = stmv_logfile(p=p, flag= paste("Exhaustive Interpolation correlation basis phase", p$runmode, "completed ...") ) # final update before continuing
        invisible( stmv_db(p=p, DS="save_current_state", runmode="interpolate_predictions", datasubset="P" ) )
        invisible( stmv_db(p=p, DS="save_current_state", runmode="interpolate_predictions", datasubset="Psd" ) )
        invisible( stmv_db(p=p, DS="save_current_state", runmode="interpolate_predictions", datasubset="Pn" ) )

      }
      message( paste( "Time used for <interpolations", ">: ", format(difftime(  Sys.time(), p$time_start_runmode )), "\n" ) )
      p = p0
  }


  # -----------------------------------------------------

  if ("save_completed_data" %in% runmode) {
    message( "\n||| Saving final results: ", format(Sys.time()),  "\n" )
    stmv_db( p=p, DS="stmv.results" ) # save to disk for use outside stmv*, returning to user scale
  }


  message( "||| Temporary files exist in case you need to examine them or restart a process. ")
  message( "||| You can delete them by running: stmv_db( p=p, DS='cleanup.all' ), once you are done. ")
  message( "||| Or by adding 'cleanup.all' to the runmodes"  )


  if ("cleanup.all" %in% runmode) {
    if ( p$storage_backend !="bigmemory.ram" ) {
      resp = readline( "||| Delete temporary files? Type to confirm <YES>:  ")
      if (resp=="YES")  stmv_db( p=p, DS="cleanup.all" )
    }
  }

  message( "||| NOTE:: parameter 'p' has been updated in case you need to re-run something \n" )

  p$time_total = round( difftime( Sys.time(), p$time_start, units="hours" ), 3)
  message( "||| Total time taken (hours):", p$time_total  )

}
