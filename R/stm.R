

stm = function( p, runmode="default", DATA=NULL, storage.backend="bigmemory.ram",  debug_plot_variable_index=1 ) {

  #\\ localized modelling of space and time data to predict/interpolate upon a grid 
  #\\ speed ratings: bigmemory.ram (1), ff (2), bigmemory.filebacked (3)

  # TODO: splancs::kernel3d as a method ? .. for count data?
  # TODO: gaussian process // see stm_interpolate_xy_simple
  #       .. the convoSPAT seems almost fast enough
  # TODO: MBA mba.surf method? ... seems very fast

  if (!exists("time.start", p) ) p$time.start = Sys.time()

  p$savedir = file.path(p$data_root, "modelled", p$variables$Y, p$spatial.domain )
  if ( !file.exists(p$savedir)) dir.create( p$savedir, recursive=TRUE, showWarnings=FALSE )
  message( "||| In case something should go wrong, intermediary outputs will be placed at:" )
  message( "|||",  p$savedir  )

  # determine storage format
  p$libs = unique( c( p$libs, "sp", "rgdal", "parallel", "RandomFields", "geoR" ) )
  if (!exists("storage.backend", p)) p$storage.backend = storage.backend
  if (any( grepl ("ff", p$storage.backend)))         p$libs = c( p$libs, "ff", "ffbase" )
  if (any( grepl ("bigmemory", p$storage.backend)))  p$libs = c( p$libs, "bigmemory" )
  if (p$storage.backend=="bigmemory.ram") {
    if ( length( unique(p$clusters)) > 1 ) {
      stop( "||| More than one unique cluster server was specified .. the RAM-based method only works within one server." )
    }
  }
  
  # other libs
  if (exists("stm_local_modelengine", p)) {
    if (p$stm_local_modelengine=="bayesx")  p$libs = c( p$libs, "R2BayesX" )
    if (p$stm_local_modelengine %in% c("gam", "mgcv", "habitat") )  p$libs = c( p$libs, "mgcv" )
    if (p$stm_local_modelengine %in% c("LaplacesDemon") )  p$libs = c( p$libs, "LaplacesDemonCpp" )
    if (p$stm_local_modelengine %in% c("inla") )  p$libs = c( p$libs, "INLA" )
    if (p$stm_local_modelengine %in% c("fft", "gaussianprocess2Dt") )  p$libs = c( p$libs, "fields" )
    if (p$stm_local_modelengine %in% c("gaussianprocess") )  p$libs = c( p$libs  )
    if (p$stm_local_modelengine %in% c("splancs") )  p$libs = c( p$libs, "splancs" )
    if (p$stm_local_modelengine %in% c("twostep") )  p$libs = c( p$libs, "mgcv", "fields" )
    if (p$stm_local_modelengine %in% c("krige") ) p$libs = c( p$libs, "fields" )
    if (p$stm_local_modelengine %in% c("gstat") ) p$libs = c( p$libs, "gstat" )
    if (p$stm_local_modelengine %in% c("stan") ) p$libs = c( p$libs, "rstan" )
    # if (p$stm_local_modelengine %in% c("spate") )  p$libs = c( p$libs, "spate" ) # now copied directly into stm
  }

  if (exists("stm_global_modelengine", p)) {
    if (p$stm_global_modelengine %in% c("gam", "mgcv") ) p$libs = c( p$libs, "mgcv" )
    if (p$stm_global_modelengine %in% c("bigglm", "biglm") ) p$libs = c( p$libs, "biglm" )
  }


  p$libs = unique( p$libs )
  suppressMessages( RLibrary( p$libs ) )

  
  if (is.null(DATA) ) {
    # in here as it assumes that initiation of process (save of data files) is complete and "saved" objects are present
    message( "||| No DATA provided, assuming we are continuing from an interrupted start" )
    message( "|||  and using parameters from saved configuration:", file.path( p$savedir, 'p.rdata' ) )
    message( "|||  Delete this if you want to over-ride these settings.")

    p = stm_db( p=p, DS="load.parameters" )  
    stm_db(p=p, DS="statistics.status.reset" )

  }  else {
    # not a restart .. new instance:

    p$stm_current_status = file.path( p$savedir, "stm_current_status" ) 

    p = stm_parameters(p=p) # fill in parameters with defaults where possible
    p = stm_db( p=p, DS="filenames" )
    p$ptr = list() # location for data pointers

    # set up the data and problem using data objects

    tmpfiles = unlist( p$cache)
    for (tf in tmpfiles) if (file.exists( tf)) file.remove(tf)

    p$variables$all = NULL
    if (exists("stm_local_modelformula", p))  {
      p$variables$local_all = all.vars( p$stm_local_modelformula )
      p$variables$local_cov = intersect( p$variables$local_all, p$variables$COV ) 
      p$variables$all = unique( c( p$variables$all, p$variables$local_all ) )
    }
    if (exists("stm_global_modelformula", p)) {
      p$variables$global_all = all.vars( p$stm_global_modelformula )
      p$variables$global_cov = intersect( p$variables$global_all, p$variables$COV )      
      p$variables$all = unique( c( p$variables$all, p$variables$global_all ) )
    }

    # permit passing a function rather than data directly .. less RAM usage in parent call
    if (class(DATA)=="character") assign("DATA", eval(parse(text=DATA) ) )
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

    # require knowledge of size of stats output before create S, which varies with a given type of analysis
    othervars = c( )
    if (exists("stm_local_modelengine", p)) {
      if (p$stm_local_modelengine == "habitat") othervars = c( )
      if (p$stm_local_modelengine == "spate") othervars = c( 
        "rho_0",  "zeta", "rho_1", "gamma", "alpha", "mu_x", "mu_y", "sigma^2", "tau^2", 
        "rho_0.sd", "zeta.sd", "rho_1.sd", "gamma.sd", "alpha.sd", "mu_x.sd", "mu_y.sd" )
    }
    if (exists("TIME", p$variables) )  othervars = c( "ar_timerange", "ar_1" )
    p$statsvars = unique( c( "sdTotal", "rsquared", "ndata", "sdSpatial", "sdObs", "range", "phi", "nu", othervars ) )

    message("||| ")
    message( "||| Initializing temporary storage of data and outputs files... ")
    message( "||| These are large files (4 to 6 X 5GB), it will take a minute ... ")
    stm_db( p=p, DS="cleanup" )
  
    p$nloccov = 0
    if (exists("local_cov", p$variables)) p$nloccov = length(p$variables$local_cov)

    # construct prediction/output grid area ('pa')
    p$windowsize.half = floor(p$stm_distance_prediction/p$pres) # convert distance to discretized increments of row/col indices

    if (exists("stm_Y_transform", p)) {
      DATA$input[, p$variables$Y ] = p$stm_Y_transform[[1]] (DATA$input[, p$variables$Y ] )
    }


    if (  exists("stm_global_modelformula", p) ) {
      # to add global covariate model (hierarchical) 
      # .. simplistic this way but faster ~ kriging with external drift
        stm_db( p=p, DS="global_model.redo", B=DATA$input )
    }


    # NOTE:: must not sink the following memory allocation into a deeper funcion as 
    # NOTE:: bigmemory RAM seems to lose the pointers if they are not made simultaneously 
 
    # init output data objects
    # statistics storage matrix ( aggregation window, coords ) .. no inputs required
    sbox = list( 
      plats = seq( p$corners$plat[1], p$corners$plat[2], by=p$stm_distance_statsgrid ),
      plons = seq( p$corners$plon[1], p$corners$plon[2], by=p$stm_distance_statsgrid ) )

      # statistics coordinates
      Sloc = as.matrix( expand.grid( sbox$plons, sbox$plats ))
        if (p$storage.backend == "bigmemory.ram" ) {
          tmp_Sloc = big.matrix(nrow=nrow(Sloc), ncol=ncol(Sloc), type="double"  )
          tmp_Sloc[] = Sloc
          p$ptr$Sloc  = bigmemory::describe( tmp_Sloc  )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Sloc  = p$cache$Sloc
          bigmemory::as.big.matrix( Sloc, type="double", backingfile=basename(p$bm$Sloc), descriptorfile=basename(p$cache$Sloc), backingpath=p$savedir )
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
          bigmemory::as.big.matrix( S, type="double", backingfile=basename(p$bm$S), descriptorfile=basename(p$cache$S), backingpath=p$savedir )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$S = ff( S, dim=dim(S), file=p$cache$S, overwrite=TRUE )
        }


      Sflag = matrix( 0L, nrow=nrow(Sloc), ncol=1 )  # 0L is the todo flag
        if (p$storage.backend == "bigmemory.ram" ) {
          tmp_Sflag = big.matrix(nrow=nrow(Sloc), ncol=1, type="double" )
          tmp_Sflag[] = 0L # TODO flag
          p$ptr$Sflag  = bigmemory::describe( tmp_Sflag )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Sflag  = p$cache$Sflag
          bigmemory::as.big.matrix( Sflag, type="double", backingfile=basename(p$bm$Sflag), descriptorfile=basename(p$cache$Sflag), backingpath=p$savedir )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Sflag = ff( Sflag, dim=dim(Sflag), file=p$cache$Sflag, overwrite=TRUE )
        }

      rm(S, Sflag, Sloc)


      # data to be worked upon .. either the raw data or covariate-residuals

      Ydata = as.matrix(DATA$input[, p$variables$Y ])
      if (exists("stm_global_modelengine", p)) {
        covmodel = stm_db( p=p, DS="global_model")
        Ypreds = predict(covmodel, type="link", se.fit=FALSE )  ## TODO .. keep track of the SE 
        Ydata  = residuals(covmodel, type="deviance") # ie. link scale .. this is the default but make it explicit 
        covmodel =NULL; gc()
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
          bigmemory::as.big.matrix( Ydata, type="double", backingfile=basename(p$bm$Y), descriptorfile=basename(p$cache$Y), backingpath=p$savedir )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Y = ff( Ydata, dim=dim(Ydata), file=p$cache$Y, overwrite=TRUE )
        }

      rm(Ydata)

      Y = stm_attach( p$storage.backend, p$ptr$Y )

     # data coordinates
      Yloc = as.matrix( DATA$input[, p$variables$LOCS ])
        if (p$storage.backend == "bigmemory.ram" ) {
          tmp_Yloc = big.matrix( nrow=nrow(Yloc), ncol=ncol(Yloc), type="double" )
          tmp_Yloc[] = Yloc
          p$ptr$Yloc = bigmemory::describe( tmp_Yloc )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$Yloc  = p$cache$Yloc
          bigmemory::as.big.matrix( Yloc, type="double", backingfile=basename(p$bm$Yloc), descriptorfile=basename(p$cache$Yloc), backingpath=p$savedir )
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
            bigmemory::as.big.matrix( Ycov, type="double", backingfile=basename(p$bm$Ycov), descriptorfile=basename(p$cache$Ycov), backingpath=p$savedir )
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
            bigmemory::as.big.matrix( Ytime, type="double", backingfile=basename(p$bm$Ytime), descriptorfile=basename(p$cache$Ytime), backingpath=p$savedir )
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
            bigmemory::as.big.matrix( Pcovdata, type="double", backingfile=basename(p$bm$Pcov[[covname]]), descriptorfile=basename(p$cache$Pcov[[covname]]), backingpath=p$savedir )
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
          bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$P), descriptorfile=basename(p$cache$P), backingpath=p$savedir )
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
          bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$Pn), descriptorfile=basename(p$cache$Pn), backingpath=p$savedir )
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
          bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$Psd), descriptorfile=basename(p$cache$Psd), backingpath=p$savedir )
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
            bigmemory::as.big.matrix( Ploc, type="double", backingfile=basename(p$bm$Ploc), descriptorfile=basename(p$cache$Ploc), backingpath=p$savedir )
          }
          if (p$storage.backend == "ff" ) {
            p$ptr$Ploc = ff( Ploc, dim=dim(Ploc), file=p$cache$Ploc, overwrite=TRUE )
          }
      Ploc = DATA = NULL; gc()

      if (exists("stm_global_modelengine", p) ) {
      # create prediction suface with covariate-based additive offsets

        if (p$storage.backend == "bigmemory.ram" ) {
          tmp_P0= big.matrix( nrow=nrow(P), ncol=ncol(P) , type="double" )
          tmp_P0[] = P
          p$ptr$P0 = bigmemory::describe(tmp_P0 )
        }
        if (p$storage.backend == "bigmemory.filebacked" ) {
          p$ptr$P0  = p$cache$P0
          bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$P0), descriptorfile=basename(p$cache$P0), backingpath=p$savedir )
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
          bigmemory::as.big.matrix( P, type="double", backingfile=basename(p$bm$P0sd), descriptorfile=basename(p$cache$P0sd), backingpath=p$savedir )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$P0sd = ff( P, dim=dim(P), file=p$cache$P0sd, overwrite=TRUE )
        }

        P=NULL; gc()

        # test to see if all covars are static as this can speed up the initial predictions
        message("||| ")
        message( "||| Predicting global effect of covariates at each prediction location ... ")
        message( "||| depending upon the size of the prediction grid and number of cpus (~1hr?).. ")

        p$timec_covariates_0 =  Sys.time()
        nc_cov =NULL
        for (i in p$variables$COV ) {
          pu = stm_attach( p$storage.backend, p$ptr$Pcov[[i]] )
          nc_cov = c( nc_cov,  ncol(pu) )
        }
        p$all.covars.static = ifelse( any(nc_cov > 1),  FALSE, TRUE )
        pc = p # copy
        if (!pc$all.covars.static) if (exists("clusters.covars", pc) ) pc$clusters = pc$clusters.covars
        pc = make.list( list( tindex=1:pc$nt) , Y=pc ) # takes about 28 GB per run .. adjust cluster number temporarily
        suppressMessages( parallel.run( stm_db, p=pc, DS="global.prediction.surface" ) )
        p$time_covariates = round(difftime( Sys.time(), p$timec_covariates_0 , units="hours"), 3)
        message( paste( "||| Time taken to predict covariate surface (hours):", p$time_covariates ) )
      }

      P = NULL; gc() # yes, repeat in case covs are not modelled
    
      stm_db( p=p, DS="statistics.Sflag" )

      Y = stm_attach( p$storage.backend, p$ptr$Y )
      Yloc = stm_attach( p$storage.backend, p$ptr$Yloc )

      Yi = 1:length(Y) # index with useable data
      bad = which( !is.finite( Y[]))
      if (length(bad)> 0 ) Yi[bad] = NA

      # data locations
      bad = which( !is.finite( rowSums(Yloc[])))
      if (length(bad)> 0 ) Yi[bad] = NA

    # data locations
      if (exists("COV", p$variables)) {
        Ycov = stm_attach( p$storage.backend, p$ptr$Ycov )
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
        Ytime = stm_attach( p$storage.backend, p$ptr$Ytime )
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
          bigmemory::as.big.matrix( Yi, type="double", backingfile=basename(p$bm$Yi), descriptorfile=basename(p$cache$Yi), backingpath=p$savedir )
        }
        if (p$storage.backend == "ff" ) {
          p$ptr$Yi = ff( Yi, dim=dim(Yi), file=p$cache$Yi, overwrite=TRUE )
        }
      rm(Yi)

      if ( !exists("stm_distance_scale", p)) {
        Yloc = stm_attach( p$storage.backend, p$ptr$Yloc )
        p$stm_distance_scale = min( diff(range( Yloc[,1]) ), diff(range( Yloc[,2]) ) ) / 10
        message( paste( "||| Crude distance scale:", p$stm_distance_scale, "" ) )
      }

      if ( !exists("stm_distance_min", p)) p$stm_distance_min = mean( c(p$stm_distance_prediction, p$stm_distance_scale /20 ) )

      if ( !exists("stm_distance_max", p)) p$stm_distance_max = mean( c(p$stm_distance_prediction*10, p$stm_distance_scale * 2 ) )

      if ( !exists("sampling", p))  {
        # fractions of distance scale  to try in local block search
        p$sampling = c( 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.5, 1.75, 2 )
      }

      #browser()

      stm_db( p=p, DS="save.parameters" )  # save in case a restart is required .. mostly for the pointers to data objects
      message( "||| Finished. Moving onto analysis... ")
      p <<- p  # push to parent in case a manual restart is needed
      gc()
      
  }  # end of intialization of data structures


  # -------------------------------------
  # localized space-time modelling/interpolation/prediction
  message("||| to view maps from an external R session: ") 
  message("|||   stm(p=p, runmode='debug_pred_static_map', debug_plot_variable_index=1) ") 
  message("|||   stm(p=p, runmode='debug_pred_static_log_map', debug_plot_variable_index=1)") 
  message("|||   stm(p=p, runmode='debug_pred_dynamic_map', debug_plot_variable_index=1)") 
  message("|||   stm(p=p, runmode='debug_stats_map', debug_plot_variable_index=1)") 
  message("||| Monitor the status of modelling by looking at the output of the following file:")
  message("||| in linux, you can issue the following command:" )
  message("|||   watch -n 60 cat ",  p$stm_current_status  )

  if ( "debug_pred_static_map" == runmode) {  
      Ploc = stm_attach( p$storage.backend, p$ptr$Ploc )
      P = stm_attach( p$storage.backend, p$ptr$P )
      lattice::levelplot( (P[,debug_plot_variable_index])~Ploc[,1]+Ploc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso")
  }

  if ( "debug_pred_static_log_map" == runmode) {  
      Ploc = stm_attach( p$storage.backend, p$ptr$Ploc )
      P = stm_attach( p$storage.backend, p$ptr$P )
      lattice::levelplot( log(P[,debug_plot_variable_index])~Ploc[,1]+Ploc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso")
  }

  if ( "debug_pred_dynamic_map" == runmode) {  
      Ploc = stm_attach( p$storage.backend, p$ptr$Ploc )
      P = stm_attach( p$storage.backend, p$ptr$P )
      for (i in 1:p$nt) {
        print( lattice::levelplot( P[,i] ~ Ploc[,1] + Ploc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
      }
  }

  if ( "debug_stats_map" == runmode) {  
      Sloc = stm_attach( p$storage.backend, p$ptr$Sloc )
      S = stm_attach( p$storage.backend, p$ptr$S )
      lattice::levelplot(S[,debug_plot_variable_index]~Sloc[,1]+Sloc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE), aspect="iso")
  }

  if ( runmode %in% c("serial_debug", "debug") ) {
    currentstatus = stm_db( p=p, DS="statistics.status" )
    p = make.list( list( locs=sample( currentstatus$todo )) , Y=p ) # random order helps use all cpus
    p <<- p  # push to parent in case a manual restart is possible
    print( c( unlist( currentstatus[ c("n.total", "n.shallow", "n.todo", "n.skipped", "n.outside", "n.complete" ) ] ) ) )
    message( "Entering browser mode ...")
    browser()
    stm_interpolate (p=p )
  }


  if ( runmode %in% c("default", "stage2", "stage3" ) ) {  
    # this is the basic run
    timei1 =  Sys.time()
    currentstatus = stm_db( p=p, DS="statistics.status" )
    p = make.list( list( locs=sample( currentstatus$todo )) , Y=p ) # random order helps use all cpus
    p <<- p  # push to parent in case a manual restart is possible
    suppressMessages( parallel.run( stm_interpolate, p=p ) )
    p$time_default = round( difftime( Sys.time(), timei1, units="hours" ), 3 )
    message("||| ")
    message( paste( "||| Time taken for main stage 1, interpolations (hours):", p$time_default, "" ) )
    currentstatus = stm_db( p=p, DS="statistics.status" )
    print( c( unlist( currentstatus[ c("n.total", "n.shallow", "n.todo", "n.skipped", "n.outside", "n.complete" ) ] ) ) )
    gc()
  }
  

  if ( runmode %in% c("stage2", "stage3" ) ) {  
    timei2 =  Sys.time()
    message("||| ")
    message( "||| Starting stage 2: more permisssive distance settings (spatial extent) " )

    for ( mult in p$stm_multiplier_stage2 ) { 
      currentstatus = stm_db(p=p, DS="statistics.status.reset" ) 
      if (length(currentstatus$todo) > 0) {
        p$stm_distance_max = p$stm_distance_max * mult
        p$stm_distance_scale = p$stm_distance_scale*mult # km ... approx guess of 95% AC range 
        p = make.list( list( locs=sample( currentstatus$todo )) , Y=p ) # random order helps use all cpus
        p <<- p  # push to parent in case a manual restart is possible
        suppressMessages( parallel.run( stm_interpolate, p=p ) )
      }
    }
    p$time_stage2 = round( difftime( Sys.time(), timei2, units="hours" ), 3)
    message("||| ---")
    message( paste( "||| Time taken to stage 2 interpolations (hours):", p$time_stage2, "" ) )
    currentstatus = stm_db( p=p, DS="statistics.status" )
    print( c( unlist( currentstatus[ c("n.total", "n.shallow", "n.todo", "n.skipped", "n.outside", "n.complete" ) ] ) ) )
    gc()
  }


  if ( runmode %in% c( "stage3" ) ) {  
    timei3 =  Sys.time()
    message("||| ---")
    message( "||| Starting stage 3: simple TPS-based failsafe method to interpolate all the remaining locations " )
    toredo = stm_db( p=p, DS="flag.incomplete.predictions" )
    if ( !is.null(toredo) && length(toredo) > 0) { 
      Sflag = stm_attach( p$storage.backend, p$ptr$Sflag )
      Sflag[toredo]=0L
      p$stm_local_modelengine = "tps"  
      p = indicators.parameters( p=p, DS="bathymetry" )
      p = make.list( list( locs=sample( toredo )) , Y=p ) # random order helps use all cpus
      p <<- p  # push to parent in case a manual restart is possible
      parallel.run( stm_interpolate, p=p )
    }
    p$time_stage3 = round( difftime( Sys.time(), timei3, units="hours" ), 3)
    message( paste( "||| Time taken to stage3 interpolations (hours):", p$time_stage3, "" ) )
    currentstatus = stm_db( p=p, DS="statistics.status" )
    print( c( unlist( currentstatus[ c("n.total", "n.shallow", "n.todo", "n.skipped", "n.outside", "n.complete" ) ] ) ) )
    gc()
  }

  # save again, in case some timings/etc needed in a restart
  p <<- p  # push to parent in case a manual restart is possible

  stm_db( p=p, DS="save.parameters" )  # save in case a restart is required .. mostly for the pointers to data objects

  resp = readline( "||| Save predictions and statistics, overwriting previous results? If you are sure type <YES>:  ")
  if (resp=="YES") {
    # save solutions to disk (again .. overwrite)
    message("||| ")
    message( "||| Saving predictions to disk .. " )
    stm_db( p=p, DS="stm.prediction.redo" ) # save to disk for use outside stm*, returning to user scale

    message( "||| Saving statistics to disk .. " )
    stm_db( p=p, DS="stats.to.prediction.grid.redo") # save to disk for use outside stm*

    message ("||| Finished! ")
  }

  if ( p$storage.backend !="bigmemory.ram" ) {
    resp = readline( "||| Delete temporary files? Type to confirm <YES>:  ")
    if (resp=="YES") {
      stm_db( p=p, DS="cleanup" )
    } else {
      message("||| ")
      message( "||| Leaving temporary files alone in case you need to examine them or restart a process. ")
      message( "||| You can delete them by running: stm_db( p=p, DS='cleanup' ), once you are done. ")
    }
  }

  p$time_total = round( difftime( Sys.time(), p$time.start, units="hours" ),3)
  message("||| ")
  message( paste( "||| Time taken for ", runmode, " (hours):", p$time_total, "\n" ) )

  message( paste( "||| Your parameter 'p' has been updated in case you need to re-run something like, etc:\n" ) )
  message( paste( "||| stm(p=p, runmode='stage3') :\n" ) )
  p <<- p  # push to parent in case a manual restart is possible

  invisible()
}

