
  stmv_db = function( DS, p, yr=NULL, ret="mean", runmode=NULL, datasubset="" ) {
    #// usage: low level function to convert data into file-based data obects to permit parallel
    #// data access and manipulation and deletes/updates

    # --------------------------
    if (!exists("stmvSaveDir", p)) {
      p$stmvSaveDir = file.path( p$modeldir, p$stmv_model_label, p$project_class, paste(  p$stmv_global_modelengine, p$stmv_local_modelengine, sep="_"), p$stmv_variables$Y, p$spatial_domain)
    }

    if (DS %in% "filenames" ) {
      # input data stored as a bigmemory file to permit operations with min memory usage
      # split into separate components to reduce filelocking conflicts

      p$cache =list()
      p$cache$Y =     file.path( p$stloc, "input.Y.cache" ) # residuals of covar model or raw data if none
      p$cache$Yloc =  file.path( p$stloc, "input.Yloc.cache" )
      p$cache$Yi =    file.path( p$stloc, "input.Yi.cache" ) # index of useable data

      p$cache$P =     file.path( p$stloc, "predictions.cache" )
      p$cache$Psd =   file.path( p$stloc, "predictions_sd.cache" )
      p$cache$Pn =    file.path( p$stloc, "predictions_n.cache" )

      if ( exists( "COV", p$stmv_variables)) {
        p$cache$Ycov =  file.path( p$stloc, "input.Ycov.cache"  )
        p$cache$Pcov =  list()
        for (cov in p$stmv_variables$COV) p$cache$Pcov[[cov]] = file.path( p$stloc, paste("predictions_cov", cov, "cache", sep=".") )
      }

      if (exists( "TIME", p$stmv_variables)){
        p$cache$Ytime = file.path( p$stloc, "input.Ytime.cache" )
      }

      p$cache$Ploc =  file.path( p$stloc, "predictions_loc.cache" )

        if (p$stmv_global_modelengine !="none" ) {
          p$cache$P0 = file.path( p$stloc, "P0.cache" )
          p$cache$P0sd = file.path( p$stloc, "P0sd.cache" )
        }

      p$cache$S =     file.path( p$stloc, "statistics.cache" )
      p$cache$Sloc =  file.path( p$stloc, "statistics_loc.cache" )
      p$cache$Sflag =     file.path( p$stloc, "statistics_flag.cache" )

      p$saved_state_fn = list()
      p$saved_state_fn$P = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "mean", "rdz", sep="." ) )
      p$saved_state_fn$Pn = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "n", "rdz", sep="." ) )
      p$saved_state_fn$Psd = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "sd", "rdz", sep="." ) )
      p$saved_state_fn$stats = file.path( p$stmvSaveDir, paste( "tmp_stmv.statistics", "rdz", sep=".") )
      p$saved_state_fn$sflag = file.path( p$stmvSaveDir, paste( "tmp_stmv.sflag", "rdz", sep=".") )
        if (p$stmv_global_modelengine !="none" ) {
          p$saved_state_fn$P0 = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "mean0", "rdz", sep="." ) )
          p$saved_state_fn$P0sd = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "sd0", "rdz", sep="." ) )
        }

      if (p$storage_backend == "bigmemory.filebacked" ) {
        p$bm = p$cache
        for ( i in names(p$bm) ) {
          if ( i=="Pcov" ) {
            for (j in p$stmv_variables$COV) p$bm$Pcov[[j]] = gsub(".cache$", ".bigmemory", p$bm$Pcov[[j]] )
          } else {
            p$bm[[i]] = gsub(".cache$", ".bigmemory", p$bm[[i]] )
          }
        }
      }

      if (p$storage_backend == "bigmemory.ram" ) {
        p$bm=list() # initial storage of ram objects
      }

      return(p)
    }

    # --------------------------

    if (DS=="save.parameters")  {
      fns = file.path( p$stmvSaveDir, "p.rdz" )
      read_write_fast( p, file=fns )
      message( "||| Saved parameters to file:")
      message( fns )
    }

    # --------------------------

    if (DS=="load.parameters")  {
      fns = file.path( p$stmvSaveDir, "p.rdz" )
      if (file.exists( fns)) p = read_write_fast( fns )
      return(p)
    }

    # --------------------------
    if (DS %in% "cleanup" ) {
      # all but not the restart files
      for (fn in unlist(p$cache) ) if (length(fn)>0) if (file.exists(fn)) file.remove(fn)
      for (fn in unlist(p$bm) ) if (length(fn)>0)  if (file.exists(fn)) file.remove(fn)
      return( NULL )
    }

    # --------------------------
    if (DS %in% "cleanup.all" ) {
      # all including the restart files
      for (fn in unlist(p$cache) ) if (length(fn)>0) if (file.exists(fn)) file.remove(fn)
      for (fn in unlist(p$bm) ) if (length(fn)>0)  if (file.exists(fn)) file.remove(fn)
      for (fn in unlist( p$saved_state_fn) ) if (length(fn)>0) if (file.exists(fn)) file.remove(fn)

      return( NULL )
    }


    # -----------------


    if (DS %in% c( "boundary.redo", "boundary" ) )  {

      fn =  file.path(p$stmvSaveDir, "boundary.rdz" )
      if (DS=="boundary") {
        boundary = NULL
        if( file.exists(fn)) boundary = read_write_fast( fn)
        return( boundary )
      }

      # data:
      Y = stmv_attach(  p$storage_backend, p$ptr$Y )
      hasdata = 1:length(Y)
      bad = which( !is.finite( Y[]))
      if (length(bad)> 0 ) hasdata[bad] = NA

      # covariates (independent vars)
      if ( exists( "COV", p$stmv_variables) ) {
        if (length(p$stmv_variables$COV) > 0) {
          Ycov = stmv_attach(  p$storage_backend, p$ptr$Ycov )
          if ( length( p$stmv_variables$COV ) == 1 ) {
            bad = which( !is.finite( Ycov[]) )
          } else {
            bad = which( !is.finite( rowSums(Ycov[])) )
          }
          if (length(bad)> 0 ) hasdata[bad] = NA
        }
      }

      ii = na.omit(hasdata)
      Yloc = stmv_attach(  p$storage_backend, p$ptr$Yloc )
      yplon = trunc( ( Yloc[ii,1] - p$origin[1] )/p$pres) + 1L
      yplat = trunc( ( Yloc[ii,2] - p$origin[2] )/p$pres) + 1L
      uu = unique( array_map( "2->1", cbind(yplon, yplat), c(p$nplons, p$nplats) ) )
      vv = array_map( "1->2", uu, c(p$nplons, p$nplats) )

      ww = cbind( (vv[,1] - 1) * p$pres + p$origin[1], (vv[,2] - 1) * p$pres + p$origin[2] )

      if (!exists("stmv_nonconvexhull_alpha", p)) p$stmv_nonconvexhull_alpha=20
      boundary=list( polygon = non_convex_hull( ww, lengthscale=p$stmv_nonconvexhull_alpha, plot=FALSE ) )

      # statistical output locations
      Sloc = stmv_attach(  p$storage_backend, p$ptr$Sloc )
      boundary$inside.polygon = point.in.polygon( Sloc[,1], Sloc[,2],
          boundary$polygon[,1], boundary$polygon[,2], mode.checked=TRUE )

      read_write_fast( boundary, file=fn )
      plot( Yloc[ii,], pch=".", col="grey" ) # data locations
      points( Sloc[which(boundary$inside.polygon==1),], pch=".", col="orange" )
      lines( boundary$polygon[] , col="green", pch=2 )
      message( "||| Check the map of data and boundaries. ")
      message( "||| If not suitable, set another value for p$stmv_nonconvexhull_alpha value (radius; distance) ")
      message( "||| and re-run stmv() " )
      return( fn )
    }

    # -----


    if (DS %in% c("stmv.results", "stmv.prediction", "stmv.stats") )  {

      if (DS=="stmv.prediction") {
        if (! exists("TIME", p$stmv_variables)) {
          fn = file.path( p$stmvSaveDir, paste("stmv.prediction",  ret, "rdz", sep="." ) )
        } else {
          fn = file.path( p$stmvSaveDir, paste("stmv.prediction",  ret, yr, "rdz", sep="." ) )
        }
        if (file.exists(fn) ) out = read_write_fast(fn)
        return(out)
      }

      if (DS=="stmv.stats") {
        fn = file.path( p$stmvSaveDir, paste( "stmv.statistics", "rdz", sep=".") )
        stats = NULL
        if (file.exists(fn)) stats=read_write_fast(fn)
        return(stats)
      }

      shallower = NULL
      if ( exists("stmv_filter_depth_m", p) && is.finite( p$stmv_filter_depth_m) ) {
        if ( "z" %in% p$stmv_variables$COV ){
          depths = stmv_attach( p$storage_backend, p$ptr$Pcov[["z"]] )[]
          ii = which( depths[] < p$stmv_filter_depth_m )
          if (length(ii) > 0) shallower = ii
          ii= depths=NULL
        }
      }

      if ( exists("TIME", p$stmv_variables)) {

        p0 = p
        p$clusters = rep( "localhost", parallel::detectCores() )
        
        parallel_run(
          p=p,
          shallower=shallower,
          runindex=list( pny=1:p$ny ),
          FUNC= function( ip=NULL, p, shallower ) {
            if (exists( "libs", p)) RLibrary( p$libs )
            if (is.null(ip)) ip = 1:p$nruns
            PP = stmv_attach( p$storage_backend, p$ptr$P )
            PPsd = stmv_attach( p$storage_backend, p$ptr$Psd )
              if (p$stmv_global_modelengine !="none" ) {
                P0 = stmv_attach( p$storage_backend, p$ptr$P0 )
                P0sd = stmv_attach( p$storage_backend, p$ptr$P0sd )
              }
            vv = ncol(PP)
            for (it in ip) {
              y = p$yrs[ p$runs[it, "pny"] ]

              if ( vv > p$ny ) {
                ww = (it-1) * p$nw + (1:p$nw)
                P = PP  [,ww]
                V = PPsd[,ww] # simpleadditive independent errors assumed
              } else if ( vv==p$ny) {
                P = PP[,it]
                V = PPsd[,it]
              }

                if (p$stmv_global_modelengine !="none" ) {
                  ## maybe add via simulation, note: P0 and P are on link scale to this point
                  uu = which(!is.finite(P[]))
                  if (length(uu)>0) P[uu] = 0 # permit global predictions to pass through ..
                  P = P[] + P0[,it]
                  nV = which(!is.finite(V[]))
                  if (length(nV)>0) V[nV] = 0 # permit covariate-base predictions to pass through ..
                  V = sqrt( V[]^2 + P0sd[,it]^2) # simple additive independent errors assumed
                }

              if ( !is.null(shallower) ){
                if ( is.vector(P) ) {
                  P[shallower] = NA
                  V[shallower] = NA
                } else {
                  P[shallower,] = NA
                  V[shallower,] = NA
                }
              }

              Pl =  P[] - 1.96* V[]
              Pu =  P[] + 1.96* V[]
              P =   P[]

              # return to user scale (that of Y)
              if (p$stmv_global_modelengine !="none" ) {
                if (exists("linkinv", p$stmv_global_family )) {
                  Pl = p$stmv_global_family$linkinv( Pl[] )
                  Pu = p$stmv_global_family$linkinv( Pu[] )
                  P = p$stmv_global_family$linkinv( P[] )
                }
              }

              # any additional transformations
              if (exists("stmv_Y_transform", p)) {
                Pl = p$stmv_Y_transform$invers (Pl[])  # p$stmv_Y_transform[2] is the inverse transform
                Pu = p$stmv_Y_transform$invers (Pu[])
                P = p$stmv_Y_transform$invers (P[])
              }

              read_write_fast( P,  file=file.path( p$stmvSaveDir, paste("stmv.prediction", "mean", y, "rdz", sep="." ) ) )
              read_write_fast( Pl, file=file.path( p$stmvSaveDir, paste("stmv.prediction", "lb",   y, "rdz", sep="." ) ) )
              read_write_fast( Pu, file=file.path( p$stmvSaveDir, paste("stmv.prediction", "ub",   y, "rdz", sep="." ) ) )
              # print ( paste("Year:", y)  )
            }
          }
        )

      } else {
        # serial run only ...

        PP = stmv_attach( p$storage_backend, p$ptr$P )
        PPsd = stmv_attach( p$storage_backend, p$ptr$Psd )
          if (p$stmv_global_modelengine !="none" ) {
            P0 = stmv_attach( p$storage_backend, p$ptr$P0 )
            P0sd = stmv_attach( p$storage_backend, p$ptr$P0sd )
          }

        P = PP[]
        V = PPsd[]
          if (p$stmv_global_modelengine !="none" ) {
            uu = which(!is.finite(P[]))
            if (length(uu)>0) P[uu] = 0 # permit covariate-base predictions to pass through ..
            P = P[] + P0[]  # both on link scale
            nV = which(!is.finite(V[]))
            if (length(nV)>0) V[nV] = 0 # permit covariate-base predictions to pass through ..
            V = sqrt( V[]^2 + P0sd[]^2) # simple additive independent errors assumed
          }
        if ( !is.null(shallower) ){
          P[shallower] = NA
          V[shallower] = NA
        }

        Pl =  P[] - 1.96* V[]
        Pu =  P[] + 1.96* V[]
        P =  P[]

        # return to user scale (that of Y)
        if (exists("linkinv", p$stmv_global_family )) {
          Pl = p$stmv_global_family$linkinv( Pl[] )
          Pu = p$stmv_global_family$linkinv( Pu[] )
          P = p$stmv_global_family$linkinv( P[] )
        }

        if (exists("stmv_Y_transform", p)) {
          Pl = p$stmv_Y_transform$invers (Pl[])  # p$stmv_Y_transform[2] is the inverse transform
          Pu = p$stmv_Y_transform$invers (Pu[])
          P = p$stmv_Y_transform$invers (P[])
        }

        read_write_fast( P,  file=file.path( p$stmvSaveDir, paste("stmv.prediction", "mean", "rdz", sep="." ) ) )
        read_write_fast( Pl, file=file.path( p$stmvSaveDir, paste("stmv.prediction", "lb",   "rdz", sep="." ) ) )
        read_write_fast( Pu, file=file.path( p$stmvSaveDir, paste("stmv.prediction", "ub",   "rdz", sep="." ) ) )
      } # end if TIME

      message( "\n||| Saving predictions complete: ", format(Sys.time()),  "\n" )

      # prediction.stats .. warp to same resolution

      Ploc = stmv_attach( p$storage_backend, p$ptr$Ploc )
      S = stmv_attach( p$storage_backend, p$ptr$S )
      Sloc = stmv_attach( p$storage_backend, p$ptr$Sloc )

      nPloc = nrow(Ploc)
      nSloc = nrow(Sloc)

      if (nPloc == nSloc) {
        stats = S[]
        names(stats) = p$statsvars
        # nothing else to do as the dim of S and P are the same..
      }  else {
        # system size: nr = nx .. x/plon ..nr = p$nplons;  nc = ny .. y/plat .. nc = p$nplats
        u = list(
          x=seq( p$corners$plon[1], p$corners$plon[2], by=p$stmv_distance_statsgrid ),
          y=seq( p$corners$plat[1], p$corners$plat[2], by=p$stmv_distance_statsgrid )
        )
        u$z = matrix( NA, nrow=length(u$x), ncol=length(u$y) )
        stats = matrix( NaN, ncol=length( p$statsvars ), nrow=nrow( Ploc) )  # output data .. ff does not handle NA's .. using NaN for now
        colnames(stats) = p$statsvars
        Sres=c(p$stmv_distance_statsgrid, p$stmv_distance_statsgrid)
        Sind = as.matrix(array_map( "xy->2", coords=Sloc[], origin=p$origin, res=Sres ))  # map Stats Locs to Plocs
        for ( i in 1:length( p$statsvars ) ) {
          # print(i)
          u$z[] = NA # reset
          u$z[Sind] = S[,i]
          stats[,i] = as.vector( fields::interp.surface( u, loc=Ploc[] ) ) # linear interpolation
          if (all(!is.finite(stats[,i]))) next()
          rY = range( S[,i], na.rm=TRUE )
          lb = which( stats[,i] < rY[1] )
          if (length(lb) > 0) stats[lb,i] = rY[1]
          lb = NULL
          ub = which( stats[,i] > rY[2] )
          if (length(ub) > 0) stats[ub,i] = rY[2]
          ub = NULL
        }
        # lattice::levelplot( stats[,1] ~ Ploc[,1]+Ploc[,2])
        boundary = try( stmv_db( p=p, DS="boundary" ) )
        if( !("try-error" %in% class(boundary) ) ) {
          if (!is.null(boundary)) {
            inside.polygon = point.in.polygon( Ploc[,1], Ploc[,2],
            boundary$polygon[,1], boundary$polygon[,2], mode.checked=TRUE )
            o = which( inside.polygon == 0 ) # outside boundary
            if (length(o) > 0) stats[o,] = NA
          }
        }
        if (length(shallower)>0) stats[shallower,] = NA
      }

      fn = file.path( p$stmvSaveDir, paste( "stmv.statistics", "rdz", sep=".") )
      read_write_fast( stats, file=fn )
      message( "\n||| Saving statistics complete: ", format(Sys.time()),  "\n" )

      return( NULL)

      if (0){
        #         p$statsvars
        # [1] "sdTotal"    "rsquared"   "ndata"      "sdSpatial"  "sdObs"      "phi"        "nu"         "localrange"
        i = 4
        lattice::levelplot( stats[,i] ~ Ploc[,1]+Ploc[,2], aspect="iso")
        # ii = which (is.finite(stats[,i]))
        # lattice::levelplot( stats[ii,i] ~ Ploc[ii,1]+Ploc[ii,2])
      }

      if(0) {
        i = 1
        Ploc = stmv_attach( p$storage_backend, p$ptr$Ploc )
        Z = smooth.2d( Y=P[], x=Ploc[], ncol=p$nplats, nrow=p$nplons, cov.function=stationary.cov, Covariance="Matern", range=p$stmv_lowpass_phi, nu=p$stmv_lowpass_nu )
        dev.new(); image(Z)
        dev.new();lattice::levelplot( P[] ~ Ploc[,1] + Ploc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
      }

    }


    # ----------------


    if (DS %in% c("save_current_state") ) {

      if ("statistics" %in% datasubset ) datasubset = unique( c( datasubset, "S", "Sflag" ) )

      # named differently to avoid collisions
      if ( "P" %in% datasubset ) {
        sP = stmv_attach( p$storage_backend, p$ptr$P )[]
        read_write_fast( sP, file=paste(p$saved_state_fn$P, runmode, sep=".") )
        sP = NULL

      }
      if ( "P0" %in% datasubset ) {
        if (p$stmv_global_modelengine !="none" ) {
          sP0 = stmv_attach( p$storage_backend, p$ptr$P0 )[]
          read_write_fast( sP0,   file=paste(p$saved_state_fn$P0, runmode, sep=".")  )
          sP0 = NULL
        }
      }

      if ( "Psd" %in% datasubset ) {
        sPsd = stmv_attach( p$storage_backend, p$ptr$Psd )[]
        read_write_fast( sPsd, file=paste(p$saved_state_fn$Psd, runmode, sep=".") )
        sPsd = NULL
      }

      if ( "P0sd" %in% datasubset ) {
        if (p$stmv_global_modelengine !="none" ) {
            sP0sd = stmv_attach( p$storage_backend, p$ptr$P0sd )[]
            read_write_fast( sP0sd, file=paste(p$saved_state_fn$P0sd, runmode, sep=".")  )
            sP0sd = NULL
          }
      }

      if ( "Pn" %in% datasubset ) {
        sPn = stmv_attach( p$storage_backend, p$ptr$Pn )[]
        read_write_fast( sPn, file=paste(p$saved_state_fn$Pn, runmode, sep=".") )
        sPn = NULL
      }

      if ( "S" %in% datasubset ) {
        sS = stmv_attach( p$storage_backend, p$ptr$S )[]
        read_write_fast( sS, file=paste(p$saved_state_fn$stats, runmode, sep=".") )
        sS = NULL
      }

      if ( "Sflag" %in% datasubset ) {
        sSflag = stmv_attach( p$storage_backend, p$ptr$Sflag )[]
        read_write_fast( sSflag, file=paste(p$saved_state_fn$sflag, runmode, sep=".") )
        sSflag = NULL
      }
      gc()
      return(NULL)
    }


    # =--------------------

    if (DS %in% c("load_saved_state") ) {
      returnflag = TRUE

      if ("statistics" %in% datasubset ) datasubset = unique( c( datasubset, "S", "Sflag" ) )

      # named differently to avoid collisions
      if ( "P" %in% datasubset ) {
        P = stmv_attach( p$storage_backend, p$ptr$P )
        sP = matrix( NaN, nrow=nrow(P), ncol=ncol(P) )
        if (file.exists(paste( p$saved_state_fn$P, runmode, sep="."))) {
          P = read_write_fast( paste( p$saved_state_fn$P, runmode, sep=".") )
        } else {
          returnflag = FALSE
        }
        if (is.vector(sP))  sP=as.matrix(sP, nrow=nrow(P), ncol=1) # big matrix does not like vectors
        P[] = sP[]
        sP = NULL
      }

      if ( "P0" %in% datasubset ) {
        if (p$stmv_global_modelengine !="none" ) {
          P0 = stmv_attach( p$storage_backend, p$ptr$P0 )
          sP0 = matrix( NaN, nrow=nrow(P0), ncol=ncol(P0) )
          if (file.exists(paste( p$saved_state_fn$P0, runmode, sep="."))) {
            P0 = read_write_fast( paste( p$saved_state_fn$P0, runmode, sep=".") )
          } else {
            returnflag = FALSE
          }
          if (is.vector(sP0))  sP0=as.matrix(sP0, nrow=nrow(P0), ncol=1) # big matrix does not like vectors
          P0[] = sP0[]
          sP0 = NULL
        }
      }

      if ( "Psd" %in% datasubset ) {
        Psd = stmv_attach( p$storage_backend, p$ptr$Psd )
        sPsd = matrix( NaN, nrow=nrow(Psd), ncol=ncol(Psd) )
        if (file.exists(paste( p$saved_state_fn$Psd, runmode, sep="."))) {
          Psd = read_write_fast( paste( p$saved_state_fn$Psd, runmode, sep=".") )
        } else {
          returnflag = FALSE
        }
        if (is.vector(sPsd)) sPsd=as.matrix(sPsd, nrow=nrow(Psd), ncol=1)  # big matrix does not like vectors
        Psd[] = sPsd[]
        sPsd = NULL
      }

      if ( "P0sd" %in% datasubset ) {
        if (p$stmv_global_modelengine !="none" ) {
          P0sd = stmv_attach( p$storage_backend, p$ptr$P0sd )
          sP0sd = matrix( NaN, nrow=nrow(P0sd), ncol=ncol(P0sd) )
          if (file.exists(paste( p$saved_state_fn$P0sd, runmode, sep="."))) {
            P0sd = read_write_fast( paste( p$saved_state_fn$P0sd, runmode, sep=".") )
          } else {
            returnflag = FALSE
          }
          if (is.vector(sP0sd)) sP0sd=as.matrix(sP0sd, nrow=nrow(P0sd), ncol=1)  # big matrix does not like vectors
          P0sd[] = sP0sd[]
          sP0sd = NULL
        }
      }

      if ( "Pn" %in% datasubset ) {
        Pn = stmv_attach( p$storage_backend, p$ptr$Pn )
        sPn = matrix( NaN, nrow=nrow(Pn), ncol=ncol(Pn) )
        if (file.exists(paste( p$saved_state_fn$Pn, runmode, sep="."))) {
          Pn = read_write_fast( paste( p$saved_state_fn$Pn, runmode, sep=".") )
        } else {
          returnflag = FALSE
        }
        if (is.vector(sPn)) sPn=as.matrix(sPn, nrow=nrow(Pn), ncol=1)  # big matrix does not like vectors
        Pn[] = sPn[]
        sPn = NULL
      }


      if ( "S" %in% datasubset ) {
        S = stmv_attach( p$storage_backend, p$ptr$S )
        sS = matrix( NaN, nrow=nrow(S), ncol=ncol(S) )
        if (file.exists(paste( p$saved_state_fn$stats, runmode, sep="."))) {
          stats = read_write_fast( paste( p$saved_state_fn$stats, runmode, sep=".") )
        } else {
          fn = file.path( p$stmvSaveDir, paste( "stmv.statistics", "rdz", sep=".") )
          if (!file.exists(fn)) stop( "stmv.stats not found")
          stats = NULL
          stats = read_write_fast(fn)
          if (is.null(stats)) stop ("stmv.stats empty")
          Sloc = stmv_attach( p$storage_backend, p$ptr$Sloc )
          Ploc = stmv_attach( p$storage_backend, p$ptr$Ploc )
          nx = length(seq( p$corners$plon[1], p$corners$plon[2], by=p$stmv_distance_statsgrid ))
          ny = length(seq( p$corners$plat[1], p$corners$plat[2], by=p$stmv_distance_statsgrid ) )
          if (nx*ny != nrow(S) ) stop( "stmv.statistics has the wrong dimensionality/size" )

          # the saved state is a dense interpolated surface.
          # Linear interpolation to get down scaling is sufficient
          # .. could use fft, but no real value as it is already interpolated
          for ( i in 1:length( p$statsvars ) ) {
            u = as.image( stats[,i], x=Ploc[,], nx=nx, ny=ny )
            S[,i] = as.vector( fields::interp.surface( u, loc=Sloc[] ) ) # linear interpolation
          }
          nx = ny = u = stats = NULL
          returnflag = FALSE
        }
        if (is.vector(sS)) sS=as.matrix(sS, nrow=nrow(S), ncol=1)  # big matrix does not like vectors
        S[] = sS[]
        sS = NULL
      }

      if ( "Sflag" %in% datasubset ) {
        Sflag = stmv_attach( p$storage_backend, p$ptr$Sflag )
        sSflag = matrix( NaN, nrow=nrow(Sflag), ncol=ncol(Sflag) )
        if (file.exists(paste( p$saved_state_fn$sflag, runmode, sep="."))) {
          sflag = read_write_fast( paste( p$saved_state_fn$sflag, runmode, sep=".") )
        } else {
          returnflag = FALSE
        }
        if (is.vector(sSflag)) sSflag=as.matrix(sSflag, nrow=nrow(Sflag), ncol=1)  # big matrix does not like vectors
        Sflag[] = sSflag[]
        sSflag = NULL
        currentstatus = stmv_statistics_status( p=p, reset=c("insufficient_data", "variogram_failure", "variogram_range_limit", "unknown" ) )
      }

      return( returnflag )
    }


  }
