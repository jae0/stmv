
  stmv_db = function( DS, p, B=NULL, yr=NULL, ret="mean" ) {
    #// usage: low level function to convert data into file-based data obects to permit parallel
    #// data access and manipulation and deletes/updates
    #// B is the xyz or xytz data or the function to get the data to work upon

    # --------------------------
    if (!exists("stmvSaveDir", p)) {
      p$stmvSaveDir = file.path(p$data_root, "modelled", p$variables$Y, p$spatial.domain )
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

      if ( exists( "COV", p$variables)) {
        p$cache$Ycov =  file.path( p$stloc, "input.Ycov.cache"  )
        p$cache$Pcov =  list()
        for (cov in p$variables$COV) p$cache$Pcov[[cov]] = file.path( p$stloc, paste("predictions_cov", cov, "cache", sep=".") )
      }

      if (exists( "TIME", p$variables)){
        p$cache$Ytime = file.path( p$stloc, "input.Ytime.cache" )
      }

      p$cache$Ploc =  file.path( p$stloc, "predictions_loc.cache" )

      if (exists("stmv_global_modelengine", p) ) {
        if (p$stmv_global_modelengine !="none" ) {
          p$cache$P0 = file.path( p$stloc, "P0.cache" )
          p$cache$P0sd = file.path( p$stloc, "P0sd.cache" )
        }
      }

      p$cache$S =     file.path( p$stloc, "statistics.cache" )
      p$cache$Sloc =  file.path( p$stloc, "statistics_loc.cache" )
      p$cache$Sflag =     file.path( p$stloc, "statistics_flag.cache" )


      p$saved_state_fn = list()
      p$saved_state_fn$P = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "mean", "rdata", sep="." ) )
      p$saved_state_fn$Pn = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "n", "rdata", sep="." ) )
      p$saved_state_fn$Psd = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "sd", "rdata", sep="." ) )
      p$saved_state_fn$stats = file.path( p$stmvSaveDir, paste( "tmp_stmv.statistics", "rdata", sep=".") )
      p$saved_state_fn$sflag = file.path( p$stmvSaveDir, paste( "tmp_stmv.sflag", "rdata", sep=".") )
      if (exists("stmv_global_modelengine", p) ) {
        if (p$stmv_global_modelengine !="none" ) {
          p$saved_state_fn$P0 = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "mean0", "rdata", sep="." ) )
          p$saved_state_fn$P0sd = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "sd0", "rdata", sep="." ) )
        }
      }

      if (p$storage.backend == "bigmemory.filebacked" ) {
        p$bm = p$cache
        for ( i in names(p$bm) ) {
          if ( i=="Pcov" ) {
            for (j in p$variables$COV) p$bm$Pcov[[j]] = gsub(".cache$", ".bigmemory", p$bm$Pcov[[j]] )
          } else {
            p$bm[[i]] = gsub(".cache$", ".bigmemory", p$bm[[i]] )
          }
        }
      }

      if (p$storage.backend == "bigmemory.ram" ) {
        p$bm=list() # initial storage of ram objects
      }

      return(p)
    }

    # --------------------------

    if (DS=="save.parameters")  {
      fns = file.path( p$stmvSaveDir, "p.rdata" )
      save( p, file=fns, compress=TRUE )
      message( "||| Saved parameters to file:")
      message( fns )
    }

    # --------------------------

    if (DS=="load.parameters")  {
      fns = file.path( p$stmvSaveDir, "p.rdata" )
      if (file.exists( fns)) load( fns )
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

    if ( DS %in% c( "statistics.status", "statistics.status.reset") ) {

      E = stmv_error_codes()
      Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )
      if ( DS=="statistics.status.reset" ) {
        # to reset all rejected locations
        toreset = which( Sflag[] >  E[["outside_bounds"]])
        if (length(toreset) > 0) {
          Sflag[toreset] = E[["todo"]] # to reset all the problem flags to todo
        }
      }
      out = list()

      out$todo = which( Sflag[]==E[["todo"]] )       # 0 = TODO
      out$done = which( Sflag[]==E[["complete"]] )       # 1 = completed
      out$outside = which( Sflag[]==E[["outside_bounds"]] )    # 2 = oustide bounds(if any)
      out$shallow = which( Sflag[]==E[["too_shallow"]] )    # 3 = depth shallower than p$depth.filter (if it exists .. z is a covariate)
      out$predareaerror = which( Sflag[]==E[["prediction_area"]] ) # 4=predictionarea not ok,
      out$nodata = which( Sflag[]==E[["insufficient_data"]] )     # 5=skipped due to insufficient data,
      out$variogramerror = which( Sflag[]==E[["variogram_failure"]] ) # 6=skipped .. fast variogram did not work
      out$vrangeerror = which( Sflag[]==E[["variogram_range_limit"]] )     # 7=variogram estimated range not ok
      out$predictionerror = which( Sflag[]==E[["prediction_error"]] )     # 8=problem with prediction and/or modelling
      out$predictionupdateerror = which( Sflag[]==E[["prediction_update_error"]] )     # 8=problem with prediction and/or modelling
      out$statisticsupdateerror = which( Sflag[]==E[["statistics_update_error"]] ) # 4=predictionarea not ok,
      out$skipped = which( Sflag[] == E[["unknown"]] )   # 9 not completed due to a failed attempt
    #
      # do some counts
      out$n.todo = length(out$todo)
      out$n.complete = length(out$done)
      out$n.outside = length(which(is.finite(out$outside)))
      out$n.shallow = length(out$shallow)
      out$n.predareaerror = length(out$predareaerror)
      out$n.nodata = length(out$nodata)
      out$n.variogramerror = length(out$variogramerror)
      out$n.vrangeerror = length(out$vrangeerror)
      out$n.predictionerror = length(out$predictionerror)
      out$n.predictionupdateerror = length(out$predictionupdateerror)
      out$n.statisticsupdateerror = length(out$statisticsupdateerror)
      out$n.skipped = length(out$skipped)
      out$n.total = length(Sflag)

      out$prop_incomp = round( out$n.todo / ( out$n.total - out$n.outside ), 3)
      message( paste("||| Proportion to do:", out$prop_incomp, "\n" ))
      return( out )

      if (0) {
        Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )
        Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )

        plot( Yloc[], pch=".", col="grey" ) # data locations
        bnds = try( stmv_db( p=p, DS="boundary" ) )
        if ( !is.null(bnds)) {
          lines( bnds$polygon[] , col="green", pch=2 )
          points( Sloc[which(bnds$inside.polygon==1),], pch=".", col="orange", cex=5 )
        }
        points( Sloc[which( Sflag[]== E[["todo"]]),], pch=".", col="blue", cex=5 )
        points( Sloc[which( Sflag[]== E[["complete"]]),], pch=".", col="purple", cex=5 )
        points( Sloc[which( Sflag[]== E[["outside_bounds"]]),], pch=".", col="red", cex=5 )
        points( Sloc[which( Sflag[]== E[["too_shallow"]]),], pch=".", col="yellow", cex=5 )
        points( Sloc[which( Sflag[]== E[["prediction_area"]]),], pch=".", col="green", cex=5 )
        points( Sloc[which( Sflag[]== E[["insufficient_data"]]),], pch=".", col="red2", cex=5 )
        points( Sloc[which( Sflag[]== E[["variogram_failure"]]),], pch=".", col="yellow2", cex=5 )
        points( Sloc[which( Sflag[]== E[["variogram_range_limit"]]),], pch=".", col="green2", cex=5 )
        points( Sloc[which( Sflag[]== E[["prediction_error"]]),], pch=".", col="green3", cex=5 )
        points( Sloc[which( Sflag[]== E[["prediction_update_error"]]),], pch=".", col="green3", cex=5 )
        points( Sloc[which( Sflag[]== E[["statistics_update_error"]]),], pch=".", col="green3", cex=5 )
        points( Sloc[which( Sflag[]== E[["unknown"]]),], pch=".", col="magenta", cex=5 )
      }
    }


    #-------------------


    if ( DS %in% c( "statistics.Sflag" ) ) {
      # create location specific flags for analysis, etc..

      # flag areas overlapping with prediction locations:
      Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
      Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )

      pidP = array_map( "xy->1", Ploc, gridparams=p$gridparams )
      pidS = array_map( "xy->1", Sloc, gridparams=p$gridparams )
      overlap = match( pidS, pidP )
      outside = which( !is.finite( overlap ))

      E = stmv_error_codes()
      Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )
      if (length(outside)  > 0 ) Sflag[outside] = E[["outside_bounds"]]  # outside of prediction domain

      # catch data boundaries if present
      if (exists( "boundary", p) && p$boundary) {
        timeb0 =  Sys.time()
        message("\n")
        message( "||| Defining boundary polygon for data .. this reduces the number of points to analyse")
        message( "||| but takes a few minutes to set up ...")
        stmv_db( p=p, DS="boundary.redo" ) # ~ 5 min on nfs
      # last set of filters to reduce problem size
        Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )
        bnds = try( stmv_db( p=p, DS="boundary" ) )
        if (!is.null(bnds)) {
          if( !("try-error" %in% class(bnds) ) ) {
            outside = which( bnds$inside.polygon == 0 ) # outside boundary
            if (length(outside)>0) Sflag[outside] = E[["outside_bounds"]]
        }}
        bnds = NULL
        message( paste( "||| Time taken to estimate spatial bounds (mins):", round( difftime( Sys.time(), timeb0, units="mins" ),3) ) )
      }

      if ( exists("depth.filter", p) && is.finite(p$depth.filter) ) {
        # additionaldepth-based filter:
        # assuming that there is depth information in Pcov, match Sloc's and filter out locations that fall on land
        if ( "z" %in% p$variables$COV ){
          z = stmv_attach( p$storage.backend, p$ptr$Pcov[["z"]] )[]
          Pabove = which( z < p$depth.filter ) # negative = above land
          Pbelow = which( z >= p$depth.filter )

          Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
          Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )

          pidA = array_map( "xy->1", Ploc[Pabove,], gridparams=p$gridparams )
          pidB = array_map( "xy->1", Ploc[Pbelow,], gridparams=p$gridparams )
          sid  = array_map( "xy->1", Sloc[], gridparams=p$gridparams )

          below = which( is.finite( match( sid, pidB ) ))
          above = which( is.finite( match( sid, pidA ) ))

          Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )
          if (length(below) > 0 ) Sflag[below] = E[["todo"]]
          if (length(above) > 0 ) Sflag[above] = E[["too_shallow"]]

          if (0) {
            Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )
            plot( Yloc[], pch=".", col="grey" ) # data locations
            bnds = try( stmv_db( p=p, DS="boundary" ) )
            if (!is.null(bnds)) {
              if ( !("try-error" %in% class(bnds) ) ) {
                points( Sloc[which(bnds$inside.polygon==1),], pch=".", col="orange" )
                lines( bnds$polygon[] , col="green", pch=2 )
              }
            }
            points( Sloc[which( Sflag[]==E[["todo"]]),], pch=".", col="blue" )
          }
        }
      }

    }


    #---------------------


    if (DS== "flag.incomplete.predictions") {
      # statistics locations where estimations need to be redone
      E = stmv_error_codes()
      P = stmv_attach( p$storage.backend, p$ptr$P )
      if (ncol(P) == 1 ) {
        noP = which( !is.finite( P[]) )
      } else {
        noP = which( !is.finite( rowSums( P[])) )
      }
      uP = NULL
      if( length(noP)>0 ) {
        Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
        Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )

        Sloc_nplat = ceiling( diff( p$corners$plat) / p$stmv_distance_statsgrid)
        Sloc_nplon = ceiling( diff( p$corners$plon) / p$stmv_distance_statsgrid)

        Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
        uS = array_map( "2->1", round( cbind(Sloc[,1]-p$origin[1], Sloc[,2]-p$origin[2])/p$stmv_distance_statsgrid)+1, c(Sloc_nplon, Sloc_nplat) )
        uP = array_map( "2->1", round( cbind(Ploc[noP,1]-p$origin[1], Ploc[noP,2]-p$origin[2])/p$stmv_distance_statsgrid)+1, c(Sloc_nplon, Sloc_nplat) )
        inrange = which( (uP >= min(uS)) & (uP <= max(uS)) )
        if (length( inrange) > 0) uP = uP[inrange]
        uP = unique(uP)
        Sflag[uP] = E[["todo"]]  # force set to redo
      }

      return(uP)
    }

    #---------------------

    if (DS %in% c( "boundary.redo", "boundary" ) )  {

      fn =  file.path(p$stmvSaveDir, "boundary.rdata" )
      if (DS=="boundary") {
        boundary = NULL
        if( file.exists(fn)) load( fn)
        return( boundary )
      }

      # data:
      Y = stmv_attach(  p$storage.backend, p$ptr$Y )
      hasdata = 1:length(Y)
      bad = which( !is.finite( Y[]))
      if (length(bad)> 0 ) hasdata[bad] = NA

      # covariates (independent vars)
      if ( exists( "COV", p$variables) ) {
        if (length(p$variables$COV) > 0) {
          Ycov = stmv_attach(  p$storage.backend, p$ptr$Ycov )
          if ( length( p$variables$COV ) == 1 ) {
            bad = which( !is.finite( Ycov[]) )
          } else {
            bad = which( !is.finite( rowSums(Ycov[])) )
          }
          if (length(bad)> 0 ) hasdata[bad] = NA
        }
      }

      ii = na.omit(hasdata)
      Yloc = stmv_attach(  p$storage.backend, p$ptr$Yloc )
      yplon = round( ( Yloc[ii,1] - p$origin[1] )/p$pres) + 1
      yplat = round( ( Yloc[ii,2] - p$origin[2] )/p$pres) + 1
      uu = unique( array_map( "2->1", cbind(yplon, yplat), c(p$nplons, p$nplats) ) )
      vv = array_map( "1->2", uu, c(p$nplons, p$nplats) )

      ww = cbind( (vv[,1] - 1) * p$pres + p$origin[1], (vv[,2] - 1) * p$pres + p$origin[2] )

      if (!exists("stmv_nonconvexhull_alpha", p)) p$stmv_nonconvexhull_alpha=20
      boundary=list( polygon = non_convex_hull( ww, alpha=p$stmv_nonconvexhull_alpha, plot=FALSE ) )

      # statistical output locations
      Sloc = stmv_attach(  p$storage.backend, p$ptr$Sloc )
      boundary$inside.polygon = point.in.polygon( Sloc[,1], Sloc[,2],
          boundary$polygon[,1], boundary$polygon[,2], mode.checked=TRUE )

      save( boundary, file=fn, compress=TRUE )
      plot( Yloc[ii,], pch=".", col="grey" ) # data locations
      points( Sloc[which(boundary$inside.polygon==1),], pch=".", col="orange" )
      lines( boundary$polygon[] , col="green", pch=2 )
      message( "||| Check the map of data and boundaries. ")
      message( "||| If not suitable, set another value for p$stmv_nonconvexhull_alpha value (radius; distance) ")
      message( "||| and re-run stmv() " )
      return( fn )
    }

    # -----

    if (DS %in% c("global_model", "global_model.redo") ) {

      fn.global_model = file.path( p$stmvSaveDir, paste( "global_model", p$stmv_global_modelengine, "rdata", sep=".") )

      if (DS =="global_model") {
        global_model = NULL
        if (file.exists( fn.global_model ))  load(fn.global_model)
        return(global_model)
      }

      if ( file.exists( fn.global_model ) ) {
        message( "||| A global model already exists. It will be overwritten in 10 seconds.")
        message( "|||   Type <ctrl-c> or <esc> to interrupt. To reuse the saved model ")
        message( "|||   leave out 'global_model' as a runmode option ... overwriting in:")
        for (i in 9:0) {
          Sys.sleep(1)
          cat( i)
          cat(" ")
        }
      }

      if ( length( p$variables$COV ) > 0 ) {
        good = which( is.finite (rowSums(B[ , c(p$variables$Y,p$variables$COV) ])) )
      } else {
        good = which( is.finite (B[,p$variables$Y ] ) )
      }

      ngood = length(good)
      if ( ngood > 0 ) {
        if ( ngood < nrow(B) ) {
          B = B[good,]
        }
      }


      if ( p$stmv_global_modelengine == "userdefined" )  {
        # need a file with predictions a function to get/return them
        if (!(exists("run", p$stmv_global_model ))) {
          message( "Must define functions of the form: ")
          message( " p$stmv_global_model$run =' " )
          message( "   gam( formula=p$stmv_global_modelformula, data=B, ")
          message( "     optimizer= p$stmv_gam_optimizer, family=p$stmv_global_family, weights=wt )' " )
          stop()
        }
        global_model = try( eval(parse(text=p$stmv_global_model$run )) )
        if (inherits(global_model, "try-error") ) stop( "Global modelling failed")
        print( summary( global_model ) )
        save( global_model, file= fn.global_model, compress=TRUE )
        return ()
      }

      # as a first pass, model the time-independent factors as a user-defined model
      if (p$stmv_global_modelengine=="glm") {
        warning( "glm does not have a spatial random effect, model covariates could be biased if residual errors are not iid ... ")
        if (!exists("wt", B)) {
          global_model = try(
            glm( formula=p$stmv_global_modelformula, data=B, family=p$stmv_global_family ) )
        } else {
          global_model = try(
            glm( formula=p$stmv_global_modelformula, data=B, family=p$stmv_global_family, weights=wt ) )
        }
      }

      if (p$stmv_global_modelengine=="bigglm") {
        warning( "bigglm does not have a spatial random effect, model covariates could be biased if residual errors are not iid ... ")
        if (!exists("wt", B)) {
          global_model = try(
            bigglm( formula=p$stmv_global_modelformula, data=B, family=p$stmv_global_family ))
        } else {
          global_model = try(
            bigglm( formula=p$stmv_global_modelformula, data=B, family=p$stmv_global_family, weights=~wt ))
        }
      }

      if (p$stmv_global_modelengine=="gam") {
        warning( "gam does not have a spatial random effect, model covariates could be biased if residual errors are not iid, spatial splines might help but this is not guaranteed. gamm might be a better choice. You have been warned ... ")
        require(mgcv)
        if (!exists("wt", B)) {
          global_model = try(
            gam( formula=p$stmv_global_modelformula, data=B,
              optimizer= p$stmv_gam_optimizer, family=p$stmv_global_family) )
        } else {
          global_model = try(
            gam( formula=p$stmv_global_modelformula, data=B,
              optimizer= p$stmv_gam_optimizer, family=p$stmv_global_family, weights=wt ) )
        }
      }


      if (p$stmv_global_modelengine=="gamm") {
        warning( "gamm should include a spatial random effect, otherwise model covariates could be biased if residual errors are not iid ... ")
        require(mgcv)
        if (!exists("wt", B)) {
          global_model = try(
            gam( formula=p$stmv_global_modelformula, data=B,
              optimizer= p$stmv_gam_optimizer, family=p$stmv_global_family) )
        } else {
          global_model = try(
            gamm( formula=p$stmv_global_modelformula, data=B,
              optimizer= p$stmv_gam_optimizer, family=p$stmv_global_family, weights=wt ) )
        }
      }


      if (p$stmv_global_modelengine=="bayesx") {
        warning( "make sure you have a spatial random effect, otherwise model covariates/predictions could be biased if residual errors are not iid ... ")

        require(mgcv)
        global_model = try(
          bayesx( formula=p$stmv_global_modelformula, data=B, family=p$stmv_global_family ) )
      }


      if (p$stmv_global_modelengine=="inla") {
        warning( "make sure you have a spatial random effect, otherwise model covariates/predictions could be biased if residual errors are not iid ... ")

        require(INLA)
        global_model = try( eval( p$stmv_global_modelformula ) )    # p$stmv_global_modelformula must contain the whole call

          # inla( formula=p$stmv_global_modelformula, data=B, family=p$stmv_global_family ) )
      }

      if ( "try-error" %in% class(global_model) ) stop( "The covariate model was problematic" )
      print( summary( global_model ) )

      save( global_model, file= fn.global_model, compress=TRUE )
      global_model = NULL
      return()
    }


    # -----


    # -----


    if (DS %in% c("stmv.results", "stmv.prediction", "stmv.stats") )  {

      if (DS=="stmv.prediction") {
        if (! exists("TIME", p$variables)) {
          fn = file.path( p$stmvSaveDir, paste("stmv.prediction",  ret, "rdata", sep="." ) )
        } else {
          fn = file.path( p$stmvSaveDir, paste("stmv.prediction",  ret, yr, "rdata", sep="." ) )
        }
        if (file.exists(fn) ) load(fn)
        if (ret=="mean") return (P)
        if (ret=="lb") return( Pl)
        if (ret=="ub") return( Pu)
      }

      if (DS=="stmv.stats") {
        fn = file.path( p$stmvSaveDir, paste( "stmv.statistics", "rdata", sep=".") )
        stats = NULL
        if (file.exists(fn)) load(fn)
        return(stats)
      }

      shallower = NULL
      if ( exists("depth.filter", p) && is.finite( p$depth.filter) ) {
        if ( "z" %in% p$variables$COV ){
          depths = stmv_attach( p$storage.backend, p$ptr$Pcov[["z"]] )[]
          ii = which( depths[] < p$depth.filter )
          if (length(ii) > 0) shallower = ii
          ii= depths=NULL
        }
      }

      PP = stmv_attach( p$storage.backend, p$ptr$P )
      PPsd = stmv_attach( p$storage.backend, p$ptr$Psd )
      if (exists("stmv_global_modelengine", p)) {
        if (p$stmv_global_modelengine !="none" ) {
          P0 = stmv_attach( p$storage.backend, p$ptr$P0 )
          P0sd = stmv_attach( p$storage.backend, p$ptr$P0sd )
        }
      }

      if ( exists("TIME", p$variables)) {
        vv = ncol(PP)
        for (it in 1:p$ny) {
          y = p$yrs[it]
          fn_P = file.path( p$stmvSaveDir, paste("stmv.prediction", "mean",  y, "rdata", sep="." ) )
          fn_Pl = file.path( p$stmvSaveDir, paste("stmv.prediction", "lb",   y, "rdata", sep="." ) )
          fn_Pu = file.path( p$stmvSaveDir, paste("stmv.prediction", "ub",   y, "rdata", sep="." ) )
          if ( vv > p$ny ) {
            ww = (it-1) * p$nw + (1:p$nw)
            P = PP  [,ww]
            V = PPsd[,ww] # simpleadditive independent errors assumed
          } else if ( vv==p$ny) {
            P = PP[,it]
            V = PPsd[,it]
          }
          if (exists("stmv_global_modelengine", p) ) {
            if (p$stmv_global_modelengine !="none" ) {
              ## maybe add via simulation, note: P0 and P are on link scale to this point
              uu = which(!is.finite(P[]))
              if (length(uu)>0) P[uu] = 0 # permit global predictions to pass through ..
              P = P[] + P0[,it]
              vv = which(!is.finite(V[]))
              if (length(vv)>0) V[vv] = 0 # permit covariate-base predictions to pass through ..
              V = sqrt( V[]^2 + P0sd[,it]^2) # simple additive independent errors assumed
            }
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
          P =  P[]

          # return to user scale (that of Y)
          if ( exists( "stmv_global_family", p)) {
            if (p$stmv_global_family !="none" ) {
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

          save( P, file=fn_P, compress=T )
          save( Pl, file=fn_Pl, compress=T )
          save( Pu, file=fn_Pu, compress=T )
          print ( paste("Year:", y)  )
        }

      } else {
        # serial run only ...

          fn_P = file.path( p$stmvSaveDir, paste("stmv.prediction", "mean", "rdata", sep="." ) )
          fn_Pl = file.path( p$stmvSaveDir, paste("stmv.prediction", "lb", "rdata", sep="." ) )
          fn_Pu = file.path( p$stmvSaveDir, paste("stmv.prediction", "ub", "rdata", sep="." ) )

          P = PP[]
          V = PPsd[]
          if (exists("stmv_global_modelengine", p) ) {
            if (p$stmv_global_modelengine !="none" ) {
              uu = which(!is.finite(P[]))
              if (length(uu)>0) P[uu] = 0 # permit covariate-base predictions to pass through ..
              P = P[] + P0[]  # both on link scale
              vv = which(!is.finite(V[]))
              if (length(vv)>0) V[vv] = 0 # permit covariate-base predictions to pass through ..
              V = sqrt( V[]^2 + P0sd[]^2) # simple additive independent errors assumed
            }
          }
          if ( !is.null(shallower) ){
            P[shallower] = NA
            V[shallower] = NA
          }

          Pl =  P[] - 1.96* V[]
          Pu =  P[] + 1.96* V[]
          P =  P[]

          # return to user scale (that of Y)
          if ( exists( "stmv_global_family", p)) {
            if (p$stmv_global_family !="none" ) {
              Pl = p$stmv_global_family$linkinv( Pl[] )
              Pu = p$stmv_global_family$linkinv( Pu[] )
              P = p$stmv_global_family$linkinv( P[] )
            }
          }

          if (exists("stmv_Y_transform", p)) {
            Pl = p$stmv_Y_transform$invers (Pl[])  # p$stmv_Y_transform[2] is the inverse transform
            Pu = p$stmv_Y_transform$invers (Pu[])
            P = p$stmv_Y_transform$invers (P[])
          }

          save( P, file=fn_P, compress=T )
          save( Pl, file=fn_Pl, compress=T )
          save( Pu, file=fn_Pu, compress=T )
      } # end if TIME


      # prediction.stats

      Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
      S = stmv_attach( p$storage.backend, p$ptr$S )
      Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )

      Sloc_nplat = ceiling( diff( p$corners$plat) / p$stmv_distance_statsgrid)
      Sloc_nplon = ceiling( diff( p$corners$plon) / p$stmv_distance_statsgrid)

      stats = matrix( NaN, ncol=length( p$statsvars ), nrow=nrow( Ploc) )  # output data .. ff does not handle NA's .. using NaN for now
      colnames(stats)=p$statsvars

      for ( i in 1:length( p$statsvars ) ) {
        print(i)
        # linear interpolation
        u = as.image( S[,i], x=Sloc[,], na.rm=TRUE, nx=Sloc_nplon, ny=Sloc_nplat )
        stats[,i] = as.vector( fields::interp.surface( u, loc=Ploc[] ) ) # linear interpolation
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

      fn = file.path( p$stmvSaveDir, paste( "stmv.statistics", "rdata", sep=".") )
      save( stats, file=fn, compress=TRUE )

      if (0){
        i = 1
        ii = which (is.finite(stats[,i]))
        lattice::levelplot( stats[ii,i] ~ Ploc[ii,1]+Ploc[ii,2])
      }

      if(0) {
        i = 1
        Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
        Z = smooth.2d( Y=P[], x=Ploc[], ncol=p$nplats, nrow=p$nplons, cov.function=stationary.cov, Covariance="Matern", range=p$stmv_lowpass_phi, nu=p$stmv_lowpass_nu )
        dev.new(); image(Z)
        dev.new();lattice::levelplot( P[] ~ Ploc[,1] + Ploc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
      }

    }

    # ----------------


    if (DS %in% c("save_current_state") ) {

      # named differently to avoid collisions
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

      save( sP, file=p$saved_state_fn$P, compress=TRUE )
      save( sPn, file=p$saved_state_fn$Pn, compress=TRUE )
      save( sPsd, file=p$saved_state_fn$Psd, compress=TRUE )
      save( sS, file=p$saved_state_fn$stats, compress=TRUE )
      save( sSflag, file=p$saved_state_fn$sflag, compress=TRUE )

      if (exists("stmv_global_modelengine", p)) {
        if (p$stmv_global_modelengine !="none" ) {
          save( sP0,   file=p$saved_state_fn$P0,   compress=TRUE )
          save( sP0sd, file=p$saved_state_fn$P0sd, compress=TRUE )
        }
      }
      sP = sPsd = sPn = NULL
      sS = sSflag = NULL
      sP0 = sP0sd = NULL

    }


    # =--------------------


    if (DS %in% c("load_saved_state") ) {

      # named differently to avoid collisions
      P = stmv_attach( p$storage.backend, p$ptr$P )
      Pn = stmv_attach( p$storage.backend, p$ptr$Pn )
      Psd = stmv_attach( p$storage.backend, p$ptr$Psd )
      S = stmv_attach( p$storage.backend, p$ptr$S )
      Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )

      if (exists("stmv_global_modelengine", p)) {
        if (p$stmv_global_modelengine !="none" ) {
          P0 = stmv_attach( p$storage.backend, p$ptr$P0 )
          P0sd = stmv_attach( p$storage.backend, p$ptr$P0sd )
        }
      }

      sP = matrix( NaN, nrow=nrow(P), ncol=ncol(P) )
      sPn = matrix( NaN, nrow=nrow(Pn), ncol=ncol(Pn) )
      sPsd = matrix( NaN, nrow=nrow(Psd), ncol=ncol(Psd) )
      sS = matrix( NaN, nrow=nrow(S), ncol=ncol(S) )
      sSflag = matrix( NaN, nrow=nrow(Sflag), ncol=ncol(Sflag) )

      if (file.exists(p$saved_state_fn$P)) load( p$saved_state_fn$P )
      if (file.exists(p$saved_state_fn$Pn)) load( p$saved_state_fn$Pn )
      if (file.exists(p$saved_state_fn$Psd)) load( p$saved_state_fn$Psd )
      if (file.exists(p$saved_state_fn$stats)) load( p$saved_state_fn$stats )
      if (file.exists(p$saved_state_fn$sflag)) load( p$saved_state_fn$sflag )

      P[] = sP[]
      Pn[] = sPn[]
      Psd[] = sPsd[]
      sP = sPsd = sPn = NULL

      S[] = sS[]
      Sflag[] = sSflag[]
      sS = sSflag = NULL

      if (exists("stmv_global_modelengine", p)) {
        if (p$stmv_global_modelengine !="none" ) {
          sP0 = matrix( NaN, nrow=nrow(P0), ncol=ncol(P0) )
          sP0sd = matrix( NaN, nrow=nrow(P0sd), ncol=ncol(P0sd) )
          if (file.exists(p$saved_state_fn$P0)) load( p$saved_state_fn$P0 )
          if (file.exists(p$saved_state_fn$P0sd)) load( p$saved_state_fn$P0sd )
          P0[] = sP0[]
          P0sd[] = sP0sd[]
          sP0 = sP0sd = NULL
        }
      }
    }

    # =--------------------


  }
