
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
      for (fn in unlist(p$cache) ) if (length(fn)>0) if (file.exists(fn)) file.remove(fn)
      for (fn in unlist(p$bm) ) if (length(fn)>0)  if (file.exists(fn)) file.remove(fn)
      return( NULL )
    }

    # --------------------------
    if (DS %in% "cleanup.all" ) {
      stmv_db( p=p, DS="cleanup")
      snapshots = c(
        file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "mean", "rdata", sep="." ) ),
        file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "sd", "rdata", sep="." ) ),
        file.path( p$stmvSaveDir, paste( "tmp_stmv.statistics", "rdata", sep=".") ),
        file.path( p$stmvSaveDir, paste( "tmp_stmv.sflag", "rdata", sep=".") ),
        file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "mean0", "rdata", sep="." ) ),
        file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "sd0",  "rdata", sep="." ) )
      )
      for (fn in snapshots ) if (file.exists(fn)) file.remove(fn)
      return( NULL )
    }


    # -----------------

    if ( DS %in% c( "statistics.status", "statistics.status.reset") ) {
    
      Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )
      
      out = list()
      
      out$todo = which( Sflag[]==0L )       # 0 = TODO
      out$done = which( Sflag[]==1L )       # 1 = completed
      out$outside = which( Sflag[]==2L )    # 2 = oustide bounds(if any)
      out$shallow = which( Sflag[]==3L )    # 3 = depth shallower than p$depth.filter (if it exists .. z is a covariate)
      out$predareaerror = which( Sflag[]==4L ) # 4=predictionarea not ok,
      out$nodata = which( Sflag[]==5L )     # 5=skipped due to insufficient data,
      out$variogramerror = which( Sflag[]==6L ) # 6=skipped .. fast variogram did not work
      out$vrangeerror = which( Sflag[]==7L )     # 7=variogram estimated range not ok
      out$modelerror = which( Sflag[]==8L )     # 8=problem with prediction and/or modelling
      out$skipped = which( Sflag[] == 9L )   # 9 not completed due to a failed attempt
      
      if ( DS=="statistics.status.reset" ) {
        # to reset all rejected locations
        if (length(which(is.finite(out$skipped))) > 0) {
          Sflag[out$skipped] = 0L  # to reset all the problem flags to todo
            out$skipped = which( Sflag[] == 9L )
          Sflag[out$predareaerror] = 0L  # to reset all the problem flags to todo
            out$predareaerror = which( Sflag[] == 4L )
          Sflag[out$variogramerror] = 0L  # to reset all the problem flags to todo
            out$variogramerror = which( Sflag[] == 6L )
          Sflag[out$vrangeerror] = 0L  # to reset all the problem flags to todo
            out$vrangeerror = which( Sflag[] == 7L )
          Sflag[out$modelerror] = 0L  # to reset all the problem flags to todo
            out$modelerror = which( Sflag[] == 8L )
          out$todo = which( Sflag[]==0L )
        }
      }
      
      # do some counts
      out$n.todo = length(out$todo)
      out$n.complete = length(out$done) 
      out$n.outside = length(which(is.finite(out$outside))) 
      out$n.shallow = length(out$shallow)
      out$n.predareaerror = length(out$predareaerror)
      out$n.nodata = length(out$nodata)
      out$n.variogramerror = length(out$variogramerror)
      out$n.vrangeerror = length(out$vrangeerror)
      out$n.modelerror = length(out$modelerror) 
      out$n.skipped = length(out$skipped)
      out$n.total = length(Sflag) 

      out$prop_incomp = round( out$n.todo / ( out$n.todo + out$n.complete), 3)
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
        points( Sloc[which( Sflag[]== 0L),], pch=".", col="blue", cex=5 )
        points( Sloc[which( Sflag[]== 1L),], pch=".", col="purple", cex=5 )
        points( Sloc[which( Sflag[]== 2L),], pch=".", col="red", cex=5 )
        points( Sloc[which( Sflag[]== 3L),], pch=".", col="yellow", cex=5 )
        points( Sloc[which( Sflag[]== 4L),], pch=".", col="green", cex=5 )
        points( Sloc[which( Sflag[]== 5L),], pch=".", col="red2", cex=5 )
        points( Sloc[which( Sflag[]== 6L),], pch=".", col="yellow2", cex=5 )
        points( Sloc[which( Sflag[]== 7L),], pch=".", col="green2", cex=5 )
        points( Sloc[which( Sflag[]== 8L),], pch=".", col="green3", cex=5 )
        points( Sloc[which( Sflag[]== 9L),], pch=".", col="magenta", cex=5 )

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
      Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )
      if (length(outside)  > 0 ) Sflag[outside] = 2L  # outside of prediction domain

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
            if (length(outside)>0) Sflag[outside] = 2L
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
          if (length(below) > 0 ) Sflag[below] = 0L
          if (length(above) > 0 ) Sflag[above] = 3L

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
            points( Sloc[which( Sflag[]== 0L),], pch=".", col="blue" )
          }
        }
      }

    }


    #---------------------


    if (DS== "flag.incomplete.predictions") {
      # statistics locations where estimations need to be redone
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
        Sflag[uP] = 0L  # force set to redo
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
        Ycov = stmv_attach(  p$storage.backend, p$ptr$Ycov )
        if ( length( p$variables$COV ) == 1 ) {
          bad = which( !is.finite( Ycov[]) )
        } else {
          bad = which( !is.finite( rowSums(Ycov[])) )
        }
        if (length(bad)> 0 ) hasdata[bad] = NA
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
        message( "|||   leave out 'globalmodel' as a runmode option ... overwriting in:")
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

      # as a first pass, model the time-independent factors as a user-defined model
      if (p$stmv_global_modelengine=="glm") {
        if (!exists("wt", B)) {
          global_model = try(
            glm( formula=p$stmv_global_modelformula, data=B, family=p$stmv_global_family ) )
        } else {
          global_model = try(
            glm( formula=p$stmv_global_modelformula, data=B, family=p$stmv_global_family, weights=wt ) )
        }
      }

      if (p$stmv_global_modelengine=="bigglm") {
        if (!exists("wt", B)) {
          global_model = try(
            bigglm( formula=p$stmv_global_modelformula, data=B, family=p$stmv_global_family ))
        } else {
          global_model = try(
            bigglm( formula=p$stmv_global_modelformula, data=B, family=p$stmv_global_family, weights=~wt ))
        }
      }

      if (p$stmv_global_modelengine=="gam") {
        require(mgcv)
        if (!exists("wt", B)) {
          global_model = try(
            gam( formula=p$stmv_global_modelformula, data=B, optimizer=c("outer","bfgs"), family=p$stmv_global_family) )
        } else {
          global_model = try(
            gam( formula=p$stmv_global_modelformula, data=B, optimizer=c("outer","bfgs"), family=p$stmv_global_family, weights=wt ) )
        }
      }

      if (p$stmv_global_modelengine=="bayesx") {
        require(mgcv)
        global_model = try(
          bayesx( formula=p$stmv_global_modelformula, data=B, family=p$stmv_global_family ) )
      }

      if ( "try-error" %in% class(global_model) ) stop( "The covariate model was problematic" )
      print( summary( global_model ) )
      
      save( global_model, file= fn.global_model, compress=TRUE )

      return (fn.global_model)
    }

    # -----cpu

    if (DS %in% c("global.prediction.surface") ) {
  
      parallel_run(
        p=p, 
        runindex=list(tindex=1:p$nt),
        FUNC = function( ip=NULL, p ) {
          if (exists( "libs", p)) suppressMessages( RLibrary( p$libs ) )
          if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns
          global_model = stmv_db( p=p, DS="global_model")
          if (is.null(global_model)) stop("Global model not found.")
          P0 = stmv_attach( p$storage.backend, p$ptr$P0 )
          P0sd = stmv_attach( p$storage.backend, p$ptr$P0sd )
          for ( ii in ip ) {
            # downscale and warp from p(0) -> p1
            it = p$runs$tindex[ii]  # == ii btw
            pa = NULL # construct prediction surface
            for (i in p$variables$COV ) {
              pu = stmv_attach( p$storage.backend, p$ptr$Pcov[[i]] )
              nc = ncol(pu)
              if ( nc== 1 ) {
                pa = cbind( pa, pu[] ) # ie. a static variable (space)
              } else if ( nc == p$nt & nc == p$ny) {
                pa = cbind( pa, pu[,it] ) # ie. same time dimension as predictive data (space.annual.seasonal)
              } else if ( nc == p$ny & p$nt > p$ny)  {
                iy = round( (it-1) / p$nw ) + 1
                pa = cbind( pa, pu[,iy] ) # ie., annual data (space.annual)
              } else if ( nc == p$nt & p$nt > p$ny) {
                pa = cbind( pa, pu[,it] ) # ie. same time dimension as predictive data (space.annual.seasonal)
              } else {
                print(i)
                print(nc)
                stop( "Erroneous data dimension ... the dataset for the above variable looks to be incomplete?")
              }
            }
            pa = as.data.frame( pa )
            names(pa) = p$variables$COV

            if ( any( p$variables$LOCS %in%  all.vars( p$stmv_global_modelformula ) ) ) {
              Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
              pa = cbind(pa, Ploc[])
              names(pa) = c( p$variables$COV, p$variables$LOCS )
            }

            if ( "yr" %in%  all.vars( p$stmv_global_modelformula ) ) {
              npa = names(pa)
              pa = cbind(pa, p$yrs[it] )
              names(pa) = c( npa, "yr" )
            }

            if ( "dyear" %in%  all.vars( p$stmv_global_modelformula ) ) {
              npa = names(pa)
              pa = cbind(pa, p$prediction.dyear )
              names(pa) = c( npa, "dyear" )
            }

            if (p$stmv_global_modelengine %in% c("glm", "bigglm", "gam") ) {
              Pbaseline = try( predict( global_model, newdata=pa, type="response", se.fit=TRUE ) )
              pa = NULL
              gc()
              if (!inherits(Pbaseline, "try-error")) {
                P0[,it] = Pbaseline$fit
                P0sd[,it] = Pbaseline$se.fit
              }
              Pbaseline = NULL; gc()
            } else if (p$stmv_global_modelengine =="bayesx") {
              stop( "not yet tested" ) # TODO
              # Pbaseline = try( predict( global_model, newdata=pa, type="response", se.fit=TRUE ) )
              # pa = NULL
              # gc()
              # if (!inherits(Pbaseline, "try-error")) {
              #   P0[,it] = Pbaseline$fit
              #   P0sd[,it] = Pbaseline$se.fit
              # }
              # Pbaseline = NULL; gc()

            } else if (p$stmv_global_modelengine =="none") {
              # nothing to do
            } else  {
              stop ("This global model method requires a bit more work .. ")
            }

            if (p$all.covars.static) {
              # if this is true then this is a single cpu run and all predictions for each time slice is the same
              # could probably catch this and keep storage small but that would make the update math a little more complex
              # this keeps it simple with a quick copy
              if (p$nt  > 1 ) {
                for (j in ip[2:p$nruns]){
                  P0[,j] = P0[,1]
                  P0sd[,j] = P0sd[,1]
                }
              }
              global_model =NULL
              return(NULL)
            }
          } # end each timeslice
          global_model =NULL
          return(NULL)
        }
      )
      message( "||| Done ... moving onto the rest of the analysis...")
    }


    # -----


    if (DS %in% c("stmv.prediction.redo", "stmv.prediction") )  {

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


      shallower = NULL
      if ( exists("depth.filter", p) && is.finite( p$depth.filter) ) {
        if ( "z" %in% p$variables$COV ){
          depths = stmv_attach( p$storage.backend, p$ptr$Pcov[["z"]] )[]
          ii = which( depths[] < p$depth.filter )
          if (length(ii) > 0) shallower = ii
          rm(depths)
        }
      }

      if ( exists("TIME", p$variables)) {
        clusters = p$clusters
        runindex = list( tindex=1:p$ny ) # annual only
      } else {
        clusters = p$clusters[1]  # force serial mode
        runindex = list( tindex=1 )  # dummy value
      }

      parallel_run(
        p=p, 
        clusters=clusters, # override
        runindex=runindex,
        shallower=shallower,
        FUNC = function( ip=NULL, p, shallower ) {
          if (exists( "libs", p)) suppressMessages( RLibrary( p$libs ) )
          if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns

          PP = stmv_attach( p$storage.backend, p$ptr$P )
          PPsd = stmv_attach( p$storage.backend, p$ptr$Psd )
          if (exists("stmv_global_modelengine", p)) {
            if (p$stmv_global_modelengine !="none" ) {
              P0 = stmv_attach( p$storage.backend, p$ptr$P0 )
              P0sd = stmv_attach( p$storage.backend, p$ptr$P0sd )
            }
          }

          if ( exists("TIME", p$variables)) {
            
            for ( ii in ip ) { # outputs are on yearly breakdown
              it = p$runs[ii, "tindex"] # == ii btw
              y = p$yrs[it]
              
              fn_P = file.path( p$stmvSaveDir, paste("stmv.prediction", "mean",  y, "rdata", sep="." ) )
              fn_Pl = file.path( p$stmvSaveDir, paste("stmv.prediction", "lb",   y, "rdata", sep="." ) )
              fn_Pu = file.path( p$stmvSaveDir, paste("stmv.prediction", "ub",   y, "rdata", sep="." ) )
              vv = ncol(PP)
              if ( vv > p$ny ) {
                col.ranges = (it-1) * p$nw + (1:p$nw)
                P = PP  [,col.ranges]
                V = PPsd[,col.ranges] # simpleadditive independent errors assumed
              } else if ( vv==p$ny) {
                P = PP[,it]
                V = PPsd[,it]
              }

              if (exists("stmv_global_modelengine", p) ) {
                if (p$stmv_global_modelengine !="none" ) {
                  ## maybe add via simulation ? ...
                  uu = which(!is.finite(P[]))
                  if (length(uu)>0) P[uu] = 0 # permit covariate-base predictions to pass through ..
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

              # return to user scale (that of Y)
              Pl = p$stmv_global_family$linkinv( P + 1.96* V )
              Pu = p$stmv_global_family$linkinv( P - 1.96* V )
              P = p$stmv_global_family$linkinv( P )

              if (exists("stmv_Y_transform", p)) {
                Pl = p$stmv_Y_transform[[2]] (Pl)  # p$stmv_Y_transform[2] is the inverse transform
                Pu = p$stmv_Y_transform[[2]] (Pu)
                P = p$stmv_Y_transform[[2]] (P)
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
                  P = P[] + P0[]
                  vv = which(!is.finite(V[]))
                  if (length(vv)>0) V[vv] = 0 # permit covariate-base predictions to pass through ..
                  V = sqrt( V[]^2 + P0sd[]^2) # simple additive independent errors assumed
                }
              }
              if ( !is.null(shallower) ){
                P[shallower,] = NA
                V[shallower,] = NA
              }

              # return to user scale
              Pl = p$stmv_global_family$linkinv( P + 1.96* V )
              Pu = p$stmv_global_family$linkinv( P - 1.96* V )
              P = p$stmv_global_family$linkinv( P )


              if (exists("stmv_Y_transform", p)) {
                Pl = p$stmv_Y_transform[[2]] (Pl)  # p$stmv_Y_transform[2] is the inverse transform
                Pu = p$stmv_Y_transform[[2]] (Pu)
                P = p$stmv_Y_transform[[2]] (P)
              }

              save( P, file=fn_P, compress=T )
              save( Pl, file=fn_Pl, compress=T )
              save( Pu, file=fn_Pu, compress=T )
          } # end if TIME
        } # end FUNC
      ) # end parallel_run
    
      if(0) {
        i = 1
        Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )

        Z = smooth.2d( Y=P[], x=Ploc[], ncol=p$nplats, nrow=p$nplons, cov.function=stationary.cov, Covariance="Matern", range=p$stmv_lowpass_phi, nu=p$stmv_lowpass_nu )
        lattice::levelplot( P[] ~ Ploc[,1] + Ploc[,2], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
      }
    
    }

    # ----------------


    if (DS %in% c("stats.to.prediction.grid.redo", "stats.to.prediction.grid") ) {
      
      # TODO:: parallelize this
      
      fn = file.path( p$stmvSaveDir, paste( "stmv.statistics", "rdata", sep=".") )
      if (DS=="stats.to.prediction.grid") {
        stats = NULL
        if (file.exists(fn)) load(fn)
        return(stats)
      }
      
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
      if (!is.null(boundary)) {
      if( !("try-error" %in% class(boundary) ) ) {
        inside.polygon = point.in.polygon( Ploc[,1], Ploc[,2],
          boundary$polygon[,1], boundary$polygon[,2], mode.checked=TRUE )
          o = which( inside.polygon == 0 ) # outside boundary
          if (length(o) > 0) stats[o,] = NA
        }
      }
      
      
      if (0){
        i = 1
        ii = which (is.finite(stats[,i]))
        lattice::levelplot( stats[ii,i] ~ Ploc[ii,1]+Ploc[ii,2])
      }
      
      if ( exists("depth.filter", p) && is.finite( p$depth.filter) ) {
        # stats is now with the same indices as Pcov, Ploc, etc..
        if ( "z" %in% p$variables$COV ){
          depths = stmv_attach( p$storage.backend, p$ptr$Pcov[["z"]] )[]
          shallower = which( depths < p$depth.filter )
          if (length(shallower)>0) stats[shallower,] = NA
          rm(shallower); gc()
        }
      }
      
      save( stats, file=fn, compress=TRUE )
      
    }
    
    #-------------

    
    if (DS %in% c("save_current_state") ) {

      # named differently to avoid collisions 
      P = stmv_attach( p$storage.backend, p$ptr$P )[]
      Psd = stmv_attach( p$storage.backend, p$ptr$Psd )[]
      S = stmv_attach( p$storage.backend, p$ptr$S )[]
      Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )[]

      if (exists("stmv_global_modelengine", p)) {
        if (p$stmv_global_modelengine !="none" ) {
          P0 = stmv_attach( p$storage.backend, p$ptr$P0 )[]
          P0sd = stmv_attach( p$storage.backend, p$ptr$P0sd )[]
        }
      }

      fn_P = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "mean", "rdata", sep="." ) )
      fn_Psd = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "sd", "rdata", sep="." ) )
      fn_stats = file.path( p$stmvSaveDir, paste( "tmp_stmv.statistics", "rdata", sep=".") )
      fn_sflag = file.path( p$stmvSaveDir, paste( "tmp_stmv.sflag", "rdata", sep=".") )
      
      save( P, file=fn_P, compress=TRUE )
      save( Psd, file=fn_Psd, compress=TRUE )
      save( S, file=fn_stats, compress=TRUE )
      save( Sflag, file=fn_sflag, compress=TRUE )

      if (exists("stmv_global_modelengine", p)) {
        if (p$stmv_global_modelengine !="none" ) {
          fn_P0 = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "mean0", "rdata", sep="." ) )
          fn_P0sd = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "sd0", "rdata", sep="." ) )
          save( P0,   file=fn_P0,   compress=TRUE )
          save( P0sd, file=fn_P0sd, compress=TRUE )
        }
      }
      
    }


    # ---------------------------


    if (DS %in% c("load_saved_state") ) {

      # named differently to avoid collisions 

      # named differently to avoid collisions 
      PP = stmv_attach( p$storage.backend, p$ptr$P )
      PPsd = stmv_attach( p$storage.backend, p$ptr$Psd )
      SS = stmv_attach( p$storage.backend, p$ptr$S )
      SSflag = stmv_attach( p$storage.backend, p$ptr$Sflag )

      if (exists("stmv_global_modelengine", p)) {
        if (p$stmv_global_modelengine !="none" ) {
          PP0 = stmv_attach( p$storage.backend, p$ptr$P0 )
          PP0sd = stmv_attach( p$storage.backend, p$ptr$P0sd )
        }
      }

      fn_P = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "mean", "rdata", sep="." ) )
      fn_Psd = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "sd", "rdata", sep="." ) )
      fn_stats = file.path( p$stmvSaveDir, paste( "tmp_stmv.statistics", "rdata", sep=".") )
      fn_sflag = file.path( p$stmvSaveDir, paste( "tmp_stmv.sflag", "rdata", sep=".") )

      load( fn_P )
      load( fn_Psd )
      load( fn_stats )
      load( fn_sflag )

      PP[] = P[]
      PPsd[] = Psd[]
      SS[] = S[]
      SSflag[] = Sflag[]
      
      if (exists("stmv_global_modelengine", p)) {
        if (p$stmv_global_modelengine !="none" ) {
          fn_P0 = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "mean0", "rdata", sep="." ) )
          fn_P0sd = file.path( p$stmvSaveDir, paste("tmp_stmv.prediction", "sd0",  "rdata", sep="." ) )
          load( fn_P0 )
          load( fn_P0sd )
          PP0[] = P0[]
          PP0sd[] = P0sd[]
        }
      }

    }

  }
