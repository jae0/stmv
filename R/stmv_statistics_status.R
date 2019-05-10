stmv_statistics_status = function(p, plotdata=FALSE, reset=NULL  ) {

  E = stmv_error_codes()
  Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )

  out = list()

  if (!is.null(reset)) {
    if ("incomplete" %in% reset) {
      locs_to_do = stmv_predictions_incomplete( p=p )
      if ( !is.null(locs_to_do) && length(locs_to_do) > 0) {
        Sflag[locs_to_do] = stmv_error_codes()[["todo"]]
      }
      locs_to_do = NULL

    } else if ("features"  %in% reset) {

      # create location specific flags for analysis, etc..

      # flag areas overlapping with prediction locations:
      Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
      Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )

      pidP = array_map( "xy->1", Ploc, gridparams=p$gridparams )
      pidS = array_map( "xy->1", Sloc, gridparams=p$gridparams )
      overlap = match( pidS, pidP )
      outside_bounds = which( !is.finite( overlap ))

      if (length(outside_bounds)  > 0 ) Sflag[outside_bounds] = E[["outside_bounds"]]  # outside of prediction domain

      # catch data boundaries if present
      if (exists( "boundary", p) && p$boundary) {
        timeb0 =  Sys.time()
        message("\n")
        message( "||| Defining boundary polygon for data .. this reduces the number of points to analyse")
        message( "||| but takes a few minutes to set up ...")
        stmv_db( p=p, DS="boundary.redo" ) # ~ 5 min on nfs
      # last set of filters to reduce problem size
        bnds = try( stmv_db( p=p, DS="boundary" ) )
        if (!is.null(bnds)) {
          if( !("try-error" %in% class(bnds) ) ) {
            outside_bounds = which( bnds$inside.polygon == 0 ) # outside boundary
            if (length(outside_bounds)>0) Sflag[outside_bounds] = E[["outside_bounds"]]
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

      return(NULL)  # done!

    } else {
      Eflags_reset = E[ reset ]
      toreset = which( Sflag[] %in% unlist(Eflags_reset) )
      if (length(toreset) > 0) Sflag[toreset] = E[["todo"]]
      toreset = NULL
    }
  }

  out$todo = which( Sflag[]==E[["todo"]] )       # 0 = TODO
  out$complete = which( Sflag[]==E[["complete"]] )       # 1 = completed
  out$outside_bounds = which( Sflag[]==E[["outside_bounds"]] )    # 2 = oustide bounds(if any)
  out$too_shallow = which( Sflag[]==E[["too_shallow"]] )    # 3 = depth shallower than p$depth.filter (if it exists .. z is a covariate)
  out$prediction_area = which( Sflag[]==E[["prediction_area"]] ) # 4=predictionarea not ok,
  out$local_model_error = which( Sflag[]==E[["local_model_error"]] ) # 4=predictionarea not ok,
  out$insufficient_data = which( Sflag[]==E[["insufficient_data"]] )     # 5=skipped due to insufficient data,
  out$variogram_failure = which( Sflag[]==E[["variogram_failure"]] ) # 6=skipped .. fast variogram did not work
  out$variogram_range_limit = which( Sflag[]==E[["variogram_range_limit"]] )     # 7=variogram estimated range not ok
  out$prediction_error = which( Sflag[]==E[["prediction_error"]] )     # 8=problem with prediction and/or modelling
  out$prediction_update_error = which( Sflag[]==E[["prediction_update_error"]] )     # 8=problem with prediction and/or modelling
  out$statistics_update_error = which( Sflag[]==E[["statistics_update_error"]] ) # 4=predictionarea not ok,
  out$unknown = which( Sflag[] == E[["unknown"]] )   # 9 not completed due to a failed attempt

  # do some counts
  out$n.todo = length(out$todo)
  out$n.complete = length(out$complete)
  out$n.outside_bounds = length(which(is.finite(out$outside_bounds)))
  out$n.too_shallow = length(out$too_shallow)
  out$n.prediction_area = length(out$prediction_area)
  out$n.local_model_error = length(out$local_model_error)
  out$n.insufficient_data = length(out$insufficient_data)
  out$n.variogram_failure = length(out$variogram_failure)
  out$n.variogram_range_limit = length(out$variogram_range_limit)
  out$n.prediction_error = length(out$prediction_error)
  out$n.prediction_update_error = length(out$prediction_update_error)
  out$n.statistics_update_error = length(out$statistics_update_error)
  out$n.unknown = length(out$unknown)
  out$n.total = length(Sflag)

  out$prop_incomp = round( out$n.todo / ( out$n.total - out$n.outside_bounds ), 3)
  message( paste("||| Proportion to do:", out$prop_incomp, "\n" ))

  if (plotdata) {
    dev.new()
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

  return( out )

}

