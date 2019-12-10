stmv_logfile = function(p, flag="default"){

  time_current = Sys.time()
  time_runmode = p$time_start_runmode

  dtime_totalelapsed = difftime( time_current, p$time_start )
  dtime = difftime( time_current, time_runmode )
  dtimehr = difftime( time_current, time_runmode, units="hours" )

  varstoout = c(
    "n.total",
    "n.too_shallow",
    "n.todo",
    "n.unknown",
    "n.outside_bounds",
    "n.insufficient_data",
    "n.variogram_range_limit",
    "n.variogram_failure",
    "n.local_model_error",
    "n.prediction_area",
    "n.prediction_error",
    "n.prediction_update_error",
    "n.statistics_update_error",
    "n.complete",
    "proportion_complete",
    "proportion_incomplete" )
  header = paste( c( varstoout) )
  currentstatus = stmv_statistics_status( p=p )
  currentstatus = c( unlist( currentstatus[ varstoout ] ) )

  nrate = currentstatus["n.complete"]/ as.numeric(dtimehr)
  tmore = currentstatus["n.todo"] / nrate
  tall = (currentstatus["n.todo"] + currentstatus["n.complete"]) / nrate

  fn = p$stmv_current_status  # reduce text

  nclusters = 1
  if (exists("clusters", p)) nclusters = length(p$clusters)

  cat( paste( "---", p$data_root, p$stmv_variables$Y, p$spatial_domain, "--- \n\n"), file=fn, append=FALSE )
  cat( paste( "Runmode : ", flag, "\n"), file=fn, append=TRUE )
  cat( paste( "Start time :", p$time_start, "\n"), file=fn, append=TRUE )
  cat( paste( "Start time of current runmode : ", time_runmode, "\n"), file=fn, append=TRUE )
  cat( paste( "Current time :", time_current, "\n"), file=fn, append=TRUE )
  cat( paste( "Total elapsed time :", format(dtime_totalelapsed), "\n" ), file=fn, append=TRUE)
  cat( paste( "Time spent in current runmode :", format(dtime), "\n" ), file=fn, append=TRUE)
  cat( paste( "Runmode rate (no. per hour) :  ", round(nrate,3), "\n"), file=fn, append=TRUE )
  cat( paste( "Runmode rate (no. per hour per core) :  ", round(nrate/nclusters,3), "\n"), file=fn, append=TRUE )
  cat( paste( "Runmode estimated time remaining (hrs) :", round( tmore,3), "\n" ), file=fn, append=TRUE)
  cat( paste( "Runmode estimated time total (hrs) :", round( tall,3), "\n" ), file=fn, append=TRUE)
    for ( hd in varstoout ){
    cat( paste( hd, ":", currentstatus[hd], "\n" ), file=fn, append=TRUE)
  }
  # message( readLines( p$stmv_current_status ) )
  return(currentstatus)
}
