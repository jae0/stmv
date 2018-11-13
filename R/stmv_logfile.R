stmv_logfile = function(p){

  time_current = Sys.time()

  time_interpolation = time_current
  t_suffix = ""

  if (exists("time_start_interpolation", p)) {
    time_interpolation = p$time_start_interpolation
    t_suffix = "main"
  }
  if (exists("time_start_interpolation_debug", p)) {
    time_interpolation = p$time_start_interpolation_debug
    t_suffix = "debug"
  }
  if (exists("time_start_interpolation_force_complete", p)) {
    time_interpolation = p$time_start_interpolation_force_complete
    t_suffix = "force-complete"
  }

  dtime_totalelapsed = difftime( time_current, p$time_start )
  dtime = difftime( time_current, time_interpolation )
  dtimehr = difftime( time_current, time_interpolation, units="hours" )

  varstoout = c("n.total", "n.shallow", "n.todo", "n.skipped", "n.predareaerror", "n.nodata", "n.variogramerror",
    "n.vrangeerror", "n.modelerror", "n.outside", "n.complete", "prop_incomp" )
  header = paste( c( varstoout) )
  currentstatus = stmv_db( p=p, DS="statistics.status" )
  currentstatus = c( unlist( currentstatus[ varstoout ] ) )

  nrate = currentstatus["n.complete"]/ as.numeric(dtimehr)
  tmore = currentstatus["n.todo"] / nrate
  tall = (currentstatus["n.todo"] + currentstatus["n.complete"]) / nrate

  fn = p$stmv_current_status  # reduce text

  nclusters = 1
  if (exists("clusters", p)) nclusters = length(p$clusters)

  cat( paste( "---", p$data_root, p$variables$Y, p$spatial.domain, "--- \n\n"), file=fn, append=FALSE )
  cat( paste( "stmv start time :", p$time_start, "\n"), file=fn, append=TRUE )
  cat( paste( paste("Interpolation start time (", t_suffix, "):", sep=""), time_interpolation, "\n"), file=fn, append=TRUE )
  cat( paste( "Current time :", time_current, "\n"), file=fn, append=TRUE )
  cat( paste( "Total elapsed time :", format(dtime_totalelapsed), "\n" ), file=fn, append=TRUE)
  cat( paste( paste("Time spent interpolating (", t_suffix, "):", sep=""), format(dtime), "\n" ), file=fn, append=TRUE)
  cat( paste( "Rate (no. per hour) :  ", round(nrate,3), "\n"), file=fn, append=TRUE )
  cat( paste( "Rate (no. per hour per core) :  ", round(nrate/nclusters,3), "\n"), file=fn, append=TRUE )
  cat( paste( "Estimated time remaining (hrs) :", round( tmore,3), "\n" ), file=fn, append=TRUE)
  cat( paste( "Estimated time total (hrs) :", round( tall,3), "\n" ), file=fn, append=TRUE)
  for ( hd in varstoout ){
    cat( paste( hd, ":", currentstatus[hd], "\n" ), file=fn, append=TRUE)
  }
  # message( readLines( p$stmv_current_status ) )
  return(currentstatus)
}
