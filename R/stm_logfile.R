stm_logfile = function(p, stime){

  varstoout = c("n.total", "n.shallow", "n.todo", "n.skipped", "n.outsidepreddomain",
    "n.outside", "n.complete", "prop_incomp" )
  header = paste( c( varstoout) )
  currentstatus = stm_db( p=p, DS="statistics.status" )
  currentstatus = c( unlist( currentstatus[ varstoout ] ) )
  dtime = difftime( Sys.time(), stime )
  dtimehr = difftime( Sys.time(), stime, units="hours" )
  nrate = currentstatus["n.complete"]/ as.numeric(dtimehr)
  tmore = currentstatus["n.todo"] / nrate
  tall = (currentstatus["n.todo"]+currentstatus["n.complete"]) / nrate
  cat( paste( "---", p$data_root, p$variables$Y, p$spatial.domain, "--- \n\n"), file=p$stm_current_status, append=FALSE )
  cat( paste( "Core start time :  ", stime, "\n"), file=p$stm_current_status, append=TRUE )
  cat( paste( "Current time :", Sys.time(), "\n"), file=p$stm_current_status, append=TRUE )
  cat( paste( "Elapsed time :", format(dtime), "\n" ), file=p$stm_current_status, append=TRUE)
  cat( paste( "Estimated time remaining (hrs) :", round( tmore,3), "\n" ), file=p$stm_current_status, append=TRUE)
  cat( paste( "Estimated time total (hrs) :", round( tall,3), "\n" ), file=p$stm_current_status, append=TRUE)
  for ( hd in varstoout ){
    cat( paste( hd, ":", currentstatus[hd], "\n" ), file=p$stm_current_status, append=TRUE)
  }
  # message( readLines( p$stm_current_status ) )
  return(currentstatus)
}
