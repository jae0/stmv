stmv_logfile = function(p, stime, local.n.start=NULL){

  varstoout = c("n.total", "n.shallow", "n.todo", "n.skipped", "n.predareaerror", "n.nodata", "n.variogramerror",
    "n.vrangeerror", "n.modelerror", "n.outside", "n.complete", "prop_incomp" )
  header = paste( c( varstoout) )
  currentstatus = stmv_db( p=p, DS="statistics.status" )
  currentstatus = c( unlist( currentstatus[ varstoout ] ) )
  dtime = difftime( Sys.time(), stime )
  dtimehr = difftime( Sys.time(), stime, units="hours" )
  
  if (is.null(local.n.start)) {
    local.n.complete= currentstatus["n.complete"] - local.n.start
  } else {
    local.n.complete= currentstatus["n.complete"] 
  }
  nrate = local.n.complete/ as.numeric(dtimehr)
  tmore = currentstatus["n.todo"] / nrate
  tall = (currentstatus["n.todo"]+local.n.complete) / nrate
  cat( paste( "---", p$data_root, p$variables$Y, p$spatial.domain, "--- \n\n"), file=p$stmv_current_status, append=FALSE )
  cat( paste( "Rate (no. per hour) :  ", round(nrate,3), "\n"), file=p$stmv_current_status, append=TRUE )
  nclusters = 1
  if (exists("clusters", p)) nclusters = length(p$clusters)
  cat( paste( "Rate (no. per hour per core) :  ", round(nrate/nclusters,3), "\n"), file=p$stmv_current_status, append=TRUE )
  cat( paste( "Core start time :  ", stime, "\n"), file=p$stmv_current_status, append=TRUE )
  cat( paste( "Current time :", Sys.time(), "\n"), file=p$stmv_current_status, append=TRUE )
  cat( paste( "Elapsed time :", format(dtime), "\n" ), file=p$stmv_current_status, append=TRUE)
  cat( paste( "Estimated time remaining (hrs) :", round( tmore,3), "\n" ), file=p$stmv_current_status, append=TRUE)
  cat( paste( "Estimated time total (hrs) :", round( tall,3), "\n" ), file=p$stmv_current_status, append=TRUE)
  for ( hd in varstoout ){
    cat( paste( hd, ":", currentstatus[hd], "\n" ), file=p$stmv_current_status, append=TRUE)
  }
  # message( readLines( p$stmv_current_status ) )
  return(currentstatus)
}
