stmv_logfile = function(p, stime, currentstatus){

  varstoout = c("n.total", "n.shallow", "n.todo", "n.skipped", "n.rangeissue", "n.nodata", "n.variogramerror",
    "n.vrangeerror", "n.modelerror", "n.outside", "n.complete", "prop_incomp" )
  header = paste( c( varstoout) )
  currentstatus = stmv_db( p=p, DS="statistics.status" )
  currentstatus = c( unlist( currentstatus[ varstoout ] ) )
  dtime = difftime( Sys.time(), stime )
  dtimehr = difftime( Sys.time(), stime, units="hours" )
  nrate = currentstatus["n.complete"]/ as.numeric(dtimehr)
  tmore = currentstatus["n.todo"] / nrate
  tall = (currentstatus["n.todo"]+currentstatus["n.complete"]) / nrate
  cat( paste( "---", p$data_root, p$variables$Y, p$spatial.domain, "--- \n\n"), file=p$stmv_current_status, append=FALSE )
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


  # numeric codes .. see stmv_db("statistics.status")
  #     out$todo = which( Sflag[]==0L )       # 0 = TODO
  #     out$done = which( Sflag[]==1L )       # 1 = completed
  #     out$outside = which( Sflag[]==2L )    # 2 = oustide bounds(if any)
  #     out$shallow = which( Sflag[]==3L )    # 3 = depth shallower than p$depth.filter (if it exists .. z is a covariate)
  #     out$rangeissue = which( Sflag[]==4L ) # 4=range not ok,
  #     out$nodata = which( Sflag[]==5L )     # 5=skipped due to insufficient data,
  #     out$variogramerror = which( Sflag[]==6L ) # 6=skipped .. fast variogram did not work
  #     out$vrangeerror = which( Sflag[]==7L )     # 7=variogram estimated range not ok
  #     out$modelerror = which( Sflag[]==8L )     # 8=problem with prediction and/or modelling
  #     out$skipped = which( Sflag[] == 9L )   # 9 not completed due to a failed attempt
  # 
