

stmv_predictionarea_lattice = function(p, sloc, windowsize.half  ) {

  pa = stmv_predictionarea_space( p=p, sloc=sloc, windowsize.half=windowsize.half )
  if (is.null(pa)) return(NULL)

  if ( exists("TIME", p$stmv_variables) )  pa = stmv_predictionarea_time( p=p, pa=pa )

  rownames(pa) = NULL
  return(pa)

}
