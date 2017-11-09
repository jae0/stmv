
matern_find_range_fast = function( xy, z, nu=0.5 ) {
  # nu=0.5 defaults to exponential
  #\\ NOTE:: the default parameterization is Wikipedia's paramterization:
  #\\ == Rasmussen, Carl Edward (2006) Gaussian Processes for Machine Learning
  #\\  sigma^2 * (2^{1-nu} / Gamma(nu) ) * (sqrt(2*nu) * ||x|| / phi)^{nu} * K_{nu}( sqrt(2*nu) * ||x|| / phi)
  #\\   where K_{nu} is the Bessel function with smooth nu and phi is known as the range parameter  
  # -------------------------

  xr = range( xy[,1], na.rm=TRUE )
  yr = range( xy[,2], na.rm=TRUE )
  drange = sqrt( diff( xr)^2 + diff( yr)^2  )

  maxdist = drange/5   # begin with this
  nbreaks = 13

  # gives a fast stable empirical variogram using nl least squares
  require( fields ) 

  autocorrelation_range = NA
  vario = vgram( xy, z, dmax=maxdist, N=nbreaks)
  vx=vario$centers
  vg=vario$stats["mean",]
  vo = ecmei_variogram_optimization( vx=vx, vg=vg, nu=nu ) 
  if ( vo$summary$range_ok ) autocorrelation_range = vo$summary$range

  cnt = 0
  while ( cnt < 7  ) {
    maxdist = maxdist * 1.25
    if ( maxdist > drange ) {
      # message ( "Autocorrelation range greater than data range .. retrying a last time at max dist with more data")
      cnt = 7
      maxdist = drange
    }
    # message( "Range longer than distance cutoff ... retrying with a larger distance cutoff")
    cnt = cnt + 1
    vario = vgram( xy, z, dmax=maxdist, N=nbreaks)
    vx = vario$centers
    vg = vario$stats["mean",]
    vo = ecmei_variogram_optimization( vx=vx, vg=vg, nu ) 
    if ( vo$summary$range_ok ) {
      autocorrelation_range = vo$summary$range
      break()
    }
  }
  return (autocorrelation_range) # ie. not found
}


