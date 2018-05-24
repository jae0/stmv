
matern_find_range_fast = function( xy, z, nu=0.5, discretized_n=100, nbreaks = 13 ) {
  # nu=0.5 defaults to exponential
  #\\ NOTE:: the default parameterization is Wikipedia's paramterization:
  #\\ == Rasmussen, Carl Edward (2006) Gaussian Processes for Machine Learning
  #\\  sigma^2 * (2^{1-nu} / Gamma(nu) ) * (sqrt(2*nu) * ||x|| / phi)^{nu} * K_{nu}( sqrt(2*nu) * ||x|| / phi)
  #\\   where K_{nu} is the Bessel function with smooth nu and phi is known as the range parameter
  #\\ As usage sometimes is for high density data, aggregation to a coarse resolution of 'discretized_n' units along
  #\\ the smaller dimension  before computation.
  # -------------------------

  g = grid_fast(xy=xy, z=z, discretized_n=discretized_n, FUNC=mean, na.rm=TRUE)

  maxdist = g$drange/4   # begin with this (diagonal)

  # gives a fast stable empirical variogram using nl least squares

  XY = as.data.frame(g$xy / g$drange )  # drange keeps things smaller in value to avoid floating point issues
  names(XY) =  c("plon", "plat" ) # arbitrary

  vEm = gstat::variogram( g$z~1, locations=~plon+plat, data=XY, cutoff=maxdist/g$drange, width=maxdist/g$drange/nbreaks, cressie=FALSE ) # empirical variogram
  if (inherits(vEm, "try-error") ) return(NULL)
  vEm$dist0 = vEm$dist * g$drange
  vMod0 = gstat::vgm(psill=0.75, model="Mat", range=1, nugget=0.25, kappa=0.5 ) # starting model parameters
  vFitgs =  try( gstat::fit.variogram( vEm, vMod0, fit.kappa =TRUE, fit.sills=TRUE, fit.ranges=TRUE ) ) ## gstat's kappa is the Bessel function's "nu" smoothness parameter
  if (inherits(vFitgs, "try-error") )  return(NULL)

  phi = matern_phi2phi( mRange=vFitgs$range[2], mSmooth=vFitgs$kappa[2], parameterization_input="gstat", parameterization_output="stmv" ) * g$drange
  nu=vFitgs$kappa[2]
  range = matern_phi2distance( phi=phi, nu=nu  )
  # varSpatial=vFitgs$psill[2]
  # varObs=vFitgs$psill[1]

  if ( range < maxdist ) return( range )

  cnt = 0
  while ( cnt < 5  ) {
    maxdist = maxdist * 1.25
    if ( maxdist > drange ) {
      # message ( "Autocorrelation range greater than data range .. retrying a last time at max dist with more data")
      cnt = 7
      maxdist = drange
    }
    # message( "Range longer than distance cutoff ... retrying with a larger distance cutoff")
    cnt = cnt + 1

    vEm = gstat::variogram( g$z~1, locations=~plon+plat, data=XY, cutoff=maxdist/g$drange, width=maxdist/g$drange/nbreaks, cressie=FALSE ) # empirical variogram
    if (inherits(vEm, "try-error") ) return(NULL)
    vEm$dist0 = vEm$dist * g$drange
    vMod0 = gstat::vgm(psill=0.75, model="Mat", range=1, nugget=0.25, kappa=0.5 ) # starting model parameters
    vFitgs =  try( gstat::fit.variogram( vEm, vMod0, fit.kappa =TRUE, fit.sills=TRUE, fit.ranges=TRUE ) ) ## gstat's kappa is the Bessel function's "nu" smoothness parameter
    if (inherits(vFitgs, "try-error") )  return(NULL)

    phi = matern_phi2phi( mRange=vFitgs$range[2], mSmooth=vFitgs$kappa[2], parameterization_input="gstat", parameterization_output="stmv" ) * g$drange
    nu=vFitgs$kappa[2]
    range = matern_phi2distance( phi=phi, nu=nu  )

    if ( range < maxdist ) return( range )
  }

  return (range) # ie. not found
}
