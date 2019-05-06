
stmv__fft = function( p=NULL, dat=NULL, pa=NULL, nu=NULL, phi=NULL, variablelist=FALSE, ... ) {

  #\\ this is the core engine of stmv .. localised space (no-time) modelling interpolation
  #\\ note: time is not being modelled and treated independently
  #\\      .. you had better have enough data in each time slice
  #\\ first a low-pass filter as defined by p$stmv_lowpass_nu, p$stmv_lowpass_phi, then a simple covariance filter determined by nu,phi ;; fft (no lowpass)

  if (variablelist)  return( c() )

  sdTotal=sd(dat[,p$variable$Y], na.rm=T)

  x_r = range(dat[,p$variables$LOCS[1]])
  x_c = range(dat[,p$variables$LOCS[2]])

  pa_r = range(pa[,p$variables$LOCS[1]])
  pa_c = range(pa[,p$variables$LOCS[2]])

  nr = round( diff(x_r)/p$pres ) + 1
  nc = round( diff(x_c)/p$pres ) + 1

  x_plons = seq( x_r[1], x_r[2], length.out=nr )
  x_plats = seq( x_c[1], x_c[2], length.out=nc )

  x_locs = expand.grid( x_plons, x_plats ) # final output grid
  attr( x_locs , "out.attrs") = NULL
  names( x_locs ) = p$variables$LOCS

  dat$mean = NA
  pa$mean = NA
  pa$sd = sdTotal  # this is ignored with fft

  dx = dy = p$pres

  nr2 = 2 * nr
  nc2 = 2 * nc

  # constainer for spatial filters
  dgrid = make.surface.grid(list((1:nr2) * dx, (1:nc2) * dy))
  center = matrix(c((dx * nr2)/2, (dy * nc2)/2), nrow = 1, ncol = 2)

  mC = matrix(0, nrow = nr2, ncol = nc2)
  mC[nr, nc] = 1

  if (!exists("stmv_fft_filter",p) ) p$stmv_fft_filter="lowpass" # default in case of no specification

  if ( p$stmv_fft_filter == "lowpass") {
    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", range=p$stmv_lowpass_phi, nu=p$stmv_lowpass_nu )
    sp.covar.surf = as.surface(dgrid, c(sp.covar))$z
    sp.covar.kernel = fft(sp.covar.surf) / ( fft(mC) * nr2 * nc2 )
  }

  if (p$stmv_fft_filter %in% c("matern") ) {
    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", range=phi, nu=nu )
    sp.covar.surf = as.surface(dgrid, c(sp.covar))$z
    sp.covar.kernel = fft(sp.covar.surf) / ( fft(mC) * nr2 * nc2 )
  }

  if (p$stmv_fft_filter == "lowpass_matern") {
    # both ..
    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", range=p$stmv_lowpass_phi, nu=p$stmv_lowpass_nu )
    sp.covar.surf = as.surface(dgrid, c(sp.covar))$z
    sp.covar2 = stationary.cov( dgrid, center, Covariance="Matern", range=phi, nu=nu )
    sp.covar.surf2 = as.surface(dgrid, c(sp.covar2))$z
    sp.covar.kernel = {fft(sp.covar.surf)/ ( fft(mC) * nr2 * nc2 )} * {fft(sp.covar.surf2)/ ( fft(mC) * nr2 * nc2 )
  } }

  sp.covar = sp.covar2 = sp.covar.surf = sp.covar.surf2 = dgrid = center = mC = NULL

  xi =1:nrow(dat)  # all data as p$nt==1
  pa_i = 1:nrow(pa)
  origin=c(x_r[1], x_c[1])
  res=c(p$pres, p$pres)

  zz = matrix(1:(nr*nc), nrow = nr, ncol = nc)


  for ( ti in 1:p$nt ) {

    if ( exists("TIME", p$variables) ) {
      xi = which( dat[ , p$variables$TIME] == p$prediction.ts[ti] )
      pa_i = which( pa[, p$variables$TIME] == p$prediction.ts[ti])
    }

    # map of row, col indices of input data in the new (output) coordinate system
    if (length(xi) < 5 ) {
      # print( ti)
      next()
    }

    x_id = array_map( "xy->2", coords=dat[xi,p$variables$LOCS], origin=origin, res=res )

    u = as.image( dat[xi,p$variables$Y], ind=as.matrix( x_id), na.rm=TRUE, nx=nr, ny=nc )

    mN = matrix(0, nrow = nr2, ncol = nc2)
    mN[1:nr,1:nc] = u$weights
    mN[!is.finite(mN)] = 0

    mY = matrix(0, nrow = nr2, ncol = nc2)
    mY[1:nr,1:nc] = u$z
    mY[!is.finite(mY)] = 0

    u = NULL

    # low pass filter based upon a global nu,phi .. remove high freq variation
    fN = Re(fft(fft(mN) * sp.covar.kernel, inverse = TRUE))[1:nr,1:nc]
    fY = Re(fft(fft(mY) * sp.covar.kernel, inverse = TRUE))[1:nr,1:nc]

    mY = mN = NULL

    Z = fY/fN
    fY = fN = NULL

    Z_i = array_map( "xy->2", coords=pa[pa_i,p$variables$LOCS], origin=origin, res=res )

    # bounds check: make sure predictions exist
    Z_i_test = NULL
    Z_i_test = which( Z_i[,1]<1 | Z_i[,2]<1  | Z_i[,1] > nr | Z_i[,2] > nc )

    if (length(Z_i_test) > 0) {
      keep = zz[ Z_i[-Z_i_test,] ]
      pa$mean[pa_i[keep]] = Z[keep]
    } else {
      pa$mean[pa_i] = Z[Z_i]
    }

    # pa$sd[pa_i] = NA  ## fix as NA
    Z = NULL
  }


  stmv_stats = list( sdTotal=sdTotal, rsquared=NA, ndata=nrow(dat) ) # must be same order as p$statsvars

  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, stmv_stats=stmv_stats ) )
}
