
stmv__fft = function( p=NULL, dat=NULL, pa=NULL, nu=NULL, phi=NULL, variablelist=FALSE, tol=1e-12, weights = 1, ... ) {

  #\\ this is the core engine of stmv .. localised space (no-time) modelling interpolation
  #\\ note: time is not being modelled and treated independently
  #\\      .. you had better have enough data in each time slice
  #\\ first a low-pass filter as defined by p$stmv_lowpass_nu, p$stmv_lowpass_phi,
  #\\ then a simple covariance filter determined by nu,phi ;; fft (no lowpass)
  #\\ based upon fields::image.smooth

  if (variablelist)  return( c() )

  params = list(...)

  sdTotal = sd(dat[,p$variable$Y], na.rm=T)

  # nr .. x/plon
  # nc .. y/plat
  # pa_r = range(pa[,p$variables$LOCS[1]])
  # pa_c = range(pa[,p$variables$LOCS[2]])

  x_r = range(dat[,p$variables$LOCS[1]])
  x_c = range(dat[,p$variables$LOCS[2]])

  nr = round( diff(x_r)/p$pres ) + 1
  nc = round( diff(x_c)/p$pres ) + 1

  # final output grid
  x_locs = expand.grid(
    seq( x_r[1], x_r[2], length.out=nr ),
    seq( x_c[1], x_c[2], length.out=nc )
  )
  attr( x_locs , "out.attrs") = NULL
  names( x_locs ) = p$variables$LOCS

  dat$mean = NA
  pa$mean = NA
  pa$sd = sdTotal  # this is ignored with fft

  dx = dy = p$pres

  nr2 = 2 * nr
  nc2 = 2 * nc

  # constainer for spatial filters
  grid.list = list((1:nr2) * dx, (1:nc2) * dy)
  dgrid = as.matrix(expand.grid(grid.list))
  dimnames(dgrid) = list(NULL, names(grid.list))
  attr(dgrid, "grid.list") = grid.list

  center = matrix(c((dx * nr), (dy * nc)), nrow = 1, ncol = 2)

  mC = matrix(0, nrow = nr2, ncol = nc2)
  mC[nr, nc] = 1
  fft_mC = fft(mC) * nr2 * nc2
  mC = NULL

  if (!exists("stmv_fft_filter",p) ) p$stmv_fft_filter="lowpass" # default in case of no specification

  if ( p$stmv_fft_filter == "lowpass") {
    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=p$stmv_lowpass_phi, smoothness=p$stmv_lowpass_nu )
    sp.covar = as.surface(dgrid, c(sp.covar))$z
    sp.covar.kernel = fft(sp.covar) / fft_mC
  }

  if (p$stmv_fft_filter %in% c("matern") ) {
    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=phi, smoothness=nu )
    sp.covar = as.surface(dgrid, c(sp.covar))$z
    sp.covar.kernel = fft(sp.covar) / fft_mC
  }

  if (p$stmv_fft_filter == "lowpass_matern") {
    # both ..
    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=p$stmv_lowpass_phi, smoothness=p$stmv_lowpass_nu )
    sp.covar = as.surface(dgrid, c(sp.covar))$z
    sp.covar2 = stationary.cov( dgrid, center, Covariance="Matern", theta=phi, smoothness=nu )
    sp.covar2 = as.surface(dgrid, c(sp.covar2))$z
    sp.covar.kernel = {fft(sp.covar) / fft_mC } * {fft(sp.covar2)/ fft_mC }
  }


  if (p$stmv_fft_filter %in% c("matern_constant") ) {
    # both ..
    if (!exists("stmv_constant_nu", p)) {
      stmv_constant_nu = 0.5
    } else {
      stmv_constant_nu = p$stmv_constant_nu
    }

    if (!exists("stmv_constant_phi", p)) {
      stmv_constant_phi = p$pres
    } else {
      stmv_constant_phi = p$stmv_constant_phi
    }

    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=stmv_constant_phi, smoothness=stmv_constant_nu )
    sp.covar = as.surface(dgrid, c(sp.covar))$z
    sp.covar.kernel = fft(sp.covar) / fft_mC
  }

  if (p$stmv_fft_filter == "lowpass_matern_constant") {
    # both ..
    if (!exists("stmv_constant_nu", p)) {
      stmv_constant_nu = 0.5
    } else {
      stmv_constant_nu = p$stmv_constant_nu
    }

    if (!exists("stmv_constant_phi", p)) {
      stmv_constant_phi = p$pres
    } else {
      stmv_constant_phi = p$stmv_constant_phi
    }

    if (!exists("stmv_lowpass_nu", p)) {
      stmv_lowpass_nu = 0.5
    } else {
      stmv_lowpass_nu = p$stmv_lowpass_nu
    }

    if (!exists("stmv_lowpass_phi", p)) {
      stmv_lowpass_phi = p$pres / 2
    } else {
      stmv_lowpass_phi = p$stmv_lowpass_phi
    }

    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=stmv_lowpass_phi, smoothness=stmv_lowpass_nu )
    sp.covar = as.surface(dgrid, c(sp.covar))$z
    sp.covar2 = stationary.cov( dgrid, center, Covariance="Matern", theta=stmv_constant_phi, smoothness=stmv_constant_nu )
    sp.covar2 = as.surface(dgrid, c(sp.covar2))$z
    sp.covar.kernel = {fft(sp.covar) / fft_mC } * {fft(sp.covar2)/ fft_mC }
  }


  if (p$stmv_fft_filter == "normal_kernel") {
      xi = seq(-(nr - 1), nr, 1) * dx / phi
      yi = seq(-(nc - 1), nc, 1) * dy / phi
      dd = ((matrix(xi, nr2, nc2)^2 + matrix(yi, nr2, nc2, byrow = TRUE)^2))  # squared distances
      # double.exp: An R function that takes as its argument the _squared_
      # distance between two points divided by the bandwidth. The
      # default is exp( -abs(x)) yielding a normal kernel
      kk = double.exp(dd)
      mK = matrix(kk, nrow = nr2, ncol = nc2)
      W = fft(mK) / fft_mC
      sp.covar.kernel = W/(nr2 * nc2)  # kernal weights
  }

  sp.covar = sp.covar2 = dgrid = center =  NULL

  dat_i =1:nrow(dat)  # all data as p$nt==1
  pa_i = 1:nrow(pa)
  origin=c(x_r[1], x_c[1])
  res=c(p$pres, p$pres)

  zz = matrix(1:(nr*nc), nrow = nr, ncol = nc)
  mY = matrix(0, nrow = nr2, ncol = nc2)
  mN = matrix(0, nrow = nr2, ncol = nc2)

  for ( ti in 1:p$nt ) {

    if ( exists("TIME", p$variables) ) {
      dat_ = which( dat[ , p$variables$TIME] == p$prediction.ts[ti] )
      pa_i = which( pa [ , p$variables$TIME] == p$prediction.ts[ti] )
      # map of row, col indices of input data in the new (output) coordinate system
      if (length(dat_i) < 5 ) {
        # print( ti)
        next()
      }
    }

    u = as.image(
      dat[dat_i,p$variables$Y],
      ind=as.matrix(array_map( "xy->2", coords=dat[dat_i,p$variables$LOCS], origin=origin, res=res )),
      na.rm=TRUE,
      nx=nr,
      ny=nc
    )
    mY[1:nr,1:nc] = u$z * weights
    mY[!is.finite(mY)] = 0
    mN[1:nr,1:nc] = ifelse(!is.na(u$z), weights, 0)
    u =NULL

    fY = Re(fft(fft(mY) * sp.covar.kernel, inverse = TRUE))[1:nr, 1:nc]
    fN = Re(fft(fft(mN) * sp.covar.kernel, inverse = TRUE))[1:nr, 1:nc]
    Z = ifelse((fN > tol), (fY/fN), NA)
    fY = fN = NULL

    # low pass filter based upon a global nu,phi .. remove high freq variation
    Z_i = array_map( "xy->2", coords=pa[pa_i,p$variables$LOCS], origin=origin, res=res )

    # bounds check: make sure predictions exist
    Z_i_test = NULL
    Z_i_test = which( Z_i[,1]<1 | Z_i[,2]<1  | Z_i[,1] > nr | Z_i[,2] > nc )

    if (length(Z_i_test) > 0) {
      keep = zz[ Z_i[-Z_i_test,] ]
      pa$mean[pa_i[keep]] = Z[keep]
      keep = NULL
    } else {
      pa$mean[pa_i] = Z[Z_i]
    }

    # pa$sd[pa_i] = NA  ## fix as NA
    Z = Z_i = Z_i_test = NULL
  }

  stmv_stats = list( sdTotal=sdTotal, rsquared=NA, ndata=nrow(dat) ) # must be same order as p$statsvars

  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, stmv_stats=stmv_stats ) )


}
