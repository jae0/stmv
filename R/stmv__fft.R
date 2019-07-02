
stmv__fft = function( p=NULL, dat=NULL, pa=NULL, nu=NULL, phi=NULL, variablelist=FALSE, eps=1e-9, ... ) {

  #\\ this is the core engine of stmv .. localised space (no-time) modelling interpolation
  #\\ note: time is not being modelled and treated independently
  #\\      .. you had better have enough data in each time slice
  #\\ first a low-pass filter as defined by p$stmv_lowpass_nu, p$stmv_lowpass_phi,
  #\\ then a simple covariance filter determined by nu,phi ;; fft (no lowpass)
  #\\ based upon fields::image.smooth and setup.image.smooth
  #\\ now that the local area of interest is stationary, we use local convolutions of highly autocorrelated (short-range)
  #\\ determined by p
  #\\ using fftw interace fftwtools::fftw2d

# The following documentation is from fields::smooth.2d:
  # The irregular locations are first discretized to a regular grid (
  #    using as.image) then a 2d- FFT is used to compute a
  #    Nadaraya-Watson type kernel estimator. Here we take advantage of
  #    two features. The kernel estimator is a convolution and by padding
  #    the regular by zeroes where data is not obsevred one can sum the
  #    kernel over irregular sets of locations.  A second convolutions to
  #    find the normalization of the kernel weights.

  #    The kernel function is specified by an function that should
  #    evaluate with the kernel for two matrices of locations. Assume
  #    that the kernel has the form: K( u-v) for two locations u and v.
  #    The function given as the argument to cov.function should have the
  #    call myfun( x1,x2) where x1 and x2 are matrices of 2-d locations
  #    if nrow(x1)=m and nrow( x2)=n then this function should return a
  #    mXn matrix where the (i,j) element is K( x1[i,]- x2[j,]).
  #    Optional arguments that are included in the ... arguments are
  #    passed to this function when it is used. The default kernel is the
  #    Gaussian and the argument theta is the bandwidth. It is easy to
  #    write other other kernels, just use Exp.cov.simple as a template.

  if (0) {
    varObs = S[Si, i_sdObs]^2
    varSpatial = S[Si, i_sdSpatial]^2
    sloc = Sloc[Si,]
    distance = localrange
    eps = 1e-9
  }

  if (variablelist)  return( c() )

  params = list(...)

  sdTotal = sd(dat[,p$variable$Y], na.rm=T)

  # nr .. x/plon
  # nc .. y/plat
  # pa_r = range(pa[,p$variables$LOCS[1]])
  # pa_c = range(pa[,p$variables$LOCS[2]])

  dx = p$pres
  dy = p$pres

  # system size
  #nr = nx
  #nc = ny

  x_r = range(dat[,p$variables$LOCS[1]])
  x_c = range(dat[,p$variables$LOCS[2]])

  rr = diff(x_r)
  rc = diff(x_c)

  nr = floor( rr/dx ) + 1
  nc = floor( rc/dy ) + 1

  dr = rr/(nr-1)
  dc = rc/(nc-1)

  nr2 = 2 * nr
  nc2 = 2 * nc

  # approx sa associated with each datum
  sa = rr * rc
  d_sa = sa/nrow(dat) # sa associated with each datum
  d_length = sqrt( d_sa/pi )  # sa = pi*l^2  # characteristic length scale
  theta.Taper = d_length * p$stmv_fft_taper_factor


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

  # constainer for spatial filters
  grid.list = list((1:nr2) * dx, (1:nc2) * dy)
  dgrid = as.matrix(expand.grid(grid.list))
  dimnames(dgrid) = list(NULL, names(grid.list))
  attr(dgrid, "grid.list") = grid.list
  grid.list = NULL

  center = matrix(c((dx * nr), (dy * nc)), nrow = 1, ncol = 2)

  mC = matrix(0, nrow = nr2, ncol = nc2)
  mC[nr, nc] = 1
  mC_fft = fftwtools::fftw2d(mC)
  mC = NULL

  if (!exists("stmv_fft_filter",p) ) p$stmv_fft_filter="lowpass" # default in case of no specification

  origin = c(x_r[1], x_c[1])
  resolution = c(dx, dy)

  fY = matrix(0, nrow = nr2, ncol = nc2)
  fN = matrix(0, nrow = nr2, ncol = nc2)

  for ( ti in 1:p$nt ) {

    if ( exists("TIME", p$variables) ) {
      xi   = which( dat[ , p$variables$TIME] == p$prediction.ts[ti] )
      pa_i = which( pa[, p$variables$TIME] == p$prediction.ts[ti] )
      if (length(xi) < 5 ) {
        # print( ti)
        next()
      }
    } else {
      xi   = 1:nrow(dat) # all data as p$nt==1
      pa_i = 1:nrow(pa)
    }
    # bounds check: make sure predictions exist
    zmean = mean(dat[xi, p$variables$Y], na.rm=TRUE)
    zsd = sd(dat[xi, p$variables$Y], na.rm=TRUE)
    zvar = zsd^2
    z = (dat[xi, p$variables$Y] - zmean) / zsd # zscore -- making it mean 0 removes the DC component

    u = as.image(
      z,
      ind=as.matrix(array_map( "xy->2", coords=dat[xi, p$variables$LOCS], origin=origin, res=resolution )),
      na.rm=TRUE,
      nx=nr,
      ny=nc
    )
    # surface (u)

    # See explanation:  https://en.wikipedia.org/wiki/Autocorrelation#Efficient_computation
    # Robertson, C., & George, S. C. (2012). Theory and practical recommendations for autocorrelation-based image
    # correlation spectroscopy. Journal of biomedical optics, 17(8), 080801. doi:10.1117/1.JBO.17.8.080801
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3414238/

    #  Nadaraya/Watson normalization for missing values s
    fY[] = 0
    fN[] = 0

    fY[1:nr,1:nc] = u$z
    fN[1:nr,1:nc] = u$weights

    u =NULL

    fY[!is.finite(fY)] = 0
    fN[!is.finite(fN)] = 0

    fY = fftwtools::fftw2d(fY)
    fN = fftwtools::fftw2d(fN)

    # fY * Conj(fY) == power spectra
    fY = Re( fftwtools::fftw2d( fY * Conj(fY), inverse=TRUE)  ) # autocorrelation (amplitude)
    fN = Re( fftwtools::fftw2d( fN * Conj(fN), inverse=TRUE)  ) # autocorrelation (amplitude) correction

    X = ifelse(( fN > eps), (fY / fN), NA) # autocorrelation (amplitude)

    # fftshift
    X = rbind( X[((nr+1):nr2), (1:nc2)], X[(1:nr), (1:nc2)] )  # swap_up_down
    X = cbind(X[1:nr2, ((nc+1):nc2)], X[1:nr2, 1:nc])  # swap_left_right

    # radial representation
    xy = expand.grid( x = c(-(nr-1):0, 0:(nr-1)) * dr,  y = c(-(nc-1):0, 0:(nc-1)) * dc )
    distances = sqrt(xy$x^2 + xy$y^2)
    dmax = max(distances, na.rm=TRUE ) * 0.4  # approx nyquist distance (<0.5 as corners exist)
    breaks = seq( 0, dmax, length.out=nr)
    db = breaks[2] - breaks[1]
    # angles = atan2( xy$y, xy$x )  # not used

    zz = cut( distances, breaks=c(breaks, max(breaks)+db), label=breaks+(db/2) )
    distances = NULL
    xy = NULL

    vgm = as.data.frame.table(tapply( X=X, INDEX=zz, FUN=mean, na.rm=TRUE ))
    names(vgm) = c("distances", "ac")
    X = NULL
    zz = NULL
    gc()

    vgm$distances = as.numeric( as.character(vgm$distances))
    vgm$sv =  zvar * (1-vgm$ac^2) # each sv are truly orthogonal

  # plot(ac ~ distances, data=vgm   )
  # plot(sv ~ distances, data=vgm   )

   # interpolated surface
   # constainer for spatial filters
    uu = which( (vgm$distances < dmax ) & is.finite(vgm$sv) )  # dmax ~ Nyquist freq
    fit = try( stmv_variogram_optimization( vx=vgm$distances[uu], vg=vgm$sv[uu], plotvgm=FALSE,
      stmv_internal_scale=dmax*0.75, cor=p$stmv_range_correlation ))
    # out$fit = fit

    local_phi = phi
    local_nu = nu

    if ( !inherits(fit, "try-error") ) {
      if (exists("summary", fit)) {
        if ( exists("nu", fit$summary )) {
          if ( is.finite( c(fit$summary$nu)) ) {
            local_nu = fit$summary$nu
          }
        }
        if ( exists("phi", fit$summary )) {
          if ( is.finite( c(fit$summary$phi)) ) {
            local_phi = fit$summary$phi
          }
        }
      }
    }

    if ( p$stmv_fft_filter == "lowpass") {
      theta = p$stmv_lowpass_phi
      sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=theta, smoothness=p$stmv_lowpass_nu )
      sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
      sp.covar = fftwtools::fftw2d(sp.covar) / mC_fft
    }

    if (p$stmv_fft_filter %in% c("matern") ) {
      sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=local_phi, smoothness=local_nu )
      sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
      sp.covar = fftwtools::fftw2d(sp.covar) / mC_fft
    }

    if (p$stmv_fft_filter == "lowpass_matern") {
      # both ..
      sp.covar.lowpass = stationary.cov( dgrid, center, Covariance="Matern", theta=p$stmv_lowpass_phi, smoothness=p$stmv_lowpass_nu )
      sp.covar.lowpass = as.surface(dgrid, c(sp.covar.lowpass))$z / (nr2*nc2)
      sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=local_phi, smoothness=local_nu )
      sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
      sp.covar = ( fftwtools::fftw2d(sp.covar.lowpass) / mC_fft ) * (fftwtools::fftw2d(sp.covar)/ mC_fft )
      sp.covar.lowpass = NULL
    }

    if (p$stmv_fft_filter == "matern_tapered") {
      #theta.Taper = matern_phi2distance( phi=local_phi, nu=nu, cor=p$stmv_fft_taper_factor )
      sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=local_phi, smoothness=local_nu,
        Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)
      sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
      sp.covar = fftwtools::fftw2d(sp.covar) / mC_fft
  }

    if (p$stmv_fft_filter == "lowpass_matern_tapered") {
      sp.covar.lowpass = stationary.cov( dgrid, center, Covariance="Matern", theta=p$stmv_lowpass_phi, smoothness=p$stmv_lowpass_nu )
      sp.covar.lowpass = as.surface(dgrid, c(sp.covar.lowpass))$z / (nr2*nc2)
      # theta.Taper = matern_phi2distance( phi=local_phi, nu=nu, cor=p$stmv_fft_taper_factor )
      sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=local_phi, smoothness=local_nu,
        Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)
      sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
      sp.covar = (fftwtools::fftw2d(sp.covar.lowpass) / mC_fft ) * ( fftwtools::fftw2d(sp.covar)/ mC_fft )
      sp.covar.lowpass = NULL
    }


    if (p$stmv_fft_filter == "normal_kernel") {
      theta = matern_phi2distance( phi=local_phi, nu=local_nu, cor=p$stmv_range_correlation )
      xseq = seq(-(nr - 1), nr, 1) * dx / theta
      yseq = seq(-(nc - 1), nc, 1) * dy / theta
      dd = ((matrix(xseq, nr2, nc2)^2 + matrix(yseq, nr2, nc2, byrow = TRUE)^2))  # squared distances
      # double.exp: An R function that takes as its argument the _squared_
      # distance between two points divided by the bandwidth. The
      # default is exp( -abs(x)) yielding a normal kernel
      sp.covar = matrix( double.exp(dd), nrow = nr2, ncol = nc2)
      sp.covar = fftwtools::fftw2d(sp.covar) / mC_fft / (nr2*nc2)# kernal weights
    }

    fY = Re( fftwtools::fftw2d( sp.covar * fY, inverse = TRUE))[1:nr, 1:nc]
    fN = Re( fftwtools::fftw2d( sp.covar * fN, inverse = TRUE))[1:nr, 1:nc]
    sp.covar = NULL
    Z = ifelse((fN > eps), (fY/fN), NA)
    Z[!is.finite(Z)] = NA
    Z = Z * zsd + zmean # revert to input scale
    if (0) {
      dev.new()
      surface(list(x=c(1:nr)*dr, y=c(1:nc)*dc, z=Z), xaxs="r", yaxs="r")
    }
    # out$Z = Z

    Z_i = array_map( "xy->2", coords=pa[pa_i, p$variables$LOCS], origin=origin, res=resolution )
    tokeep = which( Z_i[,1] >= 1 & Z_i[,2] >= 1  & Z_i[,1] <= nr & Z_i[,2] <= nc )
    if (length(tokeep) < 1) next()
    Z_i = Z_i[tokeep,]

    pa$mean[pa_i[tokeep]] = Z[Z_i]
      # pa$sd[pa_i] = NA  ## fix as NA
    Z = Z_i = Z_i_test = NULL
  }

  stmv_stats = list( sdTotal=sdTotal, rsquared=NA, ndata=nrow(dat) ) # must be same order as p$statsvars

  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=TRUE) , aspect="iso" )

  return( list( predictions=pa, stmv_stats=stmv_stats ) )


  # -------------------------------


  if (0) {

    # copy of original .. easier to understand

    stmv__fft = function( p=NULL, dat=NULL, pa=NULL, nu=NULL, phi=NULL, variablelist=FALSE, eps=1e-9, ... ) {

      #\\ this is the core engine of stmv .. localised space (no-time) modelling interpolation
      #\\ note: time is not being modelled and treated independently
      #\\      .. you had better have enough data in each time slice
      #\\ first a low-pass filter as defined by p$stmv_lowpass_nu, p$stmv_lowpass_phi,
      #\\ then a simple covariance filter determined by nu,phi ;; fft (no lowpass)
      #\\ based upon fields::image.smooth and setup.image.smooth
      #\\ now that the local area of interest is stationary, we use local convolutions of highly autocorrelated (short-range)
      #\\ determined by p
      #\\ using fftw interace fftwtools::fftw2d

    # The following documentation is from fields::smooth.2d:
      # The irregular locations are first discretized to a regular grid (
      #    using as.image) then a 2d- FFT is used to compute a
      #    Nadaraya-Watson type kernel estimator. Here we take advantage of
      #    two features. The kernel estimator is a convolution and by padding
      #    the regular by zeroes where data is not obsevred one can sum the
      #    kernel over irregular sets of locations.  A second convolutions to
      #    find the normalization of the kernel weights.

      #    The kernel function is specified by an function that should
      #    evaluate with the kernel for two matrices of locations. Assume
      #    that the kernel has the form: K( u-v) for two locations u and v.
      #    The function given as the argument to cov.function should have the
      #    call myfun( x1,x2) where x1 and x2 are matrices of 2-d locations
      #    if nrow(x1)=m and nrow( x2)=n then this function should return a
      #    mXn matrix where the (i,j) element is K( x1[i,]- x2[j,]).
      #    Optional arguments that are included in the ... arguments are
      #    passed to this function when it is used. The default kernel is the
      #    Gaussian and the argument theta is the bandwidth. It is easy to
      #    write other other kernels, just use Exp.cov.simple as a template.

      if (0) {
        varObs = S[Si, i_sdObs]^2
        varSpatial = S[Si, i_sdSpatial]^2
        sloc = Sloc[Si,]
        distance = localrange
        eps = 1e-9
      }

      if (variablelist)  return( c() )

      params = list(...)

      sdTotal = sd(dat[,p$variable$Y], na.rm=T)

      # nr .. x/plon
      # nc .. y/plat
      # pa_r = range(pa[,p$variables$LOCS[1]])
      # pa_c = range(pa[,p$variables$LOCS[2]])

      dx = p$pres
      dy = p$pres

      # system size
      #nr = nx
      #nc = ny

      x_r = range(dat[,p$variables$LOCS[1]])
      x_c = range(dat[,p$variables$LOCS[2]])

      rr = diff(x_r)
      rc = diff(x_c)

      nr = floor( rr/dx ) + 1
      nc = floor( rc/dy ) + 1

      dr = rr/(nr-1)
      dc = rc/(nc-1)

      nr2 = 2 * nr
      nc2 = 2 * nc

      # approx sa associated with each datum
      sa = rr * rc
      d_sa = sa/nrow(dat) # sa associated with each datum
      d_length = sqrt( d_sa/pi )  # sa = pi*l^2  # characteristic length scale
      theta.Taper = d_length * p$stmv_fft_taper_factor


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

      # constainer for spatial filters
      grid.list = list((1:nr2) * dx, (1:nc2) * dy)
      dgrid = as.matrix(expand.grid(grid.list))
      dimnames(dgrid) = list(NULL, names(grid.list))
      attr(dgrid, "grid.list") = grid.list
      grid.list = NULL

      center = matrix(c((dx * nr), (dy * nc)), nrow = 1, ncol = 2)

      mC = matrix(0, nrow = nr2, ncol = nc2)
      mC[nr, nc] = 1

      if (!exists("stmv_fft_filter",p) ) p$stmv_fft_filter="lowpass" # default in case of no specification

      origin = c(x_r[1], x_c[1])
      resolution = c(dx, dy)

      mY = matrix(0, nrow = nr2, ncol = nc2)
      mN = matrix(0, nrow = nr2, ncol = nc2)

      for ( ti in 1:p$nt ) {

        if ( exists("TIME", p$variables) ) {
          xi   = which( dat[ , p$variables$TIME] == p$prediction.ts[ti] )
          pa_i = which( pa[, p$variables$TIME] == p$prediction.ts[ti] )
          if (length(xi) < 5 ) {
            # print( ti)
            next()
          }
        } else {
          xi   = 1:nrow(dat) # all data as p$nt==1
          pa_i = 1:nrow(pa)
        }
        # bounds check: make sure predictions exist
        zmean = mean(dat[xi, p$variables$Y], na.rm=TRUE)
        zsd = sd(dat[xi, p$variables$Y], na.rm=TRUE)
        zvar = zsd^2
        z = (dat[xi, p$variables$Y] - zmean) / zsd # zscore -- making it mean 0 removes the DC component

        u = as.image(
          z,
          ind=as.matrix(array_map( "xy->2", coords=dat[xi, p$variables$LOCS], origin=origin, res=resolution )),
          na.rm=TRUE,
          nx=nr,
          ny=nc
        )
        # surface (u)

        # See explanation:  https://en.wikipedia.org/wiki/Autocorrelation#Efficient_computation
        # Robertson, C., & George, S. C. (2012). Theory and practical recommendations for autocorrelation-based image
        # correlation spectroscopy. Journal of biomedical optics, 17(8), 080801. doi:10.1117/1.JBO.17.8.080801
        # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3414238/

        mY[1:nr,1:nc] = u$z
        mY[!is.finite(mY)] = 0

        #  Nadaraya/Watson normalization for missing values s
        mN[1:nr,1:nc] = u$weights
        mN[!is.finite(mN)] = 0

        u =NULL

        fY = fftwtools::fftw2d(mY)
        fN = fftwtools::fftw2d(mN)

        # fY * Conj(fY) == power spectra
        ii = Re( fftwtools::fftw2d( fY * Conj(fY), inverse=TRUE)  ) # autocorrelation (amplitude)
        jj = Re( fftwtools::fftw2d( fN * Conj(fN), inverse=TRUE)  ) # autocorrelation (amplitude) correction

        X = ifelse(( jj > eps), (ii / jj), NA) # autocorrelation (amplitude)
        ii = NULL
        jj = NULL

        # fftshift
        X = rbind( X[((nr+1):nr2), (1:nc2)], X[(1:nr), (1:nc2)] )  # swap_up_down
        X = cbind(X[1:nr2, ((nc+1):nc2)], X[1:nr2, 1:nc])  # swap_left_right

        # radial representation
        xy = expand.grid( x = c(-(nr-1):0, 0:(nr-1)) * dr,  y = c(-(nc-1):0, 0:(nc-1)) * dc )
        distances = sqrt(xy$x^2 + xy$y^2)
        dmax = max(distances, na.rm=TRUE ) * 0.4  # approx nyquist distance (<0.5 as corners exist)
        breaks = seq( 0, dmax, length.out=nr)
        db = breaks[2] - breaks[1]
        # angles = atan2( xy$y, xy$x )  # not used

        zz = cut( distances, breaks=c(breaks, max(breaks)+db), label=breaks+(db/2) )
        distances = NULL
        xy = NULL

        vgm = as.data.frame.table(tapply( X=X, INDEX=zz, FUN=mean, na.rm=TRUE ))
        names(vgm) = c("distances", "ac")
        X = NULL
        zz = NULL
        gc()

        vgm$distances = as.numeric( as.character(vgm$distances))
        vgm$sv =  zvar * (1-vgm$ac^2) # each sv are truly orthogonal

      # plot(ac ~ distances, data=vgm   )
      # plot(sv ~ distances, data=vgm   )

      # interpolated surface
      # constainer for spatial filters
        uu = which( (vgm$distances < dmax ) & is.finite(vgm$sv) )  # dmax ~ Nyquist freq
        fit = try( stmv_variogram_optimization( vx=vgm$distances[uu], vg=vgm$sv[uu], plotvgm=FALSE,
          stmv_internal_scale=dmax*0.75, cor=p$stmv_range_correlation ))
        # out$fit = fit

        local_phi = phi
        local_nu = nu

        if ( !inherits(fit, "try-error") ) {
          if (exists("summary", fit)) {
            if ( exists("nu", fit$summary )) {
              if ( is.finite( c(fit$summary$nu)) ) {
                local_nu = fit$summary$nu
              }
            }
            if ( exists("phi", fit$summary )) {
              if ( is.finite( c(fit$summary$phi)) ) {
                local_phi = fit$summary$phi
              }
            }
          }
        }

        if ( p$stmv_fft_filter == "lowpass") {
          theta = p$stmv_lowpass_phi
          sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=theta, smoothness=p$stmv_lowpass_nu )
          sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
          sp.covar.kernel = fftwtools::fftw2d(sp.covar) / fftwtools::fftw2d(mC)
          sp.covar = NULL
        }

        if (p$stmv_fft_filter %in% c("matern") ) {
          sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=local_phi, smoothness=local_nu )
          sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
          sp.covar.kernel = fftwtools::fftw2d(sp.covar) / fftwtools::fftw2d(mC)
          sp.covar = NULL
        }

        if (p$stmv_fft_filter == "lowpass_matern") {
          # both ..
          sp.covar.lowpass = stationary.cov( dgrid, center, Covariance="Matern", theta=p$stmv_lowpass_phi, smoothness=p$stmv_lowpass_nu )
          sp.covar.lowpass = as.surface(dgrid, c(sp.covar.lowpass))$z / (nr2*nc2)
          sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=local_phi, smoothness=local_nu )
          sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
          sp.covar.kernel = ( fftwtools::fftw2d(sp.covar.lowpass) / fftwtools::fftw2d(mC) ) * (fftwtools::fftw2d(sp.covar)/ fftwtools::fftw2d(mC) )
          sp.covar.lowpass = NULL
          sp.covar = NULL
        }

        if (p$stmv_fft_filter == "matern_tapered") {
          #theta.Taper = matern_phi2distance( phi=local_phi, nu=nu, cor=p$stmv_fft_taper_factor )
          sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=local_phi, smoothness=local_nu,
            Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)
          sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
          sp.covar.kernel = fftwtools::fftw2d(sp.covar) / fftwtools::fftw2d(mC)
          sp.covar = NULL
      }

        if (p$stmv_fft_filter == "lowpass_matern_tapered") {
          sp.covar.lowpass = stationary.cov( dgrid, center, Covariance="Matern", theta=p$stmv_lowpass_phi, smoothness=p$stmv_lowpass_nu )
          sp.covar.lowpass = as.surface(dgrid, c(sp.covar.lowpass))$z / (nr2*nc2)
          # theta.Taper = matern_phi2distance( phi=local_phi, nu=nu, cor=p$stmv_fft_taper_factor )
          sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=local_phi, smoothness=local_nu,
            Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)
          sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
          sp.covar.kernel = (fftwtools::fftw2d(sp.covar.lowpass) / fftwtools::fftw2d(mC) ) * ( fftwtools::fftw2d(sp.covar)/ fftwtools::fftw2d(mC) )
          sp.covar.lowpass = NULL
          sp.covar = NULL
        }


        if (p$stmv_fft_filter == "normal_kernel") {
          theta = matern_phi2distance( phi=local_phi, nu=local_nu, cor=p$stmv_range_correlation )
          xseq = seq(-(nr - 1), nr, 1) * dx / theta
          yseq = seq(-(nc - 1), nc, 1) * dy / theta
          dd = ((matrix(xseq, nr2, nc2)^2 + matrix(yseq, nr2, nc2, byrow = TRUE)^2))  # squared distances
          # double.exp: An R function that takes as its argument the _squared_
          # distance between two points divided by the bandwidth. The
          # default is exp( -abs(x)) yielding a normal kernel
          sp.covar = matrix( double.exp(dd), nrow = nr2, ncol = nc2)
          sp.covar.kernel = fftwtools::fftw2d(sp.covar) / fftwtools::fftw2d(mC) / (nr2*nc2)# kernal weights
          sp.covar = NULL

        }

        ffY = Re( fftwtools::fftw2d( sp.covar.kernel * fY, inverse = TRUE))[1:nr, 1:nc]
        ffN = Re( fftwtools::fftw2d( sp.covar.kernel * fN, inverse = TRUE))[1:nr, 1:nc]
        Z = ifelse((ffN > eps), (ffY/ffN), NA)
        Z[!is.finite(Z)] = NA
        Z = Z * zsd + zmean # revert to input scale
        if (0) {
          dev.new()
          surface(list(x=c(1:nr)*dr, y=c(1:nc)*dc, z=Z), xaxs="r", yaxs="r")
        }
        # out$Z = Z
        fY = fN = ffY = ffN = NULL

        Z_i = array_map( "xy->2", coords=pa[pa_i, p$variables$LOCS], origin=origin, res=resolution )
        tokeep = which( Z_i[,1] >= 1 & Z_i[,2] >= 1  & Z_i[,1] <= nr & Z_i[,2] <= nc )
        if (length(tokeep) < 1) next()
        Z_i = Z_i[tokeep,]

        pa$mean[pa_i[tokeep]] = Z[Z_i]
          # pa$sd[pa_i] = NA  ## fix as NA
        Z = Z_i = Z_i_test = NULL
      }

      stmv_stats = list( sdTotal=sdTotal, rsquared=NA, ndata=nrow(dat) ) # must be same order as p$statsvars

      # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=TRUE) , aspect="iso" )

      return( list( predictions=pa, stmv_stats=stmv_stats ) )
    }


    # 1 D

    # useful tools and links

      # https://www.iro.umontreal.ca/~pift6080/H09/documents/
      # https://www.iro.umontreal.ca/~pift6080/H09/documents/eck_fft.pdf
      # https://riptutorial.com/r/example/14464/fourier-series
      # http://www.di.fc.ul.pt/~jpn/r/fourier/fourier.html


      x  =  1:4

      fft(x)
      fft(fft(x), inverse = TRUE)/length(x)

      ## Slow Discrete Fourier Transform (DFT) - e.g., for checking the formula
      fft0  =  function(z, inverse=FALSE) {
        n  =  length(z)
        if(n == 0) return(z)
        k  =  0:(n-1)
        ff  =  (if(inverse) 1 else -1) * 2*pi * 1i * k/n
        vapply(1:n, function(h) sum(z * exp(ff*(h-1))), complex(1))
      }

      fft0(x)


      fftshift  =  function(input_matrix, dim = -1) {

        rows  =  dim(input_matrix)[1]
        cols  =  dim(input_matrix)[2]

        swap_up_down  =  function(input_matrix) {
            rows_half  =  ceiling(rows/2)
            return(rbind(input_matrix[((rows_half+1):rows), (1:cols)], input_matrix[(1:rows_half), (1:cols)]))
        }

        swap_left_right  =  function(input_matrix) {
            cols_half  =  ceiling(cols/2)
            return(cbind(input_matrix[1:rows, ((cols_half+1):cols)], input_matrix[1:rows, 1:cols_half]))
        }

        if (dim == -1) {
            input_matrix  =  swap_up_down(input_matrix)
            return(swap_left_right(input_matrix))
        }
        else if (dim == 1) {
            return(swap_up_down(input_matrix))
        }
        else if (dim == 2) {
            return(swap_left_right(input_matrix))
        }
        else {
            stop("Invalid dimension parameter")
        }
      }

      ifftshift  =  function(input_matrix, dim = -1) {

        rows  =  dim(input_matrix)[1]
        cols  =  dim(input_matrix)[2]

        swap_up_down  =  function(input_matrix) {
            rows_half  =  floor(rows/2)
            return(rbind(input_matrix[((rows_half+1):rows), (1:cols)], input_matrix[(1:rows_half), (1:cols)]))
        }

        swap_left_right  =  function(input_matrix) {
            cols_half  =  floor(cols/2)
            return(cbind(input_matrix[1:rows, ((cols_half+1):cols)], input_matrix[1:rows, 1:cols_half]))
        }

        if (dim == -1) {
            input_matrix  =  swap_left_right(input_matrix)
            return(swap_up_down(input_matrix))
        }
        else if (dim == 1) {
            return(swap_up_down(input_matrix))
        }
        else if (dim == 2) {
            return(swap_left_right(input_matrix))
        }
        else {
            stop("Invalid dimension parameter")
        }
      }

      plot.fourier  =  function(fourier.series, f.0, ts) {
        w  =  2*pi*f.0
        trajectory  =  sapply(ts, function(t) fourier.series(t,w))
        plot(ts, trajectory, type="l", xlab="time", ylab="f(t)"); abline(h=0,lty=3)
      }

      convert.fft  =  function(cs, sample.rate=1) {
        # cs is the vector of complex points to convert
        cs  =  cs / length(cs) # normalize
        distance.center  =  function(c)signif( Mod(c),        4)
        angle            =  function(c)signif( 180*Arg(c)/pi, 3)
        df  =  data.frame(cycle    = 0:(length(cs)-1),
                          freq     = 0:(length(cs)-1) * sample.rate / length(cs),
                          strength = sapply(cs, distance.center),
                          delay    = sapply(cs, angle))
        df
      }

      convert.fft(fft(1:4))

      plot.frequency.spectrum  =  function(X.k, xlimits=c(0,length(X.k))) {
        plot.data   =  cbind(0:(length(X.k)-1), Mod(X.k))

        # TODO: why this scaling is necessary?
        plot.data[2:length(X.k),2]  =  2*plot.data[2:length(X.k),2]

        plot(plot.data, t="h", lwd=2, main="",
              xlab="Frequency (Hz)", ylab="Strength",
              xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
      }

      plot.harmonic  =  function(Xk, i, ts, acq.freq, color="red") {
        # Plot the i-th harmonic
        # Xk: the frequencies computed by the FFt
        #  i: which harmonic
        # ts: the sampling time points
        # acq.freq: the acquisition rate
        Xk.h  =  rep(0,length(Xk))
        Xk.h[i+1]  =  Xk[i+1] # i-th harmonic
        harmonic.trajectory  =  get.trajectory(Xk.h, ts, acq.freq=acq.freq)
        points(ts, harmonic.trajectory, type="l", col=color)
      }

      get.trajectory  =  function(X.k,ts,acq.freq) {
        # returns the x.n time series for a given time sequence (ts) and
        # a vector with the amount of frequencies k in the signal (X.k)
        N    =  length(ts)
        i    =  complex(real = 0, imaginary = 1)
        x.n  =  rep(0,N)           # create vector to tokeep the trajectory
        ks   =  0:(length(X.k)-1)
        for(n in 0:(N-1)) {       # compute each time point x_n based on freqs X.k
          x.n[n+1]  =  sum(X.k * exp(i*2*pi*ks*n/N)) / N
        }
        x.n * acq.freq
      }

      X.k  =  fft(c(4,0,0,0))                   # get amount of each frequency k

      time      =  4                            # measuring time interval (seconds)
      acq.freq  =  100                          # data acquisition frequency (Hz)
      ts   =  seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s)

      x.n  =  get.trajectory(X.k,ts,acq.freq)   # create time wave

      plot(ts,x.n,type="l",ylim=c(-2,4),lwd=2)
      abline(v=0:time,h=-2:4,lty=3); abline(h=0)

      plot.harmonic(X.k,1,ts,acq.freq,"red")
      plot.harmonic(X.k,2,ts,acq.freq,"green")
      plot.harmonic(X.k,3,ts,acq.freq,"blue")


      # Nyquist frequency is the maximum frequency that can be detected for a given sampling rate.

      acq.freq  =  100                    # data acquisition (sample) frequency (Hz)
      time      =  6                      # measuring time interval (seconds)
      ts        =  seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s)
      f.0  =  1/time

      dc.component  =  1
      component.freqs  =  c(3,7,10)        # frequency of signal components (Hz)
      component.delay  =  c(0,0,0)         # delay of signal components (radians)
      component.strength  =  c(1.5,.5,.75) # strength of signal components

      f    =  function(t,w) {
        dc.component +
        sum( component.strength * sin(component.freqs*w*t + component.delay))
      }

      plot.fourier(f,f.0,ts=ts)

      w  =  2*pi*f.0
      trajectory  =  sapply(ts, function(t) f(t,w))
      head(trajectory, n=30)

      X.k  =  fft(trajectory)                   # find all harmonics with fft()
      plot.frequency.spectrum(X.k, xlimits=c(0,20))

      x.n  =  get.trajectory(X.k,ts,acq.freq) / acq.freq  # TODO: why the scaling?
      plot(ts,x.n, type="l"); abline(h=0,lty=3)
      points(ts,trajectory,col="red",type="l") # compare with original


    nyq = floor(length(trajectory) / 4 )
    fY = fft( c(trajectory, rep(0, length(trajectory)) ) )
    ff = fY * Conj(fY)  # power spectrum == S(i) == |F[.]|^2
    ii = fft( ff, inverse=TRUE)
    ac = Re(ii)/(length(ii)  )
    ac = ac/ ac[1]

    x11(); plot(ac[1:nyq]); abline(h=0)

    x11(); plot(acf(tt, nyq))

    vg = 1-ac[1:nyq]
    vx = time*(c(0:(nyq-1)))
    x11(); plot(vg ~ vx)

    fit = stmv_variogram_optimization( vx=vx, vg=vg, nu=0.5 )


    image( Re(fY[1:nr, 1:nc] ))
    xlimits=c(0, length(fY))
    F  =  cbind(0:(length(fY)-1), Mod(fY))
    F[2:length(fY),2]  =  2*F[2:length(fY),2]

    plot(F, t="h", lwd=2, main="",
          xlab="Frequency (Hz)", ylab="Strength",
          xlim=xlimits, ylim=c(0,max(Mod(F[,2]))))


    image(pp)

  }


}



