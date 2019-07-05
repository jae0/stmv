
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
  # sa = rr * rc
  # d_sa = sa/nrow(dat) # sa associated with each datum
  # d_length = sqrt( d_sa/pi )  # sa = pi*l^2  # characteristic length scale
  # theta.Taper = d_length * 5


  # final output grid
  # x_locs = expand_grid_fast(
  #   seq( x_r[1], x_r[2], length.out=nr ),
  #   seq( x_c[1], x_c[2], length.out=nc )
  # )
  # attr( x_locs , "out.attrs") = NULL
  # names( x_locs ) = p$variables$LOCS

  dat$mean = NA
  pa$mean = NA
  pa$sd = sdTotal  # this is ignored with fft

  # constainer for spatial filters
  grid.list = list((1:nr2) * dx, (1:nc2) * dy)
  # dgrid = as.matrix(expand.grid(grid.list))  # a bit slower
  dgrid = expand_grid_fast(  grid.list[[1]],  grid.list[[2]] )
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

  # precompute a few things that are used repeatedly for time-specific variograms
  if ( exists("TIME", p$variables) ) {
    xy = expand_grid_fast( c(-(nr-1):0, 0:(nr-1)) * dr,  c(-(nc-1):0, 0:(nc-1)) * dc )
    distances = sqrt(xy[,1]^2 + xy[,2]^2)
    dmax = max(distances, na.rm=TRUE ) * 0.4  # approx nyquist distance (<0.5 as corners exist)
    breaks = seq( 0, dmax, length.out=nr)
    db = breaks[2] - breaks[1]
    # angles = atan2( xy[,2], xy[,1])  # not used
    zz = cut( distances, breaks=c(breaks, max(breaks)+db), label=breaks+(db/2) )
    xy = NULL
    distances = NULL
    breaks = NULL
  }

  fY = matrix(0, nrow = nr2, ncol = nc2)
  fN = matrix(0, nrow = nr2, ncol = nc2)

  # default is no time, with time, updated in loop
  xi   = 1:nrow(dat) # all data as p$nt==1
  pa_i = 1:nrow(pa)

  for ( ti in 1:p$nt ) {

    if ( exists("TIME", p$variables) ) {
      xi   = which( dat[ , p$variables$TIME] == p$prediction.ts[ti] )
      pa_i = which(  pa[ , p$variables$TIME] == p$prediction.ts[ti] )
      if (length(xi) < 5 ) {
        # print( ti)
        next()
      }
    }

    # bounds check: make sure predictions exist
    z = dat[xi, p$variables$Y]
    zmean = mean(z, na.rm=TRUE)
    zsd = sd(z, na.rm=TRUE)
    zvar = zsd^2
    z = (z - zmean) / zsd # zscore -- making it mean 0 removes the DC component

    if (0) {
      u = as.image(
        z,
        ind=as.matrix(array_map( "xy->2", coords=dat[xi, p$variables$LOCS], origin=origin, res=resolution )),
        na.rm=TRUE,
        nx=nr,
        ny=nc
      )
      # surface (u)
      # fY[1:nr,1:nc] = u$z
      # fN[1:nr,1:nc] = u$weights
      # u =NULL
  }

    # Nadaraya/Watson normalization for missing values
    # See explanation:  https://en.wikipedia.org/wiki/Autocorrelation#Efficient_computation
    # Robertson, C., & George, S. C. (2012). Theory and practical recommendations for autocorrelation-based image
    # correlation spectroscopy. Journal of biomedical optics, 17(8), 080801. doi:10.1117/1.JBO.17.8.080801
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3414238/

    coo = as.matrix(array_map( "xy->2", coords=dat[xi, p$variables$LOCS], origin=origin, res=resolution ))
    yy = tapply( X=z, INDEX=list(coo[,1], coo[,2]), FUN = function(w) {mean(w, na.rm=TRUE)}, simplify=TRUE )
    nn = tapply( X=z, INDEX=list(coo[,1], coo[,2]), FUN = function(w) {length(w)}, simplify=TRUE )
    fY[as.numeric(dimnames(yy)[[1]]),as.numeric(dimnames(yy)[[2]])] = yy
    fN[as.numeric(dimnames(nn)[[1]]),as.numeric(dimnames(nn)[[2]])] = nn
    yy = nn = NULL
    coo = NULL


    fY[!is.finite(fY)] = 0
    fN[!is.finite(fN)] = 0

    fY = fftwtools::fftw2d(fY)
    fN = fftwtools::fftw2d(fN)

    local_phi = phi
    local_nu = nu

    # create time-specific variograms
    # fY * Conj(fY) == power spectra
    acY = Re( fftwtools::fftw2d( fY * Conj(fY), inverse=TRUE)  ) # autocorrelation (amplitude)
    acN = Re( fftwtools::fftw2d( fN * Conj(fN), inverse=TRUE)  ) # autocorrelation (amplitude) correction
    X = ifelse(( acN > eps), (acY / acN), NA) # autocorrelation (amplitude)
    acY = acN = NULL
    # fftshift
    X = rbind( X[((nr+1):nr2), (1:nc2)], X[(1:nr), (1:nc2)] )  # swap_up_down
    X = cbind( X[1:nr2, ((nc+1):nc2)], X[1:nr2, 1:nc])  # swap_left_right
    # zz (distance breaks) precomputed outside of loop
    vgm = as.data.frame.table(tapply( X=X, INDEX=zz, FUN=mean, na.rm=TRUE ))
    names(vgm) = c("distances", "ac")
    X = NULL

    theta.Taper = vgm$distances[ find_intersection( vgm$ac, threshold=p$stmv_fft_taper_correlation ) ]
    theta.Taper = theta.Taper * p$stmv_fft_taper_fraction # fraction of the distance to 0 correlation; sqrt(0.5) = ~ 70% of the variability (associated with correlation = 0.5)

    vgm$distances = as.numeric( as.character(vgm$distances))
    vgm$sv =  zvar * (1-vgm$ac^2) # each sv are truly orthogonal

  # plot(ac ~ distances, data=vgm   )
  # plot(sv ~ distances, data=vgm   )

    uu = which( (vgm$distances < dmax ) & is.finite(vgm$sv) )  # dmax ~ Nyquist freq
    fit = try( stmv_variogram_optimization( vx=vgm$distances[uu], vg=vgm$sv[uu], plotvgm=FALSE,
      stmv_internal_scale=dmax*0.75, cor=p$stmv_range_correlation ))
    uu = NULL
    vgm = NULL
    # out$fit = fit

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

    fit = NULL


    gc()

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

      sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=local_phi, smoothness=local_nu,
        Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)
      sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
      sp.covar = fftwtools::fftw2d(sp.covar) / mC_fft
  }

    if (p$stmv_fft_filter == "lowpass_matern_tapered") {
      sp.covar.lowpass = stationary.cov( dgrid, center, Covariance="Matern", theta=p$stmv_lowpass_phi, smoothness=p$stmv_lowpass_nu )
      sp.covar.lowpass = as.surface(dgrid, c(sp.covar.lowpass))$z / (nr2*nc2)

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
      yseq = xseq = NULL
      # double.exp: An R function that takes as its argument the _squared_
      # distance between two points divided by the bandwidth. The
      # default is exp( -abs(x)) yielding a normal kernel
      sp.covar = matrix( double.exp(dd), nrow = nr2, ncol = nc2)
      sp.covar = fftwtools::fftw2d(sp.covar) / mC_fft / (nr2*nc2)# kernal weights
    }

    fY = Re( fftwtools::fftw2d( sp.covar * fY, inverse = TRUE))[1:nr, 1:nc]
    fN = Re( fftwtools::fftw2d( sp.covar * fN, inverse = TRUE))[1:nr, 1:nc]
    sp.covar = NULL
    X = ifelse((fN > eps), (fY/fN), NA)
    X[!is.finite(X)] = NA
    X = X * zsd + zmean # revert to input scale
    if (0) {
      dev.new()
      surface(list(x=c(1:nr)*dr, y=c(1:nc)*dc, z=Z), xaxs="r", yaxs="r")
    }

    X_i = array_map( "xy->2", coords=pa[pa_i, p$variables$LOCS], origin=origin, res=resolution )
    tokeep = which( X_i[,1] >= 1 & X_i[,2] >= 1  & X_i[,1] <= nr & X_i[,2] <= nc )
    if (length(tokeep) < 1) next()
    X_i = X_i[tokeep,]

    pa$mean[pa_i[tokeep]] = X[X_i]
      # pa$sd[pa_i] = NA  ## fix as NA
    X = X_i = NULL
  }

  # dat_i = array_map( "xy->2", coords=dat[, p$variables$LOCS], origin=origin, res=resolution )
  # ss = lm( pa$mean[dat_i] ~ dat[,p$variables$Y], na.action=na.omit)
  # if ( "try-error" %in% class( ss ) ) return( NULL )
  # rsquared = summary(ss)$r.squared
  # if (rsquared < p$stmv_rsquared_threshold ) return(NULL)

  stmv_stats = list( sdTotal=sdTotal, rsquared=NA, ndata=nrow(dat) ) # must be same order as p$statsvars

  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=TRUE) , aspect="iso" )

  return( list( predictions=pa, stmv_stats=stmv_stats ) )


}



