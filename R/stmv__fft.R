
stmv__fft = function( p=NULL, dat=NULL, pa=NULL, nu=NULL, phi=NULL, varSpatial=NULL, varObs=NULL, eps=1e-9, ... ) {

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
    eps = 1e-9
  }


  params = list(...)

  vns = p$stmv_variables$LOCS

  pa = data.table(pa)

  sdTotal = sd(dat[[p$stmv_variables$Y]], na.rm=T)

  # system size
  #nr = nx # nr .. x/plon
  #nc = ny # nc .. y/plat

  # data area dimensionality
  dx = p$pres
  dy = p$pres

  if ( grepl("fast_predictions", p$stmv_fft_filter)) {
    # predict only where required
    x_r = range(pa[[vns[1] ]])
    x_c = range(pa[[vns[2] ]])
  } else {
    # predict on full data subset
    x_r = range(dat[[vns[1] ]] )
    x_c = range(dat[[vns[2] ]] )
  }

  rr = diff(x_r)
  rc = diff(x_c)

  nr = trunc( rr/dx ) + 1L
  nc = trunc( rc/dy ) + 1L

  dr = rr/(nr-1) # == dx
  dc = rc/(nc-1) # == dy


  nr2 = 2 * nr
  nc2 = 2 * nc

  # approx sa associated with each datum
  # sa = rr * rc
  # d_sa = sa/nrow(dat) # sa associated with each datum
  # d_length = sqrt( d_sa/pi )  # sa = pi*l^2  # characteristic length scale
  # theta.Taper_static = d_length * 5
  theta.Taper_static = matern_phi2distance( phi=phi, nu=nu, cor=p$stmv_autocorrelation_fft_taper )  # default



  dat$mean = NA
  pa$mean = NA
  pa$sd = sqrt(varSpatial)

  # constainer for spatial filters ( final output grid )
  grid.list = list((1:nr2) * dx, (1:nc2) * dy)
  dgrid = expand_grid_fast(  grid.list[[1]],  grid.list[[2]] )
  dimnames(dgrid) = list(NULL, names(grid.list))
  attr(dgrid, "grid.list") = grid.list
  grid.list = NULL

  center = matrix(c((dx * nr), (dy * nc)), nrow = 1, ncol = 2)

  mC = matrix(0, nrow = nr2, ncol = nc2)
  mC[nr, nc] = 1
  mC_fft = 1 / fftwtools::fftw2d(mC)  # multiplication is faster than division
  mC = NULL

  origin = c(x_r[1], x_c[1])
  res = c(dx, dy)

  # precompute a few things that are used repeatedly for time-specific variograms
  xy = expand_grid_fast( c(-(nr-1):0, 0:(nr-1)) * dr,  c(-(nc-1):0, 0:(nc-1)) * dc )
  distances = sqrt(xy[,1]^2 + xy[,2]^2)
  dmax = max(distances, na.rm=TRUE ) * 0.4  # approx nyquist distance (<0.5 as corners exist)
  breaks = seq( 0, dmax, length.out=nr)
  db = breaks[2] - breaks[1]
  # angles = atan2( xy[,2], xy[,1])  # not used
  acs = data.table( dst=as.numeric(as.character(cut( distances, breaks=c(breaks, max(breaks)+db), label=breaks+(db/2) )))  )

  xy = NULL
  distances = NULL
  breaks = NULL

  coo = data.table(array_map( "xy->2", coords=dat[, ..vns], origin=origin, res=res ))
  names(coo) = c("x", "y")

  good = which(coo[,1] >= 1 & coo[,1] <= nr & coo[,2] >= 1  & coo[,2] )
  if (length(good) > 0) {
    coo = coo[good,]
    dat = dat[good,]
  }
  # default is no time, with time, updated in loop
  xi   = 1:nrow(dat) # all data as p$nt==1
  pa_i = 1:nrow(pa)

  stmv_stats_by_time = NA
  if ( grepl("stmv_variogram_resolve_time", p$stmv_fft_filter)) {
      statsvars_scale_names = c( "sdTotal", "sdSpatial", "sdObs", "phi", "nu", "localrange", "ndata" )
      stmv_stats_by_time = matrix( NA, ncol=length( statsvars_scale_names ), nrow=p$nt )  # container to hold time-specific results
  }

  # things to do outside of the time loop
  if ( grepl("lowpass", p$stmv_fft_filter)) {
    sp.covar.lowpass = stationary.cov( dgrid, center, Covariance="Matern", theta=p$stmv_lowpass_phi, smoothness=p$stmv_lowpass_nu )
    sp.covar.lowpass = as.surface(dgrid, c(sp.covar.lowpass))$z / (nr2*nc2)
  }

  for ( ti in 1:p$nt ) {

    if ( exists("TIME", p$stmv_variables) ) {
      xi   = which( dat[[ p$stmv_variables$TIME ]] == p$prediction_ts[ti] )
      pa_i = which(  pa[[ p$stmv_variables$TIME ]] == p$prediction_ts[ti] )
      if (length(xi) < 5 ) {
        # print( ti)
        next()
      }
    }

    theta.Taper = NULL

    # bounds check: make sure predictions exist
    zmean = mean( dat[xi] [[p$stmv_variables$Y]], na.rm=TRUE)
    zsd = sd( dat[xi] [[p$stmv_variables$Y]], na.rm=TRUE)
    zvar = zsd^2

    X = ( dat[xi] [[p$stmv_variables$Y]] - zmean) / zsd # zscore -- making it mean 0 removes the "DC component"


    if (0) {
      u = as.image(
        X,
        ind=as.matrix(array_map( "xy->2", coords=dat[xi, ..vns], origin=origin, res=res )),
        na.rm=TRUE,
        nx=nr,
        ny=nc
      )
      surface (u)
    }

    # Nadaraya/Watson normalization for missing values
    # See explanation:  https://en.wikipedia.org/wiki/Autocorrelation#Efficient_computation
    # Robertson, C., & George, S. C. (2012). Theory and practical recommendations for autocorrelation-based image
    # correlation spectroscopy. Journal of biomedical optics, 17(8), 080801. doi:10.1117/1.JBO.17.8.080801
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3414238/

    GG = coo[xi,]
    GG$X = c(X)
    yy = GG[, mean(X, na.rm=TRUE), by=.(x, y) ]
    nn = GG[, .N, by=.(x, y) ]
    GG = NULL
    if ( nrow(yy) < 5 ) return( out )

    fY = matrix(0, nrow = nr2, ncol = nc2)
    fN = matrix(0, nrow = nr2, ncol = nc2)
    fY[ cbind(yy[[1]], yy[[2]]) ] = yy[[3]]
    fN[ cbind(nn[[1]], nn[[2]]) ] = nn[[3]]
    yy = nn = NULL

    fY[!is.finite(fY)] = 0
    fN[!is.finite(fN)] = 0

    fY = fftwtools::fftw2d(fY)
    fN = fftwtools::fftw2d(fN)

    # default to crude mean phi and nu
    local_phi = phi
    local_nu = nu

    # update to time-slice averaged values
    if ( grepl("stmv_variogram_resolve_time", p$stmv_fft_filter)) {
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
        acs$X = c(X)
        vgm = acs[, mean(X, na.rm=TRUE), by=dst ][order(dst)]
        names(vgm) = c("distances", "ac")
        vgm= vgm[is.finite(distances) & is.finite(ac)]

        X = NULL

        vgm$sv =  zvar * (1-vgm$ac^2) # each sv are truly orthogonal

      # plot(ac ~ distances, data=vgm   )
      # plot(sv ~ distances, data=vgm   )

        uu = which( (vgm$distances < dmax ) & is.finite(vgm$sv) )  # dmax ~ Nyquist freq

        fit = try( stmv_variogram_optimization( vx=vgm$distances[uu], vg=vgm$sv[uu], plotvgm=FALSE,
          cor=p$stmv_autocorrelation_localrange ))
        uu = NULL

        stmv_stats_by_time[ti, match("sdTotal", statsvars_scale_names ) ] = zsd
        stmv_stats_by_time[ti, match("sdSpatial", statsvars_scale_names ) ] = sqrt(fit$summary$varSpatial)
        stmv_stats_by_time[ti, match("sdObs", statsvars_scale_names ) ] = sqrt(fit$summary$varObs)
        stmv_stats_by_time[ti, match("phi", statsvars_scale_names ) ] = fit$summary$phi
        stmv_stats_by_time[ti, match("nu", statsvars_scale_names ) ] = fit$summary$nu
        stmv_stats_by_time[ti, match("localrange", statsvars_scale_names ) ] = fit$summary$localrange
        stmv_stats_by_time[ti, match("ndata", statsvars_scale_names ) ] = length(xi)

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

        if (grepl("tapered", p$stmv_fft_filter)) {
          # figure out which taper distance to use
          if (grepl("empirical", p$stmv_fft_filter) ) {
            if (!is.null(vgm)) {
              theta.Taper = vgm$distances[ find_intersection( vgm$ac, threshold=p$stmv_autocorrelation_fft_taper ) ]
            }
          } else if (grepl("modelled", p$stmv_fft_filter) ) {
            theta.Taper = matern_phi2distance( phi=local_phi, nu=local_nu, cor=p$stmv_autocorrelation_fft_taper )
          }
        }
        vgm = NULL
  }

    gc()


    sp.covar = 1  # init collator of the kernal

    if ( grepl("lowpass", p$stmv_fft_filter) ) {
      sp.covar = sp.covar * fftwtools::fftw2d(sp.covar.lowpass) * mC_fft  # this is created outside of the loop for speed
    }

    if ( grepl("matern", p$stmv_fft_filter) ) {
      if ( grepl("tapered", p$stmv_fft_filter) ) {
        if (is.null(theta.Taper)) theta.Taper = theta.Taper_static  # if null then not time varying so use static
        spk =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=local_phi, smoothness=local_nu,
          Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2, aRange=1), spam.format=TRUE)
        spk = as.surface(dgrid, c(spk))$z / (nr2*nc2)
        sp.covar = sp.covar * fftwtools::fftw2d(spk)  * mC_fft
      } else {
        spk = stationary.cov( dgrid, center, Covariance="Matern", theta=local_phi, smoothness=local_nu )
        spk = as.surface(dgrid, c(spk))$z / (nr2*nc2)
        sp.covar = sp.covar * fftwtools::fftw2d(spk) * mC_fft
      }
    }

    if ( grepl("normal_kernel", p$stmv_fft_filter) ) {
      theta = matern_phi2distance( phi=local_phi, nu=local_nu, cor=p$stmv_autocorrelation_localrange )
      if ( grepl("tapered", p$stmv_fft_filter) ) {
        if (is.null(theta.Taper)) theta.Taper = theta.Taper_static  # if null then not time varying so use static
        spk =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Exponential", theta=theta,
          Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2, aRange=1 ), spam.format=TRUE)
        spk = as.surface(dgrid, c(spk))$z / (nr2*nc2)
        sp.covar = sp.covar * fftwtools::fftw2d(spk)  * mC_fft
      } else {
        # written out in long form to show the math
        xseq = seq(-(nr - 1), nr, 1) * dx / theta
        yseq = seq(-(nc - 1), nc, 1) * dy / theta
        dd = ((matrix(xseq, nr2, nc2)^2 + matrix(yseq, nr2, nc2, byrow = TRUE)^2))  # squared distances
        yseq = xseq = NULL
        # double.exp: An R function that takes as its argument the _squared_
        # distance between two points divided by the bandwidth. The
        # default is exp( -abs(x)) yielding a normal kernel
        spk = matrix( double.exp(dd), nrow = nr2, ncol = nc2) / (nr2*nc2)
        sp.covar = sp.covar * fftwtools::fftw2d(spk) * mC_fft # kernal weights
      }
    }

    spk = NULL

    fY = Re( fftwtools::fftw2d( sp.covar * fY, inverse = TRUE))[1:nr, 1:nc]
    fN = Re( fftwtools::fftw2d( sp.covar * fN, inverse = TRUE))[1:nr, 1:nc]
    sp.covar = NULL

    X = ifelse((fN > eps), (fY/fN), NA)
    X[!is.finite(X)] = NA
    X = X * zsd + zmean # revert to input scale

    if (0) {
      dev.new()
      image(list(x=c(1:nr)*dr, y=c(1:nc)*dc, z=X), xaxs="r", yaxs="r")
    }

    if ( grepl("fast_predictions", p$stmv_fft_filter)) {
      pa$mean[pa_i] = X
      X = NULL
    } else {
      X_i = array_map( "xy->2", coords=pa[pa_i, ..vns], origin=origin, res=res )
      tokeep = which( X_i[,1] >= 1 & X_i[,2] >= 1  & X_i[,1] <= nr & X_i[,2] <= nc )
      if (length(tokeep) < 1) next()
      X_i = X_i[tokeep,]
      pa$mean[pa_i[tokeep]] = X[X_i]
      X = X_i = NULL
    }


    if ( grepl("stmv_variogram_resolve_time", p$stmv_fft_filter)) {
      pa$sd[pa_i] = stmv_stats_by_time[ti, match("sdSpatial", statsvars_scale_names )]
    }

    iYP = match(
      array_map( "xy->1", dat[ xi, ..vns ], gridparams=p$gridparams ),
      array_map( "xy->1", pa[ pa_i , ..vns ], gridparams=p$gridparams )
    )
    dat$mean[xi] = pa$mean[pa_i][iYP]
  }

  ss = lm( dat$mean ~ dat[[p$stmv_variables$Y]], na.action=na.omit)
  rsquared = ifelse( inherits(ss, "try-error"), NA,  summary(ss)$r.squared )

  if (exists("stmv_rsquared_threshold", p)) {
    if (! is.finite(rsquared) ) return(NULL)
    if (rsquared < p$stmv_rsquared_threshold ) return(NULL)
  }

  stmv_stats = list( sdTotal=sdTotal, rsquared=rsquared, ndata=nrow(dat) )

  if ( grepl("stmv_variogram_resolve_time", p$stmv_fft_filter)) {
    stmv_stats = as.list(colMeans(stmv_stats_by_time, na.rm=TRUE))
    stmv_stats[ match("rsquared", statsvars_scale_names )] = rsquared
  }

  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=TRUE) , aspect="iso" )

  return( list( predictions=pa, stmv_stats=stmv_stats, stmv_stats_by_time=stmv_stats_by_time ) )

}
