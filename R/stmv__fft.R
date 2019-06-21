
stmv__fft = function( p=NULL, dat=NULL, pa=NULL, nu=NULL, phi=NULL, variablelist=FALSE, tol=1e-9, ... ) {

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

  if (variablelist)  return( c() )

  params = list(...)

  sdTotal = sd(dat[,p$variable$Y], na.rm=T)

  # nr .. x/plon
  # nc .. y/plat
  # pa_r = range(pa[,p$variables$LOCS[1]])
  # pa_c = range(pa[,p$variables$LOCS[2]])

  x_r = range(dat[,p$variables$LOCS[1]])
  x_c = range(dat[,p$variables$LOCS[2]])

  nr = floor( diff(x_r)/p$pres ) + 1
  nc = floor( diff(x_c)/p$pres ) + 1

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

  if (!exists("stmv_fft_filter",p) ) p$stmv_fft_filter="lowpass" # default in case of no specification

  if ( p$stmv_fft_filter == "lowpass") {
    theta = p$stmv_lowpass_phi
    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=theta, smoothness=p$stmv_lowpass_nu )
    sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
    sp.covar.kernel = fftwtools::fftw2d(sp.covar) / fftwtools::fftw2d(mC)
  }

  if (p$stmv_fft_filter %in% c("matern") ) {
    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=phi, smoothness=nu )
    sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
    sp.covar.kernel = fftwtools::fftw2d(sp.covar) / fftwtools::fftw2d(mC)
  }

  if (p$stmv_fft_filter == "lowpass_matern") {
    # both ..
    sp.covar.lowpass = stationary.cov( dgrid, center, Covariance="Matern", theta=p$stmv_lowpass_phi, smoothness=p$stmv_lowpass_nu )
    sp.covar.lowpass = as.surface(dgrid, c(sp.covar.lowpass))$z / (nr2*nc2)
    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=phi, smoothness=nu )
    sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
    sp.covar.kernel = ( fftwtools::fftw2d(sp.covar.lowpass) / fftwtools::fftw2d(mC) ) * (fftwtools::fftw2d(sp.covar)/ fftwtools::fftw2d(mC) )
    sp.covar = sp.covar.lowpass = NULL
  }

  if (p$stmv_fft_filter == "matern_tapered") {
    theta.Taper = matern_phi2distance( phi=phi, nu=nu, cor=p$stmv_fft_taper_factor )
    sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=phi, smoothness=nu,
      Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)
    sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
    sp.covar.kernel = fftwtools::fftw2d(sp.covar) / fftwtools::fftw2d(mC)
    sp.covar = theta.Taper = NULL
  }

  if (p$stmv_fft_filter == "lowpass_matern_tapered") {
    sp.covar.lowpass = stationary.cov( dgrid, center, Covariance="Matern", theta=p$stmv_lowpass_phi, smoothness=p$stmv_lowpass_nu )
    sp.covar.lowpass = as.surface(dgrid, c(sp.covar.lowpass))$z / (nr2*nc2)
    theta.Taper = matern_phi2distance( phi=phi, nu=nu, cor=p$stmv_fft_taper_factor )
    sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=phi, smoothness=nu,
      Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)
    sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
    sp.covar.kernel = (fftwtools::fftw2d(sp.covar.lowpass) / fftwtools::fftw2d(mC) ) * ( fftwtools::fftw2d(sp.covar)/ fftwtools::fftw2d(mC) )
    sp.covar = sp.covar.lowpass = theta.Taper = NULL
  }


  if (p$stmv_fft_filter == "normal_kernel") {
      theta = matern_phi2distance( phi=phi, nu=nu, cor=p$stmv_range_correlation )
      xi = seq(-(nr - 1), nr, 1) * dx / theta
      yi = seq(-(nc - 1), nc, 1) * dy / theta
      dd = ((matrix(xi, nr2, nc2)^2 + matrix(yi, nr2, nc2, byrow = TRUE)^2))  # squared distances
      # double.exp: An R function that takes as its argument the _squared_
      # distance between two points divided by the bandwidth. The
      # default is exp( -abs(x)) yielding a normal kernel
      kk = double.exp(dd)
      mK = matrix(kk, nrow = nr2, ncol = nc2)
      sp.covar.kernel = fftwtools::fftw2d(mK) / fftwtools::fftw2d(mC) / (nr2*nc2)# kernal weights
  }

  sp.covar = sp.covar.lowpass = dgrid = center =  NULL

  origin=c(x_r[1], x_c[1])
  res=c(p$pres, p$pres)

  zz = matrix(1:(nr*nc), nrow = nr, ncol = nc)
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

    u = as.image(
      dat[xi, p$variables$Y],
      ind=as.matrix(array_map( "xy->2", coords=dat[xi, p$variables$LOCS], origin=origin, res=res )),
      na.rm=TRUE,
      nx=nr,
      ny=nc
    )

    mY[1:nr,1:nc] = u$z
    mY[!is.finite(mY)] = 0

    #  Nadaraya/Watson normalization for missing values s
    mN[1:nr,1:nc] = u$weights
    mN[!is.finite(mN)] = 0

    u =NULL

    fY = Re(fftwtools::fftw2d( sp.covar.kernel * fftwtools::fftw2d(mY), inverse = TRUE))[1:nr, 1:nc]  #real amplitudes
    fN = Re(fftwtools::fftw2d( sp.covar.kernel * fftwtools::fftw2d(mN), inverse = TRUE))[1:nr, 1:nc]
    Z = ifelse((fN > tol), (fY/fN), NA)
    fY = fN = NULL

    # low pass filter based upon a global nu,phi .. remove high freq variation
    Z_i = array_map( "xy->2", coords=pa[pa_i,p$variables$LOCS], origin=origin, res=res )

    # bounds check: make sure predictions exist
    Z_i_test = which( Z_i[,1]<1 | Z_i[,2]<1  | Z_i[,1] > nr | Z_i[,2] > nc )

    if (length(Z_i_test) > 0) {
      keep = zz[ Z_i[-Z_i_test,] ]
      if ( length(keep) > 0) pa$mean[pa_i[keep]] = Z[keep]
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


  # -------------------------------



  if (0) {

    # testing and debugging

    loadfunctions( c("aegis", "stmv"))
    RLibrary(c ("fields", "MBA", "geoR") )

    # xyz = stmv_test_data( datasource="swiss" )
    # xy = xyz[, c("x", "y")]
    # mz = log( xyz$rain )
    # mm = lm( mz ~ 1 )
    # z = residuals( mm)
    # xyz = cbind(xyz[, c("x", "y")], z)
    # gr = stmv_variogram( xy, z, methods="geoR", plotdata=TRUE ) # ml via profile likelihood


    xyz = stmv_test_data( datasource="meuse" )
    xy = xyz[, c("x", "y")]
    z = log(xyz$elev)

    xyz = cbind(xyz[, c("x", "y")], z=((z-mean(z))/sd(z) ) )
    gr = stmv_variogram( xy, z, methods="geoR", plotdata=TRUE ) # ml via profile likelihood

    nu = gr$geoR$nu
    phi = gr$geoR$phi

    fit  =  Krig(xyz[, c("x", "y")], xyz[,"z"], theta=phi)
    x11()
    surface( fit, type="C") # look at the surface

    x11()
    require(MBA)
    mba.int  =  mba.surf( xyz, 100, 100, extend=TRUE)$xyz.est
    surface(mba.int, xaxs="r", yaxs="r")


    x_r = range(xyz$x)
    x_c = range(xyz$y)

    rez = diff(x_r)/100
    nr = floor( diff(x_r)/rez ) + 1
    nc = floor( diff(x_c)/rez ) + 1

    dx = dy = rez

    nr2 = 2 * nr
    nc2 = 2 * nc

    mC = matrix(0, nrow = nr2, ncol = nc2)
    mC[nr, nc] = 1

    # constainer for spatial filters
    grid.list = list((1:nr2) * dx, (1:nc2) * dy)
    dgrid = as.matrix(expand.grid(grid.list))
    dimnames(dgrid) = list(NULL, names(grid.list))
    attr(dgrid, "grid.list") = grid.list

    center = matrix(c((dx * nr), (dy * nc)), nrow = 1, ncol = 2)

    theta.Taper = matern_phi2distance( phi=phi, nu=nu, cor=p$stmv_fft_taper_factor )
    sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=phi, smoothness=nu,
      Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)
    sp.covar = as.surface(dgrid, c(sp.covar))$z
    sp.covar.kernel = fft(sp.covar) / fft(mC)

    mY = matrix(0, nrow = nr2, ncol = nc2)
    mN = matrix(0, nrow = nr2, ncol = nc2)

    u = as.image(
      Z=xyz$z,
      x=xyz[, c("x", "y")],
      na.rm=TRUE,
      nx=nr,
      ny=nc
    )
    # surface(u)

    mY[1:nr,1:nc] = u$z
    mY[!is.finite(mY)] = 0

    #  Nadaraya/Watson normalization for missing values s
    mN[1:nr,1:nc] = u$weights
    mN[!is.finite(mN)] = 0

    fY = Re( fft( sp.covar.kernel * fft(mY) / (nr2 * nc2), inverse = TRUE) )[1:nr, 1:nc]
    fN = Re( fft( sp.covar.kernel * fft(mN) / (nr2 * nc2), inverse = TRUE) )[1:nr, 1:nc]

    tol = 1e-9
    Z = fY/fN
    Z = ifelse((fN > tol), (fY/fN), NA)
    Z[!is.finite(Z)] = NA

    x11()

    surface(list(x=c(1:nr)*dx, y=c(1:nc)*dy, z=Z), xaxs="r", yaxs="r")


    # 2D
    # Robertson, C., & George, S. C. (2012). Theory and practical recommendations for autocorrelation-based image
    # correlation spectroscopy. Journal of biomedical optics, 17(8), 080801. doi:10.1117/1.JBO.17.8.080801
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3414238/


    fY = fft(mY) / (nr2*nc2)
    YY = fY * Conj(fY) # power spectra
    ii = Re( fft( YY, inverse=TRUE)  )  # autocorrelation (amplitude)

    fN = fft(mN) / (nr2*nc2)
    NN = fN * Conj(fN) # power spectra
    jj = Re( fft( NN, inverse = TRUE)  ) # autocorrelation (amplitude) correction


    tol = 1e-9
    kk = ifelse(( jj > tol), (ii / jj), NA) # autocorrelation
    # kk[kk< -1] = -1
    # kk[kk>1] = 1


    oo = fftshift(kk)

    xy = expand.grid(
      x = c(-nr:-1, 1:nr) * dx,
      y = c(-nc:-1, 1:nc) * dy
    )

    distances = sqrt(xy$x^2 + xy$y^2)
    dmax = max(distances, na.rm=TRUE ) * 0.4  # approx nyquist distance (0.9 as corners exist)
    breaks = seq( 0, dmax, length.out=20)
    db = breaks[2] - breaks[1]
    # angles = atan2( xy$y, xy$x )  # not used

    ii = cut( distances, breaks=c(breaks, max(breaks)+db), label=breaks+(db/2) )

    res = as.data.frame.table(tapply( X=oo, INDEX=ii, FUN=mean, na.rm=TRUE ))
    names(res) = c("distances", "ac")
    res$distances = as.numeric( as.character(res$distances))

    res$sv =  var(xyz$z) * (1-res$ac^2)
    uu = which(res$distances < dmax/3) & is.finite(res$sv)
    plot(ac ~ distances, data=res[uu,]   )
    plot(sv ~ distances, data=res[uu,]   )

    fit = try( stmv_variogram_optimization( vx=res$distances[uu], vg=res$sv[uu], plotvgm=TRUE,
      stmv_internal_scale=dmax/4, cor=p$stmv_fft_taper_factor ))

  }


  # -------------------------


  if (0) {

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
        x.n  =  rep(0,N)           # create vector to keep the trajectory
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



