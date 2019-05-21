
stmv__fft = function( p=NULL, dat=NULL, pa=NULL, nu=NULL, phi=NULL, distance=NULL, variablelist=FALSE, tol=1e-9, ... ) {

  #\\ this is the core engine of stmv .. localised space (no-time) modelling interpolation
  #\\ note: time is not being modelled and treated independently
  #\\      .. you had better have enough data in each time slice
  #\\ first a low-pass filter as defined by p$stmv_lowpass_nu, p$stmv_lowpass_phi,
  #\\ then a simple covariance filter determined by nu,phi ;; fft (no lowpass)
  #\\ based upon fields::image.smooth and setup.image.smooth
  #\\ now that the local area of interest is stationary, we use local convolutions of highly autocorrelated (short-range)
  #\\ determined by p

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
    theta = p$stmv_lowpass_phi
    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=theta, smoothness=p$stmv_lowpass_nu )
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
    sp.covar.lowpass = stationary.cov( dgrid, center, Covariance="Matern", theta=p$stmv_lowpass_phi, smoothness=p$stmv_lowpass_nu )
    sp.covar.lowpass = as.surface(dgrid, c(sp.covar.lowpass))$z

    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=phi, smoothness=nu )
    sp.covar = as.surface(dgrid, c(sp.covar))$z
    sp.covar.kernel = {fft(sp.covar.lowpass) / fft_mC } * {fft(sp.covar)/ fft_mC }
    sp.covar = sp.covar.lowpass = NULL
  }

  if (p$stmv_fft_filter == "matern_tapered") {
    theta.Taper = matern_phi2distance( phi=phi, nu=nu, cor=p$stmv_range_correlation_fft_taper )
    sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=phi, smoothness=nu,
      Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)
    sp.covar = as.surface(dgrid, c(sp.covar))$z
    sp.covar.kernel = fft(sp.covar) / fft_mC
    sp.covar = theta.Taper = NULL
  }

  if (p$stmv_fft_filter == "lowpass_matern_tapered") {
    sp.covar.lowpass = stationary.cov( dgrid, center, Covariance="Matern", theta=p$stmv_lowpass_phi, smoothness=p$stmv_lowpass_nu )
    sp.covar.lowpass = as.surface(dgrid, c(sp.covar.lowpass))$z

    theta.Taper = matern_phi2distance( phi=phi, nu=nu, cor=p$stmv_range_correlation_fft_taper )
    sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=phi, smoothness=nu,
      Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)
    sp.covar = as.surface(dgrid, c(sp.covar))$z
    sp.covar.kernel = {fft(sp.covar.lowpass) / fft_mC } * {fft(sp.covar)/ fft_mC }
    sp.covar = sp.covar.lowpass = theta.Taper = NULL
  }


  if (p$stmv_fft_filter == "normal_kernel") {
      theta = matern_phi2distance( phi=phi, nu=nu, cor=0.5 )
      xi = seq(-(nr - 1), nr, 1) * dx / theta
      yi = seq(-(nc - 1), nc, 1) * dy / theta
      dd = ((matrix(xi, nr2, nc2)^2 + matrix(yi, nr2, nc2, byrow = TRUE)^2))  # squared distances
      # double.exp: An R function that takes as its argument the _squared_
      # distance between two points divided by the bandwidth. The
      # default is exp( -abs(x)) yielding a normal kernel
      kk = double.exp(dd)
      mK = matrix(kk, nrow = nr2, ncol = nc2)
      W = fft(mK) / fft_mC
      sp.covar.kernel = W/(nr2 * nc2)  # kernal weights
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

    fY = Re(fft(fft(mY) * sp.covar.kernel, inverse = TRUE))[1:nr, 1:nc]  #real amplitudes
    fN = Re(fft(fft(mN) * sp.covar.kernel, inverse = TRUE))[1:nr, 1:nc]
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


}

if (0) {

  # testing and debugging

  # loadfunctions( c("aegis.env", "aegis", "stmv"))
  # RLibrary(c ("fields", "MBA", "geoR") )

  # xyz = stmv_test_data( datasource="swiss" )
  # xy = xyz[, c("x", "y")]
  # mz = log( xyz$rain )
  # mm = lm( mz ~ 1 )
  # z = residuals( mm)
  # xyz = cbind(xyz[, c("x", "y")], z)
  # gr = stmv_variogram( xy, z, methods="geoR", plotdata=TRUE ) # ml via profile likelihood


  # xyz = stmv_test_data( datasource="meuse" )
  # xy = xyz[, c("x", "y")]
  # z = log(xyz$elev)
  # xyz = cbind(xyz[, c("x", "y")], z)
  # gr = stmv_variogram( xy, z, methods="geoR", plotdata=TRUE ) # ml via profile likelihood

  # nu = gr$geoR$nu
  # phi = gr$geoR$phi

  # fit <- Krig(xyz[, c("x", "y")], xyz[,"z"], theta=phi)
  # x11()
  # surface( fit, type="C") # look at the surface

  # x11()
  # require(MBA)
  # mba.int <- mba.surf( xyz, 100, 100, extend=TRUE)$xyz.est
  # surface(mba.int, xaxs="r", yaxs="r")

  # x_r = range(xyz$x)
  # x_c = range(xyz$y)

  # rez = diff(x_r)/100
  # nr = round( diff(x_r)/rez ) + 1
  # nc = round( diff(x_c)/rez ) + 1

  # dx = dy = rez

  # nr2 = 2 * nr
  # nc2 = 2 * nc

  # mC = matrix(0, nrow = nr2, ncol = nc2)
  # mC[nr, nc] = 1
  # fft_mC = fft(mC) * nr2 * nc2

  # # constainer for spatial filters
  # grid.list = list((1:nr2) * dx, (1:nc2) * dy)
  # dgrid = as.matrix(expand.grid(grid.list))
  # dimnames(dgrid) = list(NULL, names(grid.list))
  # attr(dgrid, "grid.list") = grid.list

  # center = matrix(c((dx * nr), (dy * nc)), nrow = 1, ncol = 2)

  # theta.Taper = matern_phi2distance( phi=phi, nu=nu, cor=0.5 )
  # sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=phi, smoothness=nu,
  #   Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)
  # sp.covar = as.surface(dgrid, c(sp.covar))$z
  # sp.covar.kernel = fft(sp.covar) / fft_mC

  # mY = matrix(0, nrow = nr2, ncol = nc2)
  # mN = matrix(0, nrow = nr2, ncol = nc2)

  # u = as.image(
  #   Z=xyz$z,
  #   x=xyz[, c("x", "y")],
  #   na.rm=TRUE,
  #   nx=nr,
  #   ny=nc
  # )
  # # surface(u)

  # mY[1:nr,1:nc] = u$z
  # mY[!is.finite(mY)] = 0

  # #  Nadaraya/Watson normalization for missing values s
  # mN[1:nr,1:nc] = u$weights
  # mN[!is.finite(mN)] = 0

  # fY = Re(fft(fft(mY) * sp.covar.kernel, inverse = TRUE))[1:nr, 1:nc]
  # fN = Re(fft(fft(mN) * sp.covar.kernel, inverse = TRUE))[1:nr, 1:nc]

  # tol = 1e-12
  # Z = fY/fN
  # Z = ifelse((fN > tol), (fY/fN), NA)
  # Z[!is.finite(Z)] = NA

  # x11()

  # surface(list(x=c(1:nr)*dx, y=c(1:nc)*dy, z=Z), xaxs="r", yaxs="r")

  # -------------------------

  # useful tools and links

    # https://www.iro.umontreal.ca/~pift6080/H09/documents/
    # https://www.iro.umontreal.ca/~pift6080/H09/documents/eck_fft.pdf
    # https://riptutorial.com/r/example/14464/fourier-series
    # http://www.di.fc.ul.pt/~jpn/r/fourier/fourier.html


    # # cs is the vector of complex points to convert
    # convert.fft <- function(cs, sample.rate=1) {
    # cs <- cs / length(cs) # normalize

    # distance.center <- function(c)signif( Mod(c),        4)
    # angle           <- function(c)signif( 180*Arg(c)/pi, 3)

    # df <- data.frame(cycle    = 0:(length(cs)-1),
    #                   freq     = 0:(length(cs)-1) * sample.rate / length(cs),
    #                   strength = sapply(cs, distance.center),
    #                   delay    = sapply(cs, angle))
    # df
    # }

    # convert.fft(fft(1:4))

    # plot.frequency.spectrum <- function(X.k, xlimits=c(0,length(X.k))) {
    # plot.data  <- cbind(0:(length(X.k)-1), Mod(X.k))

    # # TODO: why this scaling is necessary?
    # plot.data[2:length(X.k),2] <- 2*plot.data[2:length(X.k),2]

    # plot(plot.data, t="h", lwd=2, main="",
    #       xlab="Frequency (Hz)", ylab="Strength",
    #       xlim=xlimits, ylim=c(0,max(Mod(plot.data[,2]))))
    # }

    # # Plot the i-th harmonic
    # # Xk: the frequencies computed by the FFt
    # #  i: which harmonic
    # # ts: the sampling time points
    # # acq.freq: the acquisition rate
    # plot.harmonic <- function(Xk, i, ts, acq.freq, color="red") {
    # Xk.h <- rep(0,length(Xk))
    # Xk.h[i+1] <- Xk[i+1] # i-th harmonic
    # harmonic.trajectory <- get.trajectory(Xk.h, ts, acq.freq=acq.freq)
    # points(ts, harmonic.trajectory, type="l", col=color)
    # }

    # # returns the x.n time series for a given time sequence (ts) and
    # # a vector with the amount of frequencies k in the signal (X.k)
    # get.trajectory <- function(X.k,ts,acq.freq) {

    # N   <- length(ts)
    # i   <- complex(real = 0, imaginary = 1)
    # x.n <- rep(0,N)           # create vector to keep the trajectory
    # ks  <- 0:(length(X.k)-1)

    # for(n in 0:(N-1)) {       # compute each time point x_n based on freqs X.k
    #   x.n[n+1] <- sum(X.k * exp(i*2*pi*ks*n/N)) / N
    # }

    # x.n * acq.freq
    # }

    # X.k <- fft(c(4,0,0,0))                   # get amount of each frequency k

    # time     <- 4                            # measuring time interval (seconds)
    # acq.freq <- 100                          # data acquisition frequency (Hz)
    # ts  <- seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s)

    # x.n <- get.trajectory(X.k,ts,acq.freq)   # create time wave

    # plot(ts,x.n,type="l",ylim=c(-2,4),lwd=2)
    # abline(v=0:time,h=-2:4,lty=3); abline(h=0)

    # plot.harmonic(X.k,1,ts,acq.freq,"red")
    # plot.harmonic(X.k,2,ts,acq.freq,"green")
    # plot.harmonic(X.k,3,ts,acq.freq,"blue")


    # Nyquist frequency is the maximum frequency that can be detected for a given sampling rate.

    # acq.freq <- 100                    # data acquisition (sample) frequency (Hz)
    # time     <- 6                      # measuring time interval (seconds)
    # ts       <- seq(0,time-1/acq.freq,1/acq.freq) # vector of sampling time-points (s)
    # f.0 <- 1/time

    # dc.component <- 1
    # component.freqs <- c(3,7,10)        # frequency of signal components (Hz)
    # component.delay <- c(0,0,0)         # delay of signal components (radians)
    # component.strength <- c(1.5,.5,.75) # strength of signal components

    # f   <- function(t,w) {
    # dc.component +
    # sum( component.strength * sin(component.freqs*w*t + component.delay))
    # }

    # plot.fourier(f,f.0,ts=ts)

    # w <- 2*pi*f.0
    # trajectory <- sapply(ts, function(t) f(t,w))
    # head(trajectory,n=30)

    # X.k <- fft(trajectory)                   # find all harmonics with fft()
    # plot.frequency.spectrum(X.k, xlimits=c(0,20))

    # x.n <- get.trajectory(X.k,ts,acq.freq) / acq.freq  # TODO: why the scaling?
    # plot(ts,x.n, type="l"); abline(h=0,lty=3)
    # points(ts,trajectory,col="red",type="l") # compare with original


}