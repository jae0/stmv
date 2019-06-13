
stmv_variogram_fft = function( xyz, nx=64, ny=64, nbreaks=13, plotdata=FALSE, eps=1e-8, add.interpolation=FALSE,
  stmv_range_correlation=0.1, stmv_range_correlation_fft_taper=0.05 ) {

  names(xyz) =c("x", "y", "z")
  zmean = mean(xyz$z, na.rm=TRUE)
  zsd = sd(xyz$z, na.rm=TRUE)
  zvar = zsd^2
  Z = (xyz$z - zmean) / zsd # zscore -- making it mean 0 removes the DC component

  nr = nx
  nc = ny

  x_r = range(xyz$x)
  x_c = range(xyz$y)

  dr = diff(x_r)/(nr-1)
  dc = diff(x_c)/(nc-1)

  nr2 = 2 * nr
  nc2 = 2 * nc

  mY = matrix(0, nrow = nr2, ncol = nc2)
  mN = matrix(0, nrow = nr2, ncol = nc2)

  u = as.image( Z=Z, x=xyz[, c("x", "y")], na.rm=TRUE, nx=nr, ny=nc )
  # surface(u)

  mY[1:nr,1:nc] = u$z
  mY[!is.finite(mY)] = 0

  #  Nadaraya/Watson normalization for missing values s
  mN[1:nr,1:nc] = u$weights
  mN[!is.finite(mN)] = 0

  # See explanation:  https://en.wikipedia.org/wiki/Autocorrelation#Efficient_computation
  # Robertson, C., & George, S. C. (2012). Theory and practical recommendations for autocorrelation-based image
  # correlation spectroscopy. Journal of biomedical optics, 17(8), 080801. doi:10.1117/1.JBO.17.8.080801
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3414238/

  fY = fftwtools::fftw2d(mY)
  fN = fftwtools::fftw2d(mN)

  # fY * Conj(fY) == power spectra
  ii = Re( fftwtools::fftw2d( fY * Conj(fY), inverse=TRUE)  ) # autocorrelation (amplitude)
  jj = Re( fftwtools::fftw2d( fN * Conj(fN), inverse=TRUE)  ) # autocorrelation (amplitude) correction

  X = ifelse(( jj > eps), (ii / jj), NA) # autocorrelation (amplitude)
  ii = jj = NULL

  # fftshift
  X = rbind( X[((nr+1):nr2), (1:nc2)], X[(1:nr), (1:nc2)] )  # swap_up_down
  X = cbind(X[1:nr2, ((nc+1):nc2)], X[1:nr2, 1:nc])  # swap_left_right

  # X[c(nr: (nr+1)), nc:(nc+1)]

  # radial representation
  xy = expand.grid( x = c(-(nr-1):0, 0:(nr-1)) * dr,  y = c(-(nc-1):0, 0:(nc-1)) * dc )
  distances = sqrt(xy$x^2 + xy$y^2)
  dmax = max(distances, na.rm=TRUE ) * 0.4  # approx nyquist distance (<0.5 as corners exist)
  breaks = seq( 0, dmax, length.out=nbreaks)
  db = breaks[2] - breaks[1]
  # angles = atan2( xy$y, xy$x )  # not used

  zz = cut( distances, breaks=c(breaks, max(breaks)+db), label=breaks+(db/2) )
  distances = NULL
  xy = NULL

  res = as.data.frame.table(tapply( X=X, INDEX=zz, FUN=mean, na.rm=TRUE ))
  names(res) = c("distances", "ac")
  res$distances = as.numeric( as.character(res$distances))

  res$sv =  zvar * (1-res$ac^2) # each sv are truly orthogonal

  # plot(ac ~ distances, data=res[which(res$distances < median(res$distances)),]   )
  # plot(sv ~ distances, data=res[which(res$distances < median(res$distances)),]   )

  out = list(res=res )

  if (add.interpolation) {
    # interpolated surface
  # constainer for spatial filters
    uu = which( (res$distances < dmax*0.75 ) & is.finite(res$sv) )  # dmax ~ Nyquist freq
    fit = try( stmv_variogram_optimization( vx=res$distances[uu], vg=res$sv[uu], plotvgm=plotdata,
      stmv_internal_scale=dmax*0.5, cor=stmv_range_correlation ))
    out$fit = fit

    if (any(!is.finite( c(fit$summary$phi, fit$summary$nu) ))) return(out)

    phi = fit$summary$phi
    nu = fit$summary$nu

    grid.list = list((1:nr2) * dr, (1:nc2) * dc)
    dgrid = as.matrix(expand.grid(grid.list))
    dimnames(dgrid) = list(NULL, names(grid.list))
    attr(dgrid, "grid.list") = grid.list
    grid.list = NULL
    mC = matrix(0, nrow = nr2, ncol = nc2)
    mC[nr, nc] = 1
    center = matrix(c((dr * nr), (dc * nc)), nrow = 1, ncol = 2)
    theta.Taper = matern_phi2distance( phi=phi, nu=nu, cor=stmv_range_correlation_fft_taper )
    sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=phi, smoothness=nu,
      Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)

    sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
    sp.covar.kernel = fftwtools::fftw2d(sp.covar) / fftwtools::fftw2d(mC)
    ffY = Re( fftwtools::fftw2d( sp.covar.kernel * fY, inverse = TRUE))[1:nr, 1:nc]
    ffN = Re( fftwtools::fftw2d( sp.covar.kernel * fN, inverse = TRUE))[1:nr, 1:nc]
    Z = ifelse((ffN > eps), (ffY/ffN), NA)
    Z[!is.finite(Z)] = NA
    Z = Z * zsd + zmean # revert to input scale
    if (plotdata) {
      dev.new()
      surface(list(x=c(1:nr)*dr, y=c(1:nc)*dc, z=Z), xaxs="r", yaxs="r")
    }
    out$Z = Z
  }

  return(out)



  if (0) {
    # testing and debugging

    loadfunctions( c("aegis", "stmv"))
    RLibrary(c ("fields", "MBA", "geoR") )

    if (0) {
      XYZ = stmv_test_data( datasource="swiss" )
      mz = log( XYZ$rain )
      mm = lm( mz ~ 1 )
      XYZ$z = residuals( mm)
      XYZ=XYZ[c("x","y","z")]
    }

    if (0) {
      XYZ = stmv_test_data( datasource="meuse" )
      XYZ$z = log(XYZ$elev)
      XYZ=XYZ[, c("x","y","z") ]
    }

    gr = stmv_variogram( XYZ[, c("x", "y")], XYZ[,"z"], methods="geoR", plotdata=TRUE ) # ml via profile likelihood

    nu = gr$geoR$nu
    phi = gr$geoR$phi

    fit  =  Krig(XYZ[, c("x", "y")], XYZ[,"z"], theta=phi)
    x11(); surface( fit, type="C") # look at the surface

    mba.int  =  mba.surf( XYZ, 64, 64, extend=TRUE)$xyz.est
    x11(); surface(mba.int, xaxs="r", yaxs="r")

    oo = stmv_variogram_fft( XYZ[c("x","y","z")], nx=64, ny=64, nbreaks=32, plotdata=TRUE,  add.interpolation=TRUE, stmv_range_correlation_fft_taper=0.01 )

  }

}
