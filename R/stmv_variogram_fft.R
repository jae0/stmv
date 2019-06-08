
stmv_variogram_fft = function( xyz, nx=101, ny=101, nbreaks=20, plotvgm=FALSE, eps=1e-9, add.interpolation=FALSE, cor=0.1 ) {

  names(xyz) =c("x", "y", "z")
  zmean = mean(xyz$z, na.rm=TRUE)
  zsd = sd(xyz$z, na.rm=TRUE)
  zvar = zsd^2
  Z = (xyz$z - zmean) / zsd # zscore

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

  fY = fft(mY)/(nr2*nc2)
  fN = fft(mN)/(nr2*nc2)

  # fY * Conj(fY) == power spectra
  ii = Re( fft( fY * Conj(fY), inverse=TRUE)  ) # autocorrelation (amplitude)
  jj = Re( fft( fN * Conj(fN), inverse=TRUE)  ) # autocorrelation (amplitude) correction

  X = ifelse(( jj > eps), (ii / jj), NA) # autocorrelation
  ii = jj = NULL

  # fftshift
  X = rbind(X[((nr+1):nr2), (1:nc2)], X[(1:nr), (1:nc2)])  # swap_up_down
  X = cbind(X[1:nr2, ((nc+1):nc2)], X[1:nr2, 1:nc])  # swap_left_right

  # radial representation
  xy = expand.grid( x = c(-nr:-1, 1:nr) * dr,  y = c(-nc:-1, 1:nc) * dc )
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
  # plot(ac ~ distances, data=res[uu,]   )
  # plot(sv ~ distances, data=res[uu,]   )

  out = list(res=res )

  if (add.interpolation) {
    # interpolated surface
  # constainer for spatial filters
    uu = which( (res$distances < dmax/3) & is.finite(res$sv) )
    fit = try( stmv_variogram_optimization( vx=res$distances[uu], vg=res$sv[uu], plotvgm=plotvgm, stmv_internal_scale=dmax/4, cor=cor ))
    phi = fit$summary$phi
    nu = fit$summary$nu

    grid.list = list((1:nr2) * dr, (1:nc2) * dc)
    dgrid = as.matrix(expand.grid(grid.list))
    dimnames(dgrid) = list(NULL, names(grid.list))
    attr(dgrid, "grid.list") = grid.list
    grid.list = NULL
    mC = matrix(0, nrow = nr2, ncol = nc2)
    mC[nr, nc] = 1
    mC =NULL
    center = matrix(c((dr * nr), (dc * nc)), nrow = 1, ncol = 2)
    theta.Taper = matern_phi2distance( phi=phi, nu=nu, cor=cor )
    sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=phi, smoothness=nu,
      Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)
    sp.covar = as.surface(dgrid, c(sp.covar))$z
    sp.covar.kernel = fft(sp.covar) / fft(mC)
    ffY = Re(fft( sp.covar.kernel * fY, inverse = TRUE))[1:nr, 1:nc]
    ffN = Re(fft( sp.covar.kernel * fN, inverse = TRUE))[1:nr, 1:nc]
    Z = ifelse((ffN > eps), (ffY/ffN), NA)
    Z[!is.finite(Z)] = NA
    Z = Z * zsd + zmean # revert to input scale
    if (plotdata) {
      dev.new()
      surface(list(x=c(1:nr)*dr, y=c(1:nc)*dc, z=Z), xaxs="r", yaxs="r")
    }
    out = c(out, fit=fit, Z=Z)
  }

  return(out)
}
