
stmv_variogram_fft = function( xyz, nx=NULL, ny=NULL, nbreaks=30, plotdata=FALSE, eps=1e-9, add.interpolation=FALSE,
  stmv_fft_taper_method="modelled", stmv_autocorrelation_localrange=0.1, stmv_autocorrelation_fft_taper=0, stmv_fft_taper_fraction=sqrt(0.5) ) {

  if (0) {
    require(fields)
    loadfunctions( c("aegis", "stmv"))
    RLibrary(c ("fields", "MBA", "geoR") )

    nx = ny = 128

    XYZ = as.data.frame( cbind( RMprecip$x, RMprecip$y ) )
    names(XYZ) =c( "x","y","z")

    xyz = XYZ[, c("x", "y", "z")]
    nbreaks=30
    plotdata=FALSE
    eps=1e-9
    add.interpolation=TRUE
    stmv_autocorrelation_localrange=0.1
    stmv_autocorrelation_fft_taper=0

  }

  out = list()

  names(xyz) =c("x", "y", "z")
  zmean = mean(xyz$z, na.rm=TRUE)
  zsd = sd(xyz$z, na.rm=TRUE)
  zvar = zsd^2
  z = (xyz$z - zmean) / zsd # zscore -- making it mean 0 removes the DC component

  if (is.null(nx)) {
    nx = ny = floor( nbreaks * 2.35 )
  }

  # system size
  nr = nx
  nc = ny

  x_r = range(xyz$x)
  x_c = range(xyz$y)

  # system length scale
  rr = diff(x_r)
  rc = diff(x_c)

  # no of elements
  dr = rr/(nr-1)
  dc = rc/(nc-1)

  # # approx sa associated with each datum
  # sa = rr * rc
  # d_sa = sa/nrow(xyz) # sa associated with each datum
  # d_length = sqrt( d_sa/pi )  # sa = pi*radius^2  # characteristic length scale

  nr2 = 2 * nr
  nc2 = 2 * nc

  fY = matrix(0, nrow = nr2, ncol = nc2)
  fN = matrix(0, nrow = nr2, ncol = nc2)

  if (0) {
    u = as.image( Z=z, x=xyz[, c("x", "y")], na.rm=TRUE, nx=nr, ny=nc )
    # surface(u)
    u = NULL
  }

  out$nr = nr
  out$nc = nc
  out$origin = origin = c(x_r[1], x_c[1])
  out$resolution  = resolution = c(dr, dc)

  #  Nadaraya/Watson normalization for missing values s
  coo = as.matrix(array_map( "xy->2", coords=xyz[,c("x", "y")], origin=origin, res=resolution ))
  yy = tapply( X=z, INDEX=list(coo[,1], coo[,2]),  FUN = function(w) {mean(w, na.rm=TRUE)}, simplify=TRUE )
  nn = tapply( X=z, INDEX=list(coo[,1], coo[,2]), FUN = function(w) {length(w)}, simplify=TRUE )

  if (length (dimnames(yy)[[1]]) < 5 ) return( out )

  fY[as.numeric(dimnames(yy)[[1]]),as.numeric(dimnames(yy)[[2]])] = yy
  fN[as.numeric(dimnames(nn)[[1]]),as.numeric(dimnames(nn)[[2]])] = nn
  yy = nn = NULL
  coo = NULL

  #  Nadaraya/Watson normalization for missing valuess
  fY[!is.finite(fY)] = 0
  fN[!is.finite(fN)] = 0


  # See explanation:  https://en.wikipedia.org/wiki/Autocorrelation#Efficient_computation
  # Robertson, C., & George, S. C. (2012). Theory and practical recommendations for autocorrelation-based image
  # correlation spectroscopy. Journal of biomedical optics, 17(8), 080801. doi:10.1117/1.JBO.17.8.080801
  # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3414238/

  fY = fftwtools::fftw2d(fY)
  fN = fftwtools::fftw2d(fN)

  # fY * Conj(fY) == power spectra
  acY = Re( fftwtools::fftw2d( fY * Conj(fY), inverse=TRUE)  ) # autocorrelation (amplitude)
  acN = Re( fftwtools::fftw2d( fN * Conj(fN), inverse=TRUE)  ) # autocorrelation (amplitude) correction
  X = ifelse(( acN > eps), (acY / acN), NA) # autocorrelation (amplitude)
  acN = NULL
  acY = NULL

  # fftshift
  X = rbind( X[((nr+1):nr2), (1:nc2)], X[(1:nr), (1:nc2)] )  # swap_up_down
  X = cbind(X[1:nr2, ((nc+1):nc2)], X[1:nr2, 1:nc])  # swap_left_right


  # radial representation
  xy = expand_grid_fast( c(-(nr-1):0, 0:(nr-1)) * dr,  c(-(nc-1):0, 0:(nc-1)) * dc )
  distances = sqrt(xy[,1]^2 + xy[,2]^2)
  dmax = max(distances, na.rm=TRUE ) * 0.4  # approx nyquist distance (<0.5 as corners exist)
  breaks = seq( 0, dmax, length.out=nr)
  db = breaks[2] - breaks[1]
  # angles = atan2( xy[,2], xy[,1])  # not used


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

  out$vgm = vgm

  if (add.interpolation) {
    # interpolated surface
  # constainer for spatial filters
    uu = which( (vgm$distances < dmax ) & is.finite(vgm$sv) )  # dmax ~ Nyquist freq

    if (plotdata) dev.new()
    out$fit = try( stmv_variogram_optimization( vx=vgm$distances[uu], vg=vgm$sv[uu], plotvgm=plotdata,
      stmv_internal_scale=dmax*0.75, cor=stmv_autocorrelation_localrange ) )

    if (any(!is.finite( c(out$fit$summary$phi, out$fit$summary$nu) ))) return(out)

    phi = out$fit$summary$phi
    nu = out$fit$summary$nu

    grid.list = list((1:nr2) * dr, (1:nc2) * dc)
    # dgrid = as.matrix(expand.grid(grid.list))  # a bit slower
    dgrid = expand_grid_fast(  grid.list[[1]],  grid.list[[2]] )
    dimnames(dgrid) = list(NULL, names(grid.list))
    attr(dgrid, "grid.list") = grid.list
    grid.list = NULL

    mC = matrix(0, nrow = nr2, ncol = nc2)
    mC[nr, nc] = 1  # equal weights
    center = matrix(c((dr * nr), (dc * nc)), nrow = 1, ncol = 2)


    if (stmv_fft_taper_method == "empirical") {
      theta.Taper = vgm$distances[ find_intersection( vgm$ac, threshold=stmv_autocorrelation_fft_taper ) ]
      theta.Taper = theta.Taper * stmv_fft_taper_fraction # fraction of the distance to 0 correlation; sqrt(0.5) = ~ 70% of the variability (associated with correlation = 0.5)
    } else if (stmv_fft_taper_method == "modelled") {
      theta.Taper = matern_phi2distance( phi=phi, nu=nu, cor=stmv_autocorrelation_fft_taper )
    }

    sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=phi, smoothness=nu,
      Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)
    sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
    sp.covar = fftwtools::fftw2d(sp.covar)  / fftwtools::fftw2d(mC)
    mC = NULL
    fY = Re( fftwtools::fftw2d( sp.covar * fY, inverse = TRUE))[1:nr, 1:nc]
    fN = Re( fftwtools::fftw2d( sp.covar * fN, inverse = TRUE))[1:nr, 1:nc]
    X = ifelse((fN > eps), (fY/fN), NA)
    fN = NULL
    fY = NULL
    X[!is.finite(X)] = NA
    X = X * zsd + zmean # revert to input scale
    if (plotdata) {
      dev.new()
      plot(sv ~ distances, data=out$vgm   )

      dev.new()
      surface(list(x=c(1:nr)*dr, y=c(1:nc)*dc, z=X), xaxs="r", yaxs="r")

    }
    out$Z = X
  }

  return(out)



  if (0) {
    # testing and debugging

    loadfunctions( c("aegis", "stmv"))
    RLibrary(c ("fields", "MBA", "geoR") )

    nx = ny = 128
    nx = ny = 64

    if (0) {
      XYZ = stmv_test_data( datasource="swiss" )
      mz = log( XYZ$rain )
      mm = lm( mz ~ 1 )
      XYZ$z = residuals( mm)
      XYZ=XYZ[c("x","y","z")]
      x11()
      oo = stmv_variogram_fft( XYZ[c("x","y","z")], nx=nx, ny=ny, nbreaks=nx, plotdata=TRUE, add.interpolation=TRUE,
        stmv_fft_taper_method="modelled", stmv_autocorrelation_fft_taper=0.05 )

      oo = stmv_variogram_fft( XYZ[c("x","y","z")], nx=nx, ny=ny, nbreaks=nx, plotdata=TRUE, add.interpolation=TRUE,
        stmv_fft_taper_method="empirical", stmv_autocorrelation_fft_taper=0.1, stmv_fft_taper_fraction=0.75 )

      gr = stmv_variogram( XYZ[, c("x", "y")], XYZ[,"z"], methods="fft", plotdata=TRUE ) # fft/nl least squares
    }


    if (0) {
      XYZ = stmv_test_data( datasource="meuse" )
      XYZ$z = log(XYZ$elev)
      XYZ=XYZ[, c("x","y","z") ]
      x11()
      oo = stmv_variogram_fft( XYZ[c("x","y","z")], nx=nx, ny=ny, nbreaks=nx, plotdata=TRUE,  add.interpolation=TRUE,
        stmv_fft_taper_method="modelled", stmv_autocorrelation_fft_taper=0.25  )

      oo = stmv_variogram_fft( XYZ[c("x","y","z")], nx=nx, ny=ny, nbreaks=nx, plotdata=TRUE,  add.interpolation=TRUE,
        stmv_fft_taper_method="empirical",  stmv_autocorrelation_fft_taper=0.1, stmv_fft_taper_fraction=0.5 )

      gr = stmv_variogram( XYZ[, c("x", "y")], XYZ[,"z"], methods="fft", plotdata=TRUE )
    }

    if (0) {
      str(fields::RMprecip)
      fit <- smooth.2d( RMprecip$y, x=RMprecip$x, theta=.5)
      image( fit )
      points( RMprecip$x, pch=".")

      XYZ = as.data.frame( cbind( RMprecip$x, RMprecip$y ) )
      names(XYZ) =c( "x","y","z")
      oo = stmv_variogram_fft( XYZ[c("x","y","z")], nx=nx, ny=ny, nbreaks=nx, plotdata=TRUE,  add.interpolation=TRUE,
          stmv_fft_taper_method="modelled", stmv_autocorrelation_fft_taper=0.2  )

      oo = stmv_variogram_fft( XYZ[c("x","y","z")], nx=nx, ny=ny, nbreaks=nx, plotdata=TRUE,  add.interpolation=TRUE,
          stmv_fft_taper_method="empirical",  stmv_autocorrelation_fft_taper=0.1, stmv_fft_taper_fraction=sqrt(.5) )

      gr = stmv_variogram( XYZ[, c("x", "y")], XYZ[,"z"], methods="fft", plotdata=TRUE )
    }

    variomethod="geoR"
    variomethod="gstat"
    variomethod="fft"
    variomethod="inla"

    gr = stmv_variogram( XYZ[, c("x", "y")], XYZ[,"z"], methods=variomethod, plotdata=TRUE ) # ml via profile likelihood
    nu = gr[[variomethod]]$nu
    phi = gr[[variomethod]]$phi
    fit  =  Krig(XYZ[, c("x", "y")], XYZ[,"z"], theta=phi)
    x11(); surface( fit, type="C") # look at the surface

    mba.int  =  mba.surf( XYZ[, c("x", "y", "z")], nx, ny, extend=TRUE)$xyz.est
    x11(); surface(mba.int, xaxs="r", yaxs="r")

    str(fields::RMprecip)
    fit <- smooth.2d( RMprecip$y, x=RMprecip$x, theta=.25)
    x11()
    image( fit )
    points( RMprecip$x, pch=".")

  }

}
