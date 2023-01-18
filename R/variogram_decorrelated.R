 
 variogram_decorrelated = function( xy, z, 
    nx=NULL, ny=NULL, nbreaks=30, 
    plotdata=FALSE, internal_scale=NULL
  ) {

  # copy and modification of stmv_variogram_fft to be stand-alone 
  
  #\\ estimate empirical variograms (actually correlation functions)
  #\\ and then model them as Matern basis
  #\\ returns empirical variogram and parameter estimates, 
  #\\ expect xy = c(plon, plat), z=variable

  #\\ NOTE:: the default parameterization is Wikipedia's paramterization:
  #\\ == Rasmussen, Carl Edward (2006) Gaussian Processes for Machine Learning
  #\\ == RandomFields and geostatsp parameterization ss
  #\\  sigma^2 * (2^{1-nu} / Gamma(nu) ) * (sqrt(2*nu) * ||x|| / phi)^{nu} * K_{nu}( sqrt(2*nu) * ||x|| / phi )
  #\\   where K_{nu} is the Bessel function with smooth nu and phi is known as the range parameter
  # -------------------------
  
  zvar = var( z, na.rm=TRUE )
  if (!is.finite(zvar)) zvar = 0 
  if (zvar == 0) return(NULL)

  zmean = mean( z, na.rm=TRUE)
  zsd = sqrt( zvar )
  z = ( z - zmean) / zsd # zscore -- making it mean 0 removes the DC component

  out = list(
    Ndata = length(z),
    zmean = zmean,
    zsd = zsd,
    zvar = zvar,  # this is the scaling factor for semivariance .. diving by sd, below reduces numerical floating point issues
    range_crude = sqrt( diff(range(xy[,1]))^2 + diff(range(xy[,2]))^2) / 4  # initial scaling distance beyond which domain edge effects become important
  )
 
 # A crude GUESS AT PHI:
  out$internal_scale  = ifelse( is.null(internal_scale), matern_distance2phi( out$range_crude, nu=0.5 ), 1)  # the presumed scaling distance to make calcs use smaller numbers
  
  xy = xy /  out$internal_scale

 
  if (is.null(nx)) {
    nx = ny = trunc( nbreaks * 2.35 )
  }

  # system size
  nr = nx
  nc = ny

  nr2 = 2 * nr
  nc2 = 2 * nc

  x_r = range(xy$x)
  x_c = range(xy$y)

  # system length scale
  rr = diff(x_r)
  rc = diff(x_c)

  # no of elements
  dr = rr/(nr-1)
  dc = rc/(nc-1)

  fY = matrix(0, nrow = nr2, ncol = nc2)
  fN = matrix(0, nrow = nr2, ncol = nc2)
 
  out$nr = nr
  out$nc = nc
  out$origin = origin = c(x_r[1], x_c[1])
  out$res  = res = c(dr, dc)


  #  Nadaraya/Watson normalization for missing values s
  coo = data.table(
    x = trunc( (xy$x - origin[1]) / dr ) + 1L ,
    y = trunc( (xy$y - origin[2]) / dc ) + 1L ,
    z = z
  ) 
  yy = coo[, mean(z, na.rm=TRUE), by=.(x, y) ]
  nn = coo[, .N, by=.(x, y) ]
  coo = NULL

  if ( nrow(yy) < 5 ) return( out )
  fY[ cbind(yy[[1]], yy[[2]]) ] = yy[[3]]
  fN[ cbind(nn[[1]], nn[[2]]) ] = nn[[3]]
  yy = nn = NULL

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
  X = ifelse(( acN > 0), (acY / acN), NA) # autocorrelation (amplitude)

  acN = NULL
  acY = NULL

  # fftshift
  X = rbind( X[((nr+1):nr2), (1:nc2)], X[(1:nr), (1:nc2)] )  # swap_up_down
  X = cbind( X[1:nr2, ((nc+1):nc2)], X[1:nr2, 1:nc])  # swap_left_right
 
  # radial representation

  # expand_grid_fast 
  # copied from http://stackoverflow.com/questions/10405637/use-outer-instead-of-expand-grid
  seq1 = c(-(nr-1):0, 0:(nr-1)) * dr
  seq2 = c(-(nc-1):0, 0:(nc-1)) * dc
  rc = cbind( 
    Var1 = rep.int(seq1, length(seq2)),
    Var2 = rep.int(seq2, rep.int(length(seq1),length(seq2)))
  )
  seq1 = seq2 = NULL

  distances = sqrt(rc[,1]^2 + rc[,2]^2) #  on internal scale
  # angles = atan2( rc[,2], rc[,1])  # not used
  rc = NULL

  dmax = max( distances, na.rm=TRUE ) * 0.4  # approx nyquist distance (<0.5 as corners exist)
  breaks = seq( 0, dmax, length.out=nr)
  db = breaks[2] - breaks[1]

  acs = data.table( 
    dst=as.numeric(as.character( 
      cut( distances, breaks=c(breaks, max(breaks)+db), label=breaks+(db/2) ) 
    )), 
    X=c(X) 
  )
  
  vgm = acs[, mean(X, na.rm=TRUE), by=dst ][order(dst)]
  names(vgm) = c("distances", "ac")
  vgm= vgm[is.finite(distances) & is.finite(ac)]

  distances = NULL
  X = NULL
  acs =NULL
  gc()

  vgm$distances = as.numeric( as.character(vgm$distances)) * out$internal_scale 
  vgm$sv =  out$zvar * (1-vgm$ac^2) # each sv are truly orthogonal
  
  uu = which( (vgm$distances < (dmax* out$internal_scale) ) & is.finite(vgm$sv) )  # dmax ~ Nyquist freq
  vgm$keep = FALSE
  if (length(uu)>0) vgm$keep[uu] = TRUE

  out$vgm = vgm
  
  if (plotdata) {
    uu = which( out$vgm$keep )
    plot( out$vgm$distances[uu], out$vgm$sv[uu], col="red", pch=19  )
  }

  return(out)


  if (0) {

    loadfunctions( c( "aegis", "stmv" ))

    xyz = stmv_test_data( datasource="meuse" )  
    xyz$log_zinc = log( xyz$zinc )
    xyz$z = residuals( lm( log_zinc ~ 1, xyz ) )
 
    # size_diagonal = sqrt( diff(range(xyz$x))^2 + diff(range(xyz$y))^2) # crude scale of problem (domain size)
    out = variogram_decorrelated( 
      xy=xyz[, c("x", "y")], 
      z=xyz$z, 
      nbreaks=30 
    )

    # first a gstat solution as reference:
    library(sp)
    library(gstat)
    xyz2 = xyz 
    coordinates(xyz2) = ~x+y
    vgm1 <- variogram( z~1, xyz2)
    plot(vgm1)

    toplot = which(out$vgm$keep)
    points( out$vgm$distances[toplot], out$vgm$sv[toplot], col="red", pch=19  )
   
    if (0) {
      nbreaks=50
      plotdata=TRUE
      internal_scale=NULL
      nx = NULL
    }

  }

}


