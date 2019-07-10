
matern_phi2distance = function( phi, nu=0.5, cor=0.05, dmax=phi*100000, nX=2000, eps=1e-6, plotdata=FALSE ) {

  #\\ estimate distance at which spatial (matern) correlations drops to a given threshold cor
  #\\ NOTE:: the default parameterization is Wikipedia's paramterization:
  #\\ == Rasmussen, Carl Edward (2006) Gaussian Processes for Machine Learning
  #\\  sigma^2 * (2^{1-nu} / Gamma(nu) ) * (sqrt(2*nu) * ||x|| / phi)^{nu} * K_{nu}( sqrt(2*nu) * ||x|| / phi)
  #\\   where K_{nu} is the Bessel function with smooth nu and phi is known as the range parameter
  # -------------------------

  phi = max( eps, phi, na.rm=TRUE )
  nu  = max( eps, nu, na.rm=TRUE )
  dmax = max( eps, dmax, na.rm=TRUE )

  u = matrix(0, ncol=2, nrow=nX*2 )
  u[,2] = c(seq(0, 1, length.out=nX ), exp(seq( 1+eps, log(dmax), length.out=nX )) )  # distances
  u[,1] = stmv_matern( u[,2], mRange=phi, mSmooth=nu ) # autocorrelation
  distance = approx( u, xout=cor )$y

  if(plotdata) {
    plot( u[,1] ~ (u[,2]), xlim = c(0, distance*1.25) )
    abline( v=(distance) )
  }


  return(distance)


  if (0) {
    loadfunctions("stmv")
    phi = 10
    nu = 1
    cor = 0.05

    phi = 80
    nu=0.11
    cor=0.9

    matern_phi2distance( phi, nu, cor, plotdata=TRUE)


  }
}

