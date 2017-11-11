
matern_phi2distance = function( phi, nu=0.5, cor=0.95, dmax=phi*13, nX=5000, eps=1e-3 ) {   

  #\\ estimate distance at which spatial (matern) correlations drops to a given threshold cor
  #\\ NOTE:: the default parameterization is Wikipedia's paramterization:
  #\\ == Rasmussen, Carl Edward (2006) Gaussian Processes for Machine Learning
  #\\  sigma^2 * (2^{1-nu} / Gamma(nu) ) * (sqrt(2*nu) * ||x|| / phi)^{nu} * K_{nu}( sqrt(2*nu) * ||x|| / phi)
  #\\   where K_{nu} is the Bessel function with smooth nu and phi is known as the range parameter  
  # -------------------------

    phi = max( eps, phi, na.rm=TRUE )
    nu  = max( eps, nu, na.rm=TRUE )
    dmax = max( eps, dmax, na.rm=TRUE )
    u = matrix(0, ncol=2, nrow=nX )
    u[,2] = seq(0, dmax, length.out=nX )
    u[,1] = 1 - stm_matern( u[,2], mRange=phi, mSmooth=nu ) # covariance
    distance = approx( u, xout=cor )$y
    return(distance)
}

