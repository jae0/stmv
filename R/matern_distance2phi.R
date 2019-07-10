

matern_distance2phi = function( distance, nu=0.5, cor=0.05, nX=1000, eps=1e-6 ) {

  #\\ NOTE:: the default parameterization is Wikipedia's paramterization:
  #\\ == Rasmussen, Carl Edward (2006) Gaussian Processes for Machine Learning
  #\\  sigma^2 * (2^{1-nu} / Gamma(nu) ) * (sqrt(2*nu) * ||x|| / phi)^{nu} * K_{nu}( sqrt(2*nu) * ||x|| / phi)
  #\\   where K_{nu} is the Bessel function with smooth nu and phi is known as the range parameter
  # -------------------------

    #\\ estimate phi from a range (matern) correlations drops to a given threshold cor
    nu  = max( eps, nu, na.rm=TRUE )
    phi_max = max( eps, 2*distance/sqrt(2*nu), na.rm=TRUE )
    u = matrix(0, ncol=2, nrow=nX )
    u[,2] = seq(0, phi_max, length.out=nX )  # distances
    u[,1] = stmv_matern( distance, mRange=u[,2], mSmooth=nu ) # autocorrel
    phi = approx( u, xout=cor, ties=mean )$y
    return(phi)

}


