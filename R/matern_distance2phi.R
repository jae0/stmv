

matern_distance2phi = function( distance, nu=0.5, cor=0.05, nX=2000, eps=1e-6, plotdata=FALSE ) {

  #\\ NOTE:: the default parameterization is Wikipedia's paramterization:
  #\\ == Rasmussen, Carl Edward (2006) Gaussian Processes for Machine Learning
  #\\  sigma^2 * (2^{1-nu} / Gamma(nu) ) * (sqrt(2*nu) * ||x|| / phi)^{nu} * K_{nu}( sqrt(2*nu) * ||x|| / phi)
  #\\   where K_{nu} is the Bessel function with smooth nu and phi is known as the range parameter
  # -------------------------

  #\\ estimate phi from a range (matern) correlations drops to a given threshold cor
  nu  = max( eps, nu, na.rm=TRUE )
  phi_max = max( eps, distance*100000, na.rm=TRUE )
  u = matrix(0, ncol=2, nrow=2*nX )
  u[,2] = c(seq(0, 1, length.out=nX ), exp(seq( 1+eps, log(phi_max), length.out=nX )))   # distances
  u[,1] = stmv_matern( distance, mRange=u[,2], mSmooth=nu ) # autocorrel
  v = which(duplicated(u[,1]))
  if (length(v) > 0) u=u[-v,]
  phi = approx( u, xout=cor, ties=mean )$y

  if(plotdata) {
    plot( u[,1] ~ u[,2], xlim = c(0, phi*1.25))
    abline( v=phi )
  }

  return(phi)

  if (0) {
    loadfunctions("stmv")
    dis = 100
    nu = 1
    cor = 0.05

    dis = 0.2
    nu=0.11
    cor=0.1

    matern_distance2phi( dis, nu, cor, plotdata=TRUE)


  }
}


