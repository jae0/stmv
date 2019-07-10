


stmv_matern= function( distance=NULL, mRange, mSmooth=0.5, parameterization="stmv", returntype="autocorrelation" ) {
  # matern autocorrelation function
  # mSmooth =  Bessl smoothness parameter (aks, "nu")
  # mRange = range parameter (aka "phi", .. not the spatial range)

  #\\ NOTE:: the default parameterization is Wikipedia's paramterization:
  #\\ == Rasmussen, Carl Edward (2006) Gaussian Processes for Machine Learning
  #\\  sigma^2 * (2^{1-nu} / Gamma(nu) ) * (sqrt(2*nu) * ||x|| / phi)^{nu} * K_{nu}( sqrt(2*nu) * ||x|| / phi)
  #\\   where K_{nu} is the Bessel function with smooth nu and phi is known as the range parameter
  #\\   ... sigma and tau are not added here ..
  # -------------------------

  if ( mSmooth <= 0 ) warning("mSmooth should be positive to be meaningful")

  mRange_stmv = matern_phi2phi(mRange, mSmooth, parameterization_input=parameterization, parameterization_output="stmv")
  r = distance / mRange_stmv  # distance scaled by range parameter
  r[r<1e-9] = 1e-9

  if (returntype=="covariance") {
    if (mSmooth==0.5) {
      # simplifies to exponential:
      covariance = exp( -r )
    } else if (mSmooth==1.5) {
      u = sqrt(3)*r
      covariance = (1+u) * exp(-u)
    } else if (mSmooth==2.5) {
      u = sqrt(5)*r
      covariance = (1+u+{u^2}/3) * exp(-u)
    } else {
      u = sqrt(2*mSmooth)*r
      covariance = 2^(1-mSmooth) / (gamma(mSmooth)) * u^mSmooth * besselK(u, mSmooth)
    }
    zerodist = which( distance == 0 )
    if (length(zerodist) > 0 ) covariance[zerodist] = 1
    return(covariance)
  }

  if (returntype=="autocorrelation") {
    if (mSmooth==0.5) {
      # simplifies to exponential:
      covariance = exp( -r )
      ac = 1-covariance
    } else if (mSmooth==1.5) {
      u = sqrt(3)*r
      covariance = (1+u) * exp(-u)
      ac =  1-covariance
    } else if (mSmooth==2.5) {
      u = sqrt(5)*r
      covariance = (1+u+(u^2)/3) * exp(-u)
      ac = 1-covariance
    } else {
      u = sqrt(2*mSmooth)*r
      # covariance = 2^(1-mSmooth) / (gamma(mSmooth)) * u^mSmooth * besselK(u, mSmooth)
      ac = 2^(1-mSmooth)/gamma(mSmooth) * (u)^mSmooth *  besselK(u, mSmooth)
    }
    zerodist = which( distance == 0 )
    if (length(zerodist) > 0 ) ac[zerodist] = 0
    # ac = zapsmall(ac)
    return(ac)
  }


}

