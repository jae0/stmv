
matern_phi2phi = function( mRange, mSmooth, parameterization_input="stmv", parameterization_output="stmv" ) {
  # convert range parameter to/from stmv form:

  #\\ NOTE:: the default parameterization is Wikipedia's paramterization:
  #\\ == Rasmussen, Carl Edward (2006) Gaussian Processes for Machine Learning
  #\\  sigma^2 * (2^{1-nu} / Gamma(nu) ) * (sqrt(2*nu) * ||x|| / phi)^{nu} * K_{nu}( sqrt(2*nu) * ||x|| / phi)
  #\\   where K_{nu} is the Bessel function with smooth nu and phi is known as the range parameter
  # -------------------------

  need.more.testing = c("geostatsp", "CompRandFld")
  of = grepl(parameterization_input, need.more.testing)
  ot = grepl(parameterization_output, need.more.testing)
  oo = of | ot
  if ( any(oo) ) {
    io = need.more.testing[ which( oo ) ]
    warning( "Need to check parameterization of ", io)
  }

  mRange_stmv = switch(parameterization_input,
      stmv = mRange ,
      wikipedia = mRange,
      RandomFields = mRange ,
      # CompRandFld = sqrt(2*mSmooth )/mRange,  # R(h) = 2^{1-nu} Gamma(nu)^{-1} x^nu  * K_nu(h) ; ||x||==h;  -- same as wikipedia except x=sqrt( 2nu*d/phi)
      geoR = mRange * sqrt(2*mSmooth ),
      fields = mRange * sqrt(2*mSmooth ),
      gstat = mRange * sqrt(2*mSmooth ),
      spBayes =  sqrt(2*mSmooth )/mRange ,
      bayesx =  sqrt(2*mSmooth )/mRange   ,
      inla = sqrt(2*mSmooth)/mRange,  # kappa_inla = sqrt(2*nu)/phi
      geostatsp = mRange*2,   #   sqrt(8*mSmooth) = 2 * sqrt(2*mSmooth)
      "not found"
  )

  mRange_output = switch( parameterization_output,
      stmv = mRange_stmv,
      wikipedia = mRange_stmv,
      RandomFields = mRange_stmv ,
      # CompRandFld = sqrt(2*mSmooth ) / mRange_stmv ,
      geoR = mRange_stmv * sqrt(2*mSmooth ),
      fields = mRange_stmv * sqrt(2*mSmooth ),
      gstat = mRange_stmv * sqrt(2*mSmooth ),
      spBayes = sqrt(2*mSmooth) /mRange_stmv  ,
      bayesx = sqrt(2*mSmooth ) / mRange_stmv ,
      inla = sqrt(2*mSmooth) / mRange_stmv  ,
      geostatsp = 2*mRange_stmv,  #   sqrt(8*mSmooth) = 2 * sqrt(2*mSmooth)
     "not found"
  )

  if ("not found" %in% c(mRange_stmv, mRange_output) ) warning("Parameterization not found:", parameterization_input, parameterization_output)

  return(mRange_output)
}
