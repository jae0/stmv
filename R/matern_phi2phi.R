
matern_phi2phi = function( mRange, mSmooth, parameterization_input="ecmei", parameterization_output="ecmei" ) {
  # convert range parameter to/from ecmei form:

  #\\ NOTE:: the default parameterization is Wikipedia's paramterization:
  #\\ == Rasmussen, Carl Edward (2006) Gaussian Processes for Machine Learning
  #\\  sigma^2 * (2^{1-nu} / Gamma(nu) ) * (sqrt(2*nu) * ||x|| / phi)^{nu} * K_{nu}( sqrt(2*nu) * ||x|| / phi)
  #\\   where K_{nu} is the Bessel function with smooth nu and phi is known as the range parameter  
  # -------------------------

  need.more.testing = c("geostatsp", "spate") 
  of = grepl(parameterization_input, need.more.testing)
  ot = grepl(parameterization_output, need.more.testing)
  oo = of | ot
  if ( any(oo) ) {
    io = need.more.testing[ which( oo ) ]
    warning( "Need to check parameterization of ", io)
  }
  
  mRange_ecmei = switch(parameterization_input,
      ecmei = mRange , 
      wikipedia = mRange,  
      RandomFields = mRange ,
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
      ecmei = mRange_ecmei,
      wikipedia = mRange_ecmei,
      RandomFields = mRange_ecmei ,
      geoR = mRange_ecmei * sqrt(2*mSmooth ),
      fields = mRange_ecmei * sqrt(2*mSmooth ), 
      gstat = mRange_ecmei * sqrt(2*mSmooth ), 
      spBayes = sqrt(2*mSmooth) /mRange_ecmei  ,
      bayesx = sqrt(2*mSmooth ) / mRange_ecmei ,
      inla = sqrt(2*mSmooth) / mRange_ecmei  ,   
      geostatsp = 2*mRange_ecmei,  #   sqrt(8*mSmooth) = 2 * sqrt(2*mSmooth) 
     "not found"
  )

  if ("not found" %in% c(mRange_ecmei, mRange_output) ) warning("Parameterization not found:", parameterization_input, parameterization_output)

  return(mRange_output)
}
