
matern_phi2phi = function( mRange, mSmooth, parameterization_input="emei", parameterization_output="emei" ) {
  # convert range parameter to/from emei form:

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
  
  mRange_emei = switch(parameterization_input,
      emei = mRange , 
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
      emei = mRange_emei,
      wikipedia = mRange_emei,
      RandomFields = mRange_emei ,
      geoR = mRange_emei * sqrt(2*mSmooth ),
      fields = mRange_emei * sqrt(2*mSmooth ), 
      gstat = mRange_emei * sqrt(2*mSmooth ), 
      spBayes = sqrt(2*mSmooth) /mRange_emei  ,
      bayesx = sqrt(2*mSmooth ) / mRange_emei ,
      inla = sqrt(2*mSmooth) / mRange_emei  ,   
      geostatsp = 2*mRange_emei,  #   sqrt(8*mSmooth) = 2 * sqrt(2*mSmooth) 
     "not found"
  )

  if ("not found" %in% c(mRange_emei, mRange_output) ) warning("Parameterization not found:", parameterization_input, parameterization_output)

  return(mRange_output)
}
