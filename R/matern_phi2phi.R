
matern_phi2phi = function( mRange, mSmooth, parameterization_input="stm", parameterization_output="stm" ) {
  # convert range parameter to/from stm form:

  #\\ NOTE:: the default parameterization is Wikipedia's paramterization:
  #\\ == Rasmussen, Carl Edward (2006) Gaussian Processes for Machine Learning
  #\\  sigma^2 * (2^{1-nu} / Gamma(nu) ) * (sqrt(2*nu) * ||x|| / phi)^{nu} * K_{nu}( sqrt(2*nu) * ||x|| / phi)
  #\\   where K_{nu} is the Bessel function with smooth nu and phi is known as the range parameter  
  # -------------------------

  need.more.testing = c("geostatsp") 
  of = grepl(parameterization_input, need.more.testing)
  ot = grepl(parameterization_output, need.more.testing)
  oo = of | ot
  if ( any(oo) ) {
    io = need.more.testing[ which( oo ) ]
    warning( "Need to check parameterization of ", io)
  }
  
  mRange_stm = switch(parameterization_input,
      stm = mRange , 
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
      stm = mRange_stm,
      wikipedia = mRange_stm,
      RandomFields = mRange_stm ,
      geoR = mRange_stm * sqrt(2*mSmooth ),
      fields = mRange_stm * sqrt(2*mSmooth ), 
      gstat = mRange_stm * sqrt(2*mSmooth ), 
      spBayes = sqrt(2*mSmooth) /mRange_stm  ,
      bayesx = sqrt(2*mSmooth ) / mRange_stm ,
      inla = sqrt(2*mSmooth) / mRange_stm  ,   
      geostatsp = 2*mRange_stm,  #   sqrt(8*mSmooth) = 2 * sqrt(2*mSmooth) 
     "not found"
  )

  if ("not found" %in% c(mRange_stm, mRange_output) ) warning("Parameterization not found:", parameterization_input, parameterization_output)

  return(mRange_output)
}
