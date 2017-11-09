
timeseries_simulator = function( DS="default", n=600, ar.cor=0.9, missing=0.1 ) {
  #\\ data for development/testing
  #\\ ar.cor is the autocorrelation coefficient

  if (DS == "default") {
    x = 1:n
    y = y0 = arima.sim(n=n, model=list( ar=ar.cor ))
    missing.values = sample( x, trunc(n * missing ))  
    y[ missing.values ]  = NA 
    y  = y + rnorm( n, mean=0 , sd=sd( y0,na.rm=TRUE ) )  # add noise
    z = data.frame( x=x, y=y, y0=y0 ) 
    z$decimal_date = 1900:(1900+n-1)
    z$timestamp = lubridate::date_decimal( z$decimal_date  ) 
  }

  if( DS=="sunspots") {
    z = data.frame( cbind( y=c(sunspot.year), decimal_date=time(sunspot.year) ) )
    z$timestamp = lubridate::date_decimal( z$decimal_date  ) 
  }

  if( DS=="sunspots.seasonal") {
    z = data.frame( cbind( y=c(sunspots), decimal_date=time(sunspots)))
    z$timestamp = lubridate::date_decimal( z$decimal_date  ) 
  }

  return (z)
}

