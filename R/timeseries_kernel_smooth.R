  
  timeseries_kernel_smooth = function( OP, vn, smoothing.kernel=kernel( "modified.daniell", c(2,1)), truncation.quants=c(0.005, 0.995) ) {
    
    # constrain data to be smoother and within empirical range (e.g., 99.9% quantiles)
    if ( any( is.finite( OP[,vn] ) ) ) {
        # default is a 2X smoothing kernel: 2 adjacent values and then 1 adjacent on the smoothed series .. i.e. a high-pass filter
        nd = length(smoothing.kernel$coef) - 1 # data points will be trimmed from tails so need to recenter:
        nd0 = length(OP[,vn] ) 
        si = (nd+1):(nd0-nd)
        # plot( eval(vn) ~ time, OP, type="l" )
        OP[si,vn]  = kernapply( as.vector(OP[,vn] ), smoothing.kernel ) 
        # lines( eval(vn) ~ time, OP, col="green" )
      
      # constrain range of predicted data to the input data range
      TR =  quantile( x$t, probs=truncation.quants, na.rm=TRUE  )
      # TR[1] = max( TR[1], -3)
      # TR[2] = min( TR[2], 30)
      toolow = which( OP[,vn]  < TR[1] )
      if ( length(toolow) > 0 )  OP[,vn] [toolow] = TR[1]
      toohigh = which( OP[,vn]  > TR[2] )
      if ( length(toohigh) > 0 ) OP[,vn] [toohigh] = TR[2]
    }

    debug = FALSE
    if(debug) {

      n=100
      ar.cor=0.9
      missing=0.1 
      x = 1:n
      y = y0 = arima.sim(n=n, model=list( ar=ar.cor ))
      missing.values = sample( x, trunc(n * missing ))  
      y[ missing.values ]  = NA 
      y  = y + rnorm( n, mean=0 , sd=sd( y0,na.rm=TRUE ) )  # add noise
      z = data.frame( x=x, y=y, y0=y0 ) 

      z$decimal_date = 1900:(1900+n-1)
      z$timestamp = lubridate::date_decimal( z$decimal_date  ) 

      OP = z
      OP$date = OP$timestamp
      OP$time = lubridate::decimal_date( OP$timestamp) 
      OP = OP[ order( OP$time ) ,]
      
      # STL method
      OPts = ts( oH1[,vn] , start=c( min(p$tyears), 1), frequency=12 )  # testing oH1's fit
      plot.default(  OPts, type="l" )

      OPstl = stl( OPts, s.window=4) # seasonal decomposition using loess 
      plot(OPstl)

      # StructTS method
      OPstr  = StructTS( OPts, type = "BSM") ### much slower ... used only as a daignostic tool checking stability
      OPstrsm = tsSmooth( OPstr )
      plot( OPstrsm )
   
      spectrum( OPts )
      kn = kernel( "modified.daniell", c(2,1))

      opop = kernapply( as.vector(oH1[,vn] ), kn)
      plot(oH1[,vn], type="l")
      lines(opop, type="l", col="orange")

      ppp = kernapply( OPts, kn)
      plot.default(  OPts, type="l", col="green" )
      lines( ppp, col = "red")  

      ooo = spectrum( OPts, kn, plot=FALSE )
      plot( ooo, plot.type = "marginal") # the default type
      plot( ooo, plot.type = "coherency")
      plot( ooo, plot.type = "phase")


      require(forecast)
      ppp = Arima( OPts )
    }

  }

