

stmv_test_data = function( datasource="swiss" ) {

  if ( datasource == "swiss" ) {
    require(geostatsp)
    data(swissRain) 
    # head(swissBorder) # SpatialPolygonsDataFrame
    xy = as.data.frame( swissRain )
    xy$swissLandType = sp::over( as(swissRain, "SpatialPoints"), as(swissLandType, "SpatialGridDataFrame"), fn = mean)[,1]
    xy$swissLandType  = as.factor( xy$swissLandType )
    xy$altitude = sp::over( as(swissRain, "SpatialPoints"), as(swissAltitude, "SpatialGridDataFrame"), fn = mean)[,1]
    xy$rain_alt = residuals( lm(rain~altitude, xy))
    return(xy)
  }


  if ( datasource == "meuse" ) {
    require(sp)
    data(meuse)
    xy = as.data.frame( meuse )
    xy$z = residuals( lm( log( meuse$zinc ) ~ sqrt( meuse$dist ), data=xy ) )
    return(xy)
  }


  if ( datasource == "binomial" ) {
    xy = stmv_test_data( "swiss")
    # mimic binomial data from rain values
    xy = as.list(xy)
    xy$N=length(xy$rain)
    xy$Np= 0  
    xy$ncases=floor(rnorm(length(xy$rain), mean=0.5, sd=0.1) * xy$rain) 
    xy$ntot=floor( xy$rain) 
    xy$X=rep(0,xy$N)
    xy$dist=as.matrix( dist( as.data.frame( xy[c("x","y")] )) )
    xy$eps = 1e-6
    xy$COVFN=1
    return(xy)
  }


  if ( datasource == "poisson" ) {
    xy = stmv_test_data( "swiss")
    # mimic poisson data from rain values
    xy = as.list(xy)
    xy$N=length(xy$rain)
    xy$Np=0
    xy$ncases=floor(rnorm(length(xy$rain), mean=0.5, sd=0.1) * xy$rain) 
    xy$ntot=floor( xy$rain) 
    xy$X=rep(0,xy$N)
    xy$dist=as.matrix( dist( as.data.frame( xy[c("x","y")] )) )
    xy$eps = 1e-6
    xy$COVFN=1
    return(xy)
  }



  if ( datasource == "gaussian" ) {
    xy = stmv_test_data( "swiss")
    # mimic binomial data from rain values
    xy = as.list(xy)
    xy$N=length(xy$rain)
    xy$Np=0
    xy$ncases=floor(rnorm(length(xy$rain), mean=0.5, sd=0.1) * xy$rain) 
    xy$ntot=floor( xy$rain) 
    xy$Y=xy$rain_alt
    xy$X=rep(0,xy$N)
    xy$dist=as.matrix( dist( as.data.frame( xy[c("x","y")] )) )
    xy$eps = 1e-6
    xy$COVFN=1
    return(xy)


  }


}
