

stmv_test_data = function( datasource="swiss", redo=FALSE, p=NULL ) {

  if ( datasource == "swiss" ) {
    require(geostatsp)
    data(swissRain)
    # head(swissBorder) # SpatialPolygonsDataFrame
    out = as.data.frame( swissRain )
    out$swissLandType = sp::over( as(swissRain, "SpatialPoints"), as(swissLandType, "SpatialGridDataFrame"), fn = mean)[,1]
    out$swissLandType  = as.factor( out$swissLandType )
    out$altitude = sp::over( as(swissRain, "SpatialPoints"), as(swissAltitude, "SpatialGridDataFrame"), fn = mean)[,1]
    out$rain_alt = residuals( lm(rain~altitude, out))
    return(out)
  }


  if ( datasource == "meuse" ) {
    require(sp)
    data(meuse)
    out = as.data.frame( meuse )
    out$z = residuals( lm( log( meuse$zinc ) ~ sqrt( meuse$dist ), data=out ) )
    return(out)
  }


  if ( datasource == "binomial" ) {
    out = stmv_test_data( "swiss")
    # mimic binomial data from rain values
    out = as.list(out)
    out$N=length(out$rain)
    out$Np= 0
    out$ncases=floor(rnorm(length(out$rain), mean=0.5, sd=0.1) * out$rain)
    out$ntot=floor( out$rain)
    out$X=rep(0,out$N)
    out$dist=as.matrix( dist( as.data.frame( out[c("x","y")] )) )
    out$eps = 1e-6
    out$COVFN=1
    return(out)
  }


  if ( datasource == "poisson" ) {
    out = stmv_test_data( "swiss")
    # mimic poisson data from rain values
    out = as.list(out)
    out$N=length(out$rain)
    out$Np=0
    out$ncases=floor(rnorm(length(out$rain), mean=0.5, sd=0.1) * out$rain)
    out$ntot=floor( out$rain)
    out$X=rep(0,out$N)
    out$dist=as.matrix( dist( as.data.frame( out[c("x","y")] )) )
    out$eps = 1e-6
    out$COVFN=1
    return(out)
  }



  if ( datasource == "gaussian" ) {
    out = stmv_test_data( "swiss")
    # mimic binomial data from rain values
    out = as.list(out)
    out$N=length(out$rain)
    out$Np=0
    out$ncases=floor(rnorm(length(out$rain), mean=0.5, sd=0.1) * out$rain)
    out$ntot=floor( out$rain)
    out$Y=out$rain_alt
    out$X=rep(0,out$N)
    out$dist=as.matrix( dist( as.data.frame( out[c("x","y")] )) )
    out$eps = 1e-6
    out$COVFN=1
    return(out)
  }


  if ( datasource == "aegis.space" ) {

    RLibrary( "aegis.bathymetry", "aegis.substrate" )
    require(sp)
    # fn = system.file( "extdata", "aegis_space_test.Rdata", package="stmv" )
    out = NULL
    fn = project.codedirectory("stmv", "inst", "extdata", "aegis_space_test.Rdata")
    if (!redo) {
      if (file.exists(fn)) load(fn)
      return(out)
    }

    discret = 0.1 # km

    if (is.null(p)) p = stmv_test_data( "aegis.test.paramaters")

    pb = bathymetry_parameters( project.mode="stmv" )
    B = bathymetry.db( p=pb,  DS="stmv.inputs")$input  # this is a subset of "complete" with depth filtered
    B = planar2lonlat( B, proj.type=pb$internal.crs )
    B = B[,c( "lon", "lat", "z")]
    B = B[ which( B$lon > p$corners$lon[1] & B$lon < p$corners$lon[2]  & B$lat > p$corners$lat[1] & B$lat < p$corners$lat[2] ), ]

    ps = substrate_parameters( project.mode="stmv" )
    S = substrate.db( p=ps, DS="lonlat.highres" )

    S = S[ which( S$lon > p$corners$lon[1] & S$lon < p$corners$lon[2] & S$lat > p$corners$lat[1] & S$lat < p$corners$lat[2] ) , ]
    S$substrate.grainsize = S$grainsize
    S = S[ ,c("lon", "lat", "substrate.grainsize" )]
    S$plon = NULL
    S$plat = NULL
    S = lonlat2planar( S, proj.type=p$internal.crs )
    S$lon = NULL
    S$lat = NULL
    S$plon = floor(S$plon/discret) * discret
    S$plat = floor(S$plat/discret) * discret
    dups = duplicates.toremove( paste( S$plon, S$plat) )
    if (length(dups) > 0 ) S = S[ -dups , ]

    B = lonlat2planar( B, proj.type=p$internal.crs )
    B$lon = NULL
    B$lat = NULL
    B$plon = floor(B$plon/discret) * discret
    B$plat = floor(B$plat/discret) * discret
    dups = duplicates.toremove( paste( B$plon, B$plat) )
    if (length(dups) > 0 ) B = B[ -dups , ]

    out = merge( B[,c("plon", "plat", "z")], S[, c("plon", "plat", "substrate.grainsize")], by=c("plon", "plat"), all.x=TRUE, all.y=TRUE )
    out = planar2lonlat( out, proj.type=p$internal.crs )
    out$plon = out$plat = NULL

    save(out, file=fn, compress=TRUE)
    return (out)
  }


  if ( datasource == "aegis.spacetime" ) {

    RLibrary( "aegis.bathymetry", "aegis.substrate", "aegis.temperature" )
    require(sp)
    # fn = system.file( "extdata", "aegis_space_test.Rdata", package="stmv" )
    out = NULL
    fn = project.codedirectory("stmv", "inst", "extdata", "aegis_spacetime_test.Rdata")
    if (!redo) {
      if (file.exists(fn)) load(fn)
      return(out)
    }

    p = temperature_parameters( spatial.domain="SSE" )

    out = temperature.db( p=p, DS="bottom.all"  )
    out = out[ which( out$lon > p$corners$lon[1] & out$lon < p$corners$lon[2] & out$lat > p$corners$lat[1] & out$lat < p$corners$lat[2] ) , ]
    out = out[ which(out$yr %in% c(1980:2010)), ]
    out$tiyr = lubridate::decimal_date ( out$date )

    # globally remove all unrealistic data
    keep = which( out$t >= -3 & out$t <= 25 ) # hard limits
    if (length(keep) > 0 ) out = out[ keep, ]
    TR = quantile(out$t, probs=c(0.0005, 0.9995), na.rm=TRUE ) # this was -1.7, 21.8 in 2015
    keep = which( out$t >=  TR[1] & out$t <=  TR[2] )
    if (length(keep) > 0 ) out = out[ keep, ]
    keep = which( out$z >=  2 ) # ignore very shallow areas ..
    if (length(keep) > 0 ) out = out[ keep, ]

    out = out[, c("lon", "lat", "t", "z", "date", "tiyr")]
    save(out, file=fn, compress=TRUE)
    return (out)

  }


  if ( datasource == "aegis.test.paramaters" ) {
      # dres  is the 15 second grid from CHS  .. default use highest resolution
    p = aegis::spatial_parameters( spatial.domain="testing", internal.crs="+proj=utm +ellps=WGS84 +zone=20 +units=km", dres=1/60/4, pres=0.5, lon0=-64, lon1=-62, lat0=44, lat1=45, psignif=2 )
    return(p)
  }

  if ( datasource == "aegis.prediction.locations" ) {
      # dres  is the 15 second grid from CHS  .. default use highest resolution
      ## 1 km resolution!,
    if (is.null(p)) p = stmv_test_data( "aegis.test.paramaters")
    LOCS = spatial_grid(p)
    return(LOCS)
  }

  if (0) {

    p = aegis::spatial_parameters( spatial.domain="testing", internal.crs="+proj=utm +ellps=WGS84 +zone=20 +units=km", dres=1/60/4, pres=0.5, lon0=-64, lon1=-62, lat0=44, lat1=45, psignif=2 )
    # or:  p = stmv_test_data( "aegis.test.paramaters")

    PREDLOCS = spatial_grid(p)

  }

}
