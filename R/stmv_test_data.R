

stmv_test_data = function( datasource="swiss", redo=FALSE, p=NULL ) {

  if ( datasource == "swiss" ) {
    require(geostatsp)
    data(swissRain)
    # head(swissBorder) # SpatialPolygonsDataFrame
    out = as.data.table( swissRain )
    out$swissLandType = sp::over( as(swissRain, "SpatialPoints"), as(swissLandType, "SpatialGridDataFrame"), fn = mean)[,1]
    out$swissLandType  = as.factor( out$swissLandType )
    out$altitude = sp::over( as(swissRain, "SpatialPoints"), as(swissAltitude, "SpatialGridDataFrame"), fn = mean)[,1]
    out$rain_alt = residuals( lm(rain~altitude, out))
    return(out)
  }


  if ( datasource == "meuse" ) {
    require(sp)
    data(meuse)
    out = as.data.table( meuse )
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
    out$dist=as.matrix( dist( as.data.table( out[c("x","y")] )) )
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
    out$dist=as.matrix( dist( as.data.table( out[c("x","y")] )) )
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
    out$dist=as.matrix( dist( as.data.table( out[c("x","y")] )) )
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

    discret = 0.5 # km

    if (is.null(p)) p = stmv_test_data( "aegis.test.parameters")

    pb = bathymetry_parameters( project_class="stmv" )
    B = bathymetry_db( p=pb,  DS="stmv_inputs")$input  # this is a subset of "complete" with depths filtered
    B = planar2lonlat( B, proj.type=pb$aegis_proj4string_planar_km )
    # B = B[,c( "lon", "lat", "z")]
    B = B[ which( B$lon > p$corners$lon[1] & B$lon < p$corners$lon[2]  & B$lat > p$corners$lat[1] & B$lat < p$corners$lat[2] ), ]

    ps = substrate_parameters( project_class="stmv" )
    S = substrate_db( p=ps, DS="lonlat.highres" )

    S = S[ which( S$lon > p$corners$lon[1] & S$lon < p$corners$lon[2] & S$lat > p$corners$lat[1] & S$lat < p$corners$lat[2] ) , ]
    S$substrate.grainsize = S$grainsize
    S = S[ ,c("lon", "lat", "substrate.grainsize" )]
    S$plon = NULL
    S$plat = NULL
    S = lonlat2planar( S, proj.type=p$aegis_proj4string_planar_km )
    S$lon = NULL
    S$lat = NULL
    S$plon = floor(S$plon/discret +1L )* discret
    S$plat = floor(S$plat/discret +1L) * discret
    dups = which(duplicated( paste( S$plon, S$plat) ) )
    if (length(dups) > 0 ) S = S[ -dups , ]

    B = lonlat2planar( B, proj.type=p$aegis_proj4string_planar_km )
    B$lon = NULL
    B$lat = NULL
    B$plon = floor(B$plon/discret +1L) * discret
    B$plat = floor(B$plat/discret +1L) * discret
    dups = which( duplicated( paste( B$plon, B$plat) ) )
    if (length(dups) > 0 ) B = B[ -dups , ]

    out = merge( B[,c("plon", "plat", "z")], S[, c("plon", "plat", "substrate.grainsize")], by=c("plon", "plat"), all.x=TRUE, all.y=TRUE )
    out = planar2lonlat( out, proj.type=p$aegis_proj4string_planar_km )
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

    p = temperature_parameters( spatial_domain="SSE" )

    out = temperature_db( p=p, DS="bottom.all"  )
    out = out[ which( out$lon > p$corners$lon[1] & out$lon < p$corners$lon[2] & out$lat > p$corners$lat[1] & out$lat < p$corners$lat[2] ) , ]
    out = out[ which(out$yr %in% c(1980:2010)), ]
    out$tiyr = lubridate::decimal_date ( out$date )

    # globally remove all unrealistic data
    tokeep = which( out$t >= -3 & out$t <= 25 ) # hard limits
    if (length(tokeep) > 0 ) out = out[ tokeep, ]
    TR = quantile(out$t, probs=c(0.0005, 0.9995), na.rm=TRUE ) # this was -1.7, 21.8 in 2015
    tokeep = which( out$t >=  TR[1] & out$t <=  TR[2] )
    if (length(tokeep) > 0 ) out = out[ tokeep, ]
    tokeep = which( out$z >=  2 ) # ignore very shallow areas ..
    if (length(tokeep) > 0 ) out = out[ tokeep, ]

    out = out[, c("lon", "lat", "t", "z", "date", "tiyr")]
    save(out, file=fn, compress=TRUE)
    return (out)

  }



  if ( datasource == "aegis.bathymetry" ) {

    RLibrary( "aegis.bathymetry" )
    require(sp)
    # fn = system.file( "extdata", "aegis_space_test.Rdata", package="stmv" )
    out = NULL
    fn = project.codedirectory("stmv", "inst", "extdata", "aegis_bathymetry.Rdata")
    if (!redo) {
      if (file.exists(fn)) load(fn)
      return(out)
    }

    pb = bathymetry_parameters( spatial_domain="SSE" )
    B = bathymetry_db( p=pb, DS="complete"  )

    # output locations
    u = list(
        x=seq(min(pb$corners$plon), max(pb$corners$plon), by = pb$pres),
        y=seq(min(pb$corners$plat), max(pb$corners$plat), by = pb$pres)
    )
    u$z = matrix( NA, nrow=length(u$x), ncol=length(u$y) )

    origin=c(min(pb$corners$plon), min(pb$corners$plat) )
    i = as.matrix(array_map( "xy->2", coords=B[,c("plon", "plat")], origin=origin, res=c(pb$pres, pb$pres) ))  # map Stats Locs to Plocs
    u$z[i] = B$z
    # image(u)

    out = fields::interp.surface( u, loc=spatial_grid(p0) ) # linear interpolation
    save(out, file=fn, compress=TRUE)
    return (out)

  }



  if ( datasource == "aegis.test.parameters" ) {
      # dres  is the 15 second grid from CHS  .. default use highest resolution
    p = aegis::spatial_parameters( spatial_domain="testing", aegis_proj4string_planar_km="+proj=utm +ellps=WGS84 +zone=20 +units=km", dres=1/60/4, pres=0.5, lon0=-64, lon1=-62, lat0=44, lat1=45, psignif=2 )
    return(p)
  }

  if ( datasource == "aegis.prediction.locations" ) {
      # dres  is the 15 second grid from CHS  .. default use highest resolution
      ## 1 km resolution!,
    if (is.null(p)) p = stmv_test_data( "aegis.test.parameters")
    LOCS = spatial_grid(p)
    return(LOCS)
  }

  if (0) {

    p = aegis::spatial_parameters( spatial_domain="testing", aegis_proj4string_planar_km="+proj=utm +ellps=WGS84 +zone=20 +units=km", dres=1/60/4, pres=0.5, lon0=-64, lon1=-62, lat0=44, lat1=45, psignif=2 )
    # or:  p = stmv_test_data( "aegis.test.parameters")

    PREDLOCS = spatial_grid(p)

  }


      if (0) {
        #check that working residual is correct
        dd = glm( Sepal.Length ~ Sepal.Width + Petal.Length, data=iris, family=gaussian(link="log"))
        ddp = predict(  dd, type="link", se.fit=FALSE )
        ddr = residuals(dd, type="working") #
        inv = family(dd)$linkinv
        ddo = inv( ddp + ddr )
        plot( ddo ~ iris$Sepal.Length )

        require(mgcv)
        dd = gam( Sepal.Length ~ s(Sepal.Width) + s(Petal.Length), data=iris, family=gaussian(link="log"))
        ddp = predict(  dd, type="link", se.fit=FALSE )
        ddr = residuals(dd, type="working") #
        inv = family(dd)$linkinv
        ddo = inv( ddp + ddr )
        plot( ddo ~ iris$Sepal.Length )

        require(mgcv)
        yy = ifelse( iris$Sepal.Length < median(iris$Sepal.Length), 0, 1 )
        dd = gam( yy ~ s(Sepal.Width), data=iris, family=binomial(link="logit"))
        ddp = predict(  dd, type="link", se.fit=FALSE )
        ddr = residuals(dd, type="working") #
        inv = family(dd)$linkinv
        lnk = family(dd)$linkfun
        ddo = inv(ddp + ddr)
        plot( ddo ~ yy )


        # test a subsample
        testdat = global_model$model
        range(testdat$substrate.grainsize)
        hist(testdat$substrate.grainsize)
        dd = global_model
        ddp = predict(  dd, type="link", se.fit=FALSE )
        ddr = residuals(dd, type="working") #
        inv = family(dd)$linkinv
        ddo = inv( ddp + ddr )
        ddo[ddo > 30] = 30
        plot( ddo ~ testdat$substrate.grainsize )
        hist(ddo, "fd")


        # test a subsample
        testdat = DATA$input[sample(nrow(DATA$input),1000),]
        range(testdat$substrate.grainsize)
        hist(testdat$substrate.grainsize)
        dd = gam( substrate.grainsize ~ s(b.sdSpatial) + s(b.localrange) + s(log(z))+ s(log(dZ))+ s(log(ddZ)) , data=testdat, family=gaussian(link="log"))
        ddp = predict(  dd, type="link", se.fit=FALSE )
        ddr = residuals(dd, type="working") #
        inv = family(dd)$linkinv
        ddo = inv( ddp + ddr )
        plot( ddo ~ testdat$substrate.grainsize )
        ddo[ddo > 10] = 10
        hist(ddo, "fd")

        # test transformations within smooths
        testdat$log_z = log(testdat$z)
        testdat$log_dZ = log(testdat$dZ)
        testdat$log_ddZ = log(testdat$ddZ)
        dd = gam( substrate.grainsize ~ s(b.sdSpatial) + s(b.localrange) + s(log_z)+ s(log_dZ)+ s(log_ddZ) , data=testdat, family=gaussian(link="log"))
        ddp = predict(  dd, type="link", se.fit=FALSE )
        ddr = residuals(dd, type="working") #
        inv = family(dd)$linkinv
        ddo = inv( ddp + ddr )
        plot( ddo ~ testdat$substrate.grainsize )
        ddo[ddo > 10] = 10
        hist(ddo, "fd")
      }


}
