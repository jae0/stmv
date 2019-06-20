

stmv_test_data = function( datasource="swiss" ) {

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
    fn = system.file( "extdata", "aegis_space_test.Rdata", package="stmv" )
    if (!redo) {

    }
    RLibrary( "aegis.bathymetry", "aegis.substrate" )
    require(sp)

    p = bathymetry_parameters( project.mode="stmv" )
    Bout = bathymetry.db( p=p,  DS="stmv.inputs")$input  # this is a subset of "complete" with depth filtered
    Bout = Bout[ which( Bout$plon <60 & Bout$plon > -60 & Bout$plat < 60 & Bout$plat > -60 ), ]
    Bout = planar2lonlat( Bout, proj.type=p$internal.crs )
    Bout = Bout[,c( "lon", "lat", "z")]
    Blon = range(B$lon)
    Blat = range(B$lat)

    p = substrate_parameters( project.mode="stmv" )
    S = substrate.db( p=p, DS="lonlat.highres" )

    S = S[ which( S$lon > Blon[1] & S$lon < Blon[2] & S$lat > Blat[1] & S$lat < Blat[2] ) , ]
    S$substrate.grainsize = S$grainsize
    S = S[ ,c("lon", "lat", "substrate.grainsize" )]

    # discretize to speed up the rest
    discret = 0.01 # angular
    S$lon = floor(S$lon/discret) * discret
    S$lat = floor(S$lat/discret) * discret
    S = S[ unique( paste( S$lon, S$lat) ), ]

    Sout = substrate.db ( p=p, DS="complete" )
    Sind = which(names(Sout) %in% "substrate.grainsize")
    if (nrow(Sout) != nrow(Bout)) stop( "Row numbers between bathymetry and substrate databases differ")
    Sout = as.data.frame( Sout[,Sind])
    names(Sout) = vars_required[toadd]
    Bout = cbind( Bout, Sout)

    return (out)
  }


  if ( datasource == "aegis.spacetime" ) {
    # small peice of temperature data

    RLibrary( "aegis.bathymetry", "aegis.substrate" )

    p = spatial_parameters( p=p, spatial.domain="SSE" )

    vars_required = c( "lon", "lat", "z", "substrate.grainsize" )
    p = bathymetry_parameters( project.mode="stmv" )
    Bout = bathymetry.db( p=p,  DS="stmv.inputs")$input  # this is a subset of "complete" with depth filtered
    Bout = Bout[ which( Bout$plon <60 & Bout$plon > -60 & Bout$plat < 60 & Bout$plat > -60 ), ]
    Bout = planar2lonlat( Bout, proj.type=p$internal.crs )
    Bout = Bout[,c( "lon", "lat", "z")]
    Blon = range(B$lon)
    Blat = range(B$lat)

    p = substrate_parameters( project.mode="stmv" )
    S = substrate.db( p=p, DS="lonlat.highres" )

    S = S[ which( S$lon > Blon[1] & S$lon < Blon[2] & S$lat > Blat[1] & S$lat < Blat[2] ) , ]
    S$substrate.grainsize = S$grainsize
    S = S[ ,c("lon", "lat", "substrate.grainsize" )]

    # discretize to speed up the rest
    discret = 0.01 # angular
    S$lon = floor(S$lon/discret) * discret
    S$lat = floor(S$lat/discret) * discret
    S = S[ unique( paste( S$lon, S$lat) ), ]

    Sout = substrate.db ( p=p, DS="complete" )
    Sind = which(names(Sout) %in% "substrate.grainsize")
    if (nrow(Sout) != nrow(Bout)) stop( "Row numbers between bathymetry and substrate databases differ")
    Sout = as.data.frame( Sout[,Sind])
    names(Sout) = vars_required[toadd]
    Bout = cbind( Bout, Sout)
    Sout = NULL; gc()
    if (length(p$variables$COV)==1) {
      covs = list( Bout[,p$variables$COV] )
      names(covs) = p$variables$COV
      OUT  = list( LOCS = Bout[,p$variables$LOCS], COV=covs )
    } else {
      OUT  = list( LOCS = Bout[,p$variables$LOCS], COV=as.list( Bout[,p$variables$COV] ) )
    }

    B = temperature.db( p=p, DS="bottom.all"  )
    B = B[ which(B$yr %in% p$yrs), ]
    B$tiyr = lubridate::decimal_date ( B$date )

    # globally remove all unrealistic data
    keep = which( B$t >= -3 & B$t <= 25 ) # hard limits
    if (length(keep) > 0 ) B = B[ keep, ]
    TR = quantile(B$t, probs=c(0.0005, 0.9995), na.rm=TRUE ) # this was -1.7, 21.8 in 2015
    keep = which( B$t >=  TR[1] & B$t <=  TR[2] )
    if (length(keep) > 0 ) B = B[ keep, ]
    keep = which( B$z >=  2 ) # ignore very shallow areas ..
    if (length(keep) > 0 ) B = B[ keep, ]

    locsmap = match(
      stmv::array_map( "xy->1", B[,c("plon","plat")], gridparams=p$gridparams ),
      stmv::array_map( "xy->1", bathymetry.db(p=p, DS="baseline"), gridparams=p$gridparams ) )

    newvars = setdiff(p$variables$COV, names(B) )
    if (length(newvars) > 0) {
      sn = Bout[locsmap,newvars]
      if (ncol(sn) > 0) {
        B = cbind( B,  sn )
      }
    }

    varstokeep = unique( c( p$variables$Y, p$variables$LOCS, p$variables$TIME, p$variables$COV ) )
    B = B[,varstokeep]

    return (list(input=B, output=OUT))
  }



}
