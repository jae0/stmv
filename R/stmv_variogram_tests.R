
stmv_variogram_tests = function() {

  # debugging / comparison of results

  loadfunctions( c( "aegis", "stmv" ))

  XYZ = NULL
  
  # xyz = stmv_test_data( datasource="swiss" )  # ? geostatsp removed the data ?
  # mz = log( xyz$rain )
  
  xyz = stmv_test_data( datasource="meuse" )  
  mz = log( xyz$lead )
  
  xy = xyz[, c("x", "y")]
  mm = lm( mz ~ 1 )
  z = residuals( mm)

  plotdata=TRUE
  nbreaks = c( 16, 32, 64, 13, 17, 21, 33, 49 )
  distance_cutoff=NA
  family=gaussian(link="identity")
  range_correlation=0.1

  # scale_distances = FALSE

  scale_distances = TRUE
  discretized_n=NULL
  range_correlation=0.1
  stmv_autocorrelation_fft_taper=0

  require( RandomFields )
  RFoptions(install="no")


  out = NULL

  vns = c(  "nu", "phi", "varSpatial", "varObs", "localrange", "phi_ok" ) 

  methods = c(
    "fft", "optim", "fast.recursive", "fields", 
    "geoR", "geoR.ML", "gstat", "RandomFields", "CompRandFld",
    "bayesx", "inla", "spBayes"
  )

  # not working/incomplete: 
  #  "LaplacesDemon", "stan", "geostatsp_poisson", "jags.exponential"
  # , "julia_turing" 


  for (method in methods ) {
 
    mb = microbenchmark::microbenchmark( {
      gr = stmv_variogram( xy, z, methods=method, plotdata=FALSE )
    }, times= 10 )   

    # gr = try( stmv_variogram( xy, z, methods=method, plotdata=TRUE  ))

    if (!inherits( gr, "try-error")) {
      newdat = cbind( method, summary(mb)$mean, t(gr[[method]][vns] ) )
      out = rbind( out, newdat )
      print( newdat)
    }
  } 

  names(out) = c( "method", "time_mean", vns )

  print(out)
 

  

    # testing and debugging

    loadfunctions( c("aegis", "stmv"))
    RLibrary(c ("fields", "MBA", "geoR") )

    nx = ny = 128
    nx = ny = 64

    if (0) {
      XYZ = stmv_test_data( datasource="swiss" )
      mz = log( XYZ$rain )
      mm = lm( mz ~ 1 )
      XYZ$z = residuals( mm)
      XYZ=XYZ[c("x","y","z")]
      dev.new()
      oo = stmv_variogram_fft( XYZ[c("x","y","z")], nx=nx, ny=ny, nbreaks=nx, plotdata=TRUE, add.interpolation=TRUE,
        stmv_fft_filter="matern_tapered_modelled", stmv_autocorrelation_fft_taper=0.05 )

      oo = stmv_variogram_fft( XYZ[c("x","y","z")], nx=nx, ny=ny, nbreaks=nx, plotdata=TRUE, add.interpolation=TRUE,
        stmv_fft_filter="matern_tapered_empirical", stmv_autocorrelation_fft_taper=0.09 )

      gr = stmv_variogram( XYZ[, c("x", "y")], XYZ[,"z"], methods="fft", plotdata=TRUE ) # fft/nl least squares
    }


    if (0) {
      XYZ = stmv_test_data( datasource="meuse" )
      XYZ$z = log(XYZ$elev)
      XYZ=XYZ[, c("x","y","z") ]
      dev.new()
      oo = stmv_variogram_fft( XYZ[c("x","y","z")], nx=nx, ny=ny, nbreaks=nx, plotdata=TRUE,  add.interpolation=TRUE,
        stmv_fft_filter="matern_tapered_modelled", stmv_autocorrelation_fft_taper=0.25  )

      oo = stmv_variogram_fft( XYZ[c("x","y","z")], nx=nx, ny=ny, nbreaks=nx, plotdata=TRUE,  add.interpolation=TRUE,
        stmv_fft_filter="matern_tapered_empirical",  stmv_autocorrelation_fft_taper=0.05 )

      gr = stmv_variogram( XYZ[, c("x", "y")], XYZ[,"z"], methods="fft", plotdata=TRUE )
    }

    if (0) {
      str(fields::RMprecip)
      fit <- smooth.2d( RMprecip$y, x=RMprecip$x, theta=.5)
      image( fit )
      points( RMprecip$x, pch=".")

      XYZ = as.data.frame( cbind( RMprecip$x, RMprecip$y ) )
      names(XYZ) =c( "x","y","z")
      oo = stmv_variogram_fft( XYZ[c("x","y","z")], nx=nx, ny=ny, nbreaks=nx, plotdata=TRUE,  add.interpolation=TRUE,
          stmv_fft_filter="matern_tapered_modelled", stmv_autocorrelation_fft_taper=0.2  )

      oo = stmv_variogram_fft( XYZ[c("x","y","z")], nx=nx, ny=ny, nbreaks=nx, plotdata=TRUE,  add.interpolation=TRUE,
          stmv_fft_filter="matern_tapered_empirical",  stmv_autocorrelation_fft_taper=0.07  )

      gr = stmv_variogram( XYZ[, c("x", "y")], XYZ[,"z"], methods="fft", plotdata=TRUE )
    }

    variomethod="geoR"
    variomethod="gstat"
    variomethod="fft"
    variomethod="inla"

    gr = stmv_variogram( XYZ[, c("x", "y")], XYZ[,"z"], methods=variomethod, plotdata=TRUE ) # ml via profile likelihood
    nu = gr[[variomethod]]$nu
    phi = gr[[variomethod]]$phi
    fit  =  Krig(XYZ[, c("x", "y")], XYZ[,"z"], theta=phi)
    dev.new(); surface( fit, type="C") # look at the surface

    mba.int  =  mba.surf( XYZ[, c("x", "y", "z")], nx, ny, extend=TRUE)$xyz.est
    dev.new(); surface(mba.int, xaxs="r", yaxs="r")

    str(fields::RMprecip)
    fit <- smooth.2d( RMprecip$y, x=RMprecip$x, theta=.25)
    dev.new()
    image( fit )
    points( RMprecip$x, pch=".")


} 
