
stmv_interpolate_function_lookup = function( modelengine ) {

  # wrapper to copy interpolating function as a generic script
  local_fn = NULL
  local_fn = switch( modelengine,
    akima = stmv__akima,
    bayesx = stmv__bayesx,
    constant = stmv__constant,
    fft = stmv__fft,
    gaussianprocess2Dt = stmv__gaussianprocess2Dt,
    gam = stmv__gam,
    glm = stmv__glm,
    gstat = stmv__gstat,
    krige = stmv__krige,
    kernel = stmv__kernel,
    linear = stmv__linear,
    tps = stmv__tps,
    carstm = stmv__carstm,
    twostep = stmv__twostep
  )
  if ( is.null(local_fn) ) {
    message( "Interpolation module: ", modelengine, " not found." )
    stop ()
  }
  return( local_fn )
}
