
stmv_interpolation_function = function( modelengine ) {

  # wrapper to copy interpolating function as a generic script
  local_fn = NULL
  local_fn = switch( modelengine,
    bayesx = stmv__bayesx,
    gaussianprocess2Dt = stmv__gaussianprocess2Dt,
    gam = stmv__gam,
    glm = stmv__glm,
    gstat = stmv__gstat,
    krige = stmv__krige,
    fft = stmv__fft,
    tps = stmv__tps,
    twostep = stmv__twostep
  )
  return( local_fn )
}
