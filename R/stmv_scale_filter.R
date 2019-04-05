stmv_scale_filter = function( p, Si ) {

  S = stmv_attach( p$storage.backend, p$ptr$S )
  Sflag = stmv_attach( p$storage.backend, p$ptr$Sflag )
  E = stmv_error_codes()

  distance_limits = c( p$stmv_distance_prediction, max(p$stmv_distance_scale))  # for range estimate

  vg = list(
    range = S[Si, match("range", p$statsvars)] ,
    nu    = S[Si, match("nu",   p$statsvars)] ,
    phi   = S[Si, match("phi",  p$statsvars)] ,
    varSpatial = S[Si, match("sdSpatial", p$statsvars)]^2 ,
    varObs     = S[Si, match("sdObs",   p$statsvars)]^2 ,
    varTotal   = S[Si, match("sdTotal", p$statsvars)]^2 ,
    flag = "todo"
  )

  #range checks
  if (!is.finite(vg$range)) {
    ii = which(
      {abs( Sloc[Si,1] - Sloc[,1] ) <= min(p$stmv_distance_scale)} &
      {abs( Sloc[Si,2] - Sloc[,2] ) <= min(p$stmv_distance_scale)}
    )
    vg$range =  median( S[ii, match("range", p$statsvars)], na.rm=TRUE )
    vg$flag = "variogram_range_limit"
  } else {
    if (vg$range < distance_limits[1] | vg$range > distance_limits[2])  {
      ii = which(
        {abs( Sloc[Si,1] - Sloc[,1] ) <= min(p$stmv_distance_scale)} &
        {abs( Sloc[Si,2] - Sloc[,2] ) <= min(p$stmv_distance_scale)}
      )
      vg$range =  median( S[ii, match("range", p$statsvars)], na.rm=TRUE )
      vg$flag = "variogram_range_limit"
    }
  }

  # as range is now set, the following becomes fixed
  ii = which(
    {abs( Sloc[Si,1] - Sloc[,1] ) <= vg$range} &
    {abs( Sloc[Si,2] - Sloc[,2] ) <= vg$range}
  )

  # obtain indices of data locations withing a given spatial range, optimally determined via variogram
  # faster to take a block .. but easy enough to take circles ...


  if (!is.finite(vg$nu)) {
    vg$nu = median( S[ii, match("nu", p$statsvars)], na.rm=TRUE )
    vg$flag = "variogram_failure"
  } else {
    if (vg$nu < 0.25) {
      vg$nu = median( S[ii, match("nu", p$statsvars)], na.rm=TRUE )
      vg$flag = "variogram_failure"
    } else if (vg$nu > 5) {
      vg$nu = median( S[ii, match("nu", p$statsvars)], na.rm=TRUE )
      vg$flag = "variogram_failure"
    }
  }

  if (!is.finite(vg$phi)) {
    vg$phi  = median( S[ii, match("phi", p$statsvars)], na.rm=TRUE )
    vg$flag = "variogram_failure"
  } else {
    dl = distance_limits/sqrt(8*vg$nu)
    if (vg$phi < dl[1] ) {
      vg$phi = median( S[ii, match("phi", p$statsvars)], na.rm=TRUE )
      vg$flag = "variogram_failure"
    } else if (vg$phi > dl[2] ) {
      vg$phi = median( S[ii, match("phi", p$statsvars)], na.rm=TRUE )
      vg$flag = "variogram_failure"
    }
  }

  if (!is.finite(vg$varSpatial)) {
    vg$varSpatial = median( S[ii, match("varSpatial", p$statsvars)], na.rm=TRUE )
    vg$flag = "variogram_failure"
  }

  if (!is.finite(vg$varObs)) {
    vg$varObs = median( S[ii, match("varObs", p$statsvars)], na.rm=TRUE )
    vg$flag = "variogram_failure"
  }

  if ({vg$varSpatial + vg$varObs} < 0.01 * vg$varTotal) {
    vg$varObs = median( S[ii, match("varSpatial", p$statsvars)], na.rm=TRUE )
    vg$varSpatial = median( S[ii, match("varObs", p$statsvars)], na.rm=TRUE )
    vg$flag = "variogram_failure"
  }

  return(vg)
}
