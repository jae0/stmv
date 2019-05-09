stmv_scale_filter = function( p, Si ) {

  S = stmv_attach( p$storage.backend, p$ptr$S )
  Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )

  distance_limits = range( c(p$pres*3,  p$stmv_distance_scale ) )   # for range estimate

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
  range_redo = FALSE
  if (!is.finite(vg$range)) {
    range_redo = TRUE
  } else {
    if (vg$range < distance_limits[1] | vg$range > distance_limits[2])  {
      range_redo = TRUE
    }
  }

  ii = NULL
  if (range_redo) {
    vg$flag = "variogram_failure"
    max_dist = max(p$stmv_distance_scale)
    ii = which(
      {abs( Sloc[Si,1] - Sloc[,1] ) <= distance_limits[2]} &
      {abs( Sloc[Si,2] - Sloc[,2] ) <= distance_limits[2]}
    )
    if (length(ii) > 0) {
      range_median = median( S[ii, match("range", p$statsvars)], na.rm=TRUE )
      if (is.finite( range_median)) {
        vg$range =  range_median
      }
    }
  }

  if (!is.finite(vg$range)) vg$range = median(distance_limits)
  if ( vg$range < distance_limits[1] )  vg$range = distance_limits[1]
  if ( vg$range > distance_limits[2] )  vg$range = distance_limits[2]

  # as range is now set, the following becomes fixed
  if (is.null(ii)) {
    ii = which(
      {abs( Sloc[Si,1] - Sloc[,1] ) <= vg$range} &
      {abs( Sloc[Si,2] - Sloc[,2] ) <= vg$range}
    )
  }

  # obtain indices of data locations withing a given spatial range, optimally determined via variogram
  # faster to take a block .. but easy enough to take circles ... trim off corners ..
  nu_redo = FALSE
  if (!is.finite(vg$nu)) {
    nu_redo = TRUE
  } else {
    if (vg$nu < 0.25 | vg$nu > 4 )  {
      nu_redo = TRUE
    }
  }

  if (nu_redo) {
    vg$flag = "variogram_failure"
    nu_median =  median( S[ii, match("nu", p$statsvars)], na.rm=TRUE )
    if (!is.finite(nu_median)) {
      vg$nu = 0.5
    } else {
      if (nu_median < 0.25 | nu_median > 4 )  {
        vg$nu = 0.5
      } else {
        vg$nu = nu_median
      }
    }
  }


  phi_redo = FALSE
  dl = distance_limits/sqrt(8*vg$nu)
  if (!is.finite(vg$phi)) {
    phi_redo = TRUE
  } else {
    if (vg$phi < dl[1] | vg$phi > dl[2] )  {
      phi_redo = TRUE
    }
  }

  if (phi_redo) {
    vg$flag = "variogram_failure"
    phi_median =  median( S[ii, match("phi", p$statsvars)], na.rm=TRUE )
    if (is.finite(phi_median)) {
      vg$phi = phi_median
    } else {
      vg$phi = vg$range/sqrt(8*vg$nu)  # approx rule rule in inla for nu=1 (alpha=2)
    }
  }


  if (!is.finite(vg$varSpatial)) {
    vg$flag = "variogram_failure"
    varSP_median = median( S[ii, match("sdSpatial", p$statsvars)], na.rm=TRUE )^2
    if (is.finite(varSP_median)) vg$varSpatial = varSP_median
  }


  if (!is.finite(vg$varObs)) {
    vg$flag = "variogram_failure"
    varObs_median = median( S[ii, match("sdObs", p$statsvars)], na.rm=TRUE )^2
    if (is.finite(varObs_median)) vg$varObs = varObs_median
  }


  if (!is.finite(vg$varTotal)) {
    vg$flag = "variogram_failure"
    varTot_median = median( S[ii, match("sdTotal", p$statsvars)], na.rm=TRUE )^2
    if (is.finite(varTot_median)) vg$varTotal = varTot_median
  }

  return(vg)
}
