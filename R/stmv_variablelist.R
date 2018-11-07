stmv_variablelist = function( p ) {

  p$variables$COORDS = unique( c(p$variables$LOCS, p$variables$TIME) )

  p$variables$local_all = NULL
  p$variables$local_cov = NULL
  if (exists("stmv_local_modelformula", p)) {
    if (!is.null(p$stmv_local_modelformula)) {
      if (p$stmv_local_modelformula != "none") {
        oo = all.vars( p$stmv_local_modelformula )
        if (length(oo) > 0) {
          p$variables$local_all = unique(c(p$variables$local_all, oo))
        }
        oo = all.vars( p$stmv_local_modelformula[[3]] )
        if (length(oo) > 0) {
          pp = unique( c( grep("cos.w", oo), grep("sin.w", oo), which(oo %in% c(p$variables$LOCS, "yr") ) ) )
          if (length(pp) > 0) oo = oo[-pp]
          p$variables$local_cov = unique(c(p$variables$local_cov, oo))
        }
      }
    }
  }

  if (exists("stmv_local_modelformula_time", p)) {
    if (!is.null(p$stmv_local_modelformula_time)) {
      if (p$stmv_local_modelformula_time != "none") {
        oo = all.vars( p$stmv_local_modelformula_time )
        if (length(oo) > 0) {
          p$variables$local_all = unique(c(p$variables$local_all, oo))
        }
        oo = all.vars( p$stmv_local_modelformula_time[[3]] )
        if (length(oo) > 0) {
          pp = unique( c( grep("cos.w", oo), grep("sin.w", oo), which(oo %in% c(p$variables$LOCS, "yr") ) ) )
          if (length(pp) > 0) oo = oo[-pp]
          p$variables$local_cov = unique(c(p$variables$local_cov, oo))
        }
      }
    }
  }


  if (exists("stmv_local_modelformula_space", p)) {
    if (!is.null(p$stmv_local_modelformula_space)) {
      if (p$stmv_local_modelformula_space != "none") {
        oo = all.vars( p$stmv_local_modelformula_space )
        if (length(oo) > 0) {
          p$variables$local_all = unique(c(p$variables$local_all, oo))
        }
        oo = all.vars( p$stmv_local_modelformula_space[[3]] )
        if (length(oo) > 0) {
          pp = unique( c( grep("cos.w", oo), grep("sin.w", oo), which(oo %in% c(p$variables$LOCS, "yr") ) ) )
          if (length(pp) > 0) oo = oo[-pp]
          p$variables$local_cov = unique(c(p$variables$local_cov, oo))
        }
      }
    }
  }


  p$variables$global_all = NULL
  p$variables$global_cov = NULL
  if (exists("stmv_global_modelformula", p)) {
    if (!is.null(p$stmv_global_modelformula)) {
      if (p$stmv_global_modelformula != "none") {
        oo = all.vars( p$stmv_global_modelformula )
        if (length(oo) > 0) {
          p$variables$global_all = oo
        }
        oo = all.vars( p$stmv_global_modelformula[[3]] )
        if (length(oo) > 0) {
          pp = unique( c( grep("cos.w", oo), grep("sin.w", oo), which(oo %in% c(p$variables$LOCS, "yr") ) ) )
          if (length(pp) > 0) oo = oo[-pp]
          p$variables$global_cov = oo
        }
      }
    }
  }

  p$variables$ALL = c( p$variables$local_all, p$variables$global_all )

  # all external variables (remove harmonics)
  p$variables$ALL_REQUIRED = p$variables$ALL
  oo = unique( c( grep("cos.w", p$variables$ALL_REQUIRED), grep("sin.w", p$variables$ALL_REQUIRED), which(p$variables$ALL_REQUIRED=="yr") ) )
  if (length(oo) > 0) p$variables$ALL_REQUIRED = p$variables$ALL[-oo]
  # year is  computed from time index ... not required

  p$variables$COV = setdiff( p$variables$ALL_REQUIRED, c(p$variables$COORDS, p$variables$Y) )
  # if (length(p$variables$COV) ==0) p$variables$COV = NULL

  return (p)
}
