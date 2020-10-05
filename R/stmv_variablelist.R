stmv_variablelist = function( p ) {

  p$stmv_variables$COORDS = unique( c(p$stmv_variables$LOCS, p$stmv_variables$TIME) )

  p$stmv_variables$local_all = NULL
  p$stmv_variables$local_cov = NULL
  if (exists("stmv_local_modelformula", p)) {
    if (!is.null(p$stmv_local_modelformula)) {
      if (p$stmv_local_modelformula != "none") {
        oo = all.vars( p$stmv_local_modelformula )
        if (length(oo) > 0) {
          p$stmv_variables$local_all = unique(c(p$stmv_variables$local_all, oo))
          uu = all.vars( p$stmv_local_modelformula[[3]] )
          if (length(uu) > 0) {
            pp = unique( c( grep("cos.w", uu), grep("sin.w", uu), which(uu %in% c(p$stmv_variables$LOCS, "yr") ) ) )
            if (length(pp) > 0) uu = uu[-pp]
            p$stmv_variables$local_cov = unique(c(p$stmv_variables$local_cov, uu))
          }
        } else if (length(oo) == 0) {
          if ( p$stmv_local_modelengine=="carstm" ) {
            if (p$stmv_local_all_carstm !="") p$stmv_variables$local_all = p$stmv_local_all_carstm
            if (p$stmv_local_covariates_carstm !="") p$stmv_variables$local_cov = p$stmv_local_covariates_carstm
          }
        }
      }
    }
  }

  if (exists("stmv_local_modelformula_time", p)) {
    if (!is.null(p$stmv_local_modelformula_time)) {
      if (p$stmv_local_modelformula_time != "none") {
        oo = all.vars( p$stmv_local_modelformula_time )
        if (length(oo) > 0) {
          p$stmv_variables$local_all = unique(c(p$stmv_variables$local_all, oo))
        }
        oo = all.vars( p$stmv_local_modelformula_time[[3]] )
        if (length(oo) > 0) {
          pp = unique( c( grep("cos.w", oo), grep("sin.w", oo), which(oo %in% c(p$stmv_variables$LOCS, "yr") ) ) )
          if (length(pp) > 0) oo = oo[-pp]
          p$stmv_variables$local_cov = unique(c(p$stmv_variables$local_cov, oo))
        }
      }
    }
  }


  if (exists("stmv_local_modelformula_space", p)) {
    if (!is.null(p$stmv_local_modelformula_space)) {
      if (p$stmv_local_modelformula_space != "none") {
        oo = all.vars( p$stmv_local_modelformula_space )
        if (length(oo) > 0) {
          p$stmv_variables$local_all = unique(c(p$stmv_variables$local_all, oo))
        }
        oo = all.vars( p$stmv_local_modelformula_space[[3]] )
        if (length(oo) > 0) {
          pp = unique( c( grep("cos.w", oo), grep("sin.w", oo), which(oo %in% c(p$stmv_variables$LOCS, "yr") ) ) )
          if (length(pp) > 0) oo = oo[-pp]
          p$stmv_variables$local_cov = unique(c(p$stmv_variables$local_cov, oo))
        }
      }
    }
  }


  p$stmv_variables$global_all = NULL
  p$stmv_variables$global_cov = NULL
  if (exists("stmv_global_modelformula", p)) {
    if (!is.null(p$stmv_global_modelformula)) {
      if (p$stmv_global_modelformula != "none") {
        oo = all.vars( p$stmv_global_modelformula )
        if (length(oo) > 0) {
          p$stmv_variables$global_all = oo
        }
        oo = all.vars( p$stmv_global_modelformula[[3]] )
        if (length(oo) > 0) {
          pp = unique( c( grep("cos.w", oo), grep("sin.w", oo), which(oo %in% c(p$stmv_variables$LOCS, "yr") ) ) )
          if (length(pp) > 0) oo = oo[-pp]
          p$stmv_variables$global_cov = oo
        }
      }
    }
  }

  p$stmv_variables$ALL = unique( c( p$stmv_variables$local_all, p$stmv_variables$global_all ) )

  # all external stmv_variables (remove harmonics)
  p$stmv_variables$ALL_REQUIRED = p$stmv_variables$ALL
  oo = unique( c( grep("cos.w", p$stmv_variables$ALL_REQUIRED), grep("sin.w", p$stmv_variables$ALL_REQUIRED), which(p$stmv_variables$ALL_REQUIRED=="yr") ) )
  if (length(oo) > 0) p$stmv_variables$ALL_REQUIRED = p$stmv_variables$ALL[-oo]
  # year is  computed from time index ... not required

  p$stmv_variables$COV = setdiff( p$stmv_variables$ALL_REQUIRED, c(p$stmv_variables$COORDS, p$stmv_variables$Y) )
  # if (length(p$stmv_variables$COV) ==0) p$stmv_variables$COV = NULL

  return (p)
}
