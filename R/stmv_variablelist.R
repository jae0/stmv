stmv_variablelist = function( p ) {
  
  if (!exists("variables", p)) p$variables = list()
  if (!exists("Y", p$variables)) {
    if (exists("stmv_local_modelformula", p))  {
      if (!is.null(p$stmv_local_modelformula)) {
        if (p$stmv_local_modelformula != "none") {
          oo = all.vars( p$stmv_local_modelformula[[2]] ) 
          if (length(oo) > 0) p$variables$Y = oo 
        }
      }
    }
    if (exists("stmv_global_modelformula", p))  {
      if (!is.null(p$stmv_global_modelformula)) {
        if (p$stmv_global_modelformula != "none") {
          oo = all.vars( p$stmv_global_modelformula[[2]] )
          if (length(oo) > 0) p$variables$Y = oo
        }
      }
    }
  }
  
  if (!exists("Y", p$variables)) p$variables = "not_defined" # this can be called to get covars.. do not stop
   
  p$variables$local_all = NULL  
  p$variables$local_cov = NULL
  if (exists("stmv_local_modelformula", p)) {
    if (!is.null(p$stmv_local_modelformula)) {
      if (p$stmv_local_modelformula != "none") {
        oo = all.vars( p$stmv_local_modelformula )
        if (length(oo) > 0) {
          p$variables$local_all = oo
          p$variables$local_cov = all.vars( p$stmv_local_modelformula[[3]] )
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
          p$variables$global_cov = all.vars( p$stmv_global_modelformula[[3]] )
        }
      }
    }
  }
  p$variables$ALL = NULL  # "all" but Y
  p$variables$ALL = c( p$variables$local_all, p$variables$global_all )
  p$variables$ALL = unique( c( p$variables$ALL, p$variables$LOCS, p$variables$TIME ) )  
  oo = unique( c( grep("cos.w", p$variables$ALL), grep("sin.w", p$variables$ALL) ) )
  if (length(oo) > 0) p$variables$TSvars = p$variables$ALL[oo]  # harmonics
  p$variables$ALL = setdiff( p$variables$ALL, p$variables$TSvars)
  p$variables$COORDS = unique( c(p$variables$LOCS, p$variables$TIME) )
  p$variables$COV = setdiff( p$variables$ALL, p$variables$COORDS )  # non-location and non-time based covariates

  return (p)
}
