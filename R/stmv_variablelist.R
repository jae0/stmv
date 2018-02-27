stmv_variablelist = function( p ) {
  
   
  p$variables$local_all = NULL  
  p$variables$local_cov = NULL
  if (exists("stmv_local_modelformula", p)) {
    if (!is.null(p$stmv_local_modelformula)) {
      if (p$stmv_local_modelformula != "none") {
        oo = all.vars( p$stmv_local_modelformula )
        if (length(oo) > 0) {
          p$variables$local_all = oo
        }
        oo = all.vars( p$stmv_local_modelformula[[3]] )
        if (length(oo) > 0) {
          p$variables$local_cov = oo
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
          p$variables$global_cov = oo
        }
      }
    }
  }
  
  p$variables$ALL = c( p$variables$local_all, p$variables$global_all )
  
  # all external variables (remove harmonics)
  p$variables$ALL_REQUIRED = p$variables$ALL
  oo = unique( c( grep("cos.w", p$variables$ALL), grep("sin.w", p$variables$ALL_REQUIRED) ) )
  if (length(oo) > 0) p$variables$ALL_REQUIRED = p$variables$ALL[-oo]
  p$variables$ALL_REQUIRED = setdiff( p$variables$ALL_REQUIRED, "yr" )  # year is  computed from time index ... not required
  #  varstokeep = unique( c( p$variables$Y, p$variables$LOCS, p$variables$TIME, p$variables$COV ) )
  
  p$variables$COORDS = unique( c(p$variables$LOCS, p$variables$TIME) )
  
  p$variables$COV = setdiff( p$variables$ALL_REQUIRED, c(p$variables$COORDS, p$variables$Y) )
  
  return (p)
}
