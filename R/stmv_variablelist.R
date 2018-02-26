stmv_variablelist = function( p ) {
  
  if (!exists("variables", p)) p$variables = list()
  if (!exists("Y", p$variables)) {
    if (exists("stmv_local_modelformula", p))  p$variables$Y = all.vars( p$stmv_local_modelformula[[2]] ) 
    if (exists("stmv_global_modelformula", p)) p$variables$Y = all.vars( p$stmv_global_modelformula[[2]] )
  }
  if (!exists("Y", p$variables)) stop("p$variables$Y is not defined ... it is required")

  p$variables$local_all = NULL  
  p$variables$local_cov = NULL
  if (exists("stmv_local_modelformula", p)) {
    p$variables$local_all = all.vars( p$stmv_local_modelformula )
    p$variables$local_cov = all.vars( p$stmv_local_modelformula[[3]] )
  }
  p$variables$global_all = NULL  
  p$variables$global_cov = NULL
  if (exists("stmv_global_modelformula", p)) {
    p$variables$global_all = all.vars( p$stmv_global_modelformula )
    p$variables$global_cov = all.vars( p$stmv_global_modelformula[[3]] )
  }
  p$variables$ALL = NULL
  p$variables$ALL = c( p$variables$local_all, p$variables$global_all )
  p$variables$ALL = unique( c( p$variables$ALL, p$variables$LOCS, p$variables$TIME ) )  
  p$variables$TSvars = p$variables$ALL[ unique( c( grep("cos.w", p$variables$ALL), grep("sin.w", p$variables$ALL) )  )]  # harmonics
  p$variables$ALL = setdiff( p$variables$ALL, p$variables$TSvars)
  p$variables$coordinates = unique( c(p$variables$LOCS, p$variables$TIME) )
  p$variables$COV = setdiff( p$variables$ALL, p$variables$coordinates )  # non-location and non-time based covariates

  return (p)
}
