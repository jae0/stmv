stmv_variablelist = function( p ) {
  
  if (!exists("variables", p)) p$variables = list()
  if (!exists("Y", p$variables)) {
    if (exists("stmv_local_modelformula", p))  p$variables$Y = all.vars( p$stmv_local_modelformula[[2]] ) 
    if (exists("stmv_global_modelformula", p)) p$variables$Y = all.vars( p$stmv_global_modelformula[[2]] )
  }
  if (!exists("Y", p$variables)) p$variables$Y = "not_defined"

  if (!exists("LOCS", p$variables)) p$variables$LOCS = c("plon", "plat")

  if (grepl("space-year", p$stmv_dimensionality)) {
    if (!exists("TIME", p$variables)) {
      p$variables$TIME = "tiyr" 
    }
  }

  if (!exists("COV", p$variables)) {
    p$variables$local_all = NULL  
    p$variables$local_cov = NULL
    if (exists("stmv_local_modelformula", p)) {
      p$variables$local_all = all.vars( p$stmv_local_modelformula )
      p$variables$local_cov = all.vars( p$stmv_local_modelformula[[3]] )
    }
    p$nloccov = 0
    if (exists("local_cov", p$variables)) p$nloccov = length(p$variables$local_cov)
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
  }
  return (p)
}