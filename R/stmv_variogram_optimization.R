
stmv_variogram_optimization = function( vg, vx, nu=NULL, plotvgm=FALSE, eps=1e-6, stmv_internal_scale=NA, cor=0.1 ) {
  #\\ simple nonlinear least squares fit

  if (is.na(stmv_internal_scale)) stmv_internal_scale= max(vx)/2
  vx = vx / stmv_internal_scale

  if (any(!is.finite(vx)) | any(!is.finite(vg))) {
    fit=list()
    fit$summary = list(range_ok=FALSE )
    return(fit)
  }

  vgm_var_max = max(vg)
  vgm_dist_max = max(vx)

  if (!is.null(nu)) {  # ie. nu is fixed
    vario_function = function(par, vg, vx, nu){
      vgm = par["tau.sq"] + par["sigma.sq"]*{ 1-stmv_matern(distance=vx, mRange=par["phi"], mSmooth=nu) }
      obj = sum( {vg - vgm}^2, na.rm=TRUE) # vario normal errors, no weights , etc.. just the line
      return(obj)
    }
    par = c(tau.sq=vgm_var_max*0.05, sigma.sq=vgm_var_max*0.95, phi=1 )
    lower =c(0, 0, eps )
    upper =c(vgm_var_max*2, vgm_var_max*2, 5)

    fit = try( optim( par=par, vg=vg, vx=vx, nu=nu, method="L-BFGS-B", lower=lower, upper=upper, fn=vario_function ) )

  } else {

    vario_function = function(par, vg, vx){
      vgm = par["tau.sq"] + par["sigma.sq"]*{ 1-stmv_matern(distance=vx, mRange=par["phi"], mSmooth=par["nu"]) }
      obj = sum( {vg - vgm}^2, na.rm=TRUE) # vario normal errors, no weights , etc.. just the line
      return(obj)
    }
    par = c(tau.sq=vgm_var_max*0.2, sigma.sq=vgm_var_max*0.8, phi=1.1, nu=0.51)
    lower =c(0, 0, 0.5, 0.1 )
    upper =c(vgm_var_max*1.5, vgm_var_max*1.5, 4, 5)

    fit = try( optim( par=par, vg=vg, vx=vx, method="L-BFGS-B", lower=lower, upper=upper, fn=vario_function ) )
  }

  fit$summary = list(
    vg = vg,
    vx = vx,
    vgm_var_max= vgm_var_max,
    vgm_dist_max = vgm_dist_max,
    autocorrelation_function="matern",
    nu =NA,
    phi=NA,
    varSpatial=NA,
    varObs=NA ,
    range=NA,
    range_ok = FALSE,
    objfn = NA
  )

  if ( !inherits(fit, "try-error")) {
   if ( fit$convergence==0 ) {
      fit$summary$nu = ifelse( !is.null(nu), nu, fit$par[["nu"]] )
      fit$summary$phi=fit$par[["phi"]]
      fit$summary$varSpatial=fit$par[["sigma.sq"]]
      fit$summary$varObs=fit$par[["tau.sq"]]
      fit$summary$range=matern_phi2distance( phi=fit$summary$phi, nu=fit$summary$nu, cor=cor )
      fit$summary$range_ok = ifelse( fit$summary$range < fit$summary$vgm_dist_max*0.99, TRUE, FALSE )
      fit$summary$objfn=fit$value
    }
  }

  fit$summary$phi = fit$summary$phi * stmv_internal_scale
  fit$summary$range = fit$summary$range * stmv_internal_scale
  fit$summary$vx = fit$summary$vx * stmv_internal_scale
  fit$summary$vgm_dist_max = fit$summary$vgm_dist_max * stmv_internal_scale

  # message( " Optim flag (0==all good): ", fit$convergence, " --- ", fit$message )

  if( plotvgm ) {
    xlim= c(0, fit$summary$vgm_dist_max*1.1)
    ylim= c(0, vgm_var_max*1.1)
    plot( fit$summary$vx, fit$summary$vg, col="green", xlim=xlim, ylim=ylim )
    ds = seq( 0, fit$summary$vgm_dist_max, length.out=100 )
    ac = fit$summary$varObs + fit$summary$varSpatial*(1 - stmv_matern( ds, fit$summary$phi, fit$summary$nu ) )
    lines( ds, ac, col="orange" )
    abline( h=0, lwd=1, col="lightgrey" )
    abline( v=0 ,lwd=1, col="lightgrey" )
    abline( h=fit$summary$varObs, lty="dashed", col="grey" )
    abline( h=fit$summary$varObs + fit$summary$varSpatial, lty="dashed", col="grey" )
    abline( v=fit$summary$range, lty="dashed", col="grey")
  }

  return(fit)
}
