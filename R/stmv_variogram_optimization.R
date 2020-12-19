
stmv_variogram_optimization = function( vg, vx, nu=NULL, plotvgm=FALSE, cor=0.1, control=list(factr=1e-9, maxit = 1000L), weight_by_inverse_distance=TRUE, distance_scaling_factor=1 ) {
  #\\ simple nonlinear least squares fit

  vgm_var_max = max(vg)
  vgm_dist_max = max(vx)

  vgm_internal_scale = vgm_dist_max / 3

  vxs = vx / vgm_internal_scale
  vgs = vg / vgm_var_max

  if (any(!is.finite(vxs)) | any(!is.finite(vgs))) {
    fit=list()
    fit$summary = list(phi_ok=FALSE )
    return(fit)
  }


  if (weight_by_inverse_distance) {
    w = 1/vxs
  } else {
    w = rep(1, length(vxs ) )
  }


  if (!is.null(nu)) {  # ie. nu is fixed

    vario_function_phi = function(par, vgs, vxs, nu, w){
      # ie. nu is fixed
      if (par["phi"] < 0.001) par["phi"] = 0.001
      if (par["phi"] > 5 ) par["phi"] = 5
      vgm = par["total.var"] *( (1 -par["sigma.sq.fraction"])  +  par["sigma.sq.fraction"]*( 1-stmv_matern(distance=vxs, mRange=par["phi"], mSmooth=nu ) ) )
      # vgm = par["tau.sq"] + par["sigma.sq"]*{ 1-stmv_matern(distance=vxs, mRange=par["phi"], mSmooth=nu) }
      obj = sum( w * (vgs - vgm)^2, na.rm=TRUE) # vario normal errors, no weights , etc.. just the line
      return(obj)
    }

    par = c(total.var=1, sigma.sq.fraction=0.75, phi=0.9 )
    lower =c(0.5, 0, 0.001 )
    upper =c(1.5, 1, 5)

    fit = try( optim( par=par, vgs=vgs, vxs=vxs, nu=nu, w=w, method="Nelder-Mead", fn=vario_function_phi ) )

    if ( !inherits(fit, "try-error")) if ( fit$convergence != 0 ) class(fit) = "try-error"
    if (exists( "par", fit)) {
      if ( fit$par[["phi"]] < lower[3] | fit$par[["phi"]] > upper[3] ) class(fit) = "try-error"    # give up
    }


    if ( inherits(fit, "try-error") ) {
      par = c(total.var=1, sigma.sq.fraction=.25, phi=0.5)
      ilast = c(length(vgs)-1, length(vgs))  # last two data point .. drop it
      fit = try( optim( par=par, vvg=vgs[-ilast], vxs=vxs[-ilast], nu=nu, w=w[-ilast], method="L-BFGS-B", lower=lower, upper=upper, fn=vario_function_phi,
        control=list(factr=1e-12, maxit = 5000L) ))
      if ( !inherits(fit, "try-error"))  if ( fit$convergence != 0 ) class(fit) = "try-error"
      if (exists( "par", fit)) {
        if ( fit$par[["phi"]] < lower[3] | fit$par[["phi"]] > upper[3] ) class(fit) = "try-error"    # give up
      }
    }

  } else {


    vario_function_phi_nu = function(par, vgs, vxs, w){
      # ie. nu and phi are both estimated
      if (par["nu"] < 0.01) par["nu"] = 0.01
      if (par["phi"] <= 0.001 ) par["phi"] = 0.001
      if (par["nu"] > 5 ) par["nu"] = 5
      if (par["phi"] > 5 ) par["phi"] = 5

      vgm = par["total.var"] *( (1 -par["sigma.sq.fraction"])  +  par["sigma.sq.fraction"]*( 1-stmv_matern(distance=vxs, mRange=par["phi"], mSmooth=par["nu"]) ) )
      obj = sum( w * (vgs - vgm)^2, na.rm=TRUE) # vario normal errors, no weights , etc.. just the line
      return(obj)
    }


    par = c(total.var=1, sigma.sq.fraction=0.5, phi=1.1, nu=0.9)
    lower =c(0.5, 0, 0.001, 0.01 )
    upper =c(1.5, 1, 5,     5)

    fit = try( optim( par=par, vgs=vgs, vxs=vxs, w=w, method="L-BFGS-B", lower=lower, upper=upper, fn=vario_function_phi_nu, control=control ))
      if ( !inherits(fit, "try-error")) if ( fit$convergence != 0 ) class(fit) = "try-error"
      if (exists( "par", fit)) {
        if ( fit$par[["phi"]] < lower[3] | fit$par[["phi"]] > upper[3] ) class(fit) = "try-error"    # give up
        if ( fit$par[["nu"]]  < lower[4] | fit$par[["nu"]]  > upper[4] ) class(fit) = "try-error"    # give up
      }

    if ( inherits(fit, "try-error") ) {
      par = c(total.var=1, sigma.sq.fraction=0.25, phi=0.25, nu=0.25)
      ilast = c(length(vgs)-1, length(vgs))  # last two data point .. drop it
      fit = try( optim( par=par, vgs=vgs[-ilast], vxs=vxs[-ilast], w=w[-ilast], method="L-BFGS-B",  lower=lower, upper=upper, fn=vario_function_phi_nu,
        control=list(factr=1e-12, maxit = 5000L)  ))
      # fit = try( optim( par=par, vgs=vgs[-ilast], vxs=vxs[-ilast], w=w[-ilast], method="L-BFGS-B", fn=vario_function_phi_nu))
      if ( !inherits(fit, "try-error"))  if ( fit$convergence != 0 ) class(fit) = "try-error"
      if (exists( "par", fit)) {
        if ( fit$par[["phi"]] < lower[3] | fit$par[["phi"]] > upper[3] ) class(fit) = "try-error"    # give up
        if ( fit$par[["nu"]]  < lower[4] | fit$par[["nu"]]  > upper[4] ) class(fit) = "try-error"    # give up
      }
    }

    if ( inherits(fit, "try-error") ) {
      par = c(total.var=1, sigma.sq.fraction=0.5, phi=1, nu=0.5)
      ilast = c(length(vgs)-1, length(vgs))  # last two data point .. drop it
      fit = try( optim( par=par, vgs=vgs[-ilast], vxs=vxs[-ilast], w=w[-ilast], method="L-BFGS-B",  lower=lower, upper=upper, fn=vario_function_phi_nu,
        control=list(factr=1e-9, maxit = 5000L)  ))
      # fit = try( optim( par=par, vgs=vgs[-ilast], vxs=vxs[-ilast], w=w[-ilast], method="L-BFGS-B", fn=vario_function_phi_nu))
      if ( !inherits(fit, "try-error"))  if ( fit$convergence != 0 ) class(fit) = "try-error"
      if (exists( "par", fit)) {
        if ( fit$par[["phi"]] < lower[3] | fit$par[["phi"]] > upper[3] ) class(fit) = "try-error"    # give up
        if ( fit$par[["nu"]]  < lower[4] | fit$par[["nu"]]  > upper[4] ) class(fit) = "try-error"    # give up
      }
    }

  }

  fit$summary = list(
    vg = vg,
    vx = vx * distance_scaling_factor,
    vgm_var_max= vgm_var_max,
    vgm_dist_max = vgm_dist_max * distance_scaling_factor ,
    autocorrelation_function="matern",
    nu =NA,
    phi=NA,
    localrange = NA,
    varSpatial=NA,
    varObs=NA ,
    phi_ok = FALSE,
    objfn = NA
  )

  if ( !inherits(fit, "try-error")) {
   if ( fit$convergence==0 ) {
      fit$summary$nu = ifelse( !is.null(nu), nu, fit$par[["nu"]] )
      fit$summary$phi=fit$par[["phi"]] * vgm_internal_scale * distance_scaling_factor
      fit$summary$localrange =matern_phi2distance( phi=fit$summary$phi, nu=fit$summary$nu, cor=cor )
      fit$summary$varSpatial = fit$par[["total.var"]] * fit$par[["sigma.sq.fraction"]]  * vgm_var_max
      fit$summary$varObs = fit$par[["total.var"]] * (1- fit$par[["sigma.sq.fraction"]]) * vgm_var_max
      fit$summary$phi_ok = ifelse( fit$summary$phi < fit$summary$vgm_dist_max*0.99 , TRUE, FALSE )
      fit$summary$objfn=fit$value
   }
  }


  # message( " Optim flag (0==all good): ", fit$convergence, " --- ", fit$message )
  if (!is.finite(  fit$summary$phi *  fit$summary$nu  )) return(fit)

  if( plotvgm ) {
    dev.new()
    xlim= c(0, fit$summary$vgm_dist_max*1.1 )
    ylim= c(0, vgm_var_max*1.1)
    plot( fit$summary$vx, fit$summary$vg, col="green", xlim=xlim, ylim=ylim )
    ds = seq( 0, fit$summary$vgm_dist_max  , length.out=100 )

    ac = fit$summary$varObs + fit$summary$varSpatial*(1 - stmv_matern( ds, fit$summary$phi, fit$summary$nu ) )
    lines( ds, ac, col="orange" )
    abline( h=0, lwd=1, col="lightgrey" )
    abline( v=0 ,lwd=1, col="lightgrey" )
    abline( h=fit$summary$varObs, lty="dashed", col="grey" )
    abline( h=fit$summary$varObs + fit$summary$varSpatial, lty="dashed", col="grey" )
    abline( v=fit$summary$localrange, lty="dashed", col="grey")
  }

  return(fit)
}
