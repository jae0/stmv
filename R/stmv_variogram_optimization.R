
stmv_variogram_optimization = function( vg, vx, nu=NULL, plotvgm=FALSE, stmv_internal_scale=NA, cor=0.1, control=list(factr=1e-9, maxit = 400L), weight_by_inverse_distance=TRUE ) {
  #\\ simple nonlinear least squares fit

  if (is.na(stmv_internal_scale)) stmv_internal_scale= max(vx)/2
  vx = vx / stmv_internal_scale

  if (any(!is.finite(vx)) | any(!is.finite(vg))) {
    fit=list()
    fit$summary = list(phi_ok=FALSE )
    return(fit)
  }

  vgm_var_max = max(vg)
  vgm_dist_max = max(vx)
  if (weight_by_inverse_distance) {
    w = 1/vx
  } else {
    w = rep(1, length(vx ) )
  }

  if (!is.null(nu)) {  # ie. nu is fixed
    vario_function = function(par, vg, vx, nu, w){
      if (par["phi"] < 0.001) return(Inf)
      vgm = par["tau.sq"] + par["sigma.sq"]*{ 1-stmv_matern(distance=vx, mRange=par["phi"], mSmooth=nu) }
      obj = sum( w * (vg - vgm)^2, na.rm=TRUE) # vario normal errors, no weights , etc.. just the line
      return(obj)
    }
    par = c(tau.sq=vgm_var_max*0.25, sigma.sq=vgm_var_max*0.75, phi=0.9 )
    lower =c(0, 0, 0.01 )
    upper =c(vgm_var_max*2, vgm_var_max*2, 3)

    fit = try( optim( par=par, vg=vg, vx=vx, nu=nu, w=w, method="Nelder-Meads", fn=vario_function ) )

    if ( !inherits(fit, "try-error")) if ( fit$convergence != 0 ) class(fit) = "try-error"
    if (exists( "par", fit)) {
      if ( fit$par[["phi"]] < lower[3] | fit$par[["phi"]] > upper[3] ) class(fit) = "try-error"    # give up
    }

    if ( inherits(fit, "try-error") ) {
      ilast = length(vg)  # last data point .. drop it
      # fit = try( optim( par=par, vg=vg[-ilast], vx=vx[-ilast], nu=nu, w=w[-ilast], method="L-BFGS-B", lower=lower, upper=upper, fn=vario_function,
      #   control=list(factr=1e-9, maxit = 500L) ))
      fit = try( optim( par=par, vg=vg[-ilast], vx=vx[-ilast], nu=nu, w=w[-ilast], method="BFGS", fn=vario_function) )
      if ( !inherits(fit, "try-error"))  if ( fit$convergence != 0 ) class(fit) = "try-error"
      if (exists( "par", fit)) {
        if ( fit$par[["phi"]] < lower[3] | fit$par[["phi"]] > upper[3] ) class(fit) = "try-error"    # give up
      }

    }

    if ( inherits(fit, "try-error") ) {
      par = c(tau.sq=vgm_var_max*0.25, sigma.sq=vgm_var_max*0.25, phi=0.5)
      ilast = c(length(vg)-1, length(vg))  # last two data point .. drop it
      fit = try( optim( par=par, vvg=vg[-ilast], vx=vx[-ilast], nu=nu, w=w[-ilast], method="L-BFGS-B", lower=lower, upper=upper, fn=vario_function,
        control=list(factr=1e-11, maxit = 1000L) ))
      if ( !inherits(fit, "try-error"))  if ( fit$convergence != 0 ) class(fit) = "try-error"
      if (exists( "par", fit)) {
        if ( fit$par[["phi"]] < lower[3] | fit$par[["phi"]] > upper[3] ) class(fit) = "try-error"    # give up
      }
    }

  } else {

    vario_function = function(par, vg, vx, w){
      if (par["nu"] < 0.001) return(Inf)
      if (par["phi"] < 0 ) return(Inf)
      vgm = par["tau.sq"] + par["sigma.sq"]*{ 1-stmv_matern(distance=vx, mRange=par["phi"], mSmooth=par["nu"]) }
      obj = sum( ((vg - vgm)^2) * w, na.rm=TRUE) # vario normal errors, no weights , etc.. just the line
      return(obj)
    }
    par = c(tau.sq=vgm_var_max*0.5, sigma.sq=vgm_var_max*0.5, phi=1.1, nu=0.9)
    lower =c(0, 0, 0.01, 0.01 )
    upper =c(vgm_var_max*2, vgm_var_max*2, 3, 5)

    # fit = try( optim( par=fit$par, vg=vg, vx=vx, w=w, method="L-BFGS-B", lower=lower, upper=upper, fn=vario_function, control=control ))
    fit = try( optim( par=par, vg=vg, vx=vx, w=w, method="Nelder-Mead", fn=vario_function))
    if ( !inherits(fit, "try-error")) if ( fit$convergence != 0 ) class(fit) = "try-error"
    if (exists( "par", fit)) {
      if ( fit$par[["phi"]] < lower[3] | fit$par[["phi"]] > upper[3] ) class(fit) = "try-error"    # give up
      if ( fit$par[["nu"]]  < lower[4] | fit$par[["nu"]]  > upper[4] ) class(fit) = "try-error"    # give up
    }

    if ( inherits(fit, "try-error")) {
      ilast = length(vg)  # last data point .. drop it
      # fit = try( optim( par=par, vg=vg[-ilast], vx=vx[-ilast], w=w[-ilast], method="L-BFGS-B",  lower=lower, upper=upper, fn=vario_function,
      #   control=list(factr=1e-9, maxit = 500L)  ) )
      fit = try( optim( par=par, vg=vg[-ilast], vx=vx[-ilast], w=w[-ilast], method="BFGS", fn=vario_function))
      if ( !inherits(fit, "try-error"))  if ( fit$convergence != 0 ) class(fit) = "try-error"
      if (exists( "par", fit)) {
        if ( fit$par[["phi"]] < lower[3] | fit$par[["phi"]] > upper[3] ) class(fit) = "try-error"    # give up
        if ( fit$par[["nu"]]  < lower[4] | fit$par[["nu"]]  > upper[4] ) class(fit) = "try-error"    # give up
      }
    }


    if ( inherits(fit, "try-error") ) {
      par = c(tau.sq=vgm_var_max*0.25, sigma.sq=vgm_var_max*0.25, phi=0.25, nu=0.25)
      ilast = c(length(vg)-1, length(vg))  # last two data point .. drop it
      fit = try( optim( par=par, vg=vg[-ilast], vx=vx[-ilast], w=w[-ilast], method="L-BFGS-B",  lower=lower, upper=upper, fn=vario_function,
        control=list(factr=1e-11, maxit = 1000L)  ))
      # fit = try( optim( par=par, vg=vg[-ilast], vx=vx[-ilast], w=w[-ilast], method="L-BFGS-B", fn=vario_function))
      if ( !inherits(fit, "try-error"))  if ( fit$convergence != 0 ) class(fit) = "try-error"
      if (exists( "par", fit)) {
        if ( fit$par[["phi"]] < lower[3] | fit$par[["phi"]] > upper[3] ) class(fit) = "try-error"    # give up
        if ( fit$par[["nu"]]  < lower[4] | fit$par[["nu"]]  > upper[4] ) class(fit) = "try-error"    # give up
      }
    }

    if ( inherits(fit, "try-error") ) {
      par = c(tau.sq=vgm_var_max*0.5, sigma.sq=vgm_var_max*0.5, phi=1, nu=0.5)
      ilast = c(length(vg)-1, length(vg))  # last two data point .. drop it
      fit = try( optim( par=par, vg=vg[-ilast], vx=vx[-ilast], w=w[-ilast], method="L-BFGS-B",  lower=lower, upper=upper, fn=vario_function,
        control=list(factr=1e-9, maxit = 1000L)  ))
      # fit = try( optim( par=par, vg=vg[-ilast], vx=vx[-ilast], w=w[-ilast], method="L-BFGS-B", fn=vario_function))
      if ( !inherits(fit, "try-error"))  if ( fit$convergence != 0 ) class(fit) = "try-error"
      if (exists( "par", fit)) {
        if ( fit$par[["phi"]] < lower[3] | fit$par[["phi"]] > upper[3] ) class(fit) = "try-error"    # give up
        if ( fit$par[["nu"]]  < lower[4] | fit$par[["nu"]]  > upper[4] ) class(fit) = "try-error"    # give up
      }
    }

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
    phi_ok = FALSE,
    objfn = NA
  )

  if ( !inherits(fit, "try-error")) {
   if ( fit$convergence==0 ) {
      fit$summary$nu = ifelse( !is.null(nu), nu, fit$par[["nu"]] )
      fit$summary$phi=fit$par[["phi"]]
      fit$summary$varSpatial=fit$par[["sigma.sq"]]
      fit$summary$varObs=fit$par[["tau.sq"]]
      fit$summary$phi_ok = ifelse( fit$summary$phi < fit$summary$vgm_dist_max*0.99, TRUE, FALSE )
      fit$summary$objfn=fit$value
    }
  }

  fit$summary$phi = fit$summary$phi * stmv_internal_scale
  fit$summary$vx = fit$summary$vx * stmv_internal_scale
  fit$summary$vgm_dist_max = fit$summary$vgm_dist_max * stmv_internal_scale


  # message( " Optim flag (0==all good): ", fit$convergence, " --- ", fit$message )
  if (!is.finite(  fit$summary$phi *  fit$summary$nu  )) return(fit)

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
    localrange = matern_phi2distance( phi=fit$summary$phi, nu=fit$summary$nu, cor=cor )
    abline( v=localrange, lty="dashed", col="grey")
  }

  return(fit)
}
