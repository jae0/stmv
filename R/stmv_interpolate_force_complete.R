
stmv_interpolate_force_complete = function( ip=NULL, p, qn = c(0.005, 0.995) ) {

  #// designed to be called from stmv
  #// for the sake of speed and parallelization, the kernel density method via fft is written out again .. it is taken from fields::smooth.2d
  #// the spatial interpolation is smoother than what is expected from a kriging covariance but faster

  if (exists( "libs", p)) RLibrary( p$libs )
  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns

  S = stmv_attach( p$storage.backend, p$ptr$S )

  # do this here as the Stats are from the most reliable estimates
  nu  = median( S[,which( p$statsvars=="nu"  )], na.rm=TRUE )
  phi = median( S[,which( p$statsvars=="phi" )], na.rm=TRUE )

  P = stmv_attach( p$storage.backend, p$ptr$P )
  Psd = stmv_attach( p$storage.backend, p$ptr$Psd )
  Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )

  P_qn = quantile(P[], probs=qn, na.rm = TRUE )
  Psd_qn = quantile(Psd[], probs=qn, na.rm = TRUE )

  dx = dy = p$pres
  nr = p$nplons
  nc = p$nplats


  if (p$stmv_force_complete_method=="fft") {

    tol = 1e-9
    x_r = range(Ploc[,1])
    x_c = range(Ploc[,2])

    # dr = ( diff(x_r) + diff(x_c) ) / 2  # system size in user units
    # nr = round( diff(x_r)/p$pres ) + 1
    # nc = round( diff(x_c)/p$pres ) + 1

    nr2 = 2 * nr
    nc2 = 2 * nc

    # spatial autocorrelation filter
    mC = matrix(0, nrow = nr2, ncol = nc2)
    mC[nr, nc] = 1  # center

    dgrid = make.surface.grid(list((1:nr2) * dx, (1:nc2) * dy))
    center = matrix(c((dx * nr2)/2, (dy * nc2)/2), nrow = 1, ncol = 2)

    theta.Taper = matern_phi2distance( phi=phi, nu=nu, cor=p$stmv_range_correlation_fft_taper )
    sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=phi, smoothness=nu,
      Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)
    sp.covar = as.surface(dgrid, c(sp.covar))$z
    sp.covar.kernel = fft(sp.covar) / {fft(mC) * nr2 * nc2}  # covar centered on middle of matrix

    dgrid = sp.covar = mC = NULL

    for ( iip in ip ) {

      ww = p$runs[ iip, "time_index" ]

      tofill = which( ! is.finite( P[,ww] ) )
      if (length( tofill) > 0 ) {

        rY = range( P[,ww], na.rm=TRUE )

        mY = mN = matrix(0, nrow = nr2, ncol = nc2 ) # density

        mY[1:nr,1:nc] = P[,ww] # fill with data in correct locations
        mY[!is.finite(mY)] = 0

        mN[1:nr,1:nc] = ifelse( is.finite( mY[1:nr,1:nc] ), 1, 0)  # locations with data
        mN[!is.finite(mN)] = 0

        # estimates based upon a global nu,phi .. they will fit to the immediate area near data and so retain their structure
        fN = Re(fft(fft(mN) * sp.covar.kernel, inverse = TRUE))[1:nr,1:nc]
        fY = Re(fft(fft(mY) * sp.covar.kernel, inverse = TRUE))[1:nr,1:nc]
        Z = ifelse((fN > tol), (fY/fN), NA)
        fY = fN = mN = mY = NULL

        iZ = which( !is.finite( Z))
        if (length(iZ) > 0) Z[iZ] = NA

        lb = which( Z < rY[1] )
        if (length(lb) > 0) Z[lb] = NA
        ub = which( Z > rY[2] )
        if (length(ub) > 0) Z[ub] = NA

        # image(Z)
        Z[ Z > P_qn[2] ]=NA
        Z[ Z < P_qn[1] ]=NA
        P[,ww][tofill] = Z[tofill]
      }


      ## SD estimates
      tofill = which( ! is.finite( Psd[,ww] ) )
      if (length( tofill) > 0 ) {

        rY = range( Psd[,ww], na.rm=TRUE )

        # counts
        mY = matrix(0, nrow = nr2, ncol = nc2 ) # density
        mY[1:nr,1:nc] = Psd[,ww] # fill with data in correct locations
        mN = ifelse( is.finite( mY ), 1, 0)  # locations with data
        mY[!is.finite(mY)] = 0  # this must come after the above

        # estimates based upon a global nu,phi .. they will fit to the immediate area near data and so retain their structure
        fN = Re(fft(fft(mN) * sp.covar.kernel, inverse = TRUE))[1:nr,1:nc]
        fY = Re(fft(fft(mY) * sp.covar.kernel, inverse = TRUE))[1:nr,1:nc]
        Z = ifelse((fN > tol), (fY/fN), NA)
        fY = fN = mN = mY = NULL

        iZ = which( !is.finite( Z))
        if (length(iZ) > 0) Z[iZ] = NA

        lb = which( Z < 0 )
        if (length(lb) > 0) Z[lb] = NA
        ub = which( Z > rY[2] )
        if (length(ub) > 0) Z[ub] = NA
        # image(Z)

        Z[ Z > Psd_qn[2] ]=NA
        Z[ Z < 0 ]=NA
        Psd[,ww][tofill] = Z[ tofill]
      }

    }

    return( "complete" )
  }

  # ------------------------

  if (p$stmv_force_complete_method=="mba") {  # cubic basis splines

    for ( iip in ip ) {

      ww = p$runs[ iip, "time_index" ]

      # means
      tofill = which( ! is.finite( P[,ww] ) )
      if (length( tofill) > 0 ) {
        # counts
        rY = range( P[,ww], na.rm=TRUE )

        Z = mba.surf(cbind(Ploc[,], P[,ww]) , no.X=nr, no.Y=nc, extend=TRUE)$xyz.est$z

        iZ = which( !is.finite( Z))
        if (length(iZ) > 0) Z[iZ] = NA

        lb = which( Z < rY[1] )
        if (length(lb) > 0) Z[lb] = NA
        ub = which( Z > rY[2] )
        if (length(ub) > 0) Z[ub] = NA

        # image(Z)
        Z[ Z > P_qn[2] ]=NA
        Z[ Z < P_qn[1] ]=NA
        P[,ww][tofill] = Z[tofill]
      }

      ## SD estimates
      tofill = which( ! is.finite( Psd[,ww] ) )
      if (length( tofill) > 0 ) {

        rY = range( Psd[,ww], na.rm=TRUE )

        Z = mba.surf( cbind(Ploc[], Psd[,ww]) , no.X=nr, no.Y=nc, extend=TRUE)$z

        iZ = which( !is.finite( Z))
        if (length(iZ) > 0) Z[iZ] = NA

        lb = which( Z < 0 )
        if (length(lb) > 0) Z[lb] = NA
        ub = which( Z > rY[2] )
        if (length(ub) > 0) Z[ub] = NA
        # image(Z)

        Z[ Z > Psd_qn[2] ]=NA
        Z[ Z < 0 ]=NA
        Psd[,ww][tofill] = Z[ tofill]

      }
    }
    return( "complete" )
  }

  # ----------------

  if (p$stmv_force_complete_method=="linear") {


    for ( iip in ip ) {

      ww = p$runs[ iip, "time_index" ]

      # means
      tofill = which( ! is.finite( P[,ww] ) )
      if (length( tofill) > 0 ) {
        # counts
        rY = range( P[,ww], na.rm=TRUE )

        u = as.image( P[,ww], x=Ploc[,], na.rm=TRUE, nx=nr, ny=nc )
        Z = as.vector( fields::interp.surface( u, loc=Ploc[,] ) ) # linear interpolation
        u = NULL

        iZ = which( !is.finite( Z))
        if (length(iZ) > 0) Z[iZ] = NA
        lb = which( Z < rY[1] )
        if (length(lb) > 0) Z[lb] = NA
        ub = which( Z > rY[2] )
        if (length(ub) > 0) Z[ub] = NA

        # image(Z)
        Z[ Z > P_qn[2] ]=NA
        Z[ Z < P_qn[1] ]=NA
        P[,ww][tofill] = Z[ tofill]
      }

      ## SD estimates
      tofill = which( ! is.finite( Psd[,ww] ) )
      if (length( tofill) > 0 ) {

        rY = range( Psd[,ww], na.rm=TRUE )

        u = as.image( Psd[,ww], x=Ploc[,], na.rm=TRUE, nx=nr, ny=nc )
        Z = as.vector( fields::interp.surface( u, loc=Ploc[,] ) ) # linear interpolation
        u = NULL

        iZ = which( !is.finite( Z))
        if (length(iZ) > 0) Z[iZ] = NA
        lb = which( Z < 0 )
        if (length(lb) > 0) Z[lb] = NA
        ub = which( Z > rY[2] )
        if (length(ub) > 0) Z[ub] = NA
        # image(Z)

        Z[ Z > Psd_qn[2] ]=NA
        Z[ Z < 0 ]=NA
        Psd[,ww][tofill] = Z[ tofill]

      }
    }
    return( "complete" )
  }

}
