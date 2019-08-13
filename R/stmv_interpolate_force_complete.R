
stmv_interpolate_force_complete = function( p, qn = c(0.005, 0.995), eps=1e-9 ) {

  #// designed to be called from stmv

  if (0) {
    # for debugging  runs ..
    eps=1e-9
    qn = c(0.005, 0.995)
    debugging=TRUE
  }

  if (exists( "libs", p)) RLibrary( p$libs )

  p = parallel_run( p=p, runindex=list( time_index=1:p$nt )  )
  ip = 1:p$nruns

  S = stmv_attach( p$storage.backend, p$ptr$S )

  P = stmv_attach( p$storage.backend, p$ptr$P )
  Psd = stmv_attach( p$storage.backend, p$ptr$Psd )
  Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )

  dx = dy = p$pres
  nr = p$nplons
  nc = p$nplats
 # system size
  #nr = nx
  #nc = ny

  # nr .. x/plon
  # nc .. y/plat

  x_r = range(Ploc[,1])
  x_c = range(Ploc[,2])

  origin=c(x_r[1], x_c[1])
  res=c(p$pres, p$pres)
  ind = as.matrix(array_map( "xy->2", coords=Ploc[], origin=origin, res=res ))


  if (p$stmv_force_complete_method=="kernel") {
    # essentially gaussian

    wght = setup.image.smooth( nrow=nr, ncol=nc,  dx=p$pres, dy=p$pres, theta=p$stmv_distance_statsgrid*3, xwidth=p$pres*10, ywidth=p$pres*10)
    # theta=p$stmv_distance_statsgrid*3 ensures that smoothing occurs. in a dense matrix, pres/3 is good to keep the texture when eg warping a full matrix.
    # but there may be many missing values and so 3*pres is safer though missing the texture
    # but this is ok as these areas repeated failed more refined methods

    for ( iip in ip ) {
      ww = p$runs[ iip, "time_index" ]
      # means
      tofill = which( ! is.finite( P[,ww] ) )
      if (length( tofill) > 0 ) {
        rY = range( P[,ww], na.rm=TRUE )
        X = matrix(NA, ncol=nc, nrow=nr )
        X[ind] = P[,ww]
        X = fields::image.smooth( X, dx=dx, dy=dy, wght )$z
        lb = which( X < rY[1] )
        if (length(lb) > 0) X[lb] = rY[1]
        lb = NULL
        ub = which( X > rY[2] )
        if (length(ub) > 0) X[ub] = rY[2]
        ub = NULL
        P[tofill,ww] = X[ tofill]
      }

      ## SD
      tofill = which( ! is.finite( Psd[,ww] ) )
      if (length( tofill) > 0 ) {
        rY = range( Psd[,ww], na.rm=TRUE )
        X = matrix(NA, ncol=nc, nrow=nr )
        X[ind] = Psd[,ww]
        X = fields::image.smooth( X, dx=dx, dy=dy, wght )$z
        lb = which( X < rY[1] )
        if (length(lb) > 0) X[lb] = rY[1]
        lb = NULL
        ub = which( X > rY[2] )
        if (length(ub) > 0) X[ub] = rY[2]
        ub = NULL
        Psd[tofill,ww] = X[ tofill]
      }
    }
    return( "complete" )
  }


  if (p$stmv_force_complete_method=="linear") {
    # fastest
    for ( iip in ip ) {
      ww = p$runs[ iip, "time_index" ]
      # means
      tofill = which( ! is.finite( P[,ww] ) )
      if (length( tofill) > 0 ) {
        rY = range( P[,ww], na.rm=TRUE )
        X = matrix(NA, ncol=nc, nrow=nr )
        X[ind] = P[,ww]
        X = fields::interp.surface( X, loc=Ploc[,] )$z # linear interpolation
        lb = which( X < rY[1] )
        if (length(lb) > 0) X[lb] = rY[1]
        lb = NULL
        ub = which( X > rY[2] )
        if (length(ub) > 0) X[ub] = rY[2]
        ub = NULL
        P[,ww][tofill] = X[ tofill]
      }

      ## SD
      tofill = which( ! is.finite( Psd[,ww] ) )
      if (length( tofill) > 0 ) {
        rY = range( Psd[,ww], na.rm=TRUE )
        X = matrix(NA, ncol=nc, nrow=nr )
        X[ind] = Psd[,ww]
        X = fields::interp.surface( X, loc=Ploc[,] )$z # linear interpolation
        # image(X)
        lb = which( X < rY[1] )
        if (length(lb) > 0) X[lb] = rY[1]
        lb = NULL
        ub = which( X > rY[2] )
        if (length(ub) > 0) X[ub] = rY[2]
        ub = NULL
        Psd[,ww][tofill] = X[ tofill]
      }
    }
    return( "complete" )
  }



  if (p$stmv_force_complete_method=="akima") {
    # basic cubic interpolation, but using all predicted fields rather than only data
    # slow .. esp on large problems
    require(akima)
    for ( iip in ip ) {
      ww = p$runs[ iip, "time_index" ]
      # means
      tofill = which( ! is.finite( P[,ww] ) )
      withdata = which( is.finite( P[,ww] ) )
      if (length( tofill) > 0 ) {
        rY = range( P[,ww], na.rm=TRUE )
        X = interp( x=Ploc[withdata,1], y=Ploc[withdata,2], z=P[withdata,ww],
          xo=p$plons, yo=p$plats, linear=FALSE, extrap=TRUE, duplicate="mean" )$z  # cannot extrapolate with linear, using cubic spline
        withdata = NULL
        lb = which( X < rY[1] )
        if (length(lb) > 0) X[lb] = rY[1]
        lb = NULL
        ub = which( X > rY[2] )
        if (length(ub) > 0) X[ub] = rY[2]
        ub = NULL
        P[,ww][tofill] = X[ tofill]
        X = NULL
        tofill = NULL
        gc()
      }

      ## SD estimates
      tofill = which( ! is.finite( Psd[,ww] ) )
      withdata = which( is.finite( Psd[,ww] ) )
      if (length( tofill) > 0 ) {
        rY = range( Psd[,ww], na.rm=TRUE )
        X = interp( x=Ploc[withdata,1], y=Ploc[withdata,2], z=Psd[withdata,ww],
          xo=p$plons, yo=p$plats, linear=FALSE, extrap=TRUE, duplicate="mean" )$z  # cannot extrapolate with linear, using cubic spline
        withdata = NULL
        lb = which( X < rY[1] )
        if (length(lb) > 0) X[lb] = rY[1]
        lb = NULL
        ub = which( X > rY[2] )
        if (length(ub) > 0) X[ub] = rY[2]
        ub = NULL
        Psd[,ww][tofill] = X[ tofill]
        X = NULL
        tofill = NULL
        gc()
      }
    }
    return( "complete" )
  }



  if (p$stmv_force_complete_method=="fft") {
    #// for the sake of speed and parallelization, the kernel density method via fft is written out again .. it is taken from fields::smooth.2d
    #// the spatial interpolation is smoother than what is expected from a kriging covariance but faster

    tol = 1e-9
    P_qn = quantile(P[], probs=qn, na.rm = TRUE )
    Psd_qn = quantile(Psd[], probs=qn, na.rm = TRUE )

    x_r = range(Ploc[,1])
    x_c = range(Ploc[,2])

    dr = dx
    dc = dy

    # dr = ( diff(x_r) + diff(x_c) ) / 2  # system size in user units
    # nr = round( diff(x_r)/p$pres ) + 1
    # nc = round( diff(x_c)/p$pres ) + 1

    nr2 = 2 * nr
    nc2 = 2 * nc

    range_95 = median( S[,which( p$statsvars=="localrange" )], na.rm=TRUE  )
    nu  = 0.5 ; #    --- exponential

    phi = matern_distance2phi( distance=range_95, nu=nu  )

  # constainer for spatial filters
    grid.list = list((1:nr2) * dx, (1:nc2) * dy)
    # dgrid = as.matrix(expand.grid(grid.list))  # a bit slower
    dgrid = expand_grid_fast(  grid.list[[1]],  grid.list[[2]] )
    dimnames(dgrid) = list(NULL, names(grid.list))
    attr(dgrid, "grid.list") = grid.list
    grid.list = NULL

    center = matrix(c((dx * nr), (dy * nc)), nrow = 1, ncol = 2)

    # spatial autocorrelation filter
    mC = matrix(0, nrow = nr2, ncol = nc2)
    mC[nr, nc] = 1
    mC_fft = 1 / fftwtools::fftw2d(mC)  # multiplication is faster than division
    mC = NULL

    vgm = NA

    origin = c(x_r[1], x_c[1])
    resolution = c(dx, dy)

    # precompute a few things that are used repeatedly for time-specific variograms
    xy = expand_grid_fast( c(-(nr-1):0, 0:(nr-1)) * dr,  c(-(nc-1):0, 0:(nc-1)) * dc )
    distances = sqrt(xy[,1]^2 + xy[,2]^2)
    dmax = max(distances, na.rm=TRUE ) * 0.4  # approx nyquist distance (<0.5 as corners exist)
    breaks = seq( 0, dmax, length.out=nr)
    db = breaks[2] - breaks[1]
    # angles = atan2( xy[,2], xy[,1])  # not used
    zz = cut( distances, breaks=c(breaks, max(breaks)+db), label=breaks+(db/2) )
    xy = NULL
    distances = NULL
    breaks = NULL

    fY = matrix(0, nrow = nr2, ncol = nc2)
    fN = matrix(0, nrow = nr2, ncol = nc2)

    sp.covar = stationary.cov( dgrid, center, Covariance="Matern", theta=phi, smoothness=nu )
    sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
    sp.covar = fftwtools::fftw2d(sp.covar) * mC_fft

    for ( iip in ip ) {

      ww = p$runs[ iip, "time_index" ]

      tofill = which( ! is.finite( P[,ww] ) )
      if (length( tofill) > 0 ) {

        rY = range( P[,ww], na.rm=TRUE )
        fY = fN = matrix(0, nrow = nr2, ncol = nc2 ) # density

    # bounds check: make sure predictions exist
        zmean = mean(P[,ww], na.rm=TRUE)
        zsd = sd(P[,ww], na.rm=TRUE)
        zvar = zsd^2

        fY[1:nr,1:nc] = (P[,ww] - zmean) / zsd # zscore -- making it mean 0 removes the DC component,  # fill with data in correct locations
        fY[!is.finite(fY)] = 0

        fN[1:nr,1:nc] = ifelse( is.finite( fY[1:nr,1:nc] ), 1, 0)  # locations with data
        fN[!is.finite(fN)] = 0

        fY = fftwtools::fftw2d(fY)
        fN = fftwtools::fftw2d(fN)

        fY = Re( fftwtools::fftw2d( sp.covar * fY, inverse = TRUE))[1:nr, 1:nc]
        fN = Re( fftwtools::fftw2d( sp.covar * fN, inverse = TRUE))[1:nr, 1:nc]

        X = ifelse((fN > eps), (fY/fN), NA)
        X[!is.finite(X)] = NA
        X = X * zsd + zmean # revert to input scale
        if (0) {
          dev.new()
          surface(list(x=c(1:nr)*dr, y=c(1:nc)*dc, z=X), xaxs="r", yaxs="r")
        }

        iX = which( !is.finite( X))
        if (length(iX) > 0) X[iX] = NA

        lb = which( X < rY[1] )
        if (length(lb) > 0) X[lb] = rY[1]
        ub = which( X > rY[2] )
        if (length(ub) > 0) X[ub] = rY[2]

        # image(X)
        X[ X > P_qn[2] ]=P_qn[2]
        X[ X < P_qn[1] ]=P_qn[1]
        P[,ww][tofill] = X[tofill]
      }


      ## SD estimates
      tofill = which( ! is.finite( Psd[,ww] ) )
      if (length( tofill) > 0 ) {

        fY = fN = matrix(0, nrow = nr2, ncol = nc2 ) # density

    # bounds check: make sure predictions exist
        rY = range( Psd[,ww], na.rm=TRUE )
        zmean = mean(Psd[,ww], na.rm=TRUE)
        zsd = sd(Psd[,ww], na.rm=TRUE)
        zvar = zsd^2

        fY[1:nr,1:nc] = (Psd[,ww] - zmean) / zsd # zscore -- making it mean 0 removes the DC component, fill with data in correct locations
        fY[!is.finite(fY)] = 0

        fN[1:nr,1:nc] = ifelse( is.finite( fY[1:nr,1:nc] ), 1, 0)  # locations with data
        fN[!is.finite(fN)] = 0

        fY = fftwtools::fftw2d(fY)
        fN = fftwtools::fftw2d(fN)

        fY = Re( fftwtools::fftw2d( sp.covar * fY, inverse = TRUE))[1:nr, 1:nc]
        fN = Re( fftwtools::fftw2d( sp.covar * fN, inverse = TRUE))[1:nr, 1:nc]

        X = ifelse((fN > eps), (fY/fN), NA)
        X[!is.finite(X)] = NA
        X = X * zsd + zmean # revert to input scale
        if (0) {
          dev.new()
          surface(list(x=c(1:nr)*dr, y=c(1:nc)*dc, z=X), xaxs="r", yaxs="r")
        }


        iX = which( !is.finite( X))
        if (length(iX) > 0) X[iX] = NA

        lb = which( X < 0 )
        if (length(lb) > 0) X[lb] = NA
        ub = which( X > rY[2] )
        if (length(ub) > 0) X[ub] = NA
        # image(X)

        X[ X > Psd_qn[2] ]=NA
        X[ X < 0 ]=NA
        Psd[,ww][tofill] = X[ tofill]
      }

    }

    return( "complete" )
  }

  # ------------------------

  if (p$stmv_force_complete_method=="mba") {  # cubic basis splines

    P_qn = quantile(P[], probs=qn, na.rm = TRUE )
    Psd_qn = quantile(Psd[], probs=qn, na.rm = TRUE )

    for ( iip in ip ) {
      ww = p$runs[ iip, "time_index" ]
      # means
      tofill = which( ! is.finite( P[,ww] ) )
      if (length( tofill) > 0 ) {
        # counts
        rY = range( P[,ww], na.rm=TRUE )

        X = mba.surf(cbind(Ploc[,], P[,ww]) , no.X=nr, no.Y=nc, extend=TRUE)$xyz.est$z

        iX = which( !is.finite( X))
        if (length(iX) > 0) X[iX] = NA

        lb = which( X < rY[1] )
        if (length(lb) > 0) X[lb] = NA
        ub = which( X > rY[2] )
        if (length(ub) > 0) X[ub] = NA

        # image(X)
        X[ X > P_qn[2] ]=NA
        X[ X < P_qn[1] ]=NA
        P[,ww][tofill] = X[tofill]
      }

      ## SD estimates
      tofill = which( ! is.finite( Psd[,ww] ) )
      if (length( tofill) > 0 ) {
        rY = range( Psd[,ww], na.rm=TRUE )
        X = mba.surf( cbind(Ploc[], Psd[,ww]) , no.X=nr, no.Y=nc, extend=TRUE)$z
        iX = which( !is.finite( X))
        if (length(iX) > 0) X[iX] = NA

        lb = which( X < 0 )
        if (length(lb) > 0) X[lb] = NA
        ub = which( X > rY[2] )
        if (length(ub) > 0) X[ub] = NA
        # image(X)

        X[ X > Psd_qn[2] ]=NA
        X[ X < 0 ]=NA
        Psd[,ww][tofill] = X[ tofill]

      }
    }
    return( "complete" )
  }

  # ----------------

}
