    stmv_variogram_fft_original = function( xyz, nx=NULL, ny=NULL, nbreaks=30, plotdata=FALSE, eps=1e-9, add.interpolation=FALSE,
      stmv_localrange_correlation=0.1, stmv_fft_taper_factor=5 ) {

    ## This is a copy that is easier to follow .. but uses more RAM

      if (0) {
        XYZ = as.data.frame( cbind( RMprecip$x, RMprecip$y, RMprecip$elev ) )
        names(XYZ) =c( "x","y","z", "elev" )
        xyz = XYZ[, c("x", "y", "z")]
        nbreaks=30
        plotdata=FALSE
        eps=1e-9
        add.interpolation=FALSE
        stmv_localrange_correlation=0.1
        stmv_fft_taper_factor=5  # multiplier
      }


      names(xyz) =c("x", "y", "z")
      zmean = mean(xyz$z, na.rm=TRUE)
      zsd = sd(xyz$z, na.rm=TRUE)
      zvar = zsd^2
      Z = (xyz$z - zmean) / zsd # zscore -- making it mean 0 removes the DC component

      if (is.null(nx)) {
        nx = ny = floor( nbreaks * 2.35 )
      }

      # system size
      nr = nx
      nc = ny

      rr = diff( range(xyz$x) )  # system length scale
      rc = diff( range(xyz$y) )

      # no of elements
      dr = rr/(nr-1)
      dc = rc/(nc-1)

      # approx sa associate with each datum
      sa = rr * rc
      d_sa = sa/nrow(xyz) # sa associated with each datum
      d_length = sqrt( d_sa/pi )  # sa = pi*l^2  # characteristic length scale

      nr2 = 2 * nr
      nc2 = 2 * nc

      mY = matrix(0, nrow = nr2, ncol = nc2)
      mN = matrix(0, nrow = nr2, ncol = nc2)


      if (0) {
        u = as.image( Z=Z, x=xyz[, c("x", "y")], na.rm=TRUE, nx=nr, ny=nc )
        # surface(u)
      }


      #  Nadaraya/Watson normalization for missing values s
      coo = as.matrix(array_map( "xy->2", coords=xyz[,c("x", "y")], origin=origin, res=resolution ))
      yy = tapply( X=Z, INDEX=list(coo[,1], coo[,2]),  FUN = function(w) {mean(w, na.rm=TRUE)}, simplify=TRUE )
      nn = tapply( X=Z, INDEX=list(coo[,1], coo[,2]), FUN = function(w) {length(w)}, simplify=TRUE )
      fY[as.numeric(dimnames(yy)[[1]]),as.numeric(dimnames(yy)[[2]])] = yy
      fN[as.numeric(dimnames(nn)[[1]]),as.numeric(dimnames(nn)[[2]])] = nn
      yy = nn = NULL
      coo = NULL
      mY[!is.finite(mY)] = 0
      mN[!is.finite(mN)] = 0

      # See explanation:  https://en.wikipedia.org/wiki/Autocorrelation#Efficient_computation
      # Robertson, C., & George, S. C. (2012). Theory and practical recommendations for autocorrelation-based image
      # correlation spectroscopy. Journal of biomedical optics, 17(8), 080801. doi:10.1117/1.JBO.17.8.080801
      # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3414238/

      fY = fftwtools::fftw2d(mY)
      fN = fftwtools::fftw2d(mN)

      mY = NULL
      mN = NULL

      # fY * Conj(fY) == power spectra
      ii = Re( fftwtools::fftw2d( fY * Conj(fY), inverse=TRUE)  ) # autocorrelation (amplitude)
      jj = Re( fftwtools::fftw2d( fN * Conj(fN), inverse=TRUE)  ) # autocorrelation (amplitude) correction
      X = ifelse(( jj > eps), (ii / jj), NA) # autocorrelation (amplitude)
      ii = NULL
      jj = NULL

      # fftshift
      X = rbind( X[((nr+1):nr2), (1:nc2)], X[(1:nr), (1:nc2)] )  # swap_up_down
      X = cbind(X[1:nr2, ((nc+1):nc2)], X[1:nr2, 1:nc])  # swap_left_right

      # radial representation
      xy = expand_grid_fast( c(-(nr-1):0, 0:(nr-1)) * dr,  c(-(nc-1):0, 0:(nc-1)) * dc )
      distances = sqrt(xy[,1]^2 + xy[,2]^2)
      dmax = max(distances, na.rm=TRUE ) * 0.4  # approx nyquist distance (<0.5 as corners exist)
      breaks = seq( 0, dmax, length.out=nbreaks)
      db = breaks[2] - breaks[1]
      # angles = atan2( xy[,2], xy[,1] )  # not used

      zz = cut( distances, breaks=c(breaks, max(breaks)+db), label=breaks+(db/2) )
      distances = NULL
      xy = NULL

      vgm = as.data.frame.table(tapply( X=X, INDEX=zz, FUN=mean, na.rm=TRUE ))
      names(vgm) = c("distances", "ac")
      X = NULL
      zz = NULL
      gc()

      vgm$distances = as.numeric( as.character(vgm$distances))
      vgm$sv =  zvar * (1-vgm$ac^2) # each sv are truly orthogonal

      # plot(ac ~ distances, data=vgm   )
      # plot(sv ~ distances, data=vgm   )

      out = list(vgm=vgm )

      if (add.interpolation) {
        # interpolated surface
      # constainer for spatial filters
        uu = which( (vgm$distances < dmax ) & is.finite(vgm$sv) )  # dmax ~ Nyquist freq
        fit = try( stmv_variogram_optimization( vx=vgm$distances[uu], vg=vgm$sv[uu], plotvgm=plotdata,
          stmv_internal_scale=dmax*0.75, cor=stmv_localrange_correlation ))
        out$fit = fit

        if (any(!is.finite( c(fit$summary$phi, fit$summary$nu) ))) return(out)

        phi = fit$summary$phi
        nu = fit$summary$nu

        grid.list = list((1:nr2) * dr, (1:nc2) * dc)
        # dgrid = as.matrix(expand.grid(grid.list))  # a bit slower
        dgrid = expand_grid_fast(  grid.list[[1]],  grid.list[[2]] )
        dimnames(dgrid) = list(NULL, names(grid.list))
        attr(dgrid, "grid.list") = grid.list
        grid.list = NULL

        mC = matrix(0, nrow = nr2, ncol = nc2)
        mC[nr, nc] = 1  # equal weights
        center = matrix(c((dr * nr), (dc * nc)), nrow = 1, ncol = 2)
        # theta.Taper = matern_phi2distance( phi=phi, nu=nu, cor=stmv_fft_taper_factor )
        theta.Taper = d_length * stmv_fft_taper_factor
        sp.covar =  stationary.taper.cov( x1=dgrid, x2=center, Covariance="Matern", theta=phi, smoothness=nu,
          Taper="Wendland", Taper.args=list(theta=theta.Taper, k=2, dimension=2), spam.format=TRUE)
        sp.covar = as.surface(dgrid, c(sp.covar))$z / (nr2*nc2)
        sp.covar.kernel = fftwtools::fftw2d(sp.covar)  / fftwtools::fftw2d(mC)
        ffY = Re( fftwtools::fftw2d( sp.covar.kernel * fY, inverse = TRUE))[1:nr, 1:nc]
        ffN = Re( fftwtools::fftw2d( sp.covar.kernel * fN, inverse = TRUE))[1:nr, 1:nc]
        Z = ifelse((ffN > eps), (ffY/ffN), NA)
        Z[!is.finite(Z)] = NA
        Z = Z * zsd + zmean # revert to input scale
        if (plotdata) {
          dev.new()
          surface(list(x=c(1:nr)*dr, y=c(1:nc)*dc, z=Z), xaxs="r", yaxs="r")
        }
        out$Z = Z
      }

      return(out)
    }
