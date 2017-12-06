
stm__twostep = function( p, dat, pa, nu=NULL, phi=NULL, varObs=varObs, varSpatial=varSpatial, variablelist=FALSE, ... ) {

  #\\ twostep modelling time first as a simple ts and then spatial or spatio-temporal interpolation
  #\\ nu is the bessel smooth param

  # step 1 -- timeseries modelling
  # use all available data in 'dat' to get a time trend .. and assume it applies to the prediction area of interest 'pa' 
     # some methods require a uniform (temporal with associated covariates) prediction grid based upon all dat locations 
  if (variablelist)  return( c() )

  px = dat # only the static parts .. time has to be a uniform grid so reconstruct below

  ids = array_map( "xy->1", px[, c("plon", "plat")], gridparams=p$gridparams ) # 100X faster than paste / merge
  todrop = which(duplicated( ids) )
  if (length(todrop>0)) px = px[-todrop,]
  rm(ids, todrop)

  # static vars .. don't need to look up
  tokeep = c(p$variables$LOCS )
  if (exists("weights", dat) ) tokeep = c(tokeep, "weights")
  if (p$nloccov > 0) {
    for (ci in 1:p$nloccov) {
      vn = p$variables$local_cov[ci]
      pu = stm_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
      nts = ncol(pu)
      if ( nts==1 ) tokeep = c(tokeep, vn ) 
    }
  }
  px = px[ , tokeep ]
  px_n = nrow(px)
  nts = vn = NULL

  # add temporal grid
  if ( exists("TIME", p$variables) ) {
    px = cbind( px[ rep.int(1:px_n, p$nt), ], 
                    rep.int(p$prediction.ts, rep(px_n, p$nt )) )
    names(px)[ ncol(px) ] = p$variables$TIME 
    px = cbind( px, stm_timecovars ( vars=p$variables$local_all, ti=px[,p$variables$TIME]  ) )
  }

  if (p$nloccov > 0) {
    # add time-varying covars .. not necessary except when covars are modelled locally
    for (ci in 1:p$nloccov) {
      vn = p$variables$local_cov[ci]
      pu = stm_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
      nts = ncol(pu)
      if ( nts== 1) {
        # static vars are retained in the previous step
      } else if ( nts == p$ny )  {
        px$iy = px$yr - p$yrs[1] + 1 #yr index
        px[,vn] = pu[ cbind(px$i, px$iy) ]  
       } else if ( nts == p$nt) {
        px$it = p$nw*(px$tiyr - p$yrs[1] - p$tres/2) + 1 #ts index
        px[,vn] = pu[ cbind(px$i, px$it) ]  
      }
    } # end for loop
    nts = vn = NULL
  } # end if
  rownames(px) = NULL

  ts_gam = stm__gam( p, dat, px ) # currently only a GAM is enabled for the TS component

  if (is.null( ts_gam)) return(NULL)
  if (ts_gam$stm_stats$rsquared < p$stm_rsquared_threshold ) return(NULL)

  # range checks
  rY = range( dat[,p$variables$Y], na.rm=TRUE)
  toosmall = which( ts_gam$predictions$mean < rY[1] )
  toolarge = which( ts_gam$predictions$mean > rY[2] )
  if (length(toosmall) > 0) ts_gam$predictions$mean[toosmall] =  rY[1]  
  if (length(toolarge) > 0) ts_gam$predictions$mean[toolarge] =  rY[2]
 
  pxts = ts_gam$predictions
  rownames(pxts) = NULL
  ts_gam = NULL
  gc()
  
  names(pxts)[which(names(pxts)=="mean")] = p$variables$Y
  names(pxts)[which(names(pxts)=="sd")] = paste(p$variables$Y, "sd", sep=".")

  if(0){
      # debugging plots
      ti = 668
      xi = which( pxts[ , p$variables$TIME ] == p$prediction.ts[ti] )
      mbas = MBA::mba.surf( pxts[xi, c( p$variables$LOCS, p$variables$Y) ], 300, 300, extend=TRUE)$xyz.est
      image(mbas)
  }

  out = NULL

  # step 2 :: spatial modelling .. essentially a time-space separable solution

  if (!exists( "stm_twostep_space", p)) p$stm_twostep_space="krige" # default
  
  if ( p$stm_twostep_space == "krige" ) {
    out = stm__krige( p, dat=pxts, pa=pa, nu=nu, phi=phi, varObs=varObs, varSpatial=varSpatial ) 
    if (is.null( out)) return(NULL)
  }

  if ( p$stm_twostep_space == "gstat" ) {
    out = stm__gstat( p, dat=pxts, pa=pa, nu=nu, phi=phi, varObs=varObs, varSpatial=varSpatial ) 
    if (is.null( out)) return(NULL)
  }

  if (p$stm_twostep_space %in% c("tps") ) {
    out = stm__tps( p, dat=pxts, pa=pa, lambda=varObs/varSpatial  )  
    if (is.null( out)) return(NULL)
  }

  if (p$stm_twostep_space %in% c("fft", "lowpass", "spatial.process", "lowpass_spatial.process") ) {
    out = stm__fft( p, dat=pxts, pa=pa, nu=nu, phi=phi )  
    if (is.null( out)) return(NULL)
  }


  return( out )  
}

