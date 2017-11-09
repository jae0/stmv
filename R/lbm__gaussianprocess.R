
ecmei__gaussianprocess = function( p, dat, pa ) {
  #\\ this is the core engine of ecmei .. localised space  and time modelling/ interpolation 
  # \ as a gaussian process
  # TODO 

  sdTotal=sd(dat[,p$variable$Y], na.rm=T)

  # plot(pred ~ z , dat)
  ss = lm( dat$mean ~ dat[,p$variables$Y], na.action=na.omit)
  if ( "try-error" %in% class( ss ) ) return( NULL )
  rsquared = summary(ss)$r.squared
  if (rsquared < p$ecmei_rsquared_threshold ) return(NULL)

  ecmei_stats = list( sdTotal=sdTotal, rsquared=rsquared, ndata=nrow(dat) ) # must be same order as p$statsvars
  
 
  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )

  return( list( predictions=pa, ecmei_stats=ecmei_stats ) )  
}

