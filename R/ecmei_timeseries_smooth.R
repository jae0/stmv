
ecmei_timeseries_smooth = function(p, dat, sloc=sloc, distance=distance, datvarsout=c("id", p$variables$TIME, p$variables$LOCS, p$variables$Y)) {

  # static vars .. don't need to look up
    datgridded = dat # only the static parts .. time has to be a uniform grid so reconstruct below

    ids = array_map( "xy->1", datgridded[, c("plon", "plat")], gridparams=p$gridparams ) # 100X faster than paste / merge
    todrop = which(duplicated( ids) )
    if (length(todrop>0)) datgridded = datgridded[-todrop,]
    rm(ids, todrop)

    tokeep = c(p$variables$LOCS )
    if (exists("weights", dat) ) tokeep = c(tokeep, "weights")

    if (p$nloccov > 0) {
      for (ci in 1:p$nloccov) {
        vn = p$variables$local_cov[ci]
        pu = ecmei_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
        nts = ncol(pu)
        if ( nts==1 ) tokeep = c(tokeep, vn ) 
      }
    }

    datgridded = datgridded[ , tokeep ]
    datgridded_n = nrow(datgridded)
    nts = vn = NULL

    # add temporal grid
    if ( exists("TIME", p$variables) ) {
      datgridded = cbind( datgridded[ rep.int(1:datgridded_n, p$nt), ], 
                      rep.int(p$prediction.ts, rep(datgridded_n, p$nt )) )
      names(datgridded)[ ncol(datgridded) ] = p$variables$TIME 
      datgridded = cbind( datgridded, ecmei_timecovars ( vars=p$variables$local_all, ti=datgridded[,p$variables$TIME]  ) )
    }

    if (p$nloccov > 0) {
      # add time-varying covars .. not necessary except when covars are modelled locally
      for (ci in 1:p$nloccov) {
        vn = p$variables$local_cov[ci]
        pu = ecmei_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
        nts = ncol(pu)
        if ( nts== 1) {
          # static vars are retained in the previous step
        } else if ( nts == p$ny )  {
          datgridded$iy = datgridded$yr - p$yrs[1] + 1 #yr index
          datgridded[,vn] = pu[ cbind(datgridded$i, datgridded$iy) ]  
         } else if ( nts == p$nt) {
          datgridded$it = p$nw*(datgridded$tiyr - p$yrs[1] - p$tres/2) + 1 #ts index
          datgridded[,vn] = pu[ cbind(datgridded$i, datgridded$it) ]  
        }
      } # end for loop
      nts = vn = NULL
    } # end if

    rownames(datgridded) = NULL

    ts_gam = ecmei__gam( p, dat, datgridded ) # currently only a GAM is enabled for the TS component

    if (is.null( ts_gam)) return(NULL)
    if (ts_gam$ecmei_stats$rsquared < p$ecmei_rsquared_threshold ) return(NULL)

    # range checks
    rY = range( dat[,p$variables$Y], na.rm=TRUE)
    toosmall = which( ts_gam$predictions$mean < rY[1] )
    toolarge = which( ts_gam$predictions$mean > rY[2] )
    if (length(toosmall) > 0) ts_gam$predictions$mean[toosmall] =  NA 
    if (length(toolarge) > 0) ts_gam$predictions$mean[toolarge] =  NA
   
    # overwrite dat with interpolated predictions
    datgridded = ts_gam$predictions
    rownames(datgridded) = NULL
    ts_gam = NULL
    
    names(datgridded)[which(names(datgridded)=="mean")] = p$variables$Y
    names(datgridded)[which(names(datgridded)=="sd")] = paste(p$variables$Y, "sd", sep=".")

    windowsize.half = floor(distance/p$pres)
    pa_w = -windowsize.half : windowsize.half # default window size 
    # pa_w = pa_w[ -length(pa_w)] # must be even 
    pa_w_n = length(pa_w)
    adims = c(p$nt, pa_w_n, pa_w_n ) 
    datgridded$id = array_map( "3->1", 
      coords = round( cbind( 
        ( (datgridded[,p$variables$TIME ] - p$prediction.ts[1] ) / p$tres) + 1 ,
        ( windowsize.half + (datgridded[,p$variables$LOCS[1]] - sloc[1]) / p$pres) + 1, 
        ( windowsize.half + (datgridded[,p$variables$LOCS[2]] - sloc[2]) / p$pres) + 1)), 
      dims=adims )
    o = which( datgridded$id > prod(adims)  | datgridded$id <= 0 ) # remove the area outside the aoi
    if (length(o) > 0 ) datgridded = datgridded[-o,]
    
    ddup = which(duplicated(datgridded$id))
    
    if ( length( ddup) > 0 ) {
      dups = unique(datgridded$id[ddup])
      for ( i in dups ) {
        j = which( datgridded$id== i)
        meanvalue = mean( datgridded[j, p$variable$Y], na.rm=TRUE)  
        datgridded[j, p$variable$Y] = NA
        datgridded[j[1], p$variable$Y] = meanvalue
      }
      datgridded = datgridded[ which(is.finite(datgridded[, p$variable$Y]) ) , ]
    } 

    out = list()
    out$pa_w_n = pa_w_n
    out$windowsize.half = windowsize.half
    out$adims = adims
    out$xM = array( NA, dim=adims )
    out$xM[datgridded$id] = datgridded[,p$variables$Y]

    if(0){
    # ignore this for now .. as it is unlikely to be used ..
    # prediction covariates i.e., independent variables/ covariates
      pvars = c("plon", "plat", "id")
      if (p$nloccov > 0) {
        # .. not necessary except when covars are modelled locally
        for (ci in 1:p$nloccov) {
          vn = p$variables$local_cov[ci]
          pu = NULL
          pu = ecmei_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
          nts = ncol(pu)
          if ( nts== 1 ) {
            pvars = c( pvars, vn )
            datgridded[,vn] = pu[datgridded$i]  # ie. a static variable
          }
        }
      }
      datgridded = datgridded[, pvars]
      datgridded = cbind( datgridded, ecmei_timecovars ( vars=p$variables$local_all, ti=datgridded[,p$variables$TIME]  ) )

      if (p$nloccov > 0) {
        # add time-varying covars .. not necessary except when covars are modelled locally
        for (ci in 1:p$nloccov) {
          vn = p$variables$local_cov[ci]
          pu = NULL
          pu = ecmei_attach( p$storage.backend, p$ptr$Pcov[[vn]] )
          nts = ncol(pu)
          if ( nts == p$ny )  {
            datgridded$iy = datgridded$yr - p$yrs[1] + 1 #yr index
            datgridded[,vn] = pu[ cbind(datgridded$i, datgridded$iy) ]  
            message("Need to check that datgriddeda order is correct")
          } else if ( nts == p$nt ) {
            datgridded$it = p$nw*(datgridded$tiyr - p$yrs[1] - p$tres/2) + 1 #ts index
            datgridded[,vn] = pu[ cbind(datgridded$i, datgridded$it) ]  
            message("Need to check that data order is correct")
          } else if (nts==1) { } #nothing to do .. already processed above }
        }
      }
    
    }

      out$datgridded= datgridded[, datvarsout ]
    
    return(out)

}

  
