stmv_statistics_update =  function(p, res, W, Si ) {
  # compute and extract misc summary statistics from results
  S = stmv_attach( p$storage.backend, p$ptr$S )
  Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )

  sflag = "good"

  out = list()

  if (exists("stmv_stats",  res)) out = res$stmv_stats

  if (!exists("sdSpatial", out)) {
    # some methods can generate spatial stats simultaneously ..
    # it is faster to keep them all together instead of repeating here
    # field and RandomFields gaussian processes seem most promising ...
    # default to fields for speed:
    out["sdSpatial"] = NA
    out["sdObs"] = NA
    out["range"] = NA
    out["phi"] = NA
    out["nu"] = NA
    if (!is.null(W)) {
      if ( !is.na(W[["ores"]])) {
        if ( exists("varSpatial", W[["ores"]]) ) out["sdSpatial"] = sqrt( W[["ores"]][["varSpatial"]] )
        if ( exists("varObs", W[["ores"]]) )     out["sdObs"] = sqrt(W[["ores"]][["varObs"]])
        if ( exists("range", W[["ores"]]) )      out["range"] = W[["ores"]][["range"]]
        if ( exists("phi", W[["ores"]]) )        out["phi"] = W[["ores"]][["phi"]]
        if ( exists("nu", W[["ores"]]) )         out["nu"] = W[["ores"]][["nu"]]
      }
    }

  }

  if ( exists("TIME", p$variables) ){
    # annual ts, seasonally centered and spatially
    # pa_i = which( Sloc[Si,1]==Ploc[,1] & Sloc[Si,2]==Ploc[,2] )
    pac_i = which( res$predictions$plon==Sloc[Si,1] & res$predictions$plat==Sloc[Si,2] )
    # plot( mean~tiyr, res$predictions[pac_i,])
    # plot( mean~tiyr, res$predictions, pch="." )
    out["ar_timerange"] = NA
    out["ar_1"] = NA

    if (length(pac_i) > 5) {
      pac = res$predictions[ pac_i, ]
      pac$dyr = pac[, p$variables$TIME] - trunc(pac[, p$variables$TIME] )
      piid = which( zapsmall( pac$dyr - p$dyear_centre) == 0 )
      pac = pac[ piid, c(p$variables$TIME, "mean")]
      pac = pac[ order(pac[,p$variables$TIME]),]
      if (length(piid) > 5 ) {
        ts.stat = NULL
        ts.stat = try( stmv_timeseries( pac$mean, method="fft" ) )
        if (!is.null(ts.stat) && !inherits(ts.stat, "try-error") ) {
          out["ar_timerange"] = ts.stat$quantilePeriod
          if (all( is.finite(pac$mean))) {
            afin = which (is.finite(pac$mean) )
            if (length(afin) > 5 && var( pac$mean, na.rm=TRUE) > p$eps ) {
              ar1 = NULL
              ar1 = try( ar( pac$mean, order.max=1 ) )
              if (!inherits(ar1, "try-error")) {
                if ( length(ar1$ar) == 1 ) {
                  out["ar_1"] = ar1$ar
                }
              }
            }
          }
          if ( !is.finite(out[["ar_1"]]) ) {
            ar1 = try( cor( pac$mean[1:(length(piid) - 1)], pac$mean[2:(length(piid))], use="pairwise.complete.obs" ) )
            if (!inherits(ar1, "try-error")) out["ar_1"] = ar1
          }
        }

        ### Do the logistic model here ! -- if not already done ..
        if (!exists("ts_K", out)) {
          # model as a logistic with ts_r, ts_K, etc .. as stats outputs

        }

      }
      pac=piid=NULL
    }
    pac_i=NULL
  }

  # save stats
  for ( k in 1: length(p$statsvars) ) {
    if (exists( p$statsvars[k], out )) {
      S[Si,k] = out[[ p$statsvars[k] ]]
    }
  }

  return(sflag)
}
