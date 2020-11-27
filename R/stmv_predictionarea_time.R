stmv_predictionarea_time = function( p, pa ) {

  pa_n = nrow(pa) # n of space only
  pa = cbind( pa[ rep.int(1:pa_n, p$nt), ], rep.int(p$prediction_ts, rep(pa_n, p$nt )) )
  names(pa)[ncol(pa)] = p$stmv_variables$TIME
  pa = cbind( pa, stmv_timecovars( vars=p$stmv_variables$local_all, ti=pa[,p$stmv_variables$TIME]  ) )

  if (p$nloccov > 0) {
    # add time-varying covars .. not necessary except when covars are modelled locally
    for (ci in 1:p$nloccov) {
      vn = p$stmv_variables$local_cov[ci]
      pu = NULL
      pu = stmv_attach( p$storage_backend, p$ptr$Pcov[[vn]] )
      nts = ncol(pu)
      if ( nts == p$ny )  {
        pa$iy = pa$yr - p$yrs[1] + 1 #yr index
        pa[,vn] = pu[ cbind(pa$i, pa$iy) ]
        message("Need to check that data order is correct")
      } else if ( nts == p$nt ) {
        pa$it = p$nw*(pa$tiyr - p$yrs[1] - p$tres/2) + 1 #ts index
        pa[,vn] = pu[ cbind(pa$i, pa$it) ]
        message("Need to check that data order is correct")
      } else if (nts==1) { } #nothing to do .. already processed above }
    }
  }
  return(pa)
}
