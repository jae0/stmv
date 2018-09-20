    stmv_subset_distance = function( Si, upsampling, n.min=30, n.max=10000, vgm_method="fast", minresolution=NULL ) {

      E = stmv_error_codes()
      Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
      Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )
      Y = stmv_attach( p$storage.backend, p$ptr$Y )
      Yi = stmv_attach( p$storage.backend, p$ptr$Yi )

      # find data nearest Sloc[Si,] and with sufficient data
      out = NULL
      dlon = abs( Sloc[Si,1] - Yloc[Yi[],1] )
      dlat = abs( Sloc[Si,2] - Yloc[Yi[],2] )
      ndata = 0
      for ( stmv_distance_cur in upsampling )  {
        U = which( {dlon  <= stmv_distance_cur} & {dlat <= stmv_distance_cur} )  # faster to take a block
        ndata = length(U)
        if ( ndata >= n.min ) break()
      }

      yiu = Yi[U]

      out = list(flag=E[["todo"]], ndata=ndata, stmv_distance_cur=stmv_distance_cur, U=yiu, ores=NA )  # basic threshold ... now tweak it

      if (out$ndata < n.min) {
        # retain crude estimate and run with it
        out$flag = E[["insufficient_data"]]
        return( out )   #not enough data
      }

      # NOTE: this range is a crude estimate that averages across years (if any) ...
      o = NULL
      o = try( stmv_variogram( xy=Yloc[yiu,], z=Y[yiu,], methods=vgm_method, distance_cutoff=stmv_distance_cur, nbreaks=13 ) )

      if ( is.null(o)) {
        # retain crude estimate and run with it
        out$flag = E[["variogram_failure"]]
        return( out )     # fast variogram did not work
      }
      if ( inherits(o, "try-error")) {
        # retain crude estimate and run with it
        out$flag = E[["variogram_failure"]]
        return( out )     # fast variogram did not work
      }
      if ( !exists(vgm_method, o)) {
        # retain crude estimate and run with it
        out$flag =  E[["variogram_failure"]]
        return( out )     # fast variogram did not work
      }

      # if (debugging) print( paste("... range=", round(ores[['range']],3), ", ", nu=", ores$nu, ", phi=", ores$phi, ndata=", ndata ) )

      if ( !exists("range_ok", o[[vgm_method]]) ) {
        out$flag =  E[["variogram_range_limit"]]
        return( out )     #
      }

      if ( !o[[vgm_method]][["range_ok"]] ) {
        # retain crude estimate and run with it
        out$flag =  E[["variogram_range_limit"]]
        return(out)
      }

      vario_stmv_distance_cur = o[[vgm_method]][["range"]]
      vario_U = which( {dlon  <= vario_stmv_distance_cur } & {dlat <= vario_stmv_distance_cur} )  # dlon dlat indexed on Yi
      vario_ndata =length(vario_U)


      if (vario_ndata < n.min) {
        # retain crude estimate and run with it
        out$flag =  E[["variogram_failure"]]
        return(out)
      }

      if (vario_ndata <= n.max & vario_ndata >= n.min) {
        # all good .. update and return
        out$U  = Yi[vario_U]
        out$ndata = vario_ndata
        out$ores = o[[vgm_method]]
        out$stmv_distance_cur = vario_stmv_distance_cur
        return( out)
      }

      # last try, we are here because (vario_ndata > n.max)
      # try to trim
      yi_vario_u = Yi[vario_U]
      if ( exists("TIME", p$variables)) {
        Ytime = stmv_attach( p$storage.backend, p$ptr$Ytime )
        iU = stmv_discretize_coordinates( coo=cbind(Yloc[yi_vario_u,], Ytime[yi_vario_u]), ntarget=n.max, minresolution=minresolution, method="thin" )
      } else {
        iU = stmv_discretize_coordinates( coo=Yloc[yi_vario_u,], ntarget=n.max, minresolution=minresolution, method="thin" )
      }

      vario_U = vario_U[iU]
      vario_ndata = length(vario_U)

      if (vario_ndata < n.min) {
        # retain crude estimate and run with it
        out$flag =  E[["insufficient_data"]]
        return(out)
      }

      if (vario_ndata > n.max) {
        # force via a random subsample
        out$U = Yi[ vario_U[ .Internal( sample( length(vario_U), n.max, replace=FALSE, prob=NULL)) ] ]# simple random
        out$ndata = n.max
        out$ores = o[[vgm_method]]
        out$stmv_distance_cur = vario_stmv_distance_cur
        return(out)
      }

      if (vario_ndata <= n.max & vario_ndata >= n.min) {
        # all good .. update and return
        out$U = Yi[vario_U]
        out$ndata = vario_ndata
        out$ores = o[[vgm_method]]
        out$stmv_distance_cur = vario_stmv_distance_cur
        return( out)
      }

      return( error("Something went wrong") )
    }
