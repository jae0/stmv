    stmv_subset_distance = function( Si, p ) {

      E = stmv_error_codes()
      Sloc = stmv_attach( p$storage.backend, p$ptr$Sloc )
      Yloc = stmv_attach( p$storage.backend, p$ptr$Yloc )
      Y = stmv_attach( p$storage.backend, p$ptr$Y )
      Yi = stmv_attach( p$storage.backend, p$ptr$Yi )  # initial indices of good data

      # find data nearest Sloc[Si,] and with sufficient data
      out = NULL
      dlon = abs( Sloc[Si,1] - Yloc[Yi[],1] )
      dlat = abs( Sloc[Si,2] - Yloc[Yi[],2] )
      ndata = 0
      for ( stmv_distance_cur in p$upsampling )  {
        U = which( {dlon  <= stmv_distance_cur} & {dlat <= stmv_distance_cur} )  # faster to take a block
        ndata = length(U)
        if ( ndata >= p$n.min ) break()
      }

      yiu = Yi[U]

      out = list(flag=E[["todo"]], ndata=ndata, stmv_distance_cur=stmv_distance_cur, U=yiu, ores=NA )  # basic threshold ... now tweak it

      if (out$ndata < p$n.min) {
        # retain crude estimate and run with it
        out$flag = E[["insufficient_data"]]
        return( out )   #not enough data
      }

      # NOTE: this range is a crude estimate that averages across years (if any) ...
      o = NULL
      o = try( stmv_variogram( xy=Yloc[yiu,], z=Y[yiu,], methods=p$stmv_variogram_method, distance_cutoff=stmv_distance_cur, nbreaks=13 ) )

      if ( is.null(o)) out$flag = E[["variogram_failure"]]
      if ( inherits(o, "try-error")) out$flag = E[["variogram_failure"]]
      if ( !exists(p$stmv_variogram_method, o)) out$flag =  E[["variogram_failure"]]

      if ( out$flag ==  E[["variogram_failure"]] ) {
        o = try( stmv_variogram( xy=Yloc[yiu,], z=Y[yiu,], methods="fast", distance_cutoff=stmv_distance_cur, nbreaks=13 ) )
        if ( is.null(o)) out$flag = E[["variogram_failure"]]
        if ( inherits(o, "try-error")) out$flag = E[["variogram_failure"]]
        if ( !exists(p$stmv_variogram_method, o)) out$flag =  E[["variogram_failure"]]
      }
      if ( out$flag ==  E[["variogram_failure"]] ) return(out)

      # if (debugging) print( paste("... range=", round(ores[['range']],3), ", ", nu=", ores$nu, ", phi=", ores$phi, ndata=", ndata ) )

      if ( !exists("range_ok", o[[p$stmv_variogram_method]]) ) {
        out$flag =  E[["variogram_range_limit"]]
        return( out )     #
      }

      if ( !o[[p$stmv_variogram_method]][["range_ok"]] ) {
        # retain crude estimate and run with it
        out$flag =  E[["variogram_range_limit"]]
        return(out)
      }

      vario_stmv_distance_cur = o[[p$stmv_variogram_method]][["range"]]
      vario_U = which( {dlon  <= vario_stmv_distance_cur } & {dlat <= vario_stmv_distance_cur} )  # dlon dlat indexed on Yi
      vario_ndata =length(vario_U)


      if (vario_ndata < p$n.min) {
        # retain crude estimate and run with it
        out$flag =  E[["variogram_failure"]]
        return(out)
      }

      if (vario_ndata <= p$n.max & vario_ndata >= p$n.min) {
        # all good .. update and return
        out$U  = Yi[vario_U]
        out$ndata = vario_ndata
        out$ores = o[[p$stmv_variogram_method]]
        out$stmv_distance_cur = vario_stmv_distance_cur
        return( out)
      }

      # last try, we are here because (vario_ndata > p$n.max)
      # try to trim
      yi_vario_u = Yi[vario_U]
      if ( exists("TIME", p$variables)) {
        Ytime = stmv_attach( p$storage.backend, p$ptr$Ytime )
        iU = stmv_discretize_coordinates( coo=cbind(Yloc[yi_vario_u,], Ytime[yi_vario_u]), ntarget=p$n.max, minresolution=p$minresolution, method="thin" )
      } else {
        iU = stmv_discretize_coordinates( coo=Yloc[yi_vario_u,], ntarget=p$n.max, minresolution=p$minresolution, method="thin" )
      }

      vario_U = vario_U[iU]
      vario_ndata = length(vario_U)

      if (vario_ndata < p$n.min) {
        # retain crude estimate and run with it
        out$flag =  E[["insufficient_data"]]
        return(out)
      }

      if (vario_ndata > p$n.max) {
        # force via a random subsample
        out$U = Yi[ vario_U[ .Internal( sample( length(vario_U), p$n.max, replace=FALSE, prob=NULL)) ] ]# simple random
        out$ndata = p$n.max
        out$ores = o[[p$stmv_variogram_method]]
        out$stmv_distance_cur = vario_stmv_distance_cur
        return(out)
      }

      if (vario_ndata <= p$n.max & vario_ndata >= p$n.min) {
        # all good .. update and return
        out$U = Yi[vario_U]
        out$ndata = vario_ndata
        out$ores = o[[p$stmv_variogram_method]]
        out$stmv_distance_cur = vario_stmv_distance_cur
        return( out)
      }

      return( error("Something unexpected went wrong in stmv_subset_distance") )
    }
