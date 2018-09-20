stmv_predictions_update = function(p, preds ) {

  # update SD estimates of predictions with those from other locations via the
  # incremental  method ("online algorithm") of mean estimation after Knuth ;
  # see https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
  # update means: inverse-variance weighting
  # see https://en.wikipedia.org/wiki/Inverse-variance_weighting

  P = stmv_attach( p$storage.backend, p$ptr$P )
  Pn = stmv_attach( p$storage.backend, p$ptr$Pn )
  Psd = stmv_attach( p$storage.backend, p$ptr$Psd )

  sflag = "good"

  npred = nrow(preds)

  if ( ! exists("TIME", p$variables) ) {

    u = which( is.finite( P[preds$i] ) )  # these have data already .. update
    if ( length( u ) > 1 ) {
      ui = preds$i[u]  # locations of P to modify
      Pn[ui] = Pn[ui] + 1 # update counts
      stdev_update =  Psd[ui] + ( preds$sd[u] -  Psd[ui] ) / Pn[ui]
      means_update = ( P[ui] / Psd[ui]^2 + preds$mean[u] / preds$sd[u]^2 ) / ( Psd[ui]^(-2) + preds$sd[u]^(-2) )
      mm = which(is.finite( means_update + stdev_update ))
      if( length(mm)> 0) {
        iumm = ui[mm]
        Psd[iumm] = stdev_update[mm]
        P  [iumm] = means_update[mm]
      }
      stdev_update = NULL
      means_update = NULL
      ui=mm=iumm=NULL
    }

    # first time # no data yet
    v = setdiff(1:npred, u)
    if ( length(v) > 0 ) {
      vi = preds$i[v]
      Pn [vi] = 1
      P  [vi] = preds$mean[v]
      Psd[vi] = preds$sd[v]
    }
    vi = NULL
  }


  if ( exists("TIME", p$variables) ) {
    u = which( is.finite( P[preds$i,1] ) )  # these have data already .. update
    u_n = length( u )
    if ( u_n > 1 ) {  # ignore if only one point .. mostly because it can cause issues with matrix form ..
      # locations of P to modify
      ui = sort(unique(preds$i[u]))
      nc = ncol(P)
      if (p$storage.backend == "ff" ) {
        add.ff(Pn, 1, ui, 1:nc ) # same as Pn[ui,] = Pn[ui]+1 but 2X faster
      } else {
        Pn[ui,] = Pn[ui,] + 1
      }
      stdev_update =  Psd[ui,] + ( preds$sd[u] -  Psd[ui,] ) / Pn[ui,]
      means_update = ( P[ui,] / Psd[ui,]^2 + preds$mean[u] / preds$sd[u]^2 ) / ( Psd[ui,]^(-2) + preds$sd[u]^(-2) )

      updates = means_update + stdev_update
      if (!is.matrix(updates)) {
        sflag = "error"
        u = u_n = ui = nc =stdev_update = means_update = NULL
        message( "update of predictions were problematic ... this should not happen, proabbaly due to NA's" )
        next()
      }

      mm = which( is.finite( rowSums(updates)))  # created when preds go outside quantile bounds .. this removes all data from a given location rather than the space-time .. severe but likely due to a poor prediction and so remove all (it is also faster this way as few manipulations)
      if( length(mm)> 0) {
        iumm = ui[mm]
        Psd[iumm,] = stdev_update[mm,]
        P  [iumm,] = means_update[mm,]
        iumm = NULL
      }
      stdev_update = means_update = updates = ui = mm=NULL
    }

    # do this as a second pass in case NA's were introduced by the update .. unlikely , but just in case
    v = which( !is.finite( P[preds$i,1] ) )  # these have data already .. update
    nv = length(v)          # no data yet
    if ( nv > 0 ) {
      vi = sort(unique(preds$i[v]))
      Pn [vi,] = 1
      P  [vi,] = preds$mean[v]
      Psd[vi,] = preds$sd[v]
      vi = NULL
    }
  }

    if (0) {
        v = preds
        if ( exists("TIME", p$variables) ){
          v = v[which( v[,p$variables$TIME]==2000.55),]
        }
        require(lattice)
        print(
          levelplot( mean ~ plon+plat, v, aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=TRUE) )
      )

      if ("time slice at 2012.05") {
        lattice::levelplot( mean ~ plon + plat, data=preds[preds[,p$variables$TIME]==2012.05,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" )
      }
      if ("all TIME time slices from latest predictions") {
        for( i in sort(unique(preds[,p$variables$TIME])))  {
          print(lattice::levelplot( mean ~ plon + plat, data=preds[preds[,p$variables$TIME]==i,], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
        }
      }
      if ("all nt time slices in stored predictions P") {
        Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
        # pa comes from stmv_interpolate ... not carried here
        for (i in 1:p$nt) {
          print( lattice::levelplot( P[pa$i,i] ~ Ploc[pa$i,1] + Ploc[ pa$i, 2], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
        }
      }
      if ("no time slices in P") {
        Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
          print( lattice::levelplot( P[pa$i] ~ Ploc[pa$i,1] + Ploc[ pa$i, 2], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
      }
    }

  return(sflag)
}
