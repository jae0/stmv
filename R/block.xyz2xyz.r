
  # ----------------------------------------------------
  # various ways of blocking data ...

  block.xyz2xyz = function( xyz,  params, outfile="xyz.block", engine="R", method="mean", variables="z", getdata=F ) {

    if (engine=="R") {
      out = with (params, {
        plon = "plon"
        plat = "plat"
        xyz[,plon] = trunc(xyz[,plon]/pres)*pres
        xyz[,plat] = trunc(xyz[,plat]/pres)*pres
        id = paste(xyz[,plon], xyz[plat])
        dups = id[duplicated(id)]
        if (!is.null(length(dups))) {
          for (di in dups) {
            indices = which(id==di)
            for (vi in variables) {xyz[indices[1],vi] = mean(xyz[indices,vi], na.rm=T) }
            xyz[indices[-1],] = NA  # remaining get set to null
          }
        }
        out = xyz
        return (out)
      })
    }

    if (engine=="R2") {
        if( dim(xyz)[2] == 2) xyz$z=1
        plon = "plon"
        plat = "plat"
        pres = params$pres

        xyz[,plon] = trunc(xyz[,plon]/pres)*pres
        xyz[,plat] = trunc(xyz[,plat]/pres)*pres
        out = NULL
        for (vi in variables) {

          m = tapply( X=xyz[, vi], INDEX=list(xyz$plon, xyz$plat),
                      FUN=function(q) { mean ( q, na.rm=T ) }, simplify=T )
          om = as.data.frame( as.table (m) )
          names(om) = c(plon, plat, vi)
          if (is.null (out))  out = om[,c(plon, plat)]
          names0 = names(out)
          out = cbind( out, om[,vi] )
          rm (m, om); gc()
          
          s = tapply( X=xyz[, vi], INDEX=list(xyz$plon, xyz$plat),
                      FUN=function(q) { sd ( q, na.rm=T ) }, simplify=T )
          os = as.data.frame( as.table (s) )
          names(os) = c(plon, plat, vi)
          out = cbind( out, os[,vi] )
          rm (s, os); gc()
            
          n = tapply( X=xyz[, vi], INDEX=list(xyz$plon, xyz$plat),
                      FUN=function(q) { length ( is.finite (q) )  }, simplify=T )
          on = as.data.frame( as.table (n) )
          names(on) = c(plon, plat, vi)
          out = cbind( out, on[,vi] )
          rm (n, on); gc()
 
          names(out) = c( names0, vi, paste(vi, c("sd", "n"), sep=".") )
        }
        
        out[,plon] = as.numeric(as.character( out[,plon] ))
        out[,plat] = as.numeric(as.character( out[,plat] ))
        out = out[ is.finite( out[, 3] ) ,]  # <<<<<----- ie. order of variables is important .. the 1st named is used
    }

    if (engine=="R.geo") {

      if( dim(xyz)[2] == 2) xyz$z=1
        out = NULL
        for (vi in variables) {
          m = tapply( X=xyz[, vi], INDEX=list(xyz$lon, xyz$lat),
                      FUN=function(q) { mean ( q, na.rm=T ) }, simplify=T )
          s = tapply( X=xyz[, vi], INDEX=list(xyz$lon, xyz$lat),
                      FUN=function(q) { sd ( q, na.rm=T ) }, simplify=T )
          n = tapply( X=xyz[, vi], INDEX=list(xyz$lon, xyz$lat),
                      FUN=function(q) { length ( is.finite (q) )  }, simplify=T )
          om = as.data.frame( as.table (m) ); names(om) = c("lon", "lat", vi)
          os = as.data.frame( as.table (s) ); names(os) = c("lon", "lat", vi)
          on = as.data.frame( as.table (n) ); names(on) = c("lon", "lat", vi)
          if (is.null (out)) {
            out = om[,c("lon", "lat")]
          }
          names0 = names(out)
          out = cbind( out, om[,vi], os[,vi], on[,vi] )
          names(out) = c( names0, vi, paste(vi, c("sd", "n"), sep=".") )
        }

        out[,"lon"] = as.numeric(as.character( out[,"lon"] ))
        out[,"lat"] = as.numeric(as.character( out[,"lat"] ))
        out = out[ is.finite( out[, 3] ) ,]  # <<<<<----- ie. order of variables is important .. the 1st named is used

    }

    if (engine=="R.geo2") {
      lon = "lon"
      lat = "lat"
      id = paste(xyz[,lon], xyz[lat])
      dups = id[duplicated(id)]
      if (!is.null(length(dups))) {
        for (di in dups) {
          indices = which(id==di)
          for (vi in variables) xyz[indices[1],vi] = mean(xyz[indices,vi], na.rm=T)
            xyz[indices[-1],] = NA  # remaining get set to null
        }
      }
      out = xyz
    }

    if (engine=="R.geo3") {

      if( dim(xyz)[2] == 2) xyz$z=1
        out = NULL
        for (vi in variables) {
          m = tapply( X=xyz[, vi], INDEX=list(xyz$lon, xyz$lat),
                      FUN=function(q) { mean ( q, na.rm=T ) }, simplify=T )
          out = as.data.frame( as.table (m) )
          names(out) = c("lon", "lat", vi)
        }

        out[,"lon"] = as.numeric(as.character( out[,"lon"] ))
        out[,"lat"] = as.numeric(as.character( out[,"lat"] ))
        out = out[ is.finite( out[, 3] ) ,]  # <<<<<----- ie. order of variables is important .. the 1st named is used
    }

    if (engine=="R.geo.median") {

      if( dim(xyz)[2] == 2) xyz$z=1
        out = NULL
        for (vi in variables) {
          m = tapply( X=xyz[, vi], INDEX=list(xyz$lon, xyz$lat),
                      FUN=function(q) { median ( q, na.rm=T ) }, simplify=T )
          out = as.data.frame( as.table (m) )
          names(out) = c("lon", "lat", vi)
        }

        out[,"lon"] = as.numeric(as.character( out[,"lon"] ))
        out[,"lat"] = as.numeric(as.character( out[,"lat"] ))
        out = out[ is.finite( out[, 3] ) ,]  # <<<<<----- ie. order of variables is important .. the 1st named is used
    }


    return (out)

  }


