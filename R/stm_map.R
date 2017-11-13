
  stm_map = function( xyz, cfa.regions=T, depthcontours=T, pts=NULL, colpts=F, annot=NULL, annot.cex=2.2,
                 leg = NULL, projection = "utm20", col.regions=F, at=0:1,
                 fn=paste("map", trunc(runif(1)*1e8), sep=""), loc=tempdir(),
                 corners=NULL, rez=c(1,1), spatial.domain="SSE", display=F, save=T, pt.cex=0.5, pt.pch=16, pt.col='black',
                 colorkey=NULL, fill=T, scalebar=NULL, ... ) {

    # map using levelplot ... no GMT dependency

		require( lattice )

    xlim =ylim = NULL
    if(is.null(colorkey)) colorkey=list(space="right", labels=list(cex=3)) # these are lattice options

    if (ncol( xyz) == 2) { # assuming points
      xyz = cbind( xyz, 1)
      colorkey=F
    }

    if ( is.null( xyz$plon) ) {
      names( xyz ) = c( "lon", "lat", "z" )
      xyz = lonlat2planar( xyz,   proj.type=projection)
      xyz = xyz[ , c( "plon", "plat", "z" ) ]
      if ( !is.null(corners) ) {
        if ( is.null(corners$plon) ) corners = lonlat2planar( corners,  proj.type= projection )
      }

    } else {
      names( xyz ) = c( "plon", "plat", "z" )
    }

    if ( is.logical( colorkey) )  pts = xyz      # flag from above when only XY data are passed,.. after conversion to planar coords

    if ( !is.null(corners) ) {
      xlim = range( corners$plon)
      ylim = range( corners$plat)
    } else {
      xlim = range( xyz$plon )
      ylim = range( xyz$plat )
    }

    if ( ! is.null(pts) ) {
      if ( is.null( pts$plon) ) {
        pts = lonlat2planar(xyz,  proj.type=p$internal.projection)
        pts = pts[, c("plon", "plat")]
      }
    }

    if(fill){
      xyz$z[xyz$z>max(at)]=max(at)
      xyz$z[xyz$z<min(at)]=min(at)
    }


    lp = levelplot( z ~ plon+plat, data=xyz, aspect="iso", pts=pts, colpts=colpts, annot=annot, spatial.domain=spatial.domain,
      annot.cex=annot.cex, xlab="", ylab="", scales=list(draw=F), col.regions=col.regions, at=at, xlim=xlim, ylim=ylim,
      colorkey=colorkey , rez=rez, leg=leg,  cfa.regions=cfa.regions,
      panel = function(x, y, z, rez=rez,  ...) {

        panel.levelplot (x, y, z, aspect="iso", rez=rez, ...)

        if ( !is.null(pts) ) {
          if (colpts) {  # overlay above with larger sized points, esp if there is colour
            colb = findInterval( z, at)
            for ( ii in 1:length(z) ) {
              panel.xyplot( pts$plon[ii], pts$plat[ii],  aspect="iso", col=col.regions[colb[ii]],
                  panel = panel.rect, height = rez[1], width = rez[2], cex=pt.cex,... )
            }
          } else {
              panel.xyplot( pts$plon, pts$plat,  aspect="iso", col = pt.col, pch=pt.pch,
                  panel = panel.rect, height = rez[1], width = rez[2], cex=pt.cex, ... )
          }
        }

        pp = spatial_parameters( type=spatial.domain )

        if (depthcontours) {
          isobs = isobath.db( p=pp, depths=c( 100, 200, 300, 400, 500, 600, 700 ), crs=pp$internal.crs )
          depths1 = c(100, 300, 500, 700 )
          depths2 = c(200, 400, 600)
          for ( i in depths1 ) sp.lines( isobs[as.character(i) ] , col = rgb(0.2,0.2,0.2,0.5), cex=0.6 )
          for ( i in depths2 ) sp.lines( isobs[as.character(i) ] , col = rgb(0.3,0.3,0.3,0.5), cex=0.6 )
        }

        if ( cfa.regions ) {
          # coords of boundaries .. to be moved to polygon database ...
          cfa.nens.23 = data.frame( rbind(
            c(-59.85, 46),
            c(-58.40, 46)
          ))
          cfa.23.24 = data.frame( rbind(
            c(-59.065, 43.5),
            c(-59.959007, 44.829624),
            c(-60.51667, 45.61667)
          ))
          cfa.4x.24 = data.frame( rbind(
            c( -63.333333, 42.61379),
            c( -63.333333,	44.332904),
            c( -63.50242, 44.502358)
          ))

          names( cfa.nens.23 ) = names( cfa.23.24 ) = names( cfa.4x.24 ) = c("lon", "lat")

          cfa.nens.23 = lonlat2planar( cfa.nens.23, proj.type=pp$internal.projection )
          cfa.23.24 = lonlat2planar( cfa.23.24, proj.type=pp$internal.projection )
          cfa.4x.24 = lonlat2planar( cfa.4x.24,  proj.type=pp$internal.projection )

          panel.lines( cfa.nens.23$plon, cfa.nens.23$plat, col = "black", lwd=1,lty=2 )
          panel.lines( cfa.23.24$plon, cfa.23.24$plat, col = "black", lwd=1,lty=2 )
          panel.lines( cfa.4x.24$plon, cfa.4x.24$plat, col = "black", lwd=1,lty=2 )

        }

        #coastline
        coast = emaf::coastline.db(p=pp, crs=pp$internal.crs, DS="gshhg coastline highres" )
        sp.polygons( coast, col = "black", cex=1 ,fill='grey')

         
        #legend
        lx = xlim[2]-xlim[1]
        ly = ylim[2]-ylim[1]
        if (is.null(leg) )leg = c( xlim[2]-0.1*lx, ylim[1] + 0.1*ly )
        
        if (!is.null(scalebar)) {
          lx = xlim[2]-xlim[1]
          ly = ylim[2]-ylim[1]
          if (is.null(leg) )leg = c( xlim[2]-0.1*lx, ylim[1] + 0.1*ly )
          panel.arrows( x0=leg[1]-scalebar, y0=leg[2], x1=leg[1], y1=leg[2], angle=90, length=0.06, ends="both", lwd=3, col="black", ...)
          panel.text( x=leg[1]-scalebar/2 , y=leg[2]+0.05*ly , paste(scalebar,"km"), cex=1.7 )
        }

        if ( !is.null( annot ) ){
          panel.text( x=leg[1]-scalebar/2, y=leg[2]-0.05*ly, annot, cex=2 )  # pos=2 is left of (right justified)
        }
    } # end panel
    ) # end levelplot

    if(display)print(lp)

    if(save){
      dir.create (loc, showWarnings=FALSE, recursive =TRUE)
      fn = file.path( loc, paste(fn, "png", sep="." ) )
      png(  filename=fn, width=3072, height=2304, pointsize=40, res=300 )
      print(lp)
      dev.off()
    }

    return( fn )

  }


