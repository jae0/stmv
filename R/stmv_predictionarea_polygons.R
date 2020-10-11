

stmv_predictionarea_polygons = function(p, sloc, global_sppoly=NULL, windowsize.half=0, stmv_au_buffer_links=0, stmv_au_distance_reference="none" ) {

  pa_coord_names = p$stmv_variables$LOCS[1:2]

  if (is.null( global_sppoly )) {
    # essentially identical to stmv_predictionarea_lattice, but for the addition of sppoly before addtion time
    pa = stmv_predictionarea_space( p=p, sloc=sloc, windowsize.half=windowsize.half )
    if (is.null(pa)) return(NULL)

    ## now add polygon representation of spatial lattice

    dx = p$pres
    dy = p$pres

    x_r = range(pa[, pa_coord_names[1]])
    x_c = range(pa[, pa_coord_names[2]])

    nr = aegis_floor( diff(x_r)/dx ) + 1L
    nc = aegis_floor( diff(x_c)/dy ) + 1L
    # check: dr = diff(x_r)/(nr-1) == dx  ;; dc = diff(x_c)/(nc-1) # == dy

    # default behaviour .. lattice grid
      # data_subset = NULL
    sppoly = sf::st_as_sf( pa, coords=pa_coord_names )
    sf::st_crs(sppoly) = sf::st_crs( p$aegis_proj4string_planar_km )
    sp_grid = raster::raster( sppoly, nrows=nr, ncols=nc )
    sppoly = as(sp_grid, "SpatialPolygonsDataFrame")

    if (exists("i", slot(sppoly, "data")) ) {
      sppoly$AUID = sppoly$i
    } else {
      sppoly$AUID = 1:nrow(sppoly)
    }
    sppoly$AUID = as.character( sppoly$AUID )
    # poly* function operate on Spatial* data
    nb = poly2nb(sppoly, row.names=sppoly$AUID, queen=TRUE)  # slow .. ~1hr?
    nb.remove = which(card(nb) == 0)

    if ( length(nb.remove) > 0 ) {
      # remove isolated locations and recreate sppoly .. alternatively add links to nb
      nb.keep = which(card(nb) > 0)
      nb = nb_remove( nb, nb.remove )
      sppoly = sppoly[nb.keep,]
      row.names(sppoly) = sppoly$AUID
      sppoly = sp::spChFIDs( sppoly, row.names(sppoly) )  #fix id's
      # sppoly = sppoly[ order(sppoly$AUID), ]
    }

    attr(sppoly, "nb") = nb  # adding neighbourhood as an attribute to sppoly
    nb =NULL
    nb.remove =NULL

    attr( pa, "sppoly" ) = sppoly

    if ( exists("TIME", p$stmv_variables) )  pa = stmv_predictionarea_time( p=p, pa=pa )

    rownames(pa) = NULL
    return(pa)

  }


  # ------------


  if ( !is.null( global_sppoly )) {

      au = as( global_sppoly, "SpatialPolygonsDataFrame")
      au$AUID = 1:nrow(au)

      sloc_sp = SpatialPoints( t(sloc), sp::CRS( p$aegis_proj4string_planar_km ) )

      nbfocal = over( sloc_sp, au )$AUID

      if (is.na(nbfocal)) {
        message( "no data" )
        return (NULL)
      }
      nb = attr(au, "nb")  # full matrix
      nbnames = attr( nb, "region.id")
      nnAUID = nbnames[nb[[which( nbnames == nbfocal )]]]  ## nearest neighbours
      tokeep = unique( c( nnAUID, nbfocal ) )

      if ( stmv_au_buffer_links > 0 ) {
        # no of additional neighbourhood links ... 0 == nearest neighbours, 1 == nn + next nearest neighbours, etc
        for (i in 1:nlinks ) {
          new = NULL
          for ( foc in nnAUID ) {
            new = c( new, nbnames[nb[[which( nbnames == foc  )]]] ) ## nearest neighbours
          }
          tokeep = unique( c(tokeep, new) )
          nnAUID = tokeep
        }
      }

      if ( windowsize.half > 0 ) {

        if (stmv_au_distance_reference=="none") {
          # nothing to do
        }

        if (stmv_au_distance_reference=="centroid") {
          # distance based filtering based on centroids
          aucoo = coordinates( au )
          inrange =  which( (abs(aucoo[,1] - sloc[1]) <= windowsize.half) &  (abs(aucoo[,2] - sloc[2]) <= windowsize.half) )
          todrop = setdiff( nbnames, nbnames[inrange] )
          tokeep = setdiff( tokeep, todrop)
        }

        if (stmv_au_distance_reference=="inside_or_touches_boundary") {
          # distance based filtering based on centroids
          ausf = as( au, "sf")
          foc = st_buffer( ausf[ which(ausf$AUID==nbfocal), ], dist= windowsize.half )
          inrange =  which( unlist( list2DF( st_intersects( ausf, foc ) ) ) ==1 )
          todrop = setdiff( nbnames, nbnames[inrange] )
          tokeep = setdiff( tokeep, todrop)
        }

        if (stmv_au_distance_reference=="completely_inside_boundary") {
          # distance based filtering based on centroids
          ausf = as( au, "sf")
          foc = st_buffer( ausf[ which(ausf$AUID==nbfocal), ], dist= windowsize.half )
          inrange =  which( unlist( list2DF( st_contains( ausf, foc ) ) ) ==1 )
          todrop = setdiff( nbnames, nbnames[inrange] )
          tokeep = setdiff( tokeep, todrop)
        }

      }

      sppoly = au[ order( match( tokeep, au$AUID ) ), ]
      row.names(sppoly) = sppoly$AUID

      pa = data.frame( coordinates( sppoly ) )
      names(pa) = pa_coord_names
      pa$AUID = sppoly$AUID
      pa$i = match( sppoly$AUID, global_sppoly$AUID )

      # prediction covariates i.e., independent stmv_variables/ covariates
      pvars = c(pa_coord_names, "i")
      if (p$nloccov > 0) {
        # .. not necessary except when covars are modelled locally
        for (ci in 1:p$nloccov) {
          vn = p$stmv_variables$local_cov[ci]
          pu = NULL
          pu = stmv_attach( p$storage_backend, p$ptr$Pcov[[vn]] )
          nts = ncol(pu)
          if ( nts== 1 ) {
            pvars = c( pvars, vn )
            pa[,vn] = pu[pa$i]  # i.e., a static variable
          }
        }
      }
      pa = as.data.frame(pa[, ..pvars])
      # completed, reconstruction of spatial vars

      # Time next
      if ( exists("TIME", p$stmv_variables) )  pa = stmv_predictionarea_time( p=p, pa=pa )

      # poly* function operate on Spatial* data
      #       nb = attr(au, "nb")  # full matrix
      # nbnames = attr( nb, "region.id")

      nnAUID = nbnames[nb[[which( nbnames == nbfocal )]]]  ## nearest neighbours
      tokeep = unique( c( nnAUID, nbfocal ) )

      nb.remove = which( ! (nbnames %in% sppoly$AUID ) )
      if ( length(nb.remove) > 0 ) {
        # remove isolated locations and recreate sppoly .. alternatively add links to nb
        nb = nb_remove( nb, nb.remove )
        sppoly = sp::spChFIDs( sppoly, row.names(sppoly) )  #fix id's
        # sppoly = sppoly[ order(sppoly$AUID), ]
      }

      attr(sppoly, "nb") = nb  # adding neighbourhood as an attribute to sppoly
      nb =NULL
      nb.remove =NULL

      attr( pa, "sppoly" ) = sppoly

      rownames(pa) = NULL
      return(pa)

  }

  if (0) {
      jj = which( card(nb) == 0)
      jj = match( tokeep, au$AUID )
      plot(sppoly)
      plot(sppoly[jj,], add=T, col="red")
      dev.new()
      edit(nb, polys=sppoly)
      card(nb) # last check if any more  isolated areas
  }

}
