

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

    nr = trunc( diff(x_r)/dx ) + 1L
    nc = trunc( diff(x_c)/dy ) + 1L
    # check: dr = diff(x_r)/(nr-1) == dx  ;; dc = diff(x_c)/(nc-1) # == dy

    # default behaviour .. lattice grid
      # data_subset = NULL
    sppoly = sf::st_as_sf( pa, coords=pa_coord_names, crs=st_crs( p$aegis_proj4string_planar_km) )

if (0) {
    # not using raster -- incomplete

    sppoly = (
      st_make_grid( sppoly, cellsize=areal_units_resolution_km,  what="polygons", square=TRUE )
      %>% st_as_sf( crs=st_crs( p$aegis_proj4string_planar_km ))
    )
    spdf0 = sf::st_as_sf( pa, coords=pa_coord_names, crs=st_crs( p$aegis_proj4string_planar_km ) )
    sppoly = st_as_sf( st_make_grid( spdf0, n=c(nr, nc),  what="polygons", square=TRUE ) )
    sppoly$internal_id = 1:nrow(sppoly)
    spdf0$internal_id = st_points_in_polygons( spdf0, sppoly, varname="internal_id" )
    o = match( sppoly$internal_id,spdf0$internal_id )

      sppoly = st_as_sf( st_make_grid( Z, cellsize=areal_units_resolution_km,  what="polygons", square=TRUE ) )
      sppoly$AUID = as.character( 1:nrow(sppoly) )  # row index
      row.names(sppoly) = sppoly$AUID

      require(raster)
      raster_template = raster::raster(extent(Z), res=areal_units_resolution_km, crs=projection(Z) ) # +1 to increase the area
      sppoly = raster::rasterize( Z, raster_template, field=Z$z )  #NOTE :: TODO : move to stars::st_rasterize 
      sppoly = as(sppoly, "SpatialPixelsDataFrame")
      sppoly = as( as(sppoly, "SpatialPolygonsDataFrame"), "sf")
      raster_template = NULL

}

    sppoly = as( as( raster::raster( sppoly, nrows=nr, ncols=nc ), "SpatialPolygonsDataFrame" ), "sf" )

    if (exists("i", sppoly) ) {
      sppoly$AUID = sppoly$i
    } else {
      sppoly$AUID = 1:nrow(sppoly)
    }
    sppoly$AUID = as.character( sppoly$AUID )
    # poly* function operate on Spatial* data
    NB_graph = poly2nb(sppoly, row.names=sppoly$AUID, queen=TRUE)  # slow .. ~1hr?
    NB_graph.remove = which(card(NB_graph) == 0)

    if ( length(NB_graph.remove) > 0 ) {
      # remove isolated locations and recreate sppoly .. alternatively add links to NB_graph
      NB_graph.keep = which(card(NB_graph) > 0)
      NB_graph = nb_remove( NB_graph, NB_graph.remove )
      sppoly = sppoly[NB_graph.keep,]
      row.names(sppoly) = sppoly$AUID
      sppoly = sp::spChFIDs( sppoly, row.names(sppoly) )  #fix id's
      # sppoly = sppoly[ order(sppoly$AUID), ]
    }

    attr(sppoly, "NB_graph") = NB_graph  # adding neighbourhood as an attribute to sppoly
    NB_graph =NULL
    NB_graph.remove =NULL

    attr( pa, "sppoly" ) = sppoly

    if ( exists("TIME", p$stmv_variables) )  pa = stmv_predictionarea_time( p=p, pa=pa )

    rownames(pa) = NULL
    return(pa)

  }


  # ------------


  if ( !is.null( global_sppoly )) {

      au = as( global_sppoly, "sf")
      au$AUID = 1:nrow(au)

      nbfocal = st_points_in_polygons(
        pts = st_transform( st_as_sf( t(sloc), coords=c(1,2), crs=st_crs(p$aegis_proj4string_planar_km) ), crs=st_crs(au) ),
        polys = au[, "AUID"],
        varname="AUID"
      )

      if (is.na(nbfocal)) {
        message( "no data" )
        return (NULL)
      }
      NB_graph = attr(au, "NB_graph")  # full matrix
      nbnames = attr( NB_graph, "region.id")
      nnAUID = nbnames[NB_graph[[which( nbnames == nbfocal )]]]  ## nearest neighbours
      tokeep = unique( c( nnAUID, nbfocal ) )

      if ( stmv_au_buffer_links > 0 ) {
        # no of additional neighbourhood links ... 0 == nearest neighbours, 1 == nn + next nearest neighbours, etc
        for (i in 1:stmv_au_buffer_links ) {
          new = NULL
          for ( foc in nnAUID ) {
            new = c( new, nbnames[NB_graph[[which( nbnames == foc  )]]] ) ## nearest neighbours
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
      #       NB_graph = attr(au, "NB_graph")  # full matrix
      # nbnames = attr( NB_graph, "region.id")

      nnAUID = nbnames[NB_graph[[which( nbnames == nbfocal )]]]  ## nearest neighbours
      tokeep = unique( c( nnAUID, nbfocal ) )

      NB_graph.remove = which( ! (nbnames %in% sppoly$AUID ) )
      if ( length(NB_graph.remove) > 0 ) {
        # remove isolated locations and recreate sppoly .. alternatively add links to NB_graph
        NB_graph = nb_remove( NB_graph, NB_graph.remove )
        sppoly = sp::spChFIDs( sppoly, row.names(sppoly) )  #fix id's
        # sppoly = sppoly[ order(sppoly$AUID), ]
      }

      attr(sppoly, "NB_graph") = NB_graph  # adding neighbourhood as an attribute to sppoly
      NB_graph =NULL
      NB_graph.remove =NULL

      attr( pa, "sppoly" ) = sppoly

      rownames(pa) = NULL
      return(pa)

  }

  if (0) {
      jj = which( card(NB_graph) == 0)
      jj = match( tokeep, au$AUID )
      plot(sppoly)
      plot(sppoly[jj,], add=T, col="red")
      dev.new()
      edit(NB_graph, polys=sppoly)
      card(NB_graph) # last check if any more  isolated areas
  }

}
