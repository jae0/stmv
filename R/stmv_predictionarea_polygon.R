  stmv_predictionarea_polygon = function( p, pa, dx, dy, pa_proj4string_planar_km, pa_coord_names=c("plon", "plat"), global_sppoly ) {


    if (is.null(global_sppoly)) {

      # system size
      # nr = nx # nr .. x/plon ;; nc = ny # nc .. y/plat
      # dx=p$pres, dy=p$pres, pa_coord_names=p$stmv_variables$LOCS[1:2]
      x_r = range(pa[, pa_coord_names[1]])
      x_c = range(pa[, pa_coord_names[2]])

      rr = diff(x_r)
      rc = diff(x_c)

      nr = floor( rr/dx ) + 1L
      nc = floor( rc/dy ) + 1L

      # check:
      #dr = rr/(nr-1) # == dx  ;; dc = rc/(nc-1) # == dy

      # default behaviour .. lattice grid
        # data_subset = NULL
      sppoly = sf::st_as_sf( pa, coords=pa_coord_names )
      sf::st_crs(sppoly) = sf::st_crs( pa_proj4string_planar_km )
      st_geometry(sppoly) = st_geometry( sf::st_make_grid( sppoly, cellsize=c(dx, dy), n=c(nr, nc)) )
      sppoly = as(sppoly, "Spatial")
      if (exists("i", slot(sppoly, "data")) ) {
        sppoly$AUID = sppoly$i
      } else {
        sppoly$AUID = 1:nrow(sppoly)
      }
      sppoly$AUID = as.character( sppoly$AUID )
      # poly* function operate on Spatial* data
      W.nb = poly2nb(sppoly, row.names=sppoly$AUID, queen=TRUE)  # slow .. ~1hr?
      W.remove = which(card(W.nb) == 0)

      if ( length(W.remove) > 0 ) {
        # remove isolated locations and recreate sppoly .. alternatively add links to W.nb
        W.keep = which(card(W.nb) > 0)
        W.nb = nb_remove( W.nb, W.remove )
        sppoly = sppoly[W.keep,]
        row.names(sppoly) = sppoly$AUID
        sppoly = sp::spChFIDs( sppoly, row.names(sppoly) )  #fix id's
        # sppoly = sppoly[ order(sppoly$AUID), ]
      }

      attr(sppoly, "nb") = W.nb  # adding neighbourhood as an attribute to sppoly
      W.nb =NULL
      W.remove =NULL

      return(sppoly)

    }  else if (is.character(global_sppoly)) {

      if (global_sppoly="fast_dummy") {
        # a large fast polygon for fast variable testing
        nr = 2
        nc = 2
        # global sppoly exists .. partition graph
        sppoly = sf::st_as_sf( pa, coords=pa_coord_names )
        sf::st_crs(sppoly) = sf::st_crs( pa_proj4string_planar_km )
        st_geometry(sppoly) = st_geometry( sf::st_make_grid( sppoly,  n=c(nr, nc)) )
        sppoly = as(sppoly, "Spatial")
        if (exists("i", slot(sppoly, "data")) ) {
          sppoly$AUID = sppoly$i
        } else {
          sppoly$AUID = 1:nrow(sppoly)
        }
        sppoly$AUID = as.character( sppoly$AUID )
        # poly* function operate on Spatial* data
        W.nb = poly2nb(sppoly, row.names=sppoly$AUID, queen=TRUE)  # slow .. ~1hr?
        W.remove = which(card(W.nb) == 0)
        return(sppoly)
      }

    } else if (class(global_sppoly) == "sf") {


    }

}
