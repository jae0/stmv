stmv_discretize_coordinates = function(coords, z=NULL, discretized_n=100, minresolution=NULL, ti.offset=NULL, ti.min=NULL, FUNC=mean, method="default", ...) {

  ncoords = ncol(coords)
  ndata =   nrow(coords)

  if (ncoords==1) {
    # assume it is time
    # not used/finished yet
    if (!is.null(ti.offset)) coords = coords + lubridate::dyears(ti.offset)
    if ( minresolution < 1 ) {
      brks = seq( 0, 1, by=minresolution )
      center = brks + minresolution / 2
      center = center[ - length(center) ]
      ti = cut( lubridate::yday( coords ) / 366, breaks=brks, labels=center , include.lowest=TRUE )
      ti = as.numeric(levels( ti ))[ as.integer( ti )]
      out = lubridate::year( coords ) + ti
    }
    if (minresolution == 1 ) {
      out = year( coords )
    }
    if (minresolution > 1 ) {
      if (is.null( ti.min)) ti.min = lubridate::year(min(coords) )
      brks = seq( ti.min -minresolution /2 , lubridate::year( max(coords) ) + minresolution /2 , by=minresolution )
      center = brks + minresolution / 2
      center = center[ - length(center) ]
      ti  = cut( lubridate::year(coords) , breaks=brks, labels=center, include.lowest=TRUE  )
      out  = levels( ti )[ as.integer( ti )]
    }
    return( as.numeric(out) )
  }


  if (ncoords==2) {
    # 2 dims assuming spatial only
    drange = NA
    if (is.null(minresolution)){
      yr = range( coords[,2], na.rm=TRUE )
      xr = range( coords[,1], na.rm=TRUE )
      drange = sqrt( diff( xr)^2 + diff( yr)^2  )
      dmin = min( diff(xr),  diff(yr) )
      minresolution = rep( dmin / discretized_n, 2)
    }

    coords[,1] = floor( coords[,1] / minresolution[1] ) * minresolution[1]
    coords[,2] = floor( coords[,2] / minresolution[2] ) * minresolution[2]

    if ( method=="default" ) {
      # basic discretization in 2D
      return(coords)
    }

    if ( method=="aggregate" ) {
      # basic aggregation in 2D
      if (is.null(z)) stop("z is required if aggregating")
      res = tapply( X=z, INDEX=list(coords[,1], coords[,2]),
        FUN = function(w) {FUNC(w, ...)}, simplify=TRUE )
      res = as.data.frame( as.table (res) )
      res[,1] = as.numeric(as.character( res[,1] ))
      res[,2] = as.numeric(as.character( res[,2] ))
      res = res[ which( is.finite( res[,3] )) ,]
      names(res) =c("x", "y", "z")
      return( res)
    }


    if (method=="thin") {
      if (is.null(minresolution)) stop("minresolution is required")
      res = tapply( X=ndata, INDEX=list(coords[,1], coords[,2]),
        FUN = function(w) {length( which(is.finite(w) ))}, simplify=TRUE )
      res = as.data.frame( as.table (res) )
      res[,1] = as.numeric(as.character( res[,1] ))
      res[,2] = as.numeric(as.character( res[,2] ))
      res = res[ which( is.finite( res[,3] )) ,]
      names(res) =c("x", "y", "n")
      ressum = sum(res$n, na.rm=TRUE)
      ntoremove = ressum - ntarget
      invcount = 1/res$n * ntarget / ressum  # proportion to remove to make each cell equal in weight
      keep = NULL
      for (o in 1:nrow(res)) {
        oo = which( coords[,1] == res[o,1] & coords[,2] == res[o,2] )
        noo = length(oo)
        if ( noo > 1) {
          okeep = floor(invcount[o] * noo)
          if (okeep > 1) {
            ss = oo[ .Internal( sample( noo, okeep, replace=FALSE, prob=NULL)) ]
            keep = c(keep, ndata[ss] )
          }
        } else {
          keep = c(keep, ndata[oo] )
        }
      }
      return( keep )
    }

  }




  if (ncoords==3) {
      # assume third dim is time .. find an more elegant solution when required
      drange = NA
      if (is.null(minresolution)){
        xr = range( coords[,1], na.rm=TRUE )
        yr = range( coords[,2], na.rm=TRUE )
        tr = range( coords[,3], na.rm=TRUE )
        drange = sqrt( diff( xr)^2 + diff( yr)^2  )
        dmin = min( diff(xr),  diff(yr) )
        tmin = diff(tr)
        minresolution = c( dmin / discretized_n, dmin / discretized_n, tmin/discretized_n )
      }

      coords[,1] = floor( coords[,1] / minresolution[1] ) * minresolution[1]
      coords[,2] = floor( coords[,2] / minresolution[2] ) * minresolution[2]
      coords[,3] = floor( coords[,3] / minresolution[3] ) * minresolution[3]

      if ( method=="default" ) {
        # basic discretization in 2D
        return(coords)
      }

      if ( method=="aggregate" ) {
        # basic aggregation in 2D
        if (is.null(z)) stop("z is required if aggregating")
        res = tapply( X=z, INDEX=list(coords[,1], coords[,2], coords[,3]),
          FUN = function(w) {FUNC(w, ...)}, simplify=TRUE )
        res = as.data.frame( as.table (res) )
        res[,1] = as.numeric(as.character( res[,1] ))
        res[,2] = as.numeric(as.character( res[,2] ))
        res[,3] = as.numeric(as.character( res[,3] ))
        res = res[ which( is.finite( res[,4] )) ,]
        names(res) =c("x", "y", "t", "z")
        return( res )
      }


      if (method=="thin") {
        res = tapply( X=ndata, INDEX=list(coords[,1], coords[,2], coords[,3]),
          FUN = function(w) {length( which(is.finite(w) ))}, simplify=TRUE )
        res = as.data.frame( as.table (res) )
        res[,1] = as.numeric(as.character( res[,1] ))
        res[,2] = as.numeric(as.character( res[,2] ))
        res[,3] = as.numeric(as.character( res[,3] ))
        res = res[ which( is.finite( res[,4] )) ,]
        names(res) =c("x", "y", "t", "n")
        ressum = sum(res$n, na.rm=TRUE)
        ntoremove = ressum - ntarget
        invcount = 1/res$n * ntarget / ressum  # proportion to remove to make each cell equal in weight
        keep = NULL
        for (o in 1:nrow(res)) {
          oo = which( coords[,1] == res[o,1] & coords[,2] == res[o,2] & coords[,3] == res[o,3] )
          noo = length(oo)
          if ( noo > 1) {
            okeep = floor(invcount[o] * noo)
            if (okeep > 1) {
              ss = oo[ .Internal( sample( noo, okeep, replace=FALSE, prob=NULL)) ]
              keep = c(keep, ndata[ss] )
            }
          } else {
            keep = c(keep, ndata[oo] )
          }
        }
        return( keep )
      }
    }

  }
