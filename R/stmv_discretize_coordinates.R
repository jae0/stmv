stmv_discretize_coordinates = function(coo, z=NULL, discretized_n=100, ntarget=NULL, minresolution=NULL, coord_is_time=FALSE, ti.offset=NULL, ti.min=NULL, FUNC=mean, method="default", ...) {

  if (is.vector(coo)) {
    ncoo = 1
    icoo = 1:length(coo)
  } else {
    ncoo = ncol(coo)
    icoo = 1:nrow(coo)
  }

  ntarget = min( ntarget, length(icoo))

  if (ncoo==1) {

    if (coord_is_time) {
      # not used  yet  and not tested
      if (!is.null(ti.offset)) coo = coo + lubridate::dyears(ti.offset)
      if ( minresolution < 1 ) {
        brks = seq( 0, 1, by=minresolution )
        center = brks + minresolution / 2
        center = center[ - length(center) ]
        ti = cut( lubridate::yday( coo ) / 366, breaks=brks, labels=center , include.lowest=TRUE )
        ti = as.numeric(levels( ti ))[ as.integer( ti )]
        out = lubridate::year( coo ) + ti
      }
      if (minresolution == 1 ) {
        out = year( coo )
      }
      if (minresolution > 1 ) {
        if (is.null( ti.min)) ti.min = lubridate::year(min(coo) )
        brks = seq( ti.min -minresolution /2 , lubridate::year( max(coo) ) + minresolution /2 , by=minresolution )
        center = brks + minresolution / 2
        center = center[ - length(center) ]
        ti  = cut( lubridate::year(coo) , breaks=brks, labels=center, include.lowest=TRUE  )
        out  = levels( ti )[ as.integer( ti )]
      }
      return( as.numeric(out) )

    } else {

      if (is.null(minresolution)){
        xr = range( coo[,1], na.rm=TRUE )
        dmin = diff(xr)
        minresolution = diff(xr) / discretized_n
      }

      coo = floor( coo  / minresolution + 1) * minresolution

      if ( method=="default" ) {
        # basic discretization in 1D
        return(coo)
      }

      if ( method=="aggregate" ) {
        # basic aggregation in 2D
        if (is.null(z)) stop("z is required if aggregating")
        res = tapply( X=z, INDEX=coo, FUN=function(w) {FUNC(w, ...)}, simplify=TRUE )
        res = as.data.frame( as.table (res) )
        res[,1] = as.numeric(as.character( res[,1] ))
        res = res[ which( is.finite( res[,2] )) ,]
        names(res) =c("x", "z")
        return( res)
      }


      if (method=="thin") {
        if (is.null(minresolution)) stop("minresolution is required")
        res = tapply( X=icoo, INDEX=coo,  FUN=function(w) {length( which(is.finite(w) ))}, simplify=TRUE )
        res = as.data.frame( as.table (res) )
        res[,1] = as.numeric(as.character( res[,1] ))
        res = res[ which( is.finite( res[,2] )) ,]
        names(res) =c("x", "n")
        ressum = sum(res$n, na.rm=TRUE)

        invcount = 1/res$n * ntarget / ressum  # proportion to remove to make each cell equal in weight
        tokeep = NULL
        for (o in 1:nrow(res)) {
          oo = which( coo == res[o,1]  )
          noo = length(oo)
          if ( noo > 1) {
            okeep = max(1, floor(invcount[o] * noo))
            if (okeep >= 1) {
              ss = oo[ .Internal( sample( noo, okeep, replace=FALSE, prob=NULL)) ]
              tokeep = c(tokeep, icoo[ss] )
            }
          } else {
            tokeep = c(tokeep, icoo[oo] )
          }
        }
        nk = length(tokeep)
        if (nk > ntarget) tokeep = tokeep[ .Internal( sample( nk, ntarget, replace=FALSE, prob=NULL)) ]

        return( tokeep )
      }

    }
  }


  if (ncoo==2) {
    # 2 dims assuming spatial only
    if (is.null(minresolution)){
      xr = range( coo[,1], na.rm=TRUE )
      yr = range( coo[,2], na.rm=TRUE )
      dmin = min( diff(xr),  diff(yr) )
      minresolution = rep( dmin / discretized_n, 2)
    }

    coo[,1] = floor( coo[,1] / minresolution[1] + 1) * minresolution[1]
    coo[,2] = floor( coo[,2] / minresolution[2] + 1) * minresolution[2]

    if ( method=="default" ) {
      # basic discretization in 2D
      return(coo)
    }

    if ( method=="aggregate" ) {
      # basic aggregation in 2D
      if (is.null(z)) stop("z is required if aggregating")
      res = tapply( X=z, INDEX=list(coo[,1], coo[,2]),
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
      res = tapply( X=icoo, INDEX=list(coo[,1], coo[,2]),
        FUN = function(w) {length( which(is.finite(w) ))}, simplify=TRUE )
      res = as.data.frame( as.table (res) )
      res[,1] = as.numeric(as.character( res[,1] ))
      res[,2] = as.numeric(as.character( res[,2] ))
      res = res[ which( is.finite( res[,3] )) ,]
      names(res) =c("x", "y", "n")
      ressum = sum(res$n, na.rm=TRUE)

      invcount = 1/res$n * ntarget / ressum  # proportion to remove to make each cell equal in weight
      tokeep = NULL
      for (o in 1:nrow(res)) {
        oo = which( coo[,1] == res[o,1] & coo[,2] == res[o,2] )
        noo = length(oo)
        if ( noo > 1) {
          okeep = max(1, floor(invcount[o] * noo))
          if (okeep >= 1) {
            ss = oo[ .Internal( sample( noo, okeep, replace=FALSE, prob=NULL)) ]
            tokeep = c(tokeep, icoo[ss] )
          }
        } else {
          tokeep = c(tokeep, icoo[oo] )
        }
      }
      nk = length(tokeep)
      if (nk > ntarget) tokeep = tokeep[ .Internal( sample( nk, ntarget, replace=FALSE, prob=NULL)) ]

      return( tokeep )
    }

  }


  if (ncoo==3) {
      # assume third dim is time .. find an more elegant solution when required
      if (is.null(minresolution)){
        xr = range( coo[,1], na.rm=TRUE )
        yr = range( coo[,2], na.rm=TRUE )
        tr = range( coo[,3], na.rm=TRUE )
        dmin = min( diff(xr),  diff(yr) )
        tmin = diff(tr)
        minresolution = c( dmin / discretized_n, dmin / discretized_n, tmin/discretized_n )
      }

      coo[,1] = floor( coo[,1] / minresolution[1] +1 ) * minresolution[1]
      coo[,2] = floor( coo[,2] / minresolution[2] +1 ) * minresolution[2]
      coo[,3] = floor( coo[,3] / minresolution[3] +1 ) * minresolution[3]

      if ( method=="default" ) {
        # basic discretization in 2D
        return(coo)
      }

      if ( method=="aggregate" ) {
        # basic aggregation in 2D
        if (is.null(z)) stop("z is required if aggregating")
        res = tapply( X=z, INDEX=list(coo[,1], coo[,2], coo[,3]),
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
        res = tapply( X=icoo, INDEX=list(coo[,1], coo[,2], coo[,3]),
          FUN = function(w) {length( which(is.finite(w) ))}, simplify=TRUE )
        res = as.data.frame( as.table (res) )
        res[,1] = as.numeric(as.character( res[,1] ))
        res[,2] = as.numeric(as.character( res[,2] ))
        res[,3] = as.numeric(as.character( res[,3] ))
        res = res[ which( is.finite( res[,4] )) ,]
        names(res) =c("x", "y", "t", "n")
        ressum = sum(res$n, na.rm=TRUE)

        invcount = 1/res$n * ntarget / ressum  # proportion to remove to make each cell equal in weight
        tokeep = NULL
        for (o in 1:nrow(res)) {
          oo = which( coo[,1] == res[o,1] & coo[,2] == res[o,2] & coo[,3] == res[o,3] )
          noo = length(oo)
          if ( noo > 1) {
            okeep = max(1, floor(invcount[o] * noo))
            if (okeep >= 1) {
              ss = oo[ .Internal( sample( noo, okeep, replace=FALSE, prob=NULL)) ]
              tokeep = c(tokeep, icoo[ss] )
            }
          } else {
            tokeep = c(tokeep, icoo[oo] )
          }
        }
        nk = length(tokeep)
        if (nk > ntarget) tokeep = tokeep[ .Internal( sample( nk, ntarget, replace=FALSE, prob=NULL)) ]

        return( tokeep )
      }
    }

  }
