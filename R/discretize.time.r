  discretize.time = function( time.stamp, ti.scale, ti.offset=NULL, ti.min=NULL  ) {
   
    if (!is.null(ti.offset)) time.stamp = time.stamp + lubridate::dyears(ti.offset) 
    
    if ( ti.scale < 1 ) {
      brks = seq( 0, 1, by=p$ti.scale )
      center = brks + p$ti.scale / 2
      center = center[ - length(center) ]
      ti = cut( lubridate::yday( time.stamp ) / 366, breaks=brks, labels=center , include.lowest=TRUE ) 
      ti = as.numeric(levels( ti ))[ as.integer( ti )]
      out = lubridate::year( time.stamp ) + ti 
    }
    if (ti.scale == 1 ) {
      out = year( time.stamp ) 
    }
    if (ti.scale > 1 ) {
      if (is.null( ti.min)) ti.min = lubridate::year(min(time.stamp) )
      brks = seq( ti.min -p$ti.scale /2 , lubridate::year( max(time.stamp) ) + p$ti.scale /2 , by=p$ti.scale )
      center = brks + p$ti.scale / 2
      center = center[ - length(center) ]
      ti  = cut( lubridate::year(time.stamp) , breaks=brks, labels=center, include.lowest=TRUE  ) 
      out  = levels( ti )[ as.integer( ti )]
    }
    return( as.numeric(out) )
  }

