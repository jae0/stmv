

  polygon.area.lon.lat = function( p, ... ) {
    # call as: 
    # r = sapply( a, polygon.area.lon.lat, simplify=T )
    # where a is a list of polygons eg:
    # a = Map2poly(ns)
    
    p = as.data.frame( p )
    names(p) = c("lon", "lat") 
    p = p[ which(is.finite( rowSums(p) ) ) ,]
    p = rbind( p , p[ nrow(p) , ] ) # pad the ending
    p = lonlat2planar( p, ... )
    res = polygon.area( x=p$plon, y=p$plat )
    
    return(res)
  }


