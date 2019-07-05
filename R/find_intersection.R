find_intersection = function( Y, threshold=0 ) {
    # find first to cross a threshold
    good = which(is.finite(Y))
    ngood = length(good)
    X = 1:ngood
    i = 1
    for (i in 1:ngood) {
      if ( Y[good[i]] > threshold ) next()
    }
    if (i > ngood) i = ngood
    index = good[i]
    return(index)
}
