find_intersection = function( Y, threshold=0 ) {
    # find first to cross a threshold
    good = which(is.finite(Y))
    ngood = length(good)
    X = 1:ngood
    i = 1
    while( Y[good[i]] > threshold ) i=i+1
    if (i > ngood) i=ngood
    index = good[i-1]
    return(index)
}
