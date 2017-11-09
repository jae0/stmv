
  polygon.area = function(x, y) {
    area = 0
    n <- length(x)
    i2 <- c(n, 1:(n - 1))
    x2 <- x[i2]
    y2 <- y[i2]
    res = 0.5 * abs(sum(x * y2 - x2 * y))
    return( res )
  }


