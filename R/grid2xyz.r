
  grid2xyz = function (M, rows, cols, method="R") {

    if (method=="brute.force") {
      # convert matrix M to xyz coords using rows and cols
      rows = as.vector(rows, "numeric")
      cols = as.vector(cols, "numeric")
      n = 0
      y = f1 = f2 = vector("numeric", dim(M)[1] * dim(M)[2])
      for (j in 1:dim(M)[2]) {
          for (i in 1:dim(M)[1]) {
              q = M[i, j]
              if (is.finite(q) && q!=0) {
                  n = n + 1
                  y[n] = q
                  f1[n] = rows[i]
                  f2[n] = cols[j]
              }
          }
      }
      y = y[1:n]
      f1 = f1[1:n]
      f2 = f2[1:n]
      D = as.data.frame(cbind(lon=f1, lat=f2, z=y))
    }

    if (method=="R") {
      dimnames(M) = list(rows, cols)
      D = as.data.frame.table(M, optional = TRUE)
      names(D) = c("lon", "lat", "z")
      for (i in c(1,2)) D[,i] = as.numeric(as.character(D[,i]))
     }

    return(D)
  }



