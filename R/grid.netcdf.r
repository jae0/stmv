
  grid.netcdf = function(ncdf.fname) {
    require(ncvar) # to handle netcdf files (e.g., binary output format of GMT)
    z = var.get.ncv(path=ncdf.fname, name="z", mode="nometa")
    dims = var.get.ncv(path=ncdf.fname, name="dim")
    out = as.matrix( z$data )
    # out = out[,ncol(out):1]
    return (out)
  }


