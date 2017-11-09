
  # ----------------------------------------------------------------------------
  change.resolution = function (x, res=10) {
    x = convert.degdec2degmin(x)
    x$lon = trunc(x$lon*(100/res)) / (100/res)
    x$lat = trunc(x$lat*(100/res)) / (100/res)
    x = convert.degmin2degdec(x)
    return(x)
  }


