stmv_predictionarea_space = function(p, sloc, windowsize.half) {
  # lattice form
  # system size
  # nr = nx # nr .. x/plon  == pa_coord_names[1] ;;
  # nc = ny # nc .. y/plat  == pa_coord_names[2]
  dx = p$pres
  dy = p$pres
  pa_coord_names = p$stmv_variables$LOCS[1:2]

  pa_w = -windowsize.half : windowsize.half # default window size
  pa_w_n = length(pa_w)

  # determine prediction locations and time slices
  iwplon = aegis_floor( (sloc[1]-p$origin[1]) /dx + 1L) + pa_w
  iwplat = aegis_floor( (sloc[2]-p$origin[2]) /dy + 1L) + pa_w

  pa = data.table( iplon = rep.int(iwplon, pa_w_n) ,
                   iplat = rep.int(iwplat, rep.int(pa_w_n, pa_w_n)) )

  tokeep = which( pa$iplon >= 1 & pa$iplon <= p$nplons & pa$iplat >= 1 & pa$iplat <= p$nplats )
  if (length(tokeep) < 1 ) return(NULL)
  pa = pa[tokeep,]

  Ploc = stmv_attach( p$storage_backend, p$ptr$Ploc )
  ploc_ids = array_map( "xy->1", Ploc[], gridparams=p$gridparams )
  pa_ids = array_map( "2->1", pa[, c("iplon", "iplat")], gridparams=p$gridparams )

  pa$i = match( pa_ids, ploc_ids )

  ploc_ids = NULL
  pa_ids = NULL

  tokeep = which( is.finite(pa$i) )
  if (length(tokeep) < 1 ) return(NULL)
  pa = pa[tokeep,]

  pa[, pa_coord_names[1]] = Ploc[ pa$i, 1 ]
  pa[, pa_coord_names[2]] = Ploc[ pa$i, 2 ]

  # prediction covariates i.e., independent stmv_variables/ covariates
  pvars = c(pa_coord_names, "i")

  if (p$nloccov > 0) {
    # .. not necessary except when covars are modelled locally
    for (ci in 1:p$nloccov) {
      vn = p$stmv_variables$local_cov[ci]
      pu = NULL
      pu = stmv_attach( p$storage_backend, p$ptr$Pcov[[vn]] )
      nts = ncol(pu)
      if ( nts== 1 ) {
        pvars = c( pvars, vn )
        pa[,vn] = pu[pa$i]  # i.e., a static variable
      }
    }
  }

  pa = as.data.frame(pa[, ..pvars])

  return(pa)

}