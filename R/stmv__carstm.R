
stmv__carstm = function( p=NULL, dat=NULL, pa=NULL, improve.hyperparam.estimates=FALSE, ... ) {


    if (0) {
      varObs = S[Si, i_sdObs]^2
      varSpatial = S[Si, i_sdSpatial]^2
      sloc = Sloc[Si,]
      eps = 1e-9
      improve.hyperparam.estimates=FALSE
      variablelist=FALSE
    }


    params = list(...)

    sdTotal = sd( dat[[p$stmv_variables$Y]], na.rm=T)

    # system size
    #nr = nx # nr .. x/plon
    #nc = ny # nc .. y/plat

    vns = p$stmv_variables$LOCS  # for data.table., this has to be a simple variable

    sppoly = attr( pa, "sppoly" )
    dat$AUID = st_points_in_polygons(
      pts = st_transform( st_as_sf( as.matrix(dat[, ..vns ]), coords=vns, crs=st_crs(p$aegis_proj4string_planar_km) ), crs=st_crs(sppoly) ),
      polys = sppoly[, "AUID"],
      varname="AUID"
    )

    dat = dat[  !is.na(dat$AUID)]
    dat$tag = "observations"

    sppoly_df = as.data.table(sppoly)

    pa = data.table(pa)
    pa = pa[ sppoly_df[, .(AUID, i)],  on="i"]  # merge of AUID is required if there is time

    pa[, p$stmv_variables$Y] = NA
    pa$tag ="predictions"

    vn = setdiff( names( pa ), c("i", p$stmv_variables$LOCS  ) )

    dat = rbind( dat[, ..vn], pa[, ..vn] )

    sppoly_df = NULL

    dat$aui  = as.numeric( factor(dat$AUID) )  # for inla

    # get hyper param scalings
    # function copied from carstm .. used below

    j = which( is.finite(dat[[ p$stmv_variables$Y ]]) )
    m = dat[[p$stmv_variables$Y ]] [j]

    H = stmv_carstm_hyperparameters( sd(m), alpha=0.5, median(m) )
    m = NULL
    # # adjust based upon RAM requirements and ncores
    # inla.setOption(num.threads= p$inla_num.threads)
    # inla.setOption(blas.num.threads=p$inla_blas.num.threads)

    gc()

    fit  = NULL

    control.inla.variations = list(
      list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive") ,
      list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive", h=0.05, cmin=0, tolerance=1e-9),
      list( optimise.strategy="smart", stupid.search=FALSE, strategy="adaptive", h=0.1,  cmin=0, tolerance=1e-6),
      list( optimise.strategy="smart", stupid.search=FALSE, strategy="laplace", fast=FALSE, step.factor=0.1),
      list( optimise.strategy="smart", stupid.search=TRUE, h=0.1 ),
      list( optimise.strategy="smart", stupid.search=TRUE, h=0.01 ),
      list( optimise.strategy="smart", stupid.search=TRUE, h=0.001 )
    )

    for (ll in 1:length(control.inla.variations)) {
      inlacall = gsub(
        "verbose=TRUE",
        paste( "control.inla=", deparse( control.inla.variations[[ll]] , width.cutoff=500), ", verbose=TRUE ", sep=" " ),
        p$stmv_local_modelcall
      )
      assign("fit", eval(parse(text=paste( "try(", inlacall, ")" ) ) ))
      if ( !is.null(fit)) if (! inherits(fit, "try-error") ) break()
    }

    if ( is.null(fit) )  return(NULL)
    if ( inherits(fit, "try-error") ) return(NULL)

    # to improve hyper param estimates..
    if (improve.hyperparam.estimates) fit = inla.hyperpar(fit, dz=0.25, diff.logdens=18 )  # get improved estimates for the hyperparameters

    # results container
    # initial prediction container
    hyp_pars =  c( fit$summary.hyperpar$mean, fit$summary.hyperpar$sd )
    names( hyp_pars ) =  gsub( " ", "_",  c( rownames(fit$summary.hyperpar), paste( rownames(fit$summary.hyperpar), "_sd", sep="" ) ) )
    res = list( dat=dat, dimensionality = p$aegis_dimensionality )
    res$predictions = data.table( i=pa$i)
    pa = NULL

    vnm = p$stmv_variables$Y
    oo = cbind( dat[ tag=="observations", ..vnm], fit$summary.fitted.values[ which(dat$tag=="observations"), "mean" ] )
    names(oo) = c("obs", "preds")
    oores = try( lm( preds~obs, oo, , na.action=na.omit) )
    rsquared = ifelse( inherits(oores, "try-error"), NA,  summary(oores)$r.squared )
    oo = oores = NULL

    mlik_p_eff = c(fit$neffp)
    names(mlik_p_eff) = gsub( " ", "_", rownames(fit$neffp) )

    res$stmv_stats = c(
      sdTotal=sdTotal,
      ndata=nrow(dat),
      rsquared =rsquared,
      fixed_mean = fit$summary.fixed$mean,
      fixed_sd = fit$summary.fixed$sd,
      dic = fit$dic$dic,
      dic_p_eff = fit$dic$p.eff,
      waic = fit$waic$waic, #1
      waic_p_eff = fit$waic$p.eff, #1
      mlik = fit$mlik[[2,1]], #1
      mlik_p_eff,
      hyp_pars
    )



    # row indices for predictions
    if ( p$aegis_dimensionality == "space") {
      AUID = sppoly[["AUID"]]
      res$i_preds = which(
        dat$tag=="predictions" &
        dat$AUID %in% AUID
      )  # filter by AUID and years in case additional data in other areas and times are used in the input data
      res$AUID = AUID
      res$matchfrom = list( AUID=dat$AUID[res$i_preds] )
      res$matchto   = list( AUID=res$AUID )
    }

    if ( p$aegis_dimensionality == "space-year") {
      AUID = sppoly[["AUID"]]
      year = as.character( p$yrs )
      dat$year = as.character(dat$year)

      res$i_preds = which(
        dat$tag=="predictions" &
        dat$AUID %in% AUID &
        dat$year %in% year
      )  # filter by AUID and years in case additional data in other areas and times are used in the input data

      res$AUID = AUID
      res$year = year
      res$matchfrom = list( AUID=dat$AUID[res$i_preds], year=dat$year[res$i_preds] )
      res$matchto   = list( AUID=res$AUID,   year=res$year )
    }

    if ( p$aegis_dimensionality == "space-year-season") {
      AUID = sppoly[["AUID"]]
      year = as.character( p$yrs )
      dyear = as.character( discretize_data( (p$dyears + diff(p$dyears)[1]/2), p$discretization[["dyear"]] ) )
      dat$year = as.character(dat$year)
      dat$dyear = as.character( discretize_data( dat$dyear, p$discretization[["dyear"]] ) )

      res$i_preds = which(
        dat$tag=="predictions" &
        dat$AUID %in% AUID &
        dat$year %in% year
      )  # filter by AUID and years in case additional data in other areas and times are used in the input data
      res$AUID = AUID
      res$year = year
      res$dyear = dyear
      res$matchfrom = list( AUID=dat$AUID[res$i_preds], year=dat$year[res$i_preds], dyear=dat$dyear[res$i_preds] )
      res$matchto   = list( AUID=res$AUID,   year=res$year, dyear=res$dyear )
    }


    # finish predictions
      res$predictions$mean = reformat_to_array( input=fit$summary.fitted.values[ res$i_preds, "mean" ], matchfrom=res$matchfrom, matchto=res$matchto )
      res$predictions$sd = reformat_to_array( input=fit$summary.fitted.values[ res$i_preds, "sd" ], matchfrom=res$matchfrom, matchto=res$matchto )




    ## --------- predictions complete ------



  # lattice::levelplot( mean ~ plon + plat, data=pa, col.regions=heat.colors(100), scale=list(draw=TRUE) , aspect="iso" )

  return(res)

  if (0) {
    # skip

    # row indices for spatial effects (if any)
    nAUID = length(res$AUID)

    if (exists("summary.random", fit)) {

      if (exists("aui", fit$summary.random)) {

        if (nrow(fit$summary.random$aui) == nAUID*2) {
          # a single nonspatial effect (no grouping across time)
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial") )
          res$i_nonspatial = which(resout$type=="nonspatial")
          res$ns_matchfrom = list( AUID=resout$AUID[res$i_nonspatial]  )
          res$ns_matchto   = list( AUID=res$AUID  )

        } else if (nrow(fit$summary.random$aui) == nAUID*2 * p$ny ) {
          #  nonspatial effects grouped by year
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial"), year=p$yrs )
          res$i_nonspatial = which(resout$type=="nonspatial")
          res$ns_matchfrom = list( AUID=resout$AUID[res$i_nonspatial], year=resout$year[res$i_nonspatial] )
          res$ns_matchto   = list( AUID=res$AUID,   year=res$year  )
        } else if (nrow(fit$summary.random$aui) == nAUID*2 * p$nt ) {
          # nonspatial at all time slices
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial"), year=p$yrs, dyear=p$dyears )
          res$i_nonspatial = which(resout$type=="nonspatial")
          res$ns_matchfrom = list( AUID=resout$AUID[res$i_nonspatial], year=as.character(resout$year[res$i_nonspatial]), dyear=as.character( discretize_data( resout$dyear[res$i_nonspatial], p$discretization[["dyear"]] ) ) )
          res$ns_matchto   = list( AUID=res$AUID,   year=res$year, dyear=res$dyear )
        }

        if (nrow(fit$summary.random$aui) == nAUID*2) {
          # a single spatial effect (no grouping across time)
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial") )
          res$i_spatial = which(resout$type=="spatial")
          res$sp_matchfrom = list( AUID=resout$AUID[res$i_spatial]  )
          res$sp_matchto   = list( AUID=res$AUID  )

        } else if (nrow(fit$summary.random$aui) == nAUID*2 * p$ny ) {
          # spatial effects grouped by year
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial"), year=p$yrs )
          res$i_spatial = which(resout$type=="spatial")
          res$sp_matchfrom = list( AUID=resout$AUID[res$i_spatial], year=resout$year[res$i_spatial] )
          res$sp_matchto   = list( AUID=res$AUID,   year=res$year  )
        } else if (nrow(fit$summary.random$aui) == nAUID*2 * p$nt ) {
          # at every time slice
          resout = expand.grid( AUID=res$AUID, type = c("nonspatial", "spatial"), year=p$yrs, dyear=p$dyears )
          res$i_spatial = which(resout$type=="spatial")
          res$sp_matchfrom = list( AUID=resout$AUID[res$i_spatial], year=as.character(resout$year[res$i_spatial]), dyear=as.character( discretize_data( resout$dyear[res$i_spatial], p$discretization[["dyear"]] ) ) )
          res$sp_matchto   = list( AUID=res$AUID,   year=res$year, dyear=res$dyear )
        }

      }
    }

    ## --------- start random effects -------

    # match conditions for random effects .. i_preds are locations of predictions in "fit"
    # random effects results ..

    if (exists("summary.random", fit)) {

      if (exists("iid_error", fit$summary.random)) {
        # IID random effects
        vn = paste( p$stmv_variables$Y, "random_sample_iid", sep=".")
        input = fit$summary.random$iid_error[res$i_preds, "mean" ]
        res[[vn]] = reformat_to_array( input=input, matchfrom=res$matchfrom, matchto=res$matchto )
        if (!is.null(NA_mask)) res[[vn]][NA_mask] = NA
      }

      if (exists("aui", fit$summary.random)) {

        input = fit$summary.random$aui[ res$i_nonspatial, "mean" ]
        vn = paste( p$stmv_variables$Y, "random_auid_nonspatial", sep=".")
        res[[vn]] = reformat_to_array( input=input, matchfrom=res$ns_matchfrom, matchto=res$ns_matchto )
        if (!is.null(NA_mask)) res[[vn]][NA_mask] = NA
        # carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2000", dyear="0.8" ) )

        input = fit$summary.random$aui[ res$i_spatial, "mean" ]  # offset structure due to bym2
        vn = paste( p$stmv_variables$Y, "random_auid_spatial", sep=".")
        res[[vn]] = reformat_to_array( input=input, matchfrom=res$sp_matchfrom, matchto=res$sp_matchto )
        if (!is.null(NA_mask)) res[[vn]][NA_mask] = NA
        # carstm_plot( p=p, res=res, vn=vn, time_match=list(year="2000", dyear="0.8" ) )

      }
    }

  }

}