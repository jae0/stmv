
# NOTE:: not finishe porting fully ... only designed for xy data right no .. no time .. needs more testing

stmv__inla_spde = function( p=NULL, dat=NULL, pa=NULL, plotmesh=FALSE, ... ) {
  #\\ generic spatial and space-time interpolator using inla
  #\\ parameter and data requirements can be seen in bathymetry\src\bathymetry.r

  sdTotal=sd(dat[,p$stmv_variables$Y], na.rm=T)

  # the following parameters are for inside and outside ... do not make them exact multiples as this seems to make things hang ..
  if ( !exists("inla.mesh.max.edge", p))  p$inla.mesh.max.edge = c(  0.025,   0.04 )    # proportion of 2*p$stmv_distance_scale or equivalent: c(inside,outside) -- must be positive valued
  if ( !exists("inla.mesh.offset", p))  p$inla.mesh.offset   = c( - 0.025,  - 0.05 )   # how much to extend inside and outside of boundary: proportion of stmv_distance_scale .. neg val = proportion
  if ( !exists("inla.mesh.cutoff", p)) p$inla.mesh.cutoff   = c( - 0.05,   - 0.5 )    ## min distance allowed between points: proportion of stmv_distance_scale ; neg val = proportion
  if ( !exists("inla.mesh.hull.radius", p)) p$inla.mesh.hull.radius = c( -0.04, - 0.08 ) ## radius of boundary finding algorythm ; neg val = proportion
  if ( !exists("inla.mesh.hull.resolution", p)) p$inla.mesh.hull.resolution = 125  ## resolution for discretization to find boundary
  if ( !exists("inla.eps", p)) p$inla.eps = p$pres / 10  # add a little noise to coordinates to prevent a race condition
  if ( !exists("inla.alpha", p)) p$inla.alpha = 1.5 # alpha-1 = nu of bessel function curviness
  if ( !exists("inla.nsamples", p)) p$inla.nsamples = 5000 # posterior similations
  if ( !exists("inla.robust.quantiles", p)) p$inla.robust.quantiles = c(0.025, 0.975 )  # posterior predictions robustified by trimming extreme values
  if ( !exists("stmv_local_modelformula", p) )  p$stmv_local_modelformula = formula( z ~ -1 + intercept + f( spatial.field, model=SPDE ) ) # SPDE is the spatial covariance model .. defined in

  # priors
  # kappa0 = sqrt(8)/p$expected.range
  # tau0 = 1/(sqrt(4*pi)*kappa0* p$expected.sigma)

  #-----------------
  # row, col indices
  if ( !exists("stmv.posterior.extract", p)) {
    # example for simplest form
    p$stmv.posterior.extract = function(s, rnm) {
      # rnm are the rownames that will contain info about the indices ..
      # optimally the grep search should only be done once but doing so would
      # make it difficult to implement in a simple structure/manner ...
      # the overhead is minimal relative to the speed of modelling and posterior sampling
      i_intercept = grep("intercept", rnm, fixed=TRUE ) # matching the model index "intercept" above .. etc
      i_spatial.field = grep("spatial.field", rnm, fixed=TRUE )
      return(  s$latent[i_intercept,1] + s$latent[ i_spatial.field,1] )
    }
  }


  # also sending direct distances rather than proportion seems to cause issues..
  locs = dat[,p$stmv_variables$LOC] + runif( nrow(dat)*2, min=-p$pres*p$inla.eps, max=p$pres*p$inla.eps ) # add  noise  to prevent a race condition



  bnd.offset=p$inla.mesh.offset

  cutoff=p$inla.mesh.cutoff

  convex=p$inla.mesh.hull.radius

  resolution=p$inla.mesh.hull.resolution


  #\\ create a mesh object for inla-based analysis
  #\\    note: "-" parameter values means to take as a proportion in inla

  require(INLA)

  lengthscale = max( diff(range( locs[,1])), diff(range( locs[,2]) )) / 10  # in absence of range estimate take 1/10 of domain size
  max.edge = c( 0.25, 1)*lengthscale #   # max size of a triange (in, out) proportion of dist.max
  if (is.null(bnd.offset) ) bnd.offset = -c(  0.025,   0.05 )  # how much to extend inside and outside of boundary: proportion of dist.max
  if (is.null(cutoff) ) cutoff = -c(  0.005, 0.05 ) # min distance allowed between points: proportion of dist.max
  if (is.null(convex) ) convex= c( -0.04, -0.08 )  # convex hull radius
  if (is.null(resolution) ) resolution = 125  # discretization resolution for boundary delineation

  MESH = NULL

  # MESH =try(
  #   inla.mesh.2d (
  #     loc=locs,
  #     max.edge = max.edge,
  #     offset = bnd.offset,
  #     cutoff = cutoff,
  #     boundary = list(
  #       inla.nonconvex.hull(locs, convex=convex[1], resolution=resolution ) ,
  #       inla.nonconvex.hull(locs, convex=convex[2], resolution=resolution ) )
  #   ), silent=TRUE
  # )

  # if ( "try-error" %in% class(MESH)) {
    MESH =try(
      inla.mesh.2d (
        loc=locs,
        max.edge = max.edge ,
        offset = bnd.offset,
        cutoff = cutoff
      ), silent=TRUE
    )
  # }

  if ( "try-error" %in% class(MESH)) {
    MESH =try(
      inla.mesh.2d (
        loc=locs,
        max.edge = max.edge * 1.5 ,  # try a lower resolution
        offset = bnd.offset ,
        cutoff = cutoff
      ), silent=TRUE
    )
  }

  if ( "try-error" %in% class(MESH)) MESH = NULL
  if (!is.null( MESH)  & (plotmesh) ) plot( MESH)

  # the commented options seem to make mesh gen slow ... should look into this if you use this method
  SPDE = inla.spde2.matern( MESH,
    alpha=p$inla.alpha #, # alpha is the Bessel smoothness factor - 1
    # B.tau=matrix(c(log(tau0),-1,+1),nrow=1,ncol=3),
    # B.kappa=matrix(c(log(kappa0),0,-1),nrow=1,ncol=3),
    # theta.prior.mean=c(0,0), # thetas are the internal representations of prior offsets for tau and kappa (i.e.,range and sigma)
    # theta.prior.prec=c(0.1, 0.1) # precision of thetas
  )

  # effects .. a list of two elements: first is for SPDE and second is for covariates
  obs_index = inla.spde.make.index(name="spatial.field", SPDE$n.spde)
  obs_eff = list()
  obs_eff[["spde"]] = c( obs_index, list(intercept=1) )

  obs_A = NULL
  if ( exists( "COV", p$stmv_variables) ) {
    if (length( p$stmv_variables$COV) > 0) {
      covar = as.data.frame( dat[, p$stmv_variables$COV ] )
      colnames( covar ) = p$stmv_variables$COV
      obs_eff[["covar"]] = as.list(covar)
      obs_A = list( inla.spde.make.A( mesh=MESH, loc=as.matrix(dat[, p$stmv_variables$LOC ]) ), 1 )
    }
  }
  if (is.null( obs_A ))  obs_A = list( inla.spde.make.A( mesh=MESH, loc=as.matrix(dat[, p$stmv_variables$LOC ]) ) ) # no effects

  obs_ydata = list()
  obs_ydata[[ p$stmv_variables$Y ]] = dat[, p$stmv_variables$Y ]

  DATA = inla.stack( tag="obs", data=obs_ydata, A=obs_A, effects=obs_eff, remove.unused=FALSE )
  obs_index=obs_eff=obs_ydata=obs_A=NULL

  # remove.unused=FALSE ensures that we can look at the estimated field effect without
  #   having to do expensive separate predictions.
  # DATA$A is projection matrix to translate from mesh nodes to data nodes

  # -----------
  # PREDICTIONS
  #      NOTE:: by default, all areas chosen to predict within the window.. but if covariates are involved,
  #        this can be done only where covariates exists .. so next step is to determine this and predict
  #        over the correct area.
  #        Once all predictions are complete, simple (kernal-based?) interpolation
  #        for areas without covariates can be completed

  # prediction stack:: check for covariates

  preds_locs = as.matrix( pa[, p$stmv_variables$LOC ] )
  preds_index = inla.spde.make.index( name="spatial.field", SPDE$n.spde)

  preds_eff = list()
  preds_eff[["spde"]] = c( preds_index, list(intercept=1) )

  preds_A =NULL
  if ( exists( "COV", p$stmv_variables) ) {
    if (length( p$stmv_variables$COV) > 0) {
      pcovars = as.data.frame(pa[,p$stmv_variables$COV])
      colnames( pcovars ) = p$stmv_variables$COV
      preds_eff[["covar"]] = as.list( pcovars )
      preds_A = list( inla.spde.make.A(MESH, loc=preds_locs ), 1)
      pcovars = NULL
    }
  }
  if (is.null( preds_A) ) preds_A = list( inla.spde.make.A(MESH, loc=preds_locs ) )

  preds_ydata = list()
  preds_ydata[[ p$stmv_variables$Y ]] = NA ## ie. to predict
  PREDS = inla.stack( tag="preds", data=preds_ydata, A=preds_A, effects=preds_eff, remove.unused=FALSE )
  DATA = inla.stack(DATA, PREDS )
  preds_stack_index = inla.stack.index( DATA, "preds")$data  # indices of predictions in stacked data
  preds_eff=preds_ydata=preds_A=PREDS=preds_index=preds_locs=NULL


  RES = NULL

  FAMILY="gaussian"  # family should be gaussian("identity") as it is operating upon the linear link scale

    inputstack = inla.stack.data(DATA)
    A = inla.stack.A(DATA)

    itmpdir=tempfile()

    RES = try( inla( p$stmv_local_modelformula, data = inputstack, family=FAMILY,
          control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config=TRUE), # return linear predictors to compute predictions quickly
          control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
          control.predictor=list (A=A, compute=TRUE, link =1 ), # compute=TRUE on each data location
          quantiles=NULL,  # no quants to save storage requirements
          control.inla = list(h =1e-2), # h=0.01 is default step length for gradient calc of hyper params
          # control.fixed = list(expand.factor.strategy='inla') ,
          working.directory= itmpdir, keep=FALSE,
          verbose=FALSE
      ), silent=TRUE )

    # clean up local files
      fns = list.files( dirname( itmpdir), all.files=TRUE, full.names=TRUE, recursive=TRUE, include.dirs=TRUE )
      oo = grep( basename(itmpdir), fns )
      if(length(oo)>0) file.remove( sort(fns[oo], decreasing=TRUE) )

      if ("try-error" %in% class(RES)) {
        RES = try( inla( p$stmv_local_modelformula, data = inputstack, family=FAMILY,
          control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config = TRUE),
          control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
          control.predictor=list (A=A, compute=TRUE, link =1 ),
          control.inla = list(h = 1e-3, tolerance=1e-8, restart=3), # restart a few times in case posteriors are poorly defined
          working.directory=itmpdir, keep=FALSE,
          verbose=FALSE
        ), silent=TRUE )

        # clean up local files

        fns = list.files( dirname( itmpdir), all.files=TRUE, full.names=TRUE, recursive=TRUE, include.dirs=TRUE )
        oo = grep( basename(itmpdir), fns )
        if(length(oo)>0) file.remove( sort(fns[oo], decreasing=TRUE) )

      }

      if ("try-error" %in% class(RES)) return( NULL ) # give up

      if ( RES$mode$mode.status > 0) {  # make sure Eignevalues of Hessian are appropriate (>0)
        RES = try( inla( p$stmv_local_modelformula, data = inputstack,  family=FAMILY,
          control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config = TRUE),
          control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
          control.predictor=list (A=A, compute=TRUE, link =1 ),
          control.inla = list( h=1e-4, tolerance=1e-10), # increase in case values are too close to zero
          control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
          working.directory=itmpdir, keep=FALSE,
          verbose=FALSE
        ), silent=TRUE )
            # clean up local files

        fns = list.files( dirname( itmpdir), all.files=TRUE, full.names=TRUE, recursive=TRUE, include.dirs=TRUE )
        oo = grep( basename(itmpdir), fns )
        if(length(oo)>0) file.remove( sort(fns[oo], decreasing=TRUE) )

      }

      if ("try-error" %in% class(RES)) return (NULL)

      if ( RES$mode$mode.status > 0) {  # make sure Eigenvalues of Hessian are appropriate (>0)
         RES = try( inla( p$stmv_local_modelformula, data = inputstack,  family=FAMILY,
          control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config = TRUE),
          control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
          control.predictor=list (A=A, compute=TRUE, link =1 ),
          control.inla = list( h=1e-6, tolerance=1e-12), # increase in case values are too close to zero
          control.mode = list( restart=TRUE, result=RES ), # restart from previous estimates
          working.directory=itmpdir, keep=FALSE,
          verbose=FALSE
        ), silent=TRUE )
          # clean up local files

        fns = list.files( dirname( itmpdir), all.files=TRUE, full.names=TRUE, recursive=TRUE, include.dirs=TRUE )
        oo = grep( basename(itmpdir), fns )
        if(length(oo)>0) file.remove( sort(fns[oo], decreasing=TRUE) )

       }

      # if still hessian problems accept the solution .. it should be close enough
      if ("try-error" %in% class(RES)) RES=NULL



  if (is.null(RES))  return(NULL)


  if (0) {
    # inla.spde2.matern creates files to disk that are not cleaned up:
    # No longer an issue? 2016-Nov JC
    spdetmpfn = SPDE$f$spde2.prefix
    fns = list.files( dirname( spdetmpfn ), all.files=TRUE, full.names=TRUE, recursive=TRUE, include.dirs=TRUE )
    oo = grep( basename(spdetmpfn), fns )
    if(length(oo)>0) file.remove( sort(fns[oo], decreasing=TRUE) )
    fns=spdetmpfn=oo=NULL
  }
  DATA=NULL

  # -----------
  # predictions

  preds = NULL
  if (! exists("stmv_inla_prediction", p) ) p$stmv_inla_prediction="direct"

  if ( p$stmv_inla_prediction=="direct" ) {
    # precomputed ... slow and expensive in RAM/CPU, just extract from tag indices
    pa$mean = RES$summary.fitted.values[ preds_stack_index, "mean"]
    pa$sd = RES$summary.fitted.values[ preds_stack_index, "sd"]
  }

  if ( p$stmv_inla_prediction=="projected") {
    #\\ note this method only works with simple additive models
    #\\ when smoothes are involved, it becomes very complicated and direct estimation is probably faster/easier
    pG = inla.mesh.projector( MESH, loc=as.matrix( pa[,p$stmv_variables$LOC]) )
    posterior.samples = inla.posterior.sample(n=p$inla.nsamples, RES)
    rnm = rownames(posterior.samples[[1]]$latent )
    posterior = sapply( posterior.samples, p$stmv.posterior.extract, rnm=rnm )
    posterior.samples=NULL
    # robustify the predictions by trimming extreme values .. will have minimal effect upon mean
    # but variance estimates should be useful/more stable as the tails are sometimes quite long
    for (ii in 1:nrow(posterior )) {
      qnt = quantile( posterior[ii,], probs=p$inla.robust.quantiles, na.rm=TRUE )
      toolow = which( posterior[ii,] < qnt[1] )
      toohigh = which (posterior[ii,] > qnt[2] )
      if (length( toolow) > 0 ) posterior[ii,toolow] = qnt[1]
      if (length( toohigh) > 0 ) posterior[ii,toohigh] = qnt[2]
    }
    pa$mean = c( inla.mesh.project( pG, field=apply( posterior, 1, mean, na.rm=TRUE )  ))
    pa$sd   = c( inla.mesh.project( pG, field=apply( posterior, 1, sd, na.rm=TRUE )  ))


  }

  if (0) {
    plot.new(); levelplot( mean ~ plon+plat, pa, aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=TRUE) )
    plot.new(); levelplot( sd   ~ plon+plat, pa, aspect="iso", labels=TRUE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
  }

  MESH=NULL

# bathymetry_figures( DS="predictions", p=p ) # to debug
  inla.summary = stmv_summary_inla_spde2 ( RES, SPDE )
  # save statistics last as this is an indicator of completion of all tasks .. restarts would be broken otherwise
  stats = list()
  stats$sdTotal=sdTotal
  stats$rsquared=NA
  stats$ndata=nrow(dat)
  stats$sdSpatial = sqrt(inla.summary[["spatial error", "mode"]])
  stats$sdObs = sqrt(inla.summary[["observation error", "mode"]])
  stats$nu = p$inla.alpha - 1
  stats$phi = matern_phi2phi( inla.summary[["kappa","mean"]], stats$nu, "inla" )
  stats$range = matern_phi2distance(phi=stats$phi, nu=stats$nu, cor=p$stmv_autocorrelation_localrange)
  # stats$range = geoR::practicalRange("matern", phi=stats$phi, kappa=stats$nu  )
  return (list(predictions=pa, stmv_stats=stats))
}
