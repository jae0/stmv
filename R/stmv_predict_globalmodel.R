stmv_predict_globalmodel = function( ip=NULL, p, Yq_link=NULL, global_model=NULL ) {

  if (0) {
    # for debugging  runs ..
    p = parallel_run( p=p, runindex=list( pnt = 1:p$nt  ) )
    ip = 1:p$nruns # == 1:p$nt
  }

  if (exists( "libs", p)) suppressMessages( RLibrary( p$libs ) )

  ooo = NULL

  if (is.null(ip)) {
    if( exists( "nruns", p ) ) {
      ip = 1:p$nruns
      ooo = p$runs[ip, "pnt"]
    }
  }
  if (is.null(ooo))  ooo = 1:p$nt


  P0 = stmv_attach( p$storage_backend, p$ptr$P0 )  # remember this is on link scale
  P0sd = stmv_attach( p$storage_backend, p$ptr$P0sd ) # and this too

  if (is.null(global_model)) global_model = stmv_global_model( p=p, DS="global_model")
  if (is.null(global_model)) stop("Global model not found.")

# main loop over each output location in S (stats output locations)
  for ( it in ooo ) {

    # downscale and warp from p(0) -> p1
    pa = NULL # construct prediction surface
    if ( length(p$stmv_variables$COV) > 0 ) {
      for (i in p$stmv_variables$COV ) {
        pu = stmv_attach( p$storage_backend, p$ptr$Pcov[[i]] )
        nc = ncol(pu)
        if ( nc== 1 ) {
          pa = cbind( pa, pu[] ) # ie. a static variable (space)
        } else if ( nc == p$nt & nc == p$ny) {
          pa = cbind( pa, pu[,it] ) # ie. same time dimension as predictive data (space.annual.seasonal)
        } else if ( nc == p$ny & p$nt > p$ny)  {
          iy = aegis_floor( (it-1) / p$nw ) + 1L
          pa = cbind( pa, pu[,iy] ) # ie., annual data (space.annual)
        } else if ( nc == p$nt & p$nt > p$ny) {
          pa = cbind( pa, pu[,it] ) # ie. same time dimension as predictive data (space.annual.seasonal)
        } else {
          print(i)
          print(nc)
          stop( "Erroneous data dimension ... the dataset for the above variable looks to be incomplete?")
        }
      }
      pa = as.data.frame( pa )
      names(pa) = p$stmv_variables$COV
    } else {
      pa = data.frame(intercept=rep(1, length(P0)))  # just to get the size right  when constant intercept model
    }

    if ( any( p$stmv_variables$LOCS %in%  p$stmv_variables$global_cov ) ) {
      Ploc = stmv_attach( p$storage_backend, p$ptr$Ploc )
      pa = cbind(pa, Ploc[])
      names(pa) = c( p$stmv_variables$COV, p$stmv_variables$LOCS )
    }

    if ( "yr" %in%  p$stmv_variables$global_cov ) {
      npa = names(pa)
      pa = cbind(pa, p$yrs[it] )
      names(pa) = c( npa, "yr" )
    }

    if ( "dyear" %in% p$stmv_variables$global_cov ) {
      npa = names(pa)
      pa = cbind(pa, p$prediction_dyear )
      names(pa) = c( npa, "dyear" )
    }

    if (0) {
      p$stmv_global_modelengine = "userdefined"
      p$stmv_global_model$run = "
        gam( formula=p$stmv_global_modelformula, data=B,
          optimizer= p$stmv_gam_optimizer, family=p$stmv_global_family, weights=wt ) )
      "
      # must be 'global_model', newdata=pa' in following
      p$stmv_global_model$predict = "
        predict( global_model, newdata=pa, type='link', se.fit=TRUE )
      "
    }

    if ( p$stmv_global_modelengine == "userdefined" )  {
      if (exists( "prediction_method", p$stmv_global_model )) {
        Pbaseline = NULL
        Pbaseline = try( eval(parse(text=p$stmv_global_model$predict )) )
        if (inherits(Pbaseline, "try-error")) Pbaseline = NULL
      } else {
        # try in case of a default predict method exits
        warning( "p$stmv_global_model$predict() method was not found trying to do a generic prediction ...")
        Pbaseline = try( predict( global_model, newdata=pa, type='link', se.fit=TRUE ) )  # must be on link scale
        if (inherits(Pbaseline, "try-error")) stop ("Prediction failed ... ")
      }
      pa = NULL

      if ( ! is.null(Pbaseline) ) {
        # extreme data can make convergence slow and erratic .. this will be used later to limit 99.9%CI
        # do not permit extrapolation
        lb = which( Pbaseline$fit < Yq_link[1])
        ub = which( Pbaseline$fit > Yq_link[2])
        if (length(lb) > 0) Pbaseline$fit[lb] = Yq_link[1]
        if (length(ub) > 0) Pbaseline$fit[ub] = Yq_link[2]
        P0[,it] = Pbaseline$fit
        P0sd[,it] = Pbaseline$se.fit
        Pbaseline = NULL
      }

    } else if (p$stmv_global_modelengine %in% c("glm", "bigglm", "gam") ) {

      Pbaseline = try( predict( global_model, newdata=pa, type="link", se.fit=TRUE ) )  # must be on link scale
      pa = NULL
      if (!inherits(Pbaseline, "try-error")) {
        # do not permit extrapolation
        lb = which( Pbaseline$fit < Yq_link[1])
        ub = which( Pbaseline$fit > Yq_link[2])
        if (length(lb) > 0) Pbaseline$fit[lb] = Yq_link[1]
        if (length(ub) > 0) Pbaseline$fit[ub] = Yq_link[2]
        P0[,it] = Pbaseline$fit
        P0sd[,it] = Pbaseline$se.fit
      }
      Pbaseline = NULL
    } else if (p$stmv_global_modelengine =="bayesx") {
      stop( "not yet tested" ) # TODO
      # Pbaseline = try( predict( global_model, newdata=pa, type="link", se.fit=TRUE ) )
      # pa = NULL
      # if (!inherits(Pbaseline, "try-error")) {
      #   P0[,it] = Pbaseline$fit
      #   P0sd[,it] = Pbaseline$se.fit
      # }
      # Pbaseline = NULL

    } else if (p$stmv_global_modelengine =="none") {
      # nothing to do
    } else  {
      stop ("This global model method requires a bit more work .. ")
    }

    if (exists("all.covars.static", p)) {
      if (p$all.covars.static) {
      # if this is true then this is a single cpu run and all predictions for each time slice is the same
      # could probably catch this and keep storage small but that would make the update math a little more complex
      # this keeps it simple with a quick copy
        if (p$nt  > 1 ) {
          for (j in ooo[-1]){
            P0[,j] = P0[,ooo[1]]
            P0sd[,j] = P0sd[,ooo[1]]
          }
        }
        break()  # escape for loop
      }
    }
  } # end each timeslice

  return(NULL)

    if (0) {
      if ("all nt time slices in stored predictions P") {
        Ploc = stmv_attach( p$storage_backend, p$ptr$Ploc )
        # pa comes from stmv_interpolate ... not carried here
        for (i in 1:p$nt) {
          i = 1
          print( lattice::levelplot( P0[,i] ~ Ploc[,1] + Ploc[ , 2], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
        }
      }
      if ("no time slices in P") {
        Ploc = stmv_attach( p$storage_backend, p$ptr$Ploc )
          print( lattice::levelplot( P0[] ~ Ploc[,1] + Ploc[ , 2], col.regions=heat.colors(100), scale=list(draw=FALSE) , aspect="iso" ) )
      }
    }

}