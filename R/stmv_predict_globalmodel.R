stmv_predict_globalmodel = function(ip=NULL, p, debugging=FALSE, ... ) {

  if (0) {
    # for debugging  runs ..
    p = parallel_run( p=p, runindex=list( it = 1:p$nt  )
    ip = 1:p$nruns # == 1:p$nt
    debugging=TRUE
  }

  # ---------------------
  # deal with additional passed parameters
  if ( is.null(p) ) p=list()
  p_add = list(...)
  if (length(p_add) > 0 ) p = c(p, p_add)
  i = which(duplicated(names(p), fromLast=TRUE))
  if ( length(i) > 0 ) p = p[-i] # give any passed parameters a higher priority, overwriting pre-existing variable

  if (exists( "libs", p)) suppressMessages( RLibrary( p$libs ) )

  if (is.null(ip)) if( exists( "nruns", p ) ) ip = 1:p$nruns

  global_model = stmv_db( p=p, DS="global_model")
  if (is.null(global_model)) stop("Global model not found.")

  P0 = stmv_attach( p$storage.backend, p$ptr$P0 )  # remember this is on link scale
  P0sd = stmv_attach( p$storage.backend, p$ptr$P0sd ) # and this too


# main loop over each output location in S (stats output locations)
  for ( iip in ip ) {

    # downscale and warp from p(0) -> p1
    pa = NULL # construct prediction surface

    it = p$runs[iip, "it"]

    if ( length(p$variables$COV) > 0 ) {
      for (i in p$variables$COV ) {
        pu = stmv_attach( p$storage.backend, p$ptr$Pcov[[i]] )
        nc = ncol(pu)
        if ( nc== 1 ) {
          pa = cbind( pa, pu[] ) # ie. a static variable (space)
        } else if ( nc == p$nt & nc == p$ny) {
          pa = cbind( pa, pu[,it] ) # ie. same time dimension as predictive data (space.annual.seasonal)
        } else if ( nc == p$ny & p$nt > p$ny)  {
          iy = round( (it-1) / p$nw ) + 1
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
      names(pa) = p$variables$COV
    } else {
      pa = data.frame(intercept=rep(1, length(P0)))  # just to get the size right  when constant intercept model
    }

    if ( any( p$variables$LOCS %in%  p$variables$global_cov ) ) {
      Ploc = stmv_attach( p$storage.backend, p$ptr$Ploc )
      pa = cbind(pa, Ploc[])
      names(pa) = c( p$variables$COV, p$variables$LOCS )
    }

    if ( "yr" %in%  p$variables$global_cov ) {
      npa = names(pa)
      pa = cbind(pa, p$yrs[it] )
      names(pa) = c( npa, "yr" )
    }

    if ( "dyear" %in% p$variables$global_cov ) {
      npa = names(pa)
      pa = cbind(pa, p$prediction.dyear )
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
        YYY = predict( global_model, type="link", se.fit=TRUE ) # determine bounds from data
        Yq = quantile( YYY$fit, probs=p$stmv_quantile_bounds )
        YYY = NULL
        Pbaseline$fit[ Pbaseline$fit < Yq[1] ] = Yq[1]  # do not permit extrapolation
        Pbaseline$fit[ Pbaseline$fit > Yq[2] ] = Yq[2]
        P0[,it] = Pbaseline$fit
        P0sd[,it] = Pbaseline$se.fit
        Pbaseline = NULL
      }

    } else if (p$stmv_global_modelengine %in% c("glm", "bigglm", "gam") ) {

      Pbaseline = try( predict( global_model, newdata=pa, type="link", se.fit=TRUE ) )  # must be on link scale
      pa = NULL
      if (!inherits(Pbaseline, "try-error")) {
        YYY = predict( global_model, type="link", se.fit=TRUE )  # determine bounds from data
        Yq = quantile( YYY$fit, probs=p$stmv_quantile_bounds )
        YYY = NULL
        Pbaseline$fit[ Pbaseline$fit < Yq[1] ] = Yq[1]  # do not permit extrapolation
        Pbaseline$fit[ Pbaseline$fit > Yq[2] ] = Yq[2]
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
          for (j in ip[2:p$nt]){
            P0[,j] = P0[,1]
            P0sd[,j] = P0sd[,1]
          }
        }
      }
    }
  } # end each timeslice
  global_model = NULL
  return(NULL)
}
