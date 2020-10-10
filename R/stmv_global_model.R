stmv_global_model = function( p, DS, B=NULL, Yq_link=NULL, Ydata=NULL, saveresults=TRUE ) {

    #// B is the xyz or xytz data or the function to get the data to work upon
    #// require Ydata if user_defined

    if (DS %in% c("global_model", "global_model.redo") ) {

      fn.global_model = file.path( p$stmvSaveDir, paste( "global_model", p$stmv_global_modelengine, "rdata", sep=".") )

      if (DS =="global_model") {
        global_model = NULL
        if (file.exists( fn.global_model ))  load(fn.global_model)
        return(global_model)
      }

      if ( file.exists( fn.global_model ) ) {
        message( "||| A global model already exists. It will be overwritten in 10 seconds.")
        message( "|||   Type <ctrl-c> or <esc> to interrupt. To reuse the saved model ")
        message( "|||   'p$stmv_runmode$globalmodel = FALSE' ... overwriting in:")
        for (i in 9:0) {
          Sys.sleep(1)
          cat( i)
          cat(" ")
        }
      }

      if ( length( p$stmv_variables$COV ) > 0 ) {
        good = which( is.finite (rowSums(B[ , c(p$stmv_variables$Y,p$stmv_variables$COV) ])) )
      } else {
        good = which( is.finite (B[,p$stmv_variables$Y ] ) )
      }

      ngood = length(good)
      if ( ngood > 0 ) {
        if ( ngood < nrow(B) ) {
          B = B[good,]
        }
      }

      if ( p$stmv_global_modelengine == "userdefined" )  {
        # need a file with predictions a function to get/return them
        if (!(exists("run", p$stmv_global_model ))) {
          message( "Must define functions of the form: ")
          message( " p$stmv_global_model$run =' " )
          message( "   gam( formula=p$stmv_global_modelformula, data=B, ")
          message( "     optimizer= p$stmv_gam_optimizer, family=p$stmv_global_family, weights=wt )' " )
          stop()
        }
        global_model = try( eval(parse(text=p$stmv_global_model$run )) )
        if (inherits(global_model, "try-error") ) stop( "Global modelling failed")
        print( summary( global_model ) )
        if (saveresults) save( global_model, file= fn.global_model, compress=TRUE )
        return ()
      }

      # as a first pass, model the time-independent factors as a user-defined model
      if (p$stmv_global_modelengine=="glm") {
        warning( "glm does not have a spatial random effect, model covariates could be biased if residual errors are not iid ... ")
        if (!exists("wt", B)) {
          global_model = try(
            glm( formula=p$stmv_global_modelformula, data=B, family=p$stmv_global_family ) )
        } else {
          global_model = try(
            glm( formula=p$stmv_global_modelformula, data=B, family=p$stmv_global_family, weights=wt ) )
        }
      }

      if (p$stmv_global_modelengine=="bigglm") {
        warning( "bigglm does not have a spatial random effect, model covariates could be biased if residual errors are not iid ... ")
        if (!exists("wt", B)) {
          global_model = try(
            bigglm( formula=p$stmv_global_modelformula, data=B, family=p$stmv_global_family ))
        } else {
          global_model = try(
            bigglm( formula=p$stmv_global_modelformula, data=B, family=p$stmv_global_family, weights=~wt ))
        }
      }

      if (p$stmv_global_modelengine=="gam") {
        warning( "gam does not have a spatial random effect, model covariates could be biased if residual errors are not iid, spatial splines might help but this is not guaranteed. gamm might be a better choice. You have been warned ... ")
        require(mgcv)
        if (!exists("wt", B)) {
          global_model = try(
            gam( formula=p$stmv_global_modelformula, data=B,
              optimizer= p$stmv_gam_optimizer, family=p$stmv_global_family) )
        } else {
          global_model = try(
            gam( formula=p$stmv_global_modelformula, data=B,
              optimizer= p$stmv_gam_optimizer, family=p$stmv_global_family, weights=wt ) )
        }
      }


      if (p$stmv_global_modelengine=="gamm") {
        warning( "gamm should include a spatial random effect, otherwise model covariates could be biased if residual errors are not iid ... ")
        require(mgcv)
        if (!exists("wt", B)) {
          global_model = try(
            gam( formula=p$stmv_global_modelformula, data=B,
              optimizer= p$stmv_gam_optimizer, family=p$stmv_global_family) )
        } else {
          global_model = try(
            gamm( formula=p$stmv_global_modelformula, data=B,
              optimizer= p$stmv_gam_optimizer, family=p$stmv_global_family, weights=wt ) )
        }
      }


      if (p$stmv_global_modelengine=="bayesx") {
        warning( "make sure you have a spatial random effect, otherwise model covariates/predictions could be biased if residual errors are not iid ... ")

        require(mgcv)
        global_model = try(
          bayesx( formula=p$stmv_global_modelformula, data=B, family=p$stmv_global_family ) )
      }


      if (p$stmv_global_modelengine=="inla") {
        warning( "make sure you have a spatial random effect, otherwise model covariates/predictions could be biased if residual errors are not iid ... ")

        require(INLA)
        global_model = try( eval( p$stmv_global_modelformula ) )    # p$stmv_global_modelformula must contain the whole call

          # inla( formula=p$stmv_global_modelformula, data=B, family=p$stmv_global_family ) )
      }

      if ( "try-error" %in% class(global_model) ) stop( "The covariate model was problematic" )
      print( summary( global_model ) )

      if (saveresults) {
        save( global_model, file= fn.global_model, compress=TRUE )
        global_model = NULL
        return()
      } else {
        return( global_model )
      }
    }


    # -----

  if (DS=="residuals") {
    global_model = stmv_global_model( p=p, DS="global_model")
    if ( p$stmv_global_modelengine == "userdefined") {
      # TODO MUST find a generic form as below
      # # Ypreds = predict(global_model, type="link", se.fit=FALSE )  ## TODO .. keep track of the SE
      if (!exists("predict", p$stmv_global_model)) {
        message( "||| p$stmv_global_model$predict =' " )
        message( "|||   predict( global_model, newdata=pa, type='link', se.fit=TRUE )' " )
        message( "||| where 'global_model', newdata=pa' are required " )
        stop()
      }
      Ypreds = try( eval(parse(text=pp$stmv_global_model$predict )) )
      Ypreds = p$stmv_global_model$predict( global_model )  # needs to be tested. .. JC
      Yresids = Ydata - Ypreds # ie. internalR (link) scale
      Ypreds = NULL
    } else {
      # at present only those that have a predict and residuals methods ...
      Yresids  = residuals(global_model, type="working") # ie. internal (link) scale
    }
    global_model =NULL
    if (!is.null(Yq_link)) {
      # could operate upon quantiles of residuals but in poor models this can hyper inflate errors and slow down the whole estimation process
      # truncating using data range as a crude approximation of overall residual and prediction scale
      lb = which( Yresids < Yq_link[1])
      ub = which( Yresids > Yq_link[2])
      if (length(lb) > 0) Yresids[lb] = Yq_link[1]
      if (length(ub) > 0) Yresids[ub] = Yq_link[2]
    }
    return( Yresids )
  }
}
