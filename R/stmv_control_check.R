stmv_control_check = function(p) {
  #// test for presence of control file and if found do some action
  #// permits user control of parallel processes from outside of the run 
  #// the action is prefaced by '#//'

  if ( file.exists( p$stmv_control_file ) )  {

    o = scan( p$stmv_control_file, "character", sep="\n", quiet=TRUE ) # read in file
    k = try( grep ( "\\#('|\\/\\/)", o, ignore.case=TRUE  ), silent=TRUE )  # line number of control
    file.remove( p$stmv_control_file ) 

    if (! inherits(k, "try-error") ) {
      if ( length(k) > 0) {
      
        res = gsub( "^[[:space:]]*\\#('|\\/\\/)", "", o[k] )

        if (grepl( "save", res ) ) {
          stmv_db(p=p, DS="save_current_state", runmode=p$current_runmode, datasubset=c("P", "Pn", "Psd", "statistics") )  
        }

        if (grepl( "status_reset_incomplete", res ) ) {
          currentstatus = stmv_statistics_status( p=p, reset=c( "incomplete" ) ) # flags/filter stats locations based upon prediction covariates. .. speed up and reduce storage
        }

        if (grepl( "status_update", res ) ) {
          stmv_statistics_status( p=p, verbose=FALSE ) # quick update before logging
        }

        if (grepl( "stop", res ) ) {

        }

        if (grepl( "run", res ) ) {
          try( source( file=p$stmv_control_file ) )
        }
      

      }
    }
    message( "Control file found and run")  
  }

}
