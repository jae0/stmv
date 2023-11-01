


project.library( "stmv")

datafile = file.path( "c:", "nurses.data.r" )
# datafile = file.path( "~", "nath", "geocoding", "nurses.data.r" )


res$addr = data table with adddresses of residences
res$lon = NA
res$lat = NA

for ( i in 1:nrow(res) ) {
	for ( wt in c(1, 5, 10, 20, 30, 60 ) ) {
		result = geocode( res$addr[i] )
		if ( result$status == "good" ) {
			break()
		} else {
			print ( "Geocode request is too rapid, waiting a little longer ..." )
			pause (wt) # wait a few seconds before retrying
		}
	}	# end for wt
	if (result$status == "good" ) {
		res$lon[i] = result$lon
		res$lat[i] = result$lat
		next()
	} else {
		print( "Problem geocoding:" )
		print( res[i, ] )
	}

} # end for res


