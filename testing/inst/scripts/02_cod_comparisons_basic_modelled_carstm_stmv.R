
# ------------------------------------------------
# Atlantic cod comparison using basic models testing the Gaussian vs Poisson assumption

# Instead of treating things as a normal,
# we treat the data as being derived from a Poisson process.
# So we focus upon counts of cod in each set with an "offset" of
# the surface area of each tow, either as a "standard tow",
# actual swept area calculations, or simply tow length adjusted
# "areas". The "standard method" is to assume "standard tow"... due to issues
# of insufficient data with sweptarea and being too naive with standard tow

# After modelling numerical distribution of counts under a Poisson assumption,
# we can back calculate the associated weights
# using the weighted average mass of each tow-set.

# These models try to show that a model formulation reasonably computes
# stratanl solutions but for this to work one must parameterize interactions-only
# with no other model components, including the intercept ...




# ------------------------------------------------
# load data common environment and parameter setting
# source( system.file( "scripts", "00_cod_comparisons_data_environment.R", package = "carstm") )
# --------------------------------
# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
# NOTE: the data selection is the same as in (01_cod_comparisons_basic_stranal.R)
p = list(
  speciesname = "Atlantic_cod",
  groundfish_species_code = 10,   #  10= cod
  yrs = 1970:2017,
  trawlable_units = "towdistance"  # <<<<<<<<<<<<<<<<<<
  # trawlable_units = "standardtow"
  # trawlable_units = "sweptarea"
)



# --------------------------------
# parameter setting used to filter data via 'survey_db( DS="filter")'
# unlike stratanl, we do not need to remove strata until the last /aggregation step
p = aegis.survey::survey_parameters(
  p=p,
  selection=list(
    biologicals=list(
      spec_bio = bio.taxonomy::taxonomy.recode( from="spec", to="parsimonious", tolookup=p$groundfish_species_code )
    ),
    survey=list(
      data.source="groundfish",
      yr = p$yrs,  # time frame for comparison specified above
      months=6:8,  # "summer"
      # dyear = c(150,250)/365, # alternate way of specifying season: summer = which( (x>150) & (x<250) ) , spring = which(  x<149 ), winter = which(  x>251 )
      settype = 1, # same as geartype in groundfish_survey_db
      gear = c("Western IIA trawl", "Yankee #36 otter trawl"),
      polygon_enforce=TRUE,  # make sure mis-classified stations or incorrectly entered positions get filtered out
      ranged_data = c("dyear")  # not used .. just to show how to use range_data
    )
  )
)


# ------------------------------------------------
## using the "standard" polygon definitions  .. see https://cran.r-project.org/web/packages/spdep/vignettes/nb.pdf
# Here we compute surface area of each polygon via projection to utm or some other appropriate planar projection.
# This adds some variabilty relative to "statanal" (which uses sa in sq nautical miles, btw)

crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
sppoly = areal_units( areal_units_type="stratanal_polygons_pre2014",areal_units_proj4string_planar_km=p$areal_units_proj4string_planar_km )
sppoly = st_transform( sppoly, crs=crs_lonlat )

sppoly$strata_to_keep = ifelse( as.character(sppoly$AUID) %in% strata_definitions( c("Gulf", "Georges_Bank", "Spring", "Deep_Water") ), FALSE,  TRUE )


# ------------------------------------------------
# neighbourhood structure --- required to do areal unit spatial modelling
# sppoly = neighbourhood_structure( sppoly=sppoly, areal_units_type="stratanal_polygons_pre2014" )  # not used here


# --------------------------------
# Get the data
p$selection$survey$strata_toremove = NULL  # emphasize that all data enters analysis initially ..

set = survey_db( p=p, DS="filter" )

set$AUID = st_points_in_polygons(
  pts = st_as_sf( set, coords=c("lon","lat"), crs=crs_lonlat ),
  polys = sppoly[, "AUID"],
  varname="AUID"
)

set = set[ which(!is.na(set$AUID)),]

set$totno[which(!is.finite(set$totno))] = NA


# --------------------------------
# ensure we have some estimate of sweptarea and choose the appropriate
# one based upon which trawlable units we are using
ft2m = 0.3048
m2km = 1/1000
nmi2mi = 1.1507794
mi2ft = 5280
standardtow_sakm2 = (41 * ft2m * m2km ) * ( 1.75 * nmi2mi * mi2ft * ft2m * m2km )  # surface area sampled by a standard tow in km^2  1.75 nm
set$data_offset = switch( p$trawlable_units,
  standardtow =  rep(standardtow_sakm2, nrow(set)) , # "standard tow"
  towdistance = set$sa_towdistance,  # "sa"=computed from tow distance and standard width, 0.011801==),
  sweptarea = set$sa  # swept area based upon stand tow width and variable lenths based upon start-end locations wherever possible
)
set$data_offset[which(!is.finite(set$data_offset))] = median(set$data_offset, na.rm=TRUE )  # just in case missing data


# ------------------------------------------------
# update set with AUID factor variables and a few other repeatedly used variables
set$AUID = factor(set$AUID, levels=levels(sppoly$AUID))
set$yr_factor = factor(set$yr)
set$iid_error = 1:nrow(set) # for inla indexing
set$tag = "observations"


## --------------------------------
# construct meanweights matrix
weight_year = meanweights_by_arealunit( set=set, AUID=as.character( sppoly$AUID ), yrs=p$yrs, fillall=TRUE, annual_breakdown=TRUE )
# weight_year = meanweights_by_arealunit_modelled( p=p, redo=TRUE )  -- note: data passing of M needs to be modularized 

# weight_year = weight_year[, match(as.character(p$yrs), colnames(weight_year) )]
# weight_year = weight_year[ match(as.character(sppoly$AUID), rownames(weight_year) )]



#  good data
ok = which(
  is.finite(set$totno) &   # INLA can impute Y-data
  is.finite(set$t) &
  is.finite(set$z) &
  is.finite(set$data_offset) &
  set$AUID %in% sppoly$AUID[sppoly$strata_to_keep]
)

if (0) {
  # RES = data.frame(yr=p$selection$survey[["yr"]]) # collect model comparisons
  fn = file.path( getwd(), "RES.rdata" )
  # save(RES, file=fn)
  # load(fn)
}

# ------------------------------------------------
# Model 1: simple factorial operating directly upon totwgt
# this replicates the stratanal estimates but in an explicit model-based framework

fit = lm(
  formula = totwgt ~ AUID:yr_factor - 1  ,
  data=set[ok,]
)
s = summary(fit)
AIC(fit)  # 106879

# means are equal to model coefficients in this simple model
means = split_vector2matrix(
  u=s$coefficients[,"Estimate"],
  vnames=c("AUID", "yr_factor"),
  sep=":",
  matchto=list( AUID=sppoly$AUID, yr_factor=levels(set$yr_factor))
)

  if (0) {
    # or more "simply" as:
    u = s$coefficients[,"Estimate"]
    u_sd = s$coefficients[,"Std. Error"]
    w = matrix( unlist(strsplit(names(u), ":")), ncol=2, byrow=TRUE)
    w = gsub("AUID", "", w)
    w = gsub("yr_factor", "", w)
    sid = match( w[,1], as.character(sppoly$AUID) )
    yid = match( w[,2], as.character(levels(set$yr_factor)))
    means = matrix( NA, ncol=max(yid), nrow=max(sid) )
    means[ cbind(sid, yid)] = u
  }


# SD's .. one would need to propagate the errors (assume something like first order, or via simulation ... )
sds = split_vector2matrix(
  u=s$coefficients[,"Std. Error"],
  vnames=c("AUID", "yr_factor"),
  sep=":",
  matchto=list( AUID=sppoly$AUID, yr_factor=levels(set$yr_factor))
)

out = means * sppoly$au_sa_km2 / standardtow_sakm2

RES$lm_gaussian_totwgt = colSums( out[sppoly$strata_to_keep,], na.rm=TRUE )

# plot( lm_gaussian_totwgt ~ yr, data=RES, lty=2, lwd=2, col="green", type="b", ylim=c(0,9e8))
lines( lm_gaussian_totwgt ~ yr, data=RES, lty=2, lwd=2, col="green", type="b")


vn = "lm_gaussian_totwgt"
sppoly[,vn] = out[,"2017"] / standardtow_sakm2
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )


# problem:
# - data distribution is not normal .. perhaps lognormal, poisson or overdispersed something
hist(set$totwgt, breaks=50)
hist(set$totno, breaks=50)



# ------------------------------------------------
# Model 3: simple factorial with totwgt using glm
# this permits other data distributions in a model-based framework ...
# here Gaussian to check methods to make sure it is the same as the lm version
# problems:
# - operating directly upon biomass as a gaussian is biased as this looks to be ~ lognormal ..
#   but lognormal requires more complex models to handles the zero-values (hurdle, delta, etc)

fit = glm(
  formula = totwgt ~ AUID:yr_factor - 1  ,
  family=gaussian(link="identity"),
  data=set[ok,]
)
AIC(fit)  # 106879
s = summary(fit)

# means are equal to model coefficients in this simple model
means = split_vector2matrix(
  u=s$coefficients[,"Estimate"],
  vnames=c("AUID", "yr_factor"),
  sep=":",
  matchto=list( AUID=sppoly$AUID, yr_factor=levels(set$yr_factor))
)

out = means * sppoly$au_sa_km2 / standardtow_sakm2
RES$glm_gaussian_totwgt = colSums( out[sppoly$strata_to_keep,], na.rm=TRUE )
lines( glm_gaussian_totwgt ~ yr, data=RES, lty=3, lwd=3, col="black")


# ------------------------------------------------
# Model 4: Poisson distribution is perhaps appropriate due to the shape of the distribution
# and the discrete nature of the counts which are mostly low in magnitude
# Use meanweight estimates for each year and stratum to convert to biomass .. of course there is error here too but ignored (for now)
# problem: operating directly upon the totno counts assumes sampling effort (swept area, catchability, substrate, boat speed, net type, etc)
# is identical ("standard tow") .. but in reality they are not all equal

fit = glm(
  formula = totno ~ AUID:yr_factor - 1  ,
  family=poisson(link="log"),
  data=set[ok,]
)  # slow
s = summary(fit)
AIC(fit)  # 355470

#  prediction surface as a data frame
APS = expand.grid( AUID=levels(set$AUID), yr_factor=levels(set$yr_factor) )

kk = which( paste(APS$AUID, APS$yr_factor) %in% unique(paste( set$AUID, set$yr_factor)[ok] ))
preds = predict( fit, newdata=APS[kk,], type="response", na.action=na.omit, se.fit=TRUE )

APS$predictions = NA
APS$predictions[kk] = preds$fit
APS$predictions.sd = NA
APS$predictions.sd[kk] = preds$se.fit

# reformat predictions into matrix form
out = reformat_to_array(
  input = APS$predictions,
  matchfrom = list( AUID=APS$AUID, yr_factor=APS$yr_factor),
  matchto   = list( AUID=sppoly$AUID, yr_factor=levels(set$yr_factor))
)

# convert numbers/set to numbers/km to biomass/strata
RES$glm_poisson_totno = colSums(
  {out/standardtow_sakm2 * weight_year * sppoly$au_sa_km2}[sppoly$strata_to_keep , ], na.rm=TRUE )

lines( glm_poisson_totno ~ yr, data=RES, lty=4, lwd=8, col="slateblue")



dev.new()
plot(glm_poisson_totno ~ lm_gaussian_totwgt, RES  )
abline(a=0,b=1)



# ------------------------------------------------
# Model 5: poisson, simple factorial with totno and meanweight and an offset
# A better way to represent the above is to use an offset term that permits differntial sampling power or effort.
# This permits incorporation of tow-based, species-based, swept area-based, etc. adjustmensts
# this uses a constant offset and so is equivalent to Model 5 .. it is a test of the method .. and it looks good..
# problem:  adding spatial and/or temporal effects is a bit of a challenge in GLM

fit = glm(
  formula = totno ~ offset(log(data_offset)) + AUID:yr_factor - 1 ,
  family=poisson(link="log"),
  data=set[ok,]
)
AIC(fit)  # 359785
s = summary(fit)

#  prediction surface as a data frame
APS = expand.grid( AUID=levels(set$AUID), yr_factor=levels(set$yr_factor) )
APS$data_offset = 1  # force predictions to be per km^2

kk = which( paste(APS$AUID, APS$yr_factor) %in% unique(paste( set$AUID, set$yr_factor)[ok] ))
preds = predict( fit, newdata=APS[kk,], type="response", na.action=na.omit, se.fit=TRUE )

APS$predictions = NA
APS$predictions[kk] = preds$fit
APS$predictions.sd = NA
APS$predictions.sd[kk] = preds$se.fit

# reformat predictions into matrix form
out = reformat_to_array(
  input = APS$predictions,
  matchfrom = list( AUID=APS$AUID, yr_factor=APS$yr_factor),
  matchto   = list( AUID=sppoly$AUID, yr_factor=levels(set$yr_factor))
)

# convert numbers/km to biomass/strata ..
RES$glm_poisson_totno_offset = colSums( {out * weight_year * sppoly$au_sa_km2}[sppoly$strata_to_keep , ], na.rm=TRUE )
lines( glm_poisson_totno_offset ~ yr, data=RES, lty=5, lwd=4, col="red")


vn = "glm_poisson_totno_offset"
sppoly[,vn] = out[,"2010"]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )
