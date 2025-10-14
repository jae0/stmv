



# ------------------------------------------------
# Atlantic cod comparison of CAR (ICAR/BYM) Poisson process models
# using sweptarea only
#
# Now we move to the use of ICAR-BYM type approach and environmental covariates.

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
      settype = 1, 
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

sppoly = areal_units(
  areal_units_type="stratanal_polygons_pre2014",
  areal_units_proj4string_planar_km=p$areal_units_proj4string_planar_km
)
sppoly$strata_to_keep = ifelse( as.character(sppoly$AUID) %in% strata_definitions( c("Gulf", "Georges_Bank", "Spring", "Deep_Water") ), FALSE,  TRUE )


# ------------------------------------------------
# neighbourhood structure --- required to do areal unit spatial modelling

sppoly = neighbourhood_structure( sppoly=sppoly, areal_units_type="stratanal_polygons_pre2014" )


# --------------------------------
# Get the data
p$selection$survey$strata_toremove = NULL  # emphasize that all data enters analysis initially ..

set = survey_db( p=p, DS="filter" )

# categorize Strata
crs_lonlat = st_crs(projection_proj4string("lonlat_wgs84"))
sppoly = st_transform(sppoly, crs=crs_lonlat )

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
# set$AUID = factor(set$AUID, levels=levels(sppoly$AUID))

space.id = slot( sppoly, "space.id" )
set$space = set$space_time = match( set$AUID, space.id )


set$yr_factor = factor(set$yr)

set$time = set$time_space = as.numeric(  set$yr_factor)

set$iid_error = 1:nrow(set) # for inla indexing
set$tag = "observations"


## --------------------------------
# construct meanweights matrix
weight_year = meanweights_by_arealunit( set=set, AUID=as.character( sppoly$AUID ), yrs=p$yrs, fillall=TRUE, annual_breakdown=TRUE )
# weight_year = meanweights_by_arealunit_modelled( p=p, redo=TRUE )  -- note: data passing of M needs to be modularized 

# weight_year = weight_year[, match(as.character(p$yrs), colnames(weight_year) )]
# weight_year = weight_year[ match(as.character(sppoly$AUID), rownames(weight_year) )]



# adjust based upon RAM requirements and ncores
ncores = floor( ram_local( "ncores", ram_main=4, ram_process=6 ) / 2 )
inla.setOption(num.threads=ncores)
inla.setOption(blas.num.threads=ncores)



# RES = data.frame(yr=p$selection$survey[["yr"]]) # collect model comparisons
if (0) {
  fn = file.path( getwd(), "RES.rdz" )
  # read_write_fast(RES, file=fn)
  # RES = read_write_fast(fn) 
}



## ----------------------------------
# covariates of interest
covars = c("t", "tsd", "tmax", "tmin", "degreedays", "z",  "dZ", "ddZ" )

  # currently supported:
  # z = depth (m)
  # dZ = bottom slope (m/km)
  # ddZ = bottom curvature (m/km^2)
  # substrate.grainsize = mean grain size of bottom substrate (mm)
  # t = temperature (C) – subannual
  # tlb = temperature lower 95% bound (C) –subannual
  # tub = temperature upper 95% bound (C) –subannual
  # tmean = mean annual temperature
  # tsd = standard deviation of the mean annual temperature
  # tmin = minimum value of temperature in a given year – annual
  # tmax= maximum value of temperature in a given year – annual
  # tamplitude = amplitude of temperature swings in a year (tmax-tmin) – annual
  # degreedays = number of degree days in a given year – annual


# extract covariate means by strata
res = aegis_db_extract_by_polygon(
  sppoly=sppoly,
  vars=covars,
  spatial_domain=p$spatial_domain,
  yrs=p$yrs,
  dyear=0.6 # 0.6*12 months = 7.2 = early July
)

# extract covariates and supplent survey data via lookups
set = aegis_db_lookup(
  X=set,
  lookupvars=covars,
  xy_vars=c("lon", "lat"),
  time_var="timestamp"
)

#  good data
ok = which(
  is.finite(set$totno) &
  is.finite(set$t) &
  is.finite(set$z) &
  is.finite(set$data_offset) &
  set$AUID %in% sppoly$AUID[sppoly$strata_to_keep]
)



# merge data into prediction surface and add tags

APS = aegis_prediction_surface( aegis_data=res$means  )
APS$yr = as.numeric( APS$year)
APS$totno = NA
APS$data_offset = 1  # force to be density n/km^2
APS$tag = "predictions"

varstokeep = c( "totno", "AUID", "yr", "t", "z", "data_offset", "tag" )

M = rbind( set[ok, varstokeep], APS[,varstokeep] )

M$t[!is.finite(M$t)] = median(M$t, na.rm=TRUE )  # missing data .. quick fix .. do something better
M$z[!is.finite(M$z)] = median(M$z, na.rm=TRUE )  # missing data .. quick fix .. do something better

M$yr_factor = factor( as.character(M$yr) )
M$AUID  = factor( M$AUID, levels=levels(sppoly$AUID ))
M$strata  = as.numeric( M$AUID)
M$year  = as.numeric( M$yr_factor)

M$ti = discretize_data( x=M$t, brks=p$discretization$t )
M$zi = discretize_data( x=M$t, brks=p$discretization$z )


M$iid_error = 1:nrow(M) # for inla indexing for set level variation


# ---------------------
# generic PC priors
m = log( {set$totno / set$data_offset}[ok] )
m[!is.finite(m)] = min(m[is.finite(m)])

H = inla_hyperparameters( sd(m), alpha=0.5, median(m) )
# H$prec$prec.intercept = 1e-9




# ------------------------------------------------
# Model 12:
# "INLA Envir AR1	iid|year"	ar1	rw2: temp+depth 	Poisson	5966	33681	5

# simple factorial with totno and poisson
# improvement upon Model 6b .. INLA imputes missing data given proper data model .. less fiddling ..
# random effects.. `fewer params`

# interaction-only model does not make sense here .. esp as there are other model components
fit = inla(
  formula =
    totno ~ 1 + offset( log( data_offset) )  # CARSTM does log-transformation internally  but this is a direct call to inla 
      + f(strata, model="iid", group=year, hyper=H$iid)
      + f(year, model="ar1", hyper=H$ar1 )
      + f(iid_error, model="iid", hyper=H$iid)
      + f(ti, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2)
      + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-6, hyper=H$rw2),
  family = "poisson", # "zeroinflatedpoisson0",
  data= M,
  control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=FALSE, link=1 ),
  # control.fixed=H$fixed,  # priors for fixed effects, generic is ok
  # control.inla=list( strategy="laplace", cutoff=1e-6, correct=TRUE, correct.verbose=FALSE ),
  verbose=FALSE
)

s = summary(fit)
s$dic$dic  # 33681
s$dic$p.eff # 5966

plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_array(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( AUID=M$AUID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( AUID=sppoly$AUID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$INLA.Envir.AR1.iid_year = colSums( {out * weight_year * sppoly$au_sa_km2}[sppoly$strata_to_keep, ], na.rm=TRUE )

lines( INLA.Envir.AR1.iid_year ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")

# map it
vn = "pred"
yr = "2017"
sppoly[,vn] = out[,yr] * weight_year[,yr]  # biomass density
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )




# ------------------------------------------------
# Model 13:
# "INLA Envir CAR" :	bym2	iid	rw2: temp+depth 	Poisson	6006	33743	6

# CAR with totno and poisson but no AR1 just iid
# in model 11, we saw that a fixed effects model does not work very well
# .. that is, the data model has additional error sources that are not well modelled

# add a single CAR structure that is constant (shared) across time (years)

#   family = "zeroinflatedpoisson0" does complete .. but NA DICs  and solutions that make no sense
#   family = "zeroinflatedpoisson1" does not complete ..
#  use "poisson" as it is the simplest ..


fit = inla(
  formula = totno ~ 1
    + offset( log(data_offset) ) # CARSTM does log-transformation internally  but this is a direct call to inla 
    + f(iid_error, model="iid", hyper=H$iid)
    + f(ti, model="rw2", scale.model=TRUE, hyper=H$rw2)
    + f(zi, model="rw2", scale.model=TRUE, hyper=H$rw2)
    + f(year, model="iid", hyper=H$iid)
    + f(strata, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2),
  family = "poisson",
  data=M,
  control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=TRUE, link =1 ), # compute=TRUE on each data location
  control.fixed=H$fixed,  # priors for fixed effects
  control.inla=list(  correct=TRUE, correct.verbose=FALSE ), # strategy="laplace", cutoff=1e-6,
  verbose=TRUE
)
s = summary(fit)
s$dic$dic  #  33748
s$dic$p.eff # 6056

plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_array(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( AUID=M$AUID[ii], yr_factor=M$yr_factor[ii] ),
  matchto   = list( AUID=sppoly$AUID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$INLA.Envir.CAR = colSums( {out * weight_year * sppoly$au_sa_km2}[sppoly$strata_to_keep, ], na.rm=TRUE )

lines( INLA.Envir.CAR ~ yr, data=RES, lty=1, lwd=2.5, col="red", type="b" )


dev.new();
plot( fit$marginals.hyperpar$"Phi for strata", type="l")  # posterior distribution of phi nonspatial dominates



# ------------------------------
# Model 14 CAR and AR1
# "INLA Envir AR1 CAR" 	bym2	ar1	rw2: temp+depth 	Poisson	6010	33769	7
# 45 min

fit = inla(
  formula = totno ~ 1
    + offset( log(data_offset) ) # CARSTM does log-transformation internally  but this is a direct call to inla 
#    + f(iid_error, model="iid", hyper=H$iid)
    + f(ti, model="rw2", scale.model=TRUE, hyper=H$rw2)
    + f(zi, model="rw2", scale.model=TRUE, hyper=H$rw2)
    + f(year, model="ar1", hyper=H$ar1)
    + f(strata, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2),
  family = "poisson",
  data=M,
  control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=TRUE, link =1 ), # compute=TRUE on each data location
  control.fixed=H$fixed,  # priors for fixed effects
  control.inla=list(  correct=TRUE, correct.verbose=FALSE ), # strategy="laplace", cutoff=1e-6,
  verbose=TRUE
)

s = summary(fit)
s$dic$dic  #   [1] 33810
s$dic$p.eff # [1] 5976


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_array(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( AUID=M$AUID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( AUID=sppoly$AUID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA
RES$INLA.Envir.AR1.CAR = colSums( {out * weight_year * sppoly$au_sa_km2}[sppoly$strata_to_keep, ], na.rm=TRUE )

lines( INLA.Envir.AR1.CAR ~ yr, data=RES, lty=1, lwd=2.5, col="red", type="b")


# map it ..mean density
vn = "pred"
sppoly[,vn] = out[,"2017"]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )


plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )


plot( fit$marginals.hyperpar$"Phi for strata", type="l")  # posterior distribution of phi nonspatial dominates
plot( fit$marginals.hyperpar$"Precision for strata", type="l")
plot( fit$marginals.hyperpar$"Precision for setno", type="l")



# ------------------------------------------------
# Model 15:
# "INLA Envir AR1 CAR|year"	 bym2|year	ar1	rw2: temp+depth 	Poisson	5943	33650	3

# CAR with totno and poisson ... a separate CAR process for every year .. and AR1

# 20+ hrs

fit = inla(
  formula = totno ~ 1
    + offset( log(data_offset) ) # CARSTM does log-transformation internally  but this is a direct call to inla 
    + f(iid_error, model="iid", hyper=H$iid)
    + f(ti, model="rw2", scale.model=TRUE, hyper=H$rw2)
    + f(zi, model="rw2", scale.model=TRUE, hyper=H$rw2)
    + f(year, model="ar1", hyper=H$ar1 )
    + f(strata, model="bym2", graph=slot(sppoly, "nb"), group=year, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
  family = "poisson",
  data=M,
  control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=TRUE, link =1 ), # compute=TRUE on each data location
  control.fixed=H$fixed,  # priors for fixed effects
  control.inla=list(  correct=TRUE, correct.verbose=FALSE ), # strategy="laplace", cutoff=1e-6,
  verbose=TRUE
)

# save as it takes so long
if (0) {
  fn15 = "~/tmp/car_annual15.rdz"  # 1.2GB
  read_write_fast( fit, file=fn15 )
  fit = read_write_fast(fn15)
}

s = summary(fit)
s$dic$dic  #   33655
s$dic$p.eff #  5945

plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_array(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( AUID=M$AUID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( AUID=sppoly$AUID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA

RES$INLA.Envir.AR1.CAR_year = colSums( {out * weight_year * sppoly$au_sa_km2}[sppoly$strata_to_keep,], na.rm=TRUE )

lines( INLA.Envir.AR1.CAR_year ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")

# map it
vn = "pred"
sppoly[,vn] = out[,"2017"]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )

plot( fit$marginals.hyperpar$"Phi for strata", type="l")  # posterior distribution of phi
plot( fit$marginals.hyperpar$"Precision for strata", type="l")
plot( fit$marginals.hyperpar$"Precision for setno", type="l")



plot( fit$marginals.hyperpar$"Phi for strata", type="l")  # posterior distribution of phi nonspatial dominates
plot( fit$marginals.hyperpar$"Precision for strata", type="l")
plot( fit$marginals.hyperpar$"Precision for setno", type="l")


# ------------------------------------------------
# Model 16:
# "INLA Envir AR1|strata CAR" 	bym2	ar1|strata	rw2: temp+depth 	Poisson	5903	33552	2

# CAR with totno and poisson ...
# -- a single CAR process for all years
# -- AR1 on year by each stratum
# 8+ hrs

fit = inla(
  formula = totno ~ 1
    + offset( log(data_offset) ) # CARSTM does log-transformation internally  but this is a direct call to inla 
    + f(iid_error, model="iid", hyper=H$iid)
    + f(ti, model="rw2", scale.model=TRUE, hyper=H$rw2)
    + f(zi, model="rw2", scale.model=TRUE, hyper=H$rw2)
    + f(year, model="ar1", hyper=H$ar1, group=strata )
    + f(strata, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2),
  family = "poisson",
  data=M,
  control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=TRUE, link =1 ), # compute=TRUE on each data location
  control.fixed=H$fixed,  # priors for fixed effects
  control.inla=list(  correct=TRUE, correct.verbose=FALSE ), # strategy="laplace", cutoff=1e-6,
  verbose=TRUE
)


# save as it takes so long
if (0) {
  fn16 = "~/tmp/car_annual16.rdz"  # 1.2GB
  read_write_fast( fit, file=fn16 )
  fit = read_write_fast(fn15)
}

s = summary(fit)
s$dic$dic  #   33552
s$dic$p.eff #  5903

plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_array(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( AUID=M$AUID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( AUID=sppoly$AUID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA

RES$INLA.Envir.AR1_strata.CAR = colSums( {out * weight_year * sppoly$au_sa_km2}[sppoly$strata_to_keep,], na.rm=TRUE )

lines( INLA.Envir.AR1_strata.CAR ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")

# map it
vn = "pred"
sppoly[,vn] = out[,"2017"]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )

plot( fit$marginals.hyperpar$"Phi for strata", type="l")  # posterior distribution of phi
plot( fit$marginals.hyperpar$"Precision for strata", type="l")
plot( fit$marginals.hyperpar$"Precision for setno", type="l")



plot( fit$marginals.hyperpar$"Phi for strata", type="l")  # posterior distribution of phi nonspatial dominates
plot( fit$marginals.hyperpar$"Precision for strata", type="l")
plot( fit$marginals.hyperpar$"Precision for setno", type="l")






# ------------------------------------------------
# Model 17:
# "INLA Envir AR1|strata CAR|year"	bym2|year	ar1|strata	rw2: temp+depth 	Poisson	5907	33552	1

# CAR with totno and poisson ...
# -- a separate CAR process for every year
# -- AR1 on year by each stratum
# 24+ hrs

fit = inla(
  formula = totno ~ 1
    + offset( log(data_offset) ) # CARSTM does log-transformation internally  but this is a direct call to inla 
    + f(iid_error, model="iid", hyper=H$iid)
    + f(ti, model="rw2", scale.model=TRUE, hyper=H$rw2)
    + f(zi, model="rw2", scale.model=TRUE, hyper=H$rw2)
    + f(year, model="ar1", hyper=H$ar1, group=strata )
    + f(strata, model="bym2", graph=slot(sppoly, "nb"), group=year, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
  family = "poisson",
  data=M,
  control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=TRUE, link =1 ), # compute=TRUE on each data location
  control.fixed=H$fixed,  # priors for fixed effects
  control.inla=list(  correct=TRUE, correct.verbose=FALSE ), # strategy="laplace", cutoff=1e-6,
  verbose=TRUE
)

# save as it takes so long  .. 24 hrs
if (0) {
  fn17 = "~/tmp/car_annual17.rdz"  # 1.2GB
  read_write_fast( fit, file=fn17 )
  fit = read_write_fast(fn17)
}

s = summary(fit)
s$dic$dic  #    33552
s$dic$p.eff #   5907

plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_array(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( AUID=M$AUID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( AUID=sppoly$AUID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA

RES$INLA.Envir.AR1_strata.CAR_year = colSums( {out * weight_year * sppoly$au_sa_km2}[sppoly$strata_to_keep,], na.rm=TRUE )

lines( INLA.Envir.AR1_strata.CAR_year ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")

# map it
vn = "pred"
sppoly[,vn] = out[,"2017"]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )

plot( fit$marginals.hyperpar$"Phi for strata", type="l")  # posterior distribution of phi
plot( fit$marginals.hyperpar$"Precision for strata", type="l")
plot( fit$marginals.hyperpar$"Precision for setno", type="l")




# ------------------------------------------------
# Model 18:
# "INLA Envir yr_iid CAR|year"	 bym2|year	yr_iid rw2: temp+depth 	Poisson	5943	33650	3

# CAR with totno and poisson ... a separate CAR process for every year .. and AR1

# 8+ hrs

fit = inla(
  formula = totno ~ 1
    + offset( log(data_offset) ) # CARSTM does log-transformation internally  but this is a direct call to inla 
    + f(iid_error, model="iid", hyper=H$iid)
    + f(ti, model="rw2", scale.model=TRUE, hyper=H$rw2)
    + f(zi, model="rw2", scale.model=TRUE, hyper=H$rw2)
    + f(year, model="iid", hyper=H$iid )
    + f(strata, model="bym2", graph=slot(sppoly, "nb"), group=year, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
  family = "poisson",
  data=M,
  control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=TRUE, link =1 ), # compute=TRUE on each data location
  control.fixed=H$fixed,  # priors for fixed effects
  control.inla=list(  correct=TRUE, correct.verbose=FALSE ), # strategy="laplace", cutoff=1e-6,
  verbose=TRUE
)

# save as it takes so long
if (0) {
  fn18 = "~/tmp/car_annual18.rdz"  # 1.2GB
  read_write_fast( fit, file=fn18 )
  fit = read_write_fast(fn18)
}

s = summary(fit)
s$dic$dic  #   33652
s$dic$p.eff #  5947

plot(fit, plot.prior=TRUE, plot.hyperparameters=TRUE, plot.fixed.effects=FALSE )


# reformat predictions into matrix form
ii = which(M$tag=="predictions")
out = reformat_to_array(
  input = fit$summary.fitted.values[ ii, "mean" ],
  matchfrom = list( AUID=M$AUID[ii], yr_factor=M$yr_factor[ii]),
  matchto   = list( AUID=sppoly$AUID, yr_factor=factor(p$yrs) )
)
# out[ out>1e10] = NA

RES$INLA.Envir.yr_iid.CAR_year = colSums( {out * weight_year * sppoly$au_sa_km2}[sppoly$strata_to_keep,], na.rm=TRUE )

lines( INLA.Envir.yr_iid.CAR_year ~ yr, data=RES, lty=1, lwd=2.5, col="blue", type="b")

# map it
vn = "pred"
sppoly[,vn] = out[,"2017"]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )

plot( fit$marginals.hyperpar$"Phi for strata", type="l")  # posterior distribution of phi
plot( fit$marginals.hyperpar$"Precision for strata", type="l")
plot( fit$marginals.hyperpar$"Precision for setno", type="l")



plot( fit$marginals.hyperpar$"Phi for strata", type="l")  # posterior distribution of phi nonspatial dominates
plot( fit$marginals.hyperpar$"Precision for strata", type="l")
plot( fit$marginals.hyperpar$"Precision for setno", type="l")






dev.new(width=11, height=7)
col = c("slategray", "turquoise", "darkorange", "green", "blue", "darkred", "cyan", "darkgreen", "purple" )
pch = c(20, 21, 22, 23, 24, 25, 26, 27, 20)
lty = c(1, 3, 4, 5, 6, 7, 1, 3, 4 )
lwd = c(4, 4, 4, 4, 4, 4, 4, 4, 4 )
type =c("l", "l", "l", "l", "l", "l", "l", "l", "l")
legend=c("Standard tow stratanal", "INLA Envir", "INLA Envir AR1", "INLA Envir CAR", "INLA Envir AR1 CAR", "INLA Envir AR1 CAR|year", "INLA Envir AR1|strata CAR", "INLA Envir AR1|strata CAR|year", "INLA Envir CAR|year")

plot( stratanal_towdistance  ~ yr, data=RES, lty=lty[1], lwd=lwd[1], col=col[1], pch=pch[1], type=type[1], ylim=c(0,0.46e9), xlab="Year", ylab="kg")
lines( INLA.Envir.1 ~ yr, data=RES, lty=lty[2], lwd=lwd[2], col=col[2], pch=pch[2], type=type[2])
lines( INLA.Envir.AR1.iid_year ~ yr, data=RES, lty=lty[3], lwd=lwd[3], col=col[3], pch=pch[3], type=type[3])
lines( INLA.Envir.CAR ~ yr, data=RES, lty=lty[4], lwd=lwd[4], col=col[4], pch=pch[4], type=type[4])  # yr_iid
lines( INLA.Envir.AR1.CAR ~ yr, data=RES, lty=lty[5], lwd=lwd[5], col=col[5], pch=pch[5], type=type[5])
lines( INLA.Envir.AR1.CAR_year ~ yr, data=RES, lty=lty[6], lwd=lwd[6], col=col[6], pch=pch[6], type=type[6])
lines( INLA.Envir.AR1_strata.CAR ~ yr, data=RES, lty=lty[7], lwd=lwd[7], col=col[7], pch=pch[7], type=type[7])
lines( INLA.Envir.AR1_strata.CAR_year ~ yr, data=RES, lty=lty[8], lwd=lwd[8], col=col[8], pch=pch[8], type=type[8])
lines( INLA.Envir.yr_iid.CAR_year ~ yr, data=RES, lty=lty[9], lwd=lwd[9], col=col[9], pch=pch[9], type=type[9])



legend("topright", legend=legend, lty=lty, col=col, lwd=lwd )



### end
