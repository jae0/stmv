

# ------------------------------------------------
# Atlantic cod habitat  --- Areal unit modelling of habitat
# Uses carstm and aegis to develop models and predict


# ------------------------------------------------
# load data common environment and parameter setting
# source( system.file( "scripts", "00_cod_comparisons_data_environment.R", package = "carstm") )


# --------------------------------
# construct basic parameter list defining the main characteristics of the study
# and some plotting parameters (bounding box, projection, bathymetry layout, coastline)
# NOTE: the data selection is the same as in (01_cod_comparisons_basic_stranal.R)
require(carstm)
require(aegis.polygons)

p =list(
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
set$time = set$time_space = as.numeric( set$yr_factor)

set$iid_error = 1:nrow(set) # for inla indexing
set$tag = "observations"


## --------------------------------
# construct meanweights matrix
weight_year = meanweights_by_arealunit( set=set, AUID=as.character( sppoly$AUID ), yrs=p$yrs, fillall=TRUE, annual_breakdown=TRUE )
# weight_year = meanweights_by_arealunit_modelled( p=p, redo=TRUE )  -- note: data passing of M needs to be modularized 
# weight_year = weight_year[, match(as.character(p$yrs), colnames(weight_year) )]
# weight_year = weight_year[ match(as.character(sppoly$AUID), rownames(weight_year) )]



# RES = data.frame(yr=p$selection$survey[["yr"]]) # collect model comparisons
if (0) {
  fn = file.path( project.datadirectory( "carstm" ), "RES.rdz" )
 # read_write_fast(RES, file=fn)
 # RES = read_write_fast(fn)
}



## ----------------------------------
# covariate estimates for prediction in strata and year

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

pa = presence.absence( X={set$totno / set$data_offset}, px=0.05 )  # determine presence absence and weighting
set[, "pa"] = pa$pa
set[, "wt"] = pa$probs
pa = NULL


# ------------------------------------------------
#  good data
ok = which(
  is.finite(set$totno) &   # INLA can impute Y-data
  is.finite(set$t) &
  is.finite(set$z) &
  is.finite(set$data_offset) &
  set$AUID %in% sppoly$AUID[sppoly$strata_to_keep]
)


# ---------------------
# generic PC priors
H = inla_hyperparameters( sd(set$pa) )



# ------------------------------------------------
# Model 1: a binomial (presence-absence) aks, habitat probability model with linear covariate effects
fit = glm(
  formula = pa ~  1 + t + z +degreedays,
  family=binomial(link="logit"),
  data=set
)

if (0) {
  fnx = file.path( project.datadirectory( "carstm" ), "fit_habitat_glm.rdz" )
  read_write_fast( fit, file=fnx )
  # fit = read_write_fast(fnx)
}

s = summary(fit)
AIC(fit)  # 10774

toplot = expand.grid( t=seq(-1,20), z=c(5,10,20,40,80,160,320,640),degreedays=seq(0, 5000, by=100) )
toplot$predictions = predict(fit, newdata=toplot, type="response", se.fit=FALSE  )

plot( predictions ~ z, toplot[ which( {toplot$t==min(toplot$t)} & {toplot$degreedays==min(toplot$degreedays)} ), ], type="b" )
plot( predictions ~ t, toplot[ which( {toplot$z==min(toplot$z)} & {toplot$degreedays==min(toplot$degreedays)} ), ], type="b" )
plot( predictions ~ degreedays, toplot[ which( {toplot$z==min(toplot$z)} & {toplot$t==min(toplot$t)} ), ], type="b")


# predicted probabilities of observing cod, given covariates (temperature, depth, etc)

APS = aegis_prediction_surface( aegis_data=res$means )
APS$data_offset=1
APS$yr = APS$year
APS$yr_factor = factor( APS$year, levels=p$yrs)
APS$iyr = match(APS$yr_factor, p$yrs)
APS$istrata = match( APS$AUID, sppoly$AUID )

predictions = predict(fit, APS, type="response", se.fit=TRUE  )
APS$predictions = predictions$fit
APS$predictions.se = predictions$se.fit

# reformat predictions into matrix form
out = matrix(NA, nrow=length(sppoly$AUID), ncol=length(p$yrs), dimnames=list( sppoly$AUID, p$yrs) )
out[ cbind(APS$istrata, APS$iyr) ] = APS$predictions
RES$habitat_glm = colSums( {out * sppoly$au_sa_km2 }[sppoly$strata_to_keep,], na.rm=TRUE ) /sum(sppoly$au_sa_km2[sppoly$strata_to_keep]) # sa weighted average prob habitat


# map it onto strata means of temperature and depth
aps = APS[ APS$year==2017,  ]
iy = match( as.character(sppoly$AUID), aps$AUID )
vn = "pred"
sppoly[,vn] = APS$predictions[iy]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )




# ------------------------------------------------
# Model 5b: habitat model with a smoothed covariate effect
fit = mgcv::gam(
  formula = pa ~  1 + s(t, k=3, bs="ts")  +  s(z, k=3, bs="ts")  +  s(degreedays, k=3, bs="ts"),
  family=binomial(link="logit"),
  data=set
)


inverse.logit = function( x ) {
  # x should be the log odds ratio
  oddsratio = exp(x)
  prob = oddsratio / (1 + oddsratio )
  return (prob)
}

if (0) {
  fnx = file.path( project.datadirectory( "carstm" ), "fit_habitat_gam.rdz" )
  read_write_fast( fit, file=fnx )
  # fit = read_write_fast(fnx)
}


plot(fit, all.terms=TRUE, trans=inverse.logit, seWithMean=TRUE, jit=TRUE, rug=TRUE )
s = summary(fit)
AIC(fit)  #  10579 .. a little better than the glm

toplot = expand.grid( t=seq(-1,20), z=c(5,10,20,40,80,160,320,640),degreedays=seq(0, 5000, by=100) )
toplot$predictions = predict(fit, newdata=toplot, type="response", se.fit=FALSE  )

plot( predictions ~ z, toplot[ which( {toplot$t==min(toplot$t)} & {toplot$degreedays==min(toplot$degreedays)} ), ], type="b" )
plot( predictions ~ t, toplot[ which( {toplot$z==min(toplot$z)} & {toplot$degreedays==min(toplot$degreedays)} ), ], type="b" )
plot( predictions ~ degreedays, toplot[ which( {toplot$z==min(toplot$z)} & {toplot$t==min(toplot$t)} ), ], type="b")

# predicted probabilities of observing cod, given temperature and depth

APS = aegis_prediction_surface( aegis_data=res$means )
APS$data_offset=1
APS$yr = APS$year
APS$yr_factor = factor( APS$year, levels=p$yrs)
APS$iyr = match(APS$yr_factor, p$yrs)
APS$istrata = match( APS$AUID, sppoly$AUID )

predictions = predict(fit, APS, type="response", se.fit=TRUE  )
APS$predictions = predictions$fit
APS$predictions.se = predictions$se.fit

# reformat predictions into matrix form
out = matrix(NA, nrow=length(sppoly$AUID), ncol=length(p$yrs), dimnames=list( sppoly$AUID, p$yrs) )
out[ cbind(APS$istrata, APS$iyr) ] = APS$predictions

RES$habitat_gam = colSums( {out * sppoly$au_sa_km2 }[sppoly$strata_to_keep,], na.rm=TRUE ) /sum(sppoly$au_sa_km2[sppoly$strata_to_keep]) # sa weighted average prob habitat


# map it
iy = which( APS$year=="2017")
it = match( as.character(sppoly$AUID), APS$AUID[iy] )
vn = "pred"
sppoly[,vn] = APS$predictions[iy][it]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )


# ------------------------------------------------
# Model 5c: as above but via INLA  ... very slow
## similar to GAM

APS = aegis_prediction_surface( aegis_data=res$means )
APS$data_offset=1
APS$pa = NA  # what we are trying to predict
APS$tag = "predictions"
APS$yr = APS$year
APS$yr_factor = factor( APS$year, levels=p$yrs)

basic_vars = unique(c( "pa", "t", "z", "degreedays", "data_offset", "tag", "yr", "AUID"))

M = rbind( set[, basic_vars], APS[,basic_vars] )

M$t[!is.finite(M$t)] = median(M$t, na.rm=TRUE )  # missing data .. quick fix .. do something better for
M$z[!is.finite(M$z)] = median(M$z, na.rm=TRUE )  # missing data .. quick fix .. do something better for

M$yr_factor = factor( M$yr, levels=p$yrs )
M$AUID  = factor( M$AUID, levels=levels(sppoly$AUID ))

M$ti = discretize_data( x=M$t, brks=p$discretization$t )
M$zi = discretize_data( x=M$z, brks=p$discretization$z )
M$di = discretize_data( x=M$degreedays, brks=p$discretization$degreedays )


fit = inla(
  formula = pa ~ 1
      + f(ti, model="rw2", scale.model=TRUE, diagonal=1e-5, hyper=H$prec)
      + f(zi, model="rw2", scale.model=TRUE, diagonal=1e-5, hyper=H$prec)
      + f(di, model="rw2", scale.model=TRUE, diagonal=1e-5, hyper=H$prec),
  family="binomial",  # alternates family="zeroinflatedbinomial0", family="zeroinflatedbinomial1",
  data=M,
  control.family=list(control.link=list(model="logit")),
  control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=TRUE, link=1 ), # compute=TRUE on each data location
  control.fixed=H$fixed,  # priors for fixed effects
  control.inla=list(  correct=TRUE, correct.verbose=FALSE ), # strategy="laplace", cutoff=1e-6,
  verbose=TRUE
)

if (0) {
  fnx = file.path( project.datadirectory( "carstm" ), "fit_habitat_inla.rdz" )
  read_write_fast( fit, file=fnx )
  # fit = read_write_fast(fnx)
}

plot(fit )
s = summary(fit)
s$dic$dic  # 10585   .. not sure why ..
s$dic$p.eff # 17.73

APS = cbind( APS, fit$summary.fitted.values[ which(M$tag=="predictions"), ] )
APS$iyr = match(APS$yr_factor, p$yrs)
APS$istrata = match( APS$AUID, sppoly$AUID )

# reformat predictions into matrix form
out = matrix(NA, nrow=length(sppoly$AUID), ncol=length(p$yrs), dimnames=list( sppoly$AUID, p$yrs) )
out[ cbind(APS$istrata, APS$iyr) ] = APS$mean
RES$habitat_inla = colSums( {out * sppoly$au_sa_km2 }[sppoly$strata_to_keep,], na.rm=TRUE ) /sum(sppoly$au_sa_km2[sppoly$strata_to_keep]) # sa weighted average prob habitat

# map it
vn = "pred"
sppoly[,vn] = out[,"2017"]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )



# NOTE most of the decline occured when "habitat" was stable !
# NOTE habitat has become more variable since 1995
# NOTE habitat has declined in 2017
# to do : compute SE's and add to the graph



## -------------------------------------------------------
# bym, iid_year
# Model XX:


APS = aegis_prediction_surface( aegis_data=res$means )
APS$data_offset=1
APS$pa = NA  # what we are trying to predict
APS$tag = "predictions"
APS$yr = APS$year
APS$yr_factor = factor( APS$year, levels=p$yrs)

basic_vars = unique(c( "pa", "t", "z", "degreedays", "data_offset", "tag", "yr", "AUID"))

M = rbind( set[, basic_vars], APS[,basic_vars] )

M$t[!is.finite(M$t)] = median(M$t, na.rm=TRUE )  # missing data .. quick fix .. do something better for
M$z[!is.finite(M$z)] = median(M$z, na.rm=TRUE )  # missing data .. quick fix .. do something better for

M$iid_error = 1:nrow(M)
M$yr_factor = factor( as.character(M$yr) )
M$AUID  = factor( M$AUID, levels=levels(sppoly$AUID ))
M$strata  = as.numeric( M$AUID)
M$year  = as.numeric( M$yr_factor)


M$ti = discretize_data( x=M$t, brks=p$discretization$t )
M$zi = discretize_data( x=M$z, brks=p$discretization$z )
M$di = discretize_data( x=M$degreedays, brks=p$discretization$degreedays )

fit = inla(
  formula = pa ~ 1
  + f(iid_error, model="iid", hyper=H$iid)
  + f(ti, model="rw2", scale.model=TRUE, hyper=H$rw2)
  + f(zi, model="rw2", scale.model=TRUE, hyper=H$rw2)
  + f(di, model="rw2", scale.model=TRUE, hyper=H$rw2)
  + f(year, model="iid", hyper=H$iid)
  + f(strata, model="bym2", graph=slot(sppoly, "nb"), scale.model=TRUE, constr=TRUE, hyper=H$bym2),
  family="binomial",  # alternates family="zeroinflatedbinomial0", family="zeroinflatedbinomial1",
  data=M,
  control.family=list(control.link=list(model="logit")),
  control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=TRUE, link=1 ), # compute=TRUE on each data location
  control.fixed=H$fixed,  # priors for fixed effects
  control.inla=list(  correct=TRUE, correct.verbose=FALSE ), # strategy="laplace", cutoff=1e-6,
  verbose=TRUE
)

if (0) {
  fnx = file.path( project.datadirectory( "carstm" ), "fit_habitat_strata_CAR.yr_iid.rdz" )
  read_write_fast( fit, file=fnx )
  # fit = read_write_fast(fnx)
}

plot(fit )
s = summary(fit)
s$dic$dic  #8885
s$dic$p.eff #163.8

APS = cbind( APS, fit$summary.fitted.values[ which(M$tag=="predictions"), ] )

APS$iyr = match(APS$yr_factor, p$yrs)
APS$istrata = match( APS$AUID, sppoly$AUID )

# reformat predictions into matrix form
out = matrix(NA, nrow=length(sppoly$AUID), ncol=length(p$yrs), dimnames=list( sppoly$AUID, p$yrs) )
out[ cbind(APS$istrata, APS$iyr) ] = APS$mean
RES$habitat_strata_CAR.yr_iid = colSums( {out * sppoly$au_sa_km2 }[sppoly$strata_to_keep,], na.rm=TRUE ) /sum(sppoly$au_sa_km2[sppoly$strata_to_keep]) # sa weighted average prob habitat

# map it
vn = "pred"
sppoly[,vn] = out[,"2017"]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )


## -------------------------------------------------------
# bym, iid_year
# Model XX:


APS = aegis_prediction_surface( aegis_data=res$means )
APS$data_offset=1
APS$pa = NA  # what we are trying to predict
APS$tag = "predictions"
APS$yr = APS$year
APS$yr_factor = factor( APS$year, levels=p$yrs)

basic_vars = unique(c( "pa", "t", "z", "degreedays", "data_offset", "tag", "yr", "AUID"))

M = rbind( set[, basic_vars], APS[,basic_vars] )

M$t[!is.finite(M$t)] = median(M$t, na.rm=TRUE )  # missing data .. quick fix .. do something better for
M$z[!is.finite(M$z)] = median(M$z, na.rm=TRUE )  # missing data .. quick fix .. do something better for

M$iid_error = 1:nrow(M)
M$yr_factor = factor( as.character(M$yr) )
M$AUID  = factor( M$AUID, levels=levels(sppoly$AUID ))
M$strata  = as.numeric( M$AUID)
M$year  = as.numeric( M$yr_factor)


M$ti = discretize_data( x=M$t, brks=p$discretization$t )
M$zi = discretize_data( x=M$z, brks=p$discretization$z )
M$di = discretize_data( x=M$degreedays, brks=p$discretization$degreedays )


# 28404.251 seconds
fit = inla(
  formula = pa ~ 1
  + f(iid_error, model="iid", hyper=H$iid)
  + f(ti, model="rw2", scale.model=TRUE, hyper=H$rw2)
  + f(zi, model="rw2", scale.model=TRUE, hyper=H$rw2)
  + f(di, model="rw2", scale.model=TRUE, hyper=H$rw2)
  + f(year, model="iid", hyper=H$iid)
  + f(strata, model="bym2", graph=slot(sppoly, "nb"), group=year, scale.model=TRUE, constr=TRUE, hyper=H$bym2),
  family="binomial",  # alternates family="zeroinflatedbinomial0", family="zeroinflatedbinomial1",
  data=M,
  control.family=list(control.link=list(model="logit")),
  control.compute=list(cpo=TRUE, waic=TRUE, dic=TRUE, config=TRUE),
  control.results=list(return.marginals.random=TRUE, return.marginals.predictor=TRUE ),
  control.predictor=list(compute=TRUE, link=1 ), # compute=TRUE on each data location
  control.fixed=H$fixed,  # priors for fixed effects
  control.inla=list(  correct=TRUE, correct.verbose=FALSE ), # strategy="laplace", cutoff=1e-6,
  verbose=TRUE
)


if (0) {
  fnx = file.path( project.datadirectory( "carstm" ), "fit_habitat_strata_CAR_yr.yr_iid.rdz" )
  read_write_fast( fit, file=fnx )
  # fit = read_write_fast(fnx)
}


plot(fit )
s = summary(fit)
s$dic$dic  # 8730
s$dic$p.eff # 399.4

APS = cbind( APS, fit$summary.fitted.values[ which(M$tag=="predictions"), ] )

APS$iyr = match(APS$yr_factor, p$yrs)
APS$istrata = match( APS$AUID, sppoly$AUID )

# reformat predictions into matrix form
out = matrix(NA, nrow=length(sppoly$AUID), ncol=length(p$yrs), dimnames=list( sppoly$AUID, p$yrs) )
out[ cbind(APS$istrata, APS$iyr) ] = APS$mean
RES$habitat_strata_CAR_yr.yr_iid = colSums( {out * sppoly$au_sa_km2 }[sppoly$strata_to_keep,], na.rm=TRUE ) /sum(sppoly$au_sa_km2[sppoly$strata_to_keep]) # sa weighted average prob habitat

# map it
vn = "pred"
sppoly[,vn] = out[,"2017"]
brks = interval_break(X= sppoly[[vn]], n=length(p$mypalette), style="quantile")
spplot( sppoly, vn, col.regions=p$mypalette, main=vn, at=brks, sp.layout=p$coastLayout, col="transparent" )

# fn2="habitat_strata_CAR_yr.yr_iid.rdz"
# read_write_fast(fit, file=fn2)
# fit = read_write_fast( fn2)

######

if (0) {
 fn = file.path( project.datadirectory( "carstm" ), "RES.rdz" )
 # read_write_fast(RES, file=fn)
 # RES = read_write_fast(fn)
}

# RES$habitat_strata_CAR.yr_iid  = RES$habitat_bym_yriid

dev.new(width=11, height=7)
# col = c("slategray", "turquoise", "darkorange", "green", "blue", "darkred", "cyan", "darkgreen", "purple" )
col = c( "darkorange", "green", "blue",  "cyan", "darkgreen" )
pch = c(20, 21, 22, 23, 24)#, 25, 26, 27, 20)
lty = c(1, 3, 4, 5, 6)#, 7, 1, 3, 4 )
lwd = c(4, 4, 4, 4, 4)#, 4, 4, 4, 4 )
type =c("l", "l", "l", "l", "l")#, "l", "l", "l", "l")
legend=c("GLM", "GAM", "INLA", "INLA CAR", "INLA CAR|year")#, "INLA Envir AR1 CAR|year", "INLA Envir AR1|strata CAR", "INLA Envir AR1|strata CAR|year", "INLA Envir CAR|year")

plot( habitat_glm  ~ yr, data=RES, lty=lty[1], lwd=lwd[1], col=col[1], pch=pch[1], type=type[1], ylim=c(0.375,0.825), xlab="Year", ylab="Probability")
lines( habitat_gam ~ yr, data=RES, lty=lty[2], lwd=lwd[2], col=col[2], pch=pch[2], type=type[2])
lines( habitat_inla ~ yr, data=RES, lty=lty[3], lwd=lwd[3], col=col[3], pch=pch[3], type=type[3])
lines( habitat_strata_CAR.yr_iid ~ yr, data=RES, lty=lty[4], lwd=lwd[4], col=col[4], pch=pch[4], type=type[4])  # yr_iid
lines( habitat_strata_CAR_yr.yr_iid ~ yr, data=RES, lty=lty[5], lwd=lwd[5], col=col[5], pch=pch[5], type=type[5])
# lines( INLA.Envir.AR1.CAR_year ~ yr, data=RES, lty=lty[6], lwd=lwd[6], col=col[6], pch=pch[6], type=type[6])
# lines( INLA.Envir.AR1_strata.CAR ~ yr, data=RES, lty=lty[7], lwd=lwd[7], col=col[7], pch=pch[7], type=type[7])
# lines( INLA.Envir.AR1_strata.CAR_year ~ yr, data=RES, lty=lty[8], lwd=lwd[8], col=col[8], pch=pch[8], type=type[8])
# lines( INLA.Envir.yr_iid.CAR_year ~ yr, data=RES, lty=lty[9], lwd=lwd[9], col=col[9], pch=pch[9], type=type[9])



legend("topright", legend=legend, lty=lty, col=col, lwd=lwd )


# end
