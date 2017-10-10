##====================================================================================================================
##	linear and simultaneous autoregressive models (SAR) fit to log(ENSPIE)
##====================================================================================================================
rm(list=ls())
##============================================================================
##	run script that creates the spatial weights matrices (and also loads the data)
source('~/Dropbox/1current/dissectingRichness/revision/2CreateSpatWeights_new.R')
ls()
##============================================================================
##	latitude and longitude (No Atlantic Ocean data, Indo-Pacific only)
##	lm first 
summary(lENSpie_lat_long_4scale_lm <- lm(lENSpie ~ abs_lat + abs_long + I(abs_long^2) + I(abs_lat^2) + scale + 
  abs_long:scale + abs_lat:scale + abs_lat:abs_long, data=cleanTaxon_NoAtlantic_resamp200))
anova(lENSpie_lat_long_4scale_lm)

##	check for autocorrelation in the residuals...all Moran tests of the linear models need to use the same spW (50km)
#lm.morantest(lENSpie_lat_long_4scale_lm, listw=d50_noAtlantic_spW, zero.policy=T)

##	SAR: distance-based neighbours (some have no neighbours) mean distance to nearest neighbour
summary(lENSpie_lat_long_4scale_SAR_errW_dmean <- errorsarlm(lENSpie ~ abs_lat + abs_long + 
  I(abs_long^2) + I(abs_lat^2) + scale + abs_long:scale + abs_lat:scale + abs_lat:abs_long, 
  data=cleanTaxon_NoAtlantic_resamp200, listw=dmean_noAtlantic_spW, method='Matrix', zero.policy=T), Nagelkerke=T)
#moran.test(resid(lENSpie_lat_long_4scale_SAR_errW_dmean),  listw= dmean_noAtlantic_spW, zero.policy=T)

##	SAR: distance-based neighbours (some have no neighbours) all neighbours within 50km
summary(lENSpie_lat_long_4scale_SAR_errW_d50 <- errorsarlm(lENSpie ~ abs_lat + abs_long + 
  I(abs_long^2) + I(abs_lat^2) + scale + abs_long:scale + abs_lat:scale + abs_lat:abs_long, 
  data=cleanTaxon_NoAtlantic_resamp200, listw=d50_noAtlantic_spW, method='Matrix', zero.policy=T), Nagelkerke=T)
#moran.test(resid(lENSpie_lat_long_4scale_SAR_errW_d50),  listw= d50_noAtlantic_spW, zero.policy=T)

##	SAR: distance-based neighbours (some have no neighbours) all neighbours within 100km
summary(lENSpie_lat_long_4scale_SAR_errW_d100 <- errorsarlm(lENSpie ~ abs_lat + abs_long + 
  I(abs_long^2) + I(abs_lat^2) + scale + abs_long:scale + abs_lat:scale + abs_lat:abs_long, 
  data=cleanTaxon_NoAtlantic_resamp200, listw=d100_noAtlantic_spW, method='Matrix', zero.policy=T), Nagelkerke=T)
#moran.test(resid(lENSpie_lat_long_4scale_SAR_errW_d100),  listw= d100_noAtlantic_spW, zero.policy=T)

##	SAR: d200
summary(lENSpie_lat_long_4scale_SAR_errW_d200 <- errorsarlm(lENSpie ~ abs_lat + abs_long + I(abs_long^2) + I(abs_lat^2) + scale + abs_long:scale + abs_lat:scale + abs_lat:abs_long, data=cleanTaxon_NoAtlantic_resamp200, listw=d200_noAtlantic_spW, method='Matrix', zero.policy=T), Nagelkerke=T)
#moran.test(resid(lENSpie_lat_long_4scale_SAR_errW_d200),  listw= d200_noAtlantic_spW, zero.policy=T)

##	what is the best neighbourhood matrix?  50km wins AGAIN!!
ENSpie_lat_long_neighbourhood_AIC <- bbmle::AICtab(lENSpie_lat_long_4scale_lm,
	lENSpie_lat_long_4scale_SAR_errW_dmean, 
	lENSpie_lat_long_4scale_SAR_errW_d50, 
	lENSpie_lat_long_4scale_SAR_errW_d100, 
	lENSpie_lat_long_4scale_SAR_errW_d200, base=T, weights=T)

##	model selection by LR tests
##	abs_lat:abs_long STAYS (p<0.001)
summary(lENSpie_lat_long_4scale_SAR_errW_d50_1 <- update(lENSpie_lat_long_4scale_SAR_errW_d50, ~.-abs_lat:abs_long))
anova(lENSpie_lat_long_4scale_SAR_errW_d50_1, lENSpie_lat_long_4scale_SAR_errW_d50)

##	abs_lat:scale can GO (p=0.16)
summary(lENSpie_lat_long_4scale_SAR_errW_d50_1a <- update(lENSpie_lat_long_4scale_SAR_errW_d50, ~.-abs_lat:scale))
anova(lENSpie_lat_long_4scale_SAR_errW_d50_1a, lENSpie_lat_long_4scale_SAR_errW_d50)

##	abs_long:scale STAYS (p<0.001)
summary(lENSpie_lat_long_4scale_SAR_errW_d50_1b <- update(lENSpie_lat_long_4scale_SAR_errW_d50, ~.-abs_long:scale))
anova(lENSpie_lat_long_4scale_SAR_errW_d50_1b, lENSpie_lat_long_4scale_SAR_errW_d50)

##	abs_long^2?
summary(lENSpie_lat_long_4scale_SAR_errW_d50_2 <- update(lENSpie_lat_long_4scale_SAR_errW_d50_1a, ~.-I(abs_long^2)))
anova(lENSpie_lat_long_4scale_SAR_errW_d50_2, lENSpie_lat_long_4scale_SAR_errW_d50_1a)

##	abs_lat^2
summary(lENSpie_lat_long_4scale_SAR_errW_d50_2a <- update(lENSpie_lat_long_4scale_SAR_errW_d50_1a, ~.-I(abs_lat^2)))
anova(lENSpie_lat_long_4scale_SAR_errW_d50_2a, lENSpie_lat_long_4scale_SAR_errW_d50_1a)

##	is there autocorrelation in the residuals of the SARs?
#moran.test(residuals(lENSpie_lat_long_4scale_SAR_errW_d50_1a), d50_noAtlantic_spW, zero.policy=T)

##====================================================================================================================
##	visual diagnostics
##	latitude and longitude
par(mfrow=c(4,2), mar=c(4,4,2,2))
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lENSpie_lat_long_4scale_SAR_errW_d50_2) ~ realm));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lENSpie_lat_long_4scale_SAR_errW_d50_2) ~ ecoregion));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lENSpie_lat_long_4scale_SAR_errW_d50_2) ~ scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lENSpie_lat_long_4scale_SAR_errW_d50_2) ~ abs_lat, col=scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lENSpie_lat_long_4scale_SAR_errW_d50_2) ~ abs_long, col=scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lENSpie_lat_long_4scale_SAR_errW_d50_2) ~ lextent, col=scale));abline(h=0, lty=2)
plot(fitted(lENSpie_lat_long_4scale_SAR_errW_d50_2), residuals(lENSpie_lat_long_4scale_SAR_errW_d50_2));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(fitted(lENSpie_lat_long_4scale_SAR_errW_d50_2), lENSpie, col=scale));abline(c(0,1), lty=2)

