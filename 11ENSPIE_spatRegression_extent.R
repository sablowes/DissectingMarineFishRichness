##====================================================================================================================
##	linear and simultaneous autoregressive models (SAR)
##====================================================================================================================
rm(list=ls())
##============================================================================
##	run script that creates the spatial weights matrices (and also loads the data)
source('~/Dropbox/1current/dissectingRichness/revision/2CreateSpatWeights.R')
##============================================================================
##	latitude and longitude (No Atlantic Ocean data, Indo-Pacific only)
##	lm first 
summary(lENSpie_lat_long_4scale_lm <- lm(lENSpie ~ abs_lat + abs_long + I(abs_long^2) + I(abs_lat^2) + 
  scale + lextent + abs_long:scale + abs_lat:scale + abs_lat:abs_long, data=cleanTaxon_NoAtlantic_resamp200))
anova(lENSpie_lat_long_4scale_lm)

##	check for autocorrelation in the residuals...
#moran.test(residuals(lENSpie_lat_long_4scale_lm), d50_noAtlantic_spW, zero.policy=T)

##  use the same distance (50km) spatial weights matrix as main text results

##	SAR: distance-based neighbours (some have no neighbours) all neighbours within 50km
summary(lENSpie_lat_long_4scale_SAR_errW_d50 <- errorsarlm(lENSpie ~ abs_lat + abs_long + I(abs_long^2) + I(abs_lat^2) + 
  scale + lextent  + abs_long:scale + abs_lat:scale + abs_lat:abs_long, 
  data=cleanTaxon_NoAtlantic_resamp200, listw=d50_noAtlantic_spW, method='Matrix', zero.policy=T), Nagelkerke=T)
#moran.test(resid(lENSpie_lat_long_4scale_SAR_errW_d50),  listw= d50_noAtlantic_spW, zero.policy=T)

##	model selection by LR tests
##	abs_lat:abs_long 
summary(lENSpie_lat_long_4scale_SAR_errW_d50_1 <- update(lENSpie_lat_long_4scale_SAR_errW_d50, ~.-abs_lat:abs_long))
anova(lENSpie_lat_long_4scale_SAR_errW_d50_1, lENSpie_lat_long_4scale_SAR_errW_d50)

##	abs_lat:scale 
summary(lENSpie_lat_long_4scale_SAR_errW_d50_1a <- update(lENSpie_lat_long_4scale_SAR_errW_d50, ~.-abs_lat:scale))
anova(lENSpie_lat_long_4scale_SAR_errW_d50_1a, lENSpie_lat_long_4scale_SAR_errW_d50)

##	abs_long:scale
summary(lENSpie_lat_long_4scale_SAR_errW_d50_1b <- update(lENSpie_lat_long_4scale_SAR_errW_d50, ~.-abs_long:scale))
anova(lENSpie_lat_long_4scale_SAR_errW_d50_1b, lENSpie_lat_long_4scale_SAR_errW_d50)

##	abs_long^2 
summary(lENSpie_lat_long_4scale_SAR_errW_d50_2 <- update(lENSpie_lat_long_4scale_SAR_errW_d50_1a, ~.-I(abs_long^2)))
anova(lENSpie_lat_long_4scale_SAR_errW_d50_2, lENSpie_lat_long_4scale_SAR_errW_d50_1a)

##	abs_lat^2 
summary(lENSpie_lat_long_4scale_SAR_errW_d50_2a <- update(lENSpie_lat_long_4scale_SAR_errW_d50_1a, ~.-I(abs_lat^2)))
anova(lENSpie_lat_long_4scale_SAR_errW_d50_2a, lENSpie_lat_long_4scale_SAR_errW_d50_1a)

##	extent
summary(lENSpie_lat_long_4scale_SAR_errW_d50_3 <- update(lENSpie_lat_long_4scale_SAR_errW_d50_2, ~.-lextent))
anova(lENSpie_lat_long_4scale_SAR_errW_d50_3, lENSpie_lat_long_4scale_SAR_errW_d50_2)

##	is there autocorrelation in the residuals of the SARs?  
#moran.test(residuals(lENSpie_lat_long_4scale_SAR_errW_d50_1a), d50_noAtlantic_spW, zero.policy=T)

##====================================================================================================================
##	diagnostics
##	latitude and longitude
par(mfrow=c(4,2), mar=c(4,4,2,2))
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lENSpie_lat_long_4scale_SAR_errW_d50_3) ~ realm));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lENSpie_lat_long_4scale_SAR_errW_d50_3) ~ ecoregion));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lENSpie_lat_long_4scale_SAR_errW_d50_3) ~ scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lENSpie_lat_long_4scale_SAR_errW_d50_3) ~ abs_lat, col=scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lENSpie_lat_long_4scale_SAR_errW_d50_3) ~ abs_long, col=scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lENSpie_lat_long_4scale_SAR_errW_d50_3) ~ lextent, col=scale));abline(h=0, lty=2)
plot(fitted(lENSpie_lat_long_4scale_SAR_errW_d50_3), residuals(lENSpie_lat_long_4scale_SAR_errW_d50_3));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(fitted(lENSpie_lat_long_4scale_SAR_errW_d50_3), lENSpie, col=scale));abline(c(0,1), lty=2)