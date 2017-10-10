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
summary(lind_lat_long_4scale_lm <- lm(lind ~ abs_lat + abs_long + I(abs_long^2) + I(abs_lat^2) + scale + lextent + 
  abs_long:scale + abs_lat:scale + abs_lat:abs_long, data=cleanTaxon_NoAtlantic_resamp200))
anova(lind_lat_long_4scale_lm)

##	check for autocorrelation in the residuals...
#lm.morantest(lind_lat_long_4scale_lm, listw=d50_noAtlantic_spW, zero.policy=T)

##  use the same distance as results in main text
##	SAR: distance-based neighbours (some have no neighbours) all neighbours within 50km
summary(lind_lat_long_4scale_SAR_errW_d50 <- errorsarlm(lind ~ abs_lat + abs_long + I(abs_long^2) + I(abs_lat^2) + 
  scale + lextent + abs_long:scale + abs_lat:scale + abs_lat:abs_long, 
  data=cleanTaxon_NoAtlantic_resamp200, listw=d50_noAtlantic_spW, method='Matrix', zero.policy=T), Nagelkerke=T)
#moran.test(resid(lind_lat_long_4scale_SAR_errW_d50),  listw= d50_noAtlantic_spW, zero.policy=T)

##	model selection by LR tests
##	abs_lat:abs_long can GO 
summary(lind_lat_long_4scale_SAR_errW_d50_1 <- update(lind_lat_long_4scale_SAR_errW_d50, ~.-abs_lat:abs_long))
anova(lind_lat_long_4scale_SAR_errW_d50_1, lind_lat_long_4scale_SAR_errW_d50)

##	abs_lat:scale can go 
summary(lind_lat_long_4scale_SAR_errW_d50_1a <- update(lind_lat_long_4scale_SAR_errW_d50, ~.-abs_lat:scale))
anova(lind_lat_long_4scale_SAR_errW_d50_1a, lind_lat_long_4scale_SAR_errW_d50)

##	abs_long:scale 
summary(lind_lat_long_4scale_SAR_errW_d50_1b <- update(lind_lat_long_4scale_SAR_errW_d50, ~.-abs_long:scale))
anova(lind_lat_long_4scale_SAR_errW_d50_1b, lind_lat_long_4scale_SAR_errW_d50)

##	drop 2-way interactions 2 at a time
##	abs_long:scale and abs_lat:abs_long? 
summary(lind_lat_long_4scale_SAR_errW_d50_2 <- update(lind_lat_long_4scale_SAR_errW_d50, ~.-abs_long:scale - abs_lat:abs_long))
anova(lind_lat_long_4scale_SAR_errW_d50_2, lind_lat_long_4scale_SAR_errW_d50)

##	abs_lat:scale and abs_lat:abs_long? 
summary(lind_lat_long_4scale_SAR_errW_d50_2a <- update(lind_lat_long_4scale_SAR_errW_d50, ~.-abs_lat:scale - abs_lat:abs_long))
anova(lind_lat_long_4scale_SAR_errW_d50_2a, lind_lat_long_4scale_SAR_errW_d50)

##	abs_lat:scale and abs_long:scale? 
summary(lind_lat_long_4scale_SAR_errW_d50_2c <- update(lind_lat_long_4scale_SAR_errW_d50, ~.-abs_lat:scale - abs_long:scale))
anova(lind_lat_long_4scale_SAR_errW_d50_2c, lind_lat_long_4scale_SAR_errW_d50)

##	drop all three interactions? 
summary(lind_lat_long_4scale_SAR_errW_d50_3 <- update(lind_lat_long_4scale_SAR_errW_d50, ~.-abs_lat:abs_long - abs_lat:scale - abs_long:scale))
anova(lind_lat_long_4scale_SAR_errW_d50_3, lind_lat_long_4scale_SAR_errW_d50)

##	drop scale? 
summary(lind_lat_long_4scale_SAR_errW_d50_4 <- update(lind_lat_long_4scale_SAR_errW_d50_3, ~.-scale))
anova(lind_lat_long_4scale_SAR_errW_d50_4, lind_lat_long_4scale_SAR_errW_d50)

##	drop long^2? 
summary(lind_lat_long_4scale_SAR_errW_d50_4a <- update(lind_lat_long_4scale_SAR_errW_d50_3, ~.-I(abs_long^2)))
anova(lind_lat_long_4scale_SAR_errW_d50_4a, lind_lat_long_4scale_SAR_errW_d50_3)

##	drop lat^2? 
summary(lind_lat_long_4scale_SAR_errW_d50_4b <- update(lind_lat_long_4scale_SAR_errW_d50_3, ~.-I(abs_lat^2)))
anova(lind_lat_long_4scale_SAR_errW_d50_4b, lind_lat_long_4scale_SAR_errW_d50_3)

##	drop extent?
summary(lind_lat_long_4scale_SAR_errW_d50_4c <- update(lind_lat_long_4scale_SAR_errW_d50_3, ~.-lextent))
anova(lind_lat_long_4scale_SAR_errW_d50_4c, lind_lat_long_4scale_SAR_errW_d50)

##	is there autocorrelation in the residuals of the simplified SAR?  NO 
#moran.test(residuals(lind_lat_long_4scale_SAR_errW_d50_2), d50_noAtlantic_spW, zero.policy=T)

##====================================================================================================================
##	VISUAL diagnostics
##	latitude and longitude
par(mfrow=c(4,2), mar=c(4,4,2,2))
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lind_lat_long_4scale_SAR_errW_d50_3) ~ realm));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lind_lat_long_4scale_SAR_errW_d50_3) ~ ecoregion));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lind_lat_long_4scale_SAR_errW_d50_3) ~ scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lind_lat_long_4scale_SAR_errW_d50_3) ~ abs_lat, col=scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lind_lat_long_4scale_SAR_errW_d50_3) ~ abs_long, col=scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lind_lat_long_4scale_SAR_errW_d50_3) ~ lextent, col=scale));abline(h=0, lty=2)
plot(fitted(lind_lat_long_4scale_SAR_errW_d50_3), residuals(lind_lat_long_4scale_SAR_errW_d50_3));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(fitted(lind_lat_long_4scale_SAR_errW_d50_3), lind, col=scale));abline(c(0,1), lty=2)
