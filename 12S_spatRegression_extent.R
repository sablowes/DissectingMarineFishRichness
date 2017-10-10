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
summary(lspp_rich_lat_long_4scale_lm <- lm(lspp_rich ~ abs_lat + abs_long + I(abs_long^2) + I(abs_lat^2) +
  scale + lextent + abs_long:scale + abs_lat:scale + abs_lat:abs_long, data=cleanTaxon_NoAtlantic_resamp200))
anova(lspp_rich_lat_long_4scale_lm)

##  use the same distance (50km) spatial weights matrix as main text results
##	SAR: distance-based neighbours (some have no neighbours) all neighbours within 50km
summary(lspp_rich_lat_long_4scale_SAR_errW_d50 <- errorsarlm(lspp_rich ~ abs_lat + abs_long + 
  I(abs_long^2) + I(abs_lat^2) + scale + lextent  + abs_long:scale + abs_lat:scale + abs_lat:abs_long, 
  data=cleanTaxon_NoAtlantic_resamp200, listw=d50_noAtlantic_spW, method='Matrix', zero.policy=T), Nagelkerke=T)
#moran.test(resid(lspp_rich_lat_long_4scale_SAR_errW_d50),  listw= d50_noAtlantic_spW, zero.policy=T)

##	LR tests for best spatial weights matrix
##	abs_lat:abs_long 
summary(lspp_rich_lat_long_4scale_SAR_errW_d50_1 <- update(lspp_rich_lat_long_4scale_SAR_errW_d50, ~.-abs_lat:abs_long))
anova(lspp_rich_lat_long_4scale_SAR_errW_d50_1, lspp_rich_lat_long_4scale_SAR_errW_d50)

##	abs_lat:scale 
summary(lspp_rich_lat_long_4scale_SAR_errW_d50_1a <- update(lspp_rich_lat_long_4scale_SAR_errW_d50, ~.-abs_lat:scale))
anova(lspp_rich_lat_long_4scale_SAR_errW_d50_1a, lspp_rich_lat_long_4scale_SAR_errW_d50)

##	abs_long:scale
summary(lspp_rich_lat_long_4scale_SAR_errW_d50_1b <- update(lspp_rich_lat_long_4scale_SAR_errW_d50, ~.-abs_long:scale))
anova(lspp_rich_lat_long_4scale_SAR_errW_d50_1b, lspp_rich_lat_long_4scale_SAR_errW_d50)

##	abs_long^2
summary(lspp_rich_lat_long_4scale_SAR_errW_d50_2 <- update(lspp_rich_lat_long_4scale_SAR_errW_d50_1a, ~.-I(abs_long^2)))
anova(lspp_rich_lat_long_4scale_SAR_errW_d50_2, lspp_rich_lat_long_4scale_SAR_errW_d50_1a)

##	abs_lat^2
summary(lspp_rich_lat_long_4scale_SAR_errW_d50_2a <- update(lspp_rich_lat_long_4scale_SAR_errW_d50_1a, ~.-I(abs_lat^2)))
anova(lspp_rich_lat_long_4scale_SAR_errW_d50_2a, lspp_rich_lat_long_4scale_SAR_errW_d50_1a)

##	extent?
summary(lspp_rich_lat_long_4scale_SAR_errW_d50_2b <- update(lspp_rich_lat_long_4scale_SAR_errW_d50_1a, ~.-lextent))
anova(lspp_rich_lat_long_4scale_SAR_errW_d50_2b, lspp_rich_lat_long_4scale_SAR_errW_d50_1a)

##	is there autocorrelation in the residuals of the simplified SAR?  NO :-)
#moran.test(residuals(lspp_rich_lat_long_4scale_SAR_errW_d50_1a), d50_noAtlantic_spW, zero.policy=T)
##====================================================================================================================
##	visual diagnostics
##	latitude and longitude
par(mfrow=c(4,2), mar=c(4,4,2,2))
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lspp_rich_lat_long_4scale_SAR_errW_d50_1a) ~ realm));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lspp_rich_lat_long_4scale_SAR_errW_d50_1a) ~ ecoregion));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lspp_rich_lat_long_4scale_SAR_errW_d50_1a) ~ scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lspp_rich_lat_long_4scale_SAR_errW_d50_1a) ~ abs_lat, col=scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lspp_rich_lat_long_4scale_SAR_errW_d50_1a) ~ abs_long, col=scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lspp_rich_lat_long_4scale_SAR_errW_d50_1a) ~ lextent, col=scale));abline(h=0, lty=2)
plot(fitted(lspp_rich_lat_long_4scale_SAR_errW_d50_1a), residuals(lspp_rich_lat_long_4scale_SAR_errW_d50_1a));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(fitted(lspp_rich_lat_long_4scale_SAR_errW_d50_1a), lspp_rich, col=scale));abline(c(0,1), lty=2)
