##====================================================================================================================
##	linear and simultaneous autoregressive models (SAR) fit to log(species richness)
##====================================================================================================================
rm(list=ls())
##============================================================================
##	run script that creates the spatial weights matrices (and also loads the data)
source('~/Dropbox/1current/dissectingRichness/revision/2CreateSpatWeights.R')
##============================================================================
##	latitude and longitude (No Atlantic Ocean data, Indo-Pacific only)
##	lm first 
summary(lspp_rich_lat_long_4scale_lm <- lm(lspp_rich ~ abs_lat + abs_long + I(abs_long^2) + I(abs_lat^2) + scale + 
  abs_long:scale + abs_lat:scale + abs_lat:abs_long, data=cleanTaxon_NoAtlantic_resamp200))
anova(lspp_rich_lat_long_4scale_lm)

##	check for autocorrelation in the residuals...
#lm.morantest(lspp_rich_lat_long_4scale_lm, d50_noAtlantic_spW, zero.policy=T)

##	SAR: distance-based neighbours (some have no neighbours) mean distance to nearest neighbour
summary(lspp_rich_lat_long_4scale_SAR_errW_dmean <- errorsarlm(lspp_rich ~ abs_lat + abs_long + 
  I(abs_long^2) + I(abs_lat^2) + scale + abs_long:scale + abs_lat:scale + abs_lat:abs_long, 
  data=cleanTaxon_NoAtlantic_resamp200, listw=dmean_noAtlantic_spW, method='Matrix', zero.policy=T), Nagelkerke=T)
#moran.test(resid(lspp_rich_lat_long_4scale_SAR_errW_dmean),  listw= dmean_noAtlantic_spW, zero.policy=T)

##	SAR: distance-based neighbours (some have no neighbours) all neighbours within 50km
summary(lspp_rich_lat_long_4scale_SAR_errW_d50 <- errorsarlm(lspp_rich ~ abs_lat + abs_long + 
  I(abs_long^2) + I(abs_lat^2) + scale + abs_long:scale + abs_lat:scale + abs_lat:abs_long, 
  data=cleanTaxon_NoAtlantic_resamp200, listw=d50_noAtlantic_spW, method='Matrix', zero.policy=T), Nagelkerke=T)
#moran.test(resid(lspp_rich_lat_long_4scale_SAR_errW_d50),  listw= d50_noAtlantic_spW, zero.policy=T)

##	SAR: distance-based neighbours (some have no neighbours) all neighbours within 100km
summary(lspp_rich_lat_long_4scale_SAR_errW_d100 <- errorsarlm(lspp_rich ~ abs_lat + abs_long + 
  I(abs_long^2) + I(abs_lat^2) + scale + abs_long:scale + abs_lat:scale + abs_lat:abs_long, 
  data=cleanTaxon_NoAtlantic_resamp200, listw=d100_noAtlantic_spW, method='Matrix', zero.policy=T), Nagelkerke=T)
#moran.test(resid(lspp_rich_lat_long_4scale_SAR_errW_d100),  listw= d100_noAtlantic_spW, zero.policy=T)

##	SAR: d200
summary(lspp_rich_lat_long_4scale_SAR_errW_d200 <- errorsarlm(lspp_rich ~ abs_lat + abs_long + 
  I(abs_long^2) + I(abs_lat^2) + scale + abs_long:scale + abs_lat:scale + abs_lat:abs_long, 
  data=cleanTaxon_NoAtlantic_resamp200, listw=d200_noAtlantic_spW, method='Matrix', zero.policy=T), Nagelkerke=T)
#moran.test(resid(lspp_rich_lat_long_4scale_SAR_errW_d200),  listw= d200_noAtlantic_spW, zero.policy=T)

##	what is the best neighbourhood matrix?  50km wins AGAIN!!
lspp_rich_lat_long_neighbourhood_AIC <- bbmle::AICtab(lspp_rich_lat_long_4scale_lm, 
	lspp_rich_lat_long_4scale_SAR_errW_dmean, 
	lspp_rich_lat_long_4scale_SAR_errW_d50, 
	lspp_rich_lat_long_4scale_SAR_errW_d100, 
	lspp_rich_lat_long_4scale_SAR_errW_d200, base=T, weights=T)


##	LR tests for best spatial weights matrix
##	abs_lat:abs_long STAYS
summary(lspp_rich_lat_long_4scale_SAR_errW_d50_1 <- update(lspp_rich_lat_long_4scale_SAR_errW_d50, ~.-abs_lat:abs_long))
anova(lspp_rich_lat_long_4scale_SAR_errW_d50_1, lspp_rich_lat_long_4scale_SAR_errW_d50)

##	abs_lat:scale GOES 
summary(lspp_rich_lat_long_4scale_SAR_errW_d50_1a <- update(lspp_rich_lat_long_4scale_SAR_errW_d50, ~.-abs_lat:scale))
anova(lspp_rich_lat_long_4scale_SAR_errW_d50_1a, lspp_rich_lat_long_4scale_SAR_errW_d50)

##	abs_long:scale STAYS 
summary(lspp_rich_lat_long_4scale_SAR_errW_d50_1b <- update(lspp_rich_lat_long_4scale_SAR_errW_d50, ~.-abs_long:scale))
anova(lspp_rich_lat_long_4scale_SAR_errW_d50_1b, lspp_rich_lat_long_4scale_SAR_errW_d50)

##	abs_long^2 STAYS
summary(lspp_rich_lat_long_4scale_SAR_errW_d50_2 <- update(lspp_rich_lat_long_4scale_SAR_errW_d50_1a, ~.-I(abs_long^2)))
anova(lspp_rich_lat_long_4scale_SAR_errW_d50_2, lspp_rich_lat_long_4scale_SAR_errW_d50_1a)

##	abs_lat^2 STAYS
summary(lspp_rich_lat_long_4scale_SAR_errW_d50_3 <- update(lspp_rich_lat_long_4scale_SAR_errW_d50_1a, ~.-I(abs_lat^2)))
anova(lspp_rich_lat_long_4scale_SAR_errW_d50_3, lspp_rich_lat_long_4scale_SAR_errW_d50_1a)


##	is there autocorrelation in the residuals of the simplified SAR? 
#moran.test(residuals(lspp_rich_lat_long_4scale_SAR_errW_d50_d50_1a), d50_noAtlantic_spW, zero.policy=T)
##====================================================================================================================
##	visual diagnostics
par(mfrow=c(4,2), mar=c(4,4,2,2))
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lspp_rich_lat_long_4scale_SAR_errW_d50_1a) ~ realm));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lspp_rich_lat_long_4scale_SAR_errW_d50_1a) ~ ecoregion));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lspp_rich_lat_long_4scale_SAR_errW_d50_1a) ~ scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lspp_rich_lat_long_4scale_SAR_errW_d50_1a) ~ abs_lat, col=scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lspp_rich_lat_long_4scale_SAR_errW_d50_1a) ~ abs_long, col=scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lspp_rich_lat_long_4scale_SAR_errW_d50_1a) ~ lextent, col=scale));abline(h=0, lty=2)
plot(fitted(lspp_rich_lat_long_4scale_SAR_errW_d50_1a), residuals(lspp_rich_lat_long_4scale_SAR_errW_d50_1a));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(fitted(lspp_rich_lat_long_4scale_SAR_errW_d50_1a), lspp_rich, col=scale));abline(c(0,1), lty=2)

