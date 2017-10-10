##====================================================================================================================
##	linear and simultaneous autoregressive models (SAR) to log(S_asymptote) i.e., Chao estimated richness
##====================================================================================================================
rm(list=ls())
##============================================================================
##	run script that creates the spatial weights matrices (and also loads the data)
source('~/Dropbox/1current/dissectingRichness/revision/2CreateSpatWeights.R')
##============================================================================
##============================================================================
##	latitude and longitude (No Atlantic Ocean data, Indo-Pacific only)
##	lm first 
summary(lChao_lat_long_4scale_lm <- lm(lChao ~ abs_lat + abs_long + 
  I(abs_long^2) + I(abs_lat^2) + scale + abs_long:scale + abs_lat:scale + abs_lat:abs_long, 
  data=cleanTaxon_NoAtlantic_resamp200))
anova(lChao_lat_long_4scale_lm)

##	check for autocorrelation in the residuals...
#lm.morantest(lChao_lat_long_4scale_lm, listw=d50_noAtlantic_spW, zero.policy=T)

##	SAR: distance-based neighbours (some have no neighbours) mean distance to nearest neighbour
summary(lChao_lat_long_4scale_SAR_errW_dmean <- errorsarlm(lChao ~ abs_lat + abs_long + I(abs_long^2) + I(abs_lat^2) + scale + 
  abs_long:scale + abs_lat:scale + abs_lat:abs_long, 
  data=cleanTaxon_NoAtlantic_resamp200, listw=dmean_noAtlantic_spW, method='Matrix', zero.policy=T), Nagelkerke=T)
#moran.test(resid(lChao_lat_long_4scale_SAR_errW_dmean),  listw= dmean_noAtlantic_spW, zero.policy=T)

##	SAR: distance-based neighbours (some have no neighbours) all neighbours within 50km
summary(lChao_lat_long_4scale_SAR_errW_d50 <- errorsarlm(lChao ~ abs_lat + abs_long + I(abs_long^2) + I(abs_lat^2) + scale + 
  abs_long:scale + abs_lat:scale + abs_lat:abs_long, 
  data=cleanTaxon_NoAtlantic_resamp200, listw=d50_noAtlantic_spW, method='Matrix', zero.policy=T), Nagelkerke=T)
#moran.test(resid(lChao_lat_long_4scale_SAR_errW_d50),  listw= d50_noAtlantic_spW, zero.policy=T)

##	SAR: distance-based neighbours (some have no neighbours) all neighbours within 100km
summary(lChao_lat_long_4scale_SAR_errW_d100 <- errorsarlm(lChao ~ abs_lat + abs_long + I(abs_long^2) + I(abs_lat^2) + scale + 
  abs_long:scale + abs_lat:scale + abs_lat:abs_long, 
  data=cleanTaxon_NoAtlantic_resamp200, listw=d100_noAtlantic_spW, method='Matrix', zero.policy=T), Nagelkerke=T)
#moran.test(resid(lChao_lat_long_4scale_SAR_errW_d100),  listw= d100_noAtlantic_spW, zero.policy=T)

##	SAR: d200
summary(lChao_lat_long_4scale_SAR_errW_d200 <- errorsarlm(lChao ~ abs_lat + abs_long + I(abs_long^2) + I(abs_lat^2) + scale + 
  abs_long:scale + abs_lat:scale + abs_lat:abs_long, 
  data=cleanTaxon_NoAtlantic_resamp200, listw=d200_noAtlantic_spW, method='Matrix', zero.policy=T), Nagelkerke=T)
#moran.test(resid(lChao_lat_long_4scale_SAR_errW_d200),  listw= d200_noAtlantic_spW, zero.policy=T)


##	what is the best neighbourhood matrix?  50km wins AGAIN!!
lChao_lat_long_neighbourhood_AIC <- AICtab(lChao_lat_long_4scale_lm, 
	lChao_lat_long_4scale_SAR_errW_dmean, 
	lChao_lat_long_4scale_SAR_errW_d50, 
	lChao_lat_long_4scale_SAR_errW_d100, 
	lChao_lat_long_4scale_SAR_errW_d200, weights=T, base=T)

#capture.output(lChao_lat_long_neighbourhood_AIC, file='lChao_lat_long_neighbourhood_AIC.csv')

##	LR tests for best spatial weights matrix
##	abs_lat:abs_long STAYS
summary(lChao_lat_long_4scale_SAR_errW_d50_1 <- update(lChao_lat_long_4scale_SAR_errW_d50, ~.-abs_lat:abs_long))
anova(lChao_lat_long_4scale_SAR_errW_d50_1, lChao_lat_long_4scale_SAR_errW_d50)

##	abs_lat:scale GOES (p=0.07)
summary(lChao_lat_long_4scale_SAR_errW_d50_1a <- update(lChao_lat_long_4scale_SAR_errW_d50, ~.-abs_lat:scale))
anova(lChao_lat_long_4scale_SAR_errW_d50_1a, lChao_lat_long_4scale_SAR_errW_d50)

##	abs_long:scale goes (p=0.17) (this differs from spp_rich)
summary(lChao_lat_long_4scale_SAR_errW_d50_1b <- update(lChao_lat_long_4scale_SAR_errW_d50, ~.-abs_long:scale))
anova(lChao_lat_long_4scale_SAR_errW_d50_1b, lChao_lat_long_4scale_SAR_errW_d50)
##	compare these two models with AIC
AIC(lChao_lat_long_4scale_SAR_errW_d50_1a, lChao_lat_long_4scale_SAR_errW_d50_1b)

##	drop both the interactions? NO p=0.002
summary(lChao_lat_long_4scale_SAR_errW_d50_2 <- update(lChao_lat_long_4scale_SAR_errW_d50, ~.-abs_lat:scale - abs_long:scale))
anova(lChao_lat_long_4scale_SAR_errW_d50_2, lChao_lat_long_4scale_SAR_errW_d50)

##	abs_long^2 (and abs_lat:scale) NO p<0.001
summary(lChao_lat_long_4scale_SAR_errW_d50_3a <- update(lChao_lat_long_4scale_SAR_errW_d50_1a, ~.-I(abs_long^2)))
anova(lChao_lat_long_4scale_SAR_errW_d50_3a, lChao_lat_long_4scale_SAR_errW_d50_1a)

##	abs_lat^2 (and abs_lat:scale) 
summary(lChao_lat_long_4scale_SAR_errW_d50_3b <- update(lChao_lat_long_4scale_SAR_errW_d50_1a, ~.-I(abs_lat^2)))
anova(lChao_lat_long_4scale_SAR_errW_d50_3b, lChao_lat_long_4scale_SAR_errW_d50_1a)

##	abs_long^2 (and abs_long:scale) NO p<0.001
summary(lChao_lat_long_4scale_SAR_errW_d50_3c <- update(lChao_lat_long_4scale_SAR_errW_d50_1b, ~.-I(abs_long^2)))
anova(lChao_lat_long_4scale_SAR_errW_d50_3c, lChao_lat_long_4scale_SAR_errW_d50_1b)

##	abs_lat^2 (and abs_long:scale) 
summary(lChao_lat_long_4scale_SAR_errW_d50_3b <- update(lChao_lat_long_4scale_SAR_errW_d50_1b, ~.-I(abs_lat^2)))
anova(lChao_lat_long_4scale_SAR_errW_d50_3b, lChao_lat_long_4scale_SAR_errW_d50_1b)


##	is there autocorrelation in the residuals of the simplified SAR?  NO :-)
#moran.test(residuals(lChao_lat_long_4scale_SAR_errW_d50_1a), d50_noAtlantic_spW, zero.policy=T)

##====================================================================================================================
##	diagnostics

##	latitude and longitude
par(mfrow=c(4,2), mar=c(4,4,2,2))
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lChao_lat_long_4scale_SAR_errW_d50) ~ realm));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lChao_lat_long_4scale_SAR_errW_d50) ~ ecoregion));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lChao_lat_long_4scale_SAR_errW_d50) ~ scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lChao_lat_long_4scale_SAR_errW_d50) ~ abs_lat, col=scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lChao_lat_long_4scale_SAR_errW_d50) ~ abs_long, col=scale));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(residuals(lChao_lat_long_4scale_SAR_errW_d50) ~ lextent, col=scale));abline(h=0, lty=2)
plot(fitted(lChao_lat_long_4scale_SAR_errW_d50), residuals(lChao_lat_long_4scale_SAR_errW_d50));abline(h=0, lty=2)
with(cleanTaxon_NoAtlantic_resamp200, plot(fitted(lChao_lat_long_4scale_SAR_errW_d50), lChao, col=scale));abline(c(0,1), lty=2)
