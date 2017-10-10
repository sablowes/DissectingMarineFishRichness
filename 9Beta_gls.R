##============================================================================  
##	200616: modified to fit models to means of 200 resamples...
##=====================================================================================================================  
##=====================================================================================================================  
##	fit models to check these patterns
library(nlme)
library(bbmle)
library(car)
library(tidyverse)
library(spdep)

rm(list=ls())
##	get the data cleaned for beta diversity (scale-dependent ratios of ENSpie)
load('~/Dropbox/1current/dissectingRichness/revision/beta_es_revision.Rdata')
##=====================================================================================================================  
str(all_es)
all_es$abs_lat <- with(all_es, abs(lat))
all_es$abs_long <- with(all_es, abs(clong_CT))

##	reduce to ENSpie data only
beta_pie <- dplyr::filter(all_es, es=='ENSpie')
beta_s <- dplyr::filter(all_es, es=='S')
##	reduce to Indo-Pacific data for fitting models with longitude only (and lat:long interaction)

beta_pie_NoAtlantic <- dplyr::filter(beta_pie, 
  (realm!='Tropical Atlantic' & realm!= 'Temperate Northern Atlantic'& realm!='Temperate South America' & realm!='Arctic' & realm!='Southern Ocean'))
beta_pie_NoAtlantic$realm <- factor(beta_pie_NoAtlantic$realm)
beta_pie_NoAtlantic$ecoregion <- factor(beta_pie_NoAtlantic$ecoregion)
beta_s_NoAtlantic <- dplyr::filter(beta_s, 
  (realm!='Tropical Atlantic' & realm!= 'Temperate Northern Atlantic'& realm!='Temperate South America' & realm!='Arctic' & realm!='Southern Ocean'))
beta_s_NoAtlantic$realm <- factor(beta_s_NoAtlantic$realm)
beta_s_NoAtlantic$ecoregion <- factor(beta_s_NoAtlantic$ecoregion)
str(beta_s_NoAtlantic)

##=========================================================
# create spatial weights matrices for checking spatial autocorrelation of residuals
##	for the Indo-Pacific only
beta_neighb_1_only_noAtlantic <- knn2nb(knearneigh(cbind(beta_pie_NoAtlantic$long, beta_pie_NoAtlantic$lat), k=1, longlat=TRUE))

##	check symmetry
sapply(list(beta_neighb_1_only_noAtlantic), function(x) is.symmetric.nb(x, verbose=F, force=T))

##	make symmetric so as 'sparse' matrix methods are available
beta_neighb_1_only_noAtlantic_sym <- make.sym.nb(beta_neighb_1_only_noAtlantic)
##	check all rows should have at least one neighbour
which(rowSums(nb2mat(beta_neighb_1_only_noAtlantic_sym))!=1)

# create spatial weights object (list): locations are either neighbours (=1) or not (=0); 
# note these values are then row-normalized (i.e., sum to 1)
beta_idwW_noAtlantic <- nb2listw(beta_neighb_1_only_noAtlantic_sym)

##=========================================================
##	lat:long models (linear, no variance covariates)
summary(ENSpie_beta_lat_long_lm <- lm(log_ratio_es ~ abs_lat + abs_long + I(abs_lat^2) + I(abs_long^2) +
  + scale + abs_lat:scale + abs_long:scale + abs_long:abs_lat, data= beta_pie_NoAtlantic))

car::Anova(ENSpie_beta_lat_long_lm)
##	residuals of lm are NOT spatially autocorrelated
#lm.morantest(ENSpie_beta_lat_long_lm, listw=beta_idwW_noAtlantic)

# visual inspection
# ENSpie beta diversity
par(mfrow=c(4,2), mar=c(4,4,2,2))
with(beta_pie_NoAtlantic, plot(residuals(ENSpie_beta_lat_long_lm) ~ realm));abline(h=0, lty=2)
with(beta_pie_NoAtlantic, plot(residuals(ENSpie_beta_lat_long_lm) ~ ecoregion));abline(h=0, lty=2)
with(beta_pie_NoAtlantic, plot(residuals(ENSpie_beta_lat_long_lm) ~ scale));abline(h=0, lty=2)
with(beta_pie_NoAtlantic, plot(residuals(ENSpie_beta_lat_long_lm) ~ abs_lat, col=realm));abline(h=0, lty=2)
with(beta_pie_NoAtlantic, plot(residuals(ENSpie_beta_lat_long_lm) ~ abs_long, col=realm));abline(h=0, lty=2)
plot(fitted(ENSpie_beta_lat_long_lm), residuals(ENSpie_beta_lat_long_lm));abline(h=0, lty=2)
with(beta_pie_NoAtlantic, plot(fitted(ENSpie_beta_lat_long_lm), log_ratio_es, col=scale));abline(c(0,1), lty=2)


# statistical tests for heteroscedascity
# model residuals
car::leveneTest(residuals(ENSpie_beta_lat_long_lm), beta_pie_NoAtlantic$realm)
car::leveneTest(residuals(ENSpie_beta_lat_long_lm), beta_pie_NoAtlantic$ecoregion)
car::leveneTest(residuals(ENSpie_beta_lat_long_lm), beta_pie_NoAtlantic$scale)			# this is the problem

# partial residual plots
library(effects)
plot(Effect(c('abs_long', 'scale'), ENSpie_beta_lat_long_lm, partial.residuals=T))
plot(Effect(c('abs_lat', 'scale'), ENSpie_beta_lat_long_lm, partial.residuals=T))

##=====================================================================================================================  
# check whether scale variance-covariate helps with heteroscedascity? start with full model
summary(ENSpie_beta_lat_long_vc <- nlme::gls(log_ratio_es ~ abs_lat + abs_long + scale +
  abs_lat:scale + abs_long:scale + abs_long:abs_lat, weights=varIdent(form=~1|scale), data= beta_pie_NoAtlantic, method='ML'))

##	type-II (Wald's?) Chi-squared test of significance
car::Anova(ENSpie_beta_lat_long_vc)

##	model residuals are no longer heteroscedastic among scales
car::leveneTest(resid(ENSpie_beta_lat_long_vc, type='normalized'), beta_pie_NoAtlantic$realm)
car::leveneTest(resid(ENSpie_beta_lat_long_vc, type='normalized'), beta_pie_NoAtlantic$ecoregion)
car::leveneTest(resid(ENSpie_beta_lat_long_vc, type='normalized'), beta_pie_NoAtlantic$scale)			# SOLVED

##  visual inspection
par(mfrow=c(4,2), mar=c(4,4,2,2))	##	ecoregion does not look great, but there are only 3 points per ecoregion
with(beta_pie_NoAtlantic, plot(resid(ENSpie_beta_lat_long_vc, type='normalized') ~ realm));abline(h=0, lty=2)
with(beta_pie_NoAtlantic, plot(resid(ENSpie_beta_lat_long_vc, type='normalized') ~ ecoregion));abline(h=0, lty=2)
with(beta_pie_NoAtlantic, plot(resid(ENSpie_beta_lat_long_vc, type='normalized') ~ scale));abline(h=0, lty=2)
with(beta_pie_NoAtlantic, plot(resid(ENSpie_beta_lat_long_vc, type='normalized') ~ abs_lat, col=realm));abline(h=0, lty=2)
with(beta_pie_NoAtlantic, plot(resid(ENSpie_beta_lat_long_vc, type='normalized') ~ abs_long, col=realm));abline(h=0, lty=2)
plot(fitted(ENSpie_beta_lat_long_vc), resid(ENSpie_beta_lat_long_vc, type='normalized'));abline(h=0, lty=2)
with(beta_pie_NoAtlantic, plot(fitted(ENSpie_beta_lat_long_vc), log_ratio_es, col=scale));abline(c(0,1), lty=2)

##	LR tests for model selection for gls model
##	abs_lat:abs_long stays 
summary(ENSpie_beta_lat_long_vc1 <- update(ENSpie_beta_lat_long_vc, ~.-abs_lat:abs_long))
anova(ENSpie_beta_lat_long_vc1, ENSpie_beta_lat_long_vc)

##	abs_lat:scale goes 
summary(ENSpie_beta_lat_long_vc1a <- update(ENSpie_beta_lat_long_vc, ~.-abs_lat:scale))
anova(ENSpie_beta_lat_long_vc1a, ENSpie_beta_lat_long_vc)	

##	abs_long:scale stay
summary(ENSpie_beta_lat_long_vc1b <- update(ENSpie_beta_lat_long_vc, ~.-abs_long:scale))
anova(ENSpie_beta_lat_long_vc1b, ENSpie_beta_lat_long_vc)

##	abs_lat:scale and lat x scale?
summary(ENSpie_beta_lat_long_vc2 <- update(ENSpie_beta_lat_long_vc, ~.-abs_lat:scale - abs_lat:abs_long))
anova(ENSpie_beta_lat_long_vc2, ENSpie_beta_lat_long_vc)

## drop all two way interactions? NO
summary(ENSpie_beta_lat_long_vc3 <- update(ENSpie_beta_lat_long_vc, ~.-abs_long:scale - abs_lat:abs_long - abs_lat:scale))
anova(ENSpie_beta_lat_long_vc3, ENSpie_beta_lat_long_vc)

##	AICc results for comparison
bbmle::AICctab(ENSpie_beta_lat_long_vc1,
  ENSpie_beta_lat_long_vc,
  ENSpie_beta_lat_long_vc1a,
  ENSpie_beta_lat_long_vc1b, weights=T, base=T)

##=========================================================
##	repeat lat:long models (linear, no variance covariates) for beta_S
summary(S_beta_lat_long_lm <- lm(log_ratio_es ~ abs_lat + abs_long + scale + I(abs_lat^2) +I(abs_long^2) +
  abs_lat:scale + abs_long:scale + abs_long:abs_lat, data= beta_s_NoAtlantic))
car::Anova(S_beta_lat_long_lm)

##	residuals of lm are NOT spatially autocorrelated
lm.morantest(S_beta_lat_long_lm, listw=beta_idwW_noAtlantic)

# check model residuals for heteroscedascity
car::leveneTest(residuals(S_beta_lat_long_lm), beta_s_NoAtlantic$realm)
car::leveneTest(residuals(S_beta_lat_long_lm), beta_s_NoAtlantic$ecoregion)
car::leveneTest(residuals(S_beta_lat_long_lm), beta_s_NoAtlantic$scale)			# this is the problem

##=====================================================================================================================  
##	repeat gls analysis for beta_s (model with variance covariate for heteroscedascity for scale variable)
summary(S_beta_lat_long_vc <- gls(log_ratio_es ~ abs_lat + abs_long + scale + I(abs_lat^2) + I(abs_long^2) +
  abs_lat:scale + abs_long:scale + abs_long:abs_lat, weights=varIdent(form=~1|scale), data= beta_s_NoAtlantic, method='ML'))

car::Anova(S_beta_lat_long_vc)

##	model residuals are no longer heteroscedastic among scales
car::leveneTest(resid(S_beta_lat_long_vc, type='normalized'), beta_s_NoAtlantic$realm)
car::leveneTest(resid(S_beta_lat_long_vc, type='normalized'), beta_s_NoAtlantic$ecoregion)
car::leveneTest(resid(S_beta_lat_long_vc, type='normalized'), beta_s_NoAtlantic$scale)			# SOLVED

##	visual diagnostics
par(mfrow=c(4,2), mar=c(4,4,2,2))	##	ecoregion does not look great, but statistical test says ok (there are only three points)
with(beta_s_NoAtlantic, plot(resid(S_beta_lat_long_vc, type='normalized') ~ realm));abline(h=0, lty=2)
with(beta_s_NoAtlantic, plot(resid(S_beta_lat_long_vc, type='normalized') ~ ecoregion));abline(h=0, lty=2)
with(beta_s_NoAtlantic, plot(resid(S_beta_lat_long_vc, type='normalized') ~ scale));abline(h=0, lty=2)
with(beta_s_NoAtlantic, plot(resid(S_beta_lat_long_vc, type='normalized') ~ abs_lat, col=realm));abline(h=0, lty=2)
with(beta_s_NoAtlantic, plot(resid(S_beta_lat_long_vc, type='normalized') ~ abs_long, col=realm));abline(h=0, lty=2)
plot(fitted(S_beta_lat_long_vc), resid(S_beta_lat_long_vc, type='normalized'));abline(h=0, lty=2)
with(beta_s_NoAtlantic, plot(fitted(S_beta_lat_long_vc), log_ratio_es, col=scale));abline(c(0,1), lty=2)


##	LR tests for model selection for gls model
##	abs_lat:abs_long can go p = 0.8
summary(S_beta_lat_long_vc1 <- update(S_beta_lat_long_vc, ~.-abs_lat:abs_long))
anova(S_beta_lat_long_vc1, S_beta_lat_long_vc)

##	abs_lat:scale goes p = 0.2
summary(S_beta_lat_long_vc1a <- update(S_beta_lat_long_vc, ~.-abs_lat:scale))
anova(S_beta_lat_long_vc1a, S_beta_lat_long_vc)	

##	abs_long:scale goes 
summary(S_beta_lat_long_vc1b <- update(S_beta_lat_long_vc, ~.-abs_long:scale))
anova(S_beta_lat_long_vc1b, S_beta_lat_long_vc)

##	drop both lat x scale and lat x long? 
summary(S_beta_lat_long_vc2 <- update(S_beta_lat_long_vc, ~.-abs_lat:scale - abs_lat:abs_long))
anova(S_beta_lat_long_vc2, S_beta_lat_long_vc)

## drop latitude^2? 
summary(S_beta_lat_long_vc4 <- update(S_beta_lat_long_vc2, ~. - I(abs_lat^2)))
anova(S_beta_lat_long_vc4, S_beta_lat_long_vc2)

## drop longitude^2 
summary(S_beta_lat_long_vc4a <- update(S_beta_lat_long_vc2, ~. - I(abs_long^2)))
anova(S_beta_lat_long_vc4a, S_beta_lat_long_vc2)


##  check simpler model with no second order terms
summary(S_beta_lat_long_vc5 <- update(S_beta_lat_long_vc2, ~. - I(abs_lat^2) - I(abs_long^2)))
anova(S_beta_lat_long_vc5, S_beta_lat_long_vc2)
car::Anova(S_beta_lat_long_vc5)  

##=====================================================================================================================  
##	predictions without the Atlantic
str(all_es)
head(all_es)

realm_es_data_NoAtlantic <- beta_pie_NoAtlantic %>%
	dplyr::group_by(scale) %>%
	dplyr::summarise(
		min_lat = min(abs_lat),
		max_lat = max(abs_lat),
		mean_lat = mean(abs_lat),
		median_lat = median(abs_lat),
		min_long = min(abs_long),
		max_long = max(abs_long),
		mean_long = mean(abs_long),
		median_long = median(abs_long))


new_dat <- data.frame()
es_newData_lat_NoAtlantic <- data.frame()
for(i in 1:nrow(realm_es_data_NoAtlantic)){
	#realm <- rep(realm_es_data $realm[i], 100)
	scale <- rep(realm_es_data_NoAtlantic $scale[i], 100)
	abs_lat <- seq(realm_es_data_NoAtlantic $min_lat[i], realm_es_data_NoAtlantic $max_lat[i], length=100)
	abs_long <- rep(realm_es_data_NoAtlantic $mean_long[i], 100)
	new_dat <- cbind.data.frame(scale, abs_lat, abs_long)
	es_newData_lat_NoAtlantic <- rbind.data.frame(es_newData_lat_NoAtlantic, new_dat)
}

new_dat <- data.frame()
es_newData_long_NoAtlantic <- data.frame()
for(i in 1:nrow(realm_es_data_NoAtlantic)){
#	realm <- rep(realm_es_data$realm[i], 100)
	scale <- rep(realm_es_data_NoAtlantic $scale[i], 100)
	abs_long <- seq(realm_es_data_NoAtlantic $min_long[i], realm_es_data_NoAtlantic $max_long[i], length=100)
	abs_lat <- rep(realm_es_data_NoAtlantic $mean_lat[i], 100)
	new_dat <- cbind.data.frame(scale, abs_lat, abs_long)
	es_newData_long_NoAtlantic <- rbind.data.frame(es_newData_long_NoAtlantic, new_dat)
}

es_newData_lat_NoAtlantic$scale <- factor(es_newData_lat_NoAtlantic$scale, levels=c('site1', 'site2', 'ecoregion'))
es_newData_long_NoAtlantic$scale <- factor(es_newData_long_NoAtlantic$scale, levels=c('site1', 'site2', 'ecoregion'))

##	gls predictions
library(AICcmodavg)
# beta_ENSPIE
lat_gls_predictions <- predictSE.gls(ENSpie_beta_lat_long_vc1a, newdata=es_newData_lat_NoAtlantic, se.fit=T)
es_newData_lat_NoAtlantic$predicted_ENSpie_gls <- lat_gls_predictions$fit
es_newData_lat_NoAtlantic$predicted_ENSpie_glsSE <- lat_gls_predictions$se.fit
# beta_S
lat_gls_predictionsS <- predictSE.gls(S_beta_lat_long_vc4, newdata=es_newData_lat_NoAtlantic, se.fit=T)
es_newData_lat_NoAtlantic$predicted_S_gls <- lat_gls_predictionsS$fit
es_newData_lat_NoAtlantic$predicted_S_glsSE <- lat_gls_predictionsS$se.fit


long_gls_predictions <- predictSE.gls(ENSpie_beta_lat_long_vc1a, newdata=es_newData_long_NoAtlantic, se.fit=T)
es_newData_long_NoAtlantic$predicted_ENSpie_gls <- long_gls_predictions$fit
es_newData_long_NoAtlantic$predicted_ENSpie_glsSE <- long_gls_predictions$se.fit

long_gls_predictionsS <- predictSE.gls(S_beta_lat_long_vc4, newdata=es_newData_long_NoAtlantic, se.fit=T)
es_newData_long_NoAtlantic$predicted_S_gls <- long_gls_predictionsS$fit
es_newData_long_NoAtlantic$predicted_S_glsSE <- long_gls_predictionsS$se.fit

# add labels for facetting
beta_pie_NoAtlantic <- beta_pie_NoAtlantic %>%
  mutate(model = 'ENSPIE ratio')
beta_s_NoAtlantic <- beta_s_NoAtlantic %>%
  mutate(model = 'Species richness ratio')


BETA_pie_lat_NoAtlantic_gls <- ggplot(filter(beta_pie_NoAtlantic, scale!='site2')) +
  #facet_grid(~model) + 
  	geom_point(aes(x=abs_lat, y=log_ratio_es, colour=scale, shape=scale), size=2) +
	geom_line(data=dplyr::filter(es_newData_lat_NoAtlantic, scale!='site2'), 
	          aes(x=abs_lat, y= predicted_ENSpie_gls, linetype=scale, group=scale, colour=scale),lwd=1.5) +
	geom_ribbon(data=dplyr::filter(es_newData_lat_NoAtlantic, scale!='site2'), 
	            aes(x=abs_lat, ymax= predicted_ENSpie_gls +2* predicted_ENSpie_glsSE, 
	                ymin= predicted_ENSpie_gls-2* predicted_ENSpie_glsSE, linetype=NA, group=scale, colour=scale, fill=scale), alpha=0.5) +
  ylab('Effect size (log-ratio)') +
  xlab('Absolute latitude (º)') +
  theme_bw() +
	theme(legend.position='none')

BETA_pie_long_NoAtlantic_gls <- ggplot(filter(beta_pie_NoAtlantic, scale!='site2')) +
  #facet_grid(~model) + 
  geom_point(aes(x=abs_long, y=log_ratio_es, shape=scale, colour=scale), size=2) +
  geom_line(data=dplyr::filter(es_newData_long_NoAtlantic, scale!='site2'), 
  	          aes(x=abs_long, y= predicted_ENSpie_gls, colour=scale, linetype=scale, group=scale),lwd=1.5) +
	geom_ribbon(data=dplyr::filter(es_newData_long_NoAtlantic, scale!='site2'), 
	            aes(x=abs_long, ymax= predicted_ENSpie_gls +2* predicted_ENSpie_glsSE, 
	                ymin= predicted_ENSpie_gls-2* predicted_ENSpie_glsSE, linetype=NA, group=scale, colour=scale, fill=scale), alpha=0.5) +
	ylab('') +
  xlab('Absolute longitude (º, centred on 120ºE)') +
  	theme_bw() +
	theme(legend.position='none')

BETA_S_lat_NoAtlantic_gls <- ggplot(filter(beta_s_NoAtlantic, scale!='site2')) +
  #facet_grid(~model) + 
  geom_point(aes(x=abs_lat, y=log_ratio_es, colour=scale, shape=scale), size=2) +
  geom_line(data=dplyr::filter(es_newData_lat_NoAtlantic, scale!='site2'), 
            aes(x=abs_lat, y= predicted_S_gls, linetype=scale, group=scale, colour=scale),lwd=1.5) +
  geom_ribbon(data=dplyr::filter(es_newData_lat_NoAtlantic, scale!='site2'), 
              aes(x=abs_lat, ymax= predicted_S_gls +2* predicted_S_glsSE, 
                  ymin= predicted_S_gls-2* predicted_S_glsSE, linetype=NA, group=scale, colour=scale, fill=scale), alpha=0.5) +
  ylab('Effect size (log-ratio)') +
  xlab('Absolute latitude (º)') +
  theme_bw() +
  theme(legend.position='none')

BETA_S_long_NoAtlantic_gls <- ggplot(filter(beta_s_NoAtlantic, scale!='site2')) +
  #facet_grid(~model) + 
  geom_point(aes(x=abs_long, y=log_ratio_es, shape=scale, colour=scale), size=2) +
  geom_line(data=dplyr::filter(es_newData_long_NoAtlantic, scale!='site2'), 
            aes(x=abs_long, y= predicted_S_gls, colour=scale, linetype=scale, group=scale),lwd=1.5) +
  geom_ribbon(data=dplyr::filter(es_newData_long_NoAtlantic, scale!='site2'), 
              aes(x=abs_long, ymax= predicted_S_gls +2* predicted_S_glsSE, 
                  ymin= predicted_S_gls-2* predicted_S_glsSE, linetype=NA, group=scale, colour=scale, fill=scale), alpha=0.5) +
  ylab('') +
  xlab('Absolute longitude (º, centred on 120ºE)') +
  theme_bw() +
  theme(legend.position='none')


cowplot::plot_grid(plotlist = list(BETA_S_lat_NoAtlantic_gls, BETA_S_long_NoAtlantic_gls,
  BETA_pie_lat_NoAtlantic_gls, BETA_pie_long_NoAtlantic_gls), ncol=2, labels = 'auto')

setwd('~/Dropbox/1current/dissectingRichness/MS/figs/')
#ggsave('Fig4_with_betaS_labelled.pdf', width = 200, height = 200, units = 'mm')

# plot beta_S only
cowplot::plot_grid(plotlist = list(BETA_S_lat_NoAtlantic_gls, BETA_S_long_NoAtlantic_gls), ncol=2, labels = 'auto')

#ggsave('Figx_betaS.pdf', width = 200, height = 100, units = 'mm')
