##====================================================================================================================
##	create spatial weight matrices for use in simultaneous autoregressive models (SAR)
##  this code creates spatial weight matrices for data from the Indo-Pacific only
##====================================================================================================================
rm(list=ls())
##	load the data:
load('~/Dropbox/1current/dissectingRichness/revision/cleanTaxon_200_revision.Rdata')
library(tidyverse)
library(spdep)
library(ncf)
##============================================================================
# add log-transformed variables for model fitting
str(cleanTaxon_200_revision)
cleanTaxon_200_revision$lind <- log(cleanTaxon_200_revision$individuals, base=10)
cleanTaxon_200_revision$lENSpie <- log(cleanTaxon_200_revision$ENSpie, base=10)
cleanTaxon_200_revision$lspp_rich <- log(cleanTaxon_200_revision$spp_rich, base=10)
cleanTaxon_200_revision$lChao <- log(cleanTaxon_200_revision$chaoRich_estimated, 10)
cleanTaxon_200_revision$abs_lat <- abs(cleanTaxon_200_revision$lat)
cleanTaxon_200_revision$abs_long <- abs(cleanTaxon_200_revision$clong_CT)
cleanTaxon_200_revision$realm <- factor(cleanTaxon_200_revision$realm)
cleanTaxon_200_revision$ecoregion <- factor(cleanTaxon_200_revision$ecoregion)

##	convert extent to km^2 on log-scale
cleanTaxon_200_revision$lextent <- with(cleanTaxon_200_revision, log(extent*1e-6, base=10))

# remove Atlantic Ocean realms, and the Arctic and Southern Oceans
cleanTaxon_NoAtlantic_resamp200 <- dplyr::filter(data.frame(cleanTaxon_200_revision),
  (realm!='Tropical Atlantic' & realm!= 'Temperate Northern Atlantic' & realm!='Temperate South America' & realm!='Arctic' & realm!='Southern Ocean'))
cleanTaxon_NoAtlantic_resamp200$realm <- factor(cleanTaxon_NoAtlantic_resamp200$realm)
cleanTaxon_NoAtlantic_resamp200$ecoregion <- factor(cleanTaxon_NoAtlantic_resamp200$ecoregion)
str(cleanTaxon_NoAtlantic_resamp200)
##============================================================================
##	create matrix of neighbours, and subsequently spatial weights matrix
##	getting one neighbour each determines the minimum distance required for each sample to have
##	at least one neighbour
neighb_1_only_noAtlantic <- knn2nb(knearneigh(cbind(cleanTaxon_NoAtlantic_resamp200$long, cleanTaxon_NoAtlantic_resamp200$lat), k=1, longlat=TRUE))

##	we can use nbdists function to calculate list of distances to nearest neighbours
distances_noAtlantic <- nbdists(neighb_1_only_noAtlantic, coords=cbind(cleanTaxon_NoAtlantic_resamp200 $long,
                                                                       cleanTaxon_NoAtlantic_resamp200$lat), longlat=T)
# summary(unlist(distances_noAtlantic)) 			# maximum value is the distance required to ensure that all transects have at least 
                      												# one neighbour

##	NB: use the dnearneigh function to find all neighbours within a given interpoint distance (d1 - d2)

##	mean distance between neighbours first
dmean_noAtlantic <- dnearneigh(cbind(cleanTaxon_NoAtlantic_resamp200$long, cleanTaxon_NoAtlantic_resamp200$lat), 
	d1 = 0, d2 = mean(unlist(distances_noAtlantic)), longlat=TRUE) 	

##	50 km	
d50_noAtlantic <- dnearneigh(cbind(cleanTaxon_NoAtlantic_resamp200$long, cleanTaxon_NoAtlantic_resamp200$lat), 
	d1 = 0, d2 = 50, longlat=TRUE)

##	100 km	
d100_noAtlantic <- dnearneigh(cbind(cleanTaxon_NoAtlantic_resamp200$long, cleanTaxon_NoAtlantic_resamp200$lat), 
	d1 = 0, d2 = 100, longlat=TRUE)

##	200 km
d200_noAtlantic <- dnearneigh(cbind(cleanTaxon_NoAtlantic_resamp200$long, cleanTaxon_NoAtlantic_resamp200$lat), 
	d1 = 0, d2 = 200, longlat=TRUE)
	
##	check symmetry
# sapply(list(dmean, d50, d100, d200,
# dmean_noAtlantic, d50_noAtlantic, d100_noAtlantic, d200_noAtlantic), function(x) is.symmetric.nb(x, verbose=F, force=T))

##	create spatial weights lists for use in SARs
d50_noAtlantic_spW <- nb2listw(d50_noAtlantic, zero.policy=T)
d100_noAtlantic_spW <- nb2listw(d100_noAtlantic, zero.policy=T)
d200_noAtlantic_spW <- nb2listw(d200_noAtlantic, zero.policy=T)
dmean_noAtlantic_spW <- nb2listw(dmean_noAtlantic, zero.policy=T)
##============================================================================