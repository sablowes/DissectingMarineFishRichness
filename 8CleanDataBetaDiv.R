##==========================================================================
##  @sablowes
##  @beta diversity

##  new code to explore patterns of beta diveristy with 
##  latitude and longitude 

##	beta diversity quantified by the ratio of higher scales over smallest scale
##	i.e., transect scale is always the denomiator

rm(list=ls())
##	beta diversity
library(tidyverse)

load('~/Dropbox/1current/dissectingRichness/revision/cleanTaxon_200_revision.Rdata')
##==========================================================================
##	pie/individual/species richness: calculate means for ecoregions 
##	first three levels are site scale (1, 2, 4 transects), final is ecoregion (20 
##	transects) 
str(transect_summary)
str(intermediate_summary2)
str(intermediate_summary)
str(ecoregion_summary)
##====================================================================================================================		
# first calculate the ecoregion means of all the resamples

transect_effect <- transect_summary %>%
  group_by(ecoregion) %>%
  summarise(
    realm = unique(realm),
    ind_mu = mean(individuals),
    pie_mu = mean(pie),
    ENSpie_mu = mean(ENSpie),
    spp_rich_mu = mean(spp_rich),
    chaoRich_mu = mean(chaoRich_estimated),
    lat = rgeos::gCentroid(sp::SpatialPoints(cbind.data.frame(x=long, y=lat)))$y,
    long = rgeos::gCentroid(sp::SpatialPoints(cbind.data.frame(x=long, y=lat)))$x,
	  clong_CT = ifelse(long>0 & long<180, 
							120 - abs(long),
							ifelse(long<0 & long>-60, 
							120 + abs(long), 240-abs(long))),
    scale = '500m2')


intermediate2_effect <- intermediate_summary2 %>%
  ungroup() %>%
  group_by(ecoregion) %>%
  summarise(
    realm = unique(realm),
    ind_mu = mean(individuals),
    pie_mu = mean(pie),
    ENSpie_mu = mean(ENSpie),
    spp_rich_mu = mean(spp_rich),
    chaoRich_mu = mean(chaoRich_estimated),
    lat = rgeos::gCentroid(sp::SpatialPoints(cbind.data.frame(x=long, y=lat)))$y,
    long = rgeos::gCentroid(sp::SpatialPoints(cbind.data.frame(x=long, y=lat)))$x,
	  clong_CT = ifelse(long>0 & long<180, 
							120 - abs(long),
							ifelse(long<0 & long>-60, 
							120 + abs(long), 240-abs(long))),

    scale = '1000m2')

intermediate4_effect <- intermediate_summary %>%
  ungroup() %>%
  group_by(ecoregion) %>%
  summarise(
    realm = unique(realm),
    ind_mu = mean(individuals),
    pie_mu = mean(pie),
    ENSpie_mu = mean(ENSpie),
    spp_rich_mu = mean(spp_rich),
    chaoRich_mu = mean(chaoRich_estimated),
    lat = rgeos::gCentroid(sp::SpatialPoints(cbind.data.frame(x=long, y=lat)))$y,
    long = rgeos::gCentroid(sp::SpatialPoints(cbind.data.frame(x=long, y=lat)))$x,
	  clong_CT = ifelse(long>0 & long<180, 
							120 - abs(long),
							ifelse(long<0 & long>-60, 
							120 + abs(long), 240-abs(long))),

    scale = '2000m2')

ecoregion_effect <- ecoregion_means1 %>%
  ungroup() %>%
  mutate(
    realm = realm,
    ind_mu = individuals,
    pie_mu = pie,
    ENSpie_mu = ENSpie,
    spp_rich_mu = spp_rich,
    chaoRich_mu = chaoRich_estimated,
    lat = lat,
    long = long,
	  clong_CT = clong_CT,
    scale = '10 000m2') %>% 
  select(-individuals,
         -spp_rich,
         -chaoRich_observed,
         -chaoRich_estimated,
         -ENSpie,
         -pie,
         -extent,
         -n_transect,
         -n_unique_transect)

##========================================================================================================================
##	combine the summaries calculated at the ecoregion scale					  							  	
combine_scales_beta <- bind_rows(transect_effect, intermediate2_effect, intermediate4_effect, ecoregion_effect)
str(combine_scales_beta)

combine_scales_beta$scale <- factor(combine_scales_beta$scale, levels=c('500m2', '1000m2', '2000m2', '10 000m2'))

##====================================================================================================================				  	
##	want to get effect sizes for ecoregions that we have estimates at the largest scale	
transect_effect$ecoregion <- factor(transect_effect$ecoregion)
intermediate2_effect$ecoregion <- factor(intermediate2_effect$ecoregion)
intermediate4_effect$ecoregion <- factor(intermediate4_effect$ecoregion)
ecoregion_effect$ecoregion <- factor(ecoregion_effect$ecoregion)

inner_join(transect_effect, ecoregion_effect, by='ecoregion')

eco_effect_size <- matrix(nrow=length(levels(ecoregion_effect$ecoregion)), ncol=25)

for(i in 1:length(levels(ecoregion_effect$ecoregion))){
  # get site name
  ecoregion <- as.character(ecoregion_effect$ecoregion[i])
  
  # not all ecoregions occur in the lower level summaries!!
  
  p_transect <- ifelse(which(as.character(transect_effect $ecoregion)==ecoregion)==0, NA,
                       transect_effect$pie_mu[which(as.character(transect_effect $ecoregion)==ecoregion)])
  p_intermediate2 <- ifelse(sum(as.character(intermediate2_effect$ecoregion)==ecoregion)==0, NA,
                            intermediate2_effect$pie_mu[which(as.character(intermediate2_effect$ecoregion)==ecoregion)])
  p_intermediate4 <- ifelse(sum(as.character(intermediate4_effect$ecoregion)==ecoregion)==0, NA, 
                            intermediate4_effect$pie_mu[which(as.character(intermediate4_effect$ecoregion)==ecoregion)])
  p_ecoregion <- ecoregion_effect$pie_mu[which(as.character(ecoregion_effect$ecoregion)==ecoregion)]
  
  ENSp_transect <- ifelse(which(as.character(transect_effect $ecoregion)==ecoregion)==0, NA,
                          transect_effect$ENSpie_mu[which(as.character(transect_effect $ecoregion)==ecoregion)])
  ENSp_intermediate2 <- ifelse(sum(as.character(intermediate2_effect$ecoregion)==ecoregion)==0, NA,
                               intermediate2_effect$ENSpie_mu[which(as.character(intermediate2_effect$ecoregion)==ecoregion)])
  ENSp_intermediate4 <- ifelse(sum(as.character(intermediate4_effect$ecoregion)==ecoregion)==0, NA, 
                               intermediate4_effect$ENSpie_mu[which(as.character(intermediate4_effect$ecoregion)==ecoregion)])
  ENSp_ecoregion <- ecoregion_effect$ENSpie_mu[which(as.character(ecoregion_effect$ecoregion)==ecoregion)]
  
  spp_transect <- ifelse(which(as.character(transect_effect $ecoregion)==ecoregion)==0, NA,
                         transect_effect$spp_rich_mu[which(as.character(transect_effect $ecoregion)==ecoregion)])
  spp_intermediate2 <- ifelse(sum(as.character(intermediate2_effect$ecoregion)==ecoregion)==0, NA,
                              intermediate2_effect$spp_rich_mu[which(as.character(intermediate2_effect$ecoregion)==ecoregion)])
  spp_intermediate4 <- ifelse(sum(as.character(intermediate4_effect$ecoregion)==ecoregion)==0, NA, 
                              intermediate4_effect$spp_rich_mu[which(as.character(intermediate4_effect$ecoregion)==ecoregion)])
  spp_ecoregion <- ecoregion_effect$spp_rich_mu[which(as.character(ecoregion_effect$ecoregion)==ecoregion)]
  
  Chao_transect <- ifelse(which(as.character(transect_effect $ecoregion)==ecoregion)==0, NA,
                          transect_effect$chaoRich_mu[which(as.character(transect_effect $ecoregion)==ecoregion)])
  Chao_intermediate2 <- ifelse(sum(as.character(intermediate2_effect$ecoregion)==ecoregion)==0, NA,
                               intermediate2_effect$chaoRich_mu[which(as.character(intermediate2_effect$ecoregion)==ecoregion)])
  Chao_intermediate4 <- ifelse(sum(as.character(intermediate4_effect$ecoregion)==ecoregion)==0, NA, 
                               intermediate4_effect$chaoRich_mu[which(as.character(intermediate4_effect$ecoregion)==ecoregion)])	
  Chao_ecoregion <- ecoregion_effect$chaoRich_mu[which(as.character(ecoregion_effect$ecoregion)==ecoregion)]
  
  ind_transect <- ifelse(which(as.character(transect_effect $ecoregion)==ecoregion)==0, NA,
                         transect_effect$ind_mu[which(as.character(transect_effect $ecoregion)==ecoregion)])
  ind_intermediate2 <- ifelse(sum(as.character(intermediate2_effect$ecoregion)==ecoregion)==0, NA,
                              intermediate2_effect$ind_mu[which(as.character(intermediate2_effect$ecoregion)==ecoregion)])
  ind_intermediate4 <- ifelse(sum(as.character(intermediate4_effect$ecoregion)==ecoregion)==0, NA, 
                              intermediate4_effect$ind_mu[which(as.character(intermediate4_effect$ecoregion)==ecoregion)])	
  ind_ecoregion <- ecoregion_effect$ind_mu[which(as.character(ecoregion_effect$ecoregion)==ecoregion)]
  
  eco_effect_size[i,] <- cbind(realm=as.character(ecoregion_effect$realm[which(ecoregion_effect$ecoregion==ecoregion)]),
                               ecoregion=as.character(ecoregion_effect$ecoregion[which(ecoregion_effect$ecoregion==ecoregion)]),
                               lat=as.character(ecoregion_effect$lat[which(ecoregion_effect$ecoregion==ecoregion)]),
                               long=as.character(ecoregion_effect$long[which(ecoregion_effect$ecoregion==ecoregion)]),
                               clong_CT=as.character(ecoregion_effect$clong_CT[which(ecoregion_effect$ecoregion==ecoregion)]),
                               p_transect = p_transect,
                               p_intermediate2 = p_intermediate2,
                               p_intermediate4 = p_intermediate4,
                               p_ecoregion = p_ecoregion,
                               ENSp_transect = ENSp_transect,
                               ENSp_intermediate2 = ENSp_intermediate2,
                               ENSp_intermediate4 = ENSp_intermediate4,
                               ENSp_ecoregion = ENSp_ecoregion,
                               spp_transect = spp_transect,
                               spp_intermediate2 = spp_intermediate2,
                               spp_intermediate4 = spp_intermediate4,
                               spp_ecoregion = spp_ecoregion,
                               Chao_transect = Chao_transect,
                               Chao_intermediate2 = Chao_intermediate2,
                               Chao_intermediate4 = Chao_intermediate4,
                               Chao_ecoregion = Chao_ecoregion,
                               ind_transect = ind_transect,
                               ind_intermediate2 = ind_intermediate2, 
                               ind_intermediate4 = ind_intermediate4,
                               ind_ecoregion = ind_ecoregion)
}

eco_effect_size <- as.data.frame(eco_effect_size)
colnames(eco_effect_size) <- c('realm', 'ecoregion', 'lat', 'long', 'clong_CT', 'p_transect',
                               'p_intermediate2', 'p_intermediate4', 'p_ecoregion', 'ENSp_transect', 'ENSp_intermediate2', 
                               'ENSp_intermediate4', 'ENSp_ecoregion', 'spp_transect', 'spp_intermediate2', 'spp_intermediate4', 
                               'spp_ecoregion', 'Chao_transect', 'Chao_intermediate2', 'Chao_intermediate4', 'Chao_ecoregion', 
                               'ind_transect', 'ind_intermediate2', 'ind_intermediate4', 'ind_ecoregion')

str(eco_effect_size)
# convert to numbers
eco_effect_size$lat <- as.numeric(as.character(eco_effect_size$lat))
eco_effect_size$long <- as.numeric(as.character(eco_effect_size$long))
eco_effect_size$clong_CT <- as.numeric(as.character(eco_effect_size$clong_CT))
eco_effect_size$p_transect <- as.numeric(as.character(eco_effect_size$p_transect))
eco_effect_size$p_intermediate2 <- as.numeric(as.character(eco_effect_size$p_intermediate2))
eco_effect_size$p_intermediate4 <- as.numeric(as.character(eco_effect_size$p_intermediate4))
eco_effect_size$p_ecoregion <- as.numeric(as.character(eco_effect_size$p_ecoregion))
eco_effect_size$ENSp_transect <- as.numeric(as.character(eco_effect_size$ENSp_transect))
eco_effect_size$ENSp_intermediate2 <- as.numeric(as.character(eco_effect_size$ENSp_intermediate2))
eco_effect_size$ENSp_intermediate4 <- as.numeric(as.character(eco_effect_size$ENSp_intermediate4))
eco_effect_size$ENSp_ecoregion <- as.numeric(as.character(eco_effect_size$ENSp_ecoregion))
eco_effect_size$spp_transect <- as.numeric(as.character(eco_effect_size$spp_transect))
eco_effect_size$spp_intermediate2 <- as.numeric(as.character(eco_effect_size$spp_intermediate2))
eco_effect_size$spp_intermediate4 <- as.numeric(as.character(eco_effect_size$spp_intermediate4))
eco_effect_size$spp_ecoregion <- as.numeric(as.character(eco_effect_size$spp_ecoregion))
eco_effect_size$ind_transect <- as.numeric(as.character(eco_effect_size$ind_transect))
eco_effect_size$ind_intermediate2 <- as.numeric(as.character(eco_effect_size$ind_intermediate2))
eco_effect_size$ind_intermediate4 <- as.numeric(as.character(eco_effect_size$ind_intermediate4))
eco_effect_size$ind_ecoregion <- as.numeric(as.character(eco_effect_size$ind_ecoregion))
eco_effect_size$Chao_transect <- as.numeric(as.character(eco_effect_size$Chao_transect))
eco_effect_size$Chao_intermediate2 <- as.numeric(as.character(eco_effect_size$Chao_intermediate2))
eco_effect_size$Chao_intermediate4 <- as.numeric(as.character(eco_effect_size$Chao_intermediate4))
eco_effect_size$Chao_ecoregion <- as.numeric(as.character(eco_effect_size$Chao_ecoregion))
##	reduce to complete cases
eco_effect_size_complete <- eco_effect_size[complete.cases(eco_effect_size),]
head(eco_effect_size_complete)
str(eco_effect_size_complete)

filter(eco_effect_size_complete, ecoregion=='Eastern Brazil')
##===================================================================================================
##  create new dataframe with effect sizes (ratios) for plotting
##	transect scale estimates remain in denominator for all scales
str(eco_effect_size_complete)

all_es <- cbind.data.frame(eco_effect_size_complete[,c(1:5)], 
                           es=rep(c('individuals', 'S', 'ENSpie', 'ChaoRichness'), each=nrow(eco_effect_size_complete)),
                           scale=rep(c('site1', 'site2', 'ecoregion'), each=nrow(eco_effect_size_complete)*4),
                           value = with(eco_effect_size_complete, c(I(ind_intermediate2/ind_transect), 
                                                                    I(spp_intermediate2/spp_transect),
                                                                    I(ENSp_intermediate2/ENSp_transect),
                                                                    I(Chao_intermediate2/Chao_transect),
                                                                    I(ind_intermediate4/ind_transect),
                                                                    I(spp_intermediate4/spp_transect),
                                                                    I(ENSp_intermediate4/ENSp_transect),
                                                                    I(Chao_intermediate4/Chao_transect),
                                                                    I(ind_ecoregion/ind_transect),
                                                                    I(spp_ecoregion/spp_transect),
                                                                    I(ENSp_ecoregion/ENSp_transect),
                                                                    I(Chao_ecoregion/Chao_transect))))

all_es$scale <- factor(all_es$scale, levels=c('site1', 'site2', 'ecoregion'))
all_es$es <- factor(all_es$es, levels=c('individuals', 'S', 'ChaoRichness', 'ENSpie'))
##	add log-ratio
all_es$log_ratio_es <- with(all_es, log(value))
str(all_es)
setwd('~/Dropbox/1current/dissectingRichness/revision/')
##===================================================================================================================== 
#save(all_es, file='beta_es_revision.Rdata')
##===================================================================================================================== 
