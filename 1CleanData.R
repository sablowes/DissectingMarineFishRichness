##====================================================================================================================
##	update code: 25092017
##  remove Arctic and Southern Oceans from summaries for Indo-Pacific
##  check and run code for revision @ ProcB
##====================================================================================================================
##  Code to prepare data at multiple scales for examination of the components of species richness across lats and longs
##	calculate 3 components of species richness at multiple scales
##	Smith et al data: 
##	local(site) levels: 1 transect, 2 transects, 4 transects
##	(eco)regional level: 20 transects
##	load library for rarefaction curves
library(tidyverse);
library(iNEXT);						# Chao et al functions for interpolation (rarefaction) and extrapolation (S_asymptote)
library(vegan)						# Simpson diversity functions (also does 1/p^2 == ENSpie)
library(rgeos)
library(sp)
library(geosphere)
rm(list=ls())

##  data are available from reeflifesurvey.com
##	RLS survey data downloaded on 15 Feb 2017
fish2016 <- read_csv('~/Dropbox/1current/dissectingRichness/revision/data/Reef_Life_Survey_Global_reef_fish_dataset.csv')

##  set seed for results in revision
set.seed(123)
##====================================================================================================================
##====================================================================================================================
# IAA hotspot from Hughes et al 2002 Ecol Lett
# calculate absolute distance (ºlongitude) from middle of the CT (120ºE)
fish2016 <- fish2016 %>% mutate(clong_CT = ifelse(SiteLong>0 & SiteLong<180, 120 - abs(SiteLong),
					ifelse(SiteLong<0 & SiteLong>-60, 120 + abs(SiteLong), 240-abs(SiteLong))))

with(fish2016, hist(abs(clong_CT)))
str(fish2016)

##	how many transects in full data set
fish2016 %>% summarise(n_distinct(SurveyID))
##====================================================================================================================
##	combine blocks (i.e., the two sides of the transect tape) for whole dataframe
fish2016_2 <- fish2016 %>%
  # remove 0 counts and invertebrates
	dplyr::filter(Total > 0 & Phylum=='Chordata') %>%
	# retain hierarchical structure for each taxon
	group_by(Realm, Ecoregion, SiteCode, SurveyID, Taxon) %>%
	summarise(
	  SiteLat = unique(SiteLat),
		SiteLong = unique(SiteLong),
		clong_CT = unique(clong_CT),
		Total = sum(Total)) %>%
    ungroup()
			
# quick look at the data (# sites, # ecoregions, # realms, latitudinal range)
fish2016_2 %>% 
	summarise(n_transect = n_distinct(SurveyID),
		n_site = n_distinct(SiteCode),
		n_ecoregion = n_distinct(Ecoregion),
		n_realm = n_distinct(Realm),
		min_lat = min(SiteLat),
		max_lat = max(SiteLat))


##	what are the unique levels of Taxon?
##	there are a few taxa that are not ID'd to species
fish2016_2 %>% distinct(Taxon) %>% data.frame


##	clean up in steps
string1 <- 'spp.'
string2 <- '?'
string3 <- 'sp.'

fish2016_cleanTaxon3 <- fish2016_2 %>% filter(!grepl(string1, Taxon, fixed=T)) %>%
	filter(!grepl(string2, Taxon, fixed=T)) %>%
	filter(!grepl(string3, Taxon, fixed=T)) 

fish2016_cleanTaxon3 %>% distinct(Taxon)
##==================================================
##	data summary of clean data set
fish2016_cleanTaxon3 %>% 
	summarise(n_individuals = sum(Total), 
		n_species = n_distinct(Taxon), 
		n_transect = n_distinct(SurveyID),
		n_site = n_distinct(SiteCode),
		n_ecoregion = n_distinct(Ecoregion),
		n_realm = n_distinct(Realm),
		min_lat = min(SiteLat),
		max_lat = max(SiteLat))

##	summary of clean data in the Indo-Pacific
fish2016_cleanTaxon3 %>% 
	dplyr::filter((Realm!='Tropical Atlantic' & Realm!= 'Temperate Northern Atlantic' & Realm!='Temperate South America' & Realm!='Arctic' & Realm!='Southern Ocean')) %>%
	summarise(n_individuals = sum(Total), 
		n_species = n_distinct(Taxon), 
		n_transect = n_distinct(SurveyID),
		n_site = n_distinct(SiteCode),
		n_ecoregion = n_distinct(Ecoregion),
		n_realm = n_distinct(Realm),
		min_lat = min(SiteLat),
		max_lat = max(SiteLat))
##==================================================

##====================================================================================================================
##	Site scale grains (1 transect, 2 transects, 4 transects)
##	transect level (smallest grain)
transect_summary <- fish2016_cleanTaxon3 %>%	
				  # group by individual transect
				  group_by(SurveyID) %>%
				  # transect level summary
			      summarise(
			      	individuals = sum(Total),
			      	spp_rich = n_distinct(Taxon),
			      	chaoRich_observed = ChaoSpecies(Total)$Observed,		# should be equal to spp_rich
			      	chaoRich_estimated = ChaoSpecies(Total)$Estimator,
			      	ENSpie = diversity(Total, index='invsimpson'),
			      	pie = diversity(Total, index='simpson'),
				  	  lat = unique(SiteLat),
				  	  long = unique(SiteLong),
				  	  clong_CT = unique(clong_CT),
				  	  site = unique(SiteCode),
				  	  ecoregion = unique(Ecoregion),
				  	  realm = unique(Realm),
				  	  n_transect = n_distinct(SurveyID))

##====================================================================================================================
##  select sites with 2 or more transects; and randomly sample 2 transects from these sites
site_filter2 <- fish2016_cleanTaxon3 %>%
		# group_by Site
		group_by(SiteCode) %>%
		dplyr::filter(n_distinct(SurveyID)>=2) %>%
    ungroup()
    

# check how many transects at each site
trans_per_site2 <- site_filter2 %>%
    group_by(SiteCode) %>%
    summarise(n_transect = n_distinct(SurveyID)) %>%
  ungroup()

hist(trans_per_site2$n_transect)

site_filter2$Realm <- factor(site_filter2$Realm)
site_filter2$Ecoregion <- factor(site_filter2$Ecoregion)
site_filter2$SiteCode <- factor(site_filter2$SiteCode)
site_filter2$SurveyID <- factor(site_filter2$SurveyID)
str(site_filter2)

## for each unique site, randomly sample 2 transects (now changed to sample 2 transects multiple (resample) times)
site_data2 <- data.frame()
resample <- 200	# do more resamples for final analysis
for(i in 1:length(unique(trans_per_site2$SiteCode))) {
	print(paste('site', i, 'in', length(unique(trans_per_site2$SiteCode)), 'sites'))
	# how many transects at this site
	n_trans <- trans_per_site2$n_transect[i]
	if(n_trans==2){
		rows <- fish2016_cleanTaxon3[which(fish2016_cleanTaxon3$SiteCode==as.character(trans_per_site2$SiteCode[i])),]
		rows$resamp <- 0
		site_data2 <- rbind.data.frame(site_data2, rows)
	}
	else{
		site <- site_filter2[which(site_filter2$SiteCode==as.character(trans_per_site2$SiteCode[i])),]
		site$SurveyID <- factor(site$SurveyID)
		# want only 2 transects, repeated resample times
		rows2 <- data.frame()
		for(j in 1:resample){
			sample_transects <- sample(unique(site$SurveyID),2)
			rows <- fish2016_cleanTaxon3[which(fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[1]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[2])),]
			samp <- j
			rows$resamp <- samp
			rows2 <- rbind.data.frame(rows2, rows)
		}	
	site_data2 <- rbind.data.frame(site_data2, rows2)	
	}
}						

str(site_data2)
head(site_data2)

##====================================================================================================================
##====================================================================================================================
##  select sites with 4 or more transects; and randomly (re)sample 4 transects from these sites
site_filter <- fish2016_cleanTaxon3 %>%
  # group_by Site
	group_by(SiteCode) %>%
  dplyr::filter(n_distinct(SurveyID)>=4) %>%
  ungroup()


# check how many transects at each site
trans_per_site <- site_filter %>%
  group_by(SiteCode) %>%
	summarise(
    Realm = unique(Realm),
		SiteLat = unique(SiteLat),
		SiteLong = unique(SiteLong),
		n_transect = n_distinct(SurveyID)) %>%
  ungroup()

site_filter$Realm <- factor(site_filter$Realm)
site_filter$Ecoregion <- factor(site_filter$Ecoregion)
site_filter$SiteCode <- factor(site_filter$SiteCode)
site_filter$SurveyID <- factor(site_filter$SurveyID)
str(site_filter)

## for each unique site, randomly sample 4 transects
site_data <- data.frame()
resamp <- 200
for(i in 1:length(unique(trans_per_site$SiteCode))) {
	print(paste('site', i, 'in', length(unique(trans_per_site$SiteCode)), 'sites'))
	# how many transects at this site
	n_trans <- trans_per_site$n_transect[i]
	if(n_trans==4){
		rows <- fish2016_cleanTaxon3[which(fish2016_cleanTaxon3$SiteCode==as.character(trans_per_site$SiteCode[i])),]
		rows$resamp <- 0
		site_data <- rbind.data.frame(site_data, rows)
	}
	else{
		site <- site_filter[which(site_filter$SiteCode==as.character(trans_per_site$SiteCode[i])),]
		site$SurveyID <- factor(site$SurveyID)
		# want only 2 transects, repeated resample times
		rows2 <- data.frame()
		for(j in 1:resample){
			sample_transects <- sample(unique(site$SurveyID),4)
			rows <- fish2016_cleanTaxon3[which(fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[1]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[2]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[3]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[4])),]
			samp <- j
			rows$resamp <- samp
			rows2 <- rbind.data.frame(rows2, rows)
		}
	site_data <- rbind.data.frame(site_data, rows2)	
	}
}	

str(site_data)
head(site_data)
##====================================================================================================================
##	resample 20 transects from each ecoregion many (200) times
trans_per_ecoR <- fish2016_cleanTaxon3 %>%
	# group by ecoregion
	group_by(Ecoregion) %>%
	# summarise
	summarise(
		n_transects = n_distinct(SurveyID))
		
#	filter to ecoregions with at 20 transects
ecoregion_filter <- ungroup(fish2016_cleanTaxon3) %>%
  # group_by Ecoregion
	group_by(Ecoregion) %>%
	# reduce to ecoregions with at least 20 transects
	dplyr::filter(n_distinct(SurveyID)>=20)

# check how many transects at each site
ecoregion_trans_per <- ecoregion_filter %>%
							summarise(n_transect = n_distinct(SurveyID))						

ecoregion_filter$Realm <- factor(ecoregion_filter$Realm)
ecoregion_filter$Ecoregion <- factor(ecoregion_filter$Ecoregion)
ecoregion_filter$SiteCode <- factor(ecoregion_filter$SiteCode)
ecoregion_filter$SurveyID <- factor(ecoregion_filter$SurveyID)
str(ecoregion_filter)

## at each unique ecoregion, randomly sample 20 transects
ecoregion_data <- data.frame()
resample <- 200

for(i in 1:length(unique(ecoregion_trans_per$Ecoregion))) {
	print(paste('ecoregion', i, 'in', length(unique(ecoregion_trans_per$Ecoregion)), 'ecoregions'))
	# how many transects at this site
	n_trans <- ecoregion_trans_per$n_transect[i]
	if(n_trans==20){
		rows <- fish2016_cleanTaxon3[which(fish2016_cleanTaxon3$Ecoregion==as.character(ecoregion_trans_per$Ecoregion[i])),]
		rows$resamp <- 0
		ecoregion_data <- rbind.data.frame(ecoregion_data, rows)	

	}
	else{
		eco <- ecoregion_filter[which(ecoregion_filter$Ecoregion==as.character(ecoregion_trans_per$Ecoregion[i])),]
		eco$SurveyID <- factor(eco$SurveyID)
		# want only 4 transects
		rows2 <- data.frame()
		for(j in 1:resample){
			sample_transects <- sample(unique(eco$SurveyID),20)
			rows <- fish2016_cleanTaxon3[which(fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[1]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[2]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[3]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[4]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[5]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[6]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[7]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[8]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[9]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[10]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[11])  | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[12]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[13]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[14]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[15]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[16]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[17]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[18]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[19]) | fish2016_cleanTaxon3$SurveyID==as.character(sample_transects[20])),]
			samp <- j
			rows$resamp <- samp
			rows2 <- rbind.data.frame(rows2, rows)
		}
	ecoregion_data <- rbind.data.frame(ecoregion_data, rows2)
	}
}						

dim(ecoregion_data)
str(ecoregion_data)
head(ecoregion_data)
##====================================================================================================================

##====================================================================================================================
##	need to collate taxa over each resampling event within each site (for Chao richness and Simpson diversity estimators)
##	2-transect scale
intermediate_summary2a <- ungroup(site_data2) %>%
  # retain hierarchy
	group_by(Realm, Ecoregion, SiteCode, resamp, Taxon) %>%
	# collate taxa within each resample (retain hierarchy)
	# retain geographic coordinates
	  summarise(
		  Total = sum(Total),
			SiteLat = unique(SiteLat),
			SiteLong = unique(SiteLong),
			clong_CT = unique(clong_CT))
		
##	calculate metrics for each resampling event (at each site)
intermediate_summary2 <- ungroup(intermediate_summary2a) %>%
  group_by(SiteCode, resamp) %>%
	# calculate summaries stats
	summarise(
	  individuals = sum(Total),
		spp_rich = n_distinct(Taxon),
		chaoRich_observed = ChaoSpecies(Total)$Observed,				
		chaoRich_estimated = ChaoSpecies(Total)$Estimator,				
		ENSpie = diversity(Total, index='invsimpson'),
    pie = diversity(Total, index='simpson'),
		lat = unique(SiteLat),
		long = unique(SiteLong),
		clong_CT = unique(clong_CT),
		site = unique(SiteCode),
		ecoregion = unique(Ecoregion),
		realm = unique(Realm))
					  	
##	calculate metrics as mean over all resamples
intermediate_summary2_siteMeans <- ungroup(intermediate_summary2) %>%
  group_by(SiteCode) %>%
	# calculate summaries stats
	summarise(
	  individuals = mean(individuals),
		spp_rich = mean(spp_rich),			
    chaoRich_observed = mean(chaoRich_observed),					
		chaoRich_estimated = mean(chaoRich_estimated),				
		ENSpie = sum(ENSpie)/n_distinct(resamp),							
		pie = mean(pie),
		lat = unique(lat),
		long = unique(long),
		clong_CT = unique(clong_CT),
		site = unique(SiteCode),
		ecoregion = unique(ecoregion),
		realm = unique(realm))

##	add the number of transects sampled (aggregated)					  	
intermediate_summary2_siteMeans$n_transect <- 2

##	depending on the resampling of transects	 different sites now have different
##	numbers of unique transects
##	may want to include this to control for different areas sampled??
n_transect <- ungroup(site_data2) %>%
  group_by(SiteCode) %>%
	summarise(n_unique_transect = n_distinct(SurveyID))					  	


##	include the number of unique transects with the site means
trans_2_siteMeans <- inner_join(intermediate_summary2_siteMeans,  n_transect)

##=======================
##  4 transect scale
##	collate taxa at each site for calculation of Chao richness estimator (for each resampling event)
intermediate_summary_a <- ungroup(site_data) %>%
  group_by(Realm, Ecoregion, SiteCode, resamp, Taxon) %>%
	summarise(
	  Total = sum(Total),
		SiteLat = unique(SiteLat),
		SiteLong = unique(SiteLong),
		clong_CT = unique(clong_CT))
		

intermediate_summary <- ungroup(intermediate_summary_a) %>%
  group_by(SiteCode, resamp) %>%
	# calculate summary stats
	summarise(
	  individuals = sum(Total),
		spp_rich = n_distinct(Taxon),
    chaoRich_observed = ChaoSpecies(Total)$Observed,				
		chaoRich_estimated = ChaoSpecies(Total)$Estimator,				
		ENSpie = diversity(Total, index='invsimpson'),
		pie = diversity(Total, index='simpson'),
		lat = unique(SiteLat),
		long = unique(SiteLong),
		clong_CT = unique(clong_CT),
		site = unique(SiteCode),
		ecoregion = unique(Ecoregion),
		realm = unique(Realm))

##	calculate metrics as mean over all resamples within each site
intermediate_summary_siteMeans <- ungroup(intermediate_summary) %>%
  group_by(SiteCode) %>%
	# calculate summaries stats
	summarise(
	  individuals = mean(individuals),
		spp_rich = mean(spp_rich),
		chaoRich_observed = mean(chaoRich_observed),
		chaoRich_estimated = mean(chaoRich_estimated),
		ENSpie = sum(ENSpie)/n_distinct(resamp),
    pie = mean(pie),
		lat = unique(lat),
		long = unique(long),
		clong_CT = unique(clong_CT),
		site = unique(SiteCode),
		ecoregion = unique(ecoregion),
		realm = unique(realm))

intermediate_summary_siteMeans$n_transect <- 4						  	

##	depending on the resampling of transects different sites now have different
##	numbers of unique transects
n_transect <- ungroup(site_data) %>%
  group_by(SiteCode) %>%
	summarise(n_unique_transect = n_distinct(SurveyID))					  	

trans_4_siteMeans <- inner_join(intermediate_summary_siteMeans,  n_transect)

###=======================
#	20 transect scale
##	group taxa within ecoregions for calculation of Chao estimators (for each resample)
ecoregion_summary_a <- dplyr::ungroup(ecoregion_data) %>%
  dplyr::group_by(Realm, Ecoregion, resamp, Taxon) %>%
	# calculate summaries stats
	dplyr::summarise(
	  Total = sum(Total))

##  calculate metrics for each resample in each ecoregion 				      		
ecoregion_summary <- data.frame(ecoregion_summary_a) %>%
	dplyr::group_by(Ecoregion, resamp) %>%
	dplyr::summarise(
	  individuals = sum(Total),
		spp_rich = n_distinct(Taxon),
		chaoRich_observed = ChaoSpecies(Total)$Observed,				
		chaoRich_estimated = ChaoSpecies(Total)$Estimator,				
		ENSpie = diversity(Total, index='invsimpson'),
		pie = diversity(Total, index='simpson'),
		ecoregion = unique(Ecoregion),
		realm = unique(Realm))

##  take the mean across all of the resamples						  	
ecoregion_summary_ecoMeans <- ungroup(ecoregion_summary) %>%
	group_by(Ecoregion) %>%
	summarise(					  	
	  individuals = mean(individuals),
		spp_rich = mean(spp_rich),
    chaoRich_observed = mean(chaoRich_observed),
		chaoRich_estimated = mean(chaoRich_estimated),
		ENSpie = sum(ENSpie)/n_distinct(resamp),
		pie = mean(pie),
		ecoregion = unique(ecoregion),
		realm = unique(realm))
						  	
ecoregion_summary_ecoMeans$n_transect <- 20

##	want to calculate extent and a single coordinate (i.e., the centre) for each ecoregion
##	extent is equal to the geographic area within a convex hull that bounds the sites within each ecoregion
eco_site_coords_nest <- ungroup(ecoregion_data) %>%
	##	select columns
	select(Realm, Ecoregion, SiteCode, SiteLat, SiteLong) %>%
	##	get unique sites 
	distinct(SiteCode, .keep_all=T) %>%
	##	within ecoregions, get min/max coords for polygon edges
	group_by(Realm, Ecoregion) %>%
	nest()

##	i can't get purrr::map and chull to play nice so a loop it is!
extent = numeric()
centre_x <- numeric()
centre_y <- numeric()
vertices_check <- data.frame()
for(i in 1:nrow(eco_site_coords_nest)){
	hull = chull(x=unlist(eco_site_coords_nest$data[[i]][,'SiteLong']), y=unlist(eco_site_coords_nest$data[[i]][,'SiteLat']))	
	vertices = eco_site_coords_nest$data[[i]][hull,c('SiteLong', 'SiteLat')]
	info = cbind.data.frame(Realm=rep(eco_site_coords_nest$Realm[i], times=nrow(vertices)), Ecoregion=rep(eco_site_coords_nest$Ecoregion[i], times=nrow(vertices)), vertices)
	vertices_check = rbind.data.frame(vertices_check, info)
	extent[i] = areaPolygon(data.frame(x=vertices$SiteLong, y=vertices$SiteLat))	
	centre_x[i] = geomean(cbind(x=vertices$SiteLong, y=vertices$SiteLat))[1]
	centre_y[i] = geomean(cbind(x=vertices$SiteLong, y=vertices$SiteLat))[2]
}

##	combine ecoregion and extent
ecoregion_extent <- cbind.data.frame(eco_site_coords_nest[,2], extent,  lat=centre_y, long=centre_x)

##	depending on the resampling of transects	 different sites now have different
##	numbers of unique transects
n_transect <- ungroup(ecoregion_data) %>%
  group_by(Ecoregion) %>%
	summarise(n_unique_transect = n_distinct(SurveyID))					  	

ecoregion_means <- inner_join(ecoregion_summary_ecoMeans,  n_transect)
ecoregion_means_extent <- inner_join(ecoregion_means,  ecoregion_extent)

##	calculate longitude centred on 120ºE, and then fix up order of columns (i.e., match to smaller scales)
ecoregion_means_extent$clong_CT <- with(ecoregion_means_extent, ifelse(long>0 & long<180, 120 - abs(long),
	ifelse(long<0 & long>-60, 120 + abs(long), 240-abs(long))))
names(ecoregion_means_extent)
ecoregion_means_extent <- ecoregion_means_extent[,c(1:7,13:15,8:12)]
##====================================================================================================================
# join scales for saving
##	first add number of unique transects to transect-scale data
transect_summary$n_unique_transect <- 1
##	need to add an extent column (set to grain for transect and site scales)
transect_summary$extent <- transect_summary$n_transect * 500
trans_2_siteMeans$extent <- trans_2_siteMeans$n_transect * 500
trans_4_siteMeans$extent <- trans_4_siteMeans$n_transect * 500


names(transect_summary)[-c(1,11)]
names(trans_2_siteMeans)[-c(1,11)]
names(trans_4_siteMeans)[-c(1,11)]
names(ecoregion_means_extent)[-c(1)]

transect_summary1 <- dplyr::select(transect_summary, -SurveyID, - site)
trans_2_siteMeans1 <- data.frame(trans_2_siteMeans)[-c(1,11)]
trans_4_siteMeans1 <- data.frame(trans_4_siteMeans)[-c(1,11)]
ecoregion_means1 <- data.frame(ecoregion_means_extent)[-c(1)]

##  need to add scale indicator
transect_summary1$scale <- '500m2'
trans_2_siteMeans1$scale <- '1000m2'
trans_4_siteMeans1$scale <- '2000m2'
ecoregion_means1$scale <- '10 000m2'

cleanTaxon_200_revision <- rbind.data.frame(transect_summary1, 
					trans_2_siteMeans1,
					trans_4_siteMeans1,
					ecoregion_means1)

str(cleanTaxon_200_revision)
# add scale indicator covariate
cleanTaxon_200_revision$scale <- factor(cleanTaxon_200_revision$scale, levels=c('500m2', '1000m2', '2000m2', '10 000m2'))
head(cleanTaxon_200_revision)

# setwd('~/Dropbox/1current/dissectingRichness/revision/')
# save(cleanTaxon_200_revision, transect_summary,
# trans_2_siteMeans,
# trans_4_siteMeans,
# ecoregion_means, 
# transect_summary1,
# trans_2_siteMeans1,
# trans_4_siteMeans1,
# ecoregion_means1, 
# site_data2,
# site_data,
# ecoregion_data,
# intermediate_summary2a,
# intermediate_summary2,
# intermediate_summary2_siteMeans,
# intermediate_summary_a,
# intermediate_summary,
# intermediate_summary_siteMeans,
# ecoregion_summary_a,
# ecoregion_summary,
# ecoregion_summary_ecoMeans, file='cleanTaxon_200_revision.Rdata')
##====================================================================================================================