##====================================================================================================================
##  1) code to plot data and predictions from the simplified models for ms
##    (code assumes the data modelled, and th model objects are in the global environment)
##  2) code to calculate total ecoregion diversity and plot map
##    (code includes a modification to Jarrett Byrnes' meowR mapping code to centre map on 120ºE)
##====================================================================================================================
load('~/Dropbox/1current/dissectingRichness/revision/results/revision_workspace.RData')
library(tidyverse)
library(meowR)
# new data for predictions from model
# with the Atlantic Ocean removed (need data for lat and long predictions at multiple scales)
newData_lat_noAtlantic <- with(cleanTaxon_NoAtlantic_resamp200, cbind.data.frame(
  scale = rep(levels(scale), each=100),
  abs_lat = c(seq(min(abs_lat[which(scale=='500m2')]), max(abs_lat[which(scale=='500m2')]), length=100),
              seq(min(abs_lat[which(scale=='1000m2')]), max(abs_lat[which(scale=='1000m2')]), length=100),
              seq(min(abs_lat[which(scale=='2000m2')]), max(abs_lat[which(scale=='2000m2')]), length=100),
              seq(min(abs_lat[which(scale=='10 000m2')]), max(abs_lat[which(scale=='10 000m2')]), length=100)),
  abs_long = c(rep(mean(abs_long[which(scale=='500m2')]), times=100), 
               rep(mean(abs_long[which(scale=='1000m2')]), times=100), 	
               rep(mean(abs_long[which(scale=='2000m2')]), times=100), 
               rep(mean(abs_long[which(scale=='10 000m2')]), times=100)),
  lextent = c(rep(mean(lextent[which(scale=='500m2')]), times=100), 
              rep(mean(lextent[which(scale=='1000m2')]), times=100), 
              rep(mean(lextent[which(scale=='2000m2')]), times=100), 
              rep(mean(lextent[which(scale=='10 000m2')]), times=100))))

newData_lat_noAtlantic$scale <- factor(newData_lat_noAtlantic$scale, levels=c('500m2', '1000m2', '2000m2', '10 000m2'))

newData_long_noAtlantic <- with(cleanTaxon_NoAtlantic_resamp200, cbind.data.frame(
  scale = rep(levels(scale), each=100),
  abs_long = c(seq(min(abs_long[which(scale=='500m2')]), max(abs_long[which(scale=='500m2')]), length=100),
               seq(min(abs_long[which(scale=='1000m2')]), max(abs_long[which(scale=='1000m2')]), length=100),
               seq(min(abs_long[which(scale=='2000m2')]), max(abs_long[which(scale=='2000m2')]), length=100),
               seq(min(abs_long[which(scale=='10 000m2')]), max(abs_long[which(scale=='10 000m2')]), length=100)),
  abs_lat = c(rep(mean(abs_lat[which(scale=='500m2')]), times=100), 
              rep(mean(abs_lat[which(scale=='1000m2')]), times=100), 
              rep(mean(abs_lat[which(scale=='2000m2')]), times=100), 
              rep(mean(abs_lat[which(scale=='10 000m2')]), times=100))),
  lextent = c(rep(mean(lextent[which(scale=='500m2')]), times=100), 
              rep(mean(lextent[which(scale=='1000m2')]), times=100), 
              rep(mean(lextent[which(scale=='2000m2')]), times=100), 
              rep(mean(lextent[which(scale=='10 000m2')]), times=100)))

newData_long_noAtlantic$scale <- factor(newData_long_noAtlantic$scale, levels=c('500m2', '1000m2', '2000m2', '10 000m2'))

##	predictions for individuals in the Indo-Pacific only
newData_lat_noAtlantic$predicted_individuals_full <- predict(lind_lat_long_4scale_SAR_errW_d50, newdata= newData_lat_noAtlantic, listw= d50_noAtlantic_spW)[,1]
newData_lat_noAtlantic$predicted_individuals <- predict(lind_lat_long_4scale_SAR_errW_d50_3, newdata= newData_lat_noAtlantic, listw= d50_noAtlantic_spW)[,1]

newData_long_noAtlantic$predicted_individuals_full <- predict(lind_lat_long_4scale_SAR_errW_d50, newdata= newData_long_noAtlantic, listw= d50_noAtlantic_spW)[,1]
newData_long_noAtlantic$predicted_individuals <- predict(lind_lat_long_4scale_SAR_errW_d50_3, newdata= newData_long_noAtlantic, listw= d50_noAtlantic_spW)[,1]

# ENSPIE
newData_lat_noAtlantic$predicted_ENSpie_full <- predict(lENSpie_lat_long_4scale_SAR_errW_d50, newdata=newData_lat_noAtlantic)[,1]
newData_lat_noAtlantic$predicted_ENSpie <- predict(lENSpie_lat_long_4scale_SAR_errW_d50_2, newdata=newData_lat_noAtlantic)[,1]

newData_long_noAtlantic$predicted_ENSpie_full <- predict(lENSpie_lat_long_4scale_SAR_errW_d50, newdata=newData_long_noAtlantic)[,1]
newData_long_noAtlantic$predicted_ENSpie <- predict(lENSpie_lat_long_4scale_SAR_errW_d50_2, newdata=newData_long_noAtlantic)[,1]

# species richness
newData_lat_noAtlantic$predicted_spp_rich_full <- predict(lspp_rich_lat_long_4scale_SAR_errW_d50, newdata=newData_lat_noAtlantic, listw= d50_noAtlantic_spW)[,1]
newData_lat_noAtlantic$predicted_spp_rich <- predict(lspp_rich_lat_long_4scale_SAR_errW_d50_1a, newdata=newData_lat_noAtlantic, listw= d50_noAtlantic_spW)[,1]

newData_long_noAtlantic$predicted_spp_rich_full <- predict(lspp_rich_lat_long_4scale_SAR_errW_d50, newdata=newData_long_noAtlantic)[,1]
newData_long_noAtlantic$predicted_spp_rich <- predict(lspp_rich_lat_long_4scale_SAR_errW_d50_1a, newdata=newData_long_noAtlantic)[,1]


##	define alpha covariate for use in plotting
cleanTaxon_NoAtlantic_resamp200$alpha <- with(cleanTaxon_NoAtlantic_resamp200, ifelse(scale=='500m2', 0.75,1))

##	two-scales only for presentation??
ind_2scale_lat_noAtlantic <- ggplot(dplyr::filter(data.frame(cleanTaxon_NoAtlantic_resamp200), scale=='500m2' | scale=='10 000m2')) +
  #	facet_wrap(~Realm) +
  geom_point(aes(x=abs_lat, y=individuals, colour=scale, shape=scale, alpha=alpha), size=2) +
  geom_line(data=dplyr::filter(newData_lat_noAtlantic, scale=='500m2' | scale=='10 000m2'), aes(x=abs_lat, y=10^(predicted_individuals), group=scale, linetype=scale), lwd=1.25) +#, colour=scale
  scale_y_log10(limits=c(1, 116240), breaks=c(1,10,100,1000,10000,100000), labels=scales::trans_format('log10', scales::math_format(10^.x))) +
  ylab('Total abundance') +
  xlab('') +
  scale_shape(solid=T) +
  theme_bw() +
  theme(legend.position='none', panel.grid = element_blank()) 

ind_2scale_long_noAtlantic <- ggplot(dplyr::filter(data.frame(cleanTaxon_NoAtlantic_resamp200), scale=='500m2' | scale=='10 000m2')) +
  geom_point(aes(x=abs_long, y=individuals, colour=scale, shape=scale, alpha=alpha), size=2) +
  geom_line(data=dplyr::filter(newData_long_noAtlantic, scale=='500m2' | scale=='10 000m2'), aes(x=abs_long, y=10^(predicted_individuals), group=scale, linetype=scale), lwd=1.25) +
  scale_y_log10(limits=c(1, 116240), breaks=c(1,10,100,1000,10000,100000), labels=scales::trans_format('log10', scales::math_format(10^.x))) +
  ylab('') +#individuals (log-scale)
  xlab('') +#
  scale_shape(solid=T) +
  theme_bw() +
  theme(legend.position='none', panel.grid = element_blank()) 

ENSpie_2scale_lat_noAtlantic <- ggplot(dplyr::filter(cleanTaxon_NoAtlantic_resamp200, scale!='1000m2' & scale!='2000m2')) +
  #	facet_wrap(~Realm) +
  geom_point(aes(x=abs_lat, y=ENSpie, colour=scale, shape=scale, alpha=alpha), size=2) +
  geom_line(data=dplyr::filter(newData_lat_noAtlantic, scale!='1000m2' & scale!='2000m2'), 
            aes(x=abs_lat, y=10^(predicted_ENSpie), group=scale, linetype=scale), lwd=1.25) +#, colour=scale
  scale_y_log10(limits=c(1,38), breaks=c(1,2,4,8,16,32)) +
  ylab('PIE') +
  xlab('Absolute latitude (º)') +
  scale_shape(solid=T) +
  theme_bw() +
  theme(legend.position='none', panel.grid = element_blank())

ENSpie_2scale_long_noAtlantic <- ggplot(dplyr::filter(cleanTaxon_NoAtlantic_resamp200, scale!='1000m2' & scale!='2000m2')) +
  geom_point(aes(x=abs_long, y=ENSpie, colour=scale, shape=scale, alpha=alpha), size=2) +
  geom_line(data=dplyr::filter(newData_long_noAtlantic, scale!='1000m2' & scale!='2000m2'), 
            aes(x=abs_long, y=10^(predicted_ENSpie), group=scale, linetype=scale), lwd=1.25) +
  scale_y_log10(limits=c(1,38), breaks=c(1,2,4,8,16,32)) +
  ylab('') +#ENSpie (log-scale)
  xlab('Absolute longitude (º, centred on 120ºE)') +#
  scale_shape(solid=T) +
  theme_bw() +
  theme(legend.position='none', panel.grid = element_blank())

spp_rich_2scale_lat_noAtlantic <- ggplot(dplyr::filter(cleanTaxon_NoAtlantic_resamp200, scale!='1000m2' & scale!='2000m2')) +
  #	facet_wrap(~Realm) +
  geom_point(aes(x=abs_lat, y=spp_rich, colour=scale, shape=scale, alpha=alpha), size=2) +
  geom_line(data=dplyr::filter(newData_lat_noAtlantic, scale!='1000m2' & scale!='2000m2'), aes(x=abs_lat, y=10^(predicted_spp_rich), group=scale, linetype=scale), lwd=1.25) +#, colour=scale
  scale_y_log10(limits=c(1,400), breaks=c(1,10,100,200,400)) +
  ylab('Species richness') +
  xlab('') +
  scale_shape(solid=T) +
  theme_bw() +
  theme(legend.position='none', panel.grid = element_blank())


spp_rich_2scale_long_noAtlantic <- ggplot(dplyr::filter(cleanTaxon_NoAtlantic_resamp200, scale!='1000m2' & scale!='2000m2')) +
  geom_point(aes(x=abs_long, y=spp_rich, colour=scale, shape=scale, alpha=alpha), size=2) +
  geom_line(data=dplyr::filter(newData_long_noAtlantic, scale!='1000m2' & scale!='2000m2'), aes(x=abs_long, y=10^(predicted_spp_rich), group=scale, linetype=scale), lwd=1.25) +
  scale_y_log10(limits=c(1,400), breaks=c(1,10,100,200,400)) +
  ylab('') +#spp_rich (log-scale)
  xlab('') +#
  scale_shape(solid=T) +
  theme_bw() +
  theme(legend.position='none', panel.grid = element_blank())

cowplot::plot_grid(plotlist = list(spp_rich_2scale_lat_noAtlantic, spp_rich_2scale_long_noAtlantic, 
                   ind_2scale_lat_noAtlantic, ind_2scale_long_noAtlantic, 
                   ENSpie_2scale_lat_noAtlantic, ENSpie_2scale_long_noAtlantic), ncol=2, labels = 'auto')
setwd('~/Dropbox/1current/dissectingRichness/MS/figs/')
#ggsave('Fig3.pdf', width = 200, height = 290, units = 'mm')

##	S_asymptote
newData_lat_noAtlantic$predicted_Chao <- predict(lChao_lat_long_4scale_SAR_errW_d50_d50, newdata=newData_lat_noAtlantic, listw= d50_noAtlantic_spW)[,1]
newData_long_noAtlantic$predicted_Chao <- predict(lChao_lat_long_4scale_SAR_errW_d50_d50, newdata=newData_long_noAtlantic)[,1]

chao_2scale_lat_noAtlantic <- ggplot(dplyr::filter(cleanTaxon_NoAtlantic_resamp200, scale!='1000m2' & scale!='2000m2')) +
  #	facet_wrap(~Realm) +
  geom_point(aes(x=abs_lat, y=spp_rich, colour=scale, shape=scale, alpha=alpha), size=1.5) +
  geom_line(data=filter(newData_lat_noAtlantic, scale!='1000m2' & scale!='2000m2'), aes(x=abs_lat, y=10^(predicted_Chao), group=scale, linetype=scale), lwd=1.25) +#, colour=scale
  scale_y_log10(limits=c(1,440), breaks=c(1,10,100,200,400)) +
  ylab('Estimated richness\n(log-scale)') +
  xlab('') +
  scale_shape(solid=T) +
  theme_bw() +
  theme(legend.position='none') 


chao_2scale_long_noAtlantic <- ggplot(dplyr::filter(cleanTaxon_NoAtlantic_resamp200, scale!='1000m2' & scale!='2000m2')) +
  geom_point(aes(x=abs_long, y=spp_rich, colour=scale, shape=scale, alpha=alpha), size=1.5) +
  geom_line(data=filter(newData_long_noAtlantic, scale!='1000m2' & scale!='2000m2'), aes(x=abs_long, y=10^(predicted_Chao), group=scale, linetype=scale), lwd=1.25) +
  scale_y_log10(limits=c(1,440), breaks=c(1,10,100,200,400)) +
  ylab('') +#spp_rich (log-scale)
  xlab('') +#
  scale_shape(solid=T) +
  theme_bw() +
  theme(legend.position='none') 

##====================================================================================================================
##	all-scales for supplementary material
##	add alpha_4scale for plotting
cleanTaxon_NoAtlantic_resamp200 <- rbind.data.frame(cleanTaxon_NoAtlantic_resamp200 %>%	filter(scale=='500m2') %>% mutate(alpha_4scale=0.6),
                                                    cleanTaxon_NoAtlantic_resamp200 %>%	filter(scale=='1000m2') %>% mutate(alpha_4scale=0.7),
                                                    cleanTaxon_NoAtlantic_resamp200 %>%	filter(scale=='2000m2') %>% mutate(alpha_4scale=0.8),
                                                    cleanTaxon_NoAtlantic_resamp200 %>%	filter(scale=='10 000m2') %>% mutate(alpha_4scale=1))

ind_4scale_lat_noAtlantic <- ggplot(dplyr::filter(data.frame(cleanTaxon_NoAtlantic_resamp200))) +
  #	facet_wrap(~Realm) +
  geom_point(aes(x=abs_lat, y=individuals, colour=scale, shape=scale, alpha=alpha_4scale), size=2) +
  geom_line(data=dplyr::filter(newData_lat_noAtlantic), aes(x=abs_lat, y=10^(predicted_individuals), group=scale, linetype=scale), lwd=1.15) +#, colour=scale
  scale_y_log10(limits=c(1, 116240), breaks=c(1,10,100,1000,10000,100000), labels=scales::trans_format('log10', scales::math_format(10^.x))) +
  ylab('Total abundance') +
  xlab('') +
  scale_shape(solid=F) +
  scale_linetype_manual(values=c('solid', 'dashed', 'dotted', 'dotdash')) +
  theme_bw() +
  theme(legend.position='none', panel.grid = element_blank()) 

ind_4scale_long_noAtlantic <- ggplot(dplyr::filter(data.frame(cleanTaxon_NoAtlantic_resamp200))) +
  geom_point(aes(x=abs_long, y=individuals, colour=scale, shape=scale, alpha=alpha_4scale), size=2) +
  geom_line(data=dplyr::filter(newData_long_noAtlantic), aes(x=abs_long, y=10^(predicted_individuals), group=scale, linetype=scale), lwd=1.15) +
  scale_y_log10(limits=c(1, 116240), breaks=c(1,10,100,1000,10000,100000), labels=scales::trans_format('log10', scales::math_format(10^.x))) +
  ylab('') +#individuals (log-scale)
  xlab('') +#
  scale_shape(solid=F) +
  scale_linetype_manual(values=c('solid', 'dashed', 'dotted', 'dotdash')) +
  theme_bw() +
  theme(legend.position='none', panel.grid = element_blank()) 

ENSpie_4scale_lat_noAtlantic <- ggplot(dplyr::filter(cleanTaxon_NoAtlantic_resamp200)) +
  #	facet_wrap(~Realm) +
  geom_point(aes(x=abs_lat, y=ENSpie, colour=scale, shape=scale, alpha=alpha_4scale), size=2) +
  geom_line(data=dplyr::filter(newData_lat_noAtlantic), aes(x=abs_lat, y=10^(predicted_ENSpie), group=scale, linetype=scale), lwd=1.15) +#, colour=scale
  scale_y_log10(limits=c(1,38), breaks=c(1,2,4,8,16,32)) +
  ylab('PIE') +
  xlab('|Latitude| (º)') +
  scale_shape(solid=F) +
  scale_linetype_manual(values=c('solid', 'dashed', 'dotted', 'dotdash')) +	
  theme_bw() +
  theme(legend.position='none', panel.grid = element_blank()) 

ENSpie_4scale_long_noAtlantic <- ggplot(dplyr::filter(cleanTaxon_NoAtlantic_resamp200)) +
  geom_point(aes(x=abs_long, y=ENSpie, colour=scale, shape=scale, alpha=alpha_4scale), size=2) +
  geom_line(data=dplyr::filter(newData_long_noAtlantic), aes(x=abs_long, y=10^(predicted_ENSpie), group=scale, linetype=scale), lwd=1.15) +
  scale_y_log10(limits=c(1,38), breaks=c(1,2,4,8,16,32)) +
  ylab('') +#ENSpie (log-scale)
  xlab('|Longitude| (º, centred on 120ºE)') +#
  scale_shape(solid=F) +
  scale_linetype_manual(values=c('solid', 'dashed', 'dotted', 'dotdash')) +	
  theme_bw() +
  theme(legend.position='none', panel.grid = element_blank()) 

spp_rich_4scale_lat_noAtlantic <- ggplot(dplyr::filter(data.frame(cleanTaxon_NoAtlantic_resamp200))) +
  #	facet_wrap(~Realm) +
  geom_point(aes(x=abs_lat, y=spp_rich, colour=scale, shape=scale, alpha=alpha_4scale), size=2) +
  geom_line(data=dplyr::filter(newData_lat_noAtlantic), aes(x=abs_lat, y=10^(predicted_spp_rich), group=scale, linetype=scale), lwd=1.15) +#, colour=scale
  scale_y_log10(limits=c(1,400), breaks=c(1,10,100,200,400)) +
  ylab('Species richness') +
  xlab('') +
  scale_shape(solid=F) +
  scale_linetype_manual(values=c('solid', 'dashed', 'dotted', 'dotdash')) +
  theme_bw() +
  theme(legend.position='none', panel.grid = element_blank())

spp_rich_4scale_long_noAtlantic <- ggplot(dplyr::filter(cleanTaxon_NoAtlantic_resamp200)) +
  geom_point(aes(x=abs_long, y=spp_rich, colour=scale, shape=scale, alpha=alpha_4scale), size=2) +
  geom_line(data=dplyr::filter(newData_long_noAtlantic), aes(x=abs_long, y=10^(predicted_spp_rich), group=scale, linetype=scale), lwd=1.15) +
  scale_y_log10(limits=c(1,400), breaks=c(1,10,100,200,400)) +
  ylab('') +#spp_rich (log-scale)
  xlab('') +#
  #	scale_shape(solid=F) +
  guides(alpha=F) +
  scale_linetype_manual(values=c('solid', 'dashed', 'dotted', 'dotdash')) +
  theme_bw() +
  theme(legend.position='none', panel.grid = element_blank()) 
#	theme(legend.position='top', legend.direction='horizontal') 

cowplot::plot_grid(plotlist = list(spp_rich_4scale_lat_noAtlantic, spp_rich_4scale_long_noAtlantic, 
                        ind_4scale_lat_noAtlantic, ind_4scale_long_noAtlantic, 
                        ENSpie_4scale_lat_noAtlantic, ENSpie_4scale_long_noAtlantic), ncol=2, labels = 'auto')
#ggsave('FigS1.pdf', width = 200, height = 290, units = 'mm')
##====================================================================================================================
##	HACK: changes to Jarrett's functions so as map is centred on the Coral Triangle
my_makeMEOWmapData <- function (newdata, fillColName, fillAlphaColName = NULL, regionColName = type, 
                                type = "ECOREGION", excludeNoDataAreas = T) 
{
  regionData <- regions.df
  if (type == "PROVINCE") 
    regionData <- provinces.df
  if (type == "REALM") 
    regionData <- realms.df
  if (regionColName != type) 
    regionData[[regionColName]] <- regionData[[type]]
  regionData <- inner_join(regionData, newdata)
  regionData$score <- regionData[[fillColName]]
  regionData$fillAlpha <- NA
  if (!is.null(fillAlphaColName)) 
    regionData$fillAlpha <- regionData[[fillAlphaColName]]
  if (excludeNoDataAreas && sum(is.na(regionData$score)) > 
      0) 
    regionData <- regionData[-which(is.na(regionData$score)), 
                             ]
  regionData2 <- regionData
  regionData2$long <- regionData2$long+360
  regionData <- rbind(regionData, regionData2)
  
  return(regionData)
}



my_MEOW_map <- function(newdata, fillColName, regionColName = type, type = "ECOREGION", 
                        fillPal = brewer.pal(11, "Spectral"), pal = "Spectral", pathCol = "black", 
                        pathAlpha = 1, guide = guide_colourbar(title = fillColName), 
                        dataOut = FALSE, na.value = NA, add.worldmap = FALSE, fillAlphaColName = NULL, 
                        excludeNoDataAreas = T, prevggplot = NULL, ...) 
{
  regionData <- my_makeMEOWmapData(newdata, fillColName, fillAlphaColName, 
                                   regionColName, type, excludeNoDataAreas)
  if (dataOut) 
    return(regionData)
  if (is.null(prevggplot)) 
    prevggplot <- ggplot()
  ret <- prevggplot + theme_bw(base_size = 17) + geom_polygon(data = regionData, 
                                                              mapping = aes(x = long, y = lat, fill = score, group = group, 
                                                                            alpha = fillAlpha)) + geom_path(data = regionData, 
                                                                                                            color = pathCol, alpha = pathAlpha, mapping = aes(x = long, 
                                                                                                                                                              y = lat, group = group)) + coord_equal(xlim=c(-30, 330), expand=F)
  if (is.numeric(regionData$score)) {
    ret <- ret + scale_fill_gradientn(colours = fillPal, 
                                      guide = guide, na.value = na.value, ...)
  }
  else {
    ret <- ret + scale_fill_manual(values = fillPal, guide = guide, 
                                   na.value = na.value, ...)
  }
  if (add.worldmap) {
    ret <- ret + 
      geom_path(data = my_worldmap.df, aes(x = long, y = lat, group = group)) +
      geom_polygon(data = my_worldmap.df, aes(x = long, y = lat, group = group), fill='white') + 
      scale_x_continuous(name='Longitude (º)', breaks=c(0, 60, 120, 180, 240, 300), labels = c('0', '60E', '120E', '180E', '120W', '60W'), limits=c(-30, 330)) +
      scale_y_continuous(name='Latitude (º)', breaks=c(-50, 0, 50), labels=c('50S', '0', '50N'))
  }
  ret
}
##	need to trick the worldmap to be centred on the CT
my_worldmap.df1 <- fortify(maps::map(fill=TRUE, plot=FALSE))
my_worldmap.df2 <- my_worldmap.df1
my_worldmap.df2$long <- my_worldmap.df2$long+360
my_worldmap.df2$group <- my_worldmap.df2$group + max(my_worldmap.df2$group) + 1
my_worldmap.df <- rbind(my_worldmap.df1, my_worldmap.df2)
str(my_worldmap.df)
##====================================================================================================================
##  load raw data and calculate total species richness for all ecoregions
##	RLS survey data (downloaded from reeflifesurvey.com on 15/02/2016)
fish2016 <- read_csv('~/Dropbox/1current/dissectingRichness/revision/data/Reef_Life_Survey_Global_reef_fish_dataset.csv')

##====================================================================================================================
##	combine blocks (i.e., the two sides of the transect tape) for whole dataframe
fish2016_2 <- fish2016 %>%
  # remove 0 counts and invertebrates
  dplyr::filter(Total > 0 & Phylum=='Chordata') %>%
  # retain hierarchical structure for each taxon
  dplyr::group_by(Realm, Ecoregion, SiteCode, SurveyID, Taxon) %>%
  dplyr::summarise(
    SiteLat = unique(SiteLat),
    SiteLong = unique(SiteLong),
    Total = sum(Total))

##	There are records of with genus only (.spp) and others that are apparently uncertain (contain ? or sp.)
##	clean up in steps
string1 <- 'spp.'
string2 <- '?'
string3 <- 'sp.'

clean_fish2016 <- fish2016_2 %>% filter(!grepl(string1, Taxon, fixed=T)) %>%
  filter(!grepl(string2, Taxon, fixed=T)) %>%
  filter(!grepl(string3, Taxon, fixed=T))

##  calculate species richness in each ecoregion
whole_ecoregion_richness <- ungroup(clean_fish2016) %>%	
  # group by ecoregion
  dplyr::group_by(Ecoregion) %>%
  # ecoregion level summary
  dplyr::summarise(
    spp_rich = n_distinct(Taxon))

##  add log-transformed species richness to data frame for plotting
whole_ecoregion_richness$log10_spp_rich <- log(whole_ecoregion_richness$spp_rich, 10)
#World <- map_data(map = 'world')		

spp_rich_map_wholeEco <- my_MEOW_map(whole_ecoregion_richness, fillColName='log10_spp_rich', 
                                     regionColName='Ecoregion', add.worldmap=T, pathAlpha = 0, 
                                     fillPal=rev(brewer.pal(11,'Spectral')),pal='Spectral', 
                                     guide=guide_colourbar(title='Species\nrichness'), 
                                     breaks=c(log(1,10), log(10,10), log(100,10), log(300,10), log(600,10)), labels=c(1,10, 100,300,600)) 

##	plot ecoregion scale results for species richness
rich_lat_ecoregionOnly <- ggplot(dplyr::filter(cleanTaxon_NoAtlantic_resamp200, scale=='10 000m2')) +
  geom_point(aes(x=abs_lat, y=spp_rich, colour=lspp_rich), size=2) +
  geom_line(data=dplyr::filter(newData_lat_noAtlantic, scale=='10 000m2'), 
            aes(x=abs_lat, y=10^(predicted_spp_rich)), lwd=1) +
  scale_y_continuous(trans='log10', limits=c(20,406), breaks=c(20,100,200,400)) +
  scale_colour_distiller(palette = 'Spectral', direction=-1, 
    breaks=c(log(1,10), log(10,10), log(100,10), log(300,10), log(600,10)), labels=c(1,10, 100,300,600)) +	
  scale_shape(solid=T) +
  ylab('Species richness') +
  xlab('Absolute latitude (º)') +
  theme_bw() +
  theme(legend.position='none', strip.background=element_rect(fill=NA), strip.text=element_text(size=0),
        axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title = element_text(size=16)) 


rich_long_ecoregionOnly <- ggplot(dplyr::filter(cleanTaxon_NoAtlantic_resamp200, scale=='10 000m2')) +
  geom_point(aes(x=abs_long, y=spp_rich, colour=lspp_rich), alpha=1, size=2) +
  geom_line(data=dplyr::filter(newData_long_noAtlantic, scale=='10 000m2'), 
            aes(x=abs_long, y=10^(predicted_spp_rich)), lwd=1) +
  scale_y_continuous(trans='log10', limits=c(20,406), breaks=c(20,100,200,400)) +
  scale_colour_distiller(palette = 'Spectral', direction=-1, 
    breaks=c(log(1,10), log(10,10), log(100,10), log(300,10), log(600,10)), labels=c(1,10, 100,300,600)) +		
  scale_shape(solid=T) +
  ylab('Species richness') +
  xlab('Absolute longitude (º, centred 120ºE)')	+
  theme_bw() +
  theme(legend.position='none', strip.background=element_rect(fill=NA), strip.text=element_text(size=0),
        axis.text.x=element_text(size=14), axis.text.y=element_text(size=14), axis.title = element_text(size=16)) 

#setwd('~/Dropbox/1current/dissectingRichness/MS/figs/')
library(gridExtra)
library(grid)
##	new device
grid.newpage()
##	set up layout
pushViewport(viewport(layout=grid.layout(2,2)))
##  helper function to define region on layout
define_region <- function(row, col){
  viewport(layout.pos.row=row, layout.pos.col=col)
}
## arrange the plots
print(spp_rich_map_wholeEco, vp=define_region(1, 1:2))
print(rich_lat_ecoregionOnly, vp=define_region(2, 1))
print(rich_long_ecoregionOnly, vp=define_region(2, 2))
ggsave('Fig1.pdf', width = 200, height = 200, units = 'mm')
#grid.arrange(grobs=list(rich_lat_ecoregionOnly, rich_long_ecoregionOnly), ncol=2)
