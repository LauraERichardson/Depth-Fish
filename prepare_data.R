require(tidyverse)

#### Uncapped, raw SPC cylinder replicate data
#### PHIL NOTE - DO NOT SET WD - SCRIPT SHOULD BE EXECUTED IN RELEVANT DIR OR CONTAIN RELATIVE PATHS TO EXECUTION PATH
fish<-read.csv("data/Fishbiomass_depth_Pacific2010-2014.csv", header = TRUE,strip.white = T) 

# remove site slope NAs 
fish <- filter(fish, !is.na(SITE_SLOPE_400m))

# R format housekeeping
fish$SITE_SLOPE_400m <- as.numeric(fish$SITE_SLOPE_400m)
fish$POP_STATUS <- factor(fish$POP_STATUS, levels = c("U","P"))
fish$DEPTH_BIN_new<-factor(fish$DEPTH_BIN_new, levels = c("Shallow","Mid","Deep"))
fish$SITE <- factor(fish$SITE)
fish$ISLAND <- factor(fish$ISLAND)
fish$ECOREGION <- factor(fish$ECOREGION)
fish$DIVER <- factor(fish$DIVER)
fish$X <- factor(fish$X)
fish$OBS_YEAR <- factor(fish$OBS_YEAR, levels = c("2010", "2011","2012","2013","2014"))

###################################################################

# Data cleaning and organization

#Remove islands: Wake (already removed Gardner, Laysan, Maro, Midway, Necker, Nihoa, South Bank)
fish <- fish[!(fish$ISLAND %in% c('Wake')),] %>% droplevels

#Remove some sites: PAL-00107, GUA-00410, PHR-00341 (already removed)
fish <- fish[!(fish$SITE %in% c('PAL-00107', 'GUA-00410','PHR-00341')),] %>% droplevels

#Also remove inner Maug sites because site slope data unreliable:
fish <- fish[!(fish$SITE %in% c('MAU-00077', 'MAU-00124', 'MAU-00162', 'MAU-00163', 'MAU-00174', 'MAU-00178', 'MAU-00179',
                                'MAU-00208', 'MAU-00218', 'MAU-00279','MAU-00281','MAU-00282','MAU-00284','MAU-00285','MAU-00294',
                                'MAU-00334','MAU-00339')),] %>% droplevels

# create subsets of populated and unpopulated sites 
POP<-subset(fish,POP_STATUS=="P") %>% droplevels
UNPOP<-subset(fish,POP_STATUS=="U") %>% droplevels

#Rescale variables:
fish$DEPTH_c <-scale(fish$DEPTH, center = TRUE, scale = TRUE)
fish$SITE_SLOPE_400m_c<-scale(fish$SITE_SLOPE_400m, center = TRUE, scale = TRUE)

save.image('intermed_data/Depth_study_fish_data.RData')
write.csv(fish,"intermed_data/Depth_study_fish_data.csv")