rm(list = ls())
library(readr)
library(dplyr)
library(metafor)  
library(ggplot2)
library(stringr)
library(tidyverse) 
library(ncmeta)
library(viridis)
library(ggthemes)
library(LSD)
library(yardstick)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(gplots)
library(tidyselect)
library(extrafont)
devtools::load_all("/Users/yunpeng/yunkepeng/rbeni/")
library(raster)
library(maps)
library(rworldmap)
library(cowplot)
library(ncdf4)
library(scales)
library(ggpubr)
library(MAd)

#response_ratio_v2: calculate response ratio of each individuals
response_ratio_v2 <- function(df){
  df$logr <- log(df$elevated/df$ambient)
  df_c_sub <- df
  return(df_c_sub)
  #-----------------------------------------------------------------------
  # Output: df_c_sub
  #-----------------------------------------------------------------------
}

#agg_meta_sen_coef: calculate sensitivity coefficient from individuals dataset to site-mean dataset
agg_meta_sen_coef <- function(df,name){
  df$sen_coef <- df$logr/log(df$co2_e/df$co2_a)
  
  df1 <- aggregate(df,by=list(df$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","sen_coef")]
  
  names (df1) <- c("exp",name)
  return(df1)
}

#####3. CO2 effect data
#read Kevin's most recent data
###old vcmax dataset was ("~/data/gcme/kevin/orig_vcmax/JunkePeng_11252021.csv"). Don't delete this old vcmax dataset, since we used this to construct forcing info (this csv has better sampling measurement year data than MESI.csv as it has been well de-bugged for measurement year). However, MESI's vcmax has provided one more plot (biforface_c), all others data don't change. So we use MESI version here.
kevin_othervars <- read.csv("~/data/gcme/kevin_20220222/MESI_2022.csv")
kevin_othervars <- rename(kevin_othervars, c(ambient = x_c, elevated=x_t, ambient_Sd=sd_c, elevated_Sd=sd_t,ambient_Se=se_c,elevated_Se=se_t,n_plots=rep_c,
                                             z=elevation, co2_a=c_c, co2_e=c_t, nfertQ_a = n_c, nfertQ_e = n_t, pfertQ_a = p_c, pfertQ_e = p_t,kfertQ_a = k_c, kfertQ_e = k_t,
                                             warmQ_e1 = w_t1, warmQ_e2 = w_t2, warmQ_e3 = w_t3, Unit=x_units))
kevin_othervars$ambient <-as.numeric(kevin_othervars$ambient)
kevin_othervars$elevated <-as.numeric(kevin_othervars$elevated)
kevin_othervars$ambient_Sd  <-as.numeric(kevin_othervars$ambient_Sd)
kevin_othervars$elevated_Sd  <-as.numeric(kevin_othervars$elevated_Sd)
kevin_othervars$ambient_Se <- as.numeric(kevin_othervars$ambient_Se)
kevin_othervars$elevated_Se <- as.numeric(kevin_othervars$elevated_Se)
kevin_othervars$n_plots  <-as.numeric(kevin_othervars$n_plots)
kevin_othervars$z <- as.numeric(kevin_othervars$z)
kevin_othervars$exp_nam <- kevin_othervars$site

kevin_vcmax <- subset(kevin_othervars,response =="vcmax"|response =="jmax")
kevin2_final <- response_ratio_v2(kevin_vcmax)

#show pft
kevin2_final$ecosystem[kevin2_final$ecosystem=="temperate_forest"] <- "forest"
kevin2_final$ecosystem[kevin2_final$ecosystem=="heathland"] <- "grassland"
kevin2_final$ecosystem[kevin2_final$ecosystem=="shrubland"] <- "forest"
kevin_ecosystem <- as.data.frame(kevin2_final %>% group_by(exp,ecosystem) %>% summarise(number=n()))
#correct exp_nam to make it consistent with GCME
kevin2_final$site[is.na(kevin2_final$site)==TRUE & kevin2_final$exp=="riceface_japan_a_2003_3938_14057_c"] <- "riceface_japan_a_2003_3938_14057"
kevin2_final$site[is.na(kevin2_final$site)==TRUE & kevin2_final$exp=="riceface_japan_a_2004_3938_14057_c"] <- "riceface_japan_a_2004_3938_14057"

#vcmax and jmax data
kevin2_c_vcmax <- subset(kevin2_final, treatment=="c" & response =="vcmax")
kevin2_f_vcmax <- subset(kevin2_final, treatment=="f" & response =="vcmax")
kevin2_cf_vcmax <- subset(kevin2_final, treatment=="cf" & response =="vcmax")

kevin2_c_jmax <- subset(kevin2_final, treatment=="c" & response =="jmax")
kevin2_f_jmax <- subset(kevin2_final, treatment=="f" & response =="jmax")
kevin2_cf_jmax <- subset(kevin2_final, treatment=="cf" & response =="jmax")

kevin_vcmax_plotmean <- agg_meta_sen_coef(kevin2_c_vcmax,"vcmax")

kevin_jmax_plotmean <- agg_meta_sen_coef(kevin2_c_jmax,"jmax")

kevin_vj <- merge(kevin_vcmax_plotmean,kevin_jmax_plotmean,by=c("exp"),all.x=TRUE)
kevin_vj_ecosystem <- merge(kevin_vj,kevin_ecosystem[,c("exp","ecosystem")],by=c("exp"),all.x=TRUE)

#only select sites that primarily measures vcmax (and then a subset measuring jmax).
#there are two plots (riceface_japan_a_2003_3938_14057_c and riceface_japan_a_2004_3938_14057_c) measuring jmax and not measuring vcmax. Quite surprising.
# by checking their original source (https://www.jstor.org/stable/42796062?seq=6#metadata_info_tab_contents)
# jmax and vcmax was available in figure but not in numbers. So not filling them here.

#Smith's data
smith_co2 <- read.csv("~/data/smith_keenan_gcb/gcb_co2/co2_data.csv")
smith_co2$exp_nam <- smith_co2$SiteID
smith_co2 <- subset(smith_co2,is.na(lat)==FALSE) # remove empty column
summary(smith_co2)
smith_final <- smith_co2%>% group_by(exp_nam,Species)  %>% summarise(number = n())

# check if first and second column ALWAYS has the same site, and species name 
# if so then we can apply a simple for-loop function (like below) to combine data 
for (i in 1:((nrow(smith_co2)/2))){  #if first and second column has the same site, and species name then calculation can run
  if (smith_co2$SiteID[2*i-1] == smith_co2$SiteID[2*i] &
      smith_co2$Genus[2*i-1] == smith_co2$Genus[2*i] &
      smith_co2$Species[2*i-1] == smith_co2$Species[2*i]
  ){ print ("yes")
  } else { print ("no")
  }
}

#first 2 columns are the same treatment in same plot (amb + elv) - so we can calculate it directly to combine them into one column
for (i in 1:((nrow(smith_co2)/2))){
  smith_co2$logr_v[2*i-1] <- NA
  smith_co2$sen_coef_v[2*i-1] <- NA
  smith_co2$logr_j[2*i-1] <- NA
  smith_co2$sen_coef_j[2*i-1] <- NA
  #smith_co2$logr_jv[2*i-1] <- NA
  #smith_co2$sen_coef_jv[2*i-1] <- NA
  
  #fill 1st line data
  smith_co2$ambient_co2[2*i-1] <- NA
  smith_co2$ambient_vcmax[2*i-1] <- NA
  smith_co2$ambient_vcmaxSE[2*i-1] <- NA
  smith_co2$ambient_jmax[2*i-1] <- NA
  smith_co2$ambient_jmaxSE[2*i-1] <-NA
  smith_co2$ambient_vcmax_tmean[2*i-1] <- NA
  smith_co2$ambient_jmax_tmean[2*i-1] <- NA
  #smith_co2$ambient_jv_tmean[2*i-1] <- NA
  
  smith_co2$logr_v[2*i] <- log(smith_co2$Vcmax[2*i]/smith_co2$Vcmax[2*i-1])
  smith_co2$sen_coef_v[2*i] <- log(smith_co2$Vcmax[2*i]/smith_co2$Vcmax[2*i-1])/log(smith_co2$Ca[2*i]/smith_co2$Ca[2*i-1])
  smith_co2$logr_j[2*i] <- log(smith_co2$Jmax[2*i]/smith_co2$Jmax[2*i-1])
  smith_co2$sen_coef_j[2*i] <- log(smith_co2$Jmax[2*i]/smith_co2$Jmax[2*i-1])/log(smith_co2$Ca[2*i]/smith_co2$Ca[2*i-1])
  #smith_co2$logr_jv[2*i] <- log((smith_co2$Jmax[2*i]/smith_co2$Vcmax[2*i])/(smith_co2$Jmax[2*i-1]/smith_co2$Vcmax[2*i-1]))
  #smith_co2$sen_coef_jv[2*i] <- log((smith_co2$Jmax[2*i]/smith_co2$Vcmax[2*i])/(smith_co2$Jmax[2*i-1]/smith_co2$Vcmax[2*i-1]))/log(smith_co2$Ca[2*i]/smith_co2$Ca[2*i-1])
  
  #fill 1st line data
  smith_co2$ambient_co2[2*i] <- smith_co2$Ca[2*i-1]
  smith_co2$ambient_vcmax[2*i] <- smith_co2$Vcmax[2*i-1]
  smith_co2$ambient_vcmaxSE[2*i] <- smith_co2$SE[2*i-1]
  smith_co2$ambient_jmax[2*i] <- smith_co2$Jmax[2*i-1]
  smith_co2$ambient_jmaxSE[2*i] <-smith_co2$Jmax_SE[2*i-1]
  smith_co2$ambient_vcmax_tmean[2*i] <- smith_co2$Vcmax_tmean[2*i-1]
  smith_co2$ambient_jmax_tmean[2*i] <- smith_co2$Jmax_tmean[2*i-1]
  #smith_co2$ambient_jv_tmean[2*i] <- smith_co2$JV_tmean[2*i-1]
}
#wheat has two co2_a - so jmax/vcmax sensitivity coefficient should not be jmax - vcmax. But not a problem since this doens't have jmax data

smith_co2_final <- subset(smith_co2,is.na(logr_v)==FALSE)

vcmax25_warmingco2_siteinfo <- rename(smith_co2_final, c(ambient = ambient_vcmax, elevated=Vcmax, ambient_Se=ambient_vcmaxSE, elevated_Se=SE,
                                                         co2_a=ambient_co2, co2_e=Ca,
                                                         ambient_j = ambient_jmax, elevated_j =Jmax, ambient_Se_j=ambient_jmaxSE, elevated_Se_j=Jmax_SE))
#check N fix species (will be used at the end)
unique(subset(vcmax25_warmingco2_siteinfo,Nac=="Nfix")$SiteID)

vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$Stature=="Tree"] <- "forest"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="ChinaRiceFACE"] <- "cropland" #Oryza	sativa
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="UAMACFACE"] <- "cropland" #Triticum aestivum
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="GlycineCE"] <- "cropland" # Glycine max
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="GossypiumGC"] <- "cropland" # Gossypium hirsutum
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="MississippiCE"] <- "cropland" # Gossypium hirsutum
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="WheatCE"] <- "cropland" # Triticum	aestivum
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="ChinaMiniFACE"] <- "cropland" # Glycine Max
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="GrassOTC"] <- "grassland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="PHACE"] <- "grassland" # Pascopyrium smithii and Linaria	dalmatica
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="FescueCE"] <- "grassland" # Festuca arundinacea
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="NZFACE"] <- "grassland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="NevadaFACE"] <- "grassland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="BioCON"] <- "grassland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="PhalarisGC"] <- "grassland" # Phalaris	arundinacea
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="SwissFACE"] <- "grassland"

#now, combine with Kevin:
smith_all_simple <- subset(vcmax25_warmingco2_siteinfo,exp_nam!="NZFACE" &exp_nam!="BioCON" & exp_nam!="ChinaRiceFACE" & exp_nam!="DukeFACE" & exp_nam!="EUROPOPFACE" & exp_nam!="NevadaFACE" & exp_nam!="SwissFACE")

smith_all_plotmean <- aggregate(smith_all_simple,by=list(smith_all_simple$SiteID,smith_all_simple$ecosystem), FUN=mean, na.rm=TRUE)[,c("Group.1","Group.2","sen_coef_v","sen_coef_j")]
names(smith_all_plotmean) <- c("exp","ecosystem","vcmax","jmax")
smith_all_plotmean$exp[smith_all_plotmean$exp=="OakOTC"] <- "mi_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="ORNL"] <- "ornerp_liqui_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="PineOTC"] <- "christchurch_pr_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="RichmondGH"] <- "australiacotton_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="UIAPine"] <- "ua_otc_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="PHACE"] <- "phace_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="AspenFACE"] <- "aspenface_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="BilyKriz"] <- "bily_kriz_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="Headley"] <- "headley_qp_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="Viesalm"] <- "vielsalm_c"

#check if smith and GCME have consistent site (but not sometimes consistent species!)
unique(subset(smith_all_simple,exp_nam=="OakOTC")[,c("lon","lat","Genus","Species")])
unique(subset(kevin_othervars,exp=="mi_c")[,c("lon","lat","dominant_species")])

unique(subset(smith_all_simple,exp_nam=="ORNL")[,c("lon","lat","Genus","Species")])
unique(subset(kevin_othervars,exp=="ornerp_liqui_c")[,c("lon","lat","dominant_species")])

unique(subset(smith_all_simple,exp_nam=="PineOTC")[,c("lon","lat","Genus","Species")])
unique(subset(kevin_othervars,exp=="christchurch_pr_c")[,c("lon","lat","dominant_species")])

unique(subset(smith_all_simple,exp_nam=="RichmondGH")[,c("lon","lat","Genus","Species")])
unique(subset(kevin_othervars,exp=="australiacotton_c")[,c("lon","lat","dominant_species")])

unique(subset(smith_all_simple,exp_nam=="UIAPine")[,c("lon","lat","Genus","Species")])
unique(subset(kevin_othervars,exp=="ua_otc_c")[,c("lon","lat","dominant_species")])

unique(subset(smith_all_simple,exp_nam=="PHACE")[,c("lon","lat","Genus","Species")])
unique(subset(kevin_othervars,exp=="phace_c")[,c("lon","lat","dominant_species")])

unique(subset(smith_all_simple,exp_nam=="AspenFACE")[,c("lon","lat","Genus","Species")])
unique(subset(kevin_othervars,exp=="aspenface_c")[,c("lon","lat","dominant_species")])

unique(subset(smith_all_simple,exp_nam=="BilyKriz")[,c("lon","lat","Genus","Species")])
unique(subset(kevin_othervars,exp=="bily_kriz_c")[,c("lon","lat","dominant_species")])

unique(subset(smith_all_simple,exp_nam=="Headley")[,c("lon","lat","Genus","Species")])
unique(subset(kevin_othervars,exp=="headley_qp_c")[,c("lon","lat","dominant_species")])

unique(subset(smith_all_simple,exp_nam=="Viesalm")[,c("lon","lat","Genus","Species")])
unique(subset(kevin_othervars,exp=="vielsalm_c")[,c("lon","lat","dominant_species")])

smith_all_plotmean$exp <- tolower(smith_all_plotmean$exp)

smith_all_plotmean <- smith_all_plotmean[,c("exp","vcmax","jmax","ecosystem")]

obs_co2 <- rbind(kevin_vj_ecosystem,smith_all_plotmean)

#N fertilization data from Walker 
walker <- read.csv("~/data/leaf_traits/Walker/LEAF_PHOTOSYNTHESIS_TRAITS_1224/data/Leaf_Photosynthesis_Traits.csv")

Nfer_vcmax <- subset(walker,Treatment_N!="not applicable" &Treatment_CO2!="not applicable" & Vcmax>0)
Nfer_vcmax2 <- Nfer_vcmax[,c("lon","lat","species","Treatment_N","Treatment_P",
                             "Treatment_CO2","Treatment_light","Vcmax","Jmax","narea","SLA")]
summary(Nfer_vcmax2)
Nfer_vcmax2$SLA[Nfer_vcmax2$SLA==-9999.999] <- NA
Nfer_vcmax2$LMA <- 1/Nfer_vcmax2$SLA
Nfer_vcmax2$nmass <- Nfer_vcmax2$narea/Nfer_vcmax2$LMA

Nfer_vcmax3 <- Nfer_vcmax2 %>% 
  group_by(lon,lat,species,Treatment_N,Treatment_P,Treatment_CO2,Treatment_light) %>%
  summarise(vcmax = mean(Vcmax,na.rm=TRUE),jmax = mean(Jmax,na.rm=TRUE),nmass = mean(nmass,na.rm=TRUE),
            narea = mean(narea,na.rm=TRUE),LMA = mean(LMA,na.rm=TRUE)) 

#keep all other treatments consistent
highN_vcmax <- merge(subset(Nfer_vcmax3,Treatment_N=="high"&Treatment_CO2=="amb"),
                     subset(Nfer_vcmax3,Treatment_N=="high"&Treatment_CO2=="ele"),
                     by=c("lon","lat","species","Treatment_N","Treatment_P","Treatment_light"),all.x=TRUE)
highN_vcmax <- subset(highN_vcmax,Treatment_light!="shade")

dim(highN_vcmax) #now, except co2 condition, all others are consistent

highN_vcmax$logr_vcmax <- log(highN_vcmax$vcmax.y/highN_vcmax$vcmax.x)
highN_vcmax$logr_jmax <- log(highN_vcmax$jmax.y/highN_vcmax$jmax.x)
highN_vcmax$logr_nmass <- log(highN_vcmax$nmass.y/highN_vcmax$nmass.x)
highN_vcmax$logr_narea <- log(highN_vcmax$narea.y/highN_vcmax$narea.x)
highN_vcmax$logr_LMA <- log(highN_vcmax$LMA.y/highN_vcmax$LMA.x)

highN_vcmax$site_species <- paste(highN_vcmax$lon,highN_vcmax$lat,highN_vcmax$species,sep="_")

highN_vcmax

lowN_vcmax <- merge(subset(Nfer_vcmax3,Treatment_N=="low"&Treatment_CO2=="amb"),
                    subset(Nfer_vcmax3,Treatment_N=="low"&Treatment_CO2=="ele"),
                    by=c("lon","lat","species","Treatment_N","Treatment_P","Treatment_light"),all.x=TRUE)

lowN_vcmax <- subset(lowN_vcmax,Treatment_light!="shade")

#select the same site-species paralled for high-N data
lowN_vcmax$lon;highN_vcmax$lon
lowN_vcmax <- subset(lowN_vcmax,lon!=-71.03 & lon!=141.00)

lowN_vcmax$logr_vcmax <- log(lowN_vcmax$vcmax.y/lowN_vcmax$vcmax.x)
lowN_vcmax$logr_jmax <- log(lowN_vcmax$jmax.y/lowN_vcmax$jmax.x)
lowN_vcmax$logr_nmass <- log(lowN_vcmax$nmass.y/lowN_vcmax$nmass.x)
lowN_vcmax$logr_narea <- log(lowN_vcmax$narea.y/lowN_vcmax$narea.x)
lowN_vcmax$logr_LMA <- log(lowN_vcmax$LMA.y/lowN_vcmax$LMA.x)
lowN_vcmax$site_species <- paste(lowN_vcmax$lon,lowN_vcmax$lat,lowN_vcmax$species,sep="_")

low_high_vcmax <-  dplyr::bind_rows(highN_vcmax,lowN_vcmax)
#except for treatment N, all other treatment are consistent
unique(low_high_vcmax[,c("site_species","Treatment_N","Treatment_P","Treatment_light","Treatment_CO2.x")])

low_high_vcmax$co2_change <- NA

#input co2 change
low_high_vcmax$co2_change[low_high_vcmax$lon==-84.04] <- log(560/360) 
#https://pubmed.ncbi.nlm.nih.gov/11960756/
low_high_vcmax$co2_change[low_high_vcmax$lon==20.00] <- log(700/350) 
#https://www.researchgate.net/publication/242459413_Nutrient_and_genotypic_effects_on_CO2-responsiveness_Photosynthetic_regulation_in_Leucadendron_species_of_a_nutrient-poor_environment
low_high_vcmax$co2_change[low_high_vcmax$lon==89.14] <- log(720/360) 

low_high_vcmax$vcmax <- low_high_vcmax$logr_vcmax/low_high_vcmax$co2_change
low_high_vcmax$jmax <- low_high_vcmax$logr_jmax/low_high_vcmax$co2_change
low_high_vcmax$nmass <- low_high_vcmax$logr_nmass/low_high_vcmax$co2_change
low_high_vcmax$narea <- low_high_vcmax$logr_narea/low_high_vcmax$co2_change
low_high_vcmax$LMA <- low_high_vcmax$logr_LMA/low_high_vcmax$co2_change

low_high_dataset <- low_high_vcmax[,c("lon","lat","site_species","Treatment_N",
                                      "vcmax","jmax","nmass","narea","LMA")]
low_high_dataset

#####2. Input light effect data (we only obtain vcmax and jmax)
walker <- read.csv("~/data/leaf_traits/Walker/LEAF_PHOTOSYNTHESIS_TRAITS_1224/data/Leaf_Photosynthesis_Traits.csv")

light_vcmax <- subset(walker,(Treatment_light=="shade"|Treatment_light=="sun"|
                                Treatment_light=="high"|Treatment_light=="low") & Vcmax>0)
light1 <- light_vcmax[,c("lon","lat","species","Treatment_N","Treatment_P",
                         "Treatment_CO2","Treatment_light","Vcmax","Jmax")]
light1$Jmax[light1$Jmax==-9999.99] <- NA

#divided into two categories:
#(1) shade vs. sun
#(2) high vs. low light

light1a <- merge(subset(light1,Treatment_light=="shade"),
                 subset(light1,Treatment_light=="sun"),
                 by=c("lon","lat","species","Treatment_N","Treatment_P","Treatment_CO2"),
                 all.x=TRUE)

#remove ele co2 and high-N
light1a <- subset(light1a,Treatment_CO2 != "ele" & Treatment_N!="high" & Treatment_N!="fertilised")
light1a

light1a$logr_vcmax <- log(light1a$Vcmax.y/light1a$Vcmax.x)
light1a$logr_jmax <- log(light1a$Jmax.y/light1a$Jmax.x)

light1a$type_name <- "shade_to_sun"

light2 <- subset(light_vcmax,Treatment_light=="high"|Treatment_light=="low")
light2 <- light2[,c("lon","lat","species","Treatment_N","Treatment_P",
                         "Treatment_CO2","Treatment_light","Vcmax","Jmax")]
light2a <- merge(subset(light2,Treatment_light=="low"),subset(light2,Treatment_light=="high"),
                 by=c("lon","lat","species","Treatment_N","Treatment_P","Treatment_CO2"),
                 all.x=TRUE)
light2a
light2a$logr_vcmax <- log(light2a$Vcmax.y/light2a$Vcmax.x)
light2a$logr_jmax <- log(light2a$Jmax.y/light2a$Jmax.x)
light2a$type_name <- "low_to_high_light"

light_vcmax_points <- dplyr::bind_rows(light1a[,c("lon","lat","species","logr_vcmax","logr_jmax","type_name")],
                                       light2a[,c("lon","lat","species","logr_vcmax","logr_jmax","type_name")]) 
names(light_vcmax_points) <- c("lon","lat","species","vcmax","jmax","type_name")

light_vcmax_points$site_species <- paste(light_vcmax_points$lon,light_vcmax_points$lat,light_vcmax_points$species,sep="_")

obs_light <- light_vcmax_points %>% group_by(lon,lat,site_species,type_name)  %>% summarise(vcmax = mean(vcmax),jmax=mean(jmax))


####Input Warming data 
#Input data, and combine dataset to get Tleaf
df3 <- read_csv("~/data/Kumarathunge_2020_newphy/kumarathunge_2020_newphy/Data/PPC-TGlob_V1.0.csv")
df3 <- subset(df3,Growth_condition=="Glasshouse"|Growth_condition=="Whole tree chambers")
# select sites - only these sites have available Temp information and matched our purpose.
df3 <- df3 %>% filter(Dataset %in% c("Black Spruce, ON, Canada",
                                     "Corymbia calophylla provs, AU-NSW",
                                     "Eucalyptus globulus, AU-NSW",
                                     "Eucalyptus parramattensis, AU-NSW",
                                     "Eucalyptus tereticornis provs AU-NSW",
                                     "Eucalyptus tereticornis, AU-NSW"))

unique(df3$Leafspan) # all in forest - evergreen

df3 %>% group_by(Dataset,Temp_Treatment)  %>% summarise(number = n())
#here one site have three temp treatments, remove the 'highest' temperature 
df3 <- subset(df3,Temp_Treatment!=6)

df3 %>% group_by(Dataset,Temp_Treatment)  %>% summarise(number = n())

#calculate Dushan's Vcmax25
df3$vcmax25 <- df3$Vcmax*exp((65330/8.314)*((1/(df3$Tleaf+273.15))-(1/298.15)))
df3$jmax25 <- df3$Jmax*exp((43540/8.314)*((1/(df3$Tleaf+273.15))-(1/298.15)))


#now, fill temp treatment information - needs to check original data
df1 <- read_csv("~/data/Kumarathunge_2020_newphy/kumarathunge_2020_newphy/Data/ACi-TGlob_V1.0.csv")
df1 <- subset(df1,Growth_condition=="Glasshouse"|Growth_condition=="Whole tree chambers")

#For "Black Spruce, ON, Canada" - temperature was given in its original A-Ci data
unique(subset(df1,Dataset=="Black Spruce, ON, Canada")$Temp_Treatment) 
df3$Temp_Treatment[df3$Dataset=="Black Spruce, ON, Canada" & df3$Temp_Treatment=="ambient"] <- 22
df3$Temp_Treatment[df3$Dataset=="Black Spruce, ON, Canada" & df3$Temp_Treatment=="elevated"] <- 30

#For "Corymbia calophylla provs, AU-NSW"- temperature was given in its original A-Ci data
unique(subset(df1,Dataset=="Corymbia calophylla provs, AU-NSW")$Temp_Treatment) 
df3$Temp_Treatment[df3$Dataset=="Corymbia calophylla provs, AU-NSW" & df3$Temp_Treatment=="ambient"] <- 26
df3$Temp_Treatment[df3$Dataset=="Corymbia calophylla provs, AU-NSW" & df3$Temp_Treatment=="elevated"] <- 32

#For "Eucalyptus tereticornis provs AU-NSW" - temperature was given in its original A-Ci data
#needs to check ID_curve number spefically, since it has multiple temp control
#tempeature information was also given in /Users/yunpeng/data/Kumarathunge_2020_newphy/kumarathunge_2020_newphy/dataset_processing_code.R
#showing that 1 = 18C, 4=28.5C
id_one <- subset(df3,Dataset=="Eucalyptus tereticornis provs AU-NSW")$Curve_Id
unique((subset(df1,Dataset=="Eucalyptus tereticornis provs AU-NSW")  %>% filter(Curve_Id %in% id_one))$Temp_Treatment)
df3$Temp_Treatment[df3$Dataset=="Eucalyptus tereticornis provs AU-NSW" & df3$Temp_Treatment=="1"] <- 18
df3$Temp_Treatment[df3$Dataset=="Eucalyptus tereticornis provs AU-NSW" & df3$Temp_Treatment=="4"] <- 28.5

#For "Eucalyptus globulus, AU-NSW"
#see below Table s1: we can assume ambient = 19.7, elevated = 24.9 (basing on growth temperature column)
#https://nph.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fnph.15668&file=nph15668-sup-0001-SupInfo.pdf
df3$Temp_Treatment[df3$Dataset=="Eucalyptus globulus, AU-NSW" & df3$Temp_Treatment=="ambient"] <- 19.7
df3$Temp_Treatment[df3$Dataset=="Eucalyptus globulus, AU-NSW" & df3$Temp_Treatment=="elevated"] <- 24.9

#for "Eucalyptus parramattensis, AU-NSW"
#see below Table s1: we can assume ambient = 15, elevated = 18 (basing on growth temperature column)
#https://nph.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fnph.15668&file=nph15668-sup-0001-SupInfo.pdf
df3$Temp_Treatment[df3$Dataset=="Eucalyptus parramattensis, AU-NSW" & df3$Temp_Treatment=="ambient"] <- 15.1
df3$Temp_Treatment[df3$Dataset=="Eucalyptus parramattensis, AU-NSW" & df3$Temp_Treatment=="elevated"] <- 18.0

#for Eucalyptus tereticornis, AU-NSW
#see below Table s1: we can assume ambient = 14.4, elevated = 22.6 (basing on growth temperature column)
#https://nph.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fnph.15668&file=nph15668-sup-0001-SupInfo.pdf
df3$Temp_Treatment[df3$Dataset=="Eucalyptus tereticornis, AU-NSW" & df3$Temp_Treatment=="ambient"] <- 14.4
df3$Temp_Treatment[df3$Dataset=="Eucalyptus tereticornis, AU-NSW" & df3$Temp_Treatment=="elevated"] <- 22.6

df3$Temp_Treatment <- as.numeric(df3$Temp_Treatment)

#summary
df3 %>% group_by(Dataset,Temp_Treatment)  %>% summarise(number = n())


#create a dataframe

#"Black Spruce, ON, Canada": measurement year: 2006.
#basing on https://onlinelibrary.wiley.com/doi/10.1111/j.1365-3040.2008.01842.x
d1 <- tibble(
  exp= "Black Spruce, ON, Canada",
  lon = subset(df3,Dataset=="Black Spruce, ON, Canada")$seed_source_longitude[1],
  lat = subset(df3,Dataset=="Black Spruce, ON, Canada")$seed_source_latitude[1],
  vcmax25_a= mean(subset(df3,Dataset=="Black Spruce, ON, Canada" & Temp_Treatment==22)$vcmax25,na.rm=TRUE),
  vcmax25_e= mean(subset(df3,Dataset=="Black Spruce, ON, Canada" & Temp_Treatment==30)$vcmax25,na.rm=TRUE),
  jmax25_a= mean(subset(df3,Dataset=="Black Spruce, ON, Canada" & Temp_Treatment==22)$jmax25,na.rm=TRUE),
  jmax25_e= mean(subset(df3,Dataset=="Black Spruce, ON, Canada" & Temp_Treatment==30)$jmax25,na.rm=TRUE),
  temp_a = 22,
  temp_e = 30,
  year_start=2006,
  year_end=2006,
  t1 = 8)

#"Corymbia calophylla provs, AU-NSW"
unique(subset(df3,Dataset=="Corymbia calophylla provs, AU-NSW")$Date)
d2 <- tibble(
  exp= "Corymbia calophylla provs, AU-NSW",
  lon = subset(df3,Dataset=="Corymbia calophylla provs, AU-NSW")$seed_source_longitude[1],
  lat = subset(df3,Dataset=="Corymbia calophylla provs, AU-NSW")$seed_source_latitude[1],
  vcmax25_a= mean(subset(df3,Dataset=="Corymbia calophylla provs, AU-NSW" & Temp_Treatment==26)$vcmax25,na.rm=TRUE),
  vcmax25_e= mean(subset(df3,Dataset=="Corymbia calophylla provs, AU-NSW" & Temp_Treatment==32)$vcmax25,na.rm=TRUE),
  jmax25_a= mean(subset(df3,Dataset=="Corymbia calophylla provs, AU-NSW" & Temp_Treatment==26)$jmax25,na.rm=TRUE),
  jmax25_e= mean(subset(df3,Dataset=="Corymbia calophylla provs, AU-NSW" & Temp_Treatment==32)$jmax25,na.rm=TRUE),
  temp_a = 26,
  temp_e = 32,
  year_start=2014,
  year_end=2014,
  t1 = 6)

#"Eucalyptus globulus, AU-NSW"
unique(subset(df3,Dataset=="Eucalyptus globulus, AU-NSW")$Date)

d3 <- tibble(
  exp= "Eucalyptus globulus, AU-NSW",
  lon = subset(df3,Dataset=="Eucalyptus globulus, AU-NSW")$seed_source_longitude[1],
  lat = subset(df3,Dataset=="Eucalyptus globulus, AU-NSW")$seed_source_latitude[1],
  vcmax25_a= mean(subset(df3,Dataset=="Eucalyptus globulus, AU-NSW" & Temp_Treatment==20)$vcmax25,na.rm=TRUE),
  vcmax25_e= mean(subset(df3,Dataset=="Eucalyptus globulus, AU-NSW" & Temp_Treatment==23)$vcmax25,na.rm=TRUE),
  jmax25_a= mean(subset(df3,Dataset=="Eucalyptus globulus, AU-NSW" & Temp_Treatment==20)$jmax25,na.rm=TRUE),
  jmax25_e= mean(subset(df3,Dataset=="Eucalyptus globulus, AU-NSW" & Temp_Treatment==23)$jmax25,na.rm=TRUE),
  temp_a = 20,
  temp_e = 23,
  year_start=2010,
  year_end=2011,
  t1 = 3)

#"Eucalyptus globulus, AU-NSW"
unique(subset(df3,Dataset=="Eucalyptus parramattensis, AU-NSW")$Date)

d4 <- tibble(
  exp= "Eucalyptus parramattensis, AU-NSW",
  lon = subset(df3,Dataset=="Eucalyptus parramattensis, AU-NSW")$seed_source_longitude[1],
  lat = subset(df3,Dataset=="Eucalyptus parramattensis, AU-NSW")$seed_source_latitude[1],
  vcmax25_a= mean(subset(df3,Dataset=="Eucalyptus parramattensis, AU-NSW" & Temp_Treatment==15)$vcmax25,na.rm=TRUE),
  vcmax25_e= mean(subset(df3,Dataset=="Eucalyptus parramattensis, AU-NSW" & Temp_Treatment==18)$vcmax25,na.rm=TRUE),
  jmax25_a= mean(subset(df3,Dataset=="Eucalyptus parramattensis, AU-NSW" & Temp_Treatment==15)$jmax25,na.rm=TRUE),
  jmax25_e= mean(subset(df3,Dataset=="Eucalyptus parramattensis, AU-NSW" & Temp_Treatment==18)$jmax25,na.rm=TRUE),
  temp_a = 15,
  temp_e = 18,
  year_start=2016,
  year_end=2016,
  t1 = 3)

#"Eucalyptus tereticornis provs AU-NSW"
unique(subset(df3,Dataset=="Eucalyptus tereticornis provs AU-NSW")$Date)
d5 <- tibble(
  exp= "Eucalyptus tereticornis provs AU-NSW",
  lon = subset(df3,Dataset=="Eucalyptus tereticornis provs AU-NSW")$seed_source_longitude[1],
  lat = subset(df3,Dataset=="Eucalyptus tereticornis provs AU-NSW")$seed_source_latitude[1],
  vcmax25_a= mean(subset(df3,Dataset=="Eucalyptus tereticornis provs AU-NSW" & Temp_Treatment==18)$vcmax25,na.rm=TRUE),
  vcmax25_e= mean(subset(df3,Dataset=="Eucalyptus tereticornis provs AU-NSW" & Temp_Treatment==28.5)$vcmax25,na.rm=TRUE),
  jmax25_a= mean(subset(df3,Dataset=="Eucalyptus tereticornis provs AU-NSW" & Temp_Treatment==18)$jmax25,na.rm=TRUE),
  jmax25_e= mean(subset(df3,Dataset=="Eucalyptus tereticornis provs AU-NSW" & Temp_Treatment==28.5)$jmax25,na.rm=TRUE),
  temp_a = 18,
  temp_e = 28.5,
  year_start=2016,
  year_end=2016,
  t1 = 10.5)

#"Eucalyptus tereticornis, AU-NSW"
unique(subset(df3,Dataset=="Eucalyptus tereticornis, AU-NSW")$Date)

d6 <- tibble(
  exp= "Eucalyptus tereticornis, AU-NSW",
  lon = subset(df3,Dataset=="Eucalyptus tereticornis, AU-NSW")$seed_source_longitude[1],
  lat = subset(df3,Dataset=="Eucalyptus tereticornis, AU-NSW")$seed_source_latitude[1],
  vcmax25_a= mean(subset(df3,Dataset=="Eucalyptus tereticornis, AU-NSW" & Temp_Treatment==14.4)$vcmax25,na.rm=TRUE),
  vcmax25_e= mean(subset(df3,Dataset=="Eucalyptus tereticornis, AU-NSW" & Temp_Treatment==22.6)$vcmax25,na.rm=TRUE),
  jmax25_a= mean(subset(df3,Dataset=="Eucalyptus tereticornis, AU-NSW" & Temp_Treatment==14.4)$jmax25,na.rm=TRUE),
  jmax25_e= mean(subset(df3,Dataset=="Eucalyptus tereticornis, AU-NSW" & Temp_Treatment==22.6)$jmax25,na.rm=TRUE),
  temp_a = 14.4,
  temp_e = 22.6,
  year_start=2013,
  year_end=2014,
  t1 = 8.2)

obs_warming <- as.data.frame(rbind(d1,d2,d3,d4,d5,d6))



######now, combine all observational data
obs_co2$treatment <- "co2"

obs_warming$vcmax <-log(obs_warming$vcmax25_e/obs_warming$vcmax25_a)/obs_warming$t1
obs_warming$jmax <-log(obs_warming$jmax25_e/obs_warming$jmax25_a)/obs_warming$t1
obs_warming_final <- obs_warming[,c("exp","lon","lat","vcmax","jmax")]
obs_warming_final$treatment <- "warming"

obs_light$exp <- obs_light$site_species
obs_light_final <- obs_light[,c("exp","lon","lat","vcmax","jmax","type_name")]
obs_light_final$treatment <- "light"

low_high_dataset$exp <- low_high_dataset$site_species
low_high_dataset$treatment[low_high_dataset$Treatment_N=="high"] <- "highN"
low_high_dataset$treatment[low_high_dataset$Treatment_N=="low"] <- "lowN"
low_high_dataset$type_name <- low_high_dataset$treatment

low_high_dataset$ref <- "Walker 2014 dataset"
low_high_dataset$comments <- "vcmax and jmax are sensitivity coefficients of co2 at N fertilization"
low_high_dataset <- low_high_dataset[ , -which(names(low_high_dataset) %in% c("Treatment_N","site_species"))]

#merged prediction data one-by-one
#now, combine with prediction data
prediction <- read.csv("~/data/gcme/prediction/prediction.csv")
names(prediction) <- c("X","exp","lon","lat","pred_vcmax","pred_jmax","pred_jmax_vcmax",
                       "treatment","ref","comments")

#first, merge to get prediction data
obs_co2_pred <- merge(obs_co2[,c("exp","vcmax","jmax","ecosystem")],prediction, by=c("exp"),all.x=TRUE)
obs_warming_pred <- merge(obs_warming[,c("exp","vcmax","jmax")],prediction, by=c("exp"),all.x=TRUE)
obs_light_pred <- merge(obs_light_final[,c("exp","lon","lat","vcmax","jmax","type_name")],subset(prediction,treatment=="light" & exp!="light1"), by=c("lon","lat"),all.x=TRUE) # remove first column (just predicted in different years) since it will cause repeated data
obs_light_pred$exp <- obs_light_pred$exp.x

obs_light_pred <- obs_light_pred[ , -which(names(obs_light_pred) %in% c("exp.x","exp.y"))]

all_obs_pred <- dplyr::bind_rows(obs_co2_pred, obs_warming_pred,obs_light_pred,low_high_dataset) 

all_obs_pred$jmax_vcmax <- all_obs_pred$jmax - all_obs_pred$vcmax

#finally, for meta-analysis
#anpp, bnpp, nmass, LAI, soil N
kevin_othervars
#filter plots only within vcmax and jmax
photo_plot <- c(unique(all_obs_pred$exp),c("duke2_f","euroface4_pa_f","euroface4_pe_f","euroface4_pn_f","new_zealand_face_f"),
                c("duke2_cf","euroface4_pa_cf","euroface4_pe_cf","euroface4_pn_cf","new_zealand_face_cf"))
kevin_othervars_cf <- kevin_othervars %>% filter(exp %in% photo_plot)

#leaf_n needs division - but leaf_p doesn't need (as it is all pmass)
unique(subset(kevin_othervars_cf,response=="leaf_n")$Unit)

nmass <- subset(kevin_othervars_cf,response=="leaf_n") %>%         
  filter(Unit %in% c("g_100g","g_g","g_kg","mg_g","mmol_g","ug_mg","umol_g"))
nmass$response <- "nmass"
logr_c_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="c")));
logr_f_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="f")));
logr_cf_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="cf")))

narea <- subset(kevin_othervars_cf,response=="leaf_n") %>%         
  filter(Unit %in% c("g_m2","mg_cm2","mg_m2","mmol_m2","ug_cm3","umol_m2"))
narea$response <- "narea"

logr_c_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="c")));
logr_f_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="f")));
logr_cf_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="cf")))

#bnpp

#root data 1 - all used
root1 <- subset(kevin_othervars_cf,response=="root_production" & Unit!="g_m2") 
unique(root1[,c("exp","Unit")]) # the unit g_m2 was in brandbjerg_c, which was already been covered by the data with g_m2y

#root data 2: alternatively, use fine root production (but remove non-useful data)
root2 <- subset(kevin_othervars_cf,response=="fine_root_production" & Unit!="km_m3"& Unit!="m_m2")

root2_rest <- subset(root2,exp!="new_zealand_face_c"&exp!="eucface_c"&exp!="mi_c"&exp!="ornerp_liqui_c"&
                       exp!="duke_c"&exp!="biocon_c"&exp!="aspenface_c"&exp!="brandbjerg_c")

unique(root2_rest[,c("exp","response","Unit","citation")])  # all unit works well (finzi_et_al_2007 has unit gc_m2, but it should be the value from the whole year i.e. can be gc_m2_yr)

bnpp_dataset <- rbind(root1,root2_rest)

#soil mineral N in dry-mass: 7+1
#7 plots: nh4 + no3
soil_nh4 <- subset(kevin_othervars_cf,response=="soil_nh4-n") %>% group_by(exp,Unit) %>% summarise(co2_a = mean(co2_a), co2_e = mean(co2_e), ambient = mean(ambient), elevated = mean(elevated))
soil_no3 <- subset(kevin_othervars_cf,response=="soil_no3-n") %>% group_by(exp,Unit) %>% summarise(co2_a = mean(co2_a), co2_e = mean(co2_e),ambient = mean(ambient), elevated = mean(elevated))
soil_nh4no3 <- na.omit(merge(soil_nh4,soil_no3,by=c("exp","Unit","co2_a","co2_e"),all.x=TRUE))
soil_nh4no3
soil_nh4no3$soil_mineral_N <- log((soil_nh4no3$elevated.x+soil_nh4no3$elevated.y)/(soil_nh4no3$ambient.x+soil_nh4no3$ambient.y))/log(soil_nh4no3$co2_e/soil_nh4no3$co2_a)
# calculate duke_2_cf
duke_2_cf_nh4no3 <- (soil_nh4no3$elevated.x[soil_nh4no3$exp=="duke2_cf"]+soil_nh4no3$elevated.y[soil_nh4no3$exp=="duke2_cf"])
duke_2_f_nh4no3 <- (soil_nh4no3$elevated.x[soil_nh4no3$exp=="duke2_f"]+soil_nh4no3$elevated.y[soil_nh4no3$exp=="duke2_f"])
co2_a <- soil_nh4no3$co2_a[soil_nh4no3$exp=="duke2_cf"]
co2_e <- soil_nh4no3$co2_e[soil_nh4no3$exp=="duke2_cf"]

duke2_cf_soil <- log(duke_2_cf_nh4no3/duke_2_f_nh4no3)/log(co2_e/co2_a) #will be added at the end

soil_nh4no3$soil_mineral_N[soil_nh4no3$exp=="duke2_cf"] <- duke2_cf_soil

soil_mineral_1 <- subset(soil_nh4no3, exp!="duke2_f")[,c("exp","soil_mineral_N")]

#combine from multiple samples to site-mean
soil_mineral_1 <- soil_mineral_1%>% group_by(exp)  %>% summarise(soil_mineral_N = mean(soil_mineral_N))


#1 plot: soil_in
soil_mineral_3 <- tibble(exp="facts_ii_face3_pt_c",
                         soil_mineral_N =agg_meta_sen_coef(response_ratio_v2(subset(kevin_othervars_cf,response=="soil_in")),
                                                           "soil_mineral_N")$soil_mineral_N)
soil_mineral_dry <- rbind(soil_mineral_1,soil_mineral_3)

#LMA
kevin_LMA <- read.csv("~/data/gcme/kevin/orig_leaf/LMA.csv")
kevin_LMA <- rename(kevin_LMA, c(ambient = x_c, elevated=x_t, ambient_Sd=sd_c, elevated_Sd=sd_t,ambient_Se=se_c,elevated_Se=se_t,n_plots=rep_c,
                                 z=elevation, co2_a=c_c, co2_e=c_t, nfertQ_a = n_c, nfertQ_e = n_t, pfertQ_a = p_c, pfertQ_e = p_t,kfertQ_a = k_c, kfertQ_e = k_t,
                                 warmQ_e1 = w_t1, warmQ_e2 = w_t2, warmQ_e3 = w_t3, Unit=x_units))
kevin_LMA$ambient_Se <- as.numeric(kevin_LMA$ambient_Se)
kevin_LMA$elevated_Se <- as.numeric(kevin_LMA$elevated_Se)

kevin_LMA <- kevin_LMA %>% filter(exp %in% photo_plot)

unique(kevin_LMA$response)

sla <- subset(kevin_LMA,response=="sla")
unique(sla$Unit)
lma <- subset(kevin_LMA,response=="lma")
unique(lma$Unit)
#convert SLA to LMA format 
sla$ambient <- 1/sla$ambient; sla$elevated <- 1/sla$elevated

LMA <- dplyr::bind_rows(lma[,c("ambient","elevated","exp","treatment","co2_a","co2_e","Unit")],sla[,c("ambient","elevated","exp","treatment","co2_a","co2_e","Unit")])
LMA$response <- "LMA"
LMA2 <- merge(LMA,unique(kevin_othervars_cf[,c("exp","exp_nam")]),by=c("exp"),all.x=TRUE)
logr_c_LMA <- as_tibble(response_ratio_v2(subset(LMA2,treatment=="c")))
logr_f_LMA <-  as_tibble(response_ratio_v2(subset(LMA2,treatment=="f")))
logr_cf_LMA <-  as_tibble(response_ratio_v2(subset(LMA2,treatment=="cf")))


#######final aggregation 
lma_plot <- agg_meta_sen_coef(logr_c_LMA,"LMA")
narea_plot <- agg_meta_sen_coef(logr_c_narea,"narea")
nmass_plot <- agg_meta_sen_coef(logr_c_nmass,"nmass")
logr_c_bnpp <- response_ratio_v2(bnpp_dataset)
bnpp_plot <- agg_meta_sen_coef(logr_c_bnpp,"bnpp")
#some data needs manually added because co2_a and co2_e is not indicated
bnpp_plot$bnpp[bnpp_plot$exp=="duke_c"] <- log(subset(bnpp_dataset,exp=="duke_c")$elevated/subset(bnpp_dataset,exp=="duke_c")$ambient)/log(563/363)
bnpp_plot$bnpp[bnpp_plot$exp=="brandbjerg_c"] <- log(subset(bnpp_dataset,exp=="brandbjerg_c")$elevated/subset(bnpp_dataset,exp=="brandbjerg_c")$ambient)/log(510/380)
bnpp_plot

#also, some other dataset
logr_c_lai <- response_ratio_v2(subset(kevin_othervars_cf,treatment=="c" & response=="lai"))
logr_f_lai <- response_ratio_v2(subset(kevin_othervars_cf,treatment=="f"& response=="lai"))
logr_cf_lai <- response_ratio_v2(subset(kevin_othervars_cf,treatment=="cf" & response=="lai"))

logr_c_root_shoot_ratio <- response_ratio_v2(subset(kevin_othervars_cf,treatment=="c" & response=="root_shoot_ratio"))
logr_f_root_shoot_ratio <- response_ratio_v2(subset(kevin_othervars_cf,treatment=="f"& response=="root_shoot_ratio"))
logr_cf_root_shoot_ratio <- response_ratio_v2(subset(kevin_othervars_cf,treatment=="cf" & response=="root_shoot_ratio"))

lai_plot <- agg_meta_sen_coef(logr_c_lai,"lai")
root_shoot_ratio_plot <- agg_meta_sen_coef(logr_c_root_shoot_ratio,"root_shoot_ratio")

#anpp
unique(subset(kevin_othervars_cf,response=="anpp")$Unit)
#firstly, only combing with gC/m2/yr data. Remove below unit first.
anpp_new1 <- subset(kevin_othervars_cf,response=="anpp"&Unit!="t_ha"&Unit!="g_m2"&Unit!="gc_m2"&Unit!="mg"&Unit!="g_plant")
unique(anpp_new1[,c("Unit","exp")])

# this data should be included, although it is g/m2 but it is actually flux data in popface
anpp_new2 <- subset(kevin_othervars_cf,response=="anpp"&Unit=="gc_m2"&citation=="finzi_et_al_2007")
anpp_new2
#combine
anpp_new <- rbind(anpp_new1,anpp_new2)

anpp_new_c <- response_ratio_v2(subset(anpp_new,treatment=="c"))
anpp_new_f <- response_ratio_v2(subset(anpp_new,treatment=="f"))
anpp_new_cf <- response_ratio_v2(subset(anpp_new,treatment=="cf"))
anpp_plot1 <- agg_meta_sen_coef(anpp_new_c,"anpp")
anpp_plot1

#second, have a look at 'other unit'
anpp_others <- subset(kevin_othervars_cf,response=="anpp"& (Unit=="t_ha"|Unit=="g_m2"|Unit=="gc_m2"|Unit=="mg"|Unit=="g_plant"))
unique(anpp_others[,c("exp","citation")])
#checked those sites one-by-one (only need to add nevada_desert_face_c and soyfacesoy2_c; while duke2_cf will be added at the end)

#nevada_desert_face_c - checked - unit is mg, these are seed mass, cannot be used
subset(anpp_others,exp=="nevada_desert_face_c") 

#soyfacesoy2_c - checked (https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2486.2005.001017.x)
#See its Fig.3. It is cumulative anpp within one year
#so the data should be at last measurment time (15_10_2002)
subset(anpp_others,exp=="soyfacesoy2_c")[,c("ambient","elevated","sampling_date")] 
soyfacesoy2_c_anpp <- log((subset(anpp_others,exp=="soyfacesoy2_c"&sampling_date=="15_10_2002")$elevated)/
                            (subset(anpp_others,exp=="soyfacesoy2_c"&sampling_date=="15_10_2002")$ambient))/
  log(subset(anpp_others,exp=="soyfacesoy2_c"&sampling_date=="15_10_2002")$co2_e/subset(anpp_others,exp=="soyfacesoy2_c"&sampling_date=="15_10_2002")$co2_a)
anpp_plot3 <- tibble(exp="soyfacesoy2_c",anpp=soyfacesoy2_c_anpp)

#combine
anpp_plot <- rbind(anpp_plot1,anpp_plot3)

all_obs_pred <- all_obs_pred %>% rename(condition = treatment)

final_mean <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),
                    list(all_obs_pred,lma_plot,narea_plot,nmass_plot,bnpp_plot,root_shoot_ratio_plot,
                         soil_mineral_dry,
                         anpp_plot,lai_plot))

#combine data and remove not-useful column
final_mean$narea.y[is.na(final_mean$narea.x)==F] <- final_mean$narea.x[is.na(final_mean$narea.x)==F]
final_mean$nmass.y[is.na(final_mean$nmass.x)==F] <- final_mean$nmass.x[is.na(final_mean$nmass.x)==F]
final_mean$LMA.y[is.na(final_mean$LMA.x)==F] <- final_mean$LMA.x[is.na(final_mean$LMA.x)==F]

final_mean$narea <- final_mean$narea.y
final_mean$nmass <- final_mean$nmass.y
final_mean$LMA <- final_mean$LMA.y

final_mean <- final_mean[ , -which(names(final_mean) %in% c("narea.y","nmass.y","LMA.y",
                                                            "narea.x","nmass.x","LMA.x","X"))]

#create some new points to add co2-effect on vcmax (with N fertilization)
t1 <- tibble(exp="euroface4_pa_cf",condition="Fertilization")
t2 <- tibble(exp="euroface4_pe_cf",condition="Fertilization")
t3 <- tibble(exp="euroface4_pn_cf",condition="Fertilization")
t4 <- tibble(exp="new_zealand_face_cf",condition="Fertilization")
t5 <- tibble(exp="duke2_cf",condition="Fertilization")

final5 <-  dplyr::bind_rows(final_mean,t1,t2,t3,t4,t5)


#create a function to calculate co2 effect on N fertilzation 
#cal_nfer - merged by individiduals, and then averaged by inviduals to calculate logr
#if sucessfully merged (then it means ambient values are different from each other), then n(row) before and after should be the same
#if not then(then it means ambient values are not different from each other), then n(row) before and after should not be the same
cal_nfer <- function(df_f,df_cf,name_f,name_cf){
  df1 <- subset(df_f,exp==name_f)[,c("ambient","elevated")]
  df2 <- subset(df_cf,exp==name_cf)[,c("ambient","elevated")]
  df3 <- merge(df1,df2,by=c("ambient"),all.x=TRUE) 
  co2_a <- subset(df_cf,exp==name_cf)$co2_a[1]
  co2_e <- subset(df_cf,exp==name_cf)$co2_e[1]
  
  if (nrow(df3)-nrow(df2)==0){
    output_value <- mean(log(df3$elevated.y/df3$elevated.x),na.rm=TRUE)/log(co2_e/co2_a)
  } else {output_value <-NA}
  return(output_value)
  #-----------------------------------------------------------------------
  # Output: df_c_sub
  #-----------------------------------------------------------------------
}
#check which variable are missing, which are not
unique(subset(kevin_othervars_cf,exp=="euroface4_pa_cf")$response)  #checked 
unique(subset(kevin_othervars_cf,exp=="euroface4_pe_cf")$response)#checked 
unique(subset(kevin_othervars_cf,exp=="euroface4_pn_cf")$response) #checked 
unique(subset(kevin_othervars_cf,exp=="new_zealand_face_cf")$response) #checked 
unique(subset(kevin_othervars_cf,exp=="duke2_cf")$response) #checked 

#fill them
final5$vcmax[final5$exp=="euroface4_pa_cf"] <- cal_nfer(kevin2_f_vcmax,kevin2_cf_vcmax,"euroface4_pa_f","euroface4_pa_cf")
final5$vcmax[final5$exp=="euroface4_pe_cf"] <- cal_nfer(kevin2_f_vcmax,kevin2_cf_vcmax,"euroface4_pe_f","euroface4_pe_cf")
final5$vcmax[final5$exp=="euroface4_pn_cf"] <- cal_nfer(kevin2_f_vcmax,kevin2_cf_vcmax,"euroface4_pn_f","euroface4_pn_cf")
final5$vcmax[final5$exp=="new_zealand_face_cf"] <- cal_nfer(kevin2_f_vcmax,kevin2_cf_vcmax,"new_zealand_face_f","new_zealand_face_cf")

final5$jmax[final5$exp=="euroface4_pa_cf"] <- cal_nfer(kevin2_f_jmax,kevin2_cf_jmax,"euroface4_pa_f","euroface4_pa_cf")
final5$jmax[final5$exp=="euroface4_pe_cf"] <- cal_nfer(kevin2_f_jmax,kevin2_cf_jmax,"euroface4_pe_f","euroface4_pe_cf")
final5$jmax[final5$exp=="euroface4_pn_cf"] <- cal_nfer(kevin2_f_jmax,kevin2_cf_jmax,"euroface4_pn_f","euroface4_pn_cf")

final5$LMA[final5$exp=="euroface4_pn_cf"] <- cal_nfer(logr_f_LMA,logr_cf_LMA,"euroface4_pn_f","euroface4_pn_cf")
final5$LMA[final5$exp=="new_zealand_face_cf"] <- cal_nfer(logr_f_LMA,logr_cf_LMA,"new_zealand_face_f","new_zealand_face_cf")

final5$narea[final5$exp=="new_zealand_face_cf"] <- cal_nfer(logr_f_narea,logr_cf_narea,"new_zealand_face_f","new_zealand_face_cf")
final5$narea[final5$exp=="euroface4_pe_cf"] <-  cal_nfer(logr_f_narea,logr_cf_narea,"euroface4_pe_f","euroface4_pe_cf")
final5$narea[final5$exp=="euroface4_pn_cf"] <-  cal_nfer(logr_f_narea,logr_cf_narea,"euroface4_pn_f","euroface4_pn_cf")

final5$nmass[final5$exp=="new_zealand_face_cf"] <- cal_nfer(logr_f_nmass,logr_cf_nmass,"new_zealand_face_f","new_zealand_face_cf")
final5$nmass[final5$exp=="duke2_cf"] <- cal_nfer(logr_f_nmass,logr_cf_nmass,"duke2_f","duke2_cf")
final5$nmass[final5$exp=="euroface4_pe_cf"] <-  cal_nfer(logr_f_nmass,logr_cf_nmass,"euroface4_pe_f","euroface4_pe_cf")
final5$nmass[final5$exp=="euroface4_pn_cf"] <-  cal_nfer(logr_f_nmass,logr_cf_nmass,"euroface4_pn_f","euroface4_pn_cf")

final5$anpp[final5$exp=="euroface4_pa_cf"] <- cal_nfer(anpp_new_f,anpp_new_cf,"euroface4_pa_f","euroface4_pa_cf")
final5$anpp[final5$exp=="euroface4_pe_cf"] <-  cal_nfer(anpp_new_f,anpp_new_cf,"euroface4_pe_f","euroface4_pe_cf")
final5$anpp[final5$exp=="euroface4_pn_cf"] <-  cal_nfer(anpp_new_f,anpp_new_cf,"euroface4_pn_f","euroface4_pn_cf")

final5$lai[final5$exp=="euroface4_pa_cf"] <- log(subset(kevin_othervars_cf,exp=="euroface4_pa_cf" &response=="lai_max")$elevated/subset(kevin_othervars_cf,exp=="euroface4_pa_f" &response=="lai_max")$elevated)/(log(subset(kevin_othervars_cf,exp=="euroface4_pa_cf" &response=="lai_max")$co2_e/subset(kevin_othervars_cf,exp=="euroface4_pa_cf" &response=="lai_max")$co2_a))
final5$lai[final5$exp=="euroface4_pe_cf"] <- log(subset(kevin_othervars_cf,exp=="euroface4_pe_cf" &response=="lai_max")$elevated/subset(kevin_othervars_cf,exp=="euroface4_pe_f" &response=="lai_max")$elevated)/(log(subset(kevin_othervars_cf,exp=="euroface4_pe_cf" &response=="lai_max")$co2_e/subset(kevin_othervars_cf,exp=="euroface4_pe_cf" &response=="lai_max")$co2_a))
final5$lai[final5$exp=="euroface4_pn_cf"] <- log(subset(kevin_othervars_cf,exp=="euroface4_pn_cf" &response=="lai_max")$elevated/subset(kevin_othervars_cf,exp=="euroface4_pn_f" &response=="lai_max")$elevated)/(log(subset(kevin_othervars_cf,exp=="euroface4_pn_cf" &response=="lai_max")$co2_e/subset(kevin_othervars_cf,exp=="euroface4_pn_cf" &response=="lai_max")$co2_a))

#fill dukeface's anpp - something looked wrong here (why ambient all has same values. Not combined yet)
subset(kevin_othervars,exp=="duke2_cf"& response=="anpp")[,c("Unit","citation","ambient","elevated")]

#soil mineral N was added (see above)
final5$soil_mineral_N[final5$exp=="duke2_cf"] <- duke2_cf_soil

#fill euroface4_pe_cf's LMA
unique(subset(kevin_othervars_cf,exp=="euroface4_pe_cf"&response=="leaf_n")$Unit)
unique(subset(kevin_othervars_cf,exp=="euroface4_pe_f"&response=="leaf_n")$Unit)

nmass_cf <- mean(subset(kevin_othervars_cf,exp=="euroface4_pe_cf"&response=="leaf_n" & Unit=="g_kg")$elevated)
nmass_f <- mean(subset(kevin_othervars_cf,exp=="euroface4_pe_f"&response=="leaf_n" & Unit=="g_kg")$elevated)
narea_cf <- mean(subset(kevin_othervars_cf,exp=="euroface4_pe_cf"&response=="leaf_n" & Unit=="g_m2")$elevated)
narea_f <- mean(subset(kevin_othervars_cf,exp=="euroface4_pe_f"&response=="leaf_n" & Unit=="g_m2")$elevated)
co2_a <- mean(subset(kevin_othervars_cf,exp=="euroface4_pe_cf"&response=="leaf_n" & Unit=="g_kg")$co2_a)
co2_e <- mean(subset(kevin_othervars_cf,exp=="euroface4_pe_cf"&response=="leaf_n" & Unit=="g_kg")$co2_e)
final5$LMA[final5$exp=="euroface4_pe_cf"] <- log((narea_cf/nmass_cf)/(narea_f/nmass_f))/log(co2_e/co2_a)
 
#check if n fertilization at eCO2 is all site-species
unique(subset(kevin_othervars,exp=="duke2_cf")$dominant_species)
unique(subset(kevin_othervars,exp=="euroface4_pa_cf")$dominant_species)
unique(subset(kevin_othervars,exp=="euroface4_pe_cf")$dominant_species)
unique(subset(kevin_othervars,exp=="euroface4_pn_cf")$dominant_species)
unique(subset(kevin_othervars,exp=="new_zealand_face_cf")$dominant_species)
#new zealand includes multiple species because it is temperate grasslands: https://www.researchgate.net/publication/236201570_Photosynthetic_responses_of_temperate_species_to_free_air_CO2_enrichment_FACE_in_a_grazed_New_Zealand_pasture

#add ecm type
final5$ecm_type <- NA
final5$ecm_type[final5$exp=="grassotc"] <- "Nfix" # 2 N-fixing vs. 4 AM
final5$ecm_type[final5$exp=="biocon_c"]<- "Nfix"
final5$ecm_type[final5$exp=="swissface_trifolium2_c"]<- "Nfix"
final5$ecm_type[final5$exp=="new_zealand_face_c"] <- "Nfix"
final5$ecm_type[final5$exp=="chinaminiface"] <- "Nfix"
final5$ecm_type[final5$exp=="glycinece"] <- "Nfix"

#further process (gap-filling, correcting bugs)

#2. remove high Tleaf measurement's vcmax (actually, two papers: darbah_et_al_2010a and darbah_et_al_2010b)
final5$vcmax[final5$exp=="facts_ii_face3_pt_c"] <- agg_meta_sen_coef(response_ratio_v2(subset(kevin2_c_vcmax,exp=="facts_ii_face3_pt_c" & citation!="darbah_et_al_2010b" & citation!="darbah_et_al_2010a")),"vcmax")$vcmax
final5$vcmax[final5$exp=="facts_ii_face4_bp_c"] <- agg_meta_sen_coef(response_ratio_v2(subset(kevin2_c_vcmax,exp=="facts_ii_face4_bp_c" & citation!="darbah_et_al_2010b")),"vcmax")$vcmax

#add popface's soil 
old_data <- read_csv("~/data/gcme/data_received_190325/NewData_wide_CORRECTED2.csv") %>%
  mutate( ambient_Sd  = as.numeric(ambient_Sd),  ambient_Se  = as.numeric(ambient_Se), 
          elevated_Sd = as.numeric(elevated_Sd), elevated_Se = as.numeric(elevated_Se),
          co2_a  = as.numeric(co2_a),  co2_e  = as.numeric(co2_e), 
          ambient  = as.numeric(ambient),  elevated  = as.numeric(elevated))

old_data$exp <- tolower(old_data$exp_nam)
popface <-subset(old_data, (exp_nam=="POPFACE_pa"|exp_nam=="POPFACE_pe"|exp_nam=="POPFACE_pn")&Data_type=="soil_mineral_N")
new_popface <- as.data.frame(agg_meta_sen_coef(response_ratio_v2(popface),"soil_mineral_N")[,c("exp","soil_mineral_N")])

final5$soil_mineral_N[final5$exp=="popface_pa_c"] <- new_popface$soil_mineral_N[new_popface$exp=="popface_pa"]
final5$soil_mineral_N[final5$exp=="popface_pe_c"] <- new_popface$soil_mineral_N[new_popface$exp=="popface_pe"]
final5$soil_mineral_N[final5$exp=="popface_pn_c"] <- new_popface$soil_mineral_N[new_popface$exp=="popface_pn"]

#Data from Cesar's email
#Cesar_anpp <- read.csv("~/data/gcme/cesar/ANPP.csv")
#Cesar's anpp was not used at the end - because (1) it can only help to add 2 sites and (2) their site-name is really confused and not reliable to be merged with our dataset e.g.they define Aspen and FACTS II FACE3_pt_c as the same site but it is different in our database!
#they also define N-high, low, AM, ECM that we have not define yet. Best not to combine
#they also didn't define unit so it is difficult to know if it is culumative or actual flux values - not combine.

###euroface4_pn_c
check <- subset(kevin_othervars,exp=="euroface4_pn_c")%>% group_by(response,Unit)  %>% summarise(number = n())
subset(kevin_othervars,exp=="euroface4_pn_c" & response=="lai_max"& citation=="liberloo_et_al_2006")[,c("co2_a","co2_e")]
#LAI
subset(kevin_othervars,exp=="euroface4_pn_c" & response=="lai_max")[,c("ambient","elevated","citation")]
unique(subset(kevin_othervars,exp=="euroface4_pn_c" & response=="lai_max"& citation=="liberloo_et_al_2006")[,c("co2_a","co2_e")])
LAI_a <- mean(subset(kevin_othervars,exp=="euroface4_pn_c" & response=="lai_max"& citation=="liberloo_et_al_2006")$ambient,na.rm=TRUE)
LAI_e <- mean(subset(kevin_othervars,exp=="euroface4_pn_c" & response=="lai_max"& citation=="liberloo_et_al_2006")$elevated,na.rm=TRUE)
euroface4_pn_lai <- log((LAI_e/LAI_a))/log(550/368)

####duke_c
check <- subset(kevin_othervars,exp=="duke_c")%>% group_by(response,Unit)  %>% summarise(number = n())

#LMA obtained from leaf_c in area and mass basis
subset(subset(kevin_othervars,exp=="duke_c" & response=="leaf_c"& Unit=="g_100g")[,c("ambient","citation")])
subset(subset(kevin_othervars,exp=="duke_c" & response=="leaf_c"& Unit=="gc_m2")[,c("ambient","citation")])
#above, merge basing on the same cirtation
Cmass_a <- mean(subset(kevin_othervars,exp=="duke_c" & response=="leaf_c" & Unit=="g_100g" & citation=="hamilton_et_al_2004")$ambient)
Cmass_e <- mean(subset(kevin_othervars,exp=="duke_c" & response=="leaf_c" & Unit=="g_100g" & citation=="hamilton_et_al_2004")$elevated)

Carea_a <- mean(subset(kevin_othervars,exp=="duke_c" & response=="leaf_c" & Unit=="gc_m2" & citation=="hamilton_et_al_2002")$ambient)
Carea_e <- mean(subset(kevin_othervars,exp=="duke_c" & response=="leaf_c" & Unit=="gc_m2" & citation=="hamilton_et_al_2002")$elevated)
unique(subset(kevin_othervars,exp=="duke_c" & response=="leaf_c" & Unit=="gc_m2" & citation=="hamilton_et_al_2002")[,c("co2_a","co2_e")])
#LMA = Carea/Cmass
duke_c_lma <- log((Carea_e/Cmass_e)/(Carea_a/Cmass_a))/log(563/363)

########euroface4_pe_c
#LAI
check <- subset(kevin_othervars,exp=="euroface4_pe_c")%>% group_by(response,Unit)  %>% summarise(number = n())
subset(kevin_othervars,exp=="euroface4_pe_c" & response=="lai_max")[,c("co2_a","co2_e","ambient","citation")]
LAI_a <- mean(subset(kevin_othervars,exp=="euroface4_pe_c" & response=="lai_max"& citation=="liberloo_et_al_2006")$ambient,na.rm=TRUE)
LAI_e <- mean(subset(kevin_othervars,exp=="euroface4_pe_c" & response=="lai_max"& citation=="liberloo_et_al_2006")$elevated,na.rm=TRUE)
euroface4_pe_lai <- log((LAI_e/LAI_a))/log(550/368)
#LMA
nmass_a <- mean(subset(kevin_othervars,exp=="euroface4_pe_c" & response=="leaf_n"&Unit=="g_kg")$ambient,na.rm=TRUE)/1000 #convert to g/g
nmass_e <- mean(subset(kevin_othervars,exp=="euroface4_pe_c" & response=="leaf_n"&Unit=="g_kg")$elevated,na.rm=TRUE)/1000
narea_a <- mean(subset(kevin_othervars,exp=="euroface4_pe_c" & response=="leaf_n"&Unit=="g_m2")$ambient,na.rm=TRUE)
narea_e <- mean(subset(kevin_othervars,exp=="euroface4_pe_c" & response=="leaf_n"&Unit=="g_m2")$elevated,na.rm=TRUE)
#LMA = Narea/Nmass
euroface4_pe_lma <- log((narea_e/nmass_e)/(narea_a/nmass_a))/log(550/368)

######soyfacesoy2_c:
check <- subset(kevin_othervars,exp=="soyfacesoy2_c")%>% group_by(response,Unit)  %>% summarise(number = n())
#narea = lma * nmass /10 (below cm_/g is actually g/cm2)
narea_a <- mean(subset(logr_c_LMA,exp=="soyfacesoy2_c"&Unit=="cm_/g")$ambient)*mean(subset(kevin_othervars,exp=="soyfacesoy2_c"&response=="leaf_n" & Unit=="g_kg")$ambient)/10
narea_e <- mean(subset(logr_c_LMA,exp=="soyfacesoy2_c"&Unit=="cm_/g")$elevated)*mean(subset(kevin_othervars,exp=="soyfacesoy2_c"&response=="leaf_n" & Unit=="g_kg")$elevated)/10
unique(subset(logr_c_LMA,exp=="soyfacesoy2_c")[,c("co2_a","co2_e")])
soyfacesoy2_narea <- log((narea_e/narea_a))/log(548/373)

#######euroface4_pa_c
check <- subset(kevin_othervars,exp=="euroface4_pa_c")%>% group_by(response,Unit)  %>% summarise(number = n())

#lai
subset(kevin_othervars,exp=="euroface4_pa_c" & response=="lai_max"& citation=="liberloo_et_al_2006")[,c("co2_a","co2_e")]
LAI_a <- mean(subset(kevin_othervars,exp=="euroface4_pa_c" & response=="lai_max"& citation=="liberloo_et_al_2006")$ambient,na.rm=TRUE)
LAI_e <- mean(subset(kevin_othervars,exp=="euroface4_pa_c" & response=="lai_max"& citation=="liberloo_et_al_2006")$elevated,na.rm=TRUE)
euroface4_pa_lai <- log((LAI_e/LAI_a))/log(550/368)

#narea = nmass * LMA
unique(subset(kevin_othervars,exp=="soyfacetobacco9_c")[,c("co2_a","co2_e")])
soyfacetobacco9_narea <- final_mean$nmass[final_mean$exp=="soyfacetobacco9_c"] + final_mean$LMA[final_mean$exp=="soyfacetobacco9_c"]

#######giface 
#Nmass and Narea were filled by below publication
#Carbon dioxide fertilisation and supressed respiration induce enhanced spring biomass production in a mixed species temperate meadow exposed to moderate carbon dioxide enrichment
#https://pubmed.ncbi.nlm.nih.gov/32480439/
# in table 5, see Nmass (mg/g), Cmass was also recorded but not put here.
#co2 changes from 400 to 480 (see method: the study site)
Nmass_a <- (34+41.2+26.4+36.3+31.8+43.5)/6/1000 # averaged them, then converted from mg/g to g/g
Nmass_e <- (32.7+40.7+24.2+33.3+31.4+42.4)/6/1000 # in mg/g 
giface_c_nmass <- log((Nmass_e/Nmass_a))/log(480/400)

#fill biforface_c's lma
aa <- subset(kevin_othervars,exp=="biforface_c")%>% group_by(response,Unit)  %>% summarise(number = n())
nmass_a <- mean(subset(kevin_othervars,exp=="biforface_c" & response=="leaf_n" & Unit=="mg_g")$ambient)
narea_a <- mean(subset(kevin_othervars,exp=="biforface_c" & response=="leaf_n" & Unit=="g_m2")$ambient)
nmass_e <- mean(subset(kevin_othervars,exp=="biforface_c" & response=="leaf_n" & Unit=="mg_g")$elevated)
narea_e <- mean(subset(kevin_othervars,exp=="biforface_c" & response=="leaf_n" & Unit=="g_m2")$elevated)

lma_e <- narea_e/nmass_e
lma_a <- narea_a/nmass_a
unique(subset(kevin_othervars,exp=="biforface_c")[,c("co2_a","co2_e")])
biforface_c_lma <- log(lma_e/lma_a)/log(558/408)

#gap-fill duke2_cf
#vcmax - duke2 and duke cannot be merged since they duke2 is used for non-fertilization vs. fertilization.
#duke2 face description: Since 2005, the rings were split into quadrants and N fertilization (11.2 g of N m−2yr−1 as ammonium nitrate) was applied to half of them. 
#From: https://www.sciencedirect.com/science/article/pii/S0168192316303823
#by looking at measurement time - duke2 has been measured for >10 years (which means, from 2005, starts with N fertilization), while duke has been measured just for 3 years (from 1996-1999, without any fertilization). So it should be separate
subset(kevin2_f_vcmax,exp=="duke2_f")[,c("ambient","elevated","logr","citation")];subset(kevin2_cf_vcmax,exp=="duke2_cf")[,c("ambient","elevated","logr","citation")];
subset(kevin2_c_vcmax,exp=="duke_c")[,c("ambient","elevated","co2_a","co2_e","logr","citation")]
head(subset(kevin2_c_vcmax,exp=="duke_c"));head(subset(kevin2_c_vcmax,exp=="duke2_c"))

#because this cannot be merged by individuals-level. We can only calculate them at site-level
#that is, sensitivity value = mean(duke2_cf) - mean(duke2_f)
duke2_cf_vcmax <-(mean(kevin2_cf_vcmax$logr[kevin2_cf_vcmax$exp=="duke2_cf"]) - mean(kevin2_f_vcmax$logr[kevin2_f_vcmax$exp=="duke2_f"]))/
  log(kevin2_cf_vcmax$co2_e[kevin2_cf_vcmax$exp=="duke2_cf"][1]/kevin2_cf_vcmax$co2_a[kevin2_cf_vcmax$exp=="duke2_cf"][1])

#jmax
subset(kevin2_f_jmax,exp=="duke2_f")[,c("ambient","elevated","logr","citation")];subset(kevin2_cf_jmax,exp=="duke2_cf")[,c("ambient","elevated","logr","citation")];
subset(kevin2_c_jmax,exp=="duke_c")[,c("ambient","elevated","co2_a","co2_e","logr","citation")]

a1 <- subset(kevin2_f_jmax,exp=="duke2_f")[,c("ambient","elevated","logr","citation")]
#remove the one with 4.2, which cannot be parallel to cf on below, and the values looked wrong
a1$logr[a1$citation=="oishi_et_al_2014"] <- NA

a2 <- subset(kevin2_cf_jmax,exp=="duke2_cf")[,c("ambient","elevated","logr","citation")]
a1;a2
#recalculate logr = logcf - cf
duke2_cf_jmax <- (mean(a2$logr) - mean(a1$logr,na.rm=TRUE))/log(563/363)

#  LMA and narea were corrected by maier_et_al_2008 -https://academic.oup.com/treephys/article/28/4/597/1720624
#Short-term effects of fertilization on photosynthesis and leaf morphology of field-grown loblolly pine following long-term exposure to elevated CO2 concentration
#corrected by their table 2 (using 1-year old data)
unique(subset(logr_c_LMA,exp=="duke2_c")[,c("co2_a","co2_e")])
duke2_cf_lma <- log(113.3/100.6)/log(563/363)

duke2_cf_narea <- log(1.34/1.44)/log(563/363)

subset(logr_f_lai,exp=="duke2_f")[,c("ambient","elevated","logr","citation")]
subset(logr_cf_lai,exp=="duke2_cf")[,c("ambient","elevated","logr","citation")]
# oishi_et_al_2014 - lai too low - below has problem - ambient should be 3.8 not 38!
#calculate it from domec_et_al_2012 directly 
duke2_cf_lai <- log(4.3/4.05)/log(563/363)

#include all!
final5$lai[final5$exp=="euroface4_pa_c"] <- euroface4_pa_lai
final5$lai[final5$exp=="euroface4_pn_c"] <- euroface4_pn_lai
final5$lai[final5$exp=="euroface4_pe_c"] <- euroface4_pe_lai

final5$LMA[final5$exp=="duke_c"] <- duke_c_lma
final5$LMA[final5$exp=="euroface4_pe_c"] <- euroface4_pe_lma
final5$LMA[final5$exp=="biforface_c"] <- biforface_c_lma
final5$narea[final5$exp=="soyfacesoy2_c"] <- soyfacesoy2_narea
final5$narea[final5$exp=="soyfacetobacco9_c"] <- soyfacetobacco9_narea
final5$nmass[final5$exp=="giface_c"] <- giface_c_nmass

#gap-fill duke2_cf
final5$jmax[final5$exp=="duke2_cf"] <- duke2_cf_jmax
final5$jmax_vcmax[final5$exp=="duke2_cf"] <- final5$jmax[final5$exp=="duke2_cf"] -final5$vcmax[final5$exp=="duke2_cf"] 
final5$LMA[final5$exp=="duke2_cf"] <- duke2_cf_lma
final5$narea[final5$exp=="duke2_cf"] <- duke2_cf_narea
final5$lai[final5$exp=="duke2_cf"] <- duke2_cf_lai
final5$vcmax[final5$exp=="duke2_cf"] <- duke2_cf_vcmax

#make type_name more clear
final5$type_name[final5$condition=="highN"] <- "highN"
final5$type_name[final5$condition=="lowN"] <- "lowN"
final5$type_name[final5$condition=="Fertilization"] <- "Fertilization"

#check exp name
final5$exp[final5$condition=="Fertilization"]

final5$type_name[final5$exp=="duke2_c"] <- "No_fertilization"
final5$type_name[final5$exp=="euroface4_pa_c"] <- "No_fertilization"
final5$type_name[final5$exp=="euroface4_pe_c"] <- "No_fertilization"
final5$type_name[final5$exp=="euroface4_pn_c"] <- "No_fertilization"
final5$type_name[final5$exp=="new_zealand_face_c"] <- "No_fertilization"

#re-name ecosystem name 
final5$ecosystem[final5$ecosystem=="forest"] <- "Forest"
final5$ecosystem[final5$ecosystem=="grassland"] <- "Grassland"
final5$ecosystem[final5$ecosystem=="cropland"] <- "Cropland"


#edit comments
final5$comments[final5$condition=="light"] <- "vcmax and jmax are sensitivity coefficients = log(vcmax-ele/vcmax-amb)"
final5$jmax_vcmax <- final5$jmax-final5$vcmax

csvfile <- paste("~/data/gcme/MS_data/plot_data.csv")
write.csv(final5, csvfile, row.names = TRUE)

#check japan_face
#https://academic.oup.com/pcp/article/55/2/381/1862434#85129527
#ko: Koshihikari
#Ta: Takanari
#Two rice (Oryza sativa L.) cultivars were used in this study. ‘Koshihikari’ is a japonica variety, and ‘Takanari’ is an indica variety.

#https://academic.oup.com/pcp/article/55/2/370/1861701?login=true#85129192：
#L: another place
#japan should be combined

#calculate site-info
final5 <- read.csv("~/data/gcme/MS_data/plot_data.csv")
final5$type_name[is.na(final5$type_name)==T]<-"others"
dim(subset(final5,condition=="co2"))
dim(subset(final5,condition=="co2" &(type_name)=="others")) # vcmax numbers

dim(subset(final5,is.na(nmass)==F))
dim(subset(final5,is.na(narea)==F))
dim(subset(final5,is.na(LMA)==F))
dim(subset(final5,is.na(anpp)==F))
dim(subset(final5,is.na(bnpp)==F))
dim(subset(final5,is.na(lai)==F))
dim(subset(final5,is.na(root_shoot_ratio)==F))
dim(subset(final5,is.na(soil_mineral_N)==F))
