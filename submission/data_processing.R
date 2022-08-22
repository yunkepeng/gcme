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
#logr_var  calculated from ROM
response_ratio_v2 <- function(df){
  df$logr <- log(df$elevated/df$ambient)
  df_c_sub <- df
  return(df_c_sub)
  #-----------------------------------------------------------------------
  # Output: df_c_sub
  #-----------------------------------------------------------------------
}

#agg_meta_sen_coef: calculate sensitivity coefficient from individuals dataset to site-mean dataset
agg_meta_sen_coef <- function(df){
  explist <- unique(df$exp)
  mylist <- list() #create an empty list
  #convert to sensitivity coef
  df$sen_coef <- df$logr/log(df$co2_e/df$co2_a)
  
  #all using mean value
  for (i in 1:length(explist)){
    mean_value <- mean(subset(df,exp==explist[i])$sen_coef,na.rm=TRUE)
    n_samples <- length(subset(df,exp==explist[i])$logr)
    df_box <- tibble(
      exp=explist[i],
      middle = mean_value,
      no=n_samples)
    mylist[[i]] <- df_box}
  output <- do.call("rbind",mylist)
  return(output)
}

#agg_meta: calculate response ratio from individuals dataset to site-mean dataset
agg_meta <- function(df){
  explist <- unique(df$exp)
  mylist <- list() 
  for (i in 1:length(explist)){
    mean_value <- mean(subset(df,exp==explist[i])$logr,na.rm=TRUE)
    n_samples <- length(subset(df,exp==explist[i])$logr)
    df_box <- tibble(
      exp=explist[i], 
      middle = mean_value,
      no=n_samples)
    mylist[[i]] <- df_box}
  output <- do.call("rbind",mylist)
  return(output)
}

#combine co2 and co2+nfertilization site-mean dataset
combine_co2_cf <- function(logr_c_var,logr_f_var,logr_cf_var,name){
  all_logr_c_lma <- agg_meta_sen_coef(logr_c_var)[,c("exp","middle")]
  all_logr_c_lma$condition <- "co2"
  #co2 + fer / fer
  all_logr_f_lma <- agg_meta(logr_f_var)[,c("exp","middle")]
  all_logr_cf_lma <- agg_meta(logr_cf_var)[,c("exp","middle")]
  all_logr_cf_lma$exp_old <- all_logr_cf_lma$exp 
  all_logr_cf_lma$exp <-  paste(substr(all_logr_cf_lma$exp,1,nchar(all_logr_cf_lma$exp)-2),"f",sep="")
  
  all_logr_f_lma2 <- all_logr_f_lma %>% filter(exp %in%all_logr_cf_lma$exp)
  new_df <- merge(all_logr_f_lma2,all_logr_cf_lma,by=c("exp"),all.x=TRUE)
  for (i in 1:nrow(new_df)) {
    new_df$co2[i] <- subset(logr_cf_var,exp==new_df$exp_old[i])$co2_e[1]/
      subset(logr_cf_var,exp==new_df$exp_old[i])$co2_a[1]
  } 
  new_df$middle <- (new_df$middle.y - new_df$middle.x)/log(new_df$co2)
  new_df$exp <- new_df$exp_old
  new_df$condition <- "Fertilization"
  
  lma_plot <- dplyr::bind_rows(all_logr_c_lma,new_df[,c("exp","middle","condition")]) 
  names(lma_plot) <- c("exp",name,"condition")
  return(lma_plot)
}

combine_co2_c <- function(logr_c_var,logr_f_var,logr_cf_var,name){
  all_logr_c_lma <- agg_meta_sen_coef(logr_c_var)[,c("exp","middle")]
  all_logr_c_lma$condition <- "co2"
  lma_plot <- all_logr_c_lma
  names(lma_plot) <- c("exp",name,"condition")
  return(lma_plot)
}

#####3. CO2 effect data
#read Kevin
kevin <- read.csv("~/data/gcme/kevin/orig_vcmax/JunkePeng_11252021.csv")

#1. correct exp_nam to make it consistent with GCME
kevin$site[is.na(kevin$site)==TRUE] <- "riceface_japan_a_2003_3938_14057"

kevin$exp_nam <- kevin$site 

kevin <- rename(kevin, c(ambient = x_c, elevated=x_t, ambient_Sd=sd_c, elevated_Sd=sd_t,ambient_Se=se_c,elevated_Se=se_t,n_plots=rep_c,
                         z=elevation, co2_a=c_c, co2_e=c_t, nfertQ_a = n_c, nfertQ_e = n_t, pfertQ_a = p_c, pfertQ_e = p_t,kfertQ_a = k_c, kfertQ_e = k_t,
                         warmQ_e1 = w_t1, warmQ_e2 = w_t2, warmQ_e3 = w_t3, Unit=x_units))

#adjust temperature response
kevin$warmQ_e2[is.na(kevin$warmQ_e2)==TRUE] <- 0

#correct a few sampling year --> when looking at org csv
kevin$sampling_year[is.na(kevin$sampling_year)==TRUE & kevin$site=="brandbjerg"] <- 2011
kevin$sampling_year[is.na(kevin$sampling_year)==TRUE & kevin$site=="popface"] <- 2002
kevin$sampling_year[is.na(kevin$sampling_year)==TRUE & kevin$site=="biocon"] <- 2005
kevin$sampling_year[kevin$sampling_year=="2005-2010"] <- 2008
kevin$sampling_year[kevin$sampling_year=="1996-2010"] <- 2003
kevin$sampling_year[kevin$sampling_year=="2003-2006"] <- 2005
kevin$sampling_year[kevin$citation=="domec_et_al_2012"] <- 1997 # by looking at their info it says +1y. then we assume it is 1996+1
kevin$sampling_year[kevin$citation=="ellsworth_et_al_2012"] <- 1997
kevin$start_year[is.na(kevin$start_year)==TRUE] <- 1992

kevin$sampling_year <- as.numeric(kevin$sampling_year)
kevin$start_year <- as.numeric(kevin$start_year)
kevin$Year <-  kevin$sampling_year - kevin$start_year
summary(kevin$Year)

#correct elevation
aaa <- aggregate(kevin,by=list(kevin$lon,kevin$lat), FUN=mean, na.rm=TRUE)[,c("lon","lat")]
aaa$sitename <- paste("c",1:length(aaa$lon),sep="")
devtools::load_all("~/yunkepeng/gcme/pmodel/ingestr/")
df_etopo <- ingest(aaa,source = "etopo1",dir = "~/data/etopo/" )
aaa$elv <- as.numeric(as.data.frame(df_etopo$data))
kevin_z <- merge(kevin,aaa[,c("lon","lat","elv")],by=c("lon","lat"),all.x=TRUE)
plot(kevin_z$z~kevin_z$elv) # looks ok - now interploate original elevation value with etopo elevation
kevin_z$z[is.na(kevin_z$z)==TRUE] <- kevin_z$elv[is.na(kevin_z$z)==TRUE] 
kevin_z <- kevin_z[, !(colnames(kevin_z) %in% c("elv"))]

summary(kevin_z$Year)
summary(kevin_z$start_year)
summary(kevin_z$sampling_year) #all from 1992 to 2016 - good news!

kevin_z$year_start <- kevin_z$start_year
kevin_z$year_end <- kevin_z$sampling_year


kevin2_final <- response_ratio_v2(kevin_z)

#show pft
kevin2_final$ecosystem[kevin2_final$ecosystem=="temperate_forest"] <- "forest"
kevin2_final$ecosystem[kevin2_final$ecosystem=="heathland"] <- "grassland"
kevin2_final$ecosystem[kevin2_final$ecosystem=="shrubland"] <- "forest"

kevin2_c_vcmax <- subset(kevin2_final, treatment=="c" & response =="vcmax")
kevin2_f_vcmax <- subset(kevin2_final, treatment=="f" & response =="vcmax")
kevin2_cf_vcmax <- subset(kevin2_final, treatment=="cf" & response =="vcmax")

kevin2_c_jmax <- subset(kevin2_final, treatment=="c" & response =="jmax")
kevin2_f_jmax <- subset(kevin2_final, treatment=="f" & response =="jmax")
kevin2_cf_jmax <- subset(kevin2_final, treatment=="cf" & response =="jmax")

kevin_vcmax_plotmean <- agg_meta_sen_coef(kevin2_c_vcmax)[,c("exp","middle","no")]
names(kevin_vcmax_plotmean) <- c("exp","vcmax","no")

kevin_jmax_plotmean <- agg_meta_sen_coef(kevin2_c_jmax)[,c("exp","middle","no")]
names(kevin_jmax_plotmean) <- c("exp","jmax","no")

kevin_ecosystem <- as.data.frame(kevin2_c_vcmax %>% group_by(exp,ecosystem) %>% summarise(number=n()))

kevin_vj <- merge(kevin_vcmax_plotmean[,c("exp","vcmax")],kevin_jmax_plotmean[,c("exp","jmax")],by=c("exp"),all.x=TRUE)

kevin_vj_ecosystem <- merge(kevin_vj,kevin_ecosystem[,c("exp","ecosystem")],by=c("exp"),all.x=TRUE)


#Smith's data
smith_co2 <- read.csv("~/data/smith_keenan_gcb/gcb_co2/co2_data.csv")
smith_co2$exp_nam <- smith_co2$SiteID
smith_co2 <- subset(smith_co2,is.na(lat)==FALSE) # remove empty column
summary(smith_co2)
smith_final <- smith_co2%>% group_by(exp_nam,Species)  %>% summarise(number = n())

#first 2 columns are the same treatment in same plot (amb + elv) -check their Tleaf --> all equal
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

smith_co2_final <- subset(smith_co2,is.na(logr_v)==FALSE)

vcmax25_warmingco2_siteinfo <- rename(smith_co2_final, c(ambient = ambient_vcmax, elevated=Vcmax, ambient_Se=ambient_vcmaxSE, elevated_Se=SE,
                                                         co2_a=ambient_co2, co2_e=Ca,
                                                         ambient_j = ambient_jmax, elevated_j =Jmax, ambient_Se_j=ambient_jmaxSE, elevated_Se_j=Jmax_SE))


vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$Stature=="Tree"] <- "forest"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="ChinaRiceFACE"] <- "cropland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="UAMACFACE"] <- "cropland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="GlycineCE"] <- "cropland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="GossypiumGC"] <- "cropland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="MississippiCE"] <- "cropland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="WheatCE"] <- "cropland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="ChinaMiniFACE"] <- "cropland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="GrassOTC"] <- "grassland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="PHACE"] <- "grassland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="FescueCE"] <- "grassland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="NZFACE"] <- "grassland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="NevadaFACE"] <- "grassland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="BioCON"] <- "grassland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="PhalarisGC"] <- "grassland"
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
smith_all_plotmean$exp[smith_all_plotmean$exp=="AspenFACE"] <- "rhine-aspenface_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="BilyKriz"] <- "bily_kriz_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="Headley"] <- "headley_qp_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="Viesalm"] <- "vielsalm_c"

#check if smith and GCME have consistent site (but not consistent species!)
kevin_othervars <- read.csv("~/data/gcme/kevin_20220222/MESI_2022.csv")

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
unique(subset(kevin_othervars,exp=="rhine-aspenface_c")[,c("lon","lat","dominant_species")])

unique(subset(smith_all_simple,exp_nam=="BilyKriz")[,c("lon","lat","Genus","Species")])
unique(subset(kevin_othervars,exp=="bily_kriz_c")[,c("lon","lat","dominant_species")])

unique(subset(smith_all_simple,exp_nam=="Headley")[,c("lon","lat","Genus","Species")])
unique(subset(kevin_othervars,exp=="headley_qp_c")[,c("lon","lat","dominant_species")])

unique(subset(smith_all_simple,exp_nam=="Viesalm")[,c("lon","lat","Genus","Species")])
unique(subset(kevin_othervars,exp=="vielsalm_c")[,c("lon","lat","dominant_species")])

smith_all_plotmean$exp <- tolower(smith_all_plotmean$exp)

smith_all_plotmean <- smith_all_plotmean[,c("exp","vcmax","jmax","ecosystem")]

obs_co2 <- rbind(kevin_vj_ecosystem,smith_all_plotmean)

#combine_co2 under fertilization (c, cf, f)
Nfer_vcmax <- na.omit(subset(combine_co2_cf(kevin2_c_vcmax,kevin2_f_vcmax,kevin2_cf_vcmax,"middle"),condition=="Fertilization"))
Nfer_vcmax$ecosystem <- "forest"
Nfer_vcmax$ecosystem[Nfer_vcmax$exp=="new_zealand_face_cf"] <- "grassland"

Nfer_jmax <- na.omit(subset(combine_co2_cf(kevin2_c_jmax,kevin2_f_jmax,kevin2_cf_jmax,"middle"),condition=="Fertilization"))
Nfer_jmax$ecosystem <- "forest"

Nfer_kevin <- merge(Nfer_vcmax,Nfer_jmax,by=c("exp","condition","ecosystem"),all.x=TRUE)
names(Nfer_kevin) <- c("exp","condition","ecosystem","vcmax","jmax")

#N fertilization data from Walker 
walker <- read.csv("~/data/leaf_traits/Walker/LEAF_PHOTOSYNTHESIS_TRAITS_1224/data/Leaf_Photosynthesis_Traits.csv")

Nfer_vcmax <- subset(walker,Treatment_N!="not applicable" &Treatment_CO2!="not applicable" & Vcmax>0) %>% 
  group_by(lon,lat,species,Treatment_N,Treatment_P,Treatment_CO2,Treatment_light)  %>% summarise(mean = mean(Vcmax))
#keep all other treatments consistent
highN_vcmax <- merge(subset(Nfer_vcmax,Treatment_N=="high"&Treatment_CO2=="amb"),
                     subset(Nfer_vcmax,Treatment_N=="high"&Treatment_CO2=="ele"),
                     by=c("lon","lat","species","Treatment_N","Treatment_P","Treatment_light"),all.x=TRUE)
highN_vcmax <- subset(highN_vcmax,Treatment_light!="shade")

highN_vcmax #now, except co2 condition, all others are consistent

highN_vcmax$middle <- log(highN_vcmax$mean.y/highN_vcmax$mean.x)

highN_vcmax$site_species <- paste(highN_vcmax$lon,highN_vcmax$lat,highN_vcmax$species,sep="_")
highN_vcmax[,c("site_species","middle")]

lowN_vcmax <- merge(subset(Nfer_vcmax,Treatment_N=="low"&Treatment_CO2=="amb"),
                    subset(Nfer_vcmax,Treatment_N=="low"&Treatment_CO2=="ele"),
                    by=c("lon","lat","species","Treatment_N","Treatment_P","Treatment_light"),all.x=TRUE)

lowN_vcmax <- subset(lowN_vcmax,Treatment_light!="shade")
lowN_vcmax

lowN_vcmax$middle <- log(lowN_vcmax$mean.y/lowN_vcmax$mean.x)
lowN_vcmax$site_species <- paste(lowN_vcmax$lon,lowN_vcmax$lat,lowN_vcmax$species,sep="_")

#by checking them - most sites do have both high N and P comparing with low nutrient sites
low_high_vcmax <- merge(highN_vcmax[,c("lon","lat","site_species","middle")],na.omit(lowN_vcmax[,c("lon","lat","site_species","middle")]),
                        by=c("lon","lat","site_species"),all.x=TRUE)
names(low_high_vcmax) <- c("lon","lat","site_species","vcmax_highN","vcmax_lowN")
low_high_vcmax$co2_change <- NA

#input co2 change
low_high_vcmax$co2_change[low_high_vcmax$lon==-84.04] <- log(560/360) 
#https://pubmed.ncbi.nlm.nih.gov/11960756/
low_high_vcmax$co2_change[low_high_vcmax$lon==20.00] <- log(700/350) 
#https://www.researchgate.net/publication/242459413_Nutrient_and_genotypic_effects_on_CO2-responsiveness_Photosynthetic_regulation_in_Leucadendron_species_of_a_nutrient-poor_environment
low_high_vcmax$co2_change[low_high_vcmax$lon==89.14] <- log(720/360) 

low_high_vcmax$vcmax_highN_coef <- low_high_vcmax$vcmax_highN/low_high_vcmax$co2_change
low_high_vcmax$vcmax_lowN_coef <- low_high_vcmax$vcmax_lowN/low_high_vcmax$co2_change

#do the same for jmax

Nfer_jmax <- subset(walker,Treatment_N!="not applicable" &Treatment_CO2!="not applicable" & Jmax>0) %>% 
  group_by(lon,lat,species,Treatment_N,Treatment_P,Treatment_CO2,Treatment_light) %>% summarise(mean = mean(Jmax))

highN_jmax <- merge(subset(Nfer_jmax,Treatment_N=="high"&Treatment_CO2=="amb"),
                    subset(Nfer_jmax,Treatment_N=="high"&Treatment_CO2=="ele"),
                    by=c("lon","lat","species","Treatment_N","Treatment_P","Treatment_light"),all.x=TRUE)

highN_jmax <- subset(highN_jmax,Treatment_light!="shade")

highN_jmax$middle <- log(highN_jmax$mean.y/highN_jmax$mean.x)

highN_jmax$type_name <- "CO2 response at high N"
highN_jmax$site_species <- paste(highN_jmax$lon,highN_jmax$lat,highN_jmax$species,sep="_")

lowN_jmax <- merge(subset(Nfer_jmax,Treatment_N=="low"&Treatment_CO2=="amb"),
                   subset(Nfer_jmax,Treatment_N=="low"&Treatment_CO2=="ele"),
                   by=c("lon","lat","species","Treatment_N","Treatment_P","Treatment_light"),all.x=TRUE)

lowN_jmax <- subset(lowN_jmax,Treatment_light!="shade")

lowN_jmax$middle <- log(lowN_jmax$mean.y/lowN_jmax$mean.x)

lowN_jmax$type_name <- "CO2 response at low N"
lowN_jmax$site_species <- paste(lowN_jmax$lon,lowN_jmax$lat,lowN_jmax$species,sep="_")

low_high_jmax <- merge(highN_jmax[,c("lon","lat","site_species","middle")],na.omit(lowN_jmax[,c("lon","lat","site_species","middle")]),
                       by=c("lon","lat","site_species"),all.x=TRUE)
names(low_high_jmax) <- c("lon","lat","site_species","jmax_highN","jmax_lowN")
low_high_jmax$co2_change <- NA

#input co2 change
low_high_jmax$co2_change[low_high_jmax$lon==-84.04] <- log(560/360) 
# https://watermark.silverchair.com/22-5-321.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAv4wggL6BgkqhkiG9w0BBwagggLrMIIC5wIBADCCAuAGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMlQzU_4AV0WQkXxioAgEQgIICsZUQMf48_Ogae5yJOUenILkeX5E_nHvKWARnNsbV4Cd5VZ38yYrQBL5duSx1kSq97qPt5PcU5fhbJvGTYY3HdCO1nm4_DwncLNNjvauVf6SSlXs6ilSIzkeclSN2eRY7uzwXfa-JaJg8RHzLq8scOxja5YhaUCRogGCZDpTXryJ97Z1jobEOPG-YCdDI0apFWymgoaqcE286a1lG6-aQK_7kbxwinsZSf1_CcxhItSBy7LTZRxwANCisRFwFr0IXbqpDnoQJH8cBepblQkkZzcEfwclOsghsN26D0uNF4K4UTTpsHUd2jifg0ga_4fzNsP_iMmht62Y0HwRC7uSc1gD587zvIQGQoMa4yUBpPeHYXOLmcGNpURkWQOA-qfRnjAC_7PS8UUMQbSQp7LPD2ZpuAKxNMONSXTT4SwqlQYhxFyHQZyaQYEoaKe-XBSaVjdt2oA1IXfz7OZ72-KsJweC-buZHrjy64iCPtaJ7R0KsAuLtbO3Wsgu5DZyr3oddoZ30fuk2NVogC-7iHPvFfYxu0uqkRwQavdcIYvxVZvt3NkI534Wce9fIgQzoWr2PSbqJoenwM3Bs0QpemX2ByyxFWtCt-7qOcRXu71Iupj5Qt0S6DkcIqlEwEkruKStMO5MUv6CEtpeUaY2VyZ6kKnlh8ffsMkFgKFdARaaRUKc4Igp2JY-ptvQrLInZKoZlebgjHiiwWVQvDB4TqhmXi33pJ_iiriuNTd2_8AKChx7b777LNYKDri48i4Lny425sdVdxrlY_t8123c3OnnVTgmpTD8lUTRKh85iAnKcrts6hSFFeGrC9nM8zuFh0vpA2p8gq0wUqO5knbwNrqIcW_ZXV8v0cTmyNJjnKK7PvwKGq7HFBYfHLk1ZBUJNy7duVXwQ5nqvddVLFmZff_Dv2_-1
low_high_jmax$co2_change[low_high_jmax$lon==20.00] <- log(700/350) 
#https://www.researchgate.net/publication/242459413_Nutrient_and_genotypic_effects_on_CO2-responsiveness_Photosynthetic_regulation_in_Leucadendron_species_of_a_nutrient-poor_environment
low_high_jmax$co2_change[low_high_jmax$lon==89.14] <- log(720/360) 

low_high_jmax$jmax_highN_coef <- low_high_jmax$jmax_highN/low_high_jmax$co2_change
low_high_jmax$jmax_lowN_coef <- low_high_jmax$jmax_lowN/low_high_jmax$co2_change

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

#For "Black Spruce, ON, Canada"
unique(subset(df1,Dataset=="Black Spruce, ON, Canada")$Temp_Treatment) 
df3$Temp_Treatment[df3$Dataset=="Black Spruce, ON, Canada" & df3$Temp_Treatment=="ambient"] <- 22
df3$Temp_Treatment[df3$Dataset=="Black Spruce, ON, Canada" & df3$Temp_Treatment=="elevated"] <- 30

#For "Corymbia calophylla provs, AU-NSW"
unique(subset(df1,Dataset=="Corymbia calophylla provs, AU-NSW")$Temp_Treatment) 
df3$Temp_Treatment[df3$Dataset=="Corymbia calophylla provs, AU-NSW" & df3$Temp_Treatment=="ambient"] <- 26
df3$Temp_Treatment[df3$Dataset=="Corymbia calophylla provs, AU-NSW" & df3$Temp_Treatment=="elevated"] <- 32

#For "Eucalyptus tereticornis provs AU-NSW"
#needs to check ID_curve number spefically, since it has multiple temp control
id_one <- subset(df3,Dataset=="Eucalyptus tereticornis provs AU-NSW")$Curve_Id
unique((subset(df1,Dataset=="Eucalyptus tereticornis provs AU-NSW")  %>% filter(Curve_Id %in% id_one))$Temp_Treatment)
df3$Temp_Treatment[df3$Dataset=="Eucalyptus tereticornis provs AU-NSW" & df3$Temp_Treatment=="1"] <- 18
df3$Temp_Treatment[df3$Dataset=="Eucalyptus tereticornis provs AU-NSW" & df3$Temp_Treatment=="4"] <- 28.5

#For "Eucalyptus globulus, AU-NSW"
#By looking at https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.12314 
# Fig.1 (mean ambient T in whole year = 20) and description: a factorial between ambient temperature (aT) and ambient +3 °C (eT) by ambient CO2 concentration (aC)
#we can assume Tair changed from 20 to 23
df3$Temp_Treatment[df3$Dataset=="Eucalyptus globulus, AU-NSW" & df3$Temp_Treatment=="ambient"] <- 20
df3$Temp_Treatment[df3$Dataset=="Eucalyptus globulus, AU-NSW" & df3$Temp_Treatment=="elevated"] <- 23

#for "Eucalyptus parramattensis, AU-NSW"
#see below Table s1: we can assume ambient = 15, elevated = 18
#https://nph.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fnph.15668&file=nph15668-sup-0001-SupInfo.pdf
df3$Temp_Treatment[df3$Dataset=="Eucalyptus parramattensis, AU-NSW" & df3$Temp_Treatment=="ambient"] <- 15
df3$Temp_Treatment[df3$Dataset=="Eucalyptus parramattensis, AU-NSW" & df3$Temp_Treatment=="elevated"] <- 18

#for Eucalyptus tereticornis, AU-NSW
#see below Table s1: we can assume ambient = 14.4, elevated = 22.6
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

#####2. Input light effect data
walker <- read.csv("~/data/leaf_traits/Walker/LEAF_PHOTOSYNTHESIS_TRAITS_1224/data/Leaf_Photosynthesis_Traits.csv")

light_vcmax <- subset(walker,Treatment_light!="not applicable" & Vcmax>0) %>% 
  group_by(lon,lat,species,Treatment_N,Treatment_P,Treatment_CO2,Treatment_light)  %>% summarise(mean = mean(Vcmax))

#divided into two categories:
#(1) shade vs. sun
#(2) high vs. low light

light1 <- subset(light_vcmax,Treatment_light=="shade"|Treatment_light=="sun")
light1a <- merge(subset(light1,Treatment_light=="shade"),subset(light1,Treatment_light=="sun"),
                 by=c("lon","lat","species","Treatment_N","Treatment_P","Treatment_CO2"),all.x=TRUE)

#remove ele co2 and high-N
light1a <- subset(light1a,Treatment_CO2 != "ele" & Treatment_N!="high" & Treatment_N!="fertilised")
light1a

light1a$middle <- log(light1a$mean.y/light1a$mean.x)
light1a$type_name <- "shade_to_sun"

light2 <- subset(light_vcmax,Treatment_light=="high"|Treatment_light=="low")
light2a <- merge(subset(light2,Treatment_light=="low"),subset(light2,Treatment_light=="high"),
                 by=c("lon","lat","species","Treatment_N","Treatment_P","Treatment_CO2"),all.x=TRUE)
light2a$middle <- log(light2a$mean.y/light2a$mean.x)
light2a$type_name <- "low_to_high_light"

light_vcmax_points <- dplyr::bind_rows(light1a[,c("lon","lat","species","middle","type_name")],
                                       light2a[,c("lon","lat","species","middle","type_name")]) 
names(light_vcmax_points) <- c("lon","lat","species","vcmax","type_name")

#jmax
light_jmax <- subset(walker,Treatment_light!="not applicable" & Jmax>0) %>% 
  group_by(lon,lat,species,Treatment_N,Treatment_P,Treatment_CO2,Treatment_light)  %>% summarise(mean = mean(Jmax))

light1 <- subset(light_jmax,Treatment_light=="shade"|Treatment_light=="sun")

light1a <- merge(subset(light1,Treatment_light=="shade"),subset(light1,Treatment_light=="sun"),
                 by=c("lon","lat","species","Treatment_N","Treatment_P","Treatment_CO2"),all.x=TRUE)

#remove ele co2 and high-N
light1a <- subset(light1a,Treatment_CO2 != "ele" & Treatment_N!="high" & Treatment_N!="fertilised")

light1a$middle <- log(light1a$mean.y/light1a$mean.x)
light1a$type_name <- "shade_to_sun"

light2 <- subset(light_jmax,Treatment_light=="high"|Treatment_light=="low")
light2a <- merge(subset(light2,Treatment_light=="low"),subset(light2,Treatment_light=="high"),
                 by=c("lon","lat","species","Treatment_N","Treatment_P","Treatment_CO2"),all.x=TRUE)
#all consistent

light2a$middle <- log(light2a$mean.y/light2a$mean.x)
light2a$type_name <- "low_to_high_light"

light_jmax_points <- dplyr::bind_rows(light1a[,c("lon","lat","species","middle","type_name")],
                                      light2a[,c("lon","lat","species","middle","type_name")]) 

names(light_jmax_points) <- c("lon","lat","species","jmax","type_name")

#ratio
obs_light <- merge(light_vcmax_points,light_jmax_points,by=c("lon","lat","species","type_name"),all.x=TRUE)


######now, combine all observational data
obs_co2$treatment <- "co2"

obs_warming$vcmax <-log(obs_warming$vcmax25_e/obs_warming$vcmax25_a)/obs_warming$t1
obs_warming$jmax <-log(obs_warming$jmax25_e/obs_warming$jmax25_a)/obs_warming$t1
obs_warming_final <- obs_warming[,c("exp","lon","lat","vcmax","jmax")]
obs_warming_final$treatment <- "warming"

obs_light$exp <-paste(obs_light$lon,obs_light$lat,obs_light$species,sep="_")
obs_light_final <- obs_light[,c("exp","lon","lat","vcmax","jmax","type_name")]
obs_light_final$treatment <- "light"

low_high_vcmax$exp <- low_high_vcmax$site_species

low_high_jmax$exp <- low_high_jmax$site_species

Nfer_vj <- merge(low_high_vcmax[,c("vcmax_highN_coef","vcmax_lowN_coef","exp","lon","lat")],
                 low_high_jmax[,c("jmax_highN_coef","jmax_lowN_coef","exp","lon","lat")],
                 by=c("exp","lon","lat"),all.x=TRUE)

obs_lowN <- Nfer_vj[,c("exp","lon","lat","vcmax_highN_coef","jmax_highN_coef")]
names(obs_lowN) <- c("exp","lon","lat","vcmax","jmax")
obs_lowN$treatment <- "lowN"

obs_highN <- Nfer_vj[,c("exp","lon","lat","vcmax_lowN_coef","jmax_lowN_coef")]
names(obs_highN) <- c("exp","lon","lat","vcmax","jmax")
obs_highN$treatment <- "highN"

names(Nfer_kevin) <- c("exp","treatment","ecosystem","vcmax","jmax")

#add reference
obs_lowN$ref <- "Walker 2014 dataset"
obs_lowN$comments <- "vcmax and jmax are sensitivity coefficients of co2 at N fertilization"

obs_highN$ref <- "Walker 2014 dataset"
obs_highN$comments <- "vcmax and jmax are sensitivity coefficients of co2 at N fertilization"

Nfer_kevin$ref <- "GCME_kevin"
Nfer_kevin$comments <- "vcmax and jmax are sensitivity coefficients of co2 at N fertilization"

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

all_obs_pred <- dplyr::bind_rows(obs_co2_pred, obs_warming_pred,obs_light_pred,
                                 obs_lowN,obs_highN,Nfer_kevin) 

all_obs_pred$jmax_vcmax <- all_obs_pred$jmax - all_obs_pred$vcmax
names(all_obs_pred)

all_obs_pred <- all_obs_pred[,c("exp","lon","lat","vcmax","jmax","jmax_vcmax","pred_vcmax","pred_jmax","pred_jmax_vcmax","treatment","ecosystem","type_name",
                                "ref","comments")]

#####before this, it was fully checked!
#finally, for meta-analysis
#anpp, bnpp, nmass, LAI, soil N
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


#filter plots only within vcmax and jmax
photo_plot <- c(unique(all_obs_pred$exp),c("duke2_f","euroface4_pa_f","euroface4_pe_f","euroface4_pn_f","new_zealand_face_f"))
kevin_othervars_cf <- kevin_othervars %>% filter(exp %in% photo_plot)

#leaf_n needs division - but leaf_p doesn't need (as it is all pmass)
nmass <- subset(kevin_othervars_cf,response=="leaf_n") %>%         
  filter(Unit %in% c("g","g_100g","g_g","g_kg","mg_g","mg_kg","mmol_g","ug_mg","umol_g"))
nmass$response <- "nmass"
logr_c_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="c")));
logr_f_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="f")));
logr_cf_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="cf")))

narea <- subset(kevin_othervars_cf,response=="leaf_n") %>%         
  filter(Unit %in% c("g_m2","g_m3","g_pot","mg_cm2","mg_m2","mmol_m2","ug_cm2","ug_cm3","umol_m2"))
narea$response <- "narea"

logr_c_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="c")));
logr_f_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="f")));
logr_cf_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="cf")))

#anpp and anpp_grain are both anpp.
kevin_othervars_cf$response[kevin_othervars_cf$response=="anpp_grain"] <- "anpp"
unique(subset(kevin_othervars_cf,response=="anpp")$Unit)

#bnpp - make sure data between them are not mixed.
kevin_othervars_cf$output[kevin_othervars_cf$response=="coarse_root_production"|
                            kevin_othervars_cf$response=="fine_root_turnover"|
                            kevin_othervars_cf$response=="root_production"|
                            (kevin_othervars_cf$response=="fine_root_biomass" & kevin_othervars_cf$Unit=="g_m2y")|
                            (kevin_othervars_cf$response=="fine_root_production")] <- "bnpp"

bnpp_dataset <- subset(kevin_othervars_cf,output=="bnpp"&(response=="root_production"|exp=="christchurch_pr_c"|exp=="nevada_desert_face_c"|exp=="popface_pa_c"|exp=="popface_pe_c"|exp=="popface_pn_c"| (exp=="facts_ii_face3_pt_c"&response=="fine_root_production")|exp=="biforface_c"))

unique(bnpp_dataset[,c("exp","response","Unit")])

#soil mineral N in dry-mass: 7+1
#7 plots: nh4 + no3
soil_nh4 <- subset(kevin_othervars_cf,response=="soil_nh4-n") %>% group_by(exp,Unit) %>% summarise(co2_a = mean(co2_a), co2_e = mean(co2_e), ambient = mean(ambient), elevated = mean(elevated))
soil_no3 <- subset(kevin_othervars_cf,response=="soil_no3-n") %>% group_by(exp,Unit) %>% summarise(co2_a = mean(co2_a), co2_e = mean(co2_e),ambient = mean(ambient), elevated = mean(elevated))
soil_nh4no3 <- na.omit(merge(soil_nh4,soil_no3,by=c("exp","Unit"),all.x=TRUE))
soil_nh4no3$logr <- log((soil_nh4no3$elevated.x+soil_nh4no3$elevated.y)/(soil_nh4no3$ambient.x+soil_nh4no3$ambient.y))
soil_nh4no3$soil_mineral_N <- log((soil_nh4no3$elevated.x+soil_nh4no3$elevated.y)/(soil_nh4no3$ambient.x+soil_nh4no3$ambient.y))/log(soil_nh4no3$co2_e.x/soil_nh4no3$co2_a.x)

duke2_cf_soil <- (soil_nh4no3$logr[soil_nh4no3$exp=="duke2_cf"]-soil_nh4no3$logr[soil_nh4no3$exp=="duke2_f"])/log(soil_nh4no3$co2_e.x[soil_nh4no3$exp=="duke2_cf"]/soil_nh4no3$co2_a.x[soil_nh4no3$exp=="duke2_cf"])
soil_mineral_1 <- subset(soil_nh4no3,exp!="duke2_cf" & exp!="duke2_f")[,c("exp","soil_mineral_N")]
soil_mineral_2 <- tibble(exp="duke2_cf",soil_mineral_N =duke2_cf_soil)

#1 plot: soil_in
soil_mineral_3 <- tibble(exp="facts_ii_face3_pt_c",soil_mineral_N =agg_meta_sen_coef(response_ratio_v2(subset(kevin_othervars_cf,response=="soil_in")))$middle)
soil_mineral_dry <- rbind(soil_mineral_1,soil_mineral_2,soil_mineral_3)
soil_mineral_dry$type <- "dry"

#soil mineral N in solution: 2+1
solution_mineral <-  agg_meta_sen_coef(response_ratio_v2(subset(kevin_othervars_cf,response=="soil_solution_mineral_n")))
soil_mineral_4 <- tibble(exp=solution_mineral$exp,soil_mineral_N=solution_mineral$middle)

soil_solution_nh4 <- subset(kevin_othervars_cf,response=="soil_solution_nh4") %>% group_by(exp,Unit) %>% summarise(co2_a = mean(co2_a), co2_e = mean(co2_e), ambient = mean(ambient), elevated = mean(elevated))
soil_solution_no3 <- subset(kevin_othervars_cf,response=="soil_solution_no3") %>% group_by(exp,Unit) %>% summarise(co2_a = mean(co2_a), co2_e = mean(co2_e),ambient = mean(ambient), elevated = mean(elevated))
soil_solution_nh4no3 <- merge(soil_solution_nh4,soil_solution_no3,by=c("exp","Unit"),all.x=TRUE)
soil_solution_nh4no3$soil_mineral_N <- log((soil_solution_nh4no3$elevated.x+soil_solution_nh4no3$elevated.y)/(soil_solution_nh4no3$ambient.x+soil_solution_nh4no3$ambient.y))/log(soil_solution_nh4no3$co2_e.x/soil_solution_nh4no3$co2_a.x)
soil_mineral_5 <- tibble(exp=soil_solution_nh4no3$exp,soil_mineral_N=soil_solution_nh4no3$soil_mineral_N)

soil_mineral_wet <- rbind(soil_mineral_4,soil_mineral_5)
soil_mineral_wet$type <- "wet"

soil_mineral <- rbind(soil_mineral_dry,soil_mineral_wet)
soil_mineral_plotmean <- soil_mineral%>% group_by(exp,type)  %>% summarise(soil_mineral_N = mean(soil_mineral_N))



#check numbers of sites and variables 
varname <- kevin_othervars_cf%>% group_by(response)  %>% summarise(number = n())

#other vars - all created now
for (i in 1:nrow(varname)) {
  tryCatch({
    varname1 <- varname$response[i]
    df_c <- subset(kevin_othervars_cf,treatment=="c" & response==varname1)
    assign(paste("logr_c_", varname1,sep=""), as_tibble(response_ratio_v2(df_c)))
    
    df_f <- subset(kevin_othervars_cf,treatment=="f"& response==varname1)
    assign(paste("logr_f_", varname1,sep=""), as_tibble(response_ratio_v2(df_f)))
    
    df_cf <- subset(kevin_othervars_cf,treatment=="cf" & response==varname1)
    assign(paste("logr_cf_", varname1,sep=""), as_tibble(response_ratio_v2(df_cf)))
  }, error=function(e){})} 

#LMA
kevin_LMA <- read.csv("~/data/gcme/kevin/orig_leaf/LMA.csv")
kevin_LMA <- rename(kevin_LMA, c(ambient = x_c, elevated=x_t, ambient_Sd=sd_c, elevated_Sd=sd_t,ambient_Se=se_c,elevated_Se=se_t,n_plots=rep_c,
                                 z=elevation, co2_a=c_c, co2_e=c_t, nfertQ_a = n_c, nfertQ_e = n_t, pfertQ_a = p_c, pfertQ_e = p_t,kfertQ_a = k_c, kfertQ_e = k_t,
                                 warmQ_e1 = w_t1, warmQ_e2 = w_t2, warmQ_e3 = w_t3, Unit=x_units))
kevin_LMA$ambient_Se <- as.numeric(kevin_LMA$ambient_Se)
kevin_LMA$elevated_Se <- as.numeric(kevin_LMA$elevated_Se)
kevin_LMA$sampling_year[kevin_LMA$sampling_year=="2001-2007"] <- 2004
kevin_LMA$sampling_year[kevin_LMA$sampling_year=="2011-2012"] <- 2011
kevin_LMA$sampling_year <- as.numeric(kevin_LMA$sampling_year)
kevin_LMA$Year <- kevin_LMA$sampling_year - kevin_LMA$start_year
kevin_LMA$Year[kevin_LMA$Year<0] <- 0
summary(kevin_LMA$Year)
kevin_LMA <- kevin_LMA %>% filter(exp %in% photo_plot)


#THIS dangerous! However, since we don't know year of measurement, we could only assume it as average
kevin_LMA$Year[is.na(kevin_LMA$Year)==TRUE] <- 2
sla <- subset(kevin_LMA,response=="sla")
lma <- subset(kevin_LMA,response=="lma")
#convert SLA to LMA format 
sla$ambient <- 1/sla$ambient; sla$elevated <- 1/sla$elevated; sla$ambient_Sd <- 1/sla$ambient_Sd
sla$elevated_Sd <- 1/sla$elevated_Sd;sla$ambient_Se <- 1/sla$ambient_Se;sla$elevated_Se <- 1/sla$elevated_Se
LMA <- dplyr::bind_rows(lma,sla)
LMA$response <- "LMA"
LMA2 <- merge(LMA,unique(kevin_othervars_cf[,c("exp","exp_nam")]),by=c("exp"),all.x=TRUE)
logr_c_LMA <- as_tibble(response_ratio_v2(subset(LMA2,treatment=="c")));
logr_f_LMA <-  as_tibble(response_ratio_v2(subset(LMA2,treatment=="f")));
logr_cf_LMA <-  as_tibble(response_ratio_v2(subset(LMA2,treatment=="cf")))


#######final aggregation 
lma_plot <- combine_co2_cf(logr_c_LMA,logr_f_LMA,logr_cf_LMA,"LMA")
narea_plot <- combine_co2_cf(logr_c_narea,logr_f_narea,logr_cf_narea,"narea")
nmass_plot <- combine_co2_cf(logr_c_nmass,logr_f_nmass,logr_cf_nmass,"nmass")

logr_c_bnpp <- response_ratio_v2(bnpp_dataset)
bnpp_plot <- combine_co2_c(logr_c_bnpp,logr_c_bnpp,logr_c_bnpp,"bnpp")
bnpp_plot$bnpp[bnpp_plot$exp=="duke_c"] <- log(subset(bnpp_dataset,exp=="duke_c")$elevated/subset(bnpp_dataset,exp=="duke_c")$ambient)/log(563/363)
#bnpp_plot_final <- merge(bnpp_plot,unique(bnpp_dataset[,c("exp","response")]),by=c("exp"),all.x=TRUE)

soil_mineral_plotmean$condition <- "co2"
soil_mineral_plotmean$condition[soil_mineral_plotmean$exp=="duke2_cf"] <- "fertilization"
soil_mineral_plotmean <- soil_mineral_plotmean[,c("exp","condition","soil_mineral_N","type")]
soil_mineral_plotmean_dry <- subset(soil_mineral_plotmean,type=="dry")

lai_plot <- combine_co2_cf(logr_c_lai,logr_f_lai,logr_cf_lai,"lai")

root_shoot_ratio_plot <- combine_co2_c(logr_c_root_shoot_ratio,logr_f_root_shoot_ratio,logr_cf_root_shoot_ratio,"root_shoot_ratio")

#anpp
#when removing g/m2 - it becomes better
anpp_new <- subset(kevin_othervars_cf,response=="anpp"&Unit!="t_ha"&Unit!="g_m2"&Unit!="gc_m2"&Unit!="mg"&Unit!="g_plant")
anpp_new_c <- response_ratio_v2(subset(anpp_new,treatment=="c"))
anpp_new_f <- response_ratio_v2(subset(anpp_new,treatment=="f"))
anpp_new_cf <- response_ratio_v2(subset(anpp_new,treatment=="cf"))
anpp_plot <- combine_co2_cf(anpp_new_c,anpp_new_f,anpp_new_cf,"anpp")

all_obs_pred <- all_obs_pred %>% rename(condition = treatment)

final_mean <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp","condition"),all.x=TRUE),
                    list(all_obs_pred,lma_plot,narea_plot,nmass_plot,bnpp_plot,root_shoot_ratio_plot,
                         soil_mineral_plotmean_dry[,c("exp","condition","soil_mineral_N")],
                         anpp_plot,lai_plot))

#further process

#2. remove high Tleaf measurement's vcmax (actually, two papers: darbah_et_al_2010a and darbah_et_al_2010b)
final_mean$vcmax[final_mean$exp=="facts_ii_face3_pt_c"] <- agg_meta_sen_coef(response_ratio_v2(subset(logr_c_vcmax,exp=="facts_ii_face3_pt_c" & citation!="darbah_et_al_2010b" & citation!="darbah_et_al_2010a")))$middle
final_mean$vcmax[final_mean$exp=="facts_ii_face4_bp_c"] <- agg_meta_sen_coef(response_ratio_v2(subset(logr_c_vcmax,exp=="facts_ii_face4_bp_c" & citation!="darbah_et_al_2010b")))$middle

#add three anpp plots in rice_face where the unit it g/m2
#only 3 samples - include them

anpp_a <- mean(subset(kevin_othervars,exp=="riceface_japan_ko_2013_3558_13960_c" & response=="anpp_grain")$ambient)
anpp_e <- mean(subset(kevin_othervars,exp=="riceface_japan_ko_2013_3558_13960_c" & response=="anpp_grain")$elevated)
co2_a <- subset(kevin_othervars,exp=="riceface_japan_ko_2013_3558_13960_c" & response=="anpp_grain")$co2_a[1]
co2_e <- subset(kevin_othervars,exp=="riceface_japan_ko_2013_3558_13960_c" & response=="anpp_grain")$co2_e[1]
anpp_r <- log((anpp_e/anpp_a))/log(co2_e/co2_a)

anpp_a <- mean(subset(kevin_othervars,exp=="soyfacesoy1_c" & response=="anpp_grain")$ambient)
anpp_e <- mean(subset(kevin_othervars,exp=="soyfacesoy1_c" & response=="anpp_grain")$elevated)
co2_a <- subset(kevin_othervars,exp=="soyfacesoy1_c" & response=="anpp_grain")$co2_a[1]
co2_e <- subset(kevin_othervars,exp=="soyfacesoy1_c" & response=="anpp_grain")$co2_e[1]
anpp_s1 <- log((anpp_e/anpp_a))/log(co2_e/co2_a)

check <- subset(kevin_othervars,exp=="soyfacesoy2_c")%>% group_by(response,Unit)  %>% summarise(number = n())
anpp_a <- mean(subset(kevin_othervars,exp=="soyfacesoy2_c" & response=="anpp_grain")$ambient)
anpp_e <- mean(subset(kevin_othervars,exp=="soyfacesoy2_c" & response=="anpp_grain")$elevated)
co2_a <- subset(kevin_othervars,exp=="soyfacesoy2_c" & response=="anpp_grain")$co2_a[1]
co2_e <- subset(kevin_othervars,exp=="soyfacesoy2_c" & response=="anpp_grain")$co2_e[1]
anpp_s2 <- log((anpp_e/anpp_a))/log(co2_e/co2_a)

final_mean$anpp[final_mean$exp=="riceface_japan_ko_2013_3558_13960_c"] <- anpp_r
final_mean$anpp[final_mean$exp=="soyfacesoy1_c"] <- anpp_s1
final_mean$anpp[final_mean$exp=="soyfacesoy2_c"] <- anpp_s2

#add popface's soil 
old_data <- read_csv("~/data/gcme/data_received_190325/NewData_wide_CORRECTED2.csv") %>%
  mutate( ambient_Sd  = as.numeric(ambient_Sd),  ambient_Se  = as.numeric(ambient_Se), 
          elevated_Sd = as.numeric(elevated_Sd), elevated_Se = as.numeric(elevated_Se),
          co2_a  = as.numeric(co2_a),  co2_e  = as.numeric(co2_e), 
          ambient  = as.numeric(ambient),  elevated  = as.numeric(elevated))

old_data$exp <- tolower(old_data$exp_nam)
popface <-subset(old_data, (exp_nam=="POPFACE_pa"|exp_nam=="POPFACE_pe"|exp_nam=="POPFACE_pn")&Data_type=="soil_mineral_N")
new_popface <- as.data.frame(agg_meta_sen_coef(response_ratio_v2(popface))[,c("exp","middle")])
names(new_popface) <- c("exp","soil_mineral_N")
final_mean$soil_mineral_N[final_mean$exp=="popface_pa_c"] <- new_popface$soil_mineral_N[new_popface$exp=="popface_pa"]
final_mean$soil_mineral_N[final_mean$exp=="popface_pe_c"] <- new_popface$soil_mineral_N[new_popface$exp=="popface_pe"]
final_mean$soil_mineral_N[final_mean$exp=="popface_pn_c"] <- new_popface$soil_mineral_N[new_popface$exp=="popface_pn"]

final_mean$ecm_type <- NA
final_mean$ecm_type[final_mean$exp=="grassotc"] <- "Nfix" # 2 N-fixing vs. 4 AM
final_mean$ecm_type[final_mean$exp=="biocon_c"]<- "Nfix"
final_mean$ecm_type[final_mean$exp=="swissface_trifolium2_c"]<- "Nfix"
final_mean$ecm_type[final_mean$exp=="new_zealand_face_c"] <- "Nfix"
final_mean$ecm_type[final_mean$exp=="chinaminiface"] <- "Nfix"
final_mean$ecm_type[final_mean$exp=="glycinece"] <- "Nfix"


#corret one outlier of bgb_coarse

#add two plots from Cesar -->after adding them --becoming weaker!
Cesar_anpp <- read.csv("~/data/gcme/cesar/ANPP.csv")
phace <- log(112.6493/111.8533)/log(600/384)
Aspen <- (log(669.6334/517.9303)/log(550/360) + log(621.7131/429.4247)/log(550/360))/2

#further look
final_mean$anpp[final_mean$exp=="phace_c"] <- phace
final_mean$anpp[final_mean$exp=="rhine-aspenface_c"] <- Aspen

###add more data: eucface - LAI
#from https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.13151?casa_token=6CKcWQ_OHHwAAAAA%3AEsLJPJXb45rz2WxE807-NvACiQmFkELScHJiV_eaRUEPd0psT7co5ZnJp8Mo7CKaPFt4H6dkKe8XqZFXJw
eucface_lai_df <- read.csv("~/data/Duursma_gcb/EucFACE_DUURSMA_GCB_LEAFAREAINDEX/data/FACE_RA_P0037_GAPFRACLAI_20121026-20150225_L2.csv")
#150 as given in paper, also consistent with our df
subset(logr_c_vcmax,exp=="eucface_c")$co2_e[1];subset(logr_c_vcmax,exp=="eucface_c")$co2_a[1]
eucface_lai <- log(mean(subset(eucface_lai_df,treatment=="elevated")$LAI,na.rm=TRUE)/mean(subset(eucface_lai_df,treatment=="ambient")$LAI,na.rm=TRUE))/log(540/394)

###biocon_c - LAI
#LMA, everything is LMA (g/cm2) now, though shown as cm2/g - converting to g/m2
#lma: g/m2
lma_a <-10000*mean(subset(logr_c_LMA,exp=="biocon_c"&Unit=="cm2/g")$ambient)
lma_e <-10000*mean(subset(logr_c_LMA,exp=="biocon_c"&Unit=="cm2/g")$elevated,na.rm=TRUE)
check <- subset(kevin_othervars,exp=="biocon_c")%>% group_by(response,Unit)  %>% summarise(number = n())
#agb: g/m2
agb_a <- mean(subset(kevin_othervars,exp=="biocon_c" & response=="agb" & Unit=="g_m2")$ambient)
agb_e <- mean(subset(kevin_othervars,exp=="biocon_c" & response=="agb" & Unit=="g_m2")$elevated)
agb_a/lma_a;agb_e/lma_e # looks ok
subset(logr_c_vcmax,exp=="biocon_c")$co2_a;subset(logr_c_vcmax,exp=="biocon_c")$co2_e
biocon_lai <- log((agb_e/lma_e)/(agb_a/lma_a))/log(570/367)

#??? biocon_c - root/shoot = bgb/agb?
bgb_a <- mean(subset(kevin_othervars,exp=="biocon_c" & response=="bgb" & Unit=="g_m2")$ambient)
bgb_e <- mean(subset(kevin_othervars,exp=="biocon_c" & response=="bgb" & Unit=="g_m2")$elevated)
biocon_root_shoot <- log((bgb_e/agb_e)/(bgb_a/agb_a))/log(570/367)

###popface_pa_c - anpp - shown in gC/m2 but should be gC/m2/yr
check <- subset(kevin_othervars,exp=="popface_pa_c")%>% group_by(response,Unit)  %>% summarise(number = n())
anpp_a <- mean(subset(kevin_othervars,exp=="popface_pa_c" & response=="anpp")$ambient)
anpp_e <- mean(subset(kevin_othervars,exp=="popface_pa_c" & response=="anpp")$elevated)
popface_pa_c_ANPP <- log((anpp_e/anpp_a))/log(550/368)

###popface_pe_c - anpp - shown in gC/m2 but should be gC/m2/yr
check <- subset(kevin_othervars,exp=="popface_pe_c")%>% group_by(response,Unit)  %>% summarise(number = n())
anpp_a <- mean(subset(kevin_othervars,exp=="popface_pe_c" & response=="anpp")$ambient)
anpp_e <- mean(subset(kevin_othervars,exp=="popface_pe_c" & response=="anpp")$elevated)
popface_pe_c_ANPP <- log((anpp_e/anpp_a))/log(550/368)

###popface_pn_c - anpp - shown in gC/m2 but should be gC/m2/yr
check <- subset(kevin_othervars,exp=="popface_pn_c")%>% group_by(response,Unit)  %>% summarise(number = n())
anpp_a <- mean(subset(kevin_othervars,exp=="popface_pn_c" & response=="anpp")$ambient)
anpp_e <- mean(subset(kevin_othervars,exp=="popface_pn_c" & response=="anpp")$elevated)
popface_pn_c_ANPP <- log((anpp_e/anpp_a))/log(550/368)

###brandbjerg_c - LAI and ANPP missing after investigation
check <- subset(kevin_othervars,exp=="brandbjerg_c")%>% group_by(response,Unit)  %>% summarise(number = n())

###euroface4_pn_c
check <- subset(kevin_othervars,exp=="euroface4_pn_c")%>% group_by(response,Unit)  %>% summarise(number = n())
#LAI
LAI_a <- mean(subset(kevin_othervars,exp=="euroface4_pn_c" & response=="lai_max"& citation=="liberloo_et_al_2006")$ambient,na.rm=TRUE)
LAI_e <- mean(subset(kevin_othervars,exp=="euroface4_pn_c" & response=="lai_max"& citation=="liberloo_et_al_2006")$elevated,na.rm=TRUE)
euroface4_pn_lai <- log((LAI_e/LAI_a))/log(550/368)

#bnpp derived from SI of: https://www.pnas.org/doi/10.1073/pnas.0706518104#supplementary-materials
#Table S3
bnpp_e <- 2.46
bnpp_a <- 1.76
euroface4_pn_bnpp <- log((2.46/1.76))/log(550/368)

#####facts_ii_face3_pt_c
#anpp still missing
check <- subset(kevin_othervars,exp=="facts_ii_face3_pt_c")%>% group_by(response,Unit)  %>% summarise(number = n())

#??? root/shoot = bgb/agb?
bgb_a <- mean(subset(kevin_othervars,exp=="facts_ii_face3_pt_c" & response=="bgb" & Unit=="g_m2")$ambient)
bgb_e <- mean(subset(kevin_othervars,exp=="facts_ii_face3_pt_c" & response=="bgb" & Unit=="g_m2")$elevated)
agb_a <- mean(subset(kevin_othervars,exp=="facts_ii_face3_pt_c" & response=="agb" & Unit=="g_m2")$ambient)
agb_e <- mean(subset(kevin_othervars,exp=="facts_ii_face3_pt_c" & response=="agb" & Unit=="g_m2")$elevated)
co2_a <- mean(c(subset(kevin_othervars,exp=="facts_ii_face3_pt_c" & response=="bgb" & Unit=="g_m2")$co2_a,subset(kevin_othervars,exp=="facts_ii_face3_pt_c" & response=="agb" & Unit=="g_m2")$co2_a))
co2_e <- mean(c(subset(kevin_othervars,exp=="facts_ii_face3_pt_c" & response=="bgb" & Unit=="g_m2")$co2_e,subset(kevin_othervars,exp=="facts_ii_face3_pt_c" & response=="agb" & Unit=="g_m2")$co2_e))

facts_ii_face3_pt_c_root_shoot <- log((bgb_e/agb_e)/(bgb_a/agb_a))/log(co2_e/co2_a)

####new_zealand
check <- subset(kevin_othervars,exp=="new_zealand_face_c")%>% group_by(response,Unit)  %>% summarise(number = n())
#??? root/shoot = bgb/agb?
#using g/g is better because it is completely mass unit with reasonable values
bgb_a <-subset(kevin_othervars,exp=="new_zealand_face_c" & response=="bgb_c")$ambient
bgb_e <- subset(kevin_othervars,exp=="new_zealand_face_c" & response=="bgb_c")$elevated
agb_a <- subset(kevin_othervars,exp=="new_zealand_face_c" & response=="agb_c")$ambient
agb_e <-subset(kevin_othervars,exp=="new_zealand_face_c" & response=="agb_c")$elevated

#bgb_a <- mean(subset(kevin_othervars,exp=="new_zealand_face_c" & response=="bgb" & Unit=="g_m2")$ambient)
#bgb_e <- mean(subset(kevin_othervars,exp=="new_zealand_face_c" & response=="bgb" & Unit=="g_m2")$elevated)
#agb_a <- mean(subset(kevin_othervars,exp=="new_zealand_face_c" & response=="agb" & Unit=="g_m2" & citation=="allard_et_al_2005")$ambient)
#agb_e <- mean(subset(kevin_othervars,exp=="new_zealand_face_c" & response=="agb" & Unit=="g_m2" & citation=="allard_et_al_2005")$elevated)
new_zealand_c_root_shoot <- log((bgb_e/agb_e)/(bgb_a/agb_a))/log(475/364)

#LAI from LMA and agb
lma_a <-mean(subset(logr_c_LMA,exp=="new_zealand_face_c"&Unit=="g/m_")$ambient)
lma_e <-mean(subset(logr_c_LMA,exp=="new_zealand_face_c"&Unit=="g/m_")$elevated)
#agb: g/m2
agb_a <- mean(subset(kevin_othervars,exp=="new_zealand_face_c" & response=="agb" & Unit=="g_m2")$ambient)
agb_e <- mean(subset(kevin_othervars,exp=="new_zealand_face_c" & response=="agb" & Unit=="g_m2")$elevated)
agb_a/lma_a;agb_e/lma_e # looks ok
new_zealand_face_c_lai <- log((agb_e/lma_e)/(agb_a/lma_a))/log(475/364)

####duke_c
check <- subset(kevin_othervars,exp=="duke_c")%>% group_by(response,Unit)  %>% summarise(number = n())

#LMA obtained from leaf_c in area and mass basis
Cmass_a <- mean(subset(kevin_othervars,exp=="duke_c" & response=="leaf_c" & Unit=="g_100g" & citation=="hamilton_et_al_2004")$ambient)
Cmass_e <- mean(subset(kevin_othervars,exp=="duke_c" & response=="leaf_c" & Unit=="g_100g" & citation=="hamilton_et_al_2004")$elevated)

Carea_a <- mean(subset(kevin_othervars,exp=="duke_c" & response=="leaf_c" & Unit=="gc_m2" & citation=="hamilton_et_al_2002")$ambient)
Carea_e <- mean(subset(kevin_othervars,exp=="duke_c" & response=="leaf_c" & Unit=="gc_m2" & citation=="hamilton_et_al_2002")$elevated)

#LMA = Carea/Cmass
duke_c_lma <- log((Carea_e/Cmass_e)/(Carea_a/Cmass_a))/log(563/363)

#??? root/shoot = bgb/agb?
bgb_a <- mean(subset(kevin_othervars,exp=="duke_c" & response=="bgb" & Unit=="g_m2" &co2_a==363)$ambient)
bgb_e <- mean(subset(kevin_othervars,exp=="duke_c" & response=="bgb" & Unit=="g_m2"&co2_e==563)$elevated)
agb_a <- mean(subset(kevin_othervars,exp=="duke_c" & response=="agb" & Unit=="gc_m2"& citation=="schlesinger_and_lichter_2001")$ambient) # too keep column constant
agb_e <- mean(subset(kevin_othervars,exp=="duke_c" & response=="agb" & Unit=="gc_m2" & citation=="schlesinger_and_lichter_2001")$elevated)
duke_c_root_shoot <- log((bgb_e/agb_e)/(bgb_a/agb_a))/log(563/363)

#LAI still missing
check <- subset(kevin_othervars,exp=="duke2_c")%>% group_by(response,Unit)  %>% summarise(number = n())

###duke2_c
#??? root/shoot = bgb/agb?
bgb_a <- mean(subset(kevin_othervars,exp=="duke2_c" & response=="bgb" & Unit=="g_m2" & co2_a==363)$ambient)
bgb_e <- mean(subset(kevin_othervars,exp=="duke2_c" & response=="bgb" & Unit=="g_m2" & co2_e==563)$elevated)
agb_a <- mean(subset(kevin_othervars,exp=="duke2_c" & response=="agb" & Unit=="g_m2")$ambient,na.rm=TRUE)
agb_e <- mean(subset(kevin_othervars,exp=="duke2_c" & response=="agb" & Unit=="g_m2")$elevated,na.rm=TRUE)
duke2_c_root_shoot <- log((bgb_e/agb_e)/(bgb_a/agb_a))/log(563/363)

#ANPP and BNPP. See: Table S1 of https://nph.onlinelibrary.wiley.com/doi/10.1111/j.1469-8137.2009.03078.x
#McCarthy HR, Oren R, Johnsen KH, Gallet‐Budynek A, Pritchard SG, Cook CW, et al. 2010. Re‐assessment of plant carbon dynamics at the Duke free‐air CO2 enrichment site: interactions of atmospheric [CO2] with nitrogen and water availability over stand development. New Phytologist, 185(2): 514-528.
#it has foliage, branch, stem and coarse root production
#we compare R1 and R2, where R1 is ambient and R2 is elevated (for 9 years)
foliage_a <- 131+144+135+177+182+194+174+173+207; foliage_e <-140+176+203+210+222+224+211+201+255
branch_a <- 83+81+65+84+99+74+36+67+72; branch_e <- 77+92+84+102+110+87+45+68+74
stem_a <- 467+468+378+492+592+451+220+409+446;stem_e <- 455+553+512+628+691+558+292+441+490
root_a <- 112+109+91+112+145+105+53+91+103;root_e <- 109+128+126+150+161+131+67+99+110

duke2_c_anpp <- log((foliage_e+branch_e+stem_e)/(foliage_a+branch_a+stem_a))/log(563/363)
duke2_c_bnpp <- log((root_e)/(root_a))/log(563/363)

######ignore duke2_cf

########euroface4_pe_c
#LAI
check <- subset(kevin_othervars,exp=="euroface4_pe_c")%>% group_by(response,Unit)  %>% summarise(number = n())
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

#bnpp not be able capture!

######ornerp_liqui_c: Narea, LMA still missing
check <- subset(kevin_othervars,exp=="ornerp_liqui_c")%>% group_by(response,Unit)  %>% summarise(number = n())

#LAI calculated from: https://data.ess-dive.lbl.gov/view/doi%3A10.15485%2F1480325#ess-dive-ddd1cfa81a329ba-20181119T143441660
lai_amb <- mean(as.numeric(subset(read.csv("~/data/ORNL_FACE/ORNL_ax.csv"),co2=="AMB")$ptrait_lai),na.rm=T)
lai_elv <- mean(as.numeric(subset(read.csv("~/data/ORNL_FACE/ORNL_ax.csv"),co2=="ELE")$ptrait_lai),na.rm=T)

ornerp_liqui_lai <- log((lai_elv/lai_amb))/log(560/380)

######soyfacesoy2_c: bnpp and root/shoot still missing
check <- subset(kevin_othervars,exp=="soyfacesoy2_c")%>% group_by(response,Unit)  %>% summarise(number = n())
#narea = lma * nmass /10
narea_a <- mean(subset(logr_c_LMA,exp=="soyfacesoy2_c"&Unit=="cm_/g")$ambient)*mean(subset(kevin_othervars,exp=="soyfacesoy2_c"&response=="leaf_n" & Unit=="g_kg")$ambient)/10
narea_e <- mean(subset(logr_c_LMA,exp=="soyfacesoy2_c"&Unit=="cm_/g")$elevated)*mean(subset(kevin_othervars,exp=="soyfacesoy2_c"&response=="leaf_n" & Unit=="g_kg")$elevated)/10
soyfacesoy2_narea <- log((narea_e/narea_a))/log(548/373)


######facts_ii_face4_bp_c missing jmax (only has J), anpp, bnpp and root/shoot missing
check <- subset(kevin_othervars,exp=="facts_ii_face4_bp_c")%>% group_by(response,Unit)  %>% summarise(number = n())

#nevada_desert_face_c: missing Narea, LMA, root/shoot
check <- subset(kevin_othervars,exp=="nevada_desert_face_c")%>% group_by(response,Unit)  %>% summarise(number = n())
# it is Scrubland, so can be considered as grassland

#anpp, unit is wrong. so not used
#anpp_a <- mean(subset(kevin_othervars,exp=="nevada_desert_face_c"&response=="anpp")$ambient)
#anpp_e <- mean(subset(kevin_othervars,exp=="nevada_desert_face_c"&response=="anpp")$elevated)
#nevada_desert_face_c_anpp <- log((anpp_e/anpp_a))/log(550/364)

#Leaf area obtained from https://link.springer.com/content/pdf/10.1007/s10021-005-0124-4.pdf Housman D. C., Naumburg E., Huxman T.E., Charlet T.N., Nowak R.S., Smith S.D. (2006) Increase in Desert Shrub Produvtivity under Elevated Carbon Dioxide Vary with Water Availability. Ecosystems 9: 374-385
#only has ratio of LA
LA_ratio <- mean(1.18,0.92,1.06,1.28,1.10,1.33)
leafCN_ratio <- mean(subset(kevin_othervars,exp=="nevada_desert_face_c"&response=="leaf_cn")$elevated/subset(kevin_othervars,exp=="nevada_desert_face_c"&response=="leaf_cn")$ambient)
leafNmass_ratio <- mean(subset(kevin_othervars,exp=="nevada_desert_face_c"&response=="leaf_n")$elevated/subset(kevin_othervars,exp=="nevada_desert_face_c"&response=="leaf_n")$ambient)

LMA_ratio <- LA_ratio/(leafNmass_ratio*leafCN_ratio)
narea_ratio <- LMA_ratio*leafNmass_ratio
nevada_desert_face_c_lma <- log(LMA_ratio)/log(550/364)
nevada_desert_face_c_narea <- log(narea_ratio)/log(550/364)

#######euroface4_pa_c
check <- subset(kevin_othervars,exp=="euroface4_pa_c")%>% group_by(response,Unit)  %>% summarise(number = n())

#lai
LAI_a <- mean(subset(kevin_othervars,exp=="euroface4_pa_c" & response=="lai_max"& citation=="liberloo_et_al_2006")$ambient,na.rm=TRUE)
LAI_e <- mean(subset(kevin_othervars,exp=="euroface4_pa_c" & response=="lai_max"& citation=="liberloo_et_al_2006")$elevated,na.rm=TRUE)
euroface4_pa_lai <- log((LAI_e/LAI_a))/log(550/368)

#anpp, narea, nmass, lma lacked 
subset(kevin_othervars,exp=="euroface4_pa_c")$dominant_species[1]
#Na, Nm and LMA cannot be filled by https://www.sciencedirect.com/science/article/pii/S0269749106005434
# Marinari S, Calfapietra C, De Angelis P, Mugnozza GS, Grego S (2007) Impact of elevated CO2 and nitrogen fertilization on foliar elemental composition in a short rotation poplar plantation. Environmental Pollution 147:507–515
# it only has controlled condition, but not ambient condition

#######mi: missing nmass, narea, lma
check <- subset(kevin_othervars,exp=="mi_c")%>% group_by(response,Unit)  %>% summarise(number = n())

#root/shoot = leaf biomass/bgb biomass
#make sure co2 condition is the same
root_shoot_a <- mean(subset(kevin_othervars,exp=="mi_c"&response=="leaf_biomass"&co2_e==760)$ambient)/mean(subset(kevin_othervars,exp=="mi_c"&response=="bgb"&co2_e==760)$ambient)
root_shoot_e <- mean(subset(kevin_othervars,exp=="mi_c"&response=="leaf_biomass"&co2_e==760)$elevated)/mean(subset(kevin_othervars,exp=="mi_c"&response=="bgb"&co2_e==760)$elevated)

mi_root_shoot <- log((root_shoot_e/root_shoot_a))/log(760/410)

######soyfacesoy1_c: Nmass, Narea, bnpp and root/shoot missing
check <- subset(kevin_othervars,exp=="soyfacesoy1_c")%>% group_by(response,Unit)  %>% summarise(number = n())

#lma = leaf biomass (g/m2) /LAI
leaf_biomass_a <- mean(subset(kevin_othervars,exp=="soyfacesoy1_c" & response=="leaf_biomass")$ambient)
leaf_biomass_e <- mean(subset(kevin_othervars,exp=="soyfacesoy1_c" & response=="leaf_biomass")$elevated)

LAI_a <- mean(subset(kevin_othervars,exp=="soyfacesoy1_c" & response=="lai")$ambient)
LAI_e <- mean(subset(kevin_othervars,exp=="soyfacesoy1_c" & response=="lai")$elevated)

soyfacesoy1_lma <- log((leaf_biomass_e/LAI_e)/(leaf_biomass_a/LAI_a))/log(552/371)

######soyfacetobacco9_c: missing LAI, anpp, bnpp, root/shoot
check <- subset(kevin_othervars,exp=="soyfacetobacco9_c")%>% group_by(response,Unit)  %>% summarise(number = n())

#narea = nmass * LMA
soyfacetobacco9_narea <- final_mean$nmass[final_mean$exp=="soyfacetobacco9_c"] + final_mean$LMA[final_mean$exp=="soyfacetobacco9_c"]

#######giface - filling LAI (m2/m2) and LMA, Nmass and Narea
LAI_a <- mean(subset(kevin_othervars,exp=="giface_c" & response=="lai_max")$ambient)
LAI_e <- mean(subset(kevin_othervars,exp=="giface_c" & response=="lai_max")$elevated)

giface_c_lai <- log((LAI_e/LAI_a))/log(450/380)

#lma = agb (g/m2) /LAI in grassland
leaf_biomass_a <- mean(subset(kevin_othervars,exp=="giface_c" & response=="agb" & Unit=="g_m2")$ambient)
leaf_biomass_e <- mean(subset(kevin_othervars,exp=="giface_c" & response=="agb"& Unit=="g_m2")$elevated)

giface_c_lma <- log((leaf_biomass_e/LAI_e)/(leaf_biomass_a/LAI_a))/log(450/380)

#Nmass and Narea were filled by below publication
#Carbon dioxide fertilisation and supressed respiration induce enhanced spring biomass production in a mixed species temperate meadow exposed to moderate carbon dioxide enrichment
# in table 5, see Nmass (mg/g), Cmass was also recorded but not put here.
#co2 changes from 400 to 480
Nmass_a <- (34+41.2+26.4+36.3+31.8+43.5)/6/1000 # averaged them, then converted from mg/g to g/g
Nmass_e <- (32.7+40.7+24.2+33.3+31.4+42.4)/6/1000 # in mg/g 
giface_c_nmass <- log((Nmass_e/Nmass_a))/log(480/400)

#since Narea = Nmass * LMA # one paper increases from 380 to 450, another increases from 400 to 480. Not different too much - here we just used one of them
giface_c_narea <- log((Nmass_e*leaf_biomass_e/LAI_e)/(Nmass_a*leaf_biomass_a/LAI_a))/log(480/400)

#######soyfacemaiz4_c cannot be filled - just one old paper without available info
aa <- subset(kevin_othervars,exp=="riceface_japan_ko_2012_3558_13960_c")%>% group_by(response,Unit)  %>% summarise(number = n())
subset(kevin_othervars,response=="lai_max")%>% group_by(exp)  %>% summarise(number = n())

#rice_face- just keeps its original data - because it marks year and coordinates too precisely, so it is dangerous for filling data.

#include all!
#final5 <- final4 
final5 <- final_mean
final5$lai[final5$exp=="eucface_c"] <- eucface_lai
final5$lai[final5$exp=="biocon_c"] <- biocon_lai
final5$lai[final5$exp=="euroface4_pa_c"] <- euroface4_pa_lai
final5$lai[final5$exp=="euroface4_pn_c"] <- euroface4_pn_lai
final5$lai[final5$exp=="euroface4_pe_c"] <- euroface4_pe_lai
final5$lai[final5$exp=="new_zealand_face_c"] <- new_zealand_face_c_lai
final5$lai[final5$exp=="ornerp_liqui_c"] <- ornerp_liqui_lai
final5$lai[final5$exp=="giface_c"] <- giface_c_lai

final5$root_shoot_ratio[final5$exp=="biocon_c"] <- biocon_root_shoot
final5$root_shoot_ratio[final5$exp=="facts_ii_face3_pt_c"] <- facts_ii_face3_pt_c_root_shoot
final5$root_shoot_ratio[final5$exp=="new_zealand_face_c"] <- new_zealand_c_root_shoot
final5$root_shoot_ratio[final5$exp=="duke_c"] <- duke_c_root_shoot
final5$root_shoot_ratio[final5$exp=="duke2_c"] <- duke2_c_root_shoot
final5$root_shoot_ratio[final5$exp=="mi_c"] <- mi_root_shoot

#this can be imputed - as it was derived from literatures well
final5$anpp[final5$exp=="duke2_c"] <- duke2_c_anpp
final5$bnpp[final5$exp=="duke2_c"] <- duke2_c_bnpp
final5$bnpp[final5$exp=="euroface4_pn_c"] <- euroface4_pn_bnpp

#anpp, bnpp not imputed now - since their unit are somewhere wrong (gc/m2)
#final5$anpp[final5$exp=="euroface4_pa_c"] <- popface_pa_c_ANPP
#final5$anpp[final5$exp=="euroface4_pn_c"] <- popface_pn_c_ANPP
#final5$anpp[final5$exp=="euroface4_pe_c"] <- popface_pe_c_ANPP
#final5$anpp[final5$exp=="nevada_desert_face_c"] <- nevada_desert_face_c_anpp

final5$LMA[final5$exp=="duke_c"] <- duke_c_lma
final5$LMA[final5$exp=="euroface4_pe_c"] <- euroface4_pe_lma
final5$LMA[final5$exp=="nevada_desert_face_c"] <- nevada_desert_face_c_lma
final5$LMA[final5$exp=="soyfacesoy1_c"] <- soyfacesoy1_lma
final5$LMA[final5$exp=="giface_c"] <- giface_c_lma

final5$narea[final5$exp=="soyfacesoy2_c"] <- soyfacesoy2_narea
final5$narea[final5$exp=="nevada_desert_face_c"] <- nevada_desert_face_c_narea
final5$narea[final5$exp=="soyfacetobacco9_c"] <- soyfacetobacco9_narea
final5$narea[final5$exp=="giface_c"] <- giface_c_narea

final5$nmass[final5$exp=="giface_c"] <- giface_c_nmass

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


#check some very high points with at Nfer condition
subset(final5,condition=="Fertilization")
#duke2_cf: jmax,lma, narea, lai.
#new_zealand_face_cf: vcmax, narea, nmass, lma

#duke2_cf
a1 <- subset(logr_f_jmax,exp=="duke2_f")[,c("ambient","elevated","logr","citation")]
#remove the one with 4.2, which cannot be parallel to cf on below, and the values looked wrong
a1$logr[a1$citation=="oishi_et_al_2014"] <- NA

a2 <- subset(logr_cf_jmax,exp=="duke2_cf")[,c("ambient","elevated","logr","citation")]

#recalculate logr = logcf - cf
duke2_cf_jmax <- (mean(a2$logr) - mean(a1$logr,na.rm=TRUE))/log(563/363)

#  LMA and narea were corrected by maier_et_al_2008 -https://academic.oup.com/treephys/article/28/4/597/1720624
#corrected by their table 2 (using 1-year old data)
duke2_cf_lma <- log(113.3/100.6)/log(563/363)

duke2_cf_narea <- log(1.34/1.44)/log(563/363)

subset(logr_c_lai,exp=="duke2_c")[,c("ambient","elevated","logr","citation")]
subset(logr_f_lai,exp=="duke2_f")[,c("ambient","elevated","logr","citation")]
subset(logr_cf_lai,exp=="duke2_cf")[,c("ambient","elevated","logr","citation")]
# oishi_et_al_2014 - lai too low - below has problem - ambient should be 3.8 not 38!
#calculate it directly 
duke2_cf_lai <- log(4.8/4.05)/log(563/363)

final5$jmax[final5$exp=="duke2_cf"] <- duke2_cf_jmax
final5$jmax_vcmax[final5$exp=="duke2_cf"] <- final5$jmax[final5$exp=="duke2_cf"] -final5$vcmax[final5$exp=="duke2_cf"] 

final5$LMA[final5$exp=="duke2_cf"] <- duke2_cf_lma
final5$narea[final5$exp=="duke2_cf"] <- duke2_cf_narea
final5$lai[final5$exp=="duke2_cf"] <- duke2_cf_lai

#new_zealand_face_cf - not wrong
#all is on below
# https://www.publish.csiro.au/fp/pdf/PP01009
#all correct

csvfile <- paste("~/data/gcme/MS_data/plot_data.csv")
write.csv(final5, csvfile, row.names = TRUE)

#check if all ambient condition at f and cf has same data - so that we can do cf -f directly
subset(final5,condition=="Fertilization")

#vcmax, jmax already checked
#1. check if both have the same ambient values (duke2_cf already been corrected)
#vcmax, jmax, lma, narea, nmass, bnpp, anpp, lai 
subset(kevin2_f_vcmax,exp=="duke2_f")[,c("ambient","elevated")];subset(kevin2_cf_vcmax,exp=="duke2_cf")[,c("ambient","elevated")]
subset(kevin2_f_jmax,exp=="duke2_f")[,c("ambient","elevated")];subset(kevin2_cf_jmax,exp=="duke2_cf")[,c("ambient","elevated")]

subset(kevin2_f_vcmax,exp=="euroface4_pa_f")[,c("ambient","elevated")];subset(kevin2_cf_vcmax,exp=="euroface4_pa_cf")[,c("ambient","elevated")]
subset(kevin2_f_jmax,exp=="euroface4_pa_f")[,c("ambient","elevated")];subset(kevin2_cf_jmax,exp=="euroface4_pa_cf")[,c("ambient","elevated")]

subset(kevin2_f_vcmax,exp=="euroface4_pe_f")[,c("ambient","elevated")];subset(kevin2_cf_vcmax,exp=="euroface4_pe_cf")[,c("ambient","elevated")]
subset(kevin2_f_jmax,exp=="euroface4_pe_f")[,c("ambient","elevated")];subset(kevin2_cf_jmax,exp=="euroface4_pe_cf")[,c("ambient","elevated")]

subset(kevin2_f_vcmax,exp=="euroface4_pn_f")[,c("ambient","elevated")];subset(kevin2_cf_vcmax,exp=="euroface4_pn_cf")[,c("ambient","elevated")]
subset(kevin2_f_jmax,exp=="euroface4_pn_f")[,c("ambient","elevated")];subset(kevin2_cf_jmax,exp=="euroface4_pn_cf")[,c("ambient","elevated")]

subset(kevin2_f_vcmax,exp=="new_zealand_face_f")[,c("ambient","elevated")];subset(kevin2_cf_vcmax,exp=="new_zealand_face_cf")[,c("ambient","elevated")]
#subset(kevin2_f_jmax,exp=="new_zealand_face_f")[,c("ambient","elevated")];subset(kevin2_cf_jmax,exp=="new_zealand_face_cf")[,c("ambient","elevated")]

subset(logr_f_LMA,exp=="euroface4_pn_f")[,c("ambient","elevated","logr","citation")]
subset(logr_cf_LMA,exp=="euroface4_pn_cf")[,c("ambient","elevated","logr","citation")]

subset(logr_f_LMA,exp=="new_zealand_face_f")[,c("ambient","elevated","logr","citation")]
subset(logr_cf_LMA,exp=="new_zealand_face_cf")[,c("ambient","elevated","logr","citation")]

subset(logr_f_narea,exp=="euroface4_pe_f")[,c("ambient","elevated","logr","citation")]
subset(logr_cf_narea,exp=="euroface4_pe_cf")[,c("ambient","elevated","logr","citation")]
subset(logr_f_narea,exp=="euroface4_pn_f")[,c("ambient","elevated","logr","citation")]
subset(logr_cf_narea,exp=="euroface4_pn_cf")[,c("ambient","elevated","logr","citation")]
subset(logr_f_narea,exp=="new_zealand_face_f")[,c("ambient","elevated","logr","citation")]
subset(logr_cf_narea,exp=="new_zealand_face_cf")[,c("ambient","elevated","logr","citation")]

subset(logr_f_nmass,exp=="euroface4_pe_f")[,c("ambient","elevated","logr","citation")]
subset(logr_cf_nmass,exp=="euroface4_pe_cf")[,c("ambient","elevated","logr","citation")]
subset(logr_f_nmass,exp=="euroface4_pn_f")[,c("ambient","elevated","logr","citation")]
subset(logr_cf_nmass,exp=="euroface4_pn_cf")[,c("ambient","elevated","logr","citation")]
subset(logr_f_nmass,exp=="new_zealand_face_f")[,c("ambient","elevated","logr","citation")]
subset(logr_cf_nmass,exp=="new_zealand_face_cf")[,c("ambient","elevated","logr","citation")]

subset(logr_f_anpp,exp=="euroface4_pa_f")[,c("ambient","elevated","logr","citation")]
subset(logr_cf_anpp,exp=="euroface4_pa_cf")[,c("ambient","elevated","logr","citation")]
subset(logr_f_anpp,exp=="euroface4_pe_f")[,c("ambient","elevated","logr","citation")]
subset(logr_cf_anpp,exp=="euroface4_pe_cf")[,c("ambient","elevated","logr","citation")]
subset(logr_f_anpp,exp=="euroface4_pn_f")[,c("ambient","elevated","logr","citation")]
subset(logr_cf_anpp,exp=="euroface4_pn_cf")[,c("ambient","elevated","logr","citation")]
