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
#calcualte logr = log (elevated/ambient). But set logr = 0 because it is very less likely that the values are 0 under treatments 
#logr_var primarily calculated from ROM
#and alternatively, if logr_var was shown NA in plot - then relying on percentage (logr_var / logr) of other samples within the same plot
response_ratio_v2 <- function(df){
  #-----------------------------------------------------------------------
  # Input: 
  # name: df..
  # df should include: c("exp_nam","ambient","ambient_Sd","ambient_Se",
  # "elevated","elevated_Se","elevated_Sd","n_plots")
  #-----------------------------------------------------------------------
  df_c_sub <- df %>%         
    ## Here only for my variables, selected as described above 
    # get standard deviation for all data
    mutate( my_ambient_sd = ambient_Sd, my_elevated_sd = elevated_Sd ) %>%
    rowwise() %>% 
    mutate( my_ambient_sd   = ifelse( is.na(my_ambient_sd),  ambient_Se  * sqrt(n_plots), my_ambient_sd ),
            my_elevated_sd  = ifelse( is.na(my_elevated_sd), elevated_Se * sqrt(n_plots), my_elevated_sd )) %>%
    
    ## Get logarithm of response ratio and its variance
    metafor::escalc( 
      measure = "ROM", 
      m1i = elevated, sd1i = my_elevated_sd, n1i = n_plots, 
      m2i = ambient,  sd2i = my_ambient_sd,  n2i = n_plots, 
      data=., 
      append = TRUE, var.names = c("logr", "logr_var") ) %>% 
    as_tibble() %>% 
    mutate( logr_se = sqrt(logr_var)/sqrt(n_plots) )
  
  df_c_sub$id <- 1:nrow(df_c_sub)
  df_c_sub$logr[df_c_sub$logr==0] <- NA
  df_c_sub$logr_var[df_c_sub$logr_var==0] <- NA
  
  explist <- unique(df_c_sub$exp)
  df_c_sub$new_logr_var_percentage <- NA
  for (i in 1:length(explist)){ # if SD shown NA in plot - then relying on percentage of other samples within the same plot
    percentage <- mean(abs(subset(df_c_sub,exp==explist[i])$logr_var/subset(df_c_sub,exp==explist[i])$logr), na.rm=TRUE)
    df_c_sub$new_logr_var_percentage[df_c_sub$exp==explist[i]] <- percentage
  }
  
  df_c_sub$logr_var[is.na(df_c_sub$logr_var)==TRUE] <- abs(df_c_sub$new_logr_var_percentage[is.na(df_c_sub$logr_var)==TRUE] * df_c_sub$logr[is.na(df_c_sub$logr_var)==TRUE])
  df_c_sub$logr_var[df_c_sub$logr_var=="NaN"] <- NA
  df_c_sub$logr_var[is.na(df_c_sub$logr)==TRUE] <- NA
  return(df_c_sub)
  #-----------------------------------------------------------------------
  # Output: df_c_sub
  #-----------------------------------------------------------------------
}

#agg_meta_sen_coef: calculate sensitivity coefficient from individuals dataset to site-mean dataset
agg_meta_sen_coef <- function(df){
  threshold <- 0.5
  explist <- unique(df$exp)
  mylist <- list() #create an empty list
  #convert to sensitivity coef
  df$sen_coef <- df$logr/log(df$co2_e/df$co2_a)
  df$sen_coef_var <- df$logr_var/log(df$co2_e/df$co2_a)
  df$sen_coef_se <- df$logr_se/log(df$co2_e/df$co2_a)
  
  df$logr <- df$sen_coef
  df$logr_var <- df$sen_coef_var
  df$logr_se <- df$sen_coef_se
  #all using mean value
  for (i in 1:length(explist)){
    
    mean_value <- mean(subset(df,exp==explist[i])$logr,na.rm=TRUE)
    n_samples <- length(subset(df,exp==explist[i])$logr_var)
    df_box <- tibble(
      exp=explist[i], middle = mean_value,
      ymin   = NA,
      ymax   = NA,
      ymin_quarter = NA,
      ymax_quarter   = NA,
      variance_info ="No",no=n_samples,logr_var=NA)
    mylist[[i]] <- df_box}
  output <- do.call("rbind",mylist)
  return(output)
}

#agg_meta: calculate response ratio from individuals dataset to site-mean dataset
agg_meta <- function(df){
  threshold <- 0.5
  explist <- unique(df$exp)
  mylist <- list() 
  for (i in 1:length(explist)){
    mean_value <- mean(subset(df,exp==explist[i])$logr,na.rm=TRUE)
    n_samples <- length(subset(df,exp==explist[i])$logr_var)
    df_box <- tibble(
      exp=explist[i], middle = mean_value,
      ymin   = NA,
      ymax   = NA,
      ymin_quarter = NA,
      ymax_quarter   = NA,
      variance_info ="No",no=n_samples,logr_var = NA)
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


#####3. CO2 effect data
#read Kevin
kevin <- read.csv("/Users/yunpeng/data/gcme/kevin/orig_vcmax/JunkePeng_11252021.csv")

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
smith_co2 <- read.csv("/Users/yunpeng/data/smith_keenan_gcb/gcb_co2/co2_data.csv")
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
smith_all_plotmean$exp[smith_all_plotmean$exp=="PineOTC"] <- "Christchurch_pr_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="RichmondGH"] <- "australiacotton_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="UIAPine"] <- "ua_otc_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="PHACE"] <- "phace_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="AspenFACE"] <- "rhine-aspenface_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="BilyKriz"] <- "bily_kriz_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="Headley"] <- "headley_qp_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="Viesalm"] <- "vielsalm_c"

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
walker <- read.csv("/Users/yunpeng/data/leaf_traits/Walker/LEAF_PHOTOSYNTHESIS_TRAITS_1224/data/Leaf_Photosynthesis_Traits.csv")

Nfer_vcmax <- subset(walker,Treatment_N!="not applicable" &Treatment_CO2!="not applicable" & Vcmax>0) %>% 
  group_by(lon,lat,species,Treatment_N,Treatment_CO2)  %>% summarise(mean = mean(Vcmax))

highN_vcmax <- merge(subset(Nfer_vcmax,Treatment_N=="high"&Treatment_CO2=="amb"),
                     subset(Nfer_vcmax,Treatment_N=="high"&Treatment_CO2=="ele"),
                     by=c("lon","lat","species","Treatment_N"),all.x=TRUE)

highN_vcmax$middle <- log(highN_vcmax$mean.y/highN_vcmax$mean.x)

highN_vcmax$site_species <- paste(highN_vcmax$lon,highN_vcmax$lat,highN_vcmax$species,sep="_")
highN_vcmax[,c("site_species","middle")]

lowN_vcmax <- merge(subset(Nfer_vcmax,Treatment_N=="low"&Treatment_CO2=="amb"),
                    subset(Nfer_vcmax,Treatment_N=="low"&Treatment_CO2=="ele"),
                    by=c("lon","lat","species","Treatment_N"),all.x=TRUE)

lowN_vcmax$middle <- log(lowN_vcmax$mean.y/lowN_vcmax$mean.x)
lowN_vcmax$site_species <- paste(lowN_vcmax$lon,lowN_vcmax$lat,lowN_vcmax$species,sep="_")


low_high_vcmax <- merge(highN_vcmax[,c("lon","lat","site_species","middle")],na.omit(lowN_vcmax[,c("lon","lat","site_species","middle")]),
                    by=c("lon","lat","site_species"),all.x=TRUE)
names(low_high_vcmax) <- c("lon","lat","site_species","vcmax_highN","vcmax_lowN")
low_high_vcmax$co2_change <- NA

#input co2 change
low_high_vcmax$co2_change[low_high_vcmax$lon==-84.04] <- log(560/360) 
# https://watermark.silverchair.com/22-5-321.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAv4wggL6BgkqhkiG9w0BBwagggLrMIIC5wIBADCCAuAGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMlQzU_4AV0WQkXxioAgEQgIICsZUQMf48_Ogae5yJOUenILkeX5E_nHvKWARnNsbV4Cd5VZ38yYrQBL5duSx1kSq97qPt5PcU5fhbJvGTYY3HdCO1nm4_DwncLNNjvauVf6SSlXs6ilSIzkeclSN2eRY7uzwXfa-JaJg8RHzLq8scOxja5YhaUCRogGCZDpTXryJ97Z1jobEOPG-YCdDI0apFWymgoaqcE286a1lG6-aQK_7kbxwinsZSf1_CcxhItSBy7LTZRxwANCisRFwFr0IXbqpDnoQJH8cBepblQkkZzcEfwclOsghsN26D0uNF4K4UTTpsHUd2jifg0ga_4fzNsP_iMmht62Y0HwRC7uSc1gD587zvIQGQoMa4yUBpPeHYXOLmcGNpURkWQOA-qfRnjAC_7PS8UUMQbSQp7LPD2ZpuAKxNMONSXTT4SwqlQYhxFyHQZyaQYEoaKe-XBSaVjdt2oA1IXfz7OZ72-KsJweC-buZHrjy64iCPtaJ7R0KsAuLtbO3Wsgu5DZyr3oddoZ30fuk2NVogC-7iHPvFfYxu0uqkRwQavdcIYvxVZvt3NkI534Wce9fIgQzoWr2PSbqJoenwM3Bs0QpemX2ByyxFWtCt-7qOcRXu71Iupj5Qt0S6DkcIqlEwEkruKStMO5MUv6CEtpeUaY2VyZ6kKnlh8ffsMkFgKFdARaaRUKc4Igp2JY-ptvQrLInZKoZlebgjHiiwWVQvDB4TqhmXi33pJ_iiriuNTd2_8AKChx7b777LNYKDri48i4Lny425sdVdxrlY_t8123c3OnnVTgmpTD8lUTRKh85iAnKcrts6hSFFeGrC9nM8zuFh0vpA2p8gq0wUqO5knbwNrqIcW_ZXV8v0cTmyNJjnKK7PvwKGq7HFBYfHLk1ZBUJNy7duVXwQ5nqvddVLFmZff_Dv2_-1
low_high_vcmax$co2_change[low_high_vcmax$lon==20.00] <- log(700/350) 
#https://www.researchgate.net/publication/242459413_Nutrient_and_genotypic_effects_on_CO2-responsiveness_Photosynthetic_regulation_in_Leucadendron_species_of_a_nutrient-poor_environment
low_high_vcmax$co2_change[low_high_vcmax$lon==89.14] <- log(720/360) 

low_high_vcmax$vcmax_highN_coef <- low_high_vcmax$vcmax_highN/low_high_vcmax$co2_change
low_high_vcmax$vcmax_lowN_coef <- low_high_vcmax$vcmax_lowN/low_high_vcmax$co2_change

#do the same for jmax

Nfer_jmax <- subset(walker,Treatment_N!="not applicable" &Treatment_CO2!="not applicable" & Jmax>0) %>% 
  group_by(lon,lat,species,Treatment_N,Treatment_CO2)  %>% summarise(mean = mean(Jmax))

highN_jmax <- merge(subset(Nfer_jmax,Treatment_N=="high"&Treatment_CO2=="amb"),
                    subset(Nfer_jmax,Treatment_N=="high"&Treatment_CO2=="ele"),
                    by=c("lon","lat","species","Treatment_N"),all.x=TRUE)

highN_jmax$middle <- log(highN_jmax$mean.y/highN_jmax$mean.x)

highN_jmax$type_name <- "CO2 response at high N"
highN_jmax$site_species <- paste(highN_jmax$lon,highN_jmax$lat,highN_jmax$species,sep="_")

lowN_jmax <- merge(subset(Nfer_jmax,Treatment_N=="low"&Treatment_CO2=="amb"),
                   subset(Nfer_jmax,Treatment_N=="low"&Treatment_CO2=="ele"),
                   by=c("lon","lat","species","Treatment_N"),all.x=TRUE)

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
df3 <- read_csv("/Users/yunpeng/data/Kumarathunge_2020_newphy/kumarathunge_2020_newphy/Data/PPC-TGlob_V1.0.csv")
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
df1 <- read_csv("/Users/yunpeng/data/Kumarathunge_2020_newphy/kumarathunge_2020_newphy/Data/ACi-TGlob_V1.0.csv")
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
walker <- read.csv("/Users/yunpeng/data/leaf_traits/Walker/LEAF_PHOTOSYNTHESIS_TRAITS_1224/data/Leaf_Photosynthesis_Traits.csv")

#vcmax
light_vcmax <- subset(walker,Treatment_light!="not applicable") %>% 
  group_by(lon,lat,species,Treatment_light)  %>% summarise(mean = mean(Vcmax))

#divided into two categories:
#(1) shade vs. sun
#(2) high vs. low light

light1 <- subset(light_vcmax,Treatment_light=="shade"|Treatment_light=="sun")
light1a <- merge(subset(light1,Treatment_light=="shade"),subset(light1,Treatment_light=="sun"),
                 by=c("lon","lat","species"),all.x=TRUE)
light1a$middle <- log(light1a$mean.y/light1a$mean.x)
light1a$type_name <- "shade_to_sun"

light2 <- subset(light_vcmax,Treatment_light=="high"|Treatment_light=="low")
light2a <- merge(subset(light2,Treatment_light=="low"),subset(light2,Treatment_light=="high"),
                 by=c("lon","lat","species"),all.x=TRUE)
light2a$middle <- log(light2a$mean.y/light2a$mean.x)
light2a$type_name <- "low_to_high_light"

light_vcmax_points <- dplyr::bind_rows(light1a[,c("lon","lat","species","middle","type_name")],
                                       light2a[,c("lon","lat","species","middle","type_name")]) 
names(light_vcmax_points) <- c("lon","lat","species","vcmax","type_name")

light_jmax <- subset(walker,Treatment_light!="not applicable" & Jmax>0) %>% 
  group_by(lon,lat,species,Treatment_light)  %>% summarise(mean = mean(Jmax))

light1 <- subset(light_jmax,Treatment_light=="shade"|Treatment_light=="sun")
light1a <- merge(subset(light1,Treatment_light=="shade"),subset(light1,Treatment_light=="sun"),
                 by=c("lon","lat","species"),all.x=TRUE)
light1a$middle <- log(light1a$mean.y/light1a$mean.x)
light1a$type_name <- "shade_to_sun"

light2 <- subset(light_jmax,Treatment_light=="high"|Treatment_light=="low")
light2a <- merge(subset(light2,Treatment_light=="low"),subset(light2,Treatment_light=="high"),
                 by=c("lon","lat","species"),all.x=TRUE)
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
obs_light_final <- obs_light[,c("exp","lon","lat","vcmax","jmax")]
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
obs_light_pred <- merge(obs_light_final[,c("exp","lon","lat","vcmax","jmax")],prediction, by=c("lon","lat"),all.x=TRUE)
obs_light_pred$exp <- obs_light_pred$exp.x


all_obs_pred <- dplyr::bind_rows(obs_co2_pred, obs_warming_pred,obs_light_pred,
                            obs_lowN,obs_highN,Nfer_kevin) 

all_obs_pred$jmax_vcmax <- all_obs_pred$jmax - all_obs_pred$vcmax
names(all_obs_pred)

all_obs_pred <- all_obs_pred[,c("exp","lon","lat","vcmax","jmax","jmax_vcmax","pred_vcmax","pred_jmax","pred_jmax_vcmax","treatment","ecosystem",
                                "ref","comments")]

csvfile <- paste("/Users/yunpeng/data/gcme/MS_data/vcmax_obs_pred.csv")
write.csv(all_obs_pred, csvfile, row.names = TRUE)