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

#Read all original data csv and downloading forcing data
#####1. GCME (CO2 + warming), sent from Kevin
##this file has better format than MESI.csv, especially about the measurement year information
kevin <- read.csv("~/data/gcme/kevin/orig_vcmax/JunkePeng_11252021.csv")
#1. correct exp_nam to make it consistent with GCME
kevin$site[is.na(kevin$site)==TRUE] <- "riceface_japan_a_2003_3938_14057"
kevin$exp_nam <- kevin$site #new plot: capo_caccia,itatinga,itatinga_experimental_station,setres
#unify sitename from japan 
kevin[grep("riceface_japan", kevin$exp_nam),]$exp_nam <- "riceface_japan"

kevin <- rename(kevin, c(ambient = x_c, elevated=x_t, ambient_Sd=sd_c, elevated_Sd=sd_t,ambient_Se=se_c,elevated_Se=se_t,n_plots=rep_c,
                         z=elevation, co2_a=c_c, co2_e=c_t, nfertQ_a = n_c, nfertQ_e = n_t, pfertQ_a = p_c, pfertQ_e = p_t,kfertQ_a = k_c, kfertQ_e = k_t,
                         warmQ_e1 = w_t1, warmQ_e2 = w_t2, warmQ_e3 = w_t3, Unit=x_units))

#correct a few sampling year --> when looking at org csv
kevin$sampling_year[is.na(kevin$sampling_year)==TRUE & kevin$site=="brandbjerg"] <- 2011
kevin$sampling_year[is.na(kevin$sampling_year)==TRUE & kevin$site=="popface"] <- 2002
kevin$sampling_year[is.na(kevin$sampling_year)==TRUE & kevin$site=="biocon"] <- 2005
kevin$sampling_year[kevin$sampling_year=="2005-2010"] <- 2008
kevin$sampling_year[kevin$sampling_year=="1996-2010"] <- 2003
kevin$sampling_year[kevin$sampling_year=="2003-2006"] <- 2005
kevin$sampling_year[kevin$citation=="domec_et_al_2012"] <- 1997 #  by looking at their info it says +1y. then we assume it is 1996+1
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
kevin_z$elv <- kevin_z$z

kevin_forcing_sitemean <- aggregate(kevin_z,by=list(kevin_z$lon,kevin_z$lat,kevin_z$z,kevin_z$year_start,kevin_z$year_end), FUN=mean, na.rm=TRUE)
kevin_forcing_sitemean <- kevin_forcing_sitemean[,c("lon","lat","z","year_start","year_end")]

kevin_forcing_sitemean$sitename <- paste("kevin",1:nrow(kevin_forcing_sitemean),sep="")
kevin_forcing_sitemean$elv <- kevin_forcing_sitemean$z
siteinfo_final <- kevin_forcing_sitemean[,c("sitename","lon","lat","elv","year_start","year_end")]

#this was universally saved - no change anymore!
csvfile <- paste("~/data/gcme/kevin/forcing/forcing_info.csv")
write.csv(siteinfo_final, csvfile, row.names = TRUE)

devtools::load_all("~/yunkepeng/gcme/pmodel/ingestr/")
devtools::load_all("~/yunkepeng/gcme/pmodel/rsofun/")

#for (i in 1:(nrow(siteinfo_final))){
#  df_watch <- ingest(
#    siteinfo  = siteinfo_final[i,],
#    source    = "watch_wfdei",
#    getvars   = c("temp", "prec", "ppfd", "vpd", "patm"), 
#    dir       = "/Volumes/My Passport/data/watch_wfdei/",
#    settings  = list(correct_bias = "worldclim", dir_bias = "/Volumes/My Passport/data/worldclim/"))
#  df_cru <- ingest(
#    siteinfo  = siteinfo_final[i,],
#    source    = "cru",
#    getvars   = "ccov",
#    dir       = "/Volumes/My Passport/data/cru/ts_4.01/")
#  df_co2 <- ingest(
#    siteinfo  = siteinfo_final[i,],
#    source  = "co2_mlo",
#    verbose = FALSE)
#  df_co2_final <- as.data.frame(df_co2$data)
#  df_co2_final2 <- df_co2_final[!(format(df_co2_final$date,"%m") == "02" & format(df_co2_final$date, "%d") == "29"), , drop = FALSE] # make columns consistent
#  
#  co2 <- df_co2_final2$co2
#  
#  ddf_meteo <- as_tibble(cbind(as.data.frame(df_watch$data),as.data.frame(df_cru$data)[,c("ccov_int","ccov")],co2))
#  
#  ddf_meteo$fapar <- 1
#  ddf_meteo$tmax <- ddf_meteo$temp
#  ddf_meteo$tmin <- ddf_meteo$temp
#  
#  ddf_meteo$sitename <- siteinfo_final[i,c("sitename")]
#  ddf_meteo$lon <- siteinfo_final[i,c("lon")]
#  ddf_meteo$lat <-  siteinfo_final[i,c("lat")]
#  ddf_meteo$elv <-  siteinfo_final[i,c("elv")]
#  ddf_meteo$year_start <-  siteinfo_final[i,c("year_start")]
#  ddf_meteo$year_end <-  siteinfo_final[i,c("year_end")]
#  
#  csvfile <- paste("~/data/gcme/kevin/forcing/climate/",siteinfo_final$sitename[i],".csv",sep = "")
#  write.csv(ddf_meteo, csvfile, row.names = TRUE)
#  print(i)    
#}

vcmax25_warmingco2_siteinfo <- merge(kevin_z,siteinfo_final,by=c("lon","lat","elv","year_start","year_end"),all.x=TRUE)

df_soiltexture <- bind_rows(
  top    = tibble(layer = "top",    fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1),
  bottom = tibble(layer = "bottom", fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1))
params_modl <- list(
  kphio           = 0.09423773,
  soilm_par_a     = 0.33349283,
  soilm_par_b     = 1.45602286,
  tau_acclim_tempstress = 10,
  par_shape_tempstress  = 0.0)

forcing_path <- "~/data/gcme/kevin/forcing/climate/"
vcmax25_warmingco2_siteinfo$elv <- vcmax25_warmingco2_siteinfo$z
#for conversion
vcmax25_warmingco2_siteinfo$warmQ_e2[vcmax25_warmingco2_siteinfo$exp=="riceface_japan_l_2008_3938_14057_w"] <- vcmax25_warmingco2_siteinfo$warmQ_e3[vcmax25_warmingco2_siteinfo$exp=="riceface_japan_l_2008_3938_14057_w"] 

for (i in 1:nrow(vcmax25_warmingco2_siteinfo)){
  siteinfo_site <- vcmax25_warmingco2_siteinfo[i,c("sitename","lon","lat","elv","year_start","year_end")]
  sitename <- siteinfo_site$sitename
  
  df1 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for ambient
  df2 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for elevated
  
  #Q1 - add warming directly
  df2$temp <- df2$temp+vcmax25_warmingco2_siteinfo$warmQ_e2[i]
  
  #Q2 - if having co2 then replacing it directly - if not then using co2 at measurement year
  if (is.na(vcmax25_warmingco2_siteinfo$co2_a[i])==TRUE|is.na(vcmax25_warmingco2_siteinfo$co2_e[i])==TRUE) { 
    df1$co2 <- df1$co2
    df2$co2 <- df2$co2
  } else {
    df1$co2 <- vcmax25_warmingco2_siteinfo$co2_a[i] # set to ambient 
    df2$co2 <- vcmax25_warmingco2_siteinfo$co2_e[i] # set to elevated 
  }
  
  modlist1 <- run_pmodel_f_bysite( 
    sitename=sitename, 
    params_siml <- list(#ambient
      spinup             = TRUE,
      spinupyears        = 10,
      recycle            = 1,
      soilmstress        = TRUE,
      tempstress         = TRUE,
      calc_aet_fapar_vpd = FALSE,
      in_ppfd            = TRUE,
      in_netrad          = FALSE,
      outdt              = 1,
      ltre               = FALSE,
      ltne               = FALSE,
      ltrd               = FALSE,
      ltnd               = FALSE,
      lgr3               = TRUE,
      lgn3               = FALSE,
      lgr4               = FALSE,
      firstyeartrend = siteinfo_site$year_start,
      nyeartrend = siteinfo_site$year_end-siteinfo_site$year_start+1), 
    siteinfo_site, 
    df1, 
    df_soiltexture, 
    params_modl = params_modl, 
    makecheck = TRUE)
  
  modlist2 <- run_pmodel_f_bysite( #elevated
    sitename=sitename, 
    params_siml <- list(
      spinup             = TRUE,
      spinupyears        = 10,
      recycle            = 1,
      soilmstress        = TRUE,
      tempstress         = TRUE,
      calc_aet_fapar_vpd = FALSE,
      in_ppfd            = TRUE,
      in_netrad          = FALSE,
      outdt              = 1,
      ltre               = FALSE,
      ltne               = FALSE,
      ltrd               = FALSE,
      ltnd               = FALSE,
      lgr3               = TRUE,
      lgn3               = FALSE,
      lgr4               = FALSE,
      firstyeartrend = siteinfo_site$year_start,
      nyeartrend = siteinfo_site$year_end-siteinfo_site$year_start+1), 
    siteinfo_site, 
    df2, 
    df_soiltexture, 
    params_modl = params_modl, 
    makecheck = TRUE)
  #mean or max?
  mean_vcmax25_ambient <- max(modlist1$vcmax25)*1000000
  mean_vcmax25_elevated <- max(modlist2$vcmax25)*1000000
  mean_vcmax_ambient <- max(modlist1$vcmax)*1000000
  mean_vcmax_elevated <- max(modlist2$vcmax)*1000000
  
  mean_jmax25_ambient <- max(modlist1$jmax25)*1000000
  mean_jmax25_elevated <- max(modlist2$jmax25)*1000000
  mean_jmax_ambient <- max(modlist1$jmax)*1000000
  mean_jmax_elevated <- max(modlist2$jmax)*1000000
  
  vcmax25_warmingco2_siteinfo$mean_vcmax25_ambient[i] <- mean_vcmax25_ambient
  vcmax25_warmingco2_siteinfo$mean_vcmax25_elevated[i] <- mean_vcmax25_elevated
  vcmax25_warmingco2_siteinfo$mean_vcmax_ambient[i] <- mean_vcmax_ambient
  vcmax25_warmingco2_siteinfo$mean_vcmax_elevated[i] <- mean_vcmax_elevated
  vcmax25_warmingco2_siteinfo$pred_vcmax_logr[i] <- log(vcmax25_warmingco2_siteinfo$mean_vcmax_elevated[i]/vcmax25_warmingco2_siteinfo$mean_vcmax_ambient[i])
  vcmax25_warmingco2_siteinfo$pred_vcmax_coef[i] <- (log(vcmax25_warmingco2_siteinfo$mean_vcmax_elevated[i]/vcmax25_warmingco2_siteinfo$mean_vcmax_ambient[i]))/(log(vcmax25_warmingco2_siteinfo$co2_e[i]/vcmax25_warmingco2_siteinfo$co2_a[i]))
  vcmax25_warmingco2_siteinfo$pred_vcmax25_logr[i] <- log(vcmax25_warmingco2_siteinfo$mean_vcmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_vcmax25_ambient[i])
  vcmax25_warmingco2_siteinfo$pred_vcmax25_coef[i] <- (log(vcmax25_warmingco2_siteinfo$mean_vcmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_vcmax25_ambient[i]))/(log(vcmax25_warmingco2_siteinfo$co2_e[i]/vcmax25_warmingco2_siteinfo$co2_a[i]))
  
  vcmax25_warmingco2_siteinfo$mean_jmax25_ambient[i] <- mean_jmax25_ambient
  vcmax25_warmingco2_siteinfo$mean_jmax25_elevated[i] <- mean_jmax25_elevated
  vcmax25_warmingco2_siteinfo$mean_jmax_ambient[i] <- mean_jmax_ambient
  vcmax25_warmingco2_siteinfo$mean_jmax_elevated[i] <- mean_jmax_elevated
  vcmax25_warmingco2_siteinfo$pred_jmax_logr[i] <- log(vcmax25_warmingco2_siteinfo$mean_jmax_elevated[i]/vcmax25_warmingco2_siteinfo$mean_jmax_ambient[i])
  vcmax25_warmingco2_siteinfo$pred_jmax_coef[i] <- (log(vcmax25_warmingco2_siteinfo$mean_jmax_elevated[i]/vcmax25_warmingco2_siteinfo$mean_jmax_ambient[i]))/(log(vcmax25_warmingco2_siteinfo$co2_e[i]/vcmax25_warmingco2_siteinfo$co2_a[i]))
  vcmax25_warmingco2_siteinfo$pred_jmax25_logr[i] <- log(vcmax25_warmingco2_siteinfo$mean_jmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_jmax25_ambient[i])
  vcmax25_warmingco2_siteinfo$pred_jmax25_coef[i] <- (log(vcmax25_warmingco2_siteinfo$mean_jmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_jmax25_ambient[i]))/(log(vcmax25_warmingco2_siteinfo$co2_e[i]/vcmax25_warmingco2_siteinfo$co2_a[i]))
  
  vcmax25_warmingco2_siteinfo$PPFD[i] <- mean(df1$ppfd,na.rm=TRUE)*1000000 # in umol/m2/s
  vcmax25_warmingco2_siteinfo$vpd[i] <- mean(df1$vpd,na.rm=TRUE)/1000 #in kPa
  vcmax25_warmingco2_siteinfo$temp[i] <- mean(df1$temp,na.rm=TRUE) #in degree celcius
}
summary(vcmax25_warmingco2_siteinfo)

kevin2_final <- vcmax25_warmingco2_siteinfo

kevin2_c_vcmax <- subset(kevin2_final, treatment=="c" & response =="vcmax")
kevin2_c_jmax <- subset(kevin2_final, treatment=="c" & response =="jmax")

#each site set the same co2_a and co2_e
unique(kevin2_c_vcmax[,c("exp","co2_a","co2_e")])

vcmax_kevin <- kevin2_c_vcmax %>% group_by(lon,lat,exp)  %>% summarise(vcmax = mean(pred_vcmax_coef))
jmax_kevin <- kevin2_c_vcmax %>% group_by(lon,lat,exp)  %>% summarise(jmax = mean(pred_jmax_coef))

#####2. Smith: co2_data
smith_co2 <- read.csv("~/data/smith_keenan_gcb/gcb_co2/co2_data.csv")
smith_co2$exp_nam <- smith_co2$SiteID
smith_co2 <- subset(smith_co2,is.na(lat)==FALSE) # remove empty column

#create year from 1992 to 2016 -- same to GCME
smith_co2$year_start <- 1992
smith_co2$year_end <- 2016

nick_forcing_sitemean <- aggregate(smith_co2,by=list(smith_co2$lon,smith_co2$lat,smith_co2$z,smith_co2$year_start,smith_co2$year_end), FUN=mean, na.rm=TRUE)
nick_forcing_sitemean <- nick_forcing_sitemean[,c("lon","lat","z","year_start","year_end")]
dim(nick_forcing_sitemean)

nick_forcing_sitemean$sitename <- paste("nick",1:nrow(nick_forcing_sitemean),sep="")
nick_forcing_sitemean$elv <- nick_forcing_sitemean$z
siteinfo_final <- nick_forcing_sitemean[,c("sitename","lon","lat","elv","year_start","year_end")]

csvfile <- paste("/Users/yunpeng/data/smith_keenan_gcb/gcb_co2/forcing/forcing_info.csv")
write.csv(siteinfo_final, csvfile, row.names = TRUE)

#for (i in 1:(nrow(siteinfo_final))){
#  df_watch <- ingest(
#    siteinfo  = siteinfo_final[i,],
#    source    = "watch_wfdei",
#    getvars   = c("temp", "prec", "ppfd", "vpd", "patm"), 
#    dir       = "/Volumes/My Passport/data/watch_wfdei/",
#    settings  = list(correct_bias = "worldclim", dir_bias = "/Volumes/My Passport/data/worldclim/"))
#  df_cru <- ingest(
#    siteinfo  = siteinfo_final[i,],
#    source    = "cru",
#    getvars   = "ccov",
#    dir       = "/Volumes/My Passport/data/cru/ts_4.01/")
#  df_co2 <- ingest(
#    siteinfo  = siteinfo_final[i,],
#    source  = "co2_mlo",
#    verbose = FALSE)
#  df_co2_final <- as.data.frame(df_co2$data)
#  df_co2_final2 <- df_co2_final[!(format(df_co2_final$date,"%m") == "02" & format(df_co2_final$date, "%d") == "29"), , drop = FALSE] # make columns consistent
#  co2 <- df_co2_final2$co2
#  ddf_meteo <- as_tibble(cbind(as.data.frame(df_watch$data),as.data.frame(df_cru$data)[,c("ccov_int","ccov")],co2))
#  ddf_meteo$fapar <- 1
#  ddf_meteo$tmax <- ddf_meteo$temp
#  ddf_meteo$tmin <- ddf_meteo$temp
#  ddf_meteo$sitename <- siteinfo_final[i,c("sitename")]
#  ddf_meteo$lon <- siteinfo_final[i,c("lon")]
#  ddf_meteo$lat <-  siteinfo_final[i,c("lat")]
#  ddf_meteo$elv <-  siteinfo_final[i,c("elv")]
#  ddf_meteo$year_start <-  siteinfo_final[i,c("year_start")]
#  ddf_meteo$year_end <-  siteinfo_final[i,c("year_end")]
#  csvfile <- paste("~/data/smith_keenan_gcb/gcb_co2/forcing/climate/",siteinfo_final$sitename[i],".csv",sep = "")
#  write.csv(ddf_meteo, csvfile, row.names = TRUE)
#  print(i)    
#}

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

smith_co2_final <- rename(smith_co2_final, c(ambient = ambient_vcmax, elevated=Vcmax, ambient_Se=ambient_vcmaxSE, elevated_Se=SE,
                                             co2_a=ambient_co2, co2_e=Ca,
                                             ambient_j = ambient_jmax, elevated_j =Jmax, ambient_Se_j=ambient_jmaxSE, elevated_Se_j=Jmax_SE))

smith_co2_final$elv <- smith_co2_final$z

vcmax25_warmingco2_siteinfo <- merge(smith_co2_final,siteinfo_final,by=c("lon","lat","elv"),all.x=TRUE)

forcing_path <- "~/data/smith_keenan_gcb/gcb_co2/forcing/climate/"

for (i in 1:nrow(vcmax25_warmingco2_siteinfo)){
  siteinfo_site <- vcmax25_warmingco2_siteinfo[i,c("sitename","lon","lat","elv","year_start","year_end")]
  sitename <- siteinfo_site$sitename
  
  df1 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for ambient
  df2 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for elevated
  
  #Q2 - if having co2 then replacing it directly - if not then using co2 at measurement year
  if (is.na(vcmax25_warmingco2_siteinfo$co2_a[i])==TRUE|is.na(vcmax25_warmingco2_siteinfo$co2_e[i])==TRUE) { 
    df1$co2 <- df1$co2
    df2$co2 <- df2$co2
  } else {
    df1$co2 <- vcmax25_warmingco2_siteinfo$co2_a[i] # set to ambient 
    df2$co2 <- vcmax25_warmingco2_siteinfo$co2_e[i] # set to elevated 
  }
  
  modlist1 <- run_pmodel_f_bysite( 
    sitename=sitename, 
    params_siml <- list(#ambient
      spinup             = TRUE,
      spinupyears        = 10,
      recycle            = 1,
      soilmstress        = TRUE,
      tempstress         = TRUE,
      calc_aet_fapar_vpd = FALSE,
      in_ppfd            = TRUE,
      in_netrad          = FALSE,
      outdt              = 1,
      ltre               = FALSE,
      ltne               = FALSE,
      ltrd               = FALSE,
      ltnd               = FALSE,
      lgr3               = TRUE,
      lgn3               = FALSE,
      lgr4               = FALSE,
      firstyeartrend = siteinfo_site$year_start,
      nyeartrend = siteinfo_site$year_end-siteinfo_site$year_start+1), 
    siteinfo_site, 
    df1, 
    df_soiltexture, 
    params_modl = params_modl, 
    makecheck = TRUE)
  
  modlist2 <- run_pmodel_f_bysite( #elevated
    sitename=sitename, 
    params_siml <- list(
      spinup             = TRUE,
      spinupyears        = 10,
      recycle            = 1,
      soilmstress        = TRUE,
      tempstress         = TRUE,
      calc_aet_fapar_vpd = FALSE,
      in_ppfd            = TRUE,
      in_netrad          = FALSE,
      outdt              = 1,
      ltre               = FALSE,
      ltne               = FALSE,
      ltrd               = FALSE,
      ltnd               = FALSE,
      lgr3               = TRUE,
      lgn3               = FALSE,
      lgr4               = FALSE,
      firstyeartrend = siteinfo_site$year_start,
      nyeartrend = siteinfo_site$year_end-siteinfo_site$year_start+1), 
    siteinfo_site, 
    df2, 
    df_soiltexture, 
    params_modl = params_modl, 
    makecheck = TRUE)
  #mean or max?
  mean_vcmax25_ambient <- max(modlist1$vcmax25)*1000000
  mean_vcmax25_elevated <- max(modlist2$vcmax25)*1000000
  mean_vcmax_ambient <- max(modlist1$vcmax)*1000000
  mean_vcmax_elevated <- max(modlist2$vcmax)*1000000
  
  mean_jmax25_ambient <- max(modlist1$jmax25)*1000000
  mean_jmax25_elevated <- max(modlist2$jmax25)*1000000
  mean_jmax_ambient <- max(modlist1$jmax)*1000000
  mean_jmax_elevated <- max(modlist2$jmax)*1000000
  
  vcmax25_warmingco2_siteinfo$mean_vcmax25_ambient[i] <- mean_vcmax25_ambient
  vcmax25_warmingco2_siteinfo$mean_vcmax25_elevated[i] <- mean_vcmax25_elevated
  vcmax25_warmingco2_siteinfo$mean_vcmax_ambient[i] <- mean_vcmax_ambient
  vcmax25_warmingco2_siteinfo$mean_vcmax_elevated[i] <- mean_vcmax_elevated
  vcmax25_warmingco2_siteinfo$pred_vcmax_logr[i] <- log(vcmax25_warmingco2_siteinfo$mean_vcmax_elevated[i]/vcmax25_warmingco2_siteinfo$mean_vcmax_ambient[i])
  vcmax25_warmingco2_siteinfo$pred_vcmax_coef[i] <- (log(vcmax25_warmingco2_siteinfo$mean_vcmax_elevated[i]/vcmax25_warmingco2_siteinfo$mean_vcmax_ambient[i]))/(log(vcmax25_warmingco2_siteinfo$co2_e[i]/vcmax25_warmingco2_siteinfo$co2_a[i]))
  vcmax25_warmingco2_siteinfo$pred_vcmax25_logr[i] <- log(vcmax25_warmingco2_siteinfo$mean_vcmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_vcmax25_ambient[i])
  vcmax25_warmingco2_siteinfo$pred_vcmax25_coef[i] <- (log(vcmax25_warmingco2_siteinfo$mean_vcmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_vcmax25_ambient[i]))/(log(vcmax25_warmingco2_siteinfo$co2_e[i]/vcmax25_warmingco2_siteinfo$co2_a[i]))
  
  vcmax25_warmingco2_siteinfo$mean_jmax25_ambient[i] <- mean_jmax25_ambient
  vcmax25_warmingco2_siteinfo$mean_jmax25_elevated[i] <- mean_jmax25_elevated
  vcmax25_warmingco2_siteinfo$mean_jmax_ambient[i] <- mean_jmax_ambient
  vcmax25_warmingco2_siteinfo$mean_jmax_elevated[i] <- mean_jmax_elevated
  vcmax25_warmingco2_siteinfo$pred_jmax_logr[i] <- log(vcmax25_warmingco2_siteinfo$mean_jmax_elevated[i]/vcmax25_warmingco2_siteinfo$mean_jmax_ambient[i])
  vcmax25_warmingco2_siteinfo$pred_jmax_coef[i] <- (log(vcmax25_warmingco2_siteinfo$mean_jmax_elevated[i]/vcmax25_warmingco2_siteinfo$mean_jmax_ambient[i]))/(log(vcmax25_warmingco2_siteinfo$co2_e[i]/vcmax25_warmingco2_siteinfo$co2_a[i]))
  vcmax25_warmingco2_siteinfo$pred_jmax25_logr[i] <- log(vcmax25_warmingco2_siteinfo$mean_jmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_jmax25_ambient[i])
  vcmax25_warmingco2_siteinfo$pred_jmax25_coef[i] <- (log(vcmax25_warmingco2_siteinfo$mean_jmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_jmax25_ambient[i]))/(log(vcmax25_warmingco2_siteinfo$co2_e[i]/vcmax25_warmingco2_siteinfo$co2_a[i]))
  
  vcmax25_warmingco2_siteinfo$PPFD[i] <- mean(df1$ppfd,na.rm=TRUE)*1000000 # in umol/m2/s
  vcmax25_warmingco2_siteinfo$vpd[i] <- mean(df1$vpd,na.rm=TRUE)/1000 #in kPa
  vcmax25_warmingco2_siteinfo$temp[i] <- mean(df1$temp,na.rm=TRUE) #in degree celcius
}

smith_prediction <- vcmax25_warmingco2_siteinfo

#remove repeated sites and make site-name consistent
smith_all_removal <- subset(smith_prediction,exp_nam!="NZFACE" &exp_nam!="BioCON" & exp_nam!="ChinaRiceFACE" & exp_nam!="DukeFACE" & exp_nam!="EUROPOPFACE" & exp_nam!="NevadaFACE" & exp_nam!="SwissFACE")
smith_all_removal$exp <- smith_all_removal$exp_nam
smith_all_removal$exp[smith_all_removal$exp=="OakOTC"] <- "mi_c";
smith_all_removal$exp[smith_all_removal$exp=="ORNL"] <- "ornerp_liqui_c";
smith_all_removal$exp[smith_all_removal$exp=="PineOTC"] <- "Christchurch_pr_c";
smith_all_removal$exp[smith_all_removal$exp=="RichmondGH"] <- "australiacotton_c";
smith_all_removal$exp[smith_all_removal$exp=="UIAPine"] <- "ua_otc_c";
smith_all_removal$exp[smith_all_removal$exp=="PHACE"] <- "phace_c";
smith_all_removal$exp[smith_all_removal$exp=="AspenFACE"] <- "aspenface_c";
smith_all_removal$exp[smith_all_removal$exp=="BilyKriz"] <- "bily_kriz_c";
smith_all_removal$exp[smith_all_removal$exp=="Headley"] <- "headley_qp_c";
smith_all_removal$exp[smith_all_removal$exp=="Viesalm"] <- "vielsalm_c"
smith_all_removal$exp <- tolower(smith_all_removal$exp)

vcmax_smith <- smith_all_removal %>% group_by(lon,lat,exp)  %>% summarise(vcmax = mean(pred_vcmax_coef))
jmax_smith <- smith_all_removal %>% group_by(lon,lat,exp)  %>% summarise(jmax = mean(pred_jmax_coef))

######3. Warming from dushan's dataset
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
#see below Table s1: we can assume ambient = 19.7, elevated = 24.9 (basing on growth temperature column)
#https://nph.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fnph.15668&file=nph15668-sup-0001-SupInfo.pdf
df3$Temp_Treatment[df3$Dataset=="Eucalyptus globulus, AU-NSW" & df3$Temp_Treatment=="ambient"] <- 19.7
df3$Temp_Treatment[df3$Dataset=="Eucalyptus globulus, AU-NSW" & df3$Temp_Treatment=="elevated"] <- 24.9

#for "Eucalyptus parramattensis, AU-NSW"
#see below Table s1: we can assume ambient = 15, elevated = 18
#https://nph.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fnph.15668&file=nph15668-sup-0001-SupInfo.pdf
df3$Temp_Treatment[df3$Dataset=="Eucalyptus parramattensis, AU-NSW" & df3$Temp_Treatment=="ambient"] <- 15.1
df3$Temp_Treatment[df3$Dataset=="Eucalyptus parramattensis, AU-NSW" & df3$Temp_Treatment=="elevated"] <- 18.0

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

obs_warming$sitename <- paste("kumarathunge",1:length(obs_warming$lon),sep="")
devtools::load_all("~/yunkepeng/gcme/pmodel/ingestr/")
df_etopo <- ingest(obs_warming,source = "etopo1",dir = "~/data/etopo/" )
obs_warming$elv <- as.numeric(as.data.frame(df_etopo$data))
siteinfo_final <- as.data.frame(obs_warming[,c("sitename","lon","lat","elv","year_start","year_end")])

csvfile <- paste("~/data/Kumarathunge_2020_newphy/forcing_info.csv")
write.csv(siteinfo_final, csvfile, row.names = TRUE)

#for (i in 1:(nrow(siteinfo_final))){
#  df_watch <- ingest(
#    siteinfo  = siteinfo_final[i,],
#    source    = "watch_wfdei",
#    getvars   = c("temp", "prec", "ppfd", "vpd", "patm"), 
#    dir       = "/Volumes/My Passport/data/watch_wfdei/",
#    settings  = list(correct_bias = "worldclim", dir_bias = "/Volumes/My Passport/data/worldclim/"))
#    df_cru <- ingest(
#    siteinfo  = siteinfo_final[i,],
#    source    = "cru",
#    getvars   = "ccov",
#    dir       = "/Volumes/My Passport/data/cru/ts_4.01/")
#  df_co2 <- ingest(
#    siteinfo  = siteinfo_final[i,],
#    source  = "co2_mlo",
#    verbose = FALSE)
#  df_co2_final <- as.data.frame(df_co2$data)
#  df_co2_final2 <- df_co2_final[!(format(df_co2_final$date,"%m") == "02" & format(df_co2_final$date, "%d") == "29"), , drop = FALSE] # make columns consistent
#  co2 <- df_co2_final2$co2
#  ddf_meteo <- as_tibble(cbind(as.data.frame(df_watch$data),as.data.frame(df_cru$data)[,c("ccov_int","ccov")],co2))
#  ddf_meteo$fapar <- 1
#  ddf_meteo$tmax <- ddf_meteo$temp
#  ddf_meteo$tmin <- ddf_meteo$temp
#  ddf_meteo$sitename <- siteinfo_final[i,c("sitename")]
#  ddf_meteo$lon <- siteinfo_final[i,c("lon")]
#  ddf_meteo$lat <-  siteinfo_final[i,c("lat")]
#  ddf_meteo$elv <-  siteinfo_final[i,c("elv")]
#  ddf_meteo$year_start <-  siteinfo_final[i,c("year_start")]
#  ddf_meteo$year_end <-  siteinfo_final[i,c("year_end")]
#  csvfile <- paste("~/data/Kumarathunge_2020_newphy/forcing_warming/",siteinfo_final$sitename[i],".csv",sep = "")
#  write.csv(ddf_meteo, csvfile, row.names = TRUE)
#  print(i)    
#}

forcing_path <- "~/data/Kumarathunge_2020_newphy/forcing_warming/"

obs_warming


for (i in 1:nrow(obs_warming)){
  siteinfo_site <- obs_warming[i,]
  sitename <- siteinfo_site$sitename
  
  df1 <-read.csv(paste(forcing_path,sitename,".csv",sep=""))  # this is for ambient
  df2 <- read.csv(paste(forcing_path,sitename,".csv",sep=""))  # this is for elevated
  
  siteinfo_site$elv <- df1$elv[1]
  
  #Q1 - add warming directly
  df2$temp <- df2$temp+obs_warming$t1[i]
  
  
  modlist1 <- run_pmodel_f_bysite( 
    sitename=sitename, 
    params_siml <- list(#ambient
      spinup             = TRUE,
      spinupyears        = 10,
      recycle            = 1,
      soilmstress        = TRUE,
      tempstress         = TRUE,
      calc_aet_fapar_vpd = FALSE,
      in_ppfd            = TRUE,
      in_netrad          = FALSE,
      outdt              = 1,
      ltre               = FALSE,
      ltne               = FALSE,
      ltrd               = FALSE,
      ltnd               = FALSE,
      lgr3               = TRUE,
      lgn3               = FALSE,
      lgr4               = FALSE,
      firstyeartrend = siteinfo_site$year_start,
      nyeartrend = siteinfo_site$year_end-siteinfo_site$year_start+1), 
    siteinfo_site, 
    df1, 
    df_soiltexture, 
    params_modl = params_modl, 
    makecheck = TRUE)
  
  modlist2 <- run_pmodel_f_bysite( #elevated
    sitename=sitename, 
    params_siml <- list(
      spinup             = TRUE,
      spinupyears        = 10,
      recycle            = 1,
      soilmstress        = TRUE,
      tempstress         = TRUE,
      calc_aet_fapar_vpd = FALSE,
      in_ppfd            = TRUE,
      in_netrad          = FALSE,
      outdt              = 1,
      ltre               = FALSE,
      ltne               = FALSE,
      ltrd               = FALSE,
      ltnd               = FALSE,
      lgr3               = TRUE,
      lgn3               = FALSE,
      lgr4               = FALSE,
      firstyeartrend = siteinfo_site$year_start,
      nyeartrend = siteinfo_site$year_end-siteinfo_site$year_start+1), 
    siteinfo_site, 
    df2, 
    df_soiltexture, 
    params_modl = params_modl, 
    makecheck = TRUE)
  #mean or max?
  mean_vcmax25_ambient <- max(modlist1$vcmax25)*1000000
  mean_vcmax25_elevated <- max(modlist2$vcmax25)*1000000
  mean_vcmax_ambient <- max(modlist1$vcmax)*1000000
  mean_vcmax_elevated <- max(modlist2$vcmax)*1000000
  
  mean_jmax25_ambient <- max(modlist1$jmax25)*1000000
  mean_jmax25_elevated <- max(modlist2$jmax25)*1000000
  mean_jmax_ambient <- max(modlist1$jmax)*1000000
  mean_jmax_elevated <- max(modlist2$jmax)*1000000
  
  obs_warming$mean_vcmax25_ambient[i] <- mean_vcmax25_ambient
  obs_warming$mean_vcmax25_elevated[i] <- mean_vcmax25_elevated
  obs_warming$mean_vcmax_ambient[i] <- mean_vcmax_ambient
  obs_warming$mean_vcmax_elevated[i] <- mean_vcmax_elevated
  
  obs_warming$mean_jmax25_ambient[i] <- mean_jmax25_ambient
  obs_warming$mean_jmax25_elevated[i] <- mean_jmax25_elevated
  obs_warming$mean_jmax_ambient[i] <- mean_jmax_ambient
  obs_warming$mean_jmax_elevated[i] <- mean_jmax_elevated
  
  obs_warming$PPFD[i] <- mean(df1$ppfd,na.rm=TRUE)*1000000 # in umol/m2/s
  obs_warming$vpd[i] <- mean(df1$vpd,na.rm=TRUE)/1000 #in kPa
  obs_warming$temp[i] <- mean(df1$temp,na.rm=TRUE) #in degree celcius
  
  obs_warming$pred_vcmax25_coef[i] <- (log(obs_warming$mean_vcmax25_elevated[i]/obs_warming$mean_vcmax25_ambient[i]))/obs_warming$t1[i]
  obs_warming$pred_jmax25_coef[i] <- (log(obs_warming$mean_jmax25_elevated[i]/obs_warming$mean_jmax25_ambient[i]))/obs_warming$t1[i]
}

vcmax_warming <- obs_warming[,c("lon","lat","exp","pred_vcmax25_coef","pred_jmax25_coef")]
names(vcmax_warming) <- c("lon","lat","exp","vcmax","jmax")

######4. light treatment from Walker 2014

walker <- read.csv("~/data/leaf_traits/Walker/LEAF_PHOTOSYNTHESIS_TRAITS_1224/data/Leaf_Photosynthesis_Traits.csv")

#vcmax
light_vcmax <- subset(walker,Treatment_light!="not applicable") %>% 
  group_by(lon,lat,z,Year,species,Treatment_light)  %>% summarise(vcmax_mean = mean(Vcmax),
                                                                  jmax_mean = mean(Jmax))

siteinfo <- unique(light_vcmax[,c("lon","lat","z","Year")])

siteinfo$Year[siteinfo$Year=="2002/2003"] <-2002 

siteinfo$year_start <- siteinfo$Year
siteinfo$year_end <- siteinfo$Year

#two literatures did not report measurement year
#https://www.sciencedirect.com/science/article/pii/S037811270500085X
#https://onlinelibrary.wiley.com/doi/full/10.1111/j.1442-9993.2005.01486.x
#both published in 2005, so we assumed that it was measured within 5 years until publication. 

siteinfo$year_start[siteinfo$Year=="not reported"] <- 2001
siteinfo$year_end[siteinfo$Year=="not reported"] <- 2005
siteinfo$year_start<- as.numeric(siteinfo$year_start)
siteinfo$year_end<- as.numeric(siteinfo$year_end)

siteinfo$elv <- siteinfo$z
siteinfo$sitename <- paste("light",1:nrow(siteinfo),sep="")

siteinfo_final <- as.data.frame(siteinfo[,c("sitename","lon","lat","elv","year_start","year_end")])

csvfile <- paste("~/data/leaf_traits/Walker/forcing_info.csv")
write.csv(siteinfo_final, csvfile, row.names = TRUE)

#for (i in 1:(nrow(siteinfo_final))){
#  df_watch <- ingest(
#    siteinfo  = siteinfo_final[i,],
#    source    = "watch_wfdei",
#    getvars   = c("temp", "prec", "ppfd", "vpd", "patm"), 
#    dir       = "/Volumes/My Passport/data/watch_wfdei/",
#    settings  = list(correct_bias = "worldclim", dir_bias = "/Volumes/My Passport/data/worldclim/"))
  
#  df_cru <- ingest(
#    siteinfo  = siteinfo_final[i,],
#    source    = "cru",
#    getvars   = "ccov",
#    dir       = "/Volumes/My Passport/data/cru/ts_4.01/")
#  df_co2 <- ingest(
#    siteinfo  = siteinfo_final[i,],
#    source  = "co2_mlo",
#    verbose = FALSE)
#  df_co2_final <- as.data.frame(df_co2$data)
#  df_co2_final2 <- df_co2_final[!(format(df_co2_final$date,"%m") == "02" & format(df_co2_final$date, "%d") == "29"), , drop = FALSE] # make columns consistent
#  co2 <- df_co2_final2$co2
#  ddf_meteo <- as_tibble(cbind(as.data.frame(df_watch$data),as.data.frame(df_cru$data)[,c("ccov_int","ccov")],co2))
#  ddf_meteo$fapar <- 1
#  ddf_meteo$tmax <- ddf_meteo$temp
#  ddf_meteo$tmin <- ddf_meteo$temp
#  ddf_meteo$sitename <- siteinfo_final[i,c("sitename")]
#  ddf_meteo$lon <- siteinfo_final[i,c("lon")]
#  ddf_meteo$lat <-  siteinfo_final[i,c("lat")]
#  ddf_meteo$elv <-  siteinfo_final[i,c("elv")]
#  ddf_meteo$year_start <-  siteinfo_final[i,c("year_start")]
#  ddf_meteo$year_end <-  siteinfo_final[i,c("year_end")]
#  csvfile <- paste("~/data/leaf_traits/Walker/forcing_light/",siteinfo_final$sitename[i],".csv",sep = "")
#  write.csv(ddf_meteo, csvfile, row.names = TRUE)
#  print(i)    
#}

forcing_path <- "~/data/leaf_traits/Walker/forcing_light/"

for (i in 1:nrow(siteinfo_final)){
  siteinfo_final_one <- siteinfo_final[i,c("sitename","lon","lat","elv","year_start","year_end")]
  sitename <- siteinfo_final_one$sitename
  
  df1 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for ambient
  df2 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for elevated
  
  #Q1 - add warming directly
  df2$ppfd <- df1$ppfd*1.5
  
  #Q2 - if having co2 then replacing it directly - if not then using co2 at measurement year
  modlist1 <- run_pmodel_f_bysite( 
    sitename=sitename, 
    params_siml <- list(#ambient
      spinup             = TRUE,
      spinupyears        = 10,
      recycle            = 1,
      soilmstress        = TRUE,
      tempstress         = TRUE,
      calc_aet_fapar_vpd = FALSE,
      in_ppfd            = TRUE,
      in_netrad          = FALSE,
      outdt              = 1,
      ltre               = FALSE,
      ltne               = FALSE,
      ltrd               = FALSE,
      ltnd               = FALSE,
      lgr3               = TRUE,
      lgn3               = FALSE,
      lgr4               = FALSE,
      firstyeartrend = siteinfo_final_one$year_start,
      nyeartrend = siteinfo_final_one$year_end-siteinfo_final_one$year_start+1), 
    siteinfo_final_one, 
    df1, 
    df_soiltexture, 
    params_modl = params_modl, 
    makecheck = TRUE)
  
  modlist2 <- run_pmodel_f_bysite( #elevated
    sitename=sitename, 
    params_siml <- list(
      spinup             = TRUE,
      spinupyears        = 10,
      recycle            = 1,
      soilmstress        = TRUE,
      tempstress         = TRUE,
      calc_aet_fapar_vpd = FALSE,
      in_ppfd            = TRUE,
      in_netrad          = FALSE,
      outdt              = 1,
      ltre               = FALSE,
      ltne               = FALSE,
      ltrd               = FALSE,
      ltnd               = FALSE,
      lgr3               = TRUE,
      lgn3               = FALSE,
      lgr4               = FALSE,
      firstyeartrend = siteinfo_final_one$year_start,
      nyeartrend = siteinfo_final_one$year_end-siteinfo_final_one$year_start+1), 
    siteinfo_final_one, 
    df2, 
    df_soiltexture, 
    params_modl = params_modl, 
    makecheck = TRUE)
  #mean or max?
  mean_vcmax25_ambient <- max(modlist1$vcmax25)*1000000
  mean_vcmax25_elevated <- max(modlist2$vcmax25)*1000000
  mean_vcmax_ambient <- max(modlist1$vcmax)*1000000
  mean_vcmax_elevated <- max(modlist2$vcmax)*1000000
  
  mean_jmax25_ambient <- max(modlist1$jmax25)*1000000
  mean_jmax25_elevated <- max(modlist2$jmax25)*1000000
  mean_jmax_ambient <- max(modlist1$jmax)*1000000
  mean_jmax_elevated <- max(modlist2$jmax)*1000000
  
  siteinfo_final$mean_vcmax25_ambient[i] <- mean_vcmax25_ambient
  siteinfo_final$mean_vcmax25_elevated[i] <- mean_vcmax25_elevated
  siteinfo_final$mean_vcmax_ambient[i] <- mean_vcmax_ambient
  siteinfo_final$mean_vcmax_elevated[i] <- mean_vcmax_elevated
  siteinfo_final$pred_vcmax_logr[i] <- log(siteinfo_final$mean_vcmax_elevated[i]/siteinfo_final$mean_vcmax_ambient[i])
  siteinfo_final$pred_vcmax25_logr[i] <- log(siteinfo_final$mean_vcmax25_elevated[i]/siteinfo_final$mean_vcmax25_ambient[i])
  
  siteinfo_final$mean_jmax25_ambient[i] <- mean_jmax25_ambient
  siteinfo_final$mean_jmax25_elevated[i] <- mean_jmax25_elevated
  siteinfo_final$mean_jmax_ambient[i] <- mean_jmax_ambient
  siteinfo_final$mean_jmax_elevated[i] <- mean_jmax_elevated
  siteinfo_final$pred_jmax_logr[i] <- log(siteinfo_final$mean_jmax_elevated[i]/siteinfo_final$mean_jmax_ambient[i])
  siteinfo_final$pred_jmax25_logr[i] <- log(siteinfo_final$mean_jmax25_elevated[i]/siteinfo_final$mean_jmax25_ambient[i])
  
  siteinfo_final$PPFD[i] <- mean(df1$ppfd,na.rm=TRUE)*1000000 # in umol/m2/s
  siteinfo_final$vpd[i] <- mean(df1$vpd,na.rm=TRUE)/1000 #in kPa
  siteinfo_final$temp[i] <- mean(df1$temp,na.rm=TRUE) #in degree celcius
}

vcmax_light <- siteinfo_final[,c("lon","lat","sitename","pred_vcmax_logr","pred_jmax_logr")]
names(vcmax_light) <- c("lon","lat","exp","vcmax","jmax")

#now, combine all dataframe
kevin_co2 <- merge(vcmax_kevin,jmax_kevin,by=c("lon","lat","exp"))
kevin_co2$ref <- "GCME_kevin"
smith_co2 <- merge(vcmax_smith,jmax_smith,by=c("lon","lat","exp"))
smith_co2$ref <- "Smith_GCB"
co2_prediction <- rbind(kevin_co2,smith_co2)
co2_prediction$treatment <- "co2"
co2_prediction$comments <- "vcmax and jmax are sensitivity coefficients = log(vcmax-ele/vcmax-amb)/log(co2-e/co2-a)"

vcmax_warming$ref <- "Kumarathunge New Phytol"
vcmax_warming$treatment <- "warming"
vcmax_warming$comments <- "vcmax and jmax are sensitivity coefficients of Vcmax25 and Jmax25 = log(vcmax25-ele/vcmax25-amb)/ delta T"

vcmax_light$ref <- "Walker 2014 dataset"
vcmax_light$treatment <- "light"
vcmax_light$comments <- "vcmax and jmax are sensitivity coefficients = log(vcmax25-ele/vcmax25-amb)"

final_prediction <- rbind(co2_prediction,vcmax_warming,vcmax_light)

final_prediction$jmax_vcmax <- final_prediction$jmax - final_prediction$vcmax

final_prediction <- final_prediction[,c("exp","lon","lat","vcmax","jmax","jmax_vcmax","treatment","ref","comments")]


###############Add one more site: biforface_c
###according to MESI_2022.xlsx (the most recent data after above vcmax data version) - one more exp was newly added for vcmax data
kevin_othervars <- read.csv("~/data/gcme/kevin_20220222/MESI_2022.csv")
unique(subset(kevin_othervars,exp=="biforface_c" & response=="vcmax")[,c("lon","lat","elevation","start_year","sampling_year","c_c","c_t")])
#add them here
siteinfo_site <- tibble(sitename="kevin61",lon=-2.3,lat=52.8,elv=108,year_start=2012,year_end=2016) # forcing data was saved in ~/data/gcme/kevin/forcing/climate/kevin61.csv
sitename <- "biforface_c"

df1 <- read.csv("~/data/gcme/kevin/forcing/climate/kevin61.csv") # this is for ambient
df2 <- read.csv("~/data/gcme/kevin/forcing/climate/kevin61.csv") # this is for elevated

df1$co2 <- 408 # set to ambient
df2$co2 <- 558 # set to elevated 

df_soiltexture <- bind_rows(
  top    = tibble(layer = "top",    fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1),
  bottom = tibble(layer = "bottom", fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1))
params_modl <- list(
  kphio           = 0.09423773,
  soilm_par_a     = 0.33349283,
  soilm_par_b     = 1.45602286,
  tau_acclim_tempstress = 10,
  par_shape_tempstress  = 0.0)

modlist1 <- run_pmodel_f_bysite( 
  sitename=sitename, 
  params_siml <- list(#ambient
    spinup             = TRUE,
    spinupyears        = 10,
    recycle            = 1,
    soilmstress        = TRUE,
    tempstress         = TRUE,
    calc_aet_fapar_vpd = FALSE,
    in_ppfd            = TRUE,
    in_netrad          = FALSE,
    outdt              = 1,
    ltre               = FALSE,
    ltne               = FALSE,
    ltrd               = FALSE,
    ltnd               = FALSE,
    lgr3               = TRUE,
    lgn3               = FALSE,
    lgr4               = FALSE,
    firstyeartrend = siteinfo_site$year_start,
    nyeartrend = siteinfo_site$year_end-siteinfo_site$year_start+1), 
  siteinfo_site, 
  df1, 
  df_soiltexture, 
  params_modl = params_modl, 
  makecheck = TRUE)

modlist2 <- run_pmodel_f_bysite( #elevated
  sitename=sitename, 
  params_siml <- list(
    spinup             = TRUE,
    spinupyears        = 10,
    recycle            = 1,
    soilmstress        = TRUE,
    tempstress         = TRUE,
    calc_aet_fapar_vpd = FALSE,
    in_ppfd            = TRUE,
    in_netrad          = FALSE,
    outdt              = 1,
    ltre               = FALSE,
    ltne               = FALSE,
    ltrd               = FALSE,
    ltnd               = FALSE,
    lgr3               = TRUE,
    lgn3               = FALSE,
    lgr4               = FALSE,
    firstyeartrend = siteinfo_site$year_start,
    nyeartrend = siteinfo_site$year_end-siteinfo_site$year_start+1), 
  siteinfo_site, 
  df2, 
  df_soiltexture, 
  params_modl = params_modl, 
  makecheck = TRUE)
#mean or max?
mean_vcmax_ambient <- max(modlist1$vcmax)*1000000
mean_vcmax_elevated <- max(modlist2$vcmax)*1000000
mean_jmax_ambient <- max(modlist1$jmax)*1000000
mean_jmax_elevated <- max(modlist2$jmax)*1000000

pred_vcmax_coef <- (log(mean_vcmax_elevated/mean_vcmax_ambient))/(log(558/408))
pred_jmax_coef <- (log(mean_jmax_elevated/mean_jmax_ambient))/(log(558/408))

biforface_c <- tibble( exp="biforface_c",lon=-2.3,lat=52.8,vcmax=pred_vcmax_coef,jmax=pred_jmax_coef,jmax_vcmax=pred_jmax_coef-pred_vcmax_coef,
        treatment="co2",ref="GCME_kevin",comments="vcmax and jmax are sensitivity coefficients = log(vcmax-ele/vcmax-amb)/log(co2-e/co2-a)")

final_prediction2 <- rbind(final_prediction[,c("exp","lon","lat","vcmax","jmax","jmax_vcmax","treatment","ref","comments")],biforface_c)

csvfile <- paste("~/data/gcme/prediction/prediction.csv")
write.csv(final_prediction2, csvfile, row.names = TRUE)
