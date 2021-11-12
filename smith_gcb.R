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


df <- read.csv("~/data/smith_keenan_gcb/gcb_data/Growth_chamber_dataset.csv")
df[df == -9999] <- NA
summary(df)
#firstly,  output measured vcmax25 at different temperature
#unit should be in K
calc_vcmax_to_vcmax25 <- function(Tleaf, Tg){
  R_constant <- 8.314 # universal gas constant J/mol/K
  Ha <- 71513 #activitation energy J/mol
  delta_s <- 641.64 #entropy term J/mol/k
  d_ha <- 200000 #deactivation energy J/mol
  ftemp <- exp( Ha * (Tleaf - 298.15) / (Tg * R_constant * 298.15) )
  fva <- (1.0 + exp( (298.15 * delta_s - d_ha)/(R_constant * 298.15) ) ) / (1.0 + exp( (Tleaf * delta_s - d_ha)/(R_constant * Tleaf) ) )
  vcmax25_conversion <- ftemp * fva
  return(vcmax25_conversion)
}
calc_jmax_to_jmax25 <- function(Tleaf, Tg){
  R_constant <- 8.314 # universal gas constant J/mol/K
  Ha <- 49884 #activitation energy J/mol
  delta_s <- 642.2 #entropy term J/mol/k
  d_ha <- 200000 #deactivation energy J/mol
  ftemp <- exp( Ha * (Tleaf - 298.15) / (Tg * R_constant * 298.15) )
  fva <- (1.0 + exp( (298.15 * delta_s - d_ha)/(R_constant * 298.15) ) ) / (1.0 + exp( (Tleaf * delta_s - d_ha)/(R_constant * Tleaf) ) )
  jmax25_conversion <- ftemp * fva
  return(jmax25_conversion)
}

df$vcmax25_14obs <- df$Vcmax14_obs/calc_vcmax_to_vcmax25(df$Tleaf14_Vcmax_obs+273.15,df$GrowthT+273.15)
df$vcmax25_23obs <- df$Vcmax23_obs/calc_vcmax_to_vcmax25(df$Tleaf23_Vcmax_obs+273.15,df$GrowthT+273.15)
df$vcmax25_32obs <- df$Vcmax32_obs/calc_vcmax_to_vcmax25(df$Tleaf32_Vcmax_obs+273.15,df$GrowthT+273.15)
df$vcmax25_41obs <- df$Vcmax41_obs/calc_vcmax_to_vcmax25(df$Tleaf41_Vcmax_obs+273.15,df$GrowthT+273.15)
df$vcmax25_50obs <- df$Vcmax50_obs/calc_vcmax_to_vcmax25(df$Tleaf50_Vcmax_obs+273.15,df$GrowthT+273.15)

df$jmax25_14obs <- df$Jmax14_obs/calc_jmax_to_jmax25(df$Tleaf14_Jmax_obs+273.15,df$GrowthT+273.15)
df$jmax25_23obs <- df$Jmax23_obs/calc_jmax_to_jmax25(df$Tleaf23_Jmax_obs+273.15,df$GrowthT+273.15)
df$jmax25_32obs <- df$Jmax32_obs/calc_jmax_to_jmax25(df$Tleaf32_Jmax_obs+273.15,df$GrowthT+273.15)
df$jmax25_41obs <- df$Jmax41_obs/calc_jmax_to_jmax25(df$Tleaf41_Jmax_obs+273.15,df$GrowthT+273.15)
df$jmax25_50obs <- df$Jmax50_obs/calc_jmax_to_jmax25(df$Tleaf50_Jmax_obs+273.15,df$GrowthT+273.15)

#prepare forcing data
# first - collect climate forcing
#devtools::load_all("~/yunkepeng/gcme/pmodel/ingestr/")
sitename <- "smith_keenan"
lon <- -86.92
lat <- 40.82
elv <- 187
year_start <- 2013 #2013-212: 212 day of 2013 --> /Users/yunpeng/data/smith_keenan_gcb/gcb_data/raw_ACi/agerardii_gt25C_r1_mt14c.csv
year_end <- 2015 # 2015-50: 50 day of 2015 --> /Users/yunpeng/data/smith_keenan_gcb/gcb_data/raw_ACi/zmays_b_gt35C_r3_mt50c.csv
siteinfo_final <- data.frame(sitename,lon,lat,elv,year_start,year_end)
#df_watch <- ingest(siteinfo  = siteinfo_final,source    = "watch_wfdei",getvars= c("temp", "prec", "ppfd", "vpd", "patm"), dir= "/Volumes/My Passport/data/watch_wfdei/",settings= list(correct_bias = "worldclim", dir_bias = "/Volumes/My Passport/data/worldclim/"))
#df_cru <- ingest(siteinfo  = siteinfo_final,source = "cru",getvars   = "ccov",dir= "/Volumes/My Passport/data/cru/ts_4.01/")
#df_co2 <- ingest(siteinfo  = siteinfo_final,source  = "co2_mlo",verbose = FALSE)
#df_co2_final <- as.data.frame(df_co2$data)
#df_co2_final2 <- df_co2_final[!(format(df_co2_final$date,"%m") == "02" & format(df_co2_final$date, "%d") == "29"), , drop = FALSE] # make columns consistent
#co2 <- df_co2_final2$co2
#ddf_meteo <- as_tibble(cbind(as.data.frame(df_watch$data),as.data.frame(df_cru$data)[,c("ccov_int","ccov")],co2))
#ddf_meteo$fapar <- 1
#ddf_meteo$tmax <- ddf_meteo$temp
#ddf_meteo$tmin <- ddf_meteo$temp
#ddf_meteo$sitename <- sitename
#ddf_meteo$lon <- lon
#ddf_meteo$lat <-  lat
#ddf_meteo$elv <-  elv
#ddf_meteo$year_start <-  year_start
#ddf_meteo$year_end <- year_end
#csvfile <- paste("/Users/yunpeng/data/smith_keenan_gcb/gcb_data/smith_keenan_forcing.csv",sep = "")
#write.csv(ddf_meteo, csvfile, row.names = TRUE)

#secondly, modify climate variables and run p-model!
forcing_site <- read.csv("/Users/yunpeng/data/smith_keenan_gcb/gcb_data/smith_keenan_forcing.csv")
df_soiltexture <- bind_rows(
  top    = tibble(layer = "top",    fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1),
  bottom = tibble(layer = "bottom", fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1))
params_modl <- list(
  kphio           = 0.09423773,
  soilm_par_a     = 0.33349283,
  soilm_par_b     = 1.45602286,
  tau_acclim_tempstress = 10,
  par_shape_tempstress  = 0.0)

#convert to vpd --> smith has set vpd just depend on temp and 50% RH (relative humidity), and also patm - but  this is very small
##' Convert specific humidity to relative humidity
##see original formula in: https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity
calc_vpd_smith <- function(patm, temp){ #patm in pa and temp in degree celcius
  press <- patm/100 #convert Pa to hPa
  rh <- 0.50
  es <-  6.112 * exp((17.67 * temp)/(temp + 243.5))
  qair <- 0.622 / ((press/rh/es)-0.378)
  kR  = 8.31446262  # universal gas constant, J/mol/K (Allen, 1973)
  kMv = 18.02                 # molecular weight of water vapor, g/mol (Tsilingiris, 2008)
  kMa = 28.963                #molecular weight of dry air, g/mol (Tsilingiris, 2008) XXX this was in SPLASH (WITH 1E-3 IN EQUATION) XXX
  wair <- qair / ( 1 - qair )
  rv <- kR / kMv
  rd <- kR / kMa
  ea <- (patm * wair * rv / (rd + wair * rv))/1000 #ea in Pa
  vpd <- 0.611*exp(17.27*(temp)/((temp)+237.3)) - ea         #vpd in KPa -> month/daily
  vpd[vpd < 0] <- NA
  vpd <- vpd*1000 #convert to pa
  return(vpd)
}


#option 1:  calculate vcmax25 and jmax25 at given temperature (Tgrowth) and vpd (dependent on Tgrowth + 50% humidity)
forcing_site_pred1 <- forcing_site
for (i in 1:nrow(df)){
  forcing_site_pred1$temp <- df$GrowthT[i]
  forcing_site_pred1$vpd <- calc_vpd_smith(forcing_site_pred1$patm,df$GrowthT[i])
  modlist <- run_pmodel_f_bysite( 
    sitename="smith_keenan", 
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
      firstyeartrend = 2013,
      nyeartrend = 3), 
    siteinfo_final, 
    forcing_site_pred1, 
    df_soiltexture, 
    params_modl = params_modl, 
    makecheck = TRUE)
  df$pred_vcmax[i] <- mean(modlist$vcmax*1000000,na.rm=TRUE)
  df$pred_jmax[i] <- mean(modlist$jmax*1000000,na.rm=TRUE)
  df$pred_vcmax25[i] <- mean(modlist$vcmax25*1000000,na.rm=TRUE)
  df$pred_jmax25[i] <- mean(modlist$jmax25*1000000,na.rm=TRUE)
}

#option 2:  calculate vcmax25 and jmax25 at given temperature (Tgrowth) and vpd (dependent on Tgrowth + 50% humidity)
#AND set below constant value
ppfd <- 800 # par = 800
ca <- 400 #ppm
forcing_site_pred2 <- forcing_site
for (i in 1:nrow(df)){
  forcing_site_pred2$temp <- df$GrowthT[i]
  forcing_site_pred2$vpd <- calc_vpd_smith(forcing_site_pred2$patm,df$GrowthT[i])
  forcing_site_pred2$ppfd <- 800/1000000
  forcing_site_pred2$co2 <- 400
  modlist <- run_pmodel_f_bysite( 
    sitename="smith_keenan", 
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
      firstyeartrend = 2013,
      nyeartrend = 3), 
    siteinfo_final, 
    forcing_site_pred2, 
    df_soiltexture, 
    params_modl = params_modl, 
    makecheck = TRUE)
  df$pred_vcmax_smith[i] <- mean(modlist$vcmax*1000000,na.rm=TRUE)
  df$pred_jmax_smith[i] <- mean(modlist$jmax*1000000,na.rm=TRUE)
  df$pred_vcmax25_smith[i] <- mean(modlist$vcmax25*1000000,na.rm=TRUE)
  df$pred_jmax25_smith[i] <- mean(modlist$jmax25*1000000,na.rm=TRUE)
}
#df$GrowthT # 15,20,25,30,35
#tleaf #14,23,32,41,50
names(df)
ggplot(df, aes(x=GrowthT)) +geom_point(aes(y=vcmax25_14obs),color="black")+geom_smooth(aes(y=vcmax25_14obs),color="black")+geom_point(aes(y=pred_vcmax25),color="red")+geom_smooth(aes(y=pred_vcmax25),color="red")
ggplot(df, aes(x=GrowthT)) +geom_point(aes(y=vcmax25_23obs),color="black")+geom_smooth(aes(y=vcmax25_23obs),color="black")+geom_point(aes(y=pred_vcmax25),color="red")+geom_smooth(aes(y=pred_vcmax25),color="red")
ggplot(df, aes(x=GrowthT)) +geom_point(aes(y=vcmax25_32obs),color="black")+geom_smooth(aes(y=vcmax25_32obs),color="black")+geom_point(aes(y=pred_vcmax25),color="red")+geom_smooth(aes(y=pred_vcmax25),color="red")
ggplot(df, aes(x=GrowthT)) +geom_point(aes(y=vcmax25_41obs),color="black")+geom_smooth(aes(y=vcmax25_41obs),color="black")+geom_point(aes(y=pred_vcmax25),color="red")+geom_smooth(aes(y=pred_vcmax25),color="red")
ggplot(df, aes(x=GrowthT)) +geom_point(aes(y=vcmax25_50obs),color="black")+geom_smooth(aes(y=vcmax25_50obs),color="black")+geom_point(aes(y=pred_vcmax25),color="red")+geom_smooth(aes(y=pred_vcmax25),color="red")

ggplot(df, aes(x=GrowthT))+geom_point(aes(y=pred_vcmax25))+geom_smooth(aes(y=pred_vcmax25,colour="prediction"))+
  geom_smooth(aes(y=vcmax25_14obs,color="Tleaf = 14"),se=F)+geom_smooth(aes(y=vcmax25_23obs,color="Tleaf = 23"),se=F)+geom_smooth(aes(y=vcmax25_32obs,color="Tleaf = 32"),se=F)+
  geom_smooth(aes(y=vcmax25_41obs,color="Tleaf = 41"),se=F)+geom_smooth(aes(y=vcmax25_50obs,color="Tleaf = 50"),se=F)+
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 15, vjust=0.3))

ggplot(df, aes(x=GrowthT))+geom_point(aes(y=pred_vcmax25))+geom_smooth(aes(y=pred_vcmax25,colour="prediction"))+
  geom_smooth(aes(y=vcmax25_23obs,color="Tleaf = 23 using kinetic"),se=F)+
  geom_smooth(aes(y=vcmax25_23obs_normal,color="Tleaf = 23 using arrhenius"),se=F)+
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 15, vjust=0.3))

#vcmax?
ggplot(df, aes(x=GrowthT))+geom_point(aes(y=pred_vcmax))+geom_smooth(aes(y=pred_vcmax,colour="prediction"))+
  geom_smooth(aes(y=Vcmax14_obs,color="Tleaf = 14"))+geom_smooth(aes(y=Vcmax23_obs,color="Tleaf = 23"))+geom_smooth(aes(y=Vcmax32_obs,color="Tleaf = 32"))+
  geom_smooth(aes(y=Vcmax41_obs,color="Tleaf = 41"))+geom_smooth(aes(y=Vcmax50_obs,color="Tleaf = 50"))+
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 15, vjust=0.3))

#another: just using Arrhenius design
df$vcmax25_14obs_normal <- df$Vcmax14_obs*exp((65330/8.314)*((1/(df$Tleaf14_Vcmax_obs+273.15))-(1/298.15)))
df$vcmax25_23obs_normal <- df$Vcmax23_obs*exp((65330/8.314)*((1/(df$Tleaf23_Vcmax_obs+273.15))-(1/298.15)))
df$vcmax25_32obs_normal <- df$Vcmax32_obs*exp((65330/8.314)*((1/(df$Tleaf32_Vcmax_obs+273.15))-(1/298.15)))
df$vcmax25_41obs_normal <- df$Vcmax41_obs*exp((65330/8.314)*((1/(df$Tleaf41_Vcmax_obs+273.15))-(1/298.15)))
df$vcmax25_50obs_normal <- df$Vcmax50_obs*exp((65330/8.314)*((1/(df$Tleaf50_Vcmax_obs+273.15))-(1/298.15)))

a1 <- as.data.frame(cbind(df$vcmax25_14obs_normal,df$vcmax25_14obs,df$Tleaf14_Vcmax_obs,df$Vcmax14_obs,df$GrowthT,df$pred_vcmax,df$pred_vcmax25,df$pred_vcmax25_smith))
a2 <- as.data.frame(cbind(df$vcmax25_23obs_normal,df$vcmax25_23obs,df$Tleaf23_Vcmax_obs,df$Vcmax23_obs,df$GrowthT,df$pred_vcmax,df$pred_vcmax25,df$pred_vcmax25_smith))
a3 <- as.data.frame(cbind(df$vcmax25_32obs_normal,df$vcmax25_32obs,df$Tleaf32_Vcmax_obs,df$Vcmax32_obs,df$GrowthT,df$pred_vcmax,df$pred_vcmax25,df$pred_vcmax25_smith))
a4 <- as.data.frame(cbind(df$vcmax25_41obs_normal,df$vcmax25_41obs,df$Tleaf41_Vcmax_obs,df$Vcmax41_obs,df$GrowthT,df$pred_vcmax,df$pred_vcmax25,df$pred_vcmax25_smith))
a5 <- as.data.frame(cbind(df$vcmax25_50obs_normal,df$vcmax25_50obs,df$Tleaf50_Vcmax_obs,df$Vcmax50_obs,df$GrowthT,df$pred_vcmax,df$pred_vcmax25,df$pred_vcmax25_smith))

a_final <- dplyr::bind_rows(a1,a2,a3,a4,a5) 
names(a_final) <- c("vcmax25_Arrhenius","vcmax25","Tleaf","vcmax","Tgrowth","pred_vcmax","pred_vcmax25","pred_vcmax25_smith")
a_final$Tgrowth_type <- as.character(a_final$Tgrowth)
ggplot(a_final)+geom_smooth(aes(x=Tleaf,y=vcmax,color=Tgrowth_type),se=F)
ggplot(a_final)+geom_smooth(aes(x=Tleaf,y=vcmax25,color=Tgrowth_type),se=F)

#only subset Tleaf near equal to Tgrowth
a_final_strict <- subset(a_final,Tgrowth - Tleaf >= -3 & Tgrowth - Tleaf <=3)
dim(a_final_strict)
#vcmax25 using peaked function
ggplot(a_final_strict)+geom_smooth(aes(x=Tgrowth,y=vcmax25))+geom_point(aes(x=Tgrowth,y=vcmax25))+
  geom_smooth(aes(x=Tgrowth,y=pred_vcmax25),color="red")

#vcmax25 using normal function
ggplot(a_final_strict)+geom_smooth(aes(x=Tgrowth,y=vcmax25_Arrhenius))+geom_point(aes(x=Tgrowth,y=vcmax25))+
  geom_smooth(aes(x=Tgrowth,y=pred_vcmax25),color="red")

ggplot(a_final_strict)+geom_smooth(aes(x=Tgrowth,y=vcmax25_Arrhenius))+
  geom_smooth(aes(x=Tgrowth,y=vcmax25),color="red")

ggplot(a_final_strict)+geom_smooth(aes(x=Tgrowth,y=vcmax))+geom_point(aes(x=Tgrowth,y=vcmax))+
  geom_smooth(aes(x=Tgrowth,y=pred_vcmax),color="red")

#now, jmax25 response
ggplot(df, aes(x=GrowthT))+geom_point(aes(y=pred_jmax25))+geom_smooth(aes(y=pred_jmax25,colour="prediction"))+
  geom_smooth(aes(y=jmax25_14obs,color="Tleaf = 14"),se=F)+geom_smooth(aes(y=jmax25_23obs,color="Tleaf = 23"),se=F)+geom_smooth(aes(y=jmax25_32obs,color="Tleaf = 32"),se=F)+
  geom_smooth(aes(y=jmax25_41obs,color="Tleaf = 41"),se=F)+geom_smooth(aes(y=jmax25_50obs,color="Tleaf = 50"),se=F)+
  theme_bw() +
  theme(axis.title.x = element_text(size = 15, vjust=-.2)) +
  theme(axis.title.y = element_text(size = 15, vjust=0.3))

#another: just using Arrhenius design
df$jmax25_14obs_normal <- df$Jmax14_obs*exp((43540/8.314)*((1/(df$Tleaf14_Jmax_obs+273.15))-(1/298.15)))
df$jmax25_23obs_normal <- df$Jmax23_obs*exp((43540/8.314)*((1/(df$Tleaf23_Jmax_obs+273.15))-(1/298.15)))
df$jmax25_32obs_normal <- df$Jmax32_obs*exp((43540/8.314)*((1/(df$Tleaf32_Jmax_obs+273.15))-(1/298.15)))
df$jmax25_41obs_normal <- df$Jmax41_obs*exp((43540/8.314)*((1/(df$Tleaf41_Jmax_obs+273.15))-(1/298.15)))
df$jmax25_50obs_normal <- df$Jmax50_obs*exp((43540/8.314)*((1/(df$Tleaf50_Jmax_obs+273.15))-(1/298.15)))

b1 <- as.data.frame(cbind(df$jmax25_14obs_normal,df$jmax25_14obs,df$Tleaf14_Jmax_obs,df$Jmax14_obs,df$GrowthT,df$pred_jmax,df$pred_jmax25,df$pred_jmax25_smith))
b2 <- as.data.frame(cbind(df$jmax25_23obs_normal,df$jmax25_23obs,df$Tleaf23_Jmax_obs,df$Jmax23_obs,df$GrowthT,df$pred_jmax,df$pred_jmax25,df$pred_jmax25_smith))
b3 <- as.data.frame(cbind(df$jmax25_32obs_normal,df$jmax25_32obs,df$Tleaf32_Jmax_obs,df$Jmax32_obs,df$GrowthT,df$pred_jmax,df$pred_jmax25,df$pred_jmax25_smith))
b4 <- as.data.frame(cbind(df$jmax25_41obs_normal,df$jmax25_41obs,df$Tleaf41_Jmax_obs,df$Jmax41_obs,df$GrowthT,df$pred_jmax,df$pred_jmax25,df$pred_jmax25_smith))
b5 <- as.data.frame(cbind(df$jmax25_50obs_normal,df$jmax25_50obs,df$Tleaf50_Jmax_obs,df$Jmax50_obs,df$GrowthT,df$pred_jmax,df$pred_jmax25,df$pred_jmax25_smith))

b_final <- dplyr::bind_rows(b1,b2,b3,b4,b5) 
names(b_final) <- c("jmax25_Arrhenius","jmax25","Tleaf","jmax","Tgrowth","pred_jmax","pred_jmax25","pred_jmax25_smith")
b_final$Tgrowth_type <- as.character(b_final$Tgrowth)
ggplot(b_final)+geom_smooth(aes(x=Tleaf,y=jmax,color=Tgrowth_type),se=F)
ggplot(b_final)+geom_smooth(aes(x=Tleaf,y=jmax25,color=Tgrowth_type),se=F)

#only subset Tleaf near equal to Tgrowth
b_final_strict <- subset(b_final,Tgrowth - Tleaf >= -3 & Tgrowth - Tleaf <=3)
dim(b_final_strict)
#jmax25 using peaked function
ggplot(b_final_strict)+geom_smooth(aes(x=Tgrowth,y=jmax25))+geom_point(aes(x=Tgrowth,y=jmax25))+
  geom_smooth(aes(x=Tgrowth,y=pred_jmax25),color="red")

#jmax25 using normal function
ggplot(b_final_strict)+geom_smooth(aes(x=Tgrowth,y=jmax25_Arrhenius))+geom_point(aes(x=Tgrowth,y=jmax25))+
  geom_smooth(aes(x=Tgrowth,y=pred_jmax25),color="red")

ggplot(b_final_strict)+geom_smooth(aes(x=Tgrowth,y=jmax25_Arrhenius))+
  geom_smooth(aes(x=Tgrowth,y=jmax25),color="red")

ggplot(b_final_strict)+geom_smooth(aes(x=Tgrowth,y=jmax))+geom_point(aes(x=Tgrowth,y=jmax))+
  geom_smooth(aes(x=Tgrowth,y=pred_jmax),color="red")

#leaf N per area
df$narea <- 10*df$Nconc_P/df$Larea_P # mg/cm2 --> g/m2 
df$LMA <- 10000*df$Lmass_P/df$Larea_P # g/cm2 --> g/m2

ggplot(df, aes(x=GrowthT))+geom_point(aes(y=narea))+geom_smooth(aes(y=narea,colour="fitted"))
ggplot(df, aes(x=GrowthT))+geom_point(aes(y=LMA))+geom_smooth(aes(y=LMA,colour="fitted"))


