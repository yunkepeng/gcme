#read all csv, as outputted from vcmax.R
library(readr)
library(dplyr)
library(metafor)  
library(ggplot2)
library(stringr)
Asat_all_final<- read.csv("/Users/yunpeng/data/gcme/yunke_coord/Asat_all_df.csv")

Amax_all_final <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/Amax_all_df.csv")

vcmax_all_final <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/vcmax_all_df.csv")

Jmax_all_final <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/jmax_all_df.csv")

vcmax25_all_final <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/vcmax25_all_df.csv")

Jmax25_all_final <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/jmax25_all_df.csv")

all_forcing <- dplyr::bind_rows(Asat_all_final,Amax_all_final,vcmax_all_final,Jmax_all_final,vcmax25_all_final,Jmax25_all_final)

all_forcing$year_start <- all_forcing$Start_Year
all_forcing$year_end <- all_forcing$Sampling_Year
all_forcing$year_start[all_forcing$year_end>2016] <- 2007
all_forcing$year_end[all_forcing$year_end>2016] <- 2016
summary(all_forcing$year_end-all_forcing$year_start)

all_forcing_sitemean <- aggregate(all_forcing,by=list(all_forcing$lon,all_forcing$lat,all_forcing$z,all_forcing$year_start,all_forcing$year_end), FUN=mean, na.rm=TRUE)
all_forcing_sitemean <- all_forcing_sitemean[,c("lon","lat","z","year_start","year_end")]

all_forcing_sitemean$sitename <- paste("a",1:nrow(all_forcing_sitemean),sep="")
all_forcing_sitemean$elv <- all_forcing_sitemean$z
siteinfo_final <- all_forcing_sitemean[,c("sitename","lon","lat","elv","year_start","year_end")]

#this was universally saved - no change anymore!
#csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/forcing_vcmax.csv")
#write.csv(siteinfo_final, csvfile, row.names = TRUE)


devtools::load_all("~/yunkepeng/gcme/pmodel/ingestr/")
devtools::load_all("~/yunkepeng/gcme/pmodel/rsofun/")



for (i in 1:(nrow(siteinfo_final))){
  df_watch <- ingest(
    siteinfo  = siteinfo_final[i,],
    source    = "watch_wfdei",
    getvars   = c("temp", "prec", "ppfd", "vpd", "patm"), 
    dir       = "/Volumes/My Passport/data/watch_wfdei/",
    settings  = list(correct_bias = "worldclim", dir_bias = "/Volumes/My Passport/data/worldclim/"))
  
  df_cru <- ingest(
    siteinfo  = siteinfo_final[i,],
    source    = "cru",
    getvars   = "ccov",
    dir       = "/Volumes/My Passport/data/cru/ts_4.01/")
  
  df_co2 <- ingest(
    siteinfo  = siteinfo_final[i,],
    source  = "co2_mlo",
    verbose = FALSE)
  
  df_co2_final <- as.data.frame(df_co2$data)
  
  df_co2_final2 <- df_co2_final[!(format(df_co2_final$date,"%m") == "02" & format(df_co2_final$date, "%d") == "29"), , drop = FALSE] # make columns consistent
  
  co2 <- df_co2_final2$co2
  
  ddf_meteo <- as_tibble(cbind(as.data.frame(df_watch$data),as.data.frame(df_cru$data)[,c("ccov_int","ccov")],co2))
  
  ddf_meteo$fapar <- 1
  ddf_meteo$tmax <- ddf_meteo$temp
  ddf_meteo$tmin <- ddf_meteo$temp
  
  ddf_meteo$sitename <- siteinfo_final[i,c("sitename")]
  ddf_meteo$lon <- siteinfo_final[i,c("lon")]
  ddf_meteo$lat <-  siteinfo_final[i,c("lat")]
  ddf_meteo$elv <-  siteinfo_final[i,c("elv")]
  ddf_meteo$year_start <-  siteinfo_final[i,c("year_start")]
  ddf_meteo$year_end <-  siteinfo_final[i,c("year_end")]
  
  csvfile <- paste("~/data/gcme/vcmax_forcing/",siteinfo_final$sitename[i],".csv",sep = "")
  write.csv(ddf_meteo, csvfile, row.names = TRUE)
  print(i)    
}
