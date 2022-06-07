#co2 effect
#~/yunkepeng/gcme/forcing_kevin_nick.R
#kevin
kevin <- read.csv("/Users/yunpeng/data/gcme/kevin/forcing/forcing_info.csv")[,c("sitename","lon","lat","elv")]

#smith
smith_co2 <- read.csv("/Users/yunpeng/data/smith_keenan_gcb/gcb_co2/co2_data.csv")
smith_co2$exp_nam <- smith_co2$SiteID
smith_co2 <- subset(smith_co2,is.na(lat)==FALSE) # remove empty column

#create year from 1992 to 2016 -- same to GCME
smith_co2$year_start <- 1982
smith_co2$year_end <- 2016

nick_forcing_sitemean <- aggregate(smith_co2,by=list(smith_co2$lon,smith_co2$lat,smith_co2$z,smith_co2$year_start,smith_co2$year_end), FUN=mean, na.rm=TRUE)
nick_forcing_sitemean <- nick_forcing_sitemean[,c("lon","lat","z","year_start","year_end")]
dim(nick_forcing_sitemean)

nick_forcing_sitemean$sitename <- paste("nick",1:nrow(nick_forcing_sitemean),sep="")
nick_forcing_sitemean$elv <- nick_forcing_sitemean$z
smith <- nick_forcing_sitemean[,c("sitename","lon","lat","elv")]

#light
#~/yunkepeng/gcme/light_forcing.R
library(dplyr)
walker <- read.csv("/Users/yunpeng/data/leaf_traits/Walker/LEAF_PHOTOSYNTHESIS_TRAITS_1224/data/Leaf_Photosynthesis_Traits.csv")

#vcmax
light_vcmax <- subset(walker,Treatment_light!="not applicable") %>% 
  group_by(lon,lat,z,Year,species,Treatment_light)  %>% summarise(mean = mean(Vcmax))

#prepare forcing
aaa <- as.data.frame(unique(light_vcmax[,c("lon","lat","z")]))

aaa$year_start <- 1991
aaa$year_end <- 2010

aaa$elv <- aaa$z
aaa$sitename <- paste("light",1:nrow(aaa),sep="")

siteinfo_final <- aaa[,c("sitename","lon","lat","elv","year_start","year_end")]
light <- siteinfo_final[,c("sitename","lon","lat","elv")]

#fertilization from Walker
Nfer_vcmax <- subset(walker,Treatment_N!="not applicable" &Treatment_CO2!="not applicable" & Vcmax>0) %>% 
  group_by(lon,lat,z,species,Treatment_N,Treatment_CO2)  %>% summarise(mean = mean(Vcmax))
bbb <- as.data.frame(unique(Nfer_vcmax[,c("lon","lat","z")]))

bbb$elv <- bbb$z
bbb$sitename <- paste("Nfer",1:nrow(bbb),sep="")
Nfer <- bbb[,c("sitename","lon","lat","elv")]

#warming, smith
sitename <- "smith_keenan1"
lon <- -86.92
lat <- 40.82
elv <- 187
smith_keenan1 <- data.frame(sitename,lon,lat,elv)

#warming, smith from BACE
sitename <- "smith_keenan2"
lon <- -71.21
lat <- 42.38
smith_keenan2 <- data.frame(sitename,lon,lat)
devtools::load_all("~/yunkepeng/gcme/pmodel/ingestr/")
df_etopo <- ingest(smith_keenan2,source = "etopo1",dir = "~/data/etopo/" )
smith_keenan2$elv <- as.numeric(as.data.frame(df_etopo$data))


#warming - dushan
library(readr)
library(dplyr)
df2 <- read_csv("/Users/yunpeng/data/Kumarathunge_2020_newphy/kumarathunge_2020_newphy/Data/PPC-TGlob_V1.0.csv")
df2 <- subset(df2,Growth_condition=="Glasshouse"|Growth_condition=="Whole tree chambers")
dff <- na.omit(as.data.frame(unique(df2[,c("seed_source_latitude","seed_source_longitude")])))
names(dff) <- c("lat","lon")

dff$sitename <- paste("warming",1:length(dff$lon),sep="")
devtools::load_all("~/yunkepeng/gcme/pmodel/ingestr/")
df_etopo <- ingest(dff,source = "etopo1",dir = "~/data/etopo/" )
dff$elv <- as.numeric(as.data.frame(df_etopo$data))


siteinfo_final <- dplyr::bind_rows(kevin,smith,light,Nfer,smith_keenan1,smith_keenan2,dff)
siteinfo_final$year_start <- 1991
siteinfo_final$year_end <- 2010

csvfile <- paste("~/data/gcme/kevin/forcing/forcing_info_allnew.csv")
write.csv(siteinfo_final, csvfile, row.names = TRUE)

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
  
  csvfile <- paste("~/data/gcme/kevin/forcing/climate_large/",siteinfo_final$sitename[i],".csv",sep = "")
  write.csv(ddf_meteo, csvfile, row.names = TRUE)
  print(i)    
}

