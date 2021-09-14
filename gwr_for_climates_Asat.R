#Asat#####
leaf_traits <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/Asat_all_df.csv")
dim(leaf_traits)
gwr_sites <- aggregate(leaf_traits,by=list(leaf_traits$lon,leaf_traits$lat,leaf_traits$z,leaf_traits$Start_Year,leaf_traits$Sampling_Year), FUN=mean, na.rm=TRUE)

gwr_sites <- gwr_sites[,c("lon","lat","z","Start_Year","Sampling_Year")]
names(gwr_sites) <- c("lon","lat","z","start_yr","end_yr")

gwr_sites$z[gwr_sites$z<0] <- 0
#convert after 2016 to 1997-2016
summary(gwr_sites)
gwr_sites$year_start <- gwr_sites$start_yr
gwr_sites$year_end <- gwr_sites$end_yr
gwr_sites$year_start[gwr_sites$end_yr>2016] <- 2007
gwr_sites$year_end[gwr_sites$end_yr>2016] <- 2016

summary(gwr_sites)
dim(gwr_sites)

#input tmx, tmn, vap, pre and radi
load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/climates_30yrs_monthly/WFDEI_CRU_1980_2016.Rdata")
#input qair
load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/climates_30yrs_monthly/WFDEI_qair.Rdata")


total_month <- (2016-1980+1) * 12
library(rbeni)
elev <- as.data.frame(nc_to_df(read_nc_onefile("~/data/watch_wfdei/WFDEI-elevation.nc"), varnam = "elevation"))

monthly_tmn$lon <- elev$lon
monthly_tmn$lat <- elev$lat
monthly_tmn$z <- elev$myvar
monthly_tmn <- monthly_tmn[,c(3: (total_month+5))] #lon, lat, z in last 3
monthly_tmn <- subset(monthly_tmn,z>=0)

monthly_tmx$lon <-  elev$lon
monthly_tmx$lat <- elev$lat
monthly_tmx$z <- elev$myvar
monthly_tmx <- monthly_tmx[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_tmx <- subset(monthly_tmx,z>=0)

monthly_vap$lon <-  elev$lon
monthly_vap$lat <- elev$lat
monthly_vap$z <- elev$myvar
monthly_vap <- monthly_vap[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_vap <- subset(monthly_vap,z>=0)

monthly_pre$lon <-  elev$lon
monthly_pre$lat <-  elev$lat
monthly_pre$z <- elev$myvar
monthly_pre <- monthly_pre[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_pre <- subset(monthly_pre,z>=0)

monthly_radi$lon <- elev$lon
monthly_radi$lat <- elev$lat
monthly_radi$z <- elev$myvar
monthly_radi <- monthly_radi[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_radi <- subset(monthly_radi,z>=0)

monthly_qair$lon <- elev$lon
monthly_qair$lat <- elev$lat
monthly_qair$z <- elev$myvar
monthly_qair <- monthly_qair[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_qair <- subset(monthly_qair,z>=0)

d <- 1.5 # range of surrounded grids (in degree) when applying each gwr

#start gwr 
library(spgwr)
library(raster)

#create a function here
gwr_methods <- function(sites,monthly_grid){
  for (i in c(1:nrow(sites))){
    print(i)
    start_yr <- sites$year_start[i]
    end_yr <- sites$year_end[i]
    lon_site <- sites$lon[i]
    lat_site <- sites$lat[i]
    z_site <- sites$z[i]
    
    climate_distance <- subset(monthly_grid,lon>(lon_site-d)&lon<(lon_site+d)&lat>(lat_site-d)&lat<(lat_site+d))
    climate_distance_yr <- climate_distance[,((start_yr-1979)*12-11):((end_yr-1979)*12)]
    number_of_month <- ncol(climate_distance_yr)
    
    climate_distance_yr_xyz <- cbind(climate_distance[,c("lon","lat","z")],climate_distance_yr)
    # re-create column name for multiple month, which is important when applying gwr within a for loop
    names(climate_distance_yr_xyz) <- c("lon","lat","z",paste("a", 1:number_of_month, sep=""))
    
    #filter the condition when without available climates data after subset (filter within a distance of 1.5)
    if (is.na(mean(climate_distance[,1],na.rm=TRUE))==TRUE){
      sites[i,8:(7+number_of_month)] <- NA
    } else {
      sites[i,8:(7+number_of_month)] <- NA
      coordinates(climate_distance_yr_xyz) <- c("lon","lat")
      gridded(climate_distance_yr_xyz) <- TRUE
      
      input_sites <- sites[i,c("lon","lat","z")]
      coordinates(input_sites) <- c("lon","lat")
      for (a in c(1:number_of_month)){
        #a +5 means, add outputted monthly tmn (within measurement year) directly after lon+lat+z+Begin_year+End_year.
        sites[i,a+7] <- (gwr( get(paste("a",a, sep="")) ~ z, climate_distance_yr_xyz, bandwidth = 1.06, fit.points = input_sites,predictions=TRUE))$SDF$pred
      }
    }
  }
  output <- sites
  return(output)
}


#alpha_output <- gwr_methods(gwr_sites,empty_alpha)
#csvfile <- paste("/Users/yunpeng/data/NPP_final/forest_climates_gwr/forest_alpha.csv")
#write_csv(alpha_output, path = csvfile)
dim(gwr_sites)
tmx_output <- gwr_methods(gwr_sites,monthly_tmx)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/Asat_tmx.csv")
write_csv(tmx_output, path = csvfile)

tmn_output <- gwr_methods(gwr_sites,monthly_tmn)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/Asat_tmn.csv")
write_csv(tmn_output, path = csvfile)

vap_output <- gwr_methods(gwr_sites,monthly_vap)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/Asat_vap.csv")
write_csv(vap_output, path = csvfile)

pre_output <- gwr_methods(gwr_sites,monthly_pre)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/Asat_pre.csv")
write_csv(pre_output, path = csvfile)

radi_output <- gwr_methods(gwr_sites,monthly_radi)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/Asat_radi.csv")
write_csv(radi_output, path = csvfile)

qair_output <- gwr_methods(gwr_sites,monthly_qair)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/Asat_qair.csv")
write_csv(qair_output, path = csvfile)

#Amax####
rm(list=ls())
leaf_traits <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/Amax_all_df.csv")
dim(leaf_traits)
gwr_sites <- aggregate(leaf_traits,by=list(leaf_traits$lon,leaf_traits$lat,leaf_traits$z,leaf_traits$Start_Year,leaf_traits$Sampling_Year), FUN=mean, na.rm=TRUE)
dim(gwr_sites)
gwr_sites <- gwr_sites[,c("lon","lat","z","Start_Year","Sampling_Year")]
names(gwr_sites) <- c("lon","lat","z","start_yr","end_yr")

gwr_sites$z[gwr_sites$z<0] <- 0
#convert after 2016 to 1997-2016
summary(gwr_sites)
gwr_sites$year_start <- gwr_sites$start_yr
gwr_sites$year_end <- gwr_sites$end_yr
gwr_sites$year_start[gwr_sites$end_yr>2016] <- 2007
gwr_sites$year_end[gwr_sites$end_yr>2016] <- 2016

summary(gwr_sites)
dim(gwr_sites)

#input tmx, tmn, vap, pre and radi
load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/climates_30yrs_monthly/WFDEI_CRU_1980_2016.Rdata")
#input qair
load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/climates_30yrs_monthly/WFDEI_qair.Rdata")


total_month <- (2016-1980+1) * 12
library(rbeni)
elev <- as.data.frame(nc_to_df(read_nc_onefile("~/data/watch_wfdei/WFDEI-elevation.nc"), varnam = "elevation"))

monthly_tmn$lon <- elev$lon
monthly_tmn$lat <- elev$lat
monthly_tmn$z <- elev$myvar
monthly_tmn <- monthly_tmn[,c(3: (total_month+5))] #lon, lat, z in last 3
monthly_tmn <- subset(monthly_tmn,z>=0)

monthly_tmx$lon <-  elev$lon
monthly_tmx$lat <- elev$lat
monthly_tmx$z <- elev$myvar
monthly_tmx <- monthly_tmx[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_tmx <- subset(monthly_tmx,z>=0)

monthly_vap$lon <-  elev$lon
monthly_vap$lat <- elev$lat
monthly_vap$z <- elev$myvar
monthly_vap <- monthly_vap[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_vap <- subset(monthly_vap,z>=0)

monthly_pre$lon <-  elev$lon
monthly_pre$lat <-  elev$lat
monthly_pre$z <- elev$myvar
monthly_pre <- monthly_pre[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_pre <- subset(monthly_pre,z>=0)

monthly_radi$lon <- elev$lon
monthly_radi$lat <- elev$lat
monthly_radi$z <- elev$myvar
monthly_radi <- monthly_radi[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_radi <- subset(monthly_radi,z>=0)

monthly_qair$lon <- elev$lon
monthly_qair$lat <- elev$lat
monthly_qair$z <- elev$myvar
monthly_qair <- monthly_qair[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_qair <- subset(monthly_qair,z>=0)

d <- 1.5 # range of surrounded grids (in degree) when applying each gwr

#start gwr 
library(spgwr)
library(raster)

#create a function here
gwr_methods <- function(sites,monthly_grid){
  for (i in c(1:nrow(sites))){
    print(i)
    start_yr <- sites$year_start[i]
    end_yr <- sites$year_end[i]
    lon_site <- sites$lon[i]
    lat_site <- sites$lat[i]
    z_site <- sites$z[i]
    
    climate_distance <- subset(monthly_grid,lon>(lon_site-d)&lon<(lon_site+d)&lat>(lat_site-d)&lat<(lat_site+d))
    climate_distance_yr <- climate_distance[,((start_yr-1979)*12-11):((end_yr-1979)*12)]
    number_of_month <- ncol(climate_distance_yr)
    
    climate_distance_yr_xyz <- cbind(climate_distance[,c("lon","lat","z")],climate_distance_yr)
    # re-create column name for multiple month, which is important when applying gwr within a for loop
    names(climate_distance_yr_xyz) <- c("lon","lat","z",paste("a", 1:number_of_month, sep=""))
    
    #filter the condition when without available climates data after subset (filter within a distance of 1.5)
    if (is.na(mean(climate_distance[,1],na.rm=TRUE))==TRUE){
      sites[i,8:(7+number_of_month)] <- NA
    } else {
      sites[i,8:(7+number_of_month)] <- NA
      coordinates(climate_distance_yr_xyz) <- c("lon","lat")
      gridded(climate_distance_yr_xyz) <- TRUE
      
      input_sites <- sites[i,c("lon","lat","z")]
      coordinates(input_sites) <- c("lon","lat")
      for (a in c(1:number_of_month)){
        #a +5 means, add outputted monthly tmn (within measurement year) directly after lon+lat+z+Begin_year+End_year.
        sites[i,a+7] <- (gwr( get(paste("a",a, sep="")) ~ z, climate_distance_yr_xyz, bandwidth = 1.06, fit.points = input_sites,predictions=TRUE))$SDF$pred
      }
    }
  }
  output <- sites
  return(output)
}


#alpha_output <- gwr_methods(gwr_sites,empty_alpha)
#csvfile <- paste("/Users/yunpeng/data/NPP_final/forest_climates_gwr/forest_alpha.csv")
#write_csv(alpha_output, path = csvfile)
dim(gwr_sites)
tmx_output <- gwr_methods(gwr_sites,monthly_tmx)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/Amax_tmx.csv")
write_csv(tmx_output, path = csvfile)

tmn_output <- gwr_methods(gwr_sites,monthly_tmn)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/Amax_tmn.csv")
write_csv(tmn_output, path = csvfile)

vap_output <- gwr_methods(gwr_sites,monthly_vap)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/Amax_vap.csv")
write_csv(vap_output, path = csvfile)

pre_output <- gwr_methods(gwr_sites,monthly_pre)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/Amax_pre.csv")
write_csv(pre_output, path = csvfile)

radi_output <- gwr_methods(gwr_sites,monthly_radi)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/Amax_radi.csv")
write_csv(radi_output, path = csvfile)

qair_output <- gwr_methods(gwr_sites,monthly_qair)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/Amax_qair.csv")
write_csv(qair_output, path = csvfile)


#vcmax####
rm(list=ls())
leaf_traits <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/vcmax_all_df.csv")
dim(leaf_traits)
gwr_sites <- aggregate(leaf_traits,by=list(leaf_traits$lon,leaf_traits$lat,leaf_traits$z,leaf_traits$Start_Year,leaf_traits$Sampling_Year), FUN=mean, na.rm=TRUE)
dim(gwr_sites)
gwr_sites <- gwr_sites[,c("lon","lat","z","Start_Year","Sampling_Year")]
names(gwr_sites) <- c("lon","lat","z","start_yr","end_yr")
summary(gwr_sites)
gwr_sites$z[gwr_sites$z<0] <- 0
#convert after 2016 to 1997-2016
summary(gwr_sites)
gwr_sites$year_start <- gwr_sites$start_yr
gwr_sites$year_end <- gwr_sites$end_yr
gwr_sites$year_start[gwr_sites$end_yr>2016] <- 2007
gwr_sites$year_end[gwr_sites$end_yr>2016] <- 2016

summary(gwr_sites)
dim(gwr_sites)

#input tmx, tmn, vap, pre and radi
load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/climates_30yrs_monthly/WFDEI_CRU_1980_2016.Rdata")
#input qair
load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/climates_30yrs_monthly/WFDEI_qair.Rdata")


total_month <- (2016-1980+1) * 12
library(rbeni)
elev <- as.data.frame(nc_to_df(read_nc_onefile("~/data/watch_wfdei/WFDEI-elevation.nc"), varnam = "elevation"))

monthly_tmn$lon <- elev$lon
monthly_tmn$lat <- elev$lat
monthly_tmn$z <- elev$myvar
monthly_tmn <- monthly_tmn[,c(3: (total_month+5))] #lon, lat, z in last 3
monthly_tmn <- subset(monthly_tmn,z>=0)

monthly_tmx$lon <-  elev$lon
monthly_tmx$lat <- elev$lat
monthly_tmx$z <- elev$myvar
monthly_tmx <- monthly_tmx[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_tmx <- subset(monthly_tmx,z>=0)

monthly_vap$lon <-  elev$lon
monthly_vap$lat <- elev$lat
monthly_vap$z <- elev$myvar
monthly_vap <- monthly_vap[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_vap <- subset(monthly_vap,z>=0)

monthly_pre$lon <-  elev$lon
monthly_pre$lat <-  elev$lat
monthly_pre$z <- elev$myvar
monthly_pre <- monthly_pre[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_pre <- subset(monthly_pre,z>=0)

monthly_radi$lon <- elev$lon
monthly_radi$lat <- elev$lat
monthly_radi$z <- elev$myvar
monthly_radi <- monthly_radi[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_radi <- subset(monthly_radi,z>=0)

monthly_qair$lon <- elev$lon
monthly_qair$lat <- elev$lat
monthly_qair$z <- elev$myvar
monthly_qair <- monthly_qair[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_qair <- subset(monthly_qair,z>=0)

d <- 1.5 # range of surrounded grids (in degree) when applying each gwr

#start gwr 
library(spgwr)
library(raster)

#create a function here
gwr_methods <- function(sites,monthly_grid){
  for (i in c(1:nrow(sites))){
    print(i)
    start_yr <- sites$year_start[i]
    end_yr <- sites$year_end[i]
    lon_site <- sites$lon[i]
    lat_site <- sites$lat[i]
    z_site <- sites$z[i]
    
    climate_distance <- subset(monthly_grid,lon>(lon_site-d)&lon<(lon_site+d)&lat>(lat_site-d)&lat<(lat_site+d))
    climate_distance_yr <- climate_distance[,((start_yr-1979)*12-11):((end_yr-1979)*12)]
    number_of_month <- ncol(climate_distance_yr)
    
    climate_distance_yr_xyz <- cbind(climate_distance[,c("lon","lat","z")],climate_distance_yr)
    # re-create column name for multiple month, which is important when applying gwr within a for loop
    names(climate_distance_yr_xyz) <- c("lon","lat","z",paste("a", 1:number_of_month, sep=""))
    
    #filter the condition when without available climates data after subset (filter within a distance of 1.5)
    if (is.na(mean(climate_distance[,1],na.rm=TRUE))==TRUE){
      sites[i,8:(7+number_of_month)] <- NA
    } else {
      sites[i,8:(7+number_of_month)] <- NA
      coordinates(climate_distance_yr_xyz) <- c("lon","lat")
      gridded(climate_distance_yr_xyz) <- TRUE
      
      input_sites <- sites[i,c("lon","lat","z")]
      coordinates(input_sites) <- c("lon","lat")
      for (a in c(1:number_of_month)){
        #a +5 means, add outputted monthly tmn (within measurement year) directly after lon+lat+z+Begin_year+End_year.
        sites[i,a+7] <- (gwr( get(paste("a",a, sep="")) ~ z, climate_distance_yr_xyz, bandwidth = 1.06, fit.points = input_sites,predictions=TRUE))$SDF$pred
      }
    }
  }
  output <- sites
  return(output)
}


#alpha_output <- gwr_methods(gwr_sites,empty_alpha)
#csvfile <- paste("/Users/yunpeng/data/NPP_final/forest_climates_gwr/forest_alpha.csv")
#write_csv(alpha_output, path = csvfile)
dim(gwr_sites)
tmx_output <- gwr_methods(gwr_sites,monthly_tmx)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/vcmax_tmx.csv")
write_csv(tmx_output, path = csvfile)

tmn_output <- gwr_methods(gwr_sites,monthly_tmn)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/vcmax_tmn.csv")
write_csv(tmn_output, path = csvfile)

vap_output <- gwr_methods(gwr_sites,monthly_vap)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/vcmax_vap.csv")
write_csv(vap_output, path = csvfile)

pre_output <- gwr_methods(gwr_sites,monthly_pre)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/vcmax_pre.csv")
write_csv(pre_output, path = csvfile)

radi_output <- gwr_methods(gwr_sites,monthly_radi)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/vcmax_radi.csv")
write_csv(radi_output, path = csvfile)

qair_output <- gwr_methods(gwr_sites,monthly_qair)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/vcmax_qair.csv")
write_csv(qair_output, path = csvfile)

#jmax####
rm(list=ls())
leaf_traits <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/jmax_all_df.csv")
dim(leaf_traits)
gwr_sites <- aggregate(leaf_traits,by=list(leaf_traits$lon,leaf_traits$lat,leaf_traits$z,leaf_traits$Start_Year,leaf_traits$Sampling_Year), FUN=mean, na.rm=TRUE)
dim(gwr_sites)
gwr_sites <- gwr_sites[,c("lon","lat","z","Start_Year","Sampling_Year")]
names(gwr_sites) <- c("lon","lat","z","start_yr","end_yr")

gwr_sites$z[gwr_sites$z<0] <- 0
#convert after 2016 to 1997-2016
summary(gwr_sites)
gwr_sites$year_start <- gwr_sites$start_yr
gwr_sites$year_end <- gwr_sites$end_yr
gwr_sites$year_start[gwr_sites$end_yr>2016] <- 2007
gwr_sites$year_end[gwr_sites$end_yr>2016] <- 2016

summary(gwr_sites)
dim(gwr_sites)

#input tmx, tmn, vap, pre and radi
load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/climates_30yrs_monthly/WFDEI_CRU_1980_2016.Rdata")
#input qair
load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/climates_30yrs_monthly/WFDEI_qair.Rdata")


total_month <- (2016-1980+1) * 12
library(rbeni)
elev <- as.data.frame(nc_to_df(read_nc_onefile("~/data/watch_wfdei/WFDEI-elevation.nc"), varnam = "elevation"))

monthly_tmn$lon <- elev$lon
monthly_tmn$lat <- elev$lat
monthly_tmn$z <- elev$myvar
monthly_tmn <- monthly_tmn[,c(3: (total_month+5))] #lon, lat, z in last 3
monthly_tmn <- subset(monthly_tmn,z>=0)

monthly_tmx$lon <-  elev$lon
monthly_tmx$lat <- elev$lat
monthly_tmx$z <- elev$myvar
monthly_tmx <- monthly_tmx[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_tmx <- subset(monthly_tmx,z>=0)

monthly_vap$lon <-  elev$lon
monthly_vap$lat <- elev$lat
monthly_vap$z <- elev$myvar
monthly_vap <- monthly_vap[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_vap <- subset(monthly_vap,z>=0)

monthly_pre$lon <-  elev$lon
monthly_pre$lat <-  elev$lat
monthly_pre$z <- elev$myvar
monthly_pre <- monthly_pre[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_pre <- subset(monthly_pre,z>=0)

monthly_radi$lon <- elev$lon
monthly_radi$lat <- elev$lat
monthly_radi$z <- elev$myvar
monthly_radi <- monthly_radi[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_radi <- subset(monthly_radi,z>=0)

monthly_qair$lon <- elev$lon
monthly_qair$lat <- elev$lat
monthly_qair$z <- elev$myvar
monthly_qair <- monthly_qair[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_qair <- subset(monthly_qair,z>=0)

d <- 1.5 # range of surrounded grids (in degree) when applying each gwr

#start gwr 
library(spgwr)
library(raster)

#create a function here
gwr_methods <- function(sites,monthly_grid){
  for (i in c(1:nrow(sites))){
    print(i)
    start_yr <- sites$year_start[i]
    end_yr <- sites$year_end[i]
    lon_site <- sites$lon[i]
    lat_site <- sites$lat[i]
    z_site <- sites$z[i]
    
    climate_distance <- subset(monthly_grid,lon>(lon_site-d)&lon<(lon_site+d)&lat>(lat_site-d)&lat<(lat_site+d))
    climate_distance_yr <- climate_distance[,((start_yr-1979)*12-11):((end_yr-1979)*12)]
    number_of_month <- ncol(climate_distance_yr)
    
    climate_distance_yr_xyz <- cbind(climate_distance[,c("lon","lat","z")],climate_distance_yr)
    # re-create column name for multiple month, which is important when applying gwr within a for loop
    names(climate_distance_yr_xyz) <- c("lon","lat","z",paste("a", 1:number_of_month, sep=""))
    
    #filter the condition when without available climates data after subset (filter within a distance of 1.5)
    if (is.na(mean(climate_distance[,1],na.rm=TRUE))==TRUE){
      sites[i,8:(7+number_of_month)] <- NA
    } else {
      sites[i,8:(7+number_of_month)] <- NA
      coordinates(climate_distance_yr_xyz) <- c("lon","lat")
      gridded(climate_distance_yr_xyz) <- TRUE
      
      input_sites <- sites[i,c("lon","lat","z")]
      coordinates(input_sites) <- c("lon","lat")
      for (a in c(1:number_of_month)){
        #a +5 means, add outputted monthly tmn (within measurement year) directly after lon+lat+z+Begin_year+End_year.
        sites[i,a+7] <- (gwr( get(paste("a",a, sep="")) ~ z, climate_distance_yr_xyz, bandwidth = 1.06, fit.points = input_sites,predictions=TRUE))$SDF$pred
      }
    }
  }
  output <- sites
  return(output)
}


#alpha_output <- gwr_methods(gwr_sites,empty_alpha)
#csvfile <- paste("/Users/yunpeng/data/NPP_final/forest_climates_gwr/forest_alpha.csv")
#write_csv(alpha_output, path = csvfile)
dim(gwr_sites)
tmx_output <- gwr_methods(gwr_sites,monthly_tmx)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/jmax_tmx.csv")
write_csv(tmx_output, path = csvfile)

tmn_output <- gwr_methods(gwr_sites,monthly_tmn)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/jmax_tmn.csv")
write_csv(tmn_output, path = csvfile)

vap_output <- gwr_methods(gwr_sites,monthly_vap)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/jmax_vap.csv")
write_csv(vap_output, path = csvfile)

pre_output <- gwr_methods(gwr_sites,monthly_pre)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/jmax_pre.csv")
write_csv(pre_output, path = csvfile)

radi_output <- gwr_methods(gwr_sites,monthly_radi)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/jmax_radi.csv")
write_csv(radi_output, path = csvfile)

qair_output <- gwr_methods(gwr_sites,monthly_qair)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/jmax_qair.csv")
write_csv(qair_output, path = csvfile)

#vcmax25####
rm(list=ls())
leaf_traits <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/vcmax25_all_df.csv")
dim(leaf_traits)
gwr_sites <- aggregate(leaf_traits,by=list(leaf_traits$lon,leaf_traits$lat,leaf_traits$z,leaf_traits$Start_Year,leaf_traits$Sampling_Year), FUN=mean, na.rm=TRUE)
dim(gwr_sites)
gwr_sites <- gwr_sites[,c("lon","lat","z","Start_Year","Sampling_Year")]
names(gwr_sites) <- c("lon","lat","z","start_yr","end_yr")

gwr_sites$z[gwr_sites$z<0] <- 0
#convert after 2016 to 1997-2016
summary(gwr_sites)
gwr_sites$year_start <- gwr_sites$start_yr
gwr_sites$year_end <- gwr_sites$end_yr
gwr_sites$year_start[gwr_sites$end_yr>2016] <- 2007
gwr_sites$year_end[gwr_sites$end_yr>2016] <- 2016

summary(gwr_sites)
dim(gwr_sites)

#input tmx, tmn, vap, pre and radi
load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/climates_30yrs_monthly/WFDEI_CRU_1980_2016.Rdata")
#input qair
load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/climates_30yrs_monthly/WFDEI_qair.Rdata")


total_month <- (2016-1980+1) * 12
library(rbeni)
elev <- as.data.frame(nc_to_df(read_nc_onefile("~/data/watch_wfdei/WFDEI-elevation.nc"), varnam = "elevation"))

monthly_tmn$lon <- elev$lon
monthly_tmn$lat <- elev$lat
monthly_tmn$z <- elev$myvar
monthly_tmn <- monthly_tmn[,c(3: (total_month+5))] #lon, lat, z in last 3
monthly_tmn <- subset(monthly_tmn,z>=0)

monthly_tmx$lon <-  elev$lon
monthly_tmx$lat <- elev$lat
monthly_tmx$z <- elev$myvar
monthly_tmx <- monthly_tmx[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_tmx <- subset(monthly_tmx,z>=0)

monthly_vap$lon <-  elev$lon
monthly_vap$lat <- elev$lat
monthly_vap$z <- elev$myvar
monthly_vap <- monthly_vap[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_vap <- subset(monthly_vap,z>=0)

monthly_pre$lon <-  elev$lon
monthly_pre$lat <-  elev$lat
monthly_pre$z <- elev$myvar
monthly_pre <- monthly_pre[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_pre <- subset(monthly_pre,z>=0)

monthly_radi$lon <- elev$lon
monthly_radi$lat <- elev$lat
monthly_radi$z <- elev$myvar
monthly_radi <- monthly_radi[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_radi <- subset(monthly_radi,z>=0)

monthly_qair$lon <- elev$lon
monthly_qair$lat <- elev$lat
monthly_qair$z <- elev$myvar
monthly_qair <- monthly_qair[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_qair <- subset(monthly_qair,z>=0)

d <- 1.5 # range of surrounded grids (in degree) when applying each gwr

#start gwr 
library(spgwr)
library(raster)

#create a function here
gwr_methods <- function(sites,monthly_grid){
  for (i in c(1:nrow(sites))){
    print(i)
    start_yr <- sites$year_start[i]
    end_yr <- sites$year_end[i]
    lon_site <- sites$lon[i]
    lat_site <- sites$lat[i]
    z_site <- sites$z[i]
    
    climate_distance <- subset(monthly_grid,lon>(lon_site-d)&lon<(lon_site+d)&lat>(lat_site-d)&lat<(lat_site+d))
    climate_distance_yr <- climate_distance[,((start_yr-1979)*12-11):((end_yr-1979)*12)]
    number_of_month <- ncol(climate_distance_yr)
    
    climate_distance_yr_xyz <- cbind(climate_distance[,c("lon","lat","z")],climate_distance_yr)
    # re-create column name for multiple month, which is important when applying gwr within a for loop
    names(climate_distance_yr_xyz) <- c("lon","lat","z",paste("a", 1:number_of_month, sep=""))
    
    #filter the condition when without available climates data after subset (filter within a distance of 1.5)
    if (is.na(mean(climate_distance[,1],na.rm=TRUE))==TRUE){
      sites[i,8:(7+number_of_month)] <- NA
    } else {
      sites[i,8:(7+number_of_month)] <- NA
      coordinates(climate_distance_yr_xyz) <- c("lon","lat")
      gridded(climate_distance_yr_xyz) <- TRUE
      
      input_sites <- sites[i,c("lon","lat","z")]
      coordinates(input_sites) <- c("lon","lat")
      for (a in c(1:number_of_month)){
        #a +5 means, add outputted monthly tmn (within measurement year) directly after lon+lat+z+Begin_year+End_year.
        sites[i,a+7] <- (gwr( get(paste("a",a, sep="")) ~ z, climate_distance_yr_xyz, bandwidth = 1.06, fit.points = input_sites,predictions=TRUE))$SDF$pred
      }
    }
  }
  output <- sites
  return(output)
}


#alpha_output <- gwr_methods(gwr_sites,empty_alpha)
#csvfile <- paste("/Users/yunpeng/data/NPP_final/forest_climates_gwr/forest_alpha.csv")
#write_csv(alpha_output, path = csvfile)
dim(gwr_sites)
tmx_output <- gwr_methods(gwr_sites,monthly_tmx)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/vcmax25_tmx.csv")
write_csv(tmx_output, path = csvfile)

tmn_output <- gwr_methods(gwr_sites,monthly_tmn)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/vcmax25_tmn.csv")
write_csv(tmn_output, path = csvfile)

vap_output <- gwr_methods(gwr_sites,monthly_vap)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/vcmax25_vap.csv")
write_csv(vap_output, path = csvfile)

pre_output <- gwr_methods(gwr_sites,monthly_pre)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/vcmax25_pre.csv")
write_csv(pre_output, path = csvfile)

radi_output <- gwr_methods(gwr_sites,monthly_radi)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/vcmax25_radi.csv")
write_csv(radi_output, path = csvfile)

qair_output <- gwr_methods(gwr_sites,monthly_qair)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/vcmax25_qair.csv")
write_csv(qair_output, path = csvfile)

#jmax25####
rm(list=ls())
leaf_traits <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/jmax25_all_df.csv")
dim(leaf_traits)
gwr_sites <- aggregate(leaf_traits,by=list(leaf_traits$lon,leaf_traits$lat,leaf_traits$z,leaf_traits$Start_Year,leaf_traits$Sampling_Year), FUN=mean, na.rm=TRUE)
dim(gwr_sites)
gwr_sites <- gwr_sites[,c("lon","lat","z","Start_Year","Sampling_Year")]
names(gwr_sites) <- c("lon","lat","z","start_yr","end_yr")
summary(gwr_sites)

gwr_sites$z[gwr_sites$z<0] <- 0
#convert after 2016 to 1997-2016
summary(gwr_sites)
gwr_sites$year_start <- gwr_sites$start_yr
gwr_sites$year_end <- gwr_sites$end_yr
gwr_sites$year_start[gwr_sites$end_yr>2016] <- 2007
gwr_sites$year_end[gwr_sites$end_yr>2016] <- 2016

summary(gwr_sites)
dim(gwr_sites)

#input tmx, tmn, vap, pre and radi
load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/climates_30yrs_monthly/WFDEI_CRU_1980_2016.Rdata")
#input qair
load(file = "/Users/yunpeng/yunkepeng/nimpl_sofun_inputs_final/climates_30yrs_monthly/WFDEI_qair.Rdata")


total_month <- (2016-1980+1) * 12
library(rbeni)
elev <- as.data.frame(nc_to_df(read_nc_onefile("~/data/watch_wfdei/WFDEI-elevation.nc"), varnam = "elevation"))

monthly_tmn$lon <- elev$lon
monthly_tmn$lat <- elev$lat
monthly_tmn$z <- elev$myvar
monthly_tmn <- monthly_tmn[,c(3: (total_month+5))] #lon, lat, z in last 3
monthly_tmn <- subset(monthly_tmn,z>=0)

monthly_tmx$lon <-  elev$lon
monthly_tmx$lat <- elev$lat
monthly_tmx$z <- elev$myvar
monthly_tmx <- monthly_tmx[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_tmx <- subset(monthly_tmx,z>=0)

monthly_vap$lon <-  elev$lon
monthly_vap$lat <- elev$lat
monthly_vap$z <- elev$myvar
monthly_vap <- monthly_vap[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_vap <- subset(monthly_vap,z>=0)

monthly_pre$lon <-  elev$lon
monthly_pre$lat <-  elev$lat
monthly_pre$z <- elev$myvar
monthly_pre <- monthly_pre[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_pre <- subset(monthly_pre,z>=0)

monthly_radi$lon <- elev$lon
monthly_radi$lat <- elev$lat
monthly_radi$z <- elev$myvar
monthly_radi <- monthly_radi[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_radi <- subset(monthly_radi,z>=0)

monthly_qair$lon <- elev$lon
monthly_qair$lat <- elev$lat
monthly_qair$z <- elev$myvar
monthly_qair <- monthly_qair[,c(3:(total_month+5))] #lon, lat, z in last 3
monthly_qair <- subset(monthly_qair,z>=0)

d <- 1.5 # range of surrounded grids (in degree) when applying each gwr

#start gwr 
library(spgwr)
library(raster)

#create a function here
gwr_methods <- function(sites,monthly_grid){
  for (i in c(1:nrow(sites))){
    print(i)
    start_yr <- sites$year_start[i]
    end_yr <- sites$year_end[i]
    lon_site <- sites$lon[i]
    lat_site <- sites$lat[i]
    z_site <- sites$z[i]
    
    climate_distance <- subset(monthly_grid,lon>(lon_site-d)&lon<(lon_site+d)&lat>(lat_site-d)&lat<(lat_site+d))
    climate_distance_yr <- climate_distance[,((start_yr-1979)*12-11):((end_yr-1979)*12)]
    number_of_month <- ncol(climate_distance_yr)
    
    climate_distance_yr_xyz <- cbind(climate_distance[,c("lon","lat","z")],climate_distance_yr)
    # re-create column name for multiple month, which is important when applying gwr within a for loop
    names(climate_distance_yr_xyz) <- c("lon","lat","z",paste("a", 1:number_of_month, sep=""))
    
    #filter the condition when without available climates data after subset (filter within a distance of 1.5)
    if (is.na(mean(climate_distance[,1],na.rm=TRUE))==TRUE){
      sites[i,8:(7+number_of_month)] <- NA
    } else {
      sites[i,8:(7+number_of_month)] <- NA
      coordinates(climate_distance_yr_xyz) <- c("lon","lat")
      gridded(climate_distance_yr_xyz) <- TRUE
      
      input_sites <- sites[i,c("lon","lat","z")]
      coordinates(input_sites) <- c("lon","lat")
      for (a in c(1:number_of_month)){
        #a +5 means, add outputted monthly tmn (within measurement year) directly after lon+lat+z+Begin_year+End_year.
        sites[i,a+7] <- (gwr( get(paste("a",a, sep="")) ~ z, climate_distance_yr_xyz, bandwidth = 1.06, fit.points = input_sites,predictions=TRUE))$SDF$pred
      }
    }
  }
  output <- sites
  return(output)
}


#alpha_output <- gwr_methods(gwr_sites,empty_alpha)
#csvfile <- paste("/Users/yunpeng/data/NPP_final/forest_climates_gwr/forest_alpha.csv")
#write_csv(alpha_output, path = csvfile)
dim(gwr_sites)
tmx_output <- gwr_methods(gwr_sites,monthly_tmx)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/jmax25_tmx.csv")
write_csv(tmx_output, path = csvfile)

tmn_output <- gwr_methods(gwr_sites,monthly_tmn)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/jmax25_tmn.csv")
write_csv(tmn_output, path = csvfile)

vap_output <- gwr_methods(gwr_sites,monthly_vap)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/jmax25_vap.csv")
write_csv(vap_output, path = csvfile)

pre_output <- gwr_methods(gwr_sites,monthly_pre)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/jmax25_pre.csv")
write_csv(pre_output, path = csvfile)

radi_output <- gwr_methods(gwr_sites,monthly_radi)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/jmax25_radi.csv")
write_csv(radi_output, path = csvfile)

qair_output <- gwr_methods(gwr_sites,monthly_qair)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/climates/jmax25_qair.csv")
write_csv(qair_output, path = csvfile)