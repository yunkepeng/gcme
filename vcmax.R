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
df <- read_csv("~/data/gcme/data_received_190325/NewData_wide_CORRECTED2.csv") %>%
  mutate( ambient_Sd  = as.numeric(ambient_Sd),  ambient_Se  = as.numeric(ambient_Se), 
          elevated_Sd = as.numeric(elevated_Sd), elevated_Se = as.numeric(elevated_Se) )

# save experiments names
df_experiments <- df %>% dplyr::select(exp_nam, prev_name) %>% distinct()

### Duplicate rows
#Some rows are simply duplicates. This is dangerous.
#keyvars <- c("exp_nam","prev_name", "factors", "treatment", "Data_type", "Unit", "Sampling_date", "Start_Year", "Year", "n_plots")
#valuevars <-  c("ambient", "ambient_Se", "ambient_Sd", "elevated", "elevated_Se", "elevated_Sd")

#This may be remediated automatically by using only columns that are distinct w.r.t. `valuevars` and `keyvars` (basically all columns in original data, except `ALIAS`).
#print(nrow(df))
#df <- df %>% distinct_at(vars(one_of(keyvars, valuevars))) %>%mutate( id=1:n() )   # create new ID key (before this was 'ALIAS')
#print(nrow(df))

## Determine all experiments that have more than 1 year data
list_exp_gt1yr <- df %>% filter(!is.na(Year)) %>% group_by(exp_nam) %>% summarise(nyears=max(Year)) %>% filter(nyears>1) %>% dplyr::select(exp_nam) %>% unlist() %>% unname()

#include coordinates and experimental details
coord <- read.csv("~/data/gcme/data_received_190325/table_var_exp_names.csv")[,c("Site.Name","System.Type","Latitude","Lat","Longitude","Lon","Elevation","AM","EcM","ErM")]
for (i in 1:nrow(coord)){
  if (coord$Lat[i] == "S"){
    coord$Latitude[i] <- -as.numeric(coord$Latitude[i])
  } else {
    coord$Latitude[i] <- as.numeric(coord$Latitude[i])
  }
  if (coord$Lon[i] == "W"){
    coord$Longitude[i] <- -as.numeric(coord$Longitude[i])
  } else {
    coord$Longitude[i] <- as.numeric(coord$Longitude[i])
  }
}
coord <- coord[,c("Site.Name","System.Type","Latitude","Longitude","Elevation","AM","EcM","ErM")]
names(coord) <- c("exp_nam","pft","lat","lon","z","AM","EcM","ErM")
coord$rep <- duplicated(coord$exp_nam)
coord <- subset(coord,rep=="FALSE")# remove 11 duplicate site

coord2 <- coord[,c("exp_nam","lat","lon","z")]
csvfile <- paste("~/data/gcme/data_received_190325/table_var_exp_names_coordinates.csv")
write.csv(coord2, csvfile, row.names = TRUE)

df_experiments_coord <- merge(df_experiments,coord,by=c("exp_nam"),all.x=TRUE)

experimental <- read.csv("~/data/gcme/data_received_190325/table_var_exp_names_experiments.csv")[,c("Site.Name","prev_name","Start_Year","End_Year","Cquantity_Control","Cquantity_Treatment",
                                                                                                    "Fumigation_type","Nutrients","Type","N.Quantity.control","N.Quantity.treatment",
                                                                                                    "N.Quantity.unit","P.Quantity.control","P.Quantity.treatment","P.Quantity.unit")]

names(experimental) <- c("exp_nam","prev_name","Start_Year","End_Year","Cquantity_Control","Cquantity_Treatment",
                         "Fumigation_type","Nutrients","Type","N.Quantity.control","N.Quantity.treatment",
                         "N.Quantity.unit","P.Quantity.control","P.Quantity.treatment","P.Quantity.unit")

df_experiments_coord_exp <- merge(df_experiments_coord,experimental,by=c("exp_nam","prev_name"),all.x=TRUE)

check_na <- subset(df_experiments_coord_exp,is.na(lon)==TRUE) %>% dplyr::select(exp_nam) %>% distinct()
dim(check_na)
check_yes <- subset(df_experiments_coord_exp,is.na(lon)==FALSE) %>% dplyr::select(exp_nam) %>% distinct()
dim(check_yes)
#output to site coordinates csv
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/coord_df.csv")
write.csv(df_experiments_coord_exp, csvfile, row.names = TRUE)

#for leaf ci/ca, leaf N, vcmax, jmax, c13
Asat_leafN_c13_all <- subset(df,Data_type=="Ci"|Data_type=="Ci/Ca"|
                     Data_type=="leaf_N"|Data_type=="Asat"|Data_type=="Amax"|
                       Data_type=="CO2_assimilation_rate_(Asat)"|Data_type=="CO2_assimilation_rate_(Amax)"|
                     Data_type=="maximal_Rubisco_activity_(Vcmax)"|Data_type=="Vcmax"|
                     Data_type=="Vcmax,_maximum_carboxylation_rate_of_Rubisco"|Data_type=="Vcmax25"|
                     Data_type=="Jmax"|Data_type=="Jmax25"|Data_type=="Jmax_"|
                     Data_type=="Jmax,_maximum_rate_of_electron_transport_for_regeneration_of_RuBP"|
                     Data_type=="leaf_[C]"|Data_type=="leaf_?13C"|Data_type=="leaf_d13C")

Asat_all2 <- df[grep("photo", df$Data_type),] # for photosynthesis capacity (Asat)

Asat_all3 <- df[grep("water_use", df$Data_type),] # for ci/ca

Asat_leafN_final <- rbind(Asat_leafN_c13_all,Asat_all2,Asat_all3)

leaf_all_coord <- merge(Asat_leafN_final,df_experiments_coord_exp,by=c("exp_nam","prev_name"),all.x=TRUE)
dim(Asat_leafN_final)
dim(leaf_all_coord)

coord_Asat <- leaf_all_coord %>% dplyr::select(exp_nam,lat,lon,z) %>% distinct()
dim(coord_Asat)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/coord_Asat_leafN.csv")
write.csv(coord_Asat, csvfile, row.names = TRUE)
#105 sites overall!
#then, after looking at sources from..., this was then filled by 
#/Users/yunpeng/data/gcme/yunke_coord/coord_Asat_leafN_filled.csv
#ingest elevation
leaf_sites <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/coord_Asat_leafN_filled.csv")
leaf_sites <- subset(leaf_sites,is.na(lon)==FALSE)# remove an AUS site without coordinates
dim(leaf_sites)
devtools::load_all("~/yunkepeng/gcme/pmodel/ingestr/")
aaa <- aggregate(leaf_sites,by=list(leaf_sites$lon,leaf_sites$lat), FUN=mean, na.rm=TRUE)
aaa <- aaa[,c("lon","lat")]
aaa$sitename <- paste("c",1:length(aaa$lon),sep="")

df_etopo <- ingest(aaa,source = "etopo1",dir = "~/data/etopo/" )
aaa$elv <- as.numeric(as.data.frame(df_etopo$data))
aaa$elv[aaa$elv<0] <- 0 #two points are negative, change to 0
aaa <- aaa[,c("lon","lat","elv")]
leaf_sites_z <- merge(leaf_sites,aaa,by=c("lon","lat"),all.x=TRUE)
plot(leaf_sites_z$z~leaf_sites_z$elv) # looks ok - now interploate original elevation value with etopo elevation
leaf_sites_z$z[is.na(leaf_sites_z$z)==TRUE] <- leaf_sites_z$elv[is.na(leaf_sites_z$z)==TRUE] 
leaf_sites_z_final <- leaf_sites_z[,c("exp_nam","lon","lat","z")]

df_leafsites <- merge(df,leaf_sites_z_final,by=c("exp_nam"),all.x=TRUE)
dim(df_leafsites)

#some years were missing, interpolating them primarily basing on max measurement year from the same plot and same variable
# secondarily basing on max measurement year from the same plot
NA_years <- as.data.frame(subset(df_leafsites,is.na(Sampling_Year)==TRUE)%>% group_by(exp_nam,Data_type)%>% summarise(number = n()))

for (i in 1:nrow(NA_years)){
  exp_name <- NA_years[i,1]
  variable_name <- NA_years[i,2]
  if (is.infinite(max(df_leafsites$Sampling_Year[is.na(df_leafsites$Sampling_Year)==FALSE & df_leafsites$exp_nam==exp_name & df_leafsites$Data_type==variable_name]))==FALSE){
    df_leafsites$Sampling_Year[is.na(df_leafsites$Sampling_Year)==TRUE & df_leafsites$exp_nam==exp_name & df_leafsites$Data_type==variable_name] <-max(df_leafsites$Sampling_Year[is.na(df_leafsites$Sampling_Year)==FALSE & df_leafsites$exp_nam==exp_name & df_leafsites$Data_type==variable_name])
  } else {
    df_leafsites$Sampling_Year[is.na(df_leafsites$Sampling_Year)==TRUE & df_leafsites$exp_nam==exp_name & df_leafsites$Data_type==variable_name] <-max(df_leafsites$Sampling_Year[is.na(df_leafsites$Sampling_Year)==FALSE & df_leafsites$exp_nam==exp_name])
  }
}

#only a few missing - no worries!

#now, start working with vcmax data

vcmax_all <- subset(df_leafsites,Data_type=="maximal_Rubisco_activity_(Vcmax)"|Data_type=="Vcmax"|
                      Data_type=="Vcmax,_maximum_carboxylation_rate_of_Rubisco")
vcmax25_all <- subset(df_leafsites,Data_type=="Vcmax25") # all unit is ok

Jmax_all <- subset(df_leafsites,Data_type=="Jmax_"|Data_type=="Jmax"|
                      Data_type=="Jmax,_maximum_rate_of_electron_transport_for_regeneration_of_RuBP")

Jmax25_all <- subset(df_leafsites,Data_type=="Jmax25" & Unit!="g/m2y") # remove one 'outlier' unit

Asat_all <- df_leafsites[grep("Asat", df_leafsites$Data_type),]

Amax_all <- df_leafsites[grep("Amax", df_leafsites$Data_type),]

#now, convert unit separately
#this three value is consistent, but its value should be umol/m2/s. Let's use this as final data of Asat
summary(subset(Asat_all,Unit=="?mol_CO2_m?2_s?1")$ambient)
summary(subset(Asat_all,Unit=="mol_CO2_/_m2_s")$ambient)
summary(subset(Asat_all,Unit=="mol/m2*s")$ambient)

Asat_all_final <- subset(Asat_all,Unit=="?mol_CO2_m?2_s?1"|
                           Unit=="mol_CO2_/_m2_s"|
                           Unit=="mol/m2*s")

Amax_all_final <- subset(Asat_all,Unit=="?mol_CO2_m?2_s?1"|
                           Unit=="mol_/_m2_s"|
                           Unit=="mol/m2*s")

vcmax_all %>% group_by(Unit)  %>% summarise(mean = mean(ambient,na.rm=TRUE))
summary(subset(vcmax_all,Unit=="mol/m2_x_s")$ambient)
summary(subset(vcmax_all,is.na(Unit)==TRUE)$ambient)

#I think all vcmax is umol/m2/s, in no doubt! So does Jmax
vcmax_all_final <- vcmax_all
vcmax25_all_final <- vcmax25_all

Jmax_all_final <- Jmax_all
Jmax25_all_final <- Jmax25_all

#one point shows start and sampling year was 1998 and 1997 separately, wrong. Convert to same year
Asat_all_final$Sampling_Year[Asat_all_final$Sampling_Year<Asat_all_final$Start_Year] <- 1998
Amax_all_final$Sampling_Year[Amax_all_final$Sampling_Year<Amax_all_final$Start_Year] <- 1998

Asat_all_final$Year <- Asat_all_final$Sampling_Year-Asat_all_final$Start_Year
summary(Asat_all_final$Year)
Amax_all_final$Year <- Amax_all_final$Sampling_Year-Amax_all_final$Start_Year
summary(Amax_all_final$Year)

vcmax_all_final$Start_Year[is.na(vcmax_all_final$Start_Year)==TRUE] <- 1995 # one place was missing - interpolate as equal to sampling_year
vcmax_all_final$Year <- vcmax_all_final$Sampling_Year-vcmax_all_final$Start_Year
summary(vcmax_all_final$Year)

Jmax_all_final$Year <- Jmax_all_final$Sampling_Year-Jmax_all_final$Start_Year
summary(Jmax_all_final$Year)

vcmax25_all_final$Year <- vcmax25_all_final$Sampling_Year-vcmax25_all_final$Start_Year
summary(vcmax25_all_final$Year)

Jmax25_all_final$Start_Year[is.na(Jmax25_all_final$Start_Year)==TRUE]  <- 1995 # one place was missing - interpolate as equal to sampling_year
Jmax25_all_final$Year <- Jmax25_all_final$Sampling_Year-Jmax25_all_final$Start_Year
summary(Jmax25_all_final$Year)


#now, all done! start to do one-point method - see gwr_for_climates_Asat.R
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/Asat_all_df.csv")
write.csv(Asat_all_final, csvfile, row.names = TRUE)

csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/Amax_all_df.csv")
write.csv(Amax_all_final, csvfile, row.names = TRUE)

csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/vcmax_all_df.csv")
write.csv(vcmax_all_final, csvfile, row.names = TRUE)

csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/jmax_all_df.csv")
write.csv(Jmax_all_final, csvfile, row.names = TRUE)

csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/vcmax25_all_df.csv")
write.csv(vcmax25_all_final, csvfile, row.names = TRUE)

csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/jmax25_all_df.csv")
write.csv(Jmax25_all_final, csvfile, row.names = TRUE)

#now, after gwr_for_climates_Asat.R by collecting all climate data for vcmax, jmax, Asat, Amax...
#basing on proposal by below priority
#1. vcmax25
#2. vcmax and then converting to vcmax25 by Tg
#3. Asat and then converting to vcmax by one-point method, and convert to vcmax25 by Tg
#so does jmax..

#1.vcmax25
summary(vcmax25_all_final)
vcmax25_all_final%>% group_by(exp_nam,treatment)  %>% summarise(number = n())
#Brandbjerg, c, d, w --> 7 types

#2.vcmax
summary(vcmax_all_final)
vcmax_all_final%>% group_by(exp_nam,treatment)  %>% summarise(number = n())
#Brandbjerg, c, d, w --> 7 types but more samples --> so still use vcmax rather than vcmax25 in Brandbjerg
#now, include Tg
tmn_df <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/climates/vcmax_tmn.csv")
tmx_df <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/climates/vcmax_tmx.csv")
lat <- tmn_df$lat
tmx <- tmx_df[,c(8:ncol(tmx_df))]
tmn <- tmn_df[,c(8:ncol(tmn_df))]
s1 <- (-23.1+ -17.3)/2
s2 <- (-17.3 + -8)/2
s3 <- (-8 + 4.1)/2
s4 <- (4.1 + 14.8)/2
s5 <- (14.8 + 21.9)/2
s6 <- (21.9 + 23.2)/2
s7 <- (23.2 + 18.3)/2
s8 <- (18.3 + 8.6)/2
s9 <-  (8.6 + -2.8)/2
s10 <- (-2.8 + -14.1)/2
s11 <- (-14.1 + -21.6)/2
s12 <- (-21.6 + -23.1)/2

s <- c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)
xx <- data.frame(matrix(nrow=nrow(tmx), ncol=ncol(tmx)))
output_Tg <- data.frame(matrix(nrow=nrow(tmx), ncol=ncol(tmx)))
#xx = acos(h), h = hour angle of the sun
for (a in 1:12){ 
  month_no <- seq(from = 1, to = ncol(tmx), by = 12)+a-1
  xx[1:nrow(tmx),month_no] <- -tan(pi*lat/180)*tan(s[a]*pi/180)
}

#check each part of Tg formula
part1 <- (0.5+((1-xx^2)^(0.5))/(2*acos(xx)))[,1:12]
part2 <- (0.5-((1-xx^2)^(0.5))/(2*acos(xx)))[,1:12]
summary(part1)
summary(part2)

#the percentage of tmx was dominated overall
Tg_site <- tmx*(0.5+((1-xx^2)^(0.5))/(2*acos(xx)))+ tmn*(0.5-((1-xx^2)^(0.5))/(2*acos(xx)))
Tg_site[Tg_site =="NaN"] <- NA
Tg_site[Tg_site < 0] <- NA

tmn_df$Tg <- rowMeans(Tg_site,na.rm=TRUE)
vcmax_Tg <- tmn_df[,c("lon","lat","z","year_start","year_end","Tg")]

#keep the design same as forcing
vcmax_all_final$year_start <- vcmax_all_final$Start_Year
vcmax_all_final$year_end <- vcmax_all_final$Sampling_Year
vcmax_all_final$year_start[vcmax_all_final$year_end>2016] <- 2007
vcmax_all_final$year_end[vcmax_all_final$year_end>2016] <- 2016

vcmax_all_final_Tg <- merge(vcmax_all_final,vcmax_Tg,by=c("lon","lat","z","year_start","year_end"),all.x=TRUE)
summary(vcmax_all_final_Tg)



#now, calculate response ratio directly
vcmax_all_final_Tg$vcmax25_a <- vcmax_all_final_Tg$ambient*exp((65330/8.314)*((1/(vcmax_all_final_Tg$Tg+273.15))-(1/298.15)))
vcmax_all_final_Tg$vcmax25_e <- vcmax_all_final_Tg$elevated*exp((65330/8.314)*((1/(vcmax_all_final_Tg$Tg+273.15))-(1/298.15)))
vcmax_all_final_Tg$response_ratio <- log(vcmax_all_final_Tg$vcmax25_e/vcmax_all_final_Tg$vcmax25_a)

#use alternative way (kattge and Knorr, 2007) - assume Tcgrowth = Tleaf = Tg
tcgrowth <- vcmax_all_final_Tg$Tg + 273.15
tkleaf <- vcmax_all_final_Tg$Tg + 273.15
dent <-  668.39-1.07*tcgrowth
fva <- exp( 71513 * (tkleaf - 298.15) / (298.15 * 8.3145 * tkleaf) )
fvb <- (1.0 + exp( (298.15 * dent - 200000)/(8.3145 * 298.15) ) ) / (1.0 + exp( (tkleaf * dent - 200000)/(8.3145 * tkleaf) ) )
fv  <- fva * fvb
vcmax_all_final_Tg$vcmax25_a_inst <- vcmax_all_final_Tg$ambient/fv
vcmax_all_final_Tg$vcmax25_e_inst <- vcmax_all_final_Tg$elevated/fv
vcmax_all_final_Tg$response_ratio_inst <- log(vcmax_all_final_Tg$vcmax25_e_inst/vcmax_all_final_Tg$vcmax25_a_inst)
vcmax_all_final_Tg$vcmax_a <- vcmax_all_final_Tg$ambient
vcmax_all_final_Tg$vcmax_e <- vcmax_all_final_Tg$elevated
#nothing changed here...just update the formula..
hist(vcmax_all_final_Tg$response_ratio)

#now, decide which rest of Asat will be used then!
vcmax_expnam <- vcmax_all_final_Tg%>% group_by(exp_nam,treatment)  %>% summarise(number = n())
names(vcmax_expnam) <- c("exp_nam","treatment","number_vcmax")
Asat_expnam <- Asat_all_final%>% group_by(exp_nam,treatment)  %>% summarise(number = n())
names(Asat_expnam) <- c("exp_nam","treatment","number_Asat")

Asat_need <- merge(Asat_expnam,vcmax_expnam,by=c("exp_nam","treatment"),all.x=TRUE)
#now, create new definition - if having vcmax, even if its number is much less than Asat.
Asat_need$method[is.na(Asat_need$number_vcmax)==TRUE] <- "Asat"
Asat_need$method[is.na(Asat_need$number_vcmax)==FALSE] <- "vcmax"

Asat_no_include_info <- Asat_need[,c("exp_nam","treatment","method")]

Asat_all_final_removal <- merge(Asat_all_final,Asat_no_include_info,by=c("exp_nam","treatment"),all.x=TRUE)
Asat_all_final_removal_yes <- subset(Asat_all_final_removal,method!="vcmax")
dim(Asat_all_final_removal)
dim(Asat_all_final_removal_yes) # 327 of 807 samples were selected, because rest of them was already included in vcmax dataset

#now, start one-point method
tmn_df <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/climates/Asat_tmn.csv")
tmx_df <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/climates/Asat_tmx.csv")
qair_df <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/climates/Asat_qair.csv")

#co2 
Ca_df <- tmn_df[,c("lon","lat","z","start_yr","end_yr")]
Ca <- read.csv("/Users/yunpeng/data/co2/cCO2_rcp85_const850-1765.csv")
for (i in 1:nrow(Ca_df)){
  Ca_df$Ca_start[i] <- Ca$co2[Ca_df$start_yr[i]==Ca$year]
  Ca_df$Ca_end[i] <- Ca$co2[Ca_df$end_yr[i]==Ca$year]
  Ca_df$Ca_Global[i] <- mean(c(Ca_df$Ca_start[i],Ca_df$Ca_end[i]))
}
summary(Ca_df)


tmx <- tmx_df[,c(8:ncol(tmx_df))]
tmn <- tmn_df[,c(8:ncol(tmn_df))]
qair <- qair_df[,c(8:ncol(qair_df))]
CaGlobal <- Ca_df$Ca_Global
lat <- Ca_df$lat # degree
elev <- Ca_df$z #m

##(1) Tg (growth temperature)
#1. Tg
#solar declination from Jan to Dec
s1 <- (-23.1+ -17.3)/2
s2 <- (-17.3 + -8)/2
s3 <- (-8 + 4.1)/2
s4 <- (4.1 + 14.8)/2
s5 <- (14.8 + 21.9)/2
s6 <- (21.9 + 23.2)/2
s7 <- (23.2 + 18.3)/2
s8 <- (18.3 + 8.6)/2
s9 <-  (8.6 + -2.8)/2
s10 <- (-2.8 + -14.1)/2
s11 <- (-14.1 + -21.6)/2
s12 <- (-21.6 + -23.1)/2

s <- c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)
xx <- data.frame(matrix(nrow=nrow(Ca_df), ncol=ncol(tmx)))
output_Tg <- data.frame(matrix(nrow=nrow(Ca_df), ncol=ncol(tmx)))
#xx = acos(h), h = hour angle of the sun
for (a in 1:12){ 
  month_no <- seq(from = 1, to = ncol(tmx), by = 12)+a-1
  xx[1:nrow(Ca_df),month_no] <- -tan(pi*lat/180)*tan(s[a]*pi/180)
}

#check each part of Tg formula
part1 <- (0.5+((1-xx^2)^(0.5))/(2*acos(xx)))[,1:12]
part2 <- (0.5-((1-xx^2)^(0.5))/(2*acos(xx)))[,1:12]
summary(part1)
summary(part2)

#the percentage of tmx was dominated overall
Tg_site <- tmx*(0.5+((1-xx^2)^(0.5))/(2*acos(xx)))+ tmn*(0.5-((1-xx^2)^(0.5))/(2*acos(xx)))
Tg_site[Tg_site =="NaN"] <- NA
Tg_site[Tg_site < 0] <- NA

Tg_K <- Tg_site +273.15

##(2) Pressure, Ca, vpd
pressure <- 101325*((1-0.0065*elev/288.15)^(9.80665*0.028963/8.31447/0.0065)) #pressure in Pa -> 1
Ca <- pressure*CaGlobal/1000000 ##Ca in Pa -> 1

# (2.2) now, adding a new function of vpd
kR  = 8.31446262  # universal gas constant, J/mol/K (Allen, 1973)
kMv = 18.02                 # molecular weight of water vapor, g/mol (Tsilingiris, 2008)
kMa = 28.963                #molecular weight of dry air, g/mol (Tsilingiris, 2008) XXX this was in SPLASH (WITH 1E-3 IN EQUATION) XXX
wair <- qair / ( 1 - qair )
rv <- kR / kMv
rd <- kR / kMa
patm <- pressure # Pa
ea <- (patm * wair * rv / (rd + wair * rv))/1000 #ea in Pa
TT_K <- ((tmx + tmn)/2) + 273.15 + Tg_K - Tg_K # averaged from tmx, tmn, then from growth season only
vpd <- 0.611*exp(17.27*(TT_K-273.15)/((TT_K-273.15)+237.3)) - ea         #vpd in KPa -> month/daily
vpd[vpd < 0] <- NA

##(3) Gstar, Kc, Oc, Ko, K, f1
Gstar <- pressure*4.22*exp((37830/8.31447)*(1/298.15 - 1/(Tg_K)))/101325    #Gstar in Pa -> month/daily
Kc <- 39.97*exp((79430/8.31447)*(1/298.15 - 1/(Tg_K)))                      #Kc in Pa -> 1
Oc <-0.21*pressure                                                              #Oc in Pa -> month/daily
Ko <-27480*exp((36380/8.31447)*(1/298.15 - 1/(Tg_K)))                       #Ko in Pa -> month/daily
K  <- Kc*(1+(Oc/Ko))                                                        #K in Pa -> month/daily
f1 <- exp((-580/((-138+Tg_K)^2))*(Tg_K-298.15))                             
E <-sqrt((146*(K+Gstar))/(1.6*f1))                                         

##(4) X, Ci, m
X <- (Gstar/Ca)+ (1-(Gstar/Ca))*E/(E+sqrt(vpd*1000))                        #predicted X -> monthly/daily
Ci <- X*Ca

Ca_df$Tg <- rowMeans(Tg_site,na.rm = TRUE)
Ca_df$Ci <- rowMeans(Ci,na.rm = TRUE)
Ca_df$Gstar <- rowMeans(Gstar,na.rm = TRUE)
Ca_df$K <- rowMeans(K,na.rm = TRUE)
Ca_df$year_start <- Ca_df$start_yr 
Ca_df$year_end <- Ca_df$end_yr 

Ca_df <- Ca_df[,c("lon","lat","z","year_start","year_end","Tg","Ci","Gstar","K")]
#now, merge into df and combine
Asat_all_final_removal_yes$year_start <- Asat_all_final_removal_yes$Start_Year
Asat_all_final_removal_yes$year_end <- Asat_all_final_removal_yes$Sampling_Year
Asat_all_final_removal_yes$year_start[Asat_all_final_removal_yes$year_end>2016] <- 2007
Asat_all_final_removal_yes$year_end[Asat_all_final_removal_yes$year_end>2016] <- 2016
summary(Asat_all_final_removal_yes)
Asat_all_final_removal_yes_onepoint <- merge(Asat_all_final_removal_yes,Ca_df,by=c("lon","lat","z","year_start","year_end"),all.x=TRUE)

Asat_all_final_removal_yes_onepoint$vcmax_a <-Asat_all_final_removal_yes_onepoint$ambient/((Asat_all_final_removal_yes_onepoint$Ci-Asat_all_final_removal_yes_onepoint$Gstar)/(Asat_all_final_removal_yes_onepoint$Ci+Asat_all_final_removal_yes_onepoint$K)-0.015)
Asat_all_final_removal_yes_onepoint$vcmax25_a <- Asat_all_final_removal_yes_onepoint$vcmax_a*exp((65330/8.314)*((1/(273.15+Asat_all_final_removal_yes_onepoint$Tg))-(1/298.15)))

Asat_all_final_removal_yes_onepoint$vcmax_e <-Asat_all_final_removal_yes_onepoint$elevated/((Asat_all_final_removal_yes_onepoint$Ci-Asat_all_final_removal_yes_onepoint$Gstar)/(Asat_all_final_removal_yes_onepoint$Ci+Asat_all_final_removal_yes_onepoint$K)-0.015)
Asat_all_final_removal_yes_onepoint$vcmax25_e <- Asat_all_final_removal_yes_onepoint$elevated*exp((65330/8.314)*((1/(273.15+Asat_all_final_removal_yes_onepoint$Tg))-(1/298.15)))

Asat_all_final_removal_yes_onepoint$response_ratio <- log(Asat_all_final_removal_yes_onepoint$vcmax25_e/Asat_all_final_removal_yes_onepoint$vcmax25_a)

#now, final combination
names(Asat_all_final_removal_yes_onepoint)
names(vcmax_all_final_Tg)

vcmax_all_final_Tg$method <- "vcmax"
Asat_all_final_removal_yes_onepoint$method <- "asat"

vcmax_all_final_Tg%>% group_by(exp_nam)  %>% summarise(number = n())
Asat_all_final_removal_yes_onepoint%>% group_by(exp_nam)  %>% summarise(number = n())

vcmax25_final <- dplyr::bind_rows(vcmax_all_final_Tg,Asat_all_final_removal_yes_onepoint)

#finally, check unique basing on exp_nam, prev_name,ambient,elevated - removing them
vcmax25_final2 <- vcmax25_final[!duplicated(vcmax25_final[c("exp_nam","prev_name","ambient","elevated")]),]

#now, firstly look at co2 compulsory effect
vcmax25_final2%>% group_by(treatment)  %>% summarise(number = n())
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="c"] <- "co2"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="cd"] <- "co2 + drought"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="cf"] <- "co2 + fertilization"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="ci"] <- "co2 + irrigation"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="cw"] <- "co2+ 1-degree warming"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="cw3"] <- "co2+ 3-degree warming"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="cwd"] <- "co2 + warming + drought"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="d"] <- "drought"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="f"] <- "fertilization"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="i"] <- "irrigation"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="w"] <- "1-degree warming"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="w2"] <- " 2-degree warming "
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="wd"] <- "warming + drought"

co2_only <- vcmax25_final2[grep("c", vcmax25_final2$treatment),]
a1 <- ggplot(co2_only, aes(x=treatment_label, y=response_ratio)) +geom_jitter()+geom_boxplot(alpha=0.5,color="red")+
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))
  
#ggsave(paste("~/data/output_gcme/test.jpg",sep=""),width = 30, height = 15)

ggplot(vcmax25_final2, aes(x=treatment_label, y=response_ratio)) +geom_jitter()+geom_boxplot(alpha=0.5,color="red")+
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of vcmax25", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()
ggsave(paste("~/data/output_gcme/test1.jpg",sep=""),width = 30, height = 15)

ggplot(subset(vcmax25_final2,treatment=="f"), aes(x=exp_nam, y=response_ratio)) +geom_jitter()+geom_boxplot(alpha=0.5,color="red")+
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of vcmax25", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()

#aggregate by site
vcmax25_final2_expnam <- aggregate(vcmax25_final2,by=list(vcmax25_final2$exp_nam,vcmax25_final2$treatment_label), FUN=mean, na.rm=TRUE)
vcmax25_final2_expnam
ggplot(vcmax25_final2_expnam, aes(x=Group.2, y=response_ratio)) +geom_jitter()+geom_boxplot(alpha=0.5,color="red")+
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of vcmax25", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()
ggsave(paste("~/data/output_gcme/test2.jpg",sep=""),width = 30, height = 15)



#now, jmax!!
#1.jmax25
summary(Jmax25_all_final)
Jmax25_all_final%>% group_by(exp_nam,treatment)  %>% summarise(number = n())
#3 sites, 12 types

#2.jmax
summary(Jmax_all_final)
Jmax_all_final%>% group_by(exp_nam,treatment)  %>% summarise(number = n())
#lots of samples than jmax25 - so use jmax at the end
#now, include Tg
tmn_df <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/climates/jmax_tmn.csv")
tmx_df <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/climates/jmax_tmx.csv")
lat <- tmn_df$lat
tmx <- tmx_df[,c(8:ncol(tmx_df))]
tmn <- tmn_df[,c(8:ncol(tmn_df))]
s1 <- (-23.1+ -17.3)/2
s2 <- (-17.3 + -8)/2
s3 <- (-8 + 4.1)/2
s4 <- (4.1 + 14.8)/2
s5 <- (14.8 + 21.9)/2
s6 <- (21.9 + 23.2)/2
s7 <- (23.2 + 18.3)/2
s8 <- (18.3 + 8.6)/2
s9 <-  (8.6 + -2.8)/2
s10 <- (-2.8 + -14.1)/2
s11 <- (-14.1 + -21.6)/2
s12 <- (-21.6 + -23.1)/2

s <- c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)
xx <- data.frame(matrix(nrow=nrow(tmx), ncol=ncol(tmx)))
output_Tg <- data.frame(matrix(nrow=nrow(tmx), ncol=ncol(tmx)))
#xx = acos(h), h = hour angle of the sun
for (a in 1:12){ 
  month_no <- seq(from = 1, to = ncol(tmx), by = 12)+a-1
  xx[1:nrow(tmx),month_no] <- -tan(pi*lat/180)*tan(s[a]*pi/180)
}

#check each part of Tg formula
part1 <- (0.5+((1-xx^2)^(0.5))/(2*acos(xx)))[,1:12]
part2 <- (0.5-((1-xx^2)^(0.5))/(2*acos(xx)))[,1:12]
summary(part1)
summary(part2)

#the percentage of tmx was dominated overall
Tg_site <- tmx*(0.5+((1-xx^2)^(0.5))/(2*acos(xx)))+ tmn*(0.5-((1-xx^2)^(0.5))/(2*acos(xx)))
Tg_site[Tg_site =="NaN"] <- NA
Tg_site[Tg_site < 0] <- NA

tmn_df$Tg <- rowMeans(Tg_site,na.rm=TRUE)
Jmax_Tg <- tmn_df[,c("lon","lat","z","year_start","year_end","Tg")]

#keep the design same as forcing
Jmax_all_final$year_start <- Jmax_all_final$Start_Year
Jmax_all_final$year_end <- Jmax_all_final$Sampling_Year
Jmax_all_final$year_start[Jmax_all_final$year_end>2016] <- 2007
Jmax_all_final$year_end[Jmax_all_final$year_end>2016] <- 2016

Jmax_all_final_Tg <- merge(Jmax_all_final,Jmax_Tg,by=c("lon","lat","z","year_start","year_end"),all.x=TRUE)
summary(Jmax_all_final_Tg)

#now, calculate response ratio directly
Jmax_all_final_Tg$Jmax25_a <- Jmax_all_final_Tg$ambient*exp((43540/8.314)*((1/(Jmax_all_final_Tg$Tg+273.15))-(1/298.15)))
Jmax_all_final_Tg$Jmax25_e <- Jmax_all_final_Tg$elevated*exp((43540/8.314)*((1/(Jmax_all_final_Tg$Tg+273.15))-(1/298.15)))
Jmax_all_final_Tg$response_ratio <- log(Jmax_all_final_Tg$Jmax25_e/Jmax_all_final_Tg$Jmax25_a)

Jmax_all_final_Tg$Jmax_a <- Jmax_all_final_Tg$ambient
Jmax_all_final_Tg$Jmax_e <- Jmax_all_final_Tg$elevated

hist(Jmax_all_final_Tg$response_ratio)


#now, decide which rest of Asat will be used then!
Jmax_expnam <- Jmax_all_final_Tg%>% group_by(exp_nam,treatment)  %>% summarise(number = n())
names(Jmax_expnam) <- c("exp_nam","treatment","number_Jmax")
Amax_expnam <- Amax_all_final%>% group_by(exp_nam,treatment)  %>% summarise(number = n())
names(Amax_expnam) <- c("exp_nam","treatment","number_Amax")

Amax_need <- merge(Amax_expnam,Jmax_expnam,by=c("exp_nam","treatment"),all.x=TRUE)
#now, create new definition - if having vcmax, even if its number is much less than Asat.
Amax_need$method[is.na(Amax_need$number_Jmax)==TRUE] <- "Amax"
Amax_need$method[is.na(Amax_need$number_Jmax)==FALSE] <- "Jmax"

Amax_no_include_info <- Amax_need[,c("exp_nam","treatment","method")]

Amax_all_final_removal <- merge(Amax_all_final,Amax_no_include_info,by=c("exp_nam","treatment"),all.x=TRUE)
Amax_all_final_removal_yes <- subset(Amax_all_final_removal,method!="Jmax")
dim(Amax_all_final_removal)
dim(Amax_all_final_removal_yes) # 314 of 771 samples were selected, because rest of them was already included in Jmax dataset
#now, include Tg
tmn_df <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/climates/Amax_tmn.csv")
tmx_df <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/climates/Amax_tmx.csv")
lat <- tmn_df$lat
tmx <- tmx_df[,c(8:ncol(tmx_df))]
tmn <- tmn_df[,c(8:ncol(tmn_df))]
s1 <- (-23.1+ -17.3)/2
s2 <- (-17.3 + -8)/2
s3 <- (-8 + 4.1)/2
s4 <- (4.1 + 14.8)/2
s5 <- (14.8 + 21.9)/2
s6 <- (21.9 + 23.2)/2
s7 <- (23.2 + 18.3)/2
s8 <- (18.3 + 8.6)/2
s9 <-  (8.6 + -2.8)/2
s10 <- (-2.8 + -14.1)/2
s11 <- (-14.1 + -21.6)/2
s12 <- (-21.6 + -23.1)/2

s <- c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)
xx <- data.frame(matrix(nrow=nrow(tmx), ncol=ncol(tmx)))
output_Tg <- data.frame(matrix(nrow=nrow(tmx), ncol=ncol(tmx)))
#xx = acos(h), h = hour angle of the sun
for (a in 1:12){ 
  month_no <- seq(from = 1, to = ncol(tmx), by = 12)+a-1
  xx[1:nrow(tmx),month_no] <- -tan(pi*lat/180)*tan(s[a]*pi/180)
}

#check each part of Tg formula
part1 <- (0.5+((1-xx^2)^(0.5))/(2*acos(xx)))[,1:12]
part2 <- (0.5-((1-xx^2)^(0.5))/(2*acos(xx)))[,1:12]
summary(part1)
summary(part2)

#the percentage of tmx was dominated overall
Tg_site <- tmx*(0.5+((1-xx^2)^(0.5))/(2*acos(xx)))+ tmn*(0.5-((1-xx^2)^(0.5))/(2*acos(xx)))
Tg_site[Tg_site =="NaN"] <- NA
Tg_site[Tg_site < 0] <- NA

tmn_df$Tg <- rowMeans(Tg_site,na.rm=TRUE)
Amax_Tg <- tmn_df[,c("lon","lat","z","year_start","year_end","Tg")]

#keep the design same as forcing
Amax_all_final_removal_yes$year_start <- Amax_all_final_removal_yes$Start_Year
Amax_all_final_removal_yes$year_end <- Amax_all_final_removal_yes$Sampling_Year
Amax_all_final_removal_yes$year_start[Amax_all_final_removal_yes$year_end>2016] <- 2007
Amax_all_final_removal_yes$year_end[Amax_all_final_removal_yes$year_end>2016] <- 2016

Amax_all_final_Tg <- merge(Amax_all_final_removal_yes,Amax_Tg,by=c("lon","lat","z","year_start","year_end"),all.x=TRUE)
summary(Amax_all_final_Tg)

#now, calculate response ratio directly Amax * 4 = Jmax!!
Amax_all_final_Tg$Jmax25_a <- Amax_all_final_Tg$ambient*4*exp((43540/8.314)*((1/(Amax_all_final_Tg$Tg+273.15))-(1/298.15)))
Amax_all_final_Tg$Jmax25_e <- Amax_all_final_Tg$elevated*4*exp((43540/8.314)*((1/(Amax_all_final_Tg$Tg+273.15))-(1/298.15)))
Amax_all_final_Tg$response_ratio <- log(Amax_all_final_Tg$Jmax25_e/Amax_all_final_Tg$Jmax25_a)

Amax_all_final_Tg$Jmax_a <- Amax_all_final_Tg$ambient*4
Amax_all_final_Tg$Jmax_e <- Amax_all_final_Tg$elevated*4

Jmax_all_final_Tg_removal <- Jmax_all_final_Tg

names(Jmax_all_final_Tg_removal)
names(Amax_all_final_Tg)
Jmax_all_final_Tg_removal$method <- "jmax"
Amax_all_final_Tg$method <- "amax"

Jmax25_final <- dplyr::bind_rows(Jmax_all_final_Tg_removal,Amax_all_final_Tg)
summary(Jmax25_final)


#finally, check unique basing on exp_nam, prev_name,ambient,elevated - removing them
Jmax25_final2 <- Jmax25_final[!duplicated(Jmax25_final[c("exp_nam","prev_name","ambient","elevated")]),]
dim(Jmax25_final2)-dim(Jmax25_final)
#now, firstly look at co2 compulsory effect
Jmax25_final2%>% group_by(treatment)  %>% summarise(number = n())
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="c"] <- "co2"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="cd"] <- "co2 + drought"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="cf"] <- "co2 + fertilization"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="ci"] <- "co2 + irrigation"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="cw"] <- "co2+ 1-degree warming"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="cw3"] <- "co2+ 3-degree warming"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="cwd"] <- "co2 + warming + drought"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="d"] <- "drought"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="f"] <- "fertilization"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="i"] <- "irrigation"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="w"] <- "1-degree warming"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="w2"] <- " 2-degree warming "
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="wd"] <- "warming + drought"

ggplot(Jmax25_final2, aes(x=treatment_label, y=response_ratio)) +geom_jitter()+geom_boxplot(alpha=0.5,color="red")+
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of Jmax25", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()
ggsave(paste("~/data/output_gcme/test3.jpg",sep=""),width = 30, height = 15)

ggplot(subset(Jmax25_final2,treatment=="f"), aes(x=exp_nam, y=response_ratio)) +geom_jitter()+geom_boxplot(alpha=0.5,color="red")+
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of vcmax25", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()

#aggregate by site
Jmax25_final2_expnam <- aggregate(Jmax25_final2,by=list(Jmax25_final2$exp_nam,Jmax25_final2$treatment_label), FUN=mean, na.rm=TRUE)
Jmax25_final2_expnam
ggplot(Jmax25_final2_expnam, aes(x=Group.2, y=response_ratio)) +geom_jitter()+geom_boxplot(alpha=0.5,color="red")+
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of Jmax25", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()
ggsave(paste("~/data/output_gcme/test4.jpg",sep=""),width = 30, height = 15)


#attempt:
#jmax only (after removal)
ggplot(subset(Jmax25_final2,method=="jmax"), aes(x=treatment_label, y=response_ratio)) +geom_jitter()+geom_boxplot(alpha=0.5,color="red")+
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of Jmax25", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()
#Amax only (after removal)
ggplot(subset(Jmax25_final2,method=="amax"), aes(x=treatment_label, y=response_ratio)) +geom_jitter()+geom_boxplot(alpha=0.5,color="red")+
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of Jmax25", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()

#include all original jmax
ggplot(Jmax_all_final_Tg, aes(x=treatment, y=response_ratio)) +geom_jitter()+geom_boxplot(alpha=0.5,color="red")+
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of Jmax25", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()


#now, start forcing
#forcing R file was available at ~/yunkepeng/gcme/vcmax_forcing.R
#forcing output was available in /Users/yunpeng/data/gcme/vcmax_forcing/
#siteinfo was available in /Users/yunpeng/data/gcme/yunke_coord/forcing_vcmax.csv
summary(vcmax25_final2)
summary(Jmax25_final2)
#c,cw,cw3,w,w2 were firstly selected to test its response in pmodel
vcmax25_warmingco2 <- subset(vcmax25_final2,treatment=="c"|treatment=="cw"|treatment=="cw3"|treatment=="w"|treatment=="w2")
Jmax25_warmingco2 <- subset(Jmax25_final2,treatment=="c"|treatment=="cw"|treatment=="cw3"|treatment=="w"|treatment=="w2")

subset(vcmax25_final2,treatment=="c"|treatment=="cw"|treatment=="cw3") %>% group_by(co2_a,co2_e)  %>% summarise(number = n())
subset(Jmax25_warmingco2,treatment=="c"|treatment=="cw"|treatment=="cw3") %>% group_by(co2_a,co2_e)  %>% summarise(number = n())

subset(vcmax25_final2,treatment=="w"|treatment=="w2"|treatment=="cw"|treatment=="cw3") %>% group_by(warm_a,warm_e,warmQ_e)  %>% summarise(number = n())
subset(Jmax25_warmingco2,treatment=="w"|treatment=="w2"|treatment=="cw"|treatment=="cw3") %>% group_by(warm_a,warm_e,warmQ_e)  %>% summarise(number = n())
#using warm_e, about how many temperature increased

#combine with siteinfo
siteinfo <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/forcing_vcmax.csv")
siteinfo <- siteinfo[,c("lon","lat","elv","year_start","year_end","sitename")]
names(siteinfo) <- c("lon","lat","z","year_start","year_end","sitename")
vcmax25_warmingco2_siteinfo <- merge(vcmax25_warmingco2,siteinfo,by=c("lon","lat","z","year_start","year_end"),all.x=TRUE)
Jmax25_warmingco2_siteinfo <- merge(Jmax25_warmingco2,siteinfo,by=c("lon","lat","z","year_start","year_end"),all.x=TRUE)

#now, start p-model
#devtools::load_all("~/yunkepeng/gcme/pmodel/ingestr/")

#default inst_vcmax
devtools::load_all("/Users/yunpeng/yunkepeng/gcme/pmodel/rsofun/")
#changing to use adjusted parameters for inst_vcmax from Smith and Keenan 2020 GCB - doesn't change anything!
#devtools::load_all("~/yunkepeng/gcme/pmodel_modified/rsofun/")

df_soiltexture <- bind_rows(
  top    = tibble(layer = "top",    fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1),
  bottom = tibble(layer = "bottom", fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1))
params_modl <- list(
  kphio           = 0.09423773,
  soilm_par_a     = 0.33349283,
  soilm_par_b     = 1.45602286,
  tau_acclim_tempstress = 10,
  par_shape_tempstress  = 0.0)

forcing_path <- "/Users/yunpeng/data/gcme/vcmax_forcing/"
vcmax25_warmingco2_siteinfo$elv <- vcmax25_warmingco2_siteinfo$z


for (i in 1:nrow(vcmax25_warmingco2_siteinfo)){
  siteinfo_site <- vcmax25_warmingco2_siteinfo[i,c("sitename","lon","lat","elv","year_start","year_end")]
  sitename <- siteinfo_site$sitename
    
  df1 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for ambient
  df2 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for elevated
    
  #Q1 - add warming directly
  df2$temp <- df2$temp+vcmax25_warmingco2_siteinfo$warm_e[i]

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
    
  mean_vcmax25_ambient <- mean(modlist1$vcmax25)*1000000
  mean_vcmax25_elevated <- mean(modlist2$vcmax25)*1000000
  mean_vcmax_ambient <- mean(modlist1$vcmax)*1000000
  mean_vcmax_elevated <- mean(modlist2$vcmax)*1000000
  
  #mean_vcmax25_ambient <- max(modlist1$vcmax25)*1000000
  #mean_vcmax25_elevated <- max(modlist2$vcmax25)*1000000
  
  vcmax25_warmingco2_siteinfo$mean_vcmax25_ambient[i] <- mean_vcmax25_ambient
  vcmax25_warmingco2_siteinfo$mean_vcmax25_elevated[i] <- mean_vcmax25_elevated
  vcmax25_warmingco2_siteinfo$mean_vcmax_ambient[i] <- mean_vcmax_ambient
  vcmax25_warmingco2_siteinfo$mean_vcmax_elevated[i] <- mean_vcmax_elevated
  vcmax25_warmingco2_siteinfo$pred_response_ratio[i] <- log(vcmax25_warmingco2_siteinfo$mean_vcmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_vcmax25_ambient[i])
}
summary(vcmax25_warmingco2_siteinfo$pred_response_ratio)

# do the same for jmax
summary(Jmax25_warmingco2_siteinfo)
Jmax25_warmingco2_siteinfo$elv <- Jmax25_warmingco2_siteinfo$z


for (i in 1:nrow(Jmax25_warmingco2_siteinfo)){
  siteinfo_site <- Jmax25_warmingco2_siteinfo[i,c("sitename","lon","lat","elv","year_start","year_end")]
  sitename <- siteinfo_site$sitename
  
  df1 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for ambient
  df2 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for elevated
  
  #Q1 - add warming directly
  df2$temp <- df2$temp+Jmax25_warmingco2_siteinfo$warm_e[i]
  
  #Q2 - if having co2 then replacing it directly - if not then using co2 at measurement year
  if (is.na(Jmax25_warmingco2_siteinfo$co2_a[i])==TRUE|is.na(Jmax25_warmingco2_siteinfo$co2_e[i])==TRUE) { 
    df1$co2 <- df1$co2
    df2$co2 <- df2$co2
  } else {
    df1$co2 <- Jmax25_warmingco2_siteinfo$co2_a[i] # set to ambient 
    df2$co2 <- Jmax25_warmingco2_siteinfo$co2_e[i] # set to elevated 
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
  
  mean_jmax25_ambient <- mean(modlist1$jmax25)*1000000
  mean_jmax25_elevated <- mean(modlist2$jmax25)*1000000
  mean_jmax_ambient <- mean(modlist1$jmax)*1000000
  mean_jmax_elevated <- mean(modlist2$jmax)*1000000
  #mean_jmax25_ambient <- max(modlist1$jmax25)*1000000
  #mean_jmax25_elevated <- max(modlist2$jmax25)*1000000
  
  Jmax25_warmingco2_siteinfo$mean_jmax25_ambient[i] <- mean_jmax25_ambient
  Jmax25_warmingco2_siteinfo$mean_jmax25_elevated[i] <- mean_jmax25_elevated
  Jmax25_warmingco2_siteinfo$mean_jmax_ambient[i] <- mean_jmax_ambient
  Jmax25_warmingco2_siteinfo$mean_jmax_elevated[i] <- mean_jmax_elevated
  Jmax25_warmingco2_siteinfo$pred_response_ratio[i] <- log(Jmax25_warmingco2_siteinfo$mean_jmax25_elevated[i]/Jmax25_warmingco2_siteinfo$mean_jmax25_ambient[i])
}

#now, merge with ECM, AM...
ECM_info <- leaf_all_coord[,c("exp_nam","prev_name","ALIAS","AM","EcM","ErM")]
vcmax25_warmingco2_siteinfo_ECM <- merge(vcmax25_warmingco2_siteinfo,ECM_info,by=c("exp_nam","prev_name","ALIAS"),all.x=TRUE)
#not too many data...

ggplot(subset(vcmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE), aes(x=treatment_label, y=response_ratio)) +geom_jitter()+
  geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=treatment_label, y=pred_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of vcmax25", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()
ggsave(paste("~/data/output_gcme/vcmax25_com.jpg",sep=""),width = 30, height = 15)

aa <- subset(vcmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE)
aaa <- aggregate(aa,by=list(aa$exp_nam,aa$treatment_label), FUN=mean, na.rm=TRUE)
aaa
ggplot(aaa, aes(x=Group.2 , y=response_ratio)) +geom_jitter()+
  geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=Group.2 , y=pred_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of vcmax25", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()


ggplot(subset(vcmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE &method=="vcmax"),
       aes(x=treatment_label, y=response_ratio)) +geom_jitter()+
  geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=treatment_label, y=pred_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="vcmax25 from measured vcmax", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=40),axis.title =element_text(size=20))+coord_flip()

ggsave(paste("~/data/output_gcme/vc25_vcmax.jpg",sep=""),width = 15, height = 15)

#now, summarise temperature  
subset(vcmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE
       &method=="vcmax"&treatment=="cw")%>% group_by(exp_nam)  %>% summarise(number = n())
subset(vcmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE
       &method=="vcmax"&treatment=="w")%>% group_by(exp_nam)  %>% summarise(number = n())


subset(vcmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE &method=="vcmax"&treatment=="w")[,c("exp_nam","year_start","year_end","Sampling_date")]

ggplot(subset(vcmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE &method=="vcmax"&treatment=="w"),
       aes(x=exp_nam, y=response_ratio)) +geom_jitter()+
  geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=exp_nam, y=pred_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="vcmax25: warming only", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()
ggsave(paste("~/data/output_gcme/vc25_vcmax_warming.jpg",sep=""),width = 15, height = 15)


ggplot(subset(vcmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE &method=="vcmax"&treatment=="cw"),
       aes(x=exp_nam, y=response_ratio)) +geom_jitter()+
  geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=exp_nam, y=pred_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="vcmax25: warming +co2", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()
ggsave(paste("~/data/output_gcme/vc25_vcmax_warming_co2.jpg",sep=""),width = 15, height = 15)

ggplot(subset(vcmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE &method=="vcmax"&treatment=="c"),
       aes(x=exp_nam, y=response_ratio)) +geom_jitter()+
  geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=exp_nam, y=pred_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="vcmax25: co2", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()
ggsave(paste("~/data/output_gcme/vc25_vcmax_co2.jpg",sep=""),width = 15, height = 15)

#now, manually merged with ECM/AM plots
naming_expnam <- vcmax25_warmingco2_siteinfo %>% group_by(exp_nam) %>% summarise(number=n()) 
naming_expnam
#csvfile <- paste("~/data/gcme/data_received_190325/co2_plot.csv")
#write.csv(naming_expnam, csvfile, row.names = TRUE)
#newly update sitename, so that it can be merged with ECM/AM info properly
naming_expnam_filled <- read.csv("/Users/yunpeng/data/gcme/data_received_190325/co2_plot_filled.csv")
naming_expnam_filled <- naming_expnam_filled[,c("exp_nam_update","exp_nam")]
names(naming_expnam_filled) <- c("Site.Name","exp_nam")
coord_org <- read.csv("~/data/gcme/data_received_190325/table_var_exp_names.csv")
df_experiments_coord <- merge(coord_org,naming_expnam_filled,by=c("Site.Name"),all.x=TRUE)
df_experiments_coord_co2 <- subset(df_experiments_coord,is.na(exp_nam)==FALSE)
dim(df_experiments_coord_co2)

#final merge
vcmax25_warmingco2_siteinfo_ecm <- merge(vcmax25_warmingco2_siteinfo,df_experiments_coord_co2,by=c("exp_nam"),all.x=TRUE)
vcmax25_warmingco2_siteinfo_ecm%>% group_by(AM)  %>% summarise(number = n())
vcmax25_warmingco2_siteinfo_ecm%>% group_by(ECM)  %>% summarise(number = n())
vcmax25_warmingco2_siteinfo_ecm$ecm_type <- "NA"
vcmax25_warmingco2_siteinfo_ecm$ecm_type[vcmax25_warmingco2_siteinfo_ecm$EcM==TRUE] <- "ECM"
vcmax25_warmingco2_siteinfo_ecm$ecm_type[vcmax25_warmingco2_siteinfo_ecm$AM==TRUE] <- "AM"
vcmax25_warmingco2_siteinfo_ecm$ecm_type[vcmax25_warmingco2_siteinfo_ecm$ErM==TRUE] <- "ErM"
vcmax25_warmingco2_siteinfo_ecm$ecm_type[vcmax25_warmingco2_siteinfo_ecm$Moss==TRUE] <- "Moss"

final_co2 <- subset(vcmax25_warmingco2_siteinfo_ecm,is.na(response_ratio)==FALSE &method=="vcmax"&treatment=="c")

ecm_info <- aggregate(final_co2,by=list(final_co2$exp_nam,final_co2$ecm_type), FUN=mean, na.rm=TRUE)
ecm_info <- ecm_info[,1:2]
names(ecm_info) <- c("exp_nam","ecm_final")
#manually add rest of points basing on Terrer et al. 2016 New Phytol- fig.1
ecm_info$ecm_final[ecm_info$exp_nam=="BioCON"] <-"AM"

final_ecm_co2 <- aggregate(final_co2,by=list(final_co2$exp_nam), FUN=mean, na.rm=TRUE)
final_ecm_co2$exp_nam <- final_ecm_co2$Group.1
final_ecm_co2_info <- merge(final_ecm_co2,ecm_info,
                                   by=c("exp_nam"),all.x=TRUE)

ggplot(final_ecm_co2_info,
       aes(x=response_ratio, y=pred_response_ratio,label=exp_nam,color=ecm_final)) +geom_point() +geom_text(aes(label=exp_nam),hjust=0, vjust=0,size=2)+
  labs(y="predicted response ratio", x = "measured response ratio") +ylim(-0.2,0.2)+xlim(-0.2,0.2) +theme_classic()+geom_abline(slope=1)+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))
ggsave(paste("~/data/output_gcme/vc25_vcmax_co2_ecm_name.jpg",sep=""),width = 15, height = 15)

ggplot(final_ecm_co2_info,
       aes(x=response_ratio, y=pred_response_ratio,color=ecm_final)) +geom_point(size=3) +
  labs(y="predicted response ratio", x = "measured response ratio") +
  ylim(-0.2,0.2)+xlim(-0.2,0.2) +theme_classic()+geom_abline(slope=1)+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))
ggsave(paste("~/data/output_gcme/vc25_vcmax_co2_ecm.jpg",sep=""),width = 15, height = 15)

ggplot(subset(vcmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE &method=="asat"),
       aes(x=treatment_label, y=response_ratio)) +geom_jitter()+
  geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=treatment_label, y=pred_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="vcmax25 from measured Asat", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=30),axis.title =element_text(size=20))+coord_flip()
ggsave(paste("~/data/output_gcme/vc25_asat.jpg",sep=""),width = 15, height = 15)


ggplot(subset(Jmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE), aes(x=treatment_label, y=response_ratio)) +geom_jitter()+
  geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=treatment_label, y=pred_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of Jmax25", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=30),axis.title =element_text(size=30))+coord_flip()
ggsave(paste("~/data/output_gcme/jmax25_com.jpg",sep=""),width = 15, height = 15)


ggplot(subset(Jmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE & method=="jmax"), aes(x=treatment_label, y=response_ratio)) +geom_jitter()+
  geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=treatment_label, y=pred_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Jmax25 from jmax", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=30),axis.title =element_text(size=30))+coord_flip()
ggsave(paste("~/data/output_gcme/jmax25_from_jmax.jpg",sep=""),width = 15, height = 15)


#divide by plots
ggplot(subset(Jmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE & method=="jmax"&treatment=="w"), aes(x=exp_nam, y=response_ratio)) +geom_jitter()+
  geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=exp_nam, y=pred_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of Jmax25", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=30))+coord_flip()
ggsave(paste("~/data/output_gcme/jmax25_from_jmax_warming.jpg",sep=""),width = 15, height = 15)

ggplot(subset(Jmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE & method=="jmax"&treatment=="cw"), aes(x=exp_nam, y=response_ratio)) +geom_jitter()+
  geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=exp_nam, y=pred_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of Jmax25", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()
ggsave(paste("~/data/output_gcme/jmax25_from_jmax_warming_co2.jpg",sep=""),width = 15, height = 15)

ggplot(subset(Jmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE & method=="jmax"&treatment=="c"), aes(x=exp_nam, y=response_ratio)) +geom_jitter()+
  geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=exp_nam, y=pred_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of Jmax25", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()
ggsave(paste("~/data/output_gcme/jmax25_from_jmax_co2.jpg",sep=""),width = 15, height = 15)

final_co2_jmax <- subset(Jmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE & method=="jmax"&treatment=="c")

final_ecm_co2_jmax <- aggregate(final_co2_jmax,by=list(final_co2_jmax$exp_nam), FUN=mean, na.rm=TRUE)
final_ecm_co2_jmax$exp_nam <- final_ecm_co2_jmax$Group.1
final_ecm_co2_jmax_info <- merge(final_ecm_co2_jmax,ecm_info,
                            by=c("exp_nam"),all.x=TRUE)
final_ecm_co2_jmax_info$ecm_final[final_ecm_co2_jmax_info$exp_nam=="RiceFACE_Japan_A_2003_39,38_140,57"] <- "AM"
final_ecm_co2_jmax_info$ecm_final[final_ecm_co2_jmax_info$exp_nam=="RiceFACE_Japan_A_2004_39,38_140,57"] <- "AM"

ggplot(final_ecm_co2_jmax_info,
       aes(x=response_ratio, y=pred_response_ratio,label=exp_nam,color=ecm_final)) +geom_point() +geom_text(aes(label=exp_nam),hjust=0, vjust=0,size=2)+
  labs(y="predicted response ratio jmax25", x = "measured response ratio jmax25") +ylim(-0.2,0.2)+xlim(-0.2,0.2) +theme_classic()+geom_abline(slope=1)+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))
ggsave(paste("~/data/output_gcme/vc25_jmax_co2_ecm_name.jpg",sep=""),width = 15, height = 15)

ggplot(final_ecm_co2_jmax_info,
       aes(x=response_ratio, y=pred_response_ratio,color=ecm_final)) +geom_point() +
  labs(y="predicted response ratio jmax25", x = "measured response ratio jmax25") +ylim(-0.2,0.2)+xlim(-0.2,0.2) +theme_classic()+geom_abline(slope=1)+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))
ggsave(paste("~/data/output_gcme/vc25_jmax_co2_ecm.jpg",sep=""),width = 15, height = 15)


ggplot(subset(Jmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE & method=="amax"), aes(x=treatment_label, y=response_ratio)) +geom_jitter()+
  geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=treatment_label, y=pred_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Jmax25 from amax", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=30),axis.title =element_text(size=30))+coord_flip()
ggsave(paste("~/data/output_gcme/jmax25_from_amax.jpg",sep=""),width = 15, height = 15)



aa <- subset(Jmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE)
aaa <- aggregate(aa,by=list(aa$exp_nam,aa$treatment_label), FUN=mean, na.rm=TRUE)
aaa
ggplot(aaa, aes(x=Group.2 , y=response_ratio)) +geom_jitter()+
  geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=Group.2 , y=pred_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of Jmax25", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()

#now, do something about fertilization exp
summary(vcmax25_final2)
#c,cw,cw3,w,w2 were firstly selected to test its response in pmodel
vcmax25_fertilization <- subset(vcmax25_final2,treatment=="f"|treatment=="cf")
fer_name <- vcmax25_fertilization%>% group_by(exp_nam,treatment)  %>% summarise(number = n())


LAI_all <- subset(df,Data_type=="LAI"|Data_type=="max_LAI"|Data_type=="maximum_LAI")
LAI_all$LAI_response_ratio <- log(LAI_all$elevated/LAI_all$ambient)
LAI_final <-LAI_all[,c("exp_nam","treatment","LAI_response_ratio")]
LAI_final_sitemean <- aggregate(LAI_final,by=list(LAI_final$exp_nam,LAI_final$treatment),
                 FUN=mean, na.rm=TRUE)
LAI_final_sitemean <-LAI_final_sitemean[,c("Group.1","Group.2","LAI_response_ratio")]
names(LAI_final_sitemean) <- c("exp_nam","treatment","LAI_response_ratio")
LAI_final_sitemean

vcmax25_fertilization_LAI <- merge(vcmax25_fertilization,LAI_final_sitemean,
                                   by=c("exp_nam","treatment"),all.x=TRUE)
vcmax25_fertilization_LAI_1 <- aggregate(vcmax25_fertilization_LAI,
                                         by=list(vcmax25_fertilization_LAI$exp_nam,vcmax25_fertilization_LAI$treatment),
                                FUN=mean, na.rm=TRUE)

subset(vcmax25_fertilization_LAI_1,is.na(LAI_response_ratio)==FALSE)$Group.2
length(subset(vcmax25_fertilization_LAI_1,is.na(LAI_response_ratio)==FALSE)$Group.1)
dim(vcmax25_fertilization_LAI_1)

#how about anpp?
ANPP_all <- subset(df,Data_type=="ANPP")
ANPP_all$ANPP_response_ratio <- log(ANPP_all$elevated/ANPP_all$ambient)
ANPP_final <-ANPP_all[,c("exp_nam","treatment","ANPP_response_ratio")]
ANPP_final_sitemean <- aggregate(ANPP_final,by=list(ANPP_final$exp_nam,ANPP_final$treatment),
                                FUN=mean, na.rm=TRUE)
ANPP_final_sitemean <-ANPP_final_sitemean[,c("Group.1","Group.2","ANPP_response_ratio")]
names(ANPP_final_sitemean) <- c("exp_nam","treatment","ANPP_response_ratio")
hist(ANPP_final_sitemean$ANPP_response_ratio)

vcmax25_fertilization_ANPP <- merge(vcmax25_fertilization,ANPP_final_sitemean,
                                   by=c("exp_nam","treatment"),all.x=TRUE)
plot(response_ratio~ANPP_response_ratio,subset(vcmax25_fertilization_ANPP,treatment=="f"))

vcmax25_fertilization_ANPP_1 <- aggregate(vcmax25_fertilization_ANPP,
                                         by=list(vcmax25_fertilization_ANPP$exp_nam,vcmax25_fertilization_ANPP$treatment),
                                         FUN=mean, na.rm=TRUE)
plot(response_ratio~ANPP_response_ratio,vcmax25_fertilization_ANPP_1)


#how about leaf N
#how about anpp?
leaf_N <- subset(df,Data_type=="leaf_N")
unit_leafN <- leaf_N%>% group_by(Unit)  %>% summarise(number = n())

leaf_N$leafN_response_ratio <- log(leaf_N$elevated/leaf_N$ambient)
leaf_N <-leaf_N[,c("exp_nam","treatment","leafN_response_ratio")]
leaf_N_sitemean <- aggregate(leaf_N,by=list(leaf_N$exp_nam,leaf_N$treatment),
                                 FUN=mean, na.rm=TRUE)
dim(leaf_N_sitemean)
leaf_N_sitemean <-leaf_N_sitemean[,c("Group.1","Group.2","leafN_response_ratio")]
names(leaf_N_sitemean) <- c("exp_nam","treatment","leafN_response_ratio")
hist(leaf_N_sitemean$leafN_response_ratio)

vcmax25_final2_leafN <- merge(vcmax25_final2,leaf_N_sitemean,
                                    by=c("exp_nam","treatment"),all.x=TRUE)

ggplot(leaf_N, aes(x=treatment, y=leafN_response_ratio)) +geom_jitter()+geom_boxplot(alpha=0.5,color="red")+
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of leaf N", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()
ggsave(paste("~/data/output_gcme/leafN.jpg",sep=""),width = 15, height = 15)


ggplot(subset(vcmax25_final2_leafN, method=="vcmax"&treatment=="c"),
       aes(x=exp_nam, y=response_ratio)) +
  geom_boxplot(aes(x=exp_nam, y=response_ratio),alpha=0.5,color="red") +
  geom_boxplot(aes(x=exp_nam, y=leafN_response_ratio),alpha=0.5,color="blue") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of vcmax25 and leaf N", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()
ggsave(paste("~/data/output_gcme/vcmax25_leafN_co2.jpg",sep=""),width = 15, height = 15)

ggplot(subset(vcmax25_final2_leafN, method=="vcmax"&treatment=="w"),
       aes(x=exp_nam, y=response_ratio)) +
  geom_boxplot(aes(x=exp_nam, y=response_ratio),alpha=0.5,color="red") +
  geom_boxplot(aes(x=exp_nam, y=leafN_response_ratio),alpha=0.5,color="blue") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of vcmax25 and leaf N", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()
ggsave(paste("~/data/output_gcme/vcmax25_leafN_warming.jpg",sep=""),width = 15, height = 15)

ggplot(subset(vcmax25_final2_leafN, treatment=="f"),
       aes(x=exp_nam, y=response_ratio)) +
  geom_boxplot(aes(x=exp_nam, y=response_ratio),alpha=0.5,color="red") +
  geom_boxplot(aes(x=exp_nam, y=leafN_response_ratio),alpha=0.5,color="blue") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of vcmax25 and leaf N", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()
ggsave(paste("~/data/output_gcme/vcmax25_all_leafN_fertilization.jpg",sep=""),width = 15, height = 15)

#remove Asat
ggplot(subset(vcmax25_final2_leafN, treatment=="f" & method=="vcmax"),
       aes(x=exp_nam, y=response_ratio)) +
  geom_boxplot(aes(x=exp_nam, y=response_ratio),alpha=0.5,color="red") +
  geom_boxplot(aes(x=exp_nam, y=leafN_response_ratio),alpha=0.5,color="blue") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of vcmax25 and leaf N", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()


Jmax25_final2_leafN <- merge(Jmax25_final2,leaf_N_sitemean,
                              by=c("exp_nam","treatment"),all.x=TRUE)

Vcmax25_final2_LAI <- merge(vcmax25_final2,LAI_final_sitemean,
                             by=c("exp_nam","treatment"),all.x=TRUE)
Vcmax25_final2_ANPP <- merge(vcmax25_final2,ANPP_final_sitemean,
                            by=c("exp_nam","treatment"),all.x=TRUE)
ggplot(subset(Vcmax25_final2_ANPP,treatment=="f"),
       aes(x=exp_nam, y=response_ratio)) +
  geom_boxplot(aes(x=exp_nam, y=response_ratio),alpha=0.5,color="red") +
  geom_boxplot(aes(x=exp_nam, y=ANPP_response_ratio),alpha=0.5,color="blue") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of vcmax25 and ANPP under N", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()
ggsave(paste("~/data/output_gcme/vcmax25_all_ANPP_fertilization.jpg",sep=""),width = 15, height = 15)


ggplot(subset(Vcmax25_final2_LAI,treatment=="f"),
       aes(x=exp_nam, y=response_ratio)) +
  geom_boxplot(aes(x=exp_nam, y=response_ratio),alpha=0.5,color="red") +
  geom_boxplot(aes(x=exp_nam, y=LAI_response_ratio),alpha=0.5,color="blue") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of vcmax25 and LAI under N ", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()
ggsave(paste("~/data/output_gcme/vcmax25_all_LAI_fertilization.jpg",sep=""),width = 15, height = 15)

#remove Asat?
ggplot(subset(Vcmax25_final2_LAI,treatment=="f" & method=="vcmax"),
       aes(x=exp_nam, y=response_ratio)) +
  geom_boxplot(aes(x=exp_nam, y=response_ratio),alpha=0.5,color="red") +
  geom_boxplot(aes(x=exp_nam, y=LAI_response_ratio),alpha=0.5,color="blue") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of vcmax25 and LAI under N ", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()


#how about leaf area
LA <- subset(df,Data_type=="leaf_area"|Data_type=="leaf_area_")
LA$LA_response_ratio <- log(LA$elevated/LA$ambient)
LA <-LA[,c("exp_nam","treatment","LA_response_ratio")]
LA_sitemean <- aggregate(LA,by=list(LA$exp_nam,LA$treatment),
                             FUN=mean, na.rm=TRUE)
dim(LA_sitemean)
LA_sitemean <-LA_sitemean[,c("Group.1","Group.2","LA_response_ratio")]
names(LA_sitemean) <- c("exp_nam","treatment","LA_response_ratio")
hist(LA_sitemean$LA_response_ratio)

vcmax25_final2_LA <- merge(vcmax25_final2,LA_sitemean,
                              by=c("exp_nam","treatment"),all.x=TRUE)
ggplot(subset(vcmax25_final2_LA,treatment=="f"),
       aes(x=exp_nam, y=response_ratio)) +
  geom_boxplot(aes(x=exp_nam, y=response_ratio),alpha=0.5,color="red") +
  geom_boxplot(aes(x=exp_nam, y=LA_response_ratio),alpha=0.5,color="blue") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of vcmax25 and LA", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()



#how about leaf N
leaf_N <- subset(df,Data_type=="leaf_N")
leaf_N$leafN_response_ratio <- log(leaf_N$elevated/leaf_N$ambient)

leaf_N$treatment_label[leaf_N$treatment=="c"] <- "co2"
leaf_N$treatment_label[leaf_N$treatment=="cf"] <- "co2 + fertilization"
leaf_N$treatment_label[leaf_N$treatment=="cw"] <- "co2+ 1-degree warming"
leaf_N$treatment_label[leaf_N$treatment=="d"] <- "drought"
leaf_N$treatment_label[leaf_N$treatment=="f"] <- "fertilization"
leaf_N$treatment_label[leaf_N$treatment=="w"] <- "1-degree warming"


unit_leafN <- leaf_N%>% group_by(Unit)  %>% summarise(number = n())
leafNmass <-  subset(leaf_N,Unit=="%"|Unit=="g/g"|Unit=="g/kg"|Unit=="mg_/_g"|Unit=="mg/g"|Unit=="g/mg"|
                            Unit=="mg/g_"|Unit=="mg/kg"|Unit=="mol/g"|Unit=="mmol/g")
dim(leafNmass)
leafNarea <-  subset(leaf_N,Unit!="%" & Unit!="g/g" & Unit!="g/kg" & Unit!="mg_/_g" & Unit!="mg/g" & Unit!="g/mg" &
                       Unit!="mg/g_" & Unit!="mg/kg" & Unit!="mol/g" & Unit!="mmol/g")
dim(leafNarea)
#already checked - well divided now

#check nmass
leafNmass_cfw <- subset(leafNmass,treatment=="f"|treatment=="c"|treatment=="w"|
                          treatment=="cf"|treatment=="cw"|treatment=="d")
p1 <- ggplot(leafNmass_cfw, aes(x=treatment_label, y=leafN_response_ratio)) +geom_jitter()+geom_boxplot(alpha=0.5,color="red")+
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of leaf Nmass", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=30))+coord_flip()

#check narea
leafNarea_cfw <- subset(leafNarea,treatment=="f"|treatment=="c"|treatment=="w"|
                          treatment=="cf"|treatment=="cw"|treatment=="d")
p2 <- ggplot(leafNarea_cfw, aes(x=treatment_label, y=leafN_response_ratio)) +geom_jitter()+geom_boxplot(alpha=0.5,color="red")+
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of leaf Narea", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=30))+coord_flip()
plot_grid(p1,p2)
ggsave(paste("~/data/output_gcme/leaf_Narea_Nmass.jpg",sep=""),width = 30, height = 15)

#now, merging vcmax25_final2 (pred+obs) and LMA into df
LMA_SLA <- subset(df,Data_type=="leaf_mass_per_area"|Data_type=="leaf_mass_per_unit_leaf_area"|Data_type=="LMA"|
                Data_type=="specific_leaf_area"|Data_type=="SLA")
LMA_SLA%>% group_by(Unit)  %>% summarise(number = n())# all reasonable unit
LMA_SLA$LMA_response_ratio <- log(LMA_SLA$elevated/LMA_SLA$ambient)
#because they are inverse!
LMA_SLA$LMA_response_ratio[LMA_SLA$Data_type=="specific_leaf_area"] <- log(LMA_SLA$ambient[LMA_SLA$Data_type=="specific_leaf_area"]/LMA_SLA$elevated[LMA_SLA$Data_type=="specific_leaf_area"])
LMA_SLA$LMA_response_ratio[LMA_SLA$Data_type=="SLA"] <- log(LMA_SLA$ambient[LMA_SLA$Data_type=="SLA"]/LMA_SLA$elevated[LMA_SLA$Data_type=="SLA"])
hist(LMA_SLA$LMA_response_ratio)

LMA_SLA_final <-LMA_SLA[,c("exp_nam","treatment","LMA_response_ratio")]
LMA_SLA_final <- aggregate(LMA_SLA_final,by=list(LMA_SLA_final$exp_nam,LMA_SLA_final$treatment),
                                 FUN=mean, na.rm=TRUE)
LMA_SLA_final <-LMA_SLA_final[,c("Group.1","Group.2","LMA_response_ratio")]
names(LMA_SLA_final) <- c("exp_nam","treatment","LMA_response_ratio")
hist(LMA_SLA_final$LMA_response_ratio)
dim(LMA_SLA_final)

dim(leafNmass)
#vcmax25_final# saves all observation but only choose vcmax data
vcmax25_obs <- subset(vcmax25_final,is.na(response_ratio)==FALSE &method=="vcmax")
dim(vcmax25_obs)
vcmax25_obs_final <-vcmax25_obs[,c("exp_nam","treatment","response_ratio")]
vcmax25_obs_final <- aggregate(vcmax25_obs_final,by=list(vcmax25_obs_final$exp_nam,vcmax25_obs_final$treatment),
                           FUN=mean, na.rm=TRUE)
vcmax25_obs_final <-vcmax25_obs_final[,c("Group.1","Group.2","response_ratio")]
names(vcmax25_obs_final) <- c("exp_nam","treatment","vcmax25_response_ratio")

#vcmax25_warmingco2_siteinfo: saves prediction only (c,w) --> filter with method=="vcmax"
vcmax25_pred <- subset(vcmax25_warmingco2_siteinfo,is.na(response_ratio)==FALSE &method=="vcmax")
dim(vcmax25_pred)
vcmax25_pred_final <-vcmax25_pred[,c("exp_nam","treatment","pred_response_ratio")]
vcmax25_pred_final <- aggregate(vcmax25_pred_final,by=list(vcmax25_pred_final$exp_nam,vcmax25_pred_final$treatment),
                               FUN=mean, na.rm=TRUE)
vcmax25_pred_final <-vcmax25_pred_final[,c("Group.1","Group.2","pred_response_ratio")]
names(vcmax25_pred_final) <- c("exp_nam","treatment","pred_vcmax25_response_ratio")

#leafNarea
leafNarea <-leafNarea[,c("exp_nam","treatment","leafN_response_ratio")]
leafNarea_final <- aggregate(leafNarea,by=list(leafNarea$exp_nam,leafNarea$treatment),
                                FUN=mean, na.rm=TRUE)
leafNarea_final <-leafNarea_final[,c("Group.1","Group.2","leafN_response_ratio")]
names(leafNarea_final) <- c("exp_nam","treatment","Narea_response_ratio")

#now, time to merged Nmass (all invidividuals) with narea,obs_vcmax25,pred_vcmax25, lma
dim(leafNmass)
leaf_traits <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp_nam","treatment"),all.x=TRUE),
                 list(leafNmass,leafNarea_final,vcmax25_obs_final,vcmax25_pred_final,LMA_SLA_final))
summary(leaf_traits)

leaf_traits_cfw <- subset(leaf_traits,treatment=="f"|treatment=="c"|treatment=="w"|
                          treatment=="cf"|treatment=="cw"|treatment=="d")
leaf_traits_c <- subset(leaf_traits,treatment=="c")

ggplot(leaf_traits_cfw, aes(x=treatment_label, y=leafN_response_ratio)) +geom_jitter()+geom_boxplot(alpha=0.5,color="red")+
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="Log response ratio of leaf Nmass", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+coord_flip()

ggplot(subset(leaf_traits,treatment=="c"),
       aes(x=exp_nam, y=leafN_response_ratio)) +geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=exp_nam, y=LMA_response_ratio),alpha=0.5,color="green") +
  geom_boxplot(aes(x=exp_nam, y=pred_vcmax25_response_ratio),alpha=0.5,color="orange") +
  geom_boxplot(aes(x=exp_nam, y=vcmax25_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="co2 effect: Log response ratio of Leaf Nmass, measured vcmax25 (red), predicted vcmax25 (orange), and lma (green)", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=30))+coord_flip()
ggsave(paste("~/data/output_gcme/Nmass_vcmax25_LMA_co2.jpg",sep=""),width = 30, height = 15)

ggplot(subset(leaf_traits,treatment=="w"),
       aes(x=exp_nam, y=leafN_response_ratio)) +geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=exp_nam, y=LMA_response_ratio),alpha=0.5,color="green") +
  geom_boxplot(aes(x=exp_nam, y=pred_vcmax25_response_ratio),alpha=0.5,color="orange") +
  geom_boxplot(aes(x=exp_nam, y=vcmax25_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="warming effect: Log response ratio of Leaf Nmass, measured vcmax25 (red), predicted vcmax25 (orange), and lma (green)", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=30))+coord_flip()
ggsave(paste("~/data/output_gcme/Nmass_vcmax25_LMA_warming.jpg",sep=""),width = 30, height = 15)

ggplot(subset(leaf_traits,treatment=="cf"),
       aes(x=exp_nam, y=leafN_response_ratio)) +geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=exp_nam, y=LMA_response_ratio),alpha=0.5,color="green") +
  geom_boxplot(aes(x=exp_nam, y=pred_vcmax25_response_ratio),alpha=0.5,color="orange") +
  geom_boxplot(aes(x=exp_nam, y=vcmax25_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="co2 + fertilization: Log response ratio of Leaf Nmass (black), vcmax25 (red) and lma (green)", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=30))+coord_flip()
ggsave(paste("~/data/output_gcme/Nmass_vcmax25_LMA_co2_fertilization.jpg",sep=""),width = 30, height = 15)

ggplot(subset(leaf_traits,treatment=="f"),
       aes(x=exp_nam, y=leafN_response_ratio)) +geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=exp_nam, y=LMA_response_ratio),alpha=0.5,color="green") +
  geom_boxplot(aes(x=exp_nam, y=pred_vcmax25_response_ratio),alpha=0.5,color="orange") +
  geom_boxplot(aes(x=exp_nam, y=vcmax25_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="fertilization: Log response ratio of Leaf Nmass (black), vcmax25 (red) and lma (green)", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=30))+coord_flip()
ggsave(paste("~/data/output_gcme/Nmass_vcmax25_LMA_fertilization.jpg",sep=""),width = 30, height = 15)

ggplot(subset(leaf_traits,treatment=="f"),
       aes(x=exp_nam, y=leafN_response_ratio)) +geom_boxplot(alpha=0.5,color="black")+
  geom_boxplot(aes(x=exp_nam, y=LMA_response_ratio),alpha=0.5,color="green") +
  geom_boxplot(aes(x=exp_nam, y=Narea_response_ratio),alpha=0.5,color="cyan") +
  geom_boxplot(aes(x=exp_nam, y=vcmax25_response_ratio),alpha=0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5 ) +labs(y="fertilization: Log response ratio of Leaf Nmass (black), Narea(cyan), vcmax25 (red) and lma (green)", x = " ") +ylim(-3,3) +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=30))+coord_flip()
ggsave(paste("~/data/output_gcme/Nmass_vcmax25_LMA_fertilization_v2.jpg",sep=""),width = 30, height = 15)

#now, filter N fertilzation plots and see which variables they have
Nfer <- subset(df,treatment=="f")
dim(Nfer)
Nfer%>% group_by(exp_nam,Data_type)  %>% summarise(number = n())# all reasonable unit
aa <- Nfer%>% group_by(Data_type)  %>% summarise(number = n())# all reasonable unit

vcmax_N_fert_sitename <- subset(vcmax25_final2,treatment=="f")%>% group_by(exp_nam)  %>% summarise(number = n())
Nfer_available <- merge(Nfer,vcmax_N_fert_sitename,by=c("exp_nam"),all.x=TRUE)
Nfer_available <- subset(Nfer_available,number>0)
aa <- Nfer_available%>% group_by(exp_nam,Data_type)  %>% summarise(number = n())# all reasonable unit

dim(vcmax_N_fert_sitename) # 9 plots including N fertilzation effect on vcmax
final_available <- aa%>% group_by(Data_type)  %>% summarise(number = n())

#check some extremely negative points
vcmax_N_fert_sitename
aa <- subset(vcmax25_final2,exp_nam=="Duke_pinus2" & treatment=="f")

#some additionals
dim(Vcmax25_final2_ANPP)
aaa <- Vcmax25_final2_ANPP[grep("f", Vcmax25_final2_ANPP$treatment),]
list_var <- aggregate(aaa,by=list(aaa$exp_nam), FUN=mean, na.rm=TRUE)$Group.1

str( df[,list_var] )
final_Nfer <- (subset(df, exp_nam %in% list_var))
fff <- final_Nfer%>% group_by(Data_type)%>% summarise(number = n())

dim(subset(final_Nfer,Data_type=="specific_leaf_area")%>% group_by(exp_nam)%>% summarise(number = n()))


