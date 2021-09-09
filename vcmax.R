library(readr)
library(dplyr)
library(metafor)  
library(ggplot2)
library(stringr)
df <- read_csv("~/data/gcme/data_received_190325/NewData_wide_CORRECTED2.csv") %>%
  mutate( ambient_Sd  = as.numeric(ambient_Sd),  ambient_Se  = as.numeric(ambient_Se), 
          elevated_Sd = as.numeric(elevated_Sd), elevated_Se = as.numeric(elevated_Se) )

# save experiments names
df_experiments <- df %>% select(exp_nam, prev_name) %>% distinct()

### Duplicate rows
#Some rows are simply duplicates. This is dangerous.
#keyvars <- c("exp_nam","prev_name", "factors", "treatment", "Data_type", "Unit", "Sampling_date", "Start_Year", "Year", "n_plots")
#valuevars <-  c("ambient", "ambient_Se", "ambient_Sd", "elevated", "elevated_Se", "elevated_Sd")

#This may be remediated automatically by using only columns that are distinct w.r.t. `valuevars` and `keyvars` (basically all columns in original data, except `ALIAS`).
#print(nrow(df))
#df <- df %>% distinct_at(vars(one_of(keyvars, valuevars))) %>%mutate( id=1:n() )   # create new ID key (before this was 'ALIAS')
#print(nrow(df))

## Determine all experiments that have more than 1 year data
list_exp_gt1yr <- df %>% filter(!is.na(Year)) %>% group_by(exp_nam) %>% summarise(nyears=max(Year)) %>% filter(nyears>1) %>% select(exp_nam) %>% unlist() %>% unname()

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

check_na <- subset(df_experiments_coord_exp,is.na(lon)==TRUE) %>% select(exp_nam) %>% distinct()
dim(check_na)
check_yes <- subset(df_experiments_coord_exp,is.na(lon)==FALSE) %>% select(exp_nam) %>% distinct()
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

coord_Asat <- leaf_all_coord %>% select(exp_nam,lat,lon,z) %>% distinct()
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

for (i in 1:nrow(NA_years)){
  exp_name <- NA_years[i,1]
  variable_name <- NA_years[i,2]
  if (is.infinite(max(df_leafsites$Sampling_Year[is.na(df_leafsites$Sampling_Year)==FALSE & df_leafsites$exp_nam==exp_name & df_leafsites$Data_type==variable_name]))==FALSE){
    df_leafsites$Sampling_Year[is.na(df_leafsites$Sampling_Year)==TRUE & df_leafsites$exp_nam==exp_name & df_leafsites$Data_type==variable_name] <-max(df_leafsites$Sampling_Year[is.na(df_leafsites$Sampling_Year)==FALSE & df_leafsites$exp_nam==exp_name & df_leafsites$Data_type==variable_name])
  } else {
    df_leafsites$Sampling_Year[is.na(df_leafsites$Sampling_Year)==TRUE & df_leafsites$exp_nam==exp_name & df_leafsites$Data_type==variable_name] <-max(df_leafsites$Sampling_Year[is.na(df_leafsites$Sampling_Year)==FALSE & df_leafsites$exp_nam==exp_name])
  }
}

NA_years <- as.data.frame(subset(df_leafsites,is.na(Sampling_Year)==TRUE)%>% group_by(exp_nam,Data_type)%>% summarise(number = n()))
#only a few missing - no worries!

#now, start working with vcmax data

vcmax_all <- subset(df_leafsites,Data_type=="maximal_Rubisco_activity_(Vcmax)"|Data_type=="Vcmax"|
                      Data_type=="Vcmax,_maximum_carboxylation_rate_of_Rubisco")
#vcmax25_all <- subset(df_leafsites,Data_type=="Vcmax25")

Jmax_all <- subset(df_leafsites,Data_type=="Jmax_"|Data_type=="Jmax"|
                      Data_type=="Jmax,_maximum_rate_of_electron_transport_for_regeneration_of_RuBP")
#Jmax25_all <- subset(df_leafsites,Data_type=="Jmax25")

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

Jmax_all%>% group_by(Unit)  %>% summarise(mean = mean(ambient,na.rm=TRUE))
Jmax_all_final <- Jmax_all

#one point shows start and sampling year was 1998 and 1997 separately, wrong. Convert to same year
Asat_all_final$Sampling_Year[Asat_all_final$Sampling_Year<Asat_all_final$Start_Year] <- Asat_all_final$Start_Year[Asat_all_final$Sampling_Year<Asat_all_final$Start_Year]
Amax_all_final$Sampling_Year[Amax_all_final$Sampling_Year<Amax_all_final$Start_Year] <- Amax_all_final$Start_Year[Amax_all_final$Sampling_Year<Amax_all_final$Start_Year]

Asat_all_final$Year <- Asat_all_final$Sampling_Year-Asat_all_final$Start_Year
Amax_all_final$Year <- Amax_all_final$Sampling_Year-Amax_all_final$Start_Year

vcmax_all_final$Year <- vcmax_all_final$Sampling_Year-vcmax_all_final$Start_Year
Jmax_all_final$Year <- Jmax_all_final$Sampling_Year-Jmax_all_final$Start_Year

#now, all done! start to do one-point method!