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
keyvars <- c("exp_nam", "prev_name","factors", "treatment", "Data_type", "Unit", "Sampling_date", "Start_Year", "Year", "n_plots")
valuevars <-  c("ambient", "ambient_Se", "ambient_Sd", "elevated", "elevated_Se", "elevated_Sd")

#This may be remediated automatically by using only columns that are distinct w.r.t. `valuevars` and `keyvars` (basically all columns in original data, except `ALIAS`).
print(nrow(df))
df <- df %>% distinct_at(vars(one_of(keyvars, valuevars))) %>%mutate( id=1:n() )   # create new ID key (before this was 'ALIAS')
print(nrow(df))

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
csvfile <- paste("~/data/gcme/yunke_coord/coord_df.csv")
write.csv(df_experiments_coord_exp, csvfile, row.names = TRUE)

#now, first attempt about data
#aa <- leafN_all %>% group_by(Unit)  %>% summarise(number = n())

#subset(leafN_all,newvar=="nmass")



#create a function to convert unit to narea (g/m2) and nmass (mg/g)
unit_convert <- function(df,unit_name_start,unit_name_end,magnitude,var_1,var_2,var_3,var_4){
  df[[var_1]][df$Unit==unit_name_start] <- df[[var_1]][df$Unit==unit_name_start]*magnitude
  df[[var_2]][df$Unit==unit_name_start] <- df[[var_2]][df$Unit==unit_name_start]*magnitude
  df[[var_3]][df$Unit==unit_name_start] <- df[[var_3]][df$Unit==unit_name_start]*magnitude
  df[[var_4]][df$Unit==unit_name_start] <- df[[var_4]][df$Unit==unit_name_start]*magnitude
  df$newvar <- NA
  df[["newvar"]][df$Unit==unit_name_start] <- unit_name_end
  return(df)}

subset(df,Data_type=="leaf_N"&Unit=="g/g")

leafN_all <- subset(df,Data_type=="leaf_N")

leafN_all <- unit_convert(leafN_all,"%","nmass",10,"ambient","ambient_Se","elevated","elevated_Se")
leafN_all <- unit_convert(leafN_all,"g/g","nmass",1000,"ambient","ambient_Se","elevated","elevated_Se")
leafN_all <- unit_convert(leafN_all,"g/kg","nmass",1,"ambient","ambient_Se","elevated","elevated_Se")
leafN_all <- unit_convert(leafN_all,"g/mg","nmass",1,"ambient","ambient_Se","elevated","elevated_Se") # it is just mg/g, unit presented wrong

leafN_all <- unit_convert(leafN_all,"g/cm2","narea",1/10000,"ambient","ambient_Se","elevated","elevated_Se")
subset(leafN_all,Unit=="g/mg")
# a test
leafN_all <- subset(df,Data_type=="leaf_N")
Asat_all <- subset(df,Data_type=="Asat")
vcmax_all <- subset(df,Data_type=="maximal_Rubisco_activity_(Vcmax)"|Data_type=="Vcmax"|
                      Data_type=="Vcmax,_maximum_carboxylation_rate_of_Rubisco"|Data_type=="Vcmax25")
Amax_all <- subset(df,Data_type=="Amax")

leafN_all_mean <- aggregate(leafN_all,by=list(leafN_all$exp_nam,leafN_all$factors,leafN_all$treatment), FUN=mean, na.rm=TRUE)[,c(1,2,3,12)]
names(leafN_all_mean) <- c("exp_nam","factors","treatment","value")
Asat_all_mean <- aggregate(Asat_all,by=list(Asat_all$exp_nam,Asat_all$factors,Asat_all$treatment), FUN=mean, na.rm=TRUE)[,c(1,2,3,12)]
names(Asat_all_mean) <- c("exp_nam","factors","treatment","value")
vcmax_all_mean <- aggregate(vcmax_all,by=list(vcmax_all$exp_nam,vcmax_all$factors,vcmax_all$treatment), FUN=mean, na.rm=TRUE)[,c(1,2,3,12)]
names(vcmax_all_mean) <- c("exp_nam","factors","treatment","value")
Amax_all_mean <- aggregate(Amax_all,by=list(Amax_all$exp_nam,Amax_all$factors,Amax_all$treatment), FUN=mean, na.rm=TRUE)[,c(1,2,3,12)]
names(Amax_all_mean) <- c("exp_nam","factors","treatment","value")

leaf_all_mean <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp_nam","factors","treatment"),all.x=TRUE),list(leafN_all_mean,Asat_all_mean,vcmax_all_mean,Amax_all_mean))
names(leaf_all_mean) <- c("exp_nam","factors","treatment","leafn","Asat","vcmax","Amax")



# a test for allocation
leafN_all <- subset(df,Data_type=="leaf_N")
Asat_all <- subset(df,Data_type=="Asat")
vcmax_all <- subset(df,Data_type=="maximal_Rubisco_activity_(Vcmax)"|Data_type=="Vcmax"|
                      Data_type=="Vcmax,_maximum_carboxylation_rate_of_Rubisco"|Data_type=="Vcmax25")
Amax_all <- subset(df,Data_type=="Amax")

leafN_all$change <- (leafN_all$elevated-leafN_all$ambient)/leafN_all$ambient
Asat_all$change <- (Asat_all$elevated-Asat_all$ambient)/Asat_all$ambient
vcmax_all$change <- (vcmax_all$elevated-vcmax_all$ambient)/vcmax_all$ambient
Amax_all$change <- (Amax_all$elevated-Amax_all$ambient)/Amax_all$ambient

leafN_all_mean <- aggregate(leafN_all,by=list(leafN_all$exp_nam,leafN_all$factors,leafN_all$treatment), FUN=mean, na.rm=TRUE)[,c(1,2,3,12,15)]
names(leafN_all_mean) <- c("exp_nam","factors","treatment","value")
Asat_all_mean <- aggregate(Asat_all,by=list(Asat_all$exp_nam,Asat_all$factors,Asat_all$treatment), FUN=mean, na.rm=TRUE)[,c(1,2,3,12,15)]
names(Asat_all_mean) <- c("exp_nam","factors","treatment","value")
vcmax_all_mean <- aggregate(vcmax_all,by=list(vcmax_all$exp_nam,vcmax_all$factors,vcmax_all$treatment), FUN=mean, na.rm=TRUE)[,c(1,2,3,12,15)]
names(vcmax_all_mean) <- c("exp_nam","factors","treatment","value")
Amax_all_mean <- aggregate(Amax_all,by=list(Amax_all$exp_nam,Amax_all$factors,Amax_all$treatment), FUN=mean, na.rm=TRUE)[,c(1,2,3,12,15)]
names(Amax_all_mean) <- c("exp_nam","factors","treatment","value")

leaf_all_mean <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp_nam","factors","treatment"),all.x=TRUE),list(leafN_all_mean,Asat_all_mean,vcmax_all_mean,Amax_all_mean))
names(leaf_all_mean) <- c("exp_nam","factors","treatment","leafn_a","leafn_e","Asat_a","Asat_e","vcmax_a","vcmax_e","Amax_a","Amax_e")

leaf_all_mean_availabile <- subset(leaf_all_mean,leafn_a>0 & Asat_a>0)
dim(leaf_all_mean_availabile)

leaf_all_mean_availabile %>% select(exp_nam) %>% distinct()
leaf_all_mean_availabile %>% group_by(treatment)  %>% summarise(number = n())

dim(subset(df,Data_type=="Asat")%>% select(exp_nam) %>% distinct())
dim(subset(df,Data_type=="GPP")%>% select(exp_nam) %>% distinct())

dim(subset(df,Data_type=="Asat")%>% select(exp_nam,factors,treatment) %>% distinct())
dim(subset(df,Data_type=="leaf_N")%>% select(exp_nam,factors,treatment) %>% distinct())

dim(subset(df,Data_type=="ANPP")%>% select(exp_nam,factors,treatment) %>% distinct())
dim(subset(df,Data_type=="NPP")%>% select(exp_nam,factors,treatment) %>% distinct())
dim(subset(df,Data_type=="ANPP_foliage")%>% select(exp_nam,factors,treatment) %>% distinct())
dim(subset(df,Data_type=="NPP_woody")%>% select(exp_nam,factors,treatment) %>% distinct())
