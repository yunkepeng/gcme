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
aaa <- df[grep("PP", df$Data_type),] # GPP,NPP
aaa <- df[grep("production", df$Data_type),] # production
aaa <- df[grep("/y", df$Data_type),] # production

aa <- aaa%>% group_by(Data_type)%>% summarise(number = n())