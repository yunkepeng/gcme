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
df_org <- read_csv("~/data/Kumarathunge_2020_newphy/kumarathunge_2020_newphy/Data/ACi-TGlob_V1.0.csv")
df_vcmax <- read_csv("~/data/Kumarathunge_2020_newphy/kumarathunge_2020_newphy/Data/PPC-TGlob_V1.0.csv")
#now, select year
df_vcmax$year <- str_sub(df_vcmax$Date,-4,-1)
df_vcmax$year <- as.numeric(df_vcmax$year) #some only have month and day but don't know year -> also convert to year


df_vcmax%>% group_by(seed_source_latitude,seed_source_longitude,year)  %>% summarise(number = n())
#around 73 sites

#now, collect Temp_Treatment, CO2_Treatment and Curve_Id in org. data
info_treatment <- df_org%>% group_by(Curve_Id,Temp_Treatment,CO2_Treatment,seed_source_latitude,seed_source_longitude)  %>% summarise(number = n())
info_treatment <- subset(info_treatment,is.na(Temp_Treatment)==FALSE)
info_treatment <- subset(info_treatment, Temp_Treatment!="ambient" & Temp_Treatment!="elevated")

info_treatment2 <- info_treatment[,c("Temp_Treatment","Curve_Id","seed_source_latitude","seed_source_longitude")]
names(info_treatment2) <- c("Temp_Treatment_filled","Curve_Id","seed_source_latitude","seed_source_longitude")
df_vcmax2 <- merge(df_vcmax,info_treatment2,by=c("Curve_Id","seed_source_latitude","seed_source_longitude"),all.x=TRUE)

#after look - it seems to be ok - at least - Temp_Treatment_filled can replace Temp_Treatment now.
#csvfile <- paste("~/data/dushan_test.csv")
#write.csv(df_vcmax2, csvfile, row.names = TRUE)

#this is our most interested plots -  7 plots overall?
subset(df_vcmax2,is.na(Temp_Treatment_filled)==F) %>% group_by(Temp_Treatment_filled,seed_source_latitude,seed_source_longitude)  %>% summarise(number = n())
#all coordinates in NA was from Nick Smith --> never minds..
subset(df_vcmax2,is.na(Temp_Treatment_filled)==F) %>% group_by(seed_source_latitude,seed_source_longitude)  %>% summarise(number = n())

df_vcmax2$Growth_T <- as.numeric(str_sub(df_vcmax2$Temp_Treatment_filled,1,2))

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
df_vcmax2$vcmax25 <- df_vcmax2$Vcmax/calc_vcmax_to_vcmax25(df_vcmax2$Tleaf+273.15,df_vcmax2$Growth_T+273.15)


#show some attempts - for unkonw growthT
final_temp <- subset(df_vcmax2,is.na(Temp_Treatment_filled)==F & is.na(seed_source_latitude)==FALSE & is.na(Growth_T)==FALSE)
#select Tleaf
final_temp2 <- subset(final_temp,Growth_T - Tleaf >= -3 & Growth_T - Tleaf <3)

ggplot(final_temp2, aes(x=Growth_T)) +
  geom_point(aes(y=vcmax25,color=Location))+geom_smooth(aes(y=vcmax25,color=Location))
#Ontario plot reference: growth chamber (same as Nick Smith): https://publish.uwo.ca/~dway4/files/way-and-sage-2008-pce.pdf
#On 7, 8 and 9 July 2014, after roughly 100 days under the treatment conditions: Richmond plot reference: https://academic.oup.com/treephys/article/37/8/1095/3778257
