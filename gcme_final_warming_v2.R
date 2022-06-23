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
library(Hmisc)
library("missMDA")
library(FactoMineR)
#Smith: BACE_2016 Gas exchange:
#MS - results: 20.1 ± 0.8 °C and 17.0 ± 0.7 °C in the high warming and no warming plots, respectively. Thus, the average warming was 3.1 °C. 
dd1 <- read.csv("/Users/yunpeng/data/BACE_ACi/bace_2016_gas_exchange/data/gas_exchange.csv")
head(dd1)

dd1$warming <- as.character(dd1$warming)
dd1$jv <- dd1$jmax/dd1$vcmax
dd1$lma <- 1/dd1$specific_leaf_area

dd_p0 <- subset(dd1,precipitation==0 & set_leaf_temperature==25)
dd_p1 <- subset(dd1,precipitation==1 & set_leaf_temperature==25)
dd_p2 <- subset(dd1,precipitation==-1 & set_leaf_temperature==25)

dd_p0 %>% ggplot( aes(x=warming, y=vcmax)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="vcmax25",x="precipitation at ambient, with warming effect") + theme_classic()+theme(axis.text=element_text(size=12))

dd_p1 %>% ggplot( aes(x=warming, y=vcmax)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="vcmax25",x="precipitation at ambient * 1.5, with warming effect") + theme_classic()+theme(axis.text=element_text(size=12))

dd_p2 %>% ggplot( aes(x=warming, y=vcmax)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="vcmax25",x="precipitation at ambient * 0.5, with warming effect") + theme_classic()+theme(axis.text=element_text(size=12))

dd_p0 %>% ggplot( aes(x=warming, y=jmax)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="jmax25",x="precipitation at ambient, with warming effect") + theme_classic()+theme(axis.text=element_text(size=12))

dd_p1 %>% ggplot( aes(x=warming, y=jmax)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="jmax25",x="precipitation at ambient * 1.5, with warming effect") + theme_classic()+theme(axis.text=element_text(size=12))

dd_p2 %>% ggplot( aes(x=warming, y=jmax)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="jmax25",x="precipitation at ambient * 0.5, with warming effect") + theme_classic()+theme(axis.text=element_text(size=12))

dd_p0 %>% ggplot( aes(x=warming, y=jv)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="jv",x="precipitation at ambient, with warming effect") + theme_classic()+theme(axis.text=element_text(size=12))

dd_p1 %>% ggplot( aes(x=warming, y=jv)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="jv",x="precipitation at ambient * 1.5, with warming effect") + theme_classic()+theme(axis.text=element_text(size=12))

dd_p2 %>% ggplot( aes(x=warming, y=jv)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="jv",x="precipitation at ambient * 0.5, with warming effect") + theme_classic()+theme(axis.text=element_text(size=12))

dd_p0 %>% ggplot( aes(x=warming, y=lma)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="lma",x="precipitation at ambient, with warming effect") + theme_classic()+theme(axis.text=element_text(size=12))

dd_p1 %>% ggplot( aes(x=warming, y=lma)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="lma",x="precipitation at ambient * 1.5, with warming effect") + theme_classic()+theme(axis.text=element_text(size=12))

dd_p2 %>% ggplot( aes(x=warming, y=lma)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="lma",x="precipitation at ambient * 0.5, with warming effect") + theme_classic()+theme(axis.text=element_text(size=12))

# warming 2 - growth chamber 
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

df$vcmax25_23obs <- df$Vcmax23_obs/calc_vcmax_to_vcmax25(df$Tleaf23_Vcmax_obs+273.15,df$GrowthT+273.15)
#another: just using Arrhenius design
df$vcmax25_23obs_normal <- df$Vcmax23_obs*exp((65330/8.314)*((1/(df$Tleaf23_Vcmax_obs+273.15))-(1/298.15)))

df$jmax25_23obs <- df$Jmax23_obs/calc_jmax_to_jmax25(df$Tleaf23_Jmax_obs+273.15,df$GrowthT+273.15)

df$GrowthT_type <- as.character(df$GrowthT)

df$jmax25_vcmax25 <- df$jmax25_23obs/df$vcmax25_23obs
df$lma <- 10000*df$Lmass_P/df$Larea_P # g/m2
df$nmass <- df$Nper_P #%
df$narea <-  df$nmass*df$lma/100

df %>% ggplot( aes(x=GrowthT_type, y=vcmax25_23obs)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="vcmax25",x="Smith - growth chamber experiment") + theme_classic()+theme(axis.text=element_text(size=12))

df %>% ggplot( aes(x=GrowthT_type, y=jmax25_23obs)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="jmax25",x="Smith - growth chamber experiment") + theme_classic()+theme(axis.text=element_text(size=12))

df %>% ggplot( aes(x=GrowthT_type, y=jmax25_vcmax25)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="jmax25/vcmax25",x="Smith - growth chamber experiment") + theme_classic()+theme(axis.text=element_text(size=12))

df %>% ggplot( aes(x=GrowthT_type, y=lma)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="lma",x="Smith - growth chamber experiment") + theme_classic()+theme(axis.text=element_text(size=12))

df %>% ggplot( aes(x=GrowthT_type, y=nmass)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="nmass",x="Smith - growth chamber experiment") + theme_classic()+theme(axis.text=element_text(size=12))

df %>% ggplot( aes(x=GrowthT_type, y=narea)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +
  labs(y="narea",x="Smith - growth chamber experiment") + theme_classic()+theme(axis.text=element_text(size=12))

#dushan
#test
library(readr)
df1 <- read_csv("/Users/yunpeng/data/Kumarathunge_2020_newphy/kumarathunge_2020_newphy/Data/ACi-TGlob_V1.0.csv")
df1 <- subset(df1,Growth_condition=="Glasshouse"|Growth_condition=="Whole tree chambers"|Growth_condition=="Growth chamber")
check <-df1 %>% group_by(Dataset,Temp_Treatment)  %>% summarise(number = n())
check

df1a_t <- aggregate(df1,by=list(df1$Dataset,df1$Curve_Id), FUN=mean, na.rm=TRUE)[,c("Group.1","Group.2","Tair","Tleaf")]
names(df1a_t) <- c("Dataset","Curve_Id","Tair","Tleaf_agg")

df2 <- read_csv("/Users/yunpeng/data/Kumarathunge_2020_newphy/kumarathunge_2020_newphy/Data/PPC-TGlob_V1.0.csv")
df2 <- subset(df2,Growth_condition=="Glasshouse"|Growth_condition=="Whole tree chambers")
check2 <-df2 %>% group_by(Dataset,Temp_Treatment)  %>% summarise(number = n())
check2

df3 <- merge(df2,df1a_t,by=c("Dataset","Curve_Id"),all.x=TRUE)
plot(df3$Tleaf~df3$Tleaf_agg)

#further merged with temp - produced some plots more!
df_temp <- unique(df1[,c("Dataset","Curve_Id","Temp_Treatment")])
names(df_temp)<- c("Dataset","Curve_Id","Temp_Treatment_agg")
df3 <- merge(df3,df_temp,by=c("Dataset","Curve_Id"),all.x=TRUE)

summary(df3)

#13 sites - here Tair is just 
unique(df3$Dataset)[2]

df3$vcmax25 <- df3$Vcmax*exp((65330/8.314)*((1/(df3$Tleaf+273.15))-(1/298.15)))
df3$jmax25 <- df3$Jmax*exp((43540/8.314)*((1/(df3$Tleaf+273.15))-(1/298.15)))

#subset(df3,Dataset==unique(df3$Dataset)[2]) %>% ggplot( aes(x=Tair, y=vcmax25)) + geom_point(aes(color=Temp_Treatment_agg),alpha = 0.6, width = 0.5) 
#subset(df3,Dataset==unique(df3$Dataset)[2]) %>% ggplot( aes(x=Temp_Treatment_agg, y=vcmax25)) + geom_boxplot()

#subset(df3,Dataset==unique(df3$Dataset)[4]) %>% ggplot( aes(x=Tair, y=vcmax25)) + geom_point(aes(color=Temp_Treatment),alpha = 0.6, width = 0.5) 

#divided into site-species level

#check experiment
check2 <- unique(df3[,c("Dataset","species","Temp_Treatment_agg")])

#For "Eucalyptus globulus, AU-NSW"
#By looking at https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.12314 
# Fig.1 (mean ambient T in whole year = 20) and description: a factorial between ambient temperature (aT) and ambient +3 °C (eT) by ambient CO2 concentration (aC)
#we can assume Tair changed from 20 to 23
df3$Temp_Treatment_agg[df3$Dataset=="Eucalyptus globulus, AU-NSW" & df3$Temp_Treatment_agg=="ambient"] <- "20C"
df3$Temp_Treatment_agg[df3$Dataset=="Eucalyptus globulus, AU-NSW" & df3$Temp_Treatment_agg=="elevated"] <- "23C"

#for "Eucalyptus parramattensis, AU-NSW"
#see below Table s1: we can assume ambient = 15, elevated = 18
#https://nph.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fnph.15668&file=nph15668-sup-0001-SupInfo.pdf
df3$Temp_Treatment_agg[df3$Dataset=="Eucalyptus parramattensis, AU-NSW" & df3$Temp_Treatment_agg=="ambient"] <- "15C"
df3$Temp_Treatment_agg[df3$Dataset=="Eucalyptus parramattensis, AU-NSW" & df3$Temp_Treatment_agg=="elevated"] <- "18C"

#for Eucalyptus tereticornis, AU-NSW
#see below Table s1: we can assume ambient = 14.4, elevated = 22.6???
#it is unpublished data - we may assume that it increases 3 degrees, as consistent with above
#https://nph.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fnph.15668&file=nph15668-sup-0001-SupInfo.pdf
df3$Temp_Treatment_agg[df3$Dataset=="Eucalyptus tereticornis, AU-NSW" & df3$Temp_Treatment_agg=="ambient"] <- "20C"
df3$Temp_Treatment_agg[df3$Dataset=="Eucalyptus tereticornis, AU-NSW" & df3$Temp_Treatment_agg=="elevated"] <- "23C"

w1 <- subset(df3,Dataset==unique(df3$Dataset)[1]) %>% ggplot( aes(x=Temp_Treatment_agg, y=vcmax25)) + geom_boxplot()+
  labs(x=unique(df3$Dataset)[1])

w2 <- subset(df3,Dataset==unique(df3$Dataset)[2]) %>% ggplot( aes(x=Temp_Treatment_agg, y=vcmax25)) + geom_boxplot()+
  labs(x=unique(df3$Dataset)[2])

w3 <- subset(df3,Dataset==unique(df3$Dataset)[4]) %>% ggplot( aes(x=Temp_Treatment_agg, y=vcmax25)) + geom_boxplot()+
  labs(x=unique(df3$Dataset)[4])

w4 <- subset(df3,Dataset==unique(df3$Dataset)[5]) %>% ggplot( aes(x=Temp_Treatment_agg, y=vcmax25)) + geom_boxplot()+
  labs(x=unique(df3$Dataset)[5])

w5 <- subset(df3,Dataset==unique(df3$Dataset)[7]) %>% ggplot( aes(x=Temp_Treatment_agg, y=vcmax25)) + geom_boxplot()+
  labs(x=unique(df3$Dataset)[7])

w6 <- subset(df3,Dataset==unique(df3$Dataset)[8]) %>% ggplot( aes(x=Temp_Treatment_agg, y=vcmax25)) + geom_boxplot()+
  labs(x=unique(df3$Dataset)[8])

#below is nick's data
subset(df3,Dataset==unique(df3$Dataset)[13]) %>% ggplot( aes(x=Temp_Treatment_agg, y=vcmax25)) + geom_boxplot()+
  labs(x=unique(df3$Dataset)[13])

#and gcme
plot_grid(w1,w2,w3,w4,w5,w6,nrow=2,label_size = 15)+theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/warming.jpg",sep=""),width = 15, height = 7.5)

#create a dataframe
d1 <- tibble(
  exp= "Black Spruce, ON, Canada",
  lon = subset(df3,Dataset=="Black Spruce, ON, Canada")$seed_source_longitude[1],
  lat = subset(df3,Dataset=="Black Spruce, ON, Canada")$seed_source_latitude[1],
  vcmax25_a= mean(subset(df3,Dataset=="Black Spruce, ON, Canada" & Temp_Treatment_agg=="22C")$vcmax25,na.rm=TRUE),
  vcmax25_e= mean(subset(df3,Dataset=="Black Spruce, ON, Canada" & Temp_Treatment_agg=="30C")$vcmax25,na.rm=TRUE),
  jmax25_a= mean(subset(df3,Dataset=="Black Spruce, ON, Canada" & Temp_Treatment_agg=="22C")$jmax25,na.rm=TRUE),
  jmax25_e= mean(subset(df3,Dataset=="Black Spruce, ON, Canada" & Temp_Treatment_agg=="30C")$jmax25,na.rm=TRUE),
  t1 = 8)
d2 <- tibble(
  exp= "Corymbia calophylla provs, AU-NSW",
  lon = subset(df3,Dataset=="Corymbia calophylla provs, AU-NSW")$seed_source_longitude[1],
  lat = subset(df3,Dataset=="Corymbia calophylla provs, AU-NSW")$seed_source_latitude[1],
  vcmax25_a= mean(subset(df3,Dataset=="Corymbia calophylla provs, AU-NSW" & Temp_Treatment_agg=="26C")$vcmax25,na.rm=TRUE),
  vcmax25_e= mean(subset(df3,Dataset=="Corymbia calophylla provs, AU-NSW" & Temp_Treatment_agg=="32C")$vcmax25,na.rm=TRUE),
  jmax25_a= mean(subset(df3,Dataset=="Corymbia calophylla provs, AU-NSW" & Temp_Treatment_agg=="26C")$jmax25,na.rm=TRUE),
  jmax25_e= mean(subset(df3,Dataset=="Corymbia calophylla provs, AU-NSW" & Temp_Treatment_agg=="32C")$jmax25,na.rm=TRUE),
  t1 = 6)
d3 <- tibble(
  exp= "Eucalyptus globulus, AU-NSW",
  lon = subset(df3,Dataset=="Eucalyptus globulus, AU-NSW")$seed_source_longitude[1],
  lat = subset(df3,Dataset=="Eucalyptus globulus, AU-NSW")$seed_source_latitude[1],
  vcmax25_a= mean(subset(df3,Dataset=="Eucalyptus globulus, AU-NSW" & Temp_Treatment_agg=="20C")$vcmax25,na.rm=TRUE),
  vcmax25_e= mean(subset(df3,Dataset=="Eucalyptus globulus, AU-NSW" & Temp_Treatment_agg=="23C")$vcmax25,na.rm=TRUE),
  jmax25_a= mean(subset(df3,Dataset=="Eucalyptus globulus, AU-NSW" & Temp_Treatment_agg=="20C")$jmax25,na.rm=TRUE),
  jmax25_e= mean(subset(df3,Dataset=="Eucalyptus globulus, AU-NSW" & Temp_Treatment_agg=="23C")$jmax25,na.rm=TRUE),
  t1 = 3)
d4 <- tibble(
  exp= "Eucalyptus parramattensis, AU-NSW",
  lon = subset(df3,Dataset=="Eucalyptus parramattensis, AU-NSW")$seed_source_longitude[1],
  lat = subset(df3,Dataset=="Eucalyptus parramattensis, AU-NSW")$seed_source_latitude[1],
  vcmax25_a= mean(subset(df3,Dataset=="Eucalyptus parramattensis, AU-NSW" & Temp_Treatment_agg=="15C")$vcmax25,na.rm=TRUE),
  vcmax25_e= mean(subset(df3,Dataset=="Eucalyptus parramattensis, AU-NSW" & Temp_Treatment_agg=="18C")$vcmax25,na.rm=TRUE),
  jmax25_a= mean(subset(df3,Dataset=="Eucalyptus parramattensis, AU-NSW" & Temp_Treatment_agg=="15C")$jmax25,na.rm=TRUE),
  jmax25_e= mean(subset(df3,Dataset=="Eucalyptus parramattensis, AU-NSW" & Temp_Treatment_agg=="18C")$jmax25,na.rm=TRUE),
  t1 = 3)
d5 <- tibble(
  exp= "Eucalyptus tereticornis provs AU-NSW",
  lon = subset(df3,Dataset=="Eucalyptus tereticornis provs AU-NSW")$seed_source_longitude[1],
  lat = subset(df3,Dataset=="Eucalyptus tereticornis provs AU-NSW")$seed_source_latitude[1],
  vcmax25_a= mean(subset(df3,Dataset=="Eucalyptus tereticornis provs AU-NSW" & Temp_Treatment_agg=="18C")$vcmax25,na.rm=TRUE),
  vcmax25_e= mean(subset(df3,Dataset=="Eucalyptus tereticornis provs AU-NSW" & Temp_Treatment_agg=="28.5C")$vcmax25,na.rm=TRUE),
  jmax25_a= mean(subset(df3,Dataset=="Eucalyptus tereticornis provs AU-NSW" & Temp_Treatment_agg=="28.5C")$jmax25,na.rm=TRUE),
  jmax25_e= mean(subset(df3,Dataset=="Eucalyptus tereticornis provs AU-NSW" & Temp_Treatment_agg=="28.5C")$jmax25,na.rm=TRUE),
  t1 = 10.5)
d6 <- tibble(
  exp= "Eucalyptus tereticornis, AU-NSW",
  lon = subset(df3,Dataset=="Eucalyptus tereticornis, AU-NSW")$seed_source_longitude[1],
  lat = subset(df3,Dataset=="Eucalyptus tereticornis, AU-NSW")$seed_source_latitude[1],
  vcmax25_a= mean(subset(df3,Dataset=="Eucalyptus tereticornis, AU-NSW" & Temp_Treatment_agg=="20C")$vcmax25,na.rm=TRUE),
  vcmax25_e= mean(subset(df3,Dataset=="Eucalyptus tereticornis, AU-NSW" & Temp_Treatment_agg=="23C")$vcmax25,na.rm=TRUE),
  jmax25_a= mean(subset(df3,Dataset=="Eucalyptus tereticornis, AU-NSW" & Temp_Treatment_agg=="20C")$jmax25,na.rm=TRUE),
  jmax25_e= mean(subset(df3,Dataset=="Eucalyptus tereticornis, AU-NSW" & Temp_Treatment_agg=="23C")$jmax25,na.rm=TRUE),
  t1 = 3)

#one plot show temp increases vcmax25 - because of Narea?
#https://academic.oup.com/treephys/article/37/8/1095/3778257?login=true#95860081
#We also found strong positive associations between Asat, Vcmax and J (at 25 °C) and Narea in the low temperature treatment, and Narea, Asat, Vcmax and J (at 25 °C) were generally higher in the warm temperature treatment compared with the cool temperature treatment

#now, for GCME (both are soil measurement)
kevin <- read.csv("/Users/yunpeng/data/gcme/kevin/orig_vcmax/JunkePeng_11252021.csv")
kevin_warming_vcmax <- subset(kevin,response=="vcmax" & (exp=="riceface_japan_l_2008_3938_14057_w" |exp=="brandbjerg_w"|exp=="riceface_japan_l_2007_3938_14057_w"))
kevin_warming_vcmax$dominant_species[kevin_warming_vcmax$exp=="riceface_japan_l_2008_3938_14057_w"] <- "rice_piddies"
kevin_warming_jmax <- subset(kevin,response=="jmax" & (exp=="riceface_japan_l_2008_3938_14057_w" |exp=="brandbjerg_w"|exp=="riceface_japan_l_2007_3938_14057_w"))
kevin_warming_jmax$dominant_species[kevin_warming_jmax$exp=="riceface_japan_l_2008_3938_14057_w"] <- "rice_piddies"

#RICE measurements were measured just at 25C:Soil and Water Warming
#Seasonal mean water temperature was 23.2°C (NT) and 25.1°C (ET) in 2007, and 22.8°C (NT) and 25.4°C (ET) in 2008 (Tokida et al. 2010).
#All measurements were conducted between 08:00 and 15:00 h on clear days; block temperature in the cuvette was fixed at 25°C, and the photosynthetic photon flux density was fixed at 1,800 µmol m−2 s−1, with a flow rate of 500 µmol s−1. 
#https://academic.oup.com/pcp/article/55/2/370/1861701#85129192

d7 <- tibble(
  exp= "riceface_japan_l_2008_3938_14057_w",
  lon = subset(kevin_warming_vcmax,exp=="riceface_japan_l_2008_3938_14057_w")$lon[1],
  lat = subset(kevin_warming_vcmax,exp=="riceface_japan_l_2008_3938_14057_w")$lat[1],
  vcmax25_a= mean(subset(kevin_warming_vcmax,exp=="riceface_japan_l_2008_3938_14057_w")$x_c,na.rm=TRUE),
  vcmax25_e= mean(subset(kevin_warming_vcmax,exp=="riceface_japan_l_2008_3938_14057_w")$x_t,na.rm=TRUE),
  jmax25_a= mean(subset(kevin_warming_jmax,exp=="riceface_japan_l_2008_3938_14057_w")$x_c,na.rm=TRUE),
  jmax25_e= mean(subset(kevin_warming_jmax,exp=="riceface_japan_l_2008_3938_14057_w")$x_t,na.rm=TRUE),
  t1 = subset(kevin_warming_vcmax,exp=="riceface_japan_l_2008_3938_14057_w")$w_t3[1])

d8 <- tibble(
  exp= "riceface_japan_l_2007_3938_14057_w",
  lon = subset(kevin_warming_vcmax,exp=="riceface_japan_l_2007_3938_14057_w")$lon[1],
  lat = subset(kevin_warming_vcmax,exp=="riceface_japan_l_2007_3938_14057_w")$lat[1],
  vcmax25_a= mean(subset(kevin_warming_vcmax,exp=="riceface_japan_l_2007_3938_14057_w")$x_c,na.rm=TRUE),
  vcmax25_e= mean(subset(kevin_warming_vcmax,exp=="riceface_japan_l_2007_3938_14057_w")$x_t,na.rm=TRUE),
  jmax25_a= mean(subset(kevin_warming_jmax,exp=="riceface_japan_l_2007_3938_14057_w")$x_c,na.rm=TRUE),
  jmax25_e= mean(subset(kevin_warming_jmax,exp=="riceface_japan_l_2007_3938_14057_w")$x_t,na.rm=TRUE),
  t1 = subset(kevin_warming_vcmax,exp=="riceface_japan_l_2007_3938_14057_w")$w_t2[1])


#brandbjerg cannot find site-species level vcmax25 - not used now
#now, for brandbjerg_w: https://backend.orbit.dtu.dk/ws/portalfiles/portal/74245683/Thesis_FINAL_kboe..PDF
#correct species by looking at Table 3 and 4
#Soil temperature in 5 cm depth was significantly higher in treatments with nighttime warming.
#kevin_warming_vcmax$dominant_species[kevin_warming_vcmax$dominant_species=="temperate_heath"] <- "Deschanpsia flexuosa"
#kevin_warming_vcmax$dominant_species[kevin_warming_vcmax$x_c==87.62] <- "Calluna vulgaris"
#kevin_warming_vcmax$dominant_species[kevin_warming_vcmax$x_c==67.21] <- "Calluna vulgaris"
#kevin_warming_vcmax$dominant_species[kevin_warming_vcmax$x_c==100.44]  <- "Calluna vulgaris"
#kevin_warming_vcmax$dominant_species[kevin_warming_vcmax$x_c==44.19]  <- "Calluna vulgaris"

#while Mean annual temperature is ~10 °C, Nighttime warming results in an increased air temperature of 1.4 °C 
#here we might just assume vcmax = vcmax25 if following paper
#see below:
#To reflect the true ambient temperature of the day, block temperature was set to the expected mean temperature of the measuring day (controlled to be ± 5 °C the target value, due to limitation of the equipment)., 
#and/or : Block temperature was set to 25 °C. Relative humidity was adjusted to 45–60% during measurements.

#we don't have vcmax25 data at species level - but have site level - see their table 1
#vcmax25_warm_am = 112.05 (13.7)
#vcmax25_cold_am = 97.05 (10.06)
#vcmax25_warm_ele = 90.45 (9.05)
#vcmax25_cold_ele = 88.25 (8.05)
#jmax25_warm_am = 141.2 (5.69)
#jmax25_cold_am = 83.98 (4.07)
#jmax25_warm_ele = 141.16 (5.08)
#jmax25_cold_ele = 83.04 (3.98)

d9 <- tibble(
  exp= "brandbjerg_w",
  lon = subset(kevin_warming_vcmax,exp=="brandbjerg_w")$lon[1],
  lat = subset(kevin_warming_vcmax,exp=="brandbjerg_w")$lat[1],
  vcmax25_a= (112.05+97.05)/2,
  vcmax25_e= (90.45+88.25)/2,
  jmax25_a= (141.2+83.98)/2,
  jmax25_e= (141.16+83.04)/2,
  t1 =1.4) # see above

#final dataset - browse
obs_warming <- as.data.frame(rbind(d1,d2,d3,d4,d5,d6,d7,d8,d9))

#siteinfo
obs_warming$year_start <- 1991
obs_warming$year_end <- 2010

forcing_path <- "/Users/yunpeng/data/gcme/kevin/forcing/climate_large/"

#warming_siteinfo <- unique(obs_warming[,c("lon","lat","t1")])
#mylist <- list()
#for (i in 1:length(list.files(forcing_path))) {lonlat<- read.csv(paste(forcing_path,list.files(forcing_path)[i],sep=""))[1,c("lon","lat","sitename")];mylist[[i]] <- lonlat} 
#siteinfo_all <- do.call("rbind",mylist)
#aa <- merge(siteinfo_all,unique(obs_warming[,c("lon","lat","t1")]),by=c("lon","lat"),all.x=TRUE)
#aa <- na.omit(aa)
#unique(aa)

obs_warming$sitename <- c("warming1","warming2","warming7","warming8",
                          "warming10","warming8","kevin37","kevin37","kevin48")

#apply warming
#default inst_vcmax
devtools::load_all("/Users/yunpeng/yunkepeng/gcme/pmodel/rsofun/")
#changing to use adjusted parameters for inst_vcmax from Smith and Keenan 2020 GCB
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
}

obs_warming$logr_obs_vcmax25 <- log(obs_warming$vcmax25_e/obs_warming$vcmax25_a)/obs_warming$t1
obs_warming$logr_obs_jmax25 <- log(obs_warming$jmax25_e/obs_warming$jmax25_a)/obs_warming$t1

obs_warming$logr_pred_vcmax25 <- log(obs_warming$mean_vcmax25_elevated/obs_warming$mean_vcmax25_ambient)/obs_warming$t1
obs_warming$logr_pred_jmax25 <- log(obs_warming$mean_jmax25_elevated/obs_warming$mean_jmax25_ambient)/obs_warming$t1

obs_warming$logr_obs_jv <- obs_warming$logr_obs_jmax25 - obs_warming$logr_obs_vcmax25 
obs_warming$logr_pred_jv <- obs_warming$logr_pred_jmax25 - obs_warming$logr_pred_vcmax25 

obs_warming$treatment <- "prediction"
obs_warming$treatment2 <- "observation"

obs_warming$points <- "air warming"
obs_warming$points[obs_warming$exp=="riceface_japan_l_2008_3938_14057_w"] <- "soil warming"
obs_warming$points[obs_warming$exp=="riceface_japan_l_2007_3938_14057_w"] <- "soil warming"
obs_warming$points[obs_warming$exp=="brandbjerg_w"] <- "soil warming"

final_prediction_v_all <- tibble(type_name="prediction",middle=median(obs_warming$logr_pred_vcmax25),ymin=quantile(obs_warming$logr_pred_vcmax25, 0.25),ymax=quantile(obs_warming$logr_pred_vcmax25, 0.75))
final_observation_v_all <- tibble(type_name="observation",middle=median(obs_warming$logr_obs_vcmax25),ymin=quantile(obs_warming$logr_obs_vcmax25, 0.25),ymax=quantile(obs_warming$logr_obs_vcmax25, 0.75))

w1 <- obs_warming %>%
  ggplot( aes(x=treatment, y=logr_pred_vcmax25)) +
  geom_boxplot(aes(x=treatment, y=logr_pred_vcmax25),width = 0.5,color="red")+
  geom_boxplot(aes(x=treatment2, y=logr_obs_vcmax25),width = 0.5,color="black")+
  geom_point(alpha = 0.6, width = 0.5,size=2,color="red" ) +
  geom_point(aes(x=treatment2, y=logr_obs_vcmax25),alpha = 0.6, width = 0.5,size=2,color="black") +
  geom_crossbar(data=final_observation_v_all,aes(x=type_name,y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="black") +
  geom_crossbar(data=final_prediction_v_all,aes(x=type_name,y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-0.15,0.05)+
  labs(x=" ", y="Vcmax25 response") +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))
w1

final_prediction_j_all <- tibble(type_name="prediction",middle=median(obs_warming$logr_pred_jmax25),ymin=quantile(obs_warming$logr_pred_jmax25, 0.25),ymax=quantile(obs_warming$logr_pred_jmax25, 0.75))
final_observation_j_all <- tibble(type_name="observation",middle=median(obs_warming$logr_obs_jmax25),ymin=quantile(obs_warming$logr_obs_jmax25, 0.25),ymax=quantile(obs_warming$logr_obs_jmax25, 0.75))

w2 <- obs_warming %>%
  ggplot( aes(x=treatment, y=logr_pred_jmax25)) +
  geom_boxplot(aes(x=treatment, y=logr_pred_jmax25),width = 0.5,color="red")+
  geom_boxplot(aes(x=treatment2, y=logr_obs_jmax25),width = 0.5,color="black")+
  geom_point(alpha = 0.6, width = 0.5,size=2,color="red" ) +
  geom_point(aes(x=treatment2, y=logr_obs_jmax25),alpha = 0.6, width = 0.5,size=2,color="black") +
  geom_crossbar(data=final_observation_j_all,aes(x=type_name,y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="black") +
  geom_crossbar(data=final_prediction_j_all,aes(x=type_name,y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-0.15,0.05)+
  labs(x=" ", y="Jmax25 response") +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))
w2

final_prediction_jv_all <- tibble(type_name="prediction",middle=median(obs_warming$logr_pred_jv),ymin=quantile(obs_warming$logr_pred_jv, 0.25),ymax=quantile(obs_warming$logr_pred_jv, 0.75))
final_observation_jv_all <- tibble(type_name="observation",middle=median(obs_warming$logr_obs_jv),ymin=quantile(obs_warming$logr_obs_jv, 0.25),ymax=quantile(obs_warming$logr_obs_jv, 0.75))
w3 <- obs_warming %>%
  ggplot( aes(x=treatment, y=logr_pred_jv)) +
  geom_boxplot(aes(x=treatment, y=logr_pred_jv),width = 0.5,color="red")+
  geom_boxplot(aes(x=treatment2, y=logr_obs_jv),width = 0.5,color="black")+
  geom_point(alpha = 0.6, width = 0.5,size=2,color="red" ) +
  geom_point(aes(x=treatment2, y=logr_obs_jv),alpha = 0.6, width = 0.5,size=2,color="black") +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-0.15,0.05)+
  labs(x=" ", y="Jmax25/Vcmax25 response") +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))
w3


plot_grid(w1,w2,w3,nrow=1,label_size = 15)+
  theme(plot.background=element_rect(fill="white", color="white"))

ggsave(paste("~/data/output_gcme/colin/warming1.jpg",sep=""),width = 15, height = 5)