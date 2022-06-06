

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


subset(df3,Dataset==unique(df3$Dataset)[1]) %>% ggplot( aes(x=Temp_Treatment_agg, y=vcmax25)) + geom_boxplot()+
  labs(x=unique(df3$Dataset)[1])

subset(df3,Dataset==unique(df3$Dataset)[2]) %>% ggplot( aes(x=Temp_Treatment_agg, y=vcmax25)) + geom_boxplot()+
  labs(x=unique(df3$Dataset)[2])

subset(df3,Dataset==unique(df3$Dataset)[4]) %>% ggplot( aes(x=Temp_Treatment_agg, y=vcmax25)) + geom_boxplot()+
  labs(x=unique(df3$Dataset)[4])

subset(df3,Dataset==unique(df3$Dataset)[5]) %>% ggplot( aes(x=Temp_Treatment_agg, y=vcmax25)) + geom_boxplot()+
  labs(x=unique(df3$Dataset)[5])

subset(df3,Dataset==unique(df3$Dataset)[7]) %>% ggplot( aes(x=Temp_Treatment_agg, y=vcmax25)) + geom_boxplot()+
  labs(x=unique(df3$Dataset)[7])

subset(df3,Dataset==unique(df3$Dataset)[8]) %>% ggplot( aes(x=Temp_Treatment_agg, y=vcmax25)) + geom_boxplot()+
  labs(x=unique(df3$Dataset)[8])

#below is nick's data
subset(df3,Dataset==unique(df3$Dataset)[13]) %>% ggplot( aes(x=Temp_Treatment_agg, y=vcmax25)) + geom_boxplot()+
  labs(x=unique(df3$Dataset)[13])

#and gcme
