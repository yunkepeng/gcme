#test 
df1 <- read_csv("/Users/yunpeng/data/Kumarathunge_2020_newphy/kumarathunge_2020_newphy/Data/ACi-TGlob_V1.0.csv")
df1 <- subset(df1,Growth_condition=="Glasshouse"|Growth_condition=="Whole tree chambers"|Growth_condition=="Growth chamber")
check <-df1 %>% group_by(Dataset,Temp_Treatment)  %>% summarise(number = n())
check

df1a <- subset(df1,Dataset=="Eucalyptus globulus, AU-NSW")
df1a_t <- aggregate(df1a,by=list(df1a$Curve_Id), FUN=mean, na.rm=TRUE)[,c("Group.1","Tair","Tleaf")]
names(df1a_t) <- c("Curve_Id","Tair","Tleaf")

df2 <- read_csv("/Users/yunpeng/data/Kumarathunge_2020_newphy/kumarathunge_2020_newphy/Data/PPC-TGlob_V1.0.csv")
df2 <- subset(df2,Growth_condition=="Glasshouse"|Growth_condition=="Whole tree chambers")
check2 <-df2 %>% group_by(Dataset,Temp_Treatment)  %>% summarise(number = n())
check2

df2a <- subset(df2,Dataset=="Eucalyptus globulus, AU-NSW")

df3 <- merge(df2a,df1a_t,by=c("Curve_Id"),all.x=TRUE)
df3$change <- (df3$Tleaf.x-df3$Tleaf.y)

#Smith: BACE_2016 Gas exchange:
#MS - results: 20.1 ± 0.8 °C and 17.0 ± 0.7 °C in the high warming and no warming plots, respectively. Thus, the average warming was 3.1 °C. 
dd1 <- read.csv("/Users/yunpeng/data/BACE_ACi/bace_2016_gas_exchange/data/gas_exchange.csv")
head(dd1)
dd_p0_t0 <- subset(dd1,precipitation==0 & set_leaf_temperature==25&warming==0)
dd_p0_t3 <- subset(dd1,precipitation==0 & set_leaf_temperature==25&warming==3)

dd_p1_t0 <- subset(dd1,precipitation==1 & set_leaf_temperature==25&warming==0)
dd_p1_t3 <- subset(dd1,precipitation==1 & set_leaf_temperature==25&warming==3)

dd_d1_t0 <- subset(dd1,precipitation==-1 & set_leaf_temperature==25&warming==0)
dd_d1_t3 <- subset(dd1,precipitation==-1 & set_leaf_temperature==25&warming==3)
