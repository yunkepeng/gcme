smith_co2 <- read.csv("/Users/yunpeng/data/smith_keenan_gcb/gcb_co2/co2_data.csv")
dim(smith_co2)
smith_co2$exp_nam <- smith_co2$SiteID
nrow(subset(smith_co2,exp_nam=="NevadaFACE"))/2

bb <- (co2_vcmax25 %>% group_by(exp_nam) %>% summarise(number=n()))
aa <- (smith_co2 %>% group_by(exp_nam) %>% summarise(number=n()))
aa$number <- aa$number/2
