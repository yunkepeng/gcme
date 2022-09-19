library(readr) 
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(Hmisc)
library(missMDA)
library(FactoMineR)
library(plotrix)


rm(list=ls())
#work on site description table
plot_output <- read.csv("~/data/gcme/MS_data/plot_data.csv")
plot_output2 <- plot_output[,c("exp","condition","lon","lat","vcmax","ecosystem","type_name","ref","comments","ecm_type")]
  
#merge with co2_a, co2_e, site-name in GCME
kevin_othervars <- read.csv("~/data/gcme/kevin_20220222/MESI_2022.csv")
kevin_othervars <- rename(kevin_othervars, c(ambient = x_c, elevated=x_t, ambient_Sd=sd_c, elevated_Sd=sd_t,ambient_Se=se_c,elevated_Se=se_t,n_plots=rep_c,
                                             z=elevation, co2_a=c_c, co2_e=c_t, nfertQ_a = n_c, nfertQ_e = n_t, pfertQ_a = p_c, pfertQ_e = p_t,kfertQ_a = k_c, kfertQ_e = k_t,
                                             warmQ_e1 = w_t1, warmQ_e2 = w_t2, warmQ_e3 = w_t3, Unit=x_units))
kevin_othervars <- subset(kevin_othervars,response=="vcmax")
siteinfo <- unique(kevin_othervars[,c("site","exp","co2_a","co2_e")])
plot1 <- merge(plot_output2,siteinfo,by=c("exp"),all.x=TRUE)

#merge with co2_a, co2_e, site-name in Smith
smith_co2 <- read.csv("~/data/smith_keenan_gcb/gcb_co2/co2_data.csv")
smith_co2$exp_nam <- smith_co2$SiteID
smith_co2 <- subset(smith_co2,is.na(lat)==FALSE) # remove empty column
summary(smith_co2)
smith_final <- smith_co2%>% group_by(exp_nam,Species)  %>% summarise(number = n())

#first 2 columns are the same treatment in same plot (amb + elv) -check their Tleaf --> all equal
for (i in 1:((nrow(smith_co2)/2))){
  smith_co2$logr_v[2*i-1] <- NA
  smith_co2$sen_coef_v[2*i-1] <- NA
  smith_co2$logr_j[2*i-1] <- NA
  smith_co2$sen_coef_j[2*i-1] <- NA
  #smith_co2$logr_jv[2*i-1] <- NA
  #smith_co2$sen_coef_jv[2*i-1] <- NA
  
  #fill 1st line data
  smith_co2$ambient_co2[2*i-1] <- NA
  smith_co2$ambient_vcmax[2*i-1] <- NA
  smith_co2$ambient_vcmaxSE[2*i-1] <- NA
  smith_co2$ambient_jmax[2*i-1] <- NA
  smith_co2$ambient_jmaxSE[2*i-1] <-NA
  smith_co2$ambient_vcmax_tmean[2*i-1] <- NA
  smith_co2$ambient_jmax_tmean[2*i-1] <- NA
  #smith_co2$ambient_jv_tmean[2*i-1] <- NA
  
  smith_co2$logr_v[2*i] <- log(smith_co2$Vcmax[2*i]/smith_co2$Vcmax[2*i-1])
  smith_co2$sen_coef_v[2*i] <- log(smith_co2$Vcmax[2*i]/smith_co2$Vcmax[2*i-1])/log(smith_co2$Ca[2*i]/smith_co2$Ca[2*i-1])
  smith_co2$logr_j[2*i] <- log(smith_co2$Jmax[2*i]/smith_co2$Jmax[2*i-1])
  smith_co2$sen_coef_j[2*i] <- log(smith_co2$Jmax[2*i]/smith_co2$Jmax[2*i-1])/log(smith_co2$Ca[2*i]/smith_co2$Ca[2*i-1])
  #smith_co2$logr_jv[2*i] <- log((smith_co2$Jmax[2*i]/smith_co2$Vcmax[2*i])/(smith_co2$Jmax[2*i-1]/smith_co2$Vcmax[2*i-1]))
  #smith_co2$sen_coef_jv[2*i] <- log((smith_co2$Jmax[2*i]/smith_co2$Vcmax[2*i])/(smith_co2$Jmax[2*i-1]/smith_co2$Vcmax[2*i-1]))/log(smith_co2$Ca[2*i]/smith_co2$Ca[2*i-1])
  
  #fill 1st line data
  smith_co2$ambient_co2[2*i] <- smith_co2$Ca[2*i-1]
  smith_co2$ambient_vcmax[2*i] <- smith_co2$Vcmax[2*i-1]
  smith_co2$ambient_vcmaxSE[2*i] <- smith_co2$SE[2*i-1]
  smith_co2$ambient_jmax[2*i] <- smith_co2$Jmax[2*i-1]
  smith_co2$ambient_jmaxSE[2*i] <-smith_co2$Jmax_SE[2*i-1]
  smith_co2$ambient_vcmax_tmean[2*i] <- smith_co2$Vcmax_tmean[2*i-1]
  smith_co2$ambient_jmax_tmean[2*i] <- smith_co2$Jmax_tmean[2*i-1]
  #smith_co2$ambient_jv_tmean[2*i] <- smith_co2$JV_tmean[2*i-1]
}

smith_co2_final <- subset(smith_co2,is.na(logr_v)==FALSE)

vcmax25_warmingco2_siteinfo <- rename(smith_co2_final, c(ambient = ambient_vcmax, elevated=Vcmax, ambient_Se=ambient_vcmaxSE, elevated_Se=SE,
                                                         co2_a=ambient_co2, co2_e=Ca,
                                                         ambient_j = ambient_jmax, elevated_j =Jmax, ambient_Se_j=ambient_jmaxSE, elevated_Se_j=Jmax_SE))


vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$Stature=="Tree"] <- "forest"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="ChinaRiceFACE"] <- "cropland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="UAMACFACE"] <- "cropland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="GlycineCE"] <- "cropland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="GossypiumGC"] <- "cropland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="MississippiCE"] <- "cropland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="WheatCE"] <- "cropland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="ChinaMiniFACE"] <- "cropland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="GrassOTC"] <- "grassland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="PHACE"] <- "grassland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="FescueCE"] <- "grassland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="NZFACE"] <- "grassland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="NevadaFACE"] <- "grassland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="BioCON"] <- "grassland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="PhalarisGC"] <- "grassland"
vcmax25_warmingco2_siteinfo$ecosystem[vcmax25_warmingco2_siteinfo$SiteID=="SwissFACE"] <- "grassland"

#now, combine with Kevin:
smith_all_plotmean <- subset(vcmax25_warmingco2_siteinfo,exp_nam!="NZFACE" &exp_nam!="BioCON" & exp_nam!="ChinaRiceFACE" & exp_nam!="DukeFACE" & exp_nam!="EUROPOPFACE" & exp_nam!="NevadaFACE" & exp_nam!="SwissFACE")

smith_all_plotmean$exp[smith_all_plotmean$exp=="OakOTC"] <- "mi_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="ORNL"] <- "ornerp_liqui_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="PineOTC"] <- "christchurch_pr_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="RichmondGH"] <- "australiacotton_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="UIAPine"] <- "ua_otc_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="PHACE"] <- "phace_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="AspenFACE"] <- "rhine-aspenface_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="BilyKriz"] <- "bily_kriz_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="Headley"] <- "headley_qp_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="Viesalm"] <- "vielsalm_c"

smith_all_plotmean$exp <- tolower(smith_all_plotmean$exp)
smith_all_plotmean2 <- unique(smith_all_plotmean[,c("exp","exp_nam","co2_a","co2_e")])

plot2 <- merge(plot1,smith_all_plotmean2,by=c("exp"),all.x=TRUE)
plot2$site[is.na(plot2$site)==TRUE] <- plot2$exp_nam[is.na(plot2$site)==TRUE]

plot2$co2_a.x[is.na(plot2$co2_a.x)==TRUE] <- plot2$co2_a.y[is.na(plot2$co2_a.x)==TRUE] 
plot2$co2_e.x[is.na(plot2$co2_e.x)==TRUE] <- plot2$co2_e.y[is.na(plot2$co2_e.x)==TRUE]

csvfile <- paste("~/data/gcme/MS_data/table_org_site.csv")
write.csv(plot2, csvfile, row.names = TRUE)
