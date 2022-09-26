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

# 4 points added - perhaps because there are more than 1 co2_a and co2_e
unique(plot1[,c("site","co2_a","co2_e")])

#they are soyface and riceface_japan_l_3938_14057


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
smith_all_plotmean$site <- smith_all_plotmean$exp_nam

smith_all_plotmean2 <- unique(smith_all_plotmean[,c("exp","site","co2_a","co2_e")])

plot2 <- merge(plot1,smith_all_plotmean2,by=c("exp"),all.x=TRUE)
# 1 points added - perhaps because there are more than 1 co2_a and co2_e 
unique(plot2[,c("site.y","co2_a.y","co2_e.y")])
#it is WheatCE

#fill info
plot2$site.x[is.na(plot2$site.x)==TRUE] <- plot2$site.y[is.na(plot2$site.x)==TRUE]
plot2$co2_a.x[is.na(plot2$co2_a.x)==TRUE] <- plot2$co2_a.y[is.na(plot2$co2_a.x)==TRUE] 
plot2$co2_e.x[is.na(plot2$co2_e.x)==TRUE] <- plot2$co2_e.y[is.na(plot2$co2_e.x)==TRUE]

plot2 <- dplyr::select(plot2, -c("site.y", "co2_a.y","co2_e.y"))

plot2$site.x[is.na(plot2$site.x)==TRUE] <- plot2$exp[is.na(plot2$site.x)==TRUE]

#add fertilization point's co2 a and e
#input co2 change
plot2$co2_a.x[plot2$lon==-84.04] <- 360
plot2$co2_e.x[plot2$lon==-84.04] <- 560 
#https://pubmed.ncbi.nlm.nih.gov/11960756/
plot2$co2_a.x[plot2$lon==20.00] <- 350
plot2$co2_e.x[plot2$lon==20.00] <- 700
#https://www.researchgate.net/publication/242459413_Nutrient_and_genotypic_effects_on_CO2-responsiveness_Photosynthetic_regulation_in_Leucadendron_species_of_a_nutrient-poor_environment
plot2$co2_a.x[plot2$lon==89.14] <-360
plot2$co2_e.x[plot2$lon==89.14] <- 720 

#clarify type-name
plot2 <- subset(plot2,condition!="highN") # remove replicated sites
plot2$condition[plot2$condition=="lowN"] <- "CO2 at high and low N"

plot2$condition[plot2$type_name=="No_fertilization"] <- "CO2 at N and without N"

#remove below since it is repliacted site-name in table (already mentioned in above)
plot2 <- subset(plot2,condition!="Fertilization")


#add temperature info
plot2$temp_a <- NA
plot2$temp_e <- NA

plot2$temp_a[plot2$site=="Black Spruce, ON, Canada"] <- 22
plot2$temp_e[plot2$site=="Black Spruce, ON, Canada"] <- 30

plot2$temp_a[plot2$site=="Corymbia calophylla provs, AU-NSW"] <- 26
plot2$temp_e[plot2$site=="Corymbia calophylla provs, AU-NSW"] <- 32

plot2$temp_a[plot2$site=="Eucalyptus tereticornis provs AU-NSW"] <- 18
plot2$temp_e[plot2$site=="Eucalyptus tereticornis provs AU-NSW"] <- 28.5

plot2$temp_a[plot2$site=="Eucalyptus globulus, AU-NSW"] <- 19.7
plot2$temp_e[plot2$site=="Eucalyptus globulus, AU-NSW"] <- 24.9

plot2$temp_a[plot2$site=="Eucalyptus parramattensis, AU-NSW"] <- 15.1
plot2$temp_e[plot2$site=="Eucalyptus parramattensis, AU-NSW"] <- 18.0

plot2$temp_a[plot2$site=="Eucalyptus tereticornis, AU-NSW"] <- 14.4
plot2$temp_e[plot2$site=="Eucalyptus tereticornis, AU-NSW"] <- 22.6

#all warming, light, N fertilization plots are in forest or shrub. Fully checked.
plot2$ecosystem[is.na(plot2$ecosystem)==T] <- "Forest"

plot2 <- plot2[order(plot2$condition, plot2$site.x),]

csvfile <- paste("~/data/gcme/MS_data/table_org_site.csv")
write.csv(plot2, csvfile, row.names = TRUE)
