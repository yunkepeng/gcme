#Nick's data

smith_co2 <- read.csv("/Users/yunpeng/data/smith_keenan_gcb/gcb_co2/co2_data.csv")
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

smith_co2_final <- rename(smith_co2_final, c(ambient = ambient_vcmax, elevated=Vcmax, ambient_Se=ambient_vcmaxSE, elevated_Se=SE,
                        co2_a=ambient_co2, co2_e=Ca,
                        ambient_j = ambient_jmax, elevated_j =Jmax, ambient_Se_j=ambient_jmaxSE, elevated_Se_j=Jmax_SE))

#start pmodel
#combine with siteinfo
#see forcing_kevin_nick.R
siteinfo <- read.csv("/Users/yunpeng/data/smith_keenan_gcb/gcb_co2/forcing/forcing_info.csv")
siteinfo <- siteinfo[,c("lon","lat","elv","year_start","year_end","sitename")]
names(siteinfo) <- c("lon","lat","z","year_start","year_end","sitename")
vcmax25_warmingco2_siteinfo <- merge(smith_co2_final,siteinfo,by=c("lon","lat","z"),all.x=TRUE)

devtools::load_all("/Users/yunpeng/yunkepeng/gcme/pmodel/rsofun/")

df_soiltexture <- bind_rows(
  top    = tibble(layer = "top",    fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1),
  bottom = tibble(layer = "bottom", fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1))
params_modl <- list(
  kphio           = 0.09423773,
  soilm_par_a     = 0.33349283,
  soilm_par_b     = 1.45602286,
  tau_acclim_tempstress = 10,
  par_shape_tempstress  = 0.0)

forcing_path <- "/Users/yunpeng/data/smith_keenan_gcb/gcb_co2/forcing/climate/"
vcmax25_warmingco2_siteinfo$elv <- vcmax25_warmingco2_siteinfo$z

for (i in 1:nrow(vcmax25_warmingco2_siteinfo)){
  siteinfo_site <- vcmax25_warmingco2_siteinfo[i,c("sitename","lon","lat","elv","year_start","year_end")]
  sitename <- siteinfo_site$sitename
  
  df1 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for ambient
  df2 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for elevated
  
  #Q2 - if having co2 then replacing it directly - if not then using co2 at measurement year
  if (is.na(vcmax25_warmingco2_siteinfo$co2_a[i])==TRUE|is.na(vcmax25_warmingco2_siteinfo$co2_e[i])==TRUE) { 
    df1$co2 <- df1$co2
    df2$co2 <- df2$co2
  } else {
    df1$co2 <- vcmax25_warmingco2_siteinfo$co2_a[i] # set to ambient 
    df2$co2 <- vcmax25_warmingco2_siteinfo$co2_e[i] # set to elevated 
  }
  
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
  
  vcmax25_warmingco2_siteinfo$mean_vcmax25_ambient[i] <- mean_vcmax25_ambient
  vcmax25_warmingco2_siteinfo$mean_vcmax25_elevated[i] <- mean_vcmax25_elevated
  vcmax25_warmingco2_siteinfo$mean_vcmax_ambient[i] <- mean_vcmax_ambient
  vcmax25_warmingco2_siteinfo$mean_vcmax_elevated[i] <- mean_vcmax_elevated
  vcmax25_warmingco2_siteinfo$pred_vcmax25_logr[i] <- log(vcmax25_warmingco2_siteinfo$mean_vcmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_vcmax25_ambient[i])
  vcmax25_warmingco2_siteinfo$pred_vcmax25_coef[i] <- (log(vcmax25_warmingco2_siteinfo$mean_vcmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_vcmax25_ambient[i]))/(log(vcmax25_warmingco2_siteinfo$co2_e[i]/vcmax25_warmingco2_siteinfo$co2_a[i]))
  
  vcmax25_warmingco2_siteinfo$mean_jmax25_ambient[i] <- mean_jmax25_ambient
  vcmax25_warmingco2_siteinfo$mean_jmax25_elevated[i] <- mean_jmax25_elevated
  vcmax25_warmingco2_siteinfo$mean_jmax_ambient[i] <- mean_jmax_ambient
  vcmax25_warmingco2_siteinfo$mean_jmax_elevated[i] <- mean_jmax_elevated
  vcmax25_warmingco2_siteinfo$pred_jmax25_logr[i] <- log(vcmax25_warmingco2_siteinfo$mean_jmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_jmax25_ambient[i])
  vcmax25_warmingco2_siteinfo$pred_jmax25_coef[i] <- (log(vcmax25_warmingco2_siteinfo$mean_jmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_jmax25_ambient[i]))/(log(vcmax25_warmingco2_siteinfo$co2_e[i]/vcmax25_warmingco2_siteinfo$co2_a[i]))
  
  vcmax25_warmingco2_siteinfo$PPFD[i] <- mean(df1$ppfd,na.rm=TRUE)*1000000 # in umol/m2/s
  vcmax25_warmingco2_siteinfo$vpd[i] <- mean(df1$vpd,na.rm=TRUE)/1000 #in kPa
  vcmax25_warmingco2_siteinfo$temp[i] <- mean(df1$temp,na.rm=TRUE) #in degree celcius
}

vcmax25_warmingco2_siteinfo$pft <- NA
vcmax25_warmingco2_siteinfo$pft[vcmax25_warmingco2_siteinfo$Stature=="Tree"] <- "forest"
vcmax25_warmingco2_siteinfo$pft[vcmax25_warmingco2_siteinfo$SiteID=="ChinaRiceFACE"] <- "cropland"
vcmax25_warmingco2_siteinfo$pft[vcmax25_warmingco2_siteinfo$SiteID=="UAMACFACE"] <- "cropland"
vcmax25_warmingco2_siteinfo$pft[vcmax25_warmingco2_siteinfo$SiteID=="GlycineCE"] <- "cropland"
vcmax25_warmingco2_siteinfo$pft[vcmax25_warmingco2_siteinfo$SiteID=="GossypiumGC"] <- "cropland"
vcmax25_warmingco2_siteinfo$pft[vcmax25_warmingco2_siteinfo$SiteID=="MississippiCE"] <- "cropland"
vcmax25_warmingco2_siteinfo$pft[vcmax25_warmingco2_siteinfo$SiteID=="WheatCE"] <- "cropland"
vcmax25_warmingco2_siteinfo$pft[vcmax25_warmingco2_siteinfo$SiteID=="ChinaMiniFACE"] <- "cropland"
vcmax25_warmingco2_siteinfo$pft[vcmax25_warmingco2_siteinfo$SiteID=="GrassOTC"] <- "grassland"
vcmax25_warmingco2_siteinfo$pft[vcmax25_warmingco2_siteinfo$SiteID=="PHACE"] <- "grassland"
vcmax25_warmingco2_siteinfo$pft[vcmax25_warmingco2_siteinfo$SiteID=="FescueCE"] <- "grassland"
vcmax25_warmingco2_siteinfo$pft[vcmax25_warmingco2_siteinfo$SiteID=="NZFACE"] <- "grassland"
vcmax25_warmingco2_siteinfo$pft[vcmax25_warmingco2_siteinfo$SiteID=="NevadaFACE"] <- "grassland"
vcmax25_warmingco2_siteinfo$pft[vcmax25_warmingco2_siteinfo$SiteID=="BioCON"] <- "grassland"
vcmax25_warmingco2_siteinfo$pft[vcmax25_warmingco2_siteinfo$SiteID=="PhalarisGC"] <- "grassland"
vcmax25_warmingco2_siteinfo$pft[vcmax25_warmingco2_siteinfo$SiteID=="SwissFACE"] <- "grassland"

csvfile <- paste("/Users/yunpeng/data/smith_keenan_gcb/gcb_co2/pred_vcmax.csv")
write.csv(vcmax25_warmingco2_siteinfo, csvfile, row.names = TRUE)



smith_co2_final_removal <- subset(vcmax25_warmingco2_siteinfo,exp_nam!="BioCON" & exp_nam!="ChinaRiceFACE" & exp_nam!="DukeFACE" & exp_nam!="EUROPOPFACE" & exp_nam!="NevadaFACE" & exp_nam!="SwissFACE")

#see number of samples within each site
ggplot(smith_co2_final_removal,aes(x=exp_nam, y=ambient)) +
  geom_point( size=2)+
  geom_errorbar(aes(ymin=ambient-ambient_Se, ymax=ambient+ambient_Se), width=.2,color = "red",
                position=position_dodge(0.05))+
  labs(x = "Ambient condition: vcmax with SE at each plot",y=" ")+theme_classic()


#now, aggregate
pred_vc_plotmean_smith <- aggregate(smith_co2_final_removal,by=list(smith_co2_final_removal$exp_nam,smith_co2_final_removal$Nac), FUN=mean, na.rm=TRUE)[,c("Group.1","Group.2","pred_vcmax25_coef","sen_coef_v")]
names(pred_vc_plotmean_smith) <- c("exp_nam","ecm_type","pred_vcmax25_coef","middle")

pred_j_plotmean_smith <- aggregate(smith_co2_final_removal,by=list(smith_co2_final_removal$exp_nam,smith_co2_final_removal$Nac), FUN=mean, na.rm=TRUE)[,c("Group.1","Group.2","pred_jmax25_coef","sen_coef_j")]
names(pred_j_plotmean_smith) <- c("exp_nam","ecm_type","pred_jmax25_coef","middle")

#combine with gcme data
vcmax25_plotmean$paper <- "GCME"
pred_vc_plotmean_smith$paper <- "Smith"

vc_twopapers <- dplyr::bind_rows(vcmax25_plotmean, pred_vc_plotmean_smith) 

ggplot(subset(vc_twopapers,paper=="GCME"),aes(x=middle, y=pred_vcmax25_coef)) +
  geom_point(aes(color=ecm_type,shape=paper),size=2)+
  geom_smooth(method="lm")+
  labs(x = "measured sensitivity coef",y="predicted sensitivity coef")+theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

smith_co2_final_removal %>%
  ggplot( aes(x=exp_nam, y=sen_coef_v,color=Nac)) +
  geom_jitter(alpha=1) +
  geom_boxplot(aes(x=exp_nam, y=pred_vcmax25_coef),color="purple",size=1)+ 
  geom_hline( yintercept=0.0, size=0.5 )+
  geom_boxplot(size=1)+ylim(-2,2)+
  labs(x="", y="sensitivity coefficient of vcmax25 - smith's data ") +
  coord_flip() 

smith_co2_final_removal %>%
  ggplot( aes(x=Nac, y=sen_coef_v)) +
  geom_jitter(alpha=1) +
  geom_hline( yintercept=0.0, size=0.5 )+
  geom_boxplot()+geom_boxplot(aes(y=pred_vcmax25_coef),color="red")+ylim(-2,2)+
  labs(x="", y="sensitivity coefficient of vcmax25 - smith's data ") +
  coord_flip() 

smith_co2_final_removal %>%
  ggplot( aes(x=exp_nam, y=sen_coef_j)) +
  geom_jitter(alpha=1) +
  geom_boxplot(aes(x=exp_nam, y=pred_jmax25_coef),color="purple",size=1)+ 
  geom_hline( yintercept=0.0, size=0.5 )+
  geom_boxplot(aes(color=Nac),size=1)+ylim(-2,2)+
  labs(x="", y="sensitivity coefficient of jmax25 - smith's data ") +
  coord_flip() 

smith_co2_final_removal %>%
  ggplot( aes(x=Nac, y=sen_coef_j)) +
  geom_jitter(alpha=1) +
  geom_hline( yintercept=0.0, size=0.5 )+
  geom_boxplot()+geom_boxplot(aes(y=pred_jmax25_coef),color="red")+geom_boxplot(aes(y=pred_jmax25_coef),color="red")+ylim(-2,2)+
  labs(x="", y="sensitivity coefficient of jmax25 - smith's data ") +
  coord_flip() 
