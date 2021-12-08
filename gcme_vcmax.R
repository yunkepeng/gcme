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
library(MAd)


#read Kevin
kevin <- read.csv("/Users/yunpeng/data/gcme/kevin/orig_vcmax/JunkePeng_11252021.csv")
#new <- subset(kevin)%>% group_by(site,lat,lon)  %>% summarise(number = n());new
#repeated sitename: itatinga,riceface_japan_ko_3558_13960,riceface_japan_ta_3558_13960
#kevin_final <- subset(kevin)%>% group_by(site,dominant_species)  %>% summarise(number = n())

#1. correct exp_nam to make it consistent with GCME
kevin$site[is.na(kevin$site)==TRUE] <- "riceface_japan_a_2003_3938_14057"
kevin$exp_nam <- kevin$site #new plot: capo_caccia,itatinga,itatinga_experimental_station,setres
#unify sitename from japan 
kevin[grep("riceface_japan", kevin$exp_nam),]$exp_nam <- "riceface_japan"

#kevin$exp_nam[kevin$site=="biocon"] <- "BioCON";kevin$exp_nam[kevin$site=="brandbjerg"] <- "Brandbjerg";kevin$exp_nam[kevin$site=="dukeface"] <- "DUKE";kevin$exp_nam[kevin$site=="eucface"] <- "EucFACE";kevin$exp_nam[kevin$site=="euroface"] <- "EUROFACE4_pe";kevin$exp_nam[kevin$site=="facts_ii_face"] <- "FACTS_II_FACE3_pt";
#kevin$exp_nam[kevin$site=="giface"] <- "GiFACE";kevin$exp_nam[kevin$site=="eucface"] <- "EucFACE";kevin$exp_nam[kevin$site=="nevada_desert_face"] <- "Nevada_Desert_FACE";kevin$exp_nam[kevin$site=="new_zealand_face"] <- "New_Zealand_FACE";kevin$exp_nam[kevin$site=="popface"] <- "POPFACE_pe";
#kevin$exp_nam[kevin$site=="riceface_china_33n_120e"] <- "RiceFACE_China_33N_120E_Or_2";kevin$exp_nam[kevin$site=="riceface_japan_ko_3558_13960"] <- "RiceFACE_Japan_Ko_2012_35,58_139,60";kevin$exp_nam[kevin$site=="riceface_japan_ko_3558_13961"] <- "RiceFACE_Japan_Ko_2012_35,58_139,60";kevin$exp_nam[kevin$site=="riceface_japan_ko_3558_13962"] <- "RiceFACE_Japan_Ko_2013_35,58_139,60";kevin$exp_nam[kevin$site=="riceface_japan_ko_3558_13963"] <- "RiceFACE_Japan_Ko_2013_35,58_139,60";
#kevin$exp_nam[kevin$site=="riceface_japan_l_3938_14057"] <- "RiceFACE_Japan_L_2007_39,38_140,57";kevin$exp_nam[kevin$site=="riceface_japan_ta_3558_13960"] <- "RiceFACE_Japan_Ta_2012_35,58_139,60";kevin$exp_nam[kevin$site=="riceface_japan_ta_3558_13961"] <- "RiceFACE_Japan_Ta_2012_35,58_139,60";kevin$exp_nam[kevin$site=="riceface_japan_ta_3558_13962"] <- "RiceFACE_Japan_Ta_2013_35,58_139,60";kevin$exp_nam[kevin$site=="riceface_japan_ta_3558_13963"] <- "RiceFACE_Japan_Ta_2013_35,58_139,60";
#kevin$exp_nam[kevin$site=="soyface"] <- "SoyFACEsoy1";kevin$exp_nam[kevin$site=="soyface"] <- "soyfaceb";kevin$exp_nam[kevin$site=="swissface_trifolium"] <- "SwissFACE_trifolium2"

#change a list of name

kevin <- rename(kevin, c(ambient = x_c, elevated=x_t, ambient_Sd=sd_c, elevated_Sd=sd_t,ambient_Se=se_c,elevated_Se=se_t,n_plots=rep_c,
                          z=elevation, co2_a=c_c, co2_e=c_t, nfertQ_a = n_c, nfertQ_e = n_t, pfertQ_a = p_c, pfertQ_e = p_t,kfertQ_a = k_c, kfertQ_e = k_t,
                          warmQ_e1 = w_t1, warmQ_e2 = w_t2, warmQ_e3 = w_t3, Unit=x_units))

#adjust temperature response
kevin$warmQ_e2[is.na(kevin$warmQ_e2)==TRUE] <- 0
kevin$warmQ_e2 # either 0, 1, or 2.2

#correct a few sampling year --> when looking at org csv
kevin$sampling_year[is.na(kevin$sampling_year)==TRUE & kevin$site=="brandbjerg"] <- 2011
kevin$sampling_year[is.na(kevin$sampling_year)==TRUE & kevin$site=="popface"] <- 2002
kevin$sampling_year[is.na(kevin$sampling_year)==TRUE & kevin$site=="biocon"] <- 2005
kevin$sampling_year[kevin$sampling_year=="2005-2010"] <- 2008
kevin$sampling_year[kevin$sampling_year=="1996-2010"] <- 2003
kevin$sampling_year[kevin$sampling_year=="2003-2006"] <- 2005
kevin$sampling_year[kevin$citation=="domec_et_al_2012"] <- 1997 # this is so weried - but by looking at their info it says +1y. then we assume it is 1996+1
kevin$sampling_year[kevin$citation=="ellsworth_et_al_2012"] <- 1997
kevin$start_year[is.na(kevin$start_year)==TRUE] <- 1992

kevin$sampling_year <- as.numeric(kevin$sampling_year)
kevin$start_year <- as.numeric(kevin$start_year)
kevin$Year <-  kevin$sampling_year - kevin$start_year
summary(kevin$Year)

#correct elevation
aaa <- aggregate(kevin,by=list(kevin$lon,kevin$lat), FUN=mean, na.rm=TRUE)[,c("lon","lat")]
aaa$sitename <- paste("c",1:length(aaa$lon),sep="")
devtools::load_all("~/yunkepeng/gcme/pmodel/ingestr/")
df_etopo <- ingest(aaa,source = "etopo1",dir = "~/data/etopo/" )
aaa$elv <- as.numeric(as.data.frame(df_etopo$data))
kevin_z <- merge(kevin,aaa[,c("lon","lat","elv")],by=c("lon","lat"),all.x=TRUE)
plot(kevin_z$z~kevin_z$elv) # looks ok - now interploate original elevation value with etopo elevation
kevin_z$z[is.na(kevin_z$z)==TRUE] <- kevin_z$elv[is.na(kevin_z$z)==TRUE] 
kevin_z <- kevin_z[, !(colnames(kevin_z) %in% c("elv"))]

summary(kevin_z$Year)
summary(kevin_z$start_year)
summary(kevin_z$sampling_year) #all from 1992 to 2016 - good news!

kevin_z$year_start <- kevin_z$start_year
kevin_z$year_end <- kevin_z$sampling_year


response_ratio <- function(df){
  #-----------------------------------------------------------------------
  # Input: 
  # name: df..
  # df should include: c("exp_nam","ambient","ambient_Sd","ambient_Se",
  # "elevated","elevated_Se","elevated_Sd","n_plots")
  #-----------------------------------------------------------------------
  df_c_sub <- df %>%         
    ## Here only for my variables, selected as described above 
    # get standard deviation for all data
    mutate( my_ambient_sd = ambient_Sd, my_elevated_sd = elevated_Sd ) %>%
    rowwise() %>% 
    mutate( my_ambient_sd   = ifelse( is.na(my_ambient_sd),  ambient_Se  * sqrt(n_plots), my_ambient_sd ),
            my_elevated_sd  = ifelse( is.na(my_elevated_sd), elevated_Se * sqrt(n_plots), my_elevated_sd )) %>%
    
    ## Get logarithm of response ratio and its variance
    metafor::escalc( 
      measure = "ROM", 
      m1i = elevated, sd1i = my_elevated_sd, n1i = n_plots, 
      m2i = ambient,  sd2i = my_ambient_sd,  n2i = n_plots, 
      data=., 
      append = TRUE, var.names = c("logr", "logr_var") ) %>% 
    as_tibble() %>% 
    mutate( logr_se = sqrt(logr_var)/sqrt(n_plots) )
  
  df_c_sub$id <- 1:nrow(df_c_sub)
  return(df_c_sub)
  #-----------------------------------------------------------------------
  # Output: df_c_sub
  #-----------------------------------------------------------------------
}
kevin2 <- response_ratio(kevin_z)
kevin2$sen_coef <- kevin2$logr/log(kevin2$co2_e/kevin2$co2_a)
kevin2$sen_coef_var <- kevin2$logr_var/log(kevin2$co2_e/kevin2$co2_a)
kevin2$sen_coef_se <- kevin2$logr_se/log(kevin2$co2_e/kevin2$co2_a)

#kevin2 - ecm type
kevin2$ecm_type <- NA
kevin2$ecm_type[kevin2$site=="biocon"] <- "AM";kevin2$ecm_type[kevin2$site=="brandbjerg"] <- "AM";
kevin2$ecm_type[kevin2$site=="giface"] <- "AM";kevin2$ecm_type[kevin2$site=="nevada_desert_face"] <- "AM";kevin2$ecm_type[kevin2$site=="new_zealand_face"] <- "AM";
kevin2$ecm_type[kevin2$site=="riceface_china_33n_120e"] <- "AM";kevin2$ecm_type[kevin2$exp_nam=="riceface_japan"] <- "AM";
kevin2$ecm_type[kevin2$site=="soyface"] <- "AM";kevin2$ecm_type[kevin2$site=="swissface_trifolium"] <- "AM"
kevin2$ecm_type[kevin2$site=="dukeface"] <- "ECM";kevin2$ecm_type[kevin2$site=="eucface"] <- "ECM";
kevin2$ecm_type[kevin2$site=="euroface"] <- "ECM";kevin2$ecm_type[kevin2$site=="facts_ii_face"] <- "ECM";kevin2$ecm_type[kevin2$site=="popface"] <- "ECM";
#now, start p-model
#combine with siteinfo
#see forcing_kevin_nick.R
siteinfo <- read.csv("/Users/yunpeng/data/gcme/kevin/forcing/forcing_info.csv")
siteinfo <- siteinfo[,c("lon","lat","elv","year_start","year_end","sitename")]
names(siteinfo) <- c("lon","lat","z","year_start","year_end","sitename")
vcmax25_warmingco2_siteinfo <- merge(kevin2,siteinfo,by=c("lon","lat","z","year_start","year_end"),all.x=TRUE)

#now, start p-model
#devtools::load_all("~/yunkepeng/gcme/pmodel/ingestr/")

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

forcing_path <- "/Users/yunpeng/data/gcme/kevin/forcing/climate/"
vcmax25_warmingco2_siteinfo$elv <- vcmax25_warmingco2_siteinfo$z

for (i in 1:nrow(vcmax25_warmingco2_siteinfo)){
  siteinfo_site <- vcmax25_warmingco2_siteinfo[i,c("sitename","lon","lat","elv","year_start","year_end")]
  sitename <- siteinfo_site$sitename
  
  df1 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for ambient
  df2 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for elevated
  
  #Q1 - add warming directly
  df2$temp <- df2$temp+vcmax25_warmingco2_siteinfo$warmQ_e2[i]
  
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
summary(vcmax25_warmingco2_siteinfo)

#plot
kevin2_final <- vcmax25_warmingco2_siteinfo

#create a function to divide them into specific cateory: cw, cf
#input:df, data_type, logr and logr_var
kevin2_c_vcmax <- subset(kevin2_final, treatment=="c" & response =="vcmax")
kevin2_w_vcmax <- subset(kevin2_final, treatment=="w" & response =="vcmax")
kevin2_f_vcmax <- subset(kevin2_final, treatment=="f" & response =="vcmax")
kevin2_cf_vcmax <- subset(kevin2_final, treatment=="cw" & response =="vcmax")
kevin2_cw_vcmax <- subset(kevin2_final, treatment=="cf" & response =="vcmax")

kevin2_c_jmax <- subset(kevin2_final, treatment=="c" & response =="jmax")
kevin2_w_jmax <- subset(kevin2_final, treatment=="w" & response =="jmax")
kevin2_f_jmax <- subset(kevin2_final, treatment=="f" & response =="jmax")
kevin2_cf_jmax <- subset(kevin2_final, treatment=="cw" & response =="jmax")
kevin2_cw_jmax <- subset(kevin2_final, treatment=="cf" & response =="jmax")

#direct compare
ggplot(kevin2_c_vcmax,aes(x=ambient, y=mean_vcmax_ambient)) +geom_point(aes(color=co2_a),size=2)+geom_smooth(method="lm")+
  labs(x = "measured vcmax at ambient co2",y="measured vcmax at elevated co2")+theme_classic()+theme(axis.text=element_text(size=20),axis.title =element_text(size=20))

#aggregate into response ratio basing exp_nam mean (where Year is random factor)
agg_meta <- function(df){
  explist <- unique(df$exp_nam)
  mylist <- list() #create an empty list
  #some are just all NA for original data's standard deviation - needs division -> in this way just calculate normal mean, lower and upper band
  for (i in 1:length(explist)){ # if over half of samples doens't have variance then we should not rely on this anyways.
    if (sum(is.na(subset(df,exp_nam==explist[i])$logr_var),na.rm=TRUE)/length(is.na(subset(df,exp_nam==explist[i])$logr_var)) >0.5){
      #uncertainty <- 1.96*sd(as.vector(subset(df,exp_nam==explist[i])$logr))
      #uncertainty_50 <- 0.67449*sd(as.vector(subset(df,exp_nam==explist[i])$logr))
      uncertainty <- 0 # calculating sd directly from plots is horrible..
      uncertainty_50 <- 0
      mean_value <- mean(subset(df,exp_nam==explist[i])$logr,na.rm=TRUE)
      n_samples <- length(subset(df,exp_nam==explist[i])$logr_var)
      df_box <- tibble(
        exp_nam=explist[i], middle = mean_value,
        ymin   = mean_value-uncertainty,
        ymax   = mean_value+uncertainty,
        ymin_quarter = mean_value-uncertainty_50,
        ymax_quarter   = mean_value+uncertainty_50,
        middle_scaled = NA, ymin_scaled   = NA, ymax_scaled   = NA,variance_info ="No_mostNA",no=n_samples)
    } else {# define if uncertainty is 0, when N = 1 with sd
      if ( length(subset(df,exp_nam==explist[i])$logr_var) ==1){
        uncertainty <- 0
        mean_value <- mean(subset(df,exp_nam==explist[i])$logr,na.rm=TRUE)
        n_samples <- length(subset(df,exp_nam==explist[i])$logr_var)
        df_box <- tibble(
          exp_nam=explist[i], middle = mean_value,
          ymin   = mean_value-uncertainty,
          ymax   = mean_value+uncertainty,
          ymin_quarter = mean_value-uncertainty,
          ymax_quarter   = mean_value+uncertainty,
          middle_scaled = NA, ymin_scaled   = NA, ymax_scaled   = NA,variance_info ="No_1_number",no=n_samples)
      } else {
        #where year is random factor
        out_meta <- metafor::rma.mv( logr, logr_var, method = "REML",
                                     random = ~ 1 | Year, slab = Year, control = list(stepadj=0.3),data = subset(df,exp_nam==explist[i]))
      
        out_meta_quarter <- metafor::rma.mv( logr, logr_var, method = "REML",level=50,
                                   random = ~ 1 | Year, slab = Year, control = list(stepadj=0.3),data = subset(df,exp_nam==explist[i]))
        
        out_meta_scaled <- predict( out_meta, transf=exp )
        n_samples <- length(subset(df,exp_nam==explist[i])$logr_var)
        
        df_box <- tibble(
        exp_nam=explist[i], middle = out_meta$b[1,1], ymin   = out_meta$ci.lb, ymax   = out_meta$ci.ub,
        ymin_quarter=out_meta_quarter$ci.lb, ymax_quarter = out_meta_quarter$ci.ub,
        middle_scaled = out_meta_scaled$pred, ymin_scaled   = out_meta_scaled$ci.lb, ymax_scaled   = out_meta_scaled$ci.ub,variance_info ="Yes",no=n_samples)}}
    mylist[[i]] <- df_box}
    output <- do.call("rbind",mylist)
    return(output)
}

plotmean_c_logrvcmax <- agg_meta(kevin2_c_vcmax)
plotmean_w_logrvcmax <- agg_meta(kevin2_w_vcmax)
plotmean_f_logrvcmax <- agg_meta(kevin2_f_vcmax)
plotmean_cf_logrvcmax <- agg_meta(kevin2_cf_vcmax)
plotmean_cw_logrvcmax <- agg_meta(kevin2_cw_vcmax)

plotmean_c_logrjmax <- agg_meta(kevin2_c_jmax)
plotmean_w_logrjmax <- agg_meta(kevin2_w_jmax)
plotmean_f_logrjmax <- agg_meta(kevin2_f_jmax)
plotmean_cf_logrjmax <- agg_meta(kevin2_cf_jmax)
plotmean_cw_logrjmax <- agg_meta(kevin2_cw_jmax)


#aggregate into sensitivity coefficient basing exp_nam mean (where Year is random factor)
agg_meta_sen_coef <- function(df){
  explist <- unique(df$exp_nam)
  mylist <- list() #create an empty list
  #convert to sensitivity coef
  df$sen_coef <- df$logr/log(df$co2_e/df$co2_a)
  df$sen_coef_var <- df$logr_var/log(df$co2_e/df$co2_a)
  df$sen_coef_se <- df$logr_se/log(df$co2_e/df$co2_a)
  
  df$logr <- df$sen_coef
  df$logr_var <- df$sen_coef_var
  df$logr_se <- df$sen_coef_se
  #some are just all NA for original data's standard deviation - needs division -> in this way just calculate normal mean, lower and upper band
  for (i in 1:length(explist)){ # if over half of samples doens't have variance then we should not rely on this anyways.
    if (sum(is.na(subset(df,exp_nam==explist[i])$logr_var),na.rm=TRUE)/length(is.na(subset(df,exp_nam==explist[i])$logr_var)) >0.5){
      #uncertainty <- 1.96*sd(as.vector(subset(df,exp_nam==explist[i])$logr))
      #uncertainty_50 <- 0.67449*sd(as.vector(subset(df,exp_nam==explist[i])$logr))
      uncertainty <- 0 # calculating sd directly from plots is horrible..
      uncertainty_50 <- 0
      mean_value <- mean(subset(df,exp_nam==explist[i])$logr,na.rm=TRUE)
      n_samples <- length(subset(df,exp_nam==explist[i])$logr_var)
      df_box <- tibble(
        exp_nam=explist[i], middle = mean_value,
        ymin   = mean_value-uncertainty,
        ymax   = mean_value+uncertainty,
        ymin_quarter = mean_value-uncertainty_50,
        ymax_quarter   = mean_value+uncertainty_50,
        middle_scaled = NA, ymin_scaled   = NA, ymax_scaled   = NA,variance_info ="No_allNA",no=n_samples)
    } else {# define if uncertainty is 0, when N = 1
      if (length(subset(df,exp_nam==explist[i])$logr_var) ==1){
        uncertainty <- 0
        mean_value <- mean(subset(df,exp_nam==explist[i])$logr,na.rm=TRUE)
        n_samples <- length(subset(df,exp_nam==explist[i])$logr_var)
        df_box <- tibble(
          exp_nam=explist[i], middle = mean_value,
          ymin   = mean_value-uncertainty,
          ymax   = mean_value+uncertainty,
          ymin_quarter = mean_value-uncertainty,
          ymax_quarter   = mean_value+uncertainty,
          middle_scaled = NA, ymin_scaled   = NA, ymax_scaled   = NA,variance_info ="No_1_number",no=n_samples)
      } else {
        #where year is random factor
        out_meta <- metafor::rma.mv( logr, logr_var, method = "REML",
                                     random = ~ 1 | Year, slab = Year, control = list(stepadj=0.3),data = subset(df,exp_nam==explist[i]))
        
        out_meta_quarter <- metafor::rma.mv( logr, logr_var, method = "REML",level=50,
                                             random = ~ 1 | Year, slab = Year, control = list(stepadj=0.3),data = subset(df,exp_nam==explist[i]))
        
        out_meta_scaled <- predict( out_meta, transf=exp )
        n_samples <- length(subset(df,exp_nam==explist[i])$logr_var)
        
        df_box <- tibble(
          exp_nam=explist[i], middle = out_meta$b[1,1], ymin   = out_meta$ci.lb, ymax   = out_meta$ci.ub,
          ymin_quarter=out_meta_quarter$ci.lb, ymax_quarter = out_meta_quarter$ci.ub,
          middle_scaled = out_meta_scaled$pred, ymin_scaled   = out_meta_scaled$ci.lb, ymax_scaled   = out_meta_scaled$ci.ub,variance_info ="Yes",no=n_samples)}}
    mylist[[i]] <- df_box}
  output <- do.call("rbind",mylist)
  return(output)
}

plotmean_c_coefvcmax <- agg_meta_sen_coef(kevin2_c_vcmax)
plotmean_c_coefjmax <- agg_meta_sen_coef(kevin2_c_jmax)


#now, shows the shift
look2_combination <- look2 %>%  filter(exp_nam %in% c("riceface_japan","brandbjerg","new_zealand_face","setres","dukeface","euroface"))

look2_cf <- look2_cf[,c("exp_nam","middle","ymin","ymax","no")]
look2_cw <- look2_cw[,c("exp_nam","middle","ymin","ymax","no")]
names(look2_cf) <- c("exp_nam","middle_cf","ymin_cf","ymax_cf","no_cf")
names(look2_cw) <- c("exp_nam","middle_cw","ymin_cw","ymax_cw","no_cw")
look2_final <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp_nam"),all.x=TRUE),
                     list(look2_combination,look2_cf,look2_cw))
look2_final %>%
  ggplot() +
  geom_crossbar(aes(x=exp_nam, y=middle, ymin=ymin, ymax=ymax), alpha = 1, width = 0.5 ) +
  geom_crossbar(aes(x=exp_nam, y=middle_cf, ymin=ymin_cf, ymax=ymax_cf), alpha = 0.2, width = 0.5, color="blue") +
  geom_crossbar(aes(x=exp_nam, y=middle_cw, ymin=ymin_cw, ymax=ymax_cw), alpha = 0.2, width = 0.5,color="red" ) +
  geom_hline( yintercept=0.0, size=0.5 )+ylim(-3,3)+
  labs(x="", y="vcmax25 from combination effects
       black (co2), blue (co2 + Nfer), red (co2 + warming)", size=expression(paste("Error"^{-1}))) +
  coord_flip() 

#later...

#direct comparsion
look2
pred_vc_plotmean <- aggregate(kevin2_c_vcmax,by=list(kevin2_c_vcmax$exp_nam), FUN=mean, na.rm=TRUE)[,c("Group.1","pred_vcmax25_coef")]
names(pred_vc_plotmean) <- c("exp_nam","pred_vcmax25_coef")
vcmax25_plotmean <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp_nam"),all.x=TRUE),
                     list(look2,pred_vc_plotmean))
vcmax25_plotmean$ecm_type <- NA
vcmax25_plotmean$ecm_type[vcmax25_plotmean$exp_nam=="dukeface"] <- "ECM";vcmax25_plotmean$ecm_type[vcmax25_plotmean$exp_nam=="eucface"] <- "ECM";vcmax25_plotmean$ecm_type[vcmax25_plotmean$exp_nam=="euroface"] <- "ECM";vcmax25_plotmean$ecm_type[vcmax25_plotmean$exp_nam=="facts_ii_face"] <- "ECM";vcmax25_plotmean$ecm_type[vcmax25_plotmean$exp_nam=="popface"] <- "ECM";
vcmax25_plotmean$ecm_type[is.na(vcmax25_plotmean$ecm_type)==TRUE] <- "AM"

ggplot(vcmax25_plotmean,aes(x=middle, y=pred_vcmax25_coef)) +
  geom_point(aes(shape=ecm_type),size=2)+
  geom_smooth(method="lm")+
  labs(x = "measured sensitivity coef",y="predicted sensitivity coef")+theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

pred_j_plotmean <- aggregate(kevin2_c_jmax,by=list(kevin2_c_jmax$exp_nam), FUN=mean, na.rm=TRUE)[,c("Group.1","pred_jmax25_coef")]
names(pred_j_plotmean) <- c("exp_nam","pred_jmax25_coef")
jmax25_plotmean <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp_nam"),all.x=TRUE),
                          list(look2j,pred_j_plotmean))
jmax25_plotmean$ecm_type <- NA
jmax25_plotmean$ecm_type[jmax25_plotmean$exp_nam=="dukeface"] <- "ECM";jmax25_plotmean$ecm_type[jmax25_plotmean$exp_nam=="eucface"] <- "ECM";jmax25_plotmean$ecm_type[jmax25_plotmean$exp_nam=="euroface"] <- "ECM";jmax25_plotmean$ecm_type[jmax25_plotmean$exp_nam=="facts_ii_face"] <- "ECM";jmax25_plotmean$ecm_type[jmax25_plotmean$exp_nam=="popface"] <- "ECM";
jmax25_plotmean$ecm_type[is.na(jmax25_plotmean$ecm_type)==TRUE] <- "AM"

ggplot(jmax25_plotmean,aes(x=middle, y=pred_jmax25_coef)) +
  geom_point(aes(shape=ecm_type),size=2)+
  geom_smooth(method="lm")+
  labs(x = "measured sensitivity coef",y="predicted sensitivity coef")+theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))



kevin2_c_vcmax %>%
  ggplot( aes(x=exp_nam, y=sen_coef)) +
  geom_jitter(alpha=0.5, aes(size = 1/sen_coef_se,color=ecm_type)) +
  geom_boxplot(aes(x=exp_nam, y=pred_vcmax25_coef),color="red",size=1)+ 
  geom_crossbar(data=look2,aes(x=exp_nam, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5 ) +
  geom_hline( yintercept=0.0, size=0.5 )+ylim(-2,2)+
  labs(x="", y="CO2 effects on vcmax25: sensitivity coefficients", size=expression(paste("Error"^{-1}))) +
  coord_flip() 

kevin2_c_vcmax %>%
  ggplot( aes(x=ecm_type, y=sen_coef)) +
  geom_jitter(alpha=1) +
  geom_hline( yintercept=0.0, size=0.5 )+
  geom_boxplot()+geom_boxplot(aes(x=ecm_type, y=pred_vcmax25_coef),color="red")+ylim(-2,2)+
  labs(x="", y="sensitivity coefficient of vcmax25 -  gcme data ") +
  coord_flip() 

look_warming <- agg_meta(kevin2_w_vcmax)
kevin2_w_vcmax %>%
  ggplot( aes(x=exp_nam, y=logr)) +
  geom_jitter(alpha=0.2, aes(size = 1/logr_se)) +
  geom_boxplot(aes(x=exp_nam, y=pred_vcmax25_logr),color="red",size=1)+ 
  geom_crossbar(data=look_warming,aes(x=exp_nam, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5 ) +
  geom_hline( yintercept=0.0, size=0.5 )+
  labs(x="", y="warming effects on vcmax25", size=expression(paste("Error"^{-1}))) +
  coord_flip() 

look_warming_j <- agg_meta(kevin2_w_jmax)
kevin2_w_jmax %>%
  ggplot( aes(x=exp_nam, y=logr)) +
  geom_jitter(alpha=0.2, aes(size = 1/logr_se)) +
  geom_boxplot(aes(x=exp_nam, y=pred_vcmax25_logr),color="red",size=1)+ 
  geom_crossbar(data=look_warming_j,aes(x=exp_nam, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5 ) +
  geom_hline( yintercept=0.0, size=0.5 )+
  labs(x="", y="warming effects on jmax25", size=expression(paste("Error"^{-1}))) +
  coord_flip() 

look_nfer <- agg_meta(kevin2_f_vcmax)
kevin2_f_vcmax %>%
  ggplot( aes(x=exp_nam, y=logr)) +
  geom_jitter(alpha=0.2, aes(size = 1/logr_se)) +
  geom_crossbar(data=look_nfer,aes(x=exp_nam, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5 ) +
  geom_hline( yintercept=0.0, size=0.5 )+
  labs(x="", y="N fertilication on vcmax25", size=expression(paste("Error"^{-1}))) +
  coord_flip() 

#jmax
look3 <- agg_meta_sen_coef(kevin2_c_jmax)
kevin2_c_jmax %>%
  ggplot( aes(x=exp_nam, y=sen_coef)) +
  geom_jitter(alpha=0.5, aes(size = 1/sen_coef_se,color=ecm_type)) +
  geom_boxplot(aes(x=exp_nam, y=pred_jmax25_coef),color="red",size=1)+ 
  geom_crossbar(data=look3,aes(x=exp_nam, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5 ) +
  geom_hline( yintercept=0.0, size=0.5 )+
  labs(x="", y="CO2 effects on jmax25: sensitivity coefficients", size=expression(paste("Error"^{-1}))) +
  coord_flip() 

kevin2_c_jmax %>%
  ggplot( aes(x=ecm_type, y=sen_coef)) +
  geom_jitter(alpha=1) +
  geom_hline( yintercept=0.0, size=0.5 )+
  geom_boxplot()+geom_boxplot(aes(x=ecm_type, y=pred_jmax25_coef),color="red")+ylim(-2,2)+
  labs(x="", y="sensitivity coefficient of jmax25 -  gcme data ") +
  coord_flip() 

look_warming2 <- agg_meta(kevin2_w_jmax)
kevin2_w_jmax %>%
  ggplot( aes(x=exp_nam, y=logr)) +
  geom_jitter(alpha=0.2, aes(size = 1/logr_se)) +
  geom_boxplot(aes(x=exp_nam, y=pred_jmax25_logr),color="red",size=1)+ 
  geom_crossbar(data=look_warming2,aes(x=exp_nam, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5 ) +
  geom_hline( yintercept=0.0, size=0.5 )+
  labs(x="", y="warming effects on jmax", size=expression(paste("Error"^{-1}))) +
  coord_flip() 

#look_nfer <- agg_meta(kevin2_f_jmax)
kevin2_f_jmax %>%
  ggplot( aes(x=exp_nam, y=logr)) +
  geom_jitter(alpha=0.2) +
  geom_boxplot()+
  geom_hline( yintercept=0.0, size=0.5 )+
  labs(x="", y="N fertilication on jmax25", size=expression(paste("Error"^{-1}))) +
  coord_flip() 
