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

#update function (if SD shown NA in plot - then relying on percentage of other samples within the same plot)
#AND --> convert logr or logr_var = 0 is NA
# first - response ratio of each
response_ratio_v2 <- function(df){
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
  df_c_sub$logr[df_c_sub$logr==0] <- NA
  df_c_sub$logr_var[df_c_sub$logr_var==0] <- NA
  
  explist <- unique(df_c_sub$exp)
  df_c_sub$new_logr_var_percentage <- NA
  for (i in 1:length(explist)){ # if SD shown NA in plot - then relying on percentage of other samples within the same plot
    percentage <- mean(abs(subset(df_c_sub,exp==explist[i])$logr_var/subset(df_c_sub,exp==explist[i])$logr), na.rm=TRUE)
    df_c_sub$new_logr_var_percentage[df_c_sub$exp==explist[i]] <- percentage
  }
  
  df_c_sub$logr_var[is.na(df_c_sub$logr_var)==TRUE] <- abs(df_c_sub$new_logr_var_percentage[is.na(df_c_sub$logr_var)==TRUE] * df_c_sub$logr[is.na(df_c_sub$logr_var)==TRUE])
  df_c_sub$logr_var[df_c_sub$logr_var=="NaN"] <- NA
  df_c_sub$logr_var[is.na(df_c_sub$logr)==TRUE] <- NA
  return(df_c_sub)
  #-----------------------------------------------------------------------
  # Output: df_c_sub
  #-----------------------------------------------------------------------
}

#for sensitivity coef
agg_meta_sen_coef <- function(df){
  threshold <- 0.5
  explist <- unique(df$exp)
  mylist <- list() #create an empty list
  #convert to sensitivity coef
  df$sen_coef <- df$logr/log(df$co2_e/df$co2_a)
  df$sen_coef_var <- df$logr_var/log(df$co2_e/df$co2_a)
  df$sen_coef_se <- df$logr_se/log(df$co2_e/df$co2_a)
  
  df$logr <- df$sen_coef
  df$logr_var <- df$sen_coef_var
  df$logr_se <- df$sen_coef_se
  #all using mean value
  for (i in 1:length(explist)){
    
    mean_value <- mean(subset(df,exp==explist[i])$logr,na.rm=TRUE)
    n_samples <- length(subset(df,exp==explist[i])$logr_var)
    df_box <- tibble(
      exp=explist[i], middle = mean_value,
      ymin   = NA,
      ymax   = NA,
      ymin_quarter = NA,
      ymax_quarter   = NA,
      variance_info ="No",no=n_samples,logr_var=NA)
    mylist[[i]] <- df_box}
  output <- do.call("rbind",mylist)
  return(output)
}

#for response ratio
agg_meta <- function(df){
  threshold <- 0.5
  explist <- unique(df$exp)
  mylist <- list() 
  for (i in 1:length(explist)){
    mean_value <- mean(subset(df,exp==explist[i])$logr,na.rm=TRUE)
    n_samples <- length(subset(df,exp==explist[i])$logr_var)
    df_box <- tibble(
      exp=explist[i], middle = mean_value,
      ymin   = NA,
      ymax   = NA,
      ymin_quarter = NA,
      ymax_quarter   = NA,
      variance_info ="No",no=n_samples,logr_var = NA)
    mylist[[i]] <- df_box}
  output <- do.call("rbind",mylist)
  return(output)
}

#combine c and cf
combine_co2_cf <- function(logr_c_var,logr_f_var,logr_cf_var,name){
  all_logr_c_lma <- agg_meta_sen_coef(logr_c_var)[,c("exp","middle")]
  all_logr_c_lma$condition <- "co2"
  #co2 + fer / fer
  all_logr_f_lma <- agg_meta(logr_f_var)[,c("exp","middle")]
  all_logr_cf_lma <- agg_meta(logr_cf_var)[,c("exp","middle")]
  all_logr_cf_lma$exp_old <- all_logr_cf_lma$exp 
  all_logr_cf_lma$exp <-  paste(substr(all_logr_cf_lma$exp,1,nchar(all_logr_cf_lma$exp)-2),"f",sep="")
  all_logr_f_lma$exp;all_logr_cf_lma$exp
  all_logr_f_lma2 <- all_logr_f_lma %>% filter(exp %in%all_logr_cf_lma$exp)
  new_df <- merge(all_logr_f_lma2,all_logr_cf_lma,by=c("exp"),all.x=TRUE)
  for (i in 1:nrow(new_df)) {
    new_df$co2[i] <- subset(logr_cf_var,exp==new_df$exp_old[i])$co2_e[1]-
      subset(logr_cf_var,exp==new_df$exp_old[i])$co2_a[1]
  } 
  new_df$middle <- (new_df$middle.y - new_df$middle.x)/log(new_df$co2)
  new_df$exp <- new_df$exp_old
  new_df$condition <- "Fertilization"
  
  lma_plot <- dplyr::bind_rows(all_logr_c_lma,new_df[,c("exp","middle","condition")]) 
  names(lma_plot) <- c("exp",name,"condition")
  return(lma_plot)
}

#aggrgate within pft
agg_meta_plots <- function(df,type_name){
  df$logr <- df$middle
  df$logr_var <- df$logr_var
  
  explist <- unique(eval(parse(text=paste("df$", type_name, sep = ""))))
  mylist <- list() #create an empty list
  
  for (i in 1:length(explist)){
    out_meta <- metafor::rma.mv( logr, logr_var, method = "REML", random = ~ 1 | exp_nam, slab = exp_nam,  control = list(stepadj=0.3), 
                                 data= subset(df,eval(parse(text=paste(type_name)))==explist[i]))
    out_meta_quarter <- metafor::rma.mv( logr, logr_var, method = "REML",random = ~ 1 | exp_nam, slab = exp_nam, level=50,control = list(stepadj=0.3),
                                         data = subset(df,eval(parse(text=paste(type_name)))==explist[i]))
    n_plots <- nrow(subset(df,eval(parse(text=paste(type_name)))==explist[i]))
    
    df_box <- tibble(
      type_name=explist[i], middle = out_meta$b[1,1], ymin   = out_meta$ci.lb, ymax   = out_meta$ci.ub,
      ymin_quarter=out_meta_quarter$ci.lb, ymax_quarter = out_meta_quarter$ci.ub,
      no_plots = n_plots)
    mylist[[i]] <- df_box
  }
  output <- do.call("rbind",mylist)
  return(output)
}

#read Kevin
kevin <- read.csv("/Users/yunpeng/data/gcme/kevin/orig_vcmax/JunkePeng_11252021.csv")
#new <- subset(kevin)%>% group_by(site,lat,lon)  %>% summarise(number = n());new
#repeated sitename: itatinga,riceface_japan_ko_3558_13960,riceface_japan_ta_3558_13960
#kevin_final <- subset(kevin)%>% group_by(site,dominant_species)  %>% summarise(number = n())

#1. correct exp_nam to make it consistent with GCME
kevin$site[is.na(kevin$site)==TRUE] <- "riceface_japan_a_2003_3938_14057"
kevin$exp_nam <- kevin$site #new plot: capo_caccia,itatinga,itatinga_experimental_station,setres
#unify sitename from japan 
#kevin[grep("riceface_japan", kevin$exp_nam),]$exp_nam <- "riceface_japan"

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

#update function (if SD shown NA in plot - then relying on percentage of other samples within the same plot)
#AND --> convert logr or logr_var = 0 is NA
# first - response ratio of each
kevin2 <- response_ratio_v2(kevin_z)

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
#for conversion
vcmax25_warmingco2_siteinfo$warmQ_e2[vcmax25_warmingco2_siteinfo$exp=="riceface_japan_l_2008_3938_14057_w"] <- vcmax25_warmingco2_siteinfo$warmQ_e3[vcmax25_warmingco2_siteinfo$exp=="riceface_japan_l_2008_3938_14057_w"] 

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
  vcmax25_warmingco2_siteinfo$pred_vcmax_logr[i] <- log(vcmax25_warmingco2_siteinfo$mean_vcmax_elevated[i]/vcmax25_warmingco2_siteinfo$mean_vcmax_ambient[i])
  vcmax25_warmingco2_siteinfo$pred_vcmax_coef[i] <- (log(vcmax25_warmingco2_siteinfo$mean_vcmax_elevated[i]/vcmax25_warmingco2_siteinfo$mean_vcmax_ambient[i]))/(log(vcmax25_warmingco2_siteinfo$co2_e[i]/vcmax25_warmingco2_siteinfo$co2_a[i]))
  vcmax25_warmingco2_siteinfo$pred_vcmax25_logr[i] <- log(vcmax25_warmingco2_siteinfo$mean_vcmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_vcmax25_ambient[i])
  vcmax25_warmingco2_siteinfo$pred_vcmax25_coef[i] <- (log(vcmax25_warmingco2_siteinfo$mean_vcmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_vcmax25_ambient[i]))/(log(vcmax25_warmingco2_siteinfo$co2_e[i]/vcmax25_warmingco2_siteinfo$co2_a[i]))
  
  vcmax25_warmingco2_siteinfo$mean_jmax25_ambient[i] <- mean_jmax25_ambient
  vcmax25_warmingco2_siteinfo$mean_jmax25_elevated[i] <- mean_jmax25_elevated
  vcmax25_warmingco2_siteinfo$mean_jmax_ambient[i] <- mean_jmax_ambient
  vcmax25_warmingco2_siteinfo$mean_jmax_elevated[i] <- mean_jmax_elevated
  vcmax25_warmingco2_siteinfo$pred_jmax_logr[i] <- log(vcmax25_warmingco2_siteinfo$mean_jmax_elevated[i]/vcmax25_warmingco2_siteinfo$mean_jmax_ambient[i])
  vcmax25_warmingco2_siteinfo$pred_jmax_coef[i] <- (log(vcmax25_warmingco2_siteinfo$mean_jmax_elevated[i]/vcmax25_warmingco2_siteinfo$mean_jmax_ambient[i]))/(log(vcmax25_warmingco2_siteinfo$co2_e[i]/vcmax25_warmingco2_siteinfo$co2_a[i]))
  vcmax25_warmingco2_siteinfo$pred_jmax25_logr[i] <- log(vcmax25_warmingco2_siteinfo$mean_jmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_jmax25_ambient[i])
  vcmax25_warmingco2_siteinfo$pred_jmax25_coef[i] <- (log(vcmax25_warmingco2_siteinfo$mean_jmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_jmax25_ambient[i]))/(log(vcmax25_warmingco2_siteinfo$co2_e[i]/vcmax25_warmingco2_siteinfo$co2_a[i]))
  
  vcmax25_warmingco2_siteinfo$PPFD[i] <- mean(df1$ppfd,na.rm=TRUE)*1000000 # in umol/m2/s
  vcmax25_warmingco2_siteinfo$vpd[i] <- mean(df1$vpd,na.rm=TRUE)/1000 #in kPa
  vcmax25_warmingco2_siteinfo$temp[i] <- mean(df1$temp,na.rm=TRUE) #in degree celcius
}
summary(vcmax25_warmingco2_siteinfo)

#plot
kevin2_final <- vcmax25_warmingco2_siteinfo

csvfile <- paste("/Users/yunpeng/data/gcme/kevin/forcing/pred_vcmax_v2.csv")
write.csv(kevin2_final, csvfile, row.names = TRUE)

#create a function to divide them into specific cateory: cw, cf
#input:df, data_type, logr and logr_var

#show pft
kevin2_final$ecosystem[kevin2_final$ecosystem=="temperate_forest"] <- "forest"
kevin2_final$ecosystem[kevin2_final$ecosystem=="heathland"] <- "grassland"
kevin2_final$ecosystem[kevin2_final$ecosystem=="shrubland"] <- "forest"

kevin2_c_vcmax <- subset(kevin2_final, treatment=="c" & response =="vcmax")
kevin2_f_vcmax <- subset(kevin2_final, treatment=="f" & response =="vcmax")
kevin2_cf_vcmax <- subset(kevin2_final, treatment=="cf" & response =="vcmax")
kevin2_w_vcmax <- subset(kevin2_final, treatment=="w" & response =="vcmax")

kevin2_c_jmax <- subset(kevin2_final, treatment=="c" & response =="jmax")
kevin2_f_jmax <- subset(kevin2_final, treatment=="f" & response =="jmax")
kevin2_cf_jmax <- subset(kevin2_final, treatment=="cf" & response =="jmax")
kevin2_w_jmax <- subset(kevin2_final, treatment=="w" & response =="jmax")

kevin_vcmax_plotmean <- agg_meta_sen_coef(kevin2_c_vcmax)

kevin_jmax_plotmean <- agg_meta_sen_coef(kevin2_c_jmax)

kevin_vcmax_w_plotmean <- agg_meta(kevin2_w_vcmax)
kevin_jmax_w_plotmean <- agg_meta(kevin2_w_jmax)

#warming - finanlise figures
pred_vc_w_kevin <- aggregate(kevin2_w_vcmax,by=list(kevin2_w_vcmax$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","pred_vcmax_logr")];names(pred_vc_w_kevin) <- c("exp","pred_vcmax")
pred_j_w_kevin <- aggregate(kevin2_w_jmax,by=list(kevin2_w_jmax$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","pred_jmax_logr")];names(pred_j_w_kevin) <- c("exp","pred_jmax")

kevin_vj_w_final <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),
                           list(pred_vc_w_kevin,pred_j_w_kevin,kevin_vcmax_w_plotmean[,c("exp","middle")],kevin_jmax_w_plotmean[,c("exp","middle","no")]))
names(kevin_vj_w_final) <- c("exp","pred_vcmax","pred_jmax","obs_vcmax","obs_jmax","no")

kevin_vj_w_final$pred_jmax_vcmax <- kevin_vj_w_final$pred_jmax- kevin_vj_w_final$pred_vcmax
kevin_vj_w_final$obs_jmax_vcmax <- kevin_vj_w_final$obs_jmax- kevin_vj_w_final$obs_vcmax
kevin_vj_w_final$ecosystem <- "grassland"

csvfile <- paste("/Users/yunpeng/data/gcme/kevin_20220222/vcmaxjmax_warming_plotmean_v2.csv")
write.csv(kevin_vj_w_final, csvfile, row.names = TRUE)




#now, show prediction
pred_vc_kevin <- aggregate(kevin2_c_vcmax,by=list(kevin2_c_vcmax$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","pred_vcmax_coef")];names(pred_vc_kevin) <- c("exp","pred_vcmax")
pred_j_kevin <- aggregate(kevin2_c_jmax,by=list(kevin2_c_jmax$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","pred_jmax_coef")];names(pred_j_kevin) <- c("exp","pred_jmax")


kevin_ecosystem <- as.data.frame(kevin2_c_vcmax %>% group_by(exp,ecosystem) %>% summarise(number=n()))

kevin_vcmax_final <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),
                     list(kevin_vcmax_plotmean,pred_vc_kevin,kevin_ecosystem[,c("exp","ecosystem")]))

kevin_jmax_final <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),
                           list(kevin_jmax_plotmean,pred_j_kevin,kevin_ecosystem[,c("exp","ecosystem")]))
kevin_jmax_final$ecosystem[is.na(kevin_jmax_final$ecosystem)==TRUE] <- "cropland" # two riceface are cropland

summary(kevin_vcmax_final)
summary(kevin_jmax_final)

#Smith's data
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

vcmax25_warmingco2_siteinfo$ecosystem <- NA
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

csvfile <- paste("/Users/yunpeng/data/smith_keenan_gcb/gcb_co2/pred_vcmax_v2.csv")
write.csv(vcmax25_warmingco2_siteinfo, csvfile, row.names = TRUE)

#now, combine with Kevin:
smith_all_removal <- subset(vcmax25_warmingco2_siteinfo,exp_nam!="NZFACE" &exp_nam!="BioCON" & exp_nam!="ChinaRiceFACE" & exp_nam!="DukeFACE" & exp_nam!="EUROPOPFACE" & exp_nam!="NevadaFACE" & exp_nam!="SwissFACE")
smith_all_simple <- smith_all_removal[,c("SiteID","ecosystem","sen_coef_v","sen_coef_j","pred_vcmax25_coef","pred_jmax25_coef")]

smith_all_plotmean <- aggregate(smith_all_simple,by=list(smith_all_simple$SiteID,smith_all_simple$ecosystem), FUN=mean, na.rm=TRUE)[,c("Group.1","Group.2","sen_coef_v","sen_coef_j","pred_vcmax25_coef","pred_jmax25_coef")]
names(smith_all_plotmean) <- c("exp","ecosystem","vcmax","jmax","pred_vcmax","pred_jmax")
smith_all_plotmean$exp[smith_all_plotmean$exp=="OakOTC"] <- "mi_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="ORNL"] <- "ornerp_liqui_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="PineOTC"] <- "Christchurch_pr_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="RichmondGH"] <- "australiacotton_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="UIAPine"] <- "ua_otc_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="PHACE"] <- "phace_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="AspenFACE"] <- "rhine-aspenface_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="BilyKriz"] <- "bily_kriz_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="Headley"] <- "headley_qp_c";
smith_all_plotmean$exp[smith_all_plotmean$exp=="Viesalm"] <- "vielsalm_c"

smith_all_plotmean$exp <- tolower(smith_all_plotmean$exp)
smith_all_plotmean$variance_info <- "smith"

smith_vcmax_final <- smith_all_plotmean[,c("exp","ecosystem","vcmax","pred_vcmax","variance_info")]
names(smith_vcmax_final) <- c("exp","ecosystem","middle","pred_vcmax","variance_info")

smith_jmax_final <- smith_all_plotmean[,c("exp","ecosystem","jmax","pred_jmax","variance_info")]
names(smith_jmax_final) <- c("exp","ecosystem","middle","pred_jmax","variance_info")

#finally, combine_co2 (c, cf, f)
Nfer_vcmax <- na.omit(subset(combine_co2_cf(kevin2_c_vcmax,kevin2_f_vcmax,kevin2_cf_vcmax,"middle"),condition=="Fertilization"))
Nfer_vcmax$ecosystem <- "forest"
Nfer_vcmax$ecosystem[Nfer_vcmax$exp=="new_zealand_face_cf"] <- "grassland"
Nfer_vcmax$variance_info <- "No"
Nfer_jmax <- na.omit(subset(combine_co2_cf(kevin2_c_jmax,kevin2_f_jmax,kevin2_cf_jmax,"middle"),condition=="Fertilization"))
Nfer_jmax$ecosystem <- "forest"
Nfer_jmax$variance_info <- "No"

#now, combine Nick and Kevin
vcmax_final <- dplyr::bind_rows(kevin_vcmax_final,smith_vcmax_final,Nfer_vcmax)
jmax_final <- dplyr::bind_rows(kevin_jmax_final,smith_jmax_final,Nfer_jmax)

csvfile <- paste("/Users/yunpeng/data/gcme/kevin_20220222/vcmax_plotmean_v2.csv")
write.csv(vcmax_final, csvfile, row.names = TRUE)

csvfile <- paste("/Users/yunpeng/data/gcme/kevin_20220222/jmax_plotmean_v2.csv")
write.csv(jmax_final, csvfile, row.names = TRUE)


#finally, merging with pft
plot_name <- unique(kevin2_c_vcmax[,c("exp","exp_nam")]);plot_name$rep <- duplicated(plot_name$exp);plot_name <- subset(plot_name,rep==FALSE)[,c("exp","exp_nam")]
plot_j_name <- unique(kevin2_c_jmax[,c("exp","exp_nam")]);plot_j_name$rep <- duplicated(plot_j_name$exp);plot_j_name <- subset(plot_j_name,rep==FALSE)[,c("exp","exp_nam")]

vcmax_final_exp <- merge(vcmax_final,plot_name,by=c("exp"),all.x=TRUE)
jmax_final_exp <- merge(jmax_final,plot_j_name,by=c("exp"),all.x=TRUE)

vcmax_final_exp$logr_se <- sqrt(vcmax_final_exp$logr_var)/sqrt(vcmax_final_exp$no)
jmax_final_exp$logr_se <- sqrt(jmax_final_exp$logr_var)/sqrt(jmax_final_exp$no)


#Kevin + Smith
vcmax_all <- subset(vcmax_final,is.na(pred_vcmax)==FALSE)
vcmax_all$source <- "kevin"
vcmax_all$source[vcmax_all$variance_info=="smith"] <- "smith"

final_prediction_vc_all <- tibble(type_name="all",middle=median(vcmax_all$pred_vcmax),ymin=quantile(vcmax_all$pred_vcmax, 0.25),ymax=quantile(vcmax_all$pred_vcmax, 0.75))
final_observation_vc_all <- tibble(type_name="all",middle=median(vcmax_all$middle),ymin=quantile(vcmax_all$middle, 0.25),ymax=quantile(vcmax_all$middle, 0.75))

b1 <- vcmax_all %>%
  ggplot( aes(x=ecosystem, y=middle)) +
  geom_boxplot()+
  geom_crossbar(data=final_observation_vc_all,aes(x=type_name,y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="black") +
  geom_crossbar(data=final_prediction_vc_all,aes(x=type_name,y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="red") +
  #geom_point(aes(color=source),alpha = 0.6, width = 0.5,size=2) +
  geom_point(alpha = 0.6, width = 0.5,size=2) +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(x="", y="Vcmax response",size=expression(paste("Standard Error"^{-1}))) +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

jmax_all <- subset(jmax_final,is.na(pred_jmax)==FALSE &is.na(middle)==FALSE )
jmax_all$source <- "kevin"
jmax_all$source[jmax_all$variance_info=="smith"] <- "smith"

final_prediction_j_all <- tibble(type_name="all",middle=median(jmax_all$pred_jmax),ymin=quantile(jmax_all$pred_jmax, 0.25),ymax=quantile(jmax_all$pred_jmax, 0.75))
final_observation_j_all <- tibble(type_name="all",middle=median(jmax_all$middle),ymin=quantile(jmax_all$middle, 0.25),ymax=quantile(jmax_all$middle, 0.75))

b2 <- jmax_all %>%
  ggplot( aes(x=ecosystem, y=middle)) +
  geom_boxplot()+
  geom_crossbar(data=final_observation_j_all,aes(x=type_name,y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="black") +
  geom_crossbar(data=final_prediction_j_all,aes(x=type_name,y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="red") +
  #geom_point(aes(color=source),alpha = 0.6, width = 0.5,size=2) +
  geom_point(alpha = 0.6, width = 0.5,size=2) +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(x="", y="Jmax response",size=expression(paste("Standard Error"^{-1}))) +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

final_vj_all <- merge(vcmax_all[,c("exp","ecosystem","middle","source","pred_vcmax")],
                  jmax_all[,c("exp","ecosystem","middle","source","pred_jmax")], by = c("exp","ecosystem","source"),all.x=TRUE)
final_vj_all <- na.omit(final_vj_all)

final_vj_all$middle <- final_vj_all$middle.y-final_vj_all$middle.x
final_vj_all$pred_vj <- final_vj_all$pred_jmax-final_vj_all$pred_vcmax

final_prediction_vj_all <- tibble(type_name="all",middle=median(final_vj_all$pred_vj),ymin=quantile(final_vj_all$pred_vj, 0.25),ymax=quantile(final_vj_all$pred_vj, 0.75))
final_observation_vj_all <- tibble(type_name="all",middle=median(final_vj_all$middle),ymin=quantile(final_vj_all$middle, 0.25),ymax=quantile(final_vj_all$middle, 0.75))

b3 <- final_vj_all %>%
  ggplot( aes(x=ecosystem, y=middle)) +
  geom_boxplot()+
  geom_crossbar(data=final_observation_vj_all,aes(x=type_name,y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="black") +
  geom_crossbar(data=final_prediction_vj_all,aes(x=type_name,y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="red") +
  #geom_point(aes(color=source),alpha = 0.6, width = 0.5,size=2) +
  geom_point(alpha = 0.6, width = 0.5,size=2) +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(x="", y="Jmax/Vcmax response",size=expression(paste("Standard Error"^{-1}))) +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))
b3

plot_grid(b1,b2,b3,nrow=1,label_size = 15)+theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/final_fig1_alternative_v2.jpg",sep=""),width = 15, height = 5)

