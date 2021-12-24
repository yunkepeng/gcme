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

# first - response ratio of each
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

#read Kevin - other vars
kevin_othervars <- read.csv("/Users/yunpeng/data/gcme/kevin/orig_leaf/YunkePeng_othervars.csv")
kevin_othervars <- rename(kevin_othervars, c(ambient = x_c, elevated=x_t, ambient_Sd=sd_c, elevated_Sd=sd_t,ambient_Se=se_c,elevated_Se=se_t,n_plots=rep_c,
                         z=elevation, co2_a=c_c, co2_e=c_t, nfertQ_a = n_c, nfertQ_e = n_t, pfertQ_a = p_c, pfertQ_e = p_t,kfertQ_a = k_c, kfertQ_e = k_t,
                         warmQ_e1 = w_t1, warmQ_e2 = w_t2, warmQ_e3 = w_t3, Unit=x_units))
kevin_othervars$ambient_Se <- as.numeric(kevin_othervars$ambient_Se)
kevin_othervars$elevated_Se <- as.numeric(kevin_othervars$elevated_Se)

kevin_othervars$sampling_year[is.na(kevin_othervars$sampling_year)==TRUE & kevin_othervars$site=="brandbjerg"] <- 2011;kevin_othervars$sampling_year[is.na(kevin_othervars$sampling_year)==TRUE & kevin_othervars$site=="popface"] <- 2002;kevin_othervars$sampling_year[is.na(kevin_othervars$sampling_year)==TRUE & kevin_othervars$site=="biocon"] <- 2005;
kevin_othervars$sampling_year[kevin_othervars$sampling_year=="2005-2010"] <- 2008;kevin_othervars$sampling_year[kevin_othervars$sampling_year=="1996-2010"] <- 2003;kevin_othervars$sampling_year[kevin_othervars$sampling_year=="2003-2006"] <- 2005;
kevin_othervars$sampling_year[kevin_othervars$sampling_year=="1993-1994"] <- 1993;kevin_othervars$sampling_year[kevin_othervars$sampling_year=="1995-1998"] <- 1996;kevin_othervars$sampling_year[kevin_othervars$sampling_year=="1997_-_2007"] <- 2002;
kevin_othervars$sampling_year[kevin_othervars$sampling_year=="1998-2001"] <- 1999;kevin_othervars$sampling_year[kevin_othervars$sampling_year=="1999-2000"] <- 1999;
kevin_othervars$sampling_year[kevin_othervars$sampling_year=="2000-2001"] <- 2000;kevin_othervars$sampling_year[kevin_othervars$sampling_year=="2000-2003"] <- 2001;
kevin_othervars$sampling_year[kevin_othervars$sampling_year=="2001-2002"] <- 2001;kevin_othervars$sampling_year[kevin_othervars$sampling_year=="2001-2004"] <- 2002;
kevin_othervars$sampling_year[kevin_othervars$sampling_year=="2002-2003"] <- 2002;kevin_othervars$sampling_year[kevin_othervars$sampling_year=="2003-2004"] <- 2003;
kevin_othervars$sampling_year[kevin_othervars$sampling_year=="2004-2005"] <- 2004;kevin_othervars$sampling_year[kevin_othervars$sampling_year=="2005-'06"] <- 2005;
kevin_othervars$sampling_year[kevin_othervars$sampling_year=="2006-2008"] <- 2007;kevin_othervars$sampling_year[kevin_othervars$sampling_year=="2007-2009"] <- 2008;
kevin_othervars$sampling_year[kevin_othervars$sampling_year=="2010_&_2011"] <- 2010;kevin_othervars$sampling_year[kevin_othervars$sampling_year=="2011-2012"] <- 2011;
kevin_othervars$sampling_year[kevin_othervars$sampling_year=="2013-2014"] <- 2013;kevin_othervars$sampling_year[kevin_othervars$sampling_year=="2016?"] <- 2016;
kevin_othervars$sampling_year[kevin_othervars$sampling_year=="average_of_1998_-2006"] <- 2002;kevin_othervars$sampling_year[kevin_othervars$sampling_year=="compiled_between_1998_-_2012"] <- 2005;kevin_othervars$sampling_year[kevin_othervars$sampling_year=="compiled_between_2003-2007"] <- 2005
kevin_othervars$sampling_year[kevin_othervars$citation=="domec_et_al_2012"] <- 1997 # this is so weried - but by looking at their info it says +1y. then we assume it is 1996+1
kevin_othervars$sampling_year[kevin_othervars$citation=="ellsworth_et_al_2012"] <- 1997

kevin_othervars$sampling_year <- as.numeric(kevin_othervars$sampling_year)

kevin_othervars$Year <- kevin_othervars$sampling_year - kevin_othervars$start_year
summary(kevin_othervars$Year)
# some are negative, which is quite weired when looking at their original data - in this way we just assume they should be equal
kevin_othervars$Year[kevin_othervars$Year<0] <- 0

#THIS dangerous! However, since we don't know year of measurement, we could only assume it as average
kevin_othervars$Year[is.na(kevin_othervars$Year)==TRUE] <- 2

#some are still NA, leave it anyways.
kevin_LMA <- read.csv("/Users/yunpeng/data/gcme/kevin/orig_leaf/LMA.csv")
kevin_LMA <- rename(kevin_LMA, c(ambient = x_c, elevated=x_t, ambient_Sd=sd_c, elevated_Sd=sd_t,ambient_Se=se_c,elevated_Se=se_t,n_plots=rep_c,
                                   z=elevation, co2_a=c_c, co2_e=c_t, nfertQ_a = n_c, nfertQ_e = n_t, pfertQ_a = p_c, pfertQ_e = p_t,kfertQ_a = k_c, kfertQ_e = k_t,
                                   warmQ_e1 = w_t1, warmQ_e2 = w_t2, warmQ_e3 = w_t3, Unit=x_units))
kevin_LMA$ambient_Se <- as.numeric(kevin_LMA$ambient_Se)
kevin_LMA$elevated_Se <- as.numeric(kevin_LMA$elevated_Se)

kevin_LMA$sampling_year[kevin_LMA$sampling_year=="2001-2007"] <- 2004
kevin_LMA$sampling_year[kevin_LMA$sampling_year=="2011-2012"] <- 2011
kevin_LMA$sampling_year <- as.numeric(kevin_LMA$sampling_year)
kevin_LMA$Year <- kevin_LMA$sampling_year - kevin_LMA$start_year
kevin_LMA$Year[kevin_LMA$Year<0] <- 0
summary(kevin_LMA$Year)

#THIS dangerous! However, since we don't know year of measurement, we could only assume it as average
kevin_LMA$Year[is.na(kevin_LMA$Year)==TRUE] <- 2

kevin_soil <- read.csv("/Users/yunpeng/data/gcme/kevin/orig_leaf/pleco_soil.csv")


sitename <- kevin_othervars%>% group_by(site)  %>% summarise(number = n())
varname <- kevin_othervars%>% group_by(response)  %>% summarise(number = n())
varname

kevin_othervars%>% group_by(treatment)  %>% summarise(number = n())

#leaf_n needs division - but leaf_p doesn't need (as it is all pmass)
nmass <- subset(kevin_othervars,response=="leaf_n") %>%         
    filter(Unit %in% c("g","g_100g","g_g","g_kg","mg_g","mg_kg","mmol_g","ug_mg","umol_g"))
nmass$response <- "nmass"
logr_c_nmass <- as_tibble(response_ratio(subset(nmass,treatment=="c")));logr_w_nmass <- as_tibble(response_ratio(subset(nmass,treatment=="w")));logr_f_nmass <- as_tibble(response_ratio(subset(nmass,treatment=="f")));
logr_cw_nmass <- as_tibble(response_ratio(subset(nmass,treatment=="cw")));logr_cf_nmass <- as_tibble(response_ratio(subset(nmass,treatment=="cf")))

narea <- subset(kevin_othervars,response=="leaf_n") %>%         
  filter(Unit %in% c("g_m2","g_m3","g_pot","mg_cm2","mg_m2","mmol_m2","ug_cm2","ug_cm3","umol_m2"))
narea$response <- "narea"

logr_c_narea <- as_tibble(response_ratio(subset(narea,treatment=="c")));logr_w_narea <- as_tibble(response_ratio(subset(narea,treatment=="w")));logr_f_narea <- as_tibble(response_ratio(subset(narea,treatment=="f")));
logr_cw_narea <- as_tibble(response_ratio(subset(narea,treatment=="cw")));logr_cf_narea <- as_tibble(response_ratio(subset(narea,treatment=="cf")))


#other vars - all created now
for (i in 1:nrow(varname)) {
  tryCatch({
    varname1 <- varname$response[i]
    df_c <- subset(kevin_othervars,treatment=="c" & response==varname1)
    assign(paste("logr_c_", varname1,sep=""), as_tibble(response_ratio(df_c)))
    
    df_f <- subset(kevin_othervars,treatment=="f"& response==varname1)
    assign(paste("logr_f_", varname1,sep=""), as_tibble(response_ratio(df_f)))
    
    df_w <- subset(kevin_othervars,treatment=="w"& response==varname1)
    assign(paste("logr_w_", varname1,sep=""), as_tibble(response_ratio(df_w)))
    
    df_cf <- subset(kevin_othervars,treatment=="cf" & response==varname1)
    assign(paste("logr_cf_", varname1,sep=""), as_tibble(response_ratio(df_cf)))
    
    df_cw <- subset(kevin_othervars,treatment=="cw" &response==varname1)
    assign(paste("logr_cw_", varname1,sep=""), as_tibble(response_ratio(df_cw)))
    
    df_fw <- subset(kevin_othervars,treatment=="fw" & response==varname1)
    assign(paste("logr_fw_", varname1,sep=""), as_tibble(response_ratio(df_fw)))
    
    df_cfw <- subset(kevin_othervars,treatment=="cfw" & response==varname1)
    assign(paste("logr_cfw_", varname1,sep=""), as_tibble(response_ratio(df_cfw)))
  }, error=function(e){})} 

#newly added
logr_c_anpp <-  logr_c_anpp %>% filter(Unit %in% c("g_m2d","g_m2gs","g_m2y","gc_m2y","kg_ha_y","kg_m2y","kgc_m2y","mg_ha_y","t_ha_y"))
logr_c_anpp%>% group_by(site)  %>% summarise(number = n())

#lma, sla
sla <- subset(kevin_LMA,response=="sla")
lma <- subset(kevin_LMA,response=="lma")
#convert SLA to LMA format 
sla$ambient <- 1/sla$ambient; sla$elevated <- 1/sla$elevated; sla$ambient_Sd <- 1/sla$ambient_Sd
sla$elevated_Sd <- 1/sla$elevated_Sd;sla$ambient_Se <- 1/sla$ambient_Se;sla$elevated_Se <- 1/sla$elevated_Se
LMA <- dplyr::bind_rows(lma,sla)
LMA$response <- "LMA"
logr_c_LMA <- as_tibble(response_ratio(subset(LMA,treatment=="c")));logr_w_LMA <-  as_tibble(response_ratio(subset(LMA,treatment=="w")));logr_f_LMA <-  as_tibble(response_ratio(subset(LMA,treatment=="f")));
logr_cw_LMA <-  as_tibble(response_ratio(subset(LMA,treatment=="cw")));logr_cf_LMA <-  as_tibble(response_ratio(subset(LMA,treatment=="cf")))


#soil plaeo
kevin_soil2 <- kevin_soil[,c("exp","pre_soil_n_.g.g.","pre_soil_c.n","pre_total_p_.g.kg.")]
names(kevin_soil2) <- c("exp","soilN","soilCN","soilP")
kevin_soil2$soilCN[kevin_soil2$soilCN=="35.7;_26.8"] <- 35.7
kevin_soil2$soilCN[kevin_soil2$soilCN=="22.2;_35"] <- 22.2
kevin_soil2$soilCN <- as.numeric(kevin_soil2$soilCN)

#final, vcmax, jmax and jmax/vcmax
kevin <- read.csv("/Users/yunpeng/data/gcme/kevin/orig_vcmax/JunkePeng_11252021.csv")
kevin$site[is.na(kevin$site)==TRUE] <- "riceface_japan_a_2003_3938_14057"
kevin$exp_nam <- kevin$site #new plot: capo_caccia,itatinga,itatinga_experimental_station,setres
kevin[grep("riceface_japan", kevin$exp_nam),]$exp_nam <- "riceface_japan"

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
kevin_z$year_start <- kevin_z$start_year
kevin_z$year_end <- kevin_z$sampling_year

logr_c_vcmax <- as_tibble(response_ratio(subset(kevin_z,treatment=="c"&response=="vcmax")));logr_w_vcmax <-  as_tibble(response_ratio(subset(kevin_z,treatment=="w"&response=="vcmax")));logr_f_vcmax <-  as_tibble(response_ratio(subset(kevin_z,treatment=="f"&response=="vcmax")));
logr_cw_vcmax <-  as_tibble(response_ratio(subset(kevin_z,treatment=="cw"&response=="vcmax")));logr_cf_vcmax <-  as_tibble(response_ratio(subset(kevin_z,treatment=="cf"&response=="vcmax")));

logr_c_jmax <- as_tibble(response_ratio(subset(kevin_z,treatment=="c"&response=="jmax")));logr_w_jmax <-  as_tibble(response_ratio(subset(kevin_z,treatment=="w"&response=="jmax")));logr_f_jmax <-  as_tibble(response_ratio(subset(kevin_z,treatment=="f"&response=="jmax")));
logr_cw_jmax <-  as_tibble(response_ratio(subset(kevin_z,treatment=="cw"&response=="jmax")));logr_cf_jmax <-  as_tibble(response_ratio(subset(kevin_z,treatment=="cf"&response=="jmax")))

#for sensitivity coef
agg_meta_sen_coef <- function(df){
  threshold <- 0.4
  #set 0.4 rather than 0.5 newly - to decide whether calculating using rma or normal
  explist <- unique(df$exp)
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
    if (sum(is.na(subset(df,exp==explist[i])$logr_var),na.rm=TRUE)/length(is.na(subset(df,exp==explist[i])$logr_var)) >=threshold){
      #uncertainty <- 1.96*sd(as.vector(subset(df,exp==explist[i])$logr))
      #uncertainty_50 <- 0.67449*sd(as.vector(subset(df,exp==explist[i])$logr))
      uncertainty <- 0 # calculating sd directly from plots is horrible..
      uncertainty_50 <- 0
      mean_value <- mean(subset(df,exp==explist[i])$logr,na.rm=TRUE)
      n_samples <- length(subset(df,exp==explist[i])$logr_var)
      df_box <- tibble(
        exp=explist[i], middle = mean_value,
        ymin   = mean_value-uncertainty,
        ymax   = mean_value+uncertainty,
        ymin_quarter = mean_value-uncertainty_50,
        ymax_quarter   = mean_value+uncertainty_50,
        middle_scaled = NA, ymin_scaled   = NA, ymax_scaled   = NA,variance_info ="No_allNA",no=n_samples)
    } else {# define if uncertainty is 0, when N = 1
      if (length(subset(df,exp==explist[i])$logr_var) ==1){
        uncertainty <- 0
        mean_value <- mean(subset(df,exp==explist[i])$logr,na.rm=TRUE)
        n_samples <- length(subset(df,exp==explist[i])$logr_var)
        df_box <- tibble(
          exp=explist[i], middle = mean_value,
          ymin   = mean_value-uncertainty,
          ymax   = mean_value+uncertainty,
          ymin_quarter = mean_value-uncertainty,
          ymax_quarter   = mean_value+uncertainty,
          middle_scaled = NA, ymin_scaled   = NA, ymax_scaled   = NA,variance_info ="No_1_number",no=n_samples)
      } else {
        #where year is random factor
        out_meta <- metafor::rma.mv( logr, logr_var, method = "REML",
                                     random = ~ 1 | Year, slab = Year, control = list(stepadj=0.3),data = subset(df,exp==explist[i]))
        
        out_meta_quarter <- metafor::rma.mv( logr, logr_var, method = "REML",level=50,
                                             random = ~ 1 | Year, slab = Year, control = list(stepadj=0.3),data = subset(df,exp==explist[i]))
        
        out_meta_scaled <- predict( out_meta, transf=exp )
        n_samples <- length(subset(df,exp==explist[i])$logr_var)
        
        df_box <- tibble(
          exp=explist[i], middle = out_meta$b[1,1], ymin   = out_meta$ci.lb, ymax   = out_meta$ci.ub,
          ymin_quarter=out_meta_quarter$ci.lb, ymax_quarter = out_meta_quarter$ci.ub,
          middle_scaled = out_meta_scaled$pred, ymin_scaled   = out_meta_scaled$ci.lb, ymax_scaled   = out_meta_scaled$ci.ub,variance_info ="Yes",no=n_samples)}}
    mylist[[i]] <- df_box}
  output <- do.call("rbind",mylist)
  return(output)
}
#for response ratio
agg_meta <- function(df){
  threshold <- 0.4
  explist <- unique(df$exp)
  mylist <- list() #create an empty list
  #some are just all NA for original data's standard deviation - needs division -> in this way just calculate normal mean, lower and upper band
  for (i in 1:length(explist)){ # if over half of samples doens't have variance then we should not rely on this anyways.
    if (sum(is.na(subset(df,exp==explist[i])$logr_var),na.rm=TRUE)/length(is.na(subset(df,exp==explist[i])$logr_var)) >=threshold){
      #uncertainty <- 1.96*sd(as.vector(subset(df,exp==explist[i])$logr))
      #uncertainty_50 <- 0.67449*sd(as.vector(subset(df,exp==explist[i])$logr))
      uncertainty <- 0 # calculating sd directly from plots is horrible..
      uncertainty_50 <- 0
      mean_value <- mean(subset(df,exp==explist[i])$logr,na.rm=TRUE)
      n_samples <- length(subset(df,exp==explist[i])$logr_var)
      df_box <- tibble(
        exp=explist[i], middle = mean_value,
        ymin   = mean_value-uncertainty,
        ymax   = mean_value+uncertainty,
        ymin_quarter = mean_value-uncertainty_50,
        ymax_quarter   = mean_value+uncertainty_50,
        middle_scaled = NA, ymin_scaled   = NA, ymax_scaled   = NA,variance_info ="No_mostNA",no=n_samples)
    } else {# define if uncertainty is 0, when N = 1 with sd
      if ( length(subset(df,exp==explist[i])$logr_var) ==1){
        uncertainty <- 0
        mean_value <- mean(subset(df,exp==explist[i])$logr,na.rm=TRUE)
        n_samples <- length(subset(df,exp==explist[i])$logr_var)
        df_box <- tibble(
          exp=explist[i], middle = mean_value,
          ymin   = mean_value-uncertainty,
          ymax   = mean_value+uncertainty,
          ymin_quarter = mean_value-uncertainty,
          ymax_quarter   = mean_value+uncertainty,
          middle_scaled = NA, ymin_scaled   = NA, ymax_scaled   = NA,variance_info ="No_1_number",no=n_samples)
      } else {
        #where year is random factor
        out_meta <- metafor::rma.mv( logr, logr_var, method = "REML",
                                     random = ~ 1 | Year, slab = Year, control = list(stepadj=0.3),data = subset(df,exp==explist[i]))
        
        out_meta_quarter <- metafor::rma.mv( logr, logr_var, method = "REML",level=50,
                                             random = ~ 1 | Year, slab = Year, control = list(stepadj=0.3),data = subset(df,exp==explist[i]))
        
        out_meta_scaled <- predict( out_meta, transf=exp )
        n_samples <- length(subset(df,exp==explist[i])$logr_var)
        
        df_box <- tibble(
          exp=explist[i], middle = out_meta$b[1,1], ymin   = out_meta$ci.lb, ymax   = out_meta$ci.ub,
          ymin_quarter=out_meta_quarter$ci.lb, ymax_quarter = out_meta_quarter$ci.ub,
          middle_scaled = out_meta_scaled$pred, ymin_scaled   = out_meta_scaled$ci.lb, ymax_scaled   = out_meta_scaled$ci.ub,variance_info ="Yes",no=n_samples)}}
    mylist[[i]] <- df_box}
  output <- do.call("rbind",mylist)
  return(output)
}
#a new function to shows 1 trait vs. 9 other traits
merge_multi <- function(var_y,var_x1,var_x2,var_x3,var_x4,var_x5,var_x6,var_x7,var_x8,var_x9){
  plot_mean_y <- agg_meta_sen_coef(var_y)[,c("exp","middle")];names(plot_mean_y) <- c("exp",var_y$response[1])
  plot_mean_x1 <- agg_meta_sen_coef(var_x1)[,c("exp","middle")];names(plot_mean_x1) <- c("exp",var_x1$response[1])
  plot_mean_x2 <- agg_meta_sen_coef(var_x2)[,c("exp","middle")];names(plot_mean_x2) <- c("exp",var_x2$response[1])
  plot_mean_x3 <- agg_meta_sen_coef(var_x3)[,c("exp","middle")];names(plot_mean_x3) <- c("exp",var_x3$response[1])
  plot_mean_x4 <- agg_meta_sen_coef(var_x4)[,c("exp","middle")];names(plot_mean_x4) <- c("exp",var_x4$response[1])
  plot_mean_x5 <- agg_meta_sen_coef(var_x5)[,c("exp","middle")];names(plot_mean_x5) <- c("exp",var_x5$response[1])
  plot_mean_x6 <- agg_meta_sen_coef(var_x6)[,c("exp","middle")];names(plot_mean_x6) <- c("exp",var_x6$response[1])
  plot_mean_x7 <- agg_meta_sen_coef(var_x7)[,c("exp","middle")];names(plot_mean_x7) <- c("exp",var_x7$response[1])
  plot_mean_x8 <- agg_meta_sen_coef(var_x8)[,c("exp","middle")];names(plot_mean_x8) <- c("exp",var_x8$response[1])
  plot_mean_x9 <- agg_meta_sen_coef(var_x9)[,c("exp","middle")];names(plot_mean_x9) <- c("exp",var_x9$response[1])
  photo_final_others <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),
                              list(plot_mean_y,plot_mean_x1,plot_mean_x2,plot_mean_x3,plot_mean_x4,plot_mean_x5,
                                   plot_mean_x6,plot_mean_x7,plot_mean_x8,plot_mean_x9))
  return(photo_final_others)
}

#merge_multi for sensitivity coefficient
merge_multi_logr <- function(var_y,var_x1,var_x2,var_x3,var_x4,var_x5,var_x6,var_x7,var_x8,var_x9){
  plot_mean_y <- agg_meta(var_y)[,c("exp","middle")];names(plot_mean_y) <- c("exp",var_y$response[1])
  plot_mean_x1 <- agg_meta(var_x1)[,c("exp","middle")];names(plot_mean_x1) <- c("exp",var_x1$response[1])
  plot_mean_x2 <- agg_meta(var_x2)[,c("exp","middle")];names(plot_mean_x2) <- c("exp",var_x2$response[1])
  plot_mean_x3 <- agg_meta(var_x3)[,c("exp","middle")];names(plot_mean_x3) <- c("exp",var_x3$response[1])
  plot_mean_x4 <- agg_meta(var_x4)[,c("exp","middle")];names(plot_mean_x4) <- c("exp",var_x4$response[1])
  plot_mean_x5 <- agg_meta(var_x5)[,c("exp","middle")];names(plot_mean_x5) <- c("exp",var_x5$response[1])
  plot_mean_x6 <- agg_meta(var_x6)[,c("exp","middle")];names(plot_mean_x6) <- c("exp",var_x6$response[1])
  plot_mean_x7 <- agg_meta(var_x7)[,c("exp","middle")];names(plot_mean_x7) <- c("exp",var_x7$response[1])
  plot_mean_x8 <- agg_meta(var_x8)[,c("exp","middle")];names(plot_mean_x8) <- c("exp",var_x8$response[1])
  plot_mean_x9 <- agg_meta(var_x9)[,c("exp","middle")];names(plot_mean_x9) <- c("exp",var_x9$response[1])
  photo_final_others <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),
                              list(plot_mean_y,plot_mean_x1,plot_mean_x2,plot_mean_x3,plot_mean_x4,plot_mean_x5,
                                   plot_mean_x6,plot_mean_x7,plot_mean_x8,plot_mean_x9))
  return(photo_final_others)
}


#a function to plot output from merge_multi
plot_multi <- function(obj){
  p <- list()
  for(i in 1:9){
    p[[i]] <- ggplot(obj,aes_string(x=names(obj)[i+2], y=names(obj)[2])) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+
      geom_point()+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
      geom_smooth(color="red",method="lm",se=F)+theme_classic()}
  pp <- plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],nrow=3,label_size = 15)
  return(pp)
}

#how about combing soil N and soil total N?
logr_c_soil_n_combined <-  rbind(logr_c_soil_n,logr_c_soil_total_n)
logr_c_soil_n_combined$response <- "soil_N"

vcmax_all <- merge_multi(logr_c_vcmax,logr_c_jmax,logr_c_narea,logr_c_nmass,logr_c_leaf_cn,logr_c_LMA,logr_c_lai,logr_c_anpp,logr_c_soil_cn,logr_c_soil_n_combined) #logr_c_soil_n? or logr_c_soil_n_combined? 
plot_multi(vcmax_all)
ggsave(paste("~/data/output_gcme/colin/multi_new_vcmax.jpg",sep=""),width = 15, height = 15)

jmax_all <- merge_multi(logr_c_jmax,logr_c_vcmax,logr_c_narea,logr_c_nmass,logr_c_leaf_cn,logr_c_LMA,logr_c_lai,logr_c_anpp,logr_c_soil_cn,logr_c_soil_n_combined)
plot_multi(jmax_all)
ggsave(paste("~/data/output_gcme/colin/multi_new_jmax.jpg",sep=""),width = 15, height = 15)

narea_all <- merge_multi(logr_c_narea,logr_c_vcmax,logr_c_jmax,logr_c_nmass,logr_c_leaf_cn,logr_c_LMA,logr_c_lai,logr_c_anpp,logr_c_soil_cn,logr_c_soil_n_combined)
plot_multi(narea_all)
ggsave(paste("~/data/output_gcme/colin/multi_new_narea.jpg",sep=""),width = 15, height = 15)

leafCN_all <- merge_multi(logr_c_leaf_cn,logr_c_narea,logr_c_vcmax,logr_c_jmax,logr_c_nmass,logr_c_LMA,logr_c_lai,logr_c_anpp,logr_c_soil_cn,logr_c_soil_n_combined)
plot_multi(leafCN_all)
ggsave(paste("~/data/output_gcme/colin/multi_new_leafCN.jpg",sep=""),width = 15, height = 15)

nmass_all <- merge_multi(logr_c_nmass,logr_c_vcmax,logr_c_jmax,logr_c_narea,logr_c_leaf_cn,logr_c_LMA,logr_c_lai,logr_c_anpp,logr_c_soil_cn,logr_c_soil_n_combined)
plot_multi(nmass_all)
ggsave(paste("~/data/output_gcme/colin/multi_new_nmass.jpg",sep=""),width = 15, height = 15)

#soil time? - all FITTED WORSE
#plot_multi(merge_multi(logr_c_vcmax,logr_c_nmass,logr_c_soil_cn_min_layer,logr_c_soil_cn_org_layer,logr_c_soil_total_p,logr_c_soil_cn,logr_c_soil_n_org_layer,logr_c_soil_n_combined,logr_c_soil_total_n,logr_c_soil_n))
#plot_multi(merge_multi(logr_c_jmax,logr_c_vcmax,logr_c_soil_cn_min_layer,logr_c_soil_cn_org_layer,logr_c_soil_total_p,logr_c_soil_cn,logr_c_soil_n_org_layer,logr_c_soil_n_combined,logr_c_soil_total_n,logr_c_soil_n))
#plot_multi(merge_multi(logr_c_narea,logr_c_vcmax,logr_c_soil_cn_min_layer,logr_c_soil_cn_org_layer,logr_c_soil_total_p,logr_c_soil_cn,logr_c_soil_n_org_layer,logr_c_soil_n_combined,logr_c_soil_total_n,logr_c_soil_n))
#plot_multi(merge_multi(logr_c_leaf_cn,logr_c_vcmax,logr_c_soil_cn_min_layer,logr_c_soil_cn_org_layer,logr_c_soil_total_p,logr_c_soil_cn,logr_c_soil_n_org_layer,logr_c_soil_n_combined,logr_c_soil_total_n,logr_c_soil_n))
#plot_multi(merge_multi(logr_c_nmass,logr_c_vcmax,logr_c_soil_cn_min_layer,logr_c_soil_cn_org_layer,logr_c_soil_total_p,logr_c_soil_cn,logr_c_soil_n_org_layer,logr_c_soil_n_combined,logr_c_soil_total_n,logr_c_soil_n))

aa <- merge_multi(logr_c_vcmax,logr_c_nmass,logr_c_soil_cn_min_layer,logr_c_soil_cn_org_layer,logr_c_soil_total_p,logr_c_soil_cn,logr_c_soil_n_org_layer,logr_c_soil_n_combined,logr_c_soil_total_n,logr_c_soil_n)

#anwser Colin's question:
summary(lm(vcmax~jmax,vcmax_all))
summary(lm(leaf_cn~nmass,leafCN_all))

#add a function for ecm/am types
ecm_types <- function(df){
  df$ecm_type <- NA
  df$ecm_type[df$exp=="biocon_c"]<- "AM";df$ecm_type[df$exp=="giface_c"]<- "AM";df$ecm_type[df$exp=="nevada_desert_face_c"]<- "AM";df$ecm_type[df$exp=="new_zealand_face_c"]<- "AM";
  df$ecm_type[df$exp=="soyfacesoy2_c"]<- "AM";df$ecm_type[df$exp=="swissface_trifolium2_c"]<- "AM";
  df$ecm_type[df$exp=="duke_c"]<- "ECM";df$ecm_type[df$exp=="euroface4_pa_c"]<- "ECM";df$ecm_type[df$exp=="euroface4_pe_c"]<- "ECM";df$ecm_type[df$exp=="euroface4_pn_c"]<- "ECM";
  df$ecm_type[df$exp=="popface_pa_c"]<- "ECM";df$ecm_type[df$exp=="popface_pe_c"]<- "ECM";df$ecm_type[df$exp=="popface_pn_c"]<- "ECM";
  return(df)
}

#or using csv from /Users/yunpeng/data/gcme/kevin/orig_vcmax/new_ecm_types.csv
#as derived from gcme_vcmax.R

vcmax_all_ecm <- ecm_types(vcmax_all)
jmax_all_ecm <- ecm_types(jmax_all)

ggplot(vcmax_all_ecm,aes(x=anpp, y=vcmax)) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=ecm_type),size=5)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="red",method="lm",se=F)+geom_text(aes(label=exp),hjust=1, vjust=0,check_overlap = T)+theme_classic()
summary(lm(vcmax~anpp,vcmax_all_ecm))

ggplot(vcmax_all_ecm,aes(x=soil_N, y=vcmax)) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=ecm_type),size=5)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="red",method="lm",se=F)+geom_text(aes(label=exp),hjust=1, vjust=0,check_overlap = T)+theme_classic()


ggplot(jmax_all_ecm,aes(x=anpp, y=jmax)) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=ecm_type),size=5)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="red",method="lm",se=F)+geom_text(aes(label=exp),hjust=0, vjust=0,check_overlap = T)+theme_classic()
summary(lm(jmax~anpp,jmax_all_ecm))

#merge with pfts?
vcmax_ecosystem <- as.data.frame(logr_c_vcmax %>% group_by(exp,ecosystem,dominant_species) %>% summarise(number=n()))
vcmax_all_ecm_ecosystem <- merge(vcmax_all_ecm,vcmax_ecosystem,by=c("exp"),all.x=TRUE)
ggplot(vcmax_all_ecm_ecosystem,aes(x=anpp, y=vcmax)) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=ecosystem),size=5)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="red",method="lm",se=F)+geom_text(aes(label=exp),hjust=0, vjust=0,check_overlap = T)+theme_classic()


#combine for final response for multi variables
#agg_meta_sen_coef(logr_c_vcmax)

#analysis soil plaeo 
kevin_soil2 <- kevin_soil[,c("exp","pre_soil_n_.g.g.","pre_soil_c.n","pre_total_p_.g.kg.")]
names(kevin_soil2) <- c("exp","soilN","soilCN","soilP")
kevin_soil2$soilCN[kevin_soil2$soilCN=="35.7;_26.8"] <- 35.7;kevin_soil2$soilCN[kevin_soil2$soilCN=="22.2;_35"] <- 22.2
kevin_soil2$soilCN <- as.numeric(kevin_soil2$soilCN)

kevin_soil2$soilN_number <- as.numeric(kevin_soil2$soilN)
as.data.frame(subset(kevin_soil2,is.na(soilN_number)==TRUE & is.na(soilN)==FALSE)%>% group_by(soilN)  %>% summarise(number = n()))$soilN #check all this and convert manually
#default unit is g/g - correct some numbers
kevin_soil2$soilN_number[kevin_soil2$soilN=="0,0009-0,0016"] <- (0.0009+0.0016)/2
kevin_soil2$soilN_number[kevin_soil2$soilN=="0,28_g_kg-1"] <- 0.28/1000
kevin_soil2$soilN_number[kevin_soil2$soilN=="0,3_g_kg-1"] <- 0.30/1000
kevin_soil2$soilN_number[kevin_soil2$soilN=="0,36_g_kg-1"] <- 0.36/1000
kevin_soil2$soilN_number[kevin_soil2$soilN=="0,5_mg_kg-1"] <- 0.5/1000000
kevin_soil2$soilN_number[kevin_soil2$soilN=="0.0065;_0.0004"] <- (0.0065+0.0004)/2
kevin_soil2$soilN_number[kevin_soil2$soilN=="0.0134;_0.0092"] <- (0.0134+0.0092)/2
kevin_soil2$soilN_number[kevin_soil2$soilN=="1,1_g_kg-1"] <- 1.1/1000
kevin_soil2$soilN_number[kevin_soil2$soilN=="1,25_g_kg-1"] <- 1.25/1000
kevin_soil2$soilN_number[kevin_soil2$soilN=="1,26_gn_kg-1"] <- 1.26/1000
kevin_soil2$soilN_number[kevin_soil2$soilN=="1,35_(g_n_m-2_y-1)"] <- NA
kevin_soil2$soilN_number[kevin_soil2$soilN=="1,59_g_kg-1"] <- 1.59/1000
kevin_soil2$soilN_number[kevin_soil2$soilN=="1,74_g_mg-1"] <- 1.74/1000
kevin_soil2$soilN_number[kevin_soil2$soilN=="12,8_mg_g-1"] <- 12.8/1000
kevin_soil2$soilN_number[kevin_soil2$soilN=="13,3_mg_g-1"] <- 13.3/1000
kevin_soil2$soilN_number[kevin_soil2$soilN=="13,3_mg_g-2"] <- 13.3/1000
kevin_soil2$soilN_number[kevin_soil2$soilN=="13,3_mg_g-3"] <- 13.3/1000
kevin_soil2$soilN_number[kevin_soil2$soilN=="2,63_(dw)"] <- NA
kevin_soil2$soilN_number[kevin_soil2$soilN=="21_kg_n_ha-1"] <- NA
kevin_soil2$soilN_number[kevin_soil2$soilN=="4,31_(t_ha-1)"] <- NA
kevin_soil2$soilN_number[kevin_soil2$soilN=="4,5_g_kg-1"] <- 4.5/1000
kevin_soil2$soilN_number[kevin_soil2$soilN=="45_kg_n/_ha_"] <- NA
kevin_soil2$soilN_number[kevin_soil2$soilN=="6,4_mg_kg-1"] <- 6.4/1000000
kevin_soil2$soilN_number[kevin_soil2$soilN=="600__mol_g-12_mg_g-1"] <- NA
kevin_soil2$soilN_number[kevin_soil2$soilN=="7,2_mg_kg-1"] <- 7.2/1000000
kevin_soil2$soilN_number[kevin_soil2$soilN=="7,39_(t_ha-1)"] <- NA
kevin_soil2$soilN_number[kevin_soil2$soilN=="9,6_mg_kg-1"] <- 9.6/1000000
as.data.frame(subset(kevin_soil2,is.na(soilN_number)==TRUE & is.na(soilN)==FALSE)%>% group_by(soilN)  %>% summarise(number = n()))$soilN #check all this and convert manually
summary(kevin_soil2$soilN_number)

soilN_final <- aggregate(kevin_soil2,by=list(kevin_soil2$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","soilN_number")]
names(soilN_final) <- c("exp","soilN")

vcmax_all_ecm_soilN <- merge(vcmax_all_ecm,soilN_final,by=c("exp"),all.x=TRUE)
#after merging...not many soil records from paleo - unfortunately - should still use logr_soilN
vcmax_all_ecm_soilN$soil_N[vcmax_all_ecm_soilN$soil_N==0] <-NA
ggplot(vcmax_all_ecm_soilN,aes(x=anpp, y=vcmax)) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(shape=ecm_type,color=(soil_N)),size=5)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="red",method="lm",se=F)+geom_text(aes(label=exp),hjust=1, vjust=0,check_overlap = T)+theme_classic()+scale_color_gradientn(colours = rainbow(5))
summary(lm(vcmax~anpp,vcmax_all_ecm))

ggplot(vcmax_all_ecm_soilN,aes(x=nmass, y=vcmax)) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=(soil_N)),size=5)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="red",method="lm",se=F)+geom_text(aes(label=exp),hjust=1, vjust=0,check_overlap = T)+theme_classic()+scale_color_gradientn(colours = rainbow(5))
summary(lm(vcmax~anpp,vcmax_all_ecm))

#now, newly adding bnpp and Asat
df <- read_csv("~/data/gcme/data_received_190325/NewData_wide_CORRECTED2.csv") %>%
  mutate( ambient_Sd  = as.numeric(ambient_Sd),  ambient_Se  = as.numeric(ambient_Se), 
          elevated_Sd = as.numeric(elevated_Sd), elevated_Se = as.numeric(elevated_Se) )

#for anpp, bnpp, Asat, c13
df$response <- NA
#bnpp
df$response[df$Data_type=="BNPP"] <- "BNPP";df$response[df$Data_type=="dead_fine_root_biomass"] <- "BNPP";df$response[df$Data_type=="DOC_exudation"] <- "BNPP";
df$response[df$Data_type=="exudation"] <- "BNPP";df$response[df$Data_type=="fine_root_biomass"] <- "BNPP";df$response[df$Data_type=="fine_root_density"] <- "BNPP";
df$response[df$Data_type=="fine_root_production"] <- "BNPP";df$response[df$Data_type=="fine_root_turnover"] <- "BNPP"
#anpp
df$response[df$Data_type=="ANPP"] <- "ANPP"#;df$response[df$Data_type=="aboveground_biomass"] <- "ANPP";df$response[df$Data_type=="aboveground_biomass_"] <- "ANPP"
#Asat
df$response[df$Data_type=="Asat"] <- "Asat";df$response[df$Data_type=="Amax"] <- "Amax";df$response[df$Data_type=="Anet"] <- "Anet"
#c13
df$response[df$Data_type=="d13C_isotope_ratio_hardwood_leaflitter"] <- "c13";df$response[df$Data_type=="d13C_isotope_ratio_leaflitter"] <- "c13";
df$response[df$Data_type=="leaf_d13C"] <- "c13";df$response[df$Data_type=="leaf_litter_d13C"] <- "c13";df$response[df$Data_type=="litter_stable_13C_isotope_signature"] <- "c13"
#N uptake
df$response[df$Data_type=="N_uptake"] <- "Nuptake";df$response[df$Data_type=="NH4+_uptake"] <- "Nuptake";df$response[df$Data_type=="NO3-_uptake"] <- "Nuptake"
#litter mass
df$response[df$Data_type=="litter_prod"] <- "litter";df$response[df$Data_type=="litter_mass"] <- "litter"

df_only <- subset(df,is.na(response)==FALSE)

#1. examine data and fill missing measurement year
empty_sampling_exp <- unique(subset(df_only,is.na(Sampling_Year==TRUE))$exp_nam)
#convert the year from missing to median value
for(i in 1:length(empty_sampling_exp)){
  sampling_value <- as.integer(median(df_only$Sampling_Year[is.na(df_only$Sampling_Year)==FALSE & df_only$exp_nam==empty_sampling_exp[i]]))
  df_only$Sampling_Year[is.na(df_only$Sampling_Year)==TRUE & df_only$exp_nam==empty_sampling_exp[i]] <- sampling_value
}
#looking at data manually and correct
df_only$Sampling_Year[is.na(df_only$Sampling_Year)==TRUE & (df_only$Start_Year)==1989] <- 2009
df_only$Sampling_Year[is.na(df_only$Sampling_Year)==TRUE & (df_only$Start_Year)==2003] <- 2003
df_only$Start_Year[is.na(df_only$Start_Year)==TRUE & (df_only$Sampling_Year)==1998] <- 1998
df_only$Start_Year[is.na(df_only$Start_Year)==TRUE & (df_only$Sampling_Year)==2014] <- 2014
#Aberdeen_lolium was missing for both Sampling_Year and Start_Year - just left here and set year =1
df_only$Year <- df_only$Sampling_Year - df_only$Start_Year
df_only$Year[df_only$Year<0] <- 0
df_only$Year[is.na(df_only$Year)==TRUE] <- 0
summary(df_only$Year)

#2. convert to lower case
varname <- df_only%>% group_by(response)  %>% summarise(number = n())
df_only$exp <- tolower(df_only$prev_name)

#3. divide into dataframe
for (i in 1:nrow(varname)) {
  tryCatch({
    varname1 <- varname$response[i]
    df_c <- subset(df_only,treatment=="c" & response==varname1)
    assign(paste("old_logr_c_", varname1,sep=""), as_tibble(response_ratio(df_c)))
    
    df_f <- subset(df_only,treatment=="f"& response==varname1)
    assign(paste("old_logr_f_", varname1,sep=""), as_tibble(response_ratio(df_f)))
    
    df_w <- subset(df_only,treatment=="w"& response==varname1)
    assign(paste("old_logr_w_", varname1,sep=""), as_tibble(response_ratio(df_w)))
    
    df_cf <- subset(df_only,treatment=="cf" & response==varname1)
    assign(paste("old_logr_cf_", varname1,sep=""), as_tibble(response_ratio(df_cf)))
    
    df_cw <- subset(df_only,treatment=="cw" &response==varname1)
    assign(paste("old_logr_cw_", varname1,sep=""), as_tibble(response_ratio(df_cw)))
    
    df_fw <- subset(df_only,treatment=="fw" & response==varname1)
    assign(paste("old_logr_fw_", varname1,sep=""), as_tibble(response_ratio(df_fw)))
    
    df_cfw <- subset(df_only,treatment=="cfw" & response==varname1)
    assign(paste("old_logr_cfw_", varname1,sep=""), as_tibble(response_ratio(df_cfw)))
  }, error=function(e){})}
# never used old_logr_c_ANPP - since it adds more data from logr_c_anpp!

#now, for vcmax
vcmax_all2_c <- merge_multi(logr_c_vcmax,logr_c_jmax,logr_c_narea,logr_c_nmass,logr_c_leaf_cn,logr_c_soil_n, logr_c_LMA,logr_c_anpp,old_logr_c_BNPP,old_logr_c_Nuptake) 
vcmax_all2_c$condition <- "co2"
vcmax_all2_cf <- merge_multi(logr_cf_vcmax,logr_cf_jmax,logr_cf_narea,logr_cf_nmass,logr_cf_leaf_cn,logr_cf_soil_n,logr_cf_LMA,logr_cf_anpp,old_logr_cf_BNPP,old_logr_cf_Nuptake) 
vcmax_all2_cf$condition <- "co2+low_N"
#one plot from cw_soil_n has weired varice - deleting them now
vcmax_all2_cw <- merge_multi(logr_cw_vcmax,logr_cw_jmax,logr_cw_narea,logr_cw_nmass,logr_cw_leaf_cn,logr_cw_soil_n[c(1:3,8:10),], logr_cw_LMA,logr_cw_anpp,old_logr_cw_BNPP,old_logr_cw_Nuptake) 
vcmax_all2_cw$condition <- "co2+T"

#now, looking at additional plots - but just logr
vcmax_all2_f_logr <- merge_multi_logr(logr_f_vcmax,logr_f_jmax,logr_f_narea,logr_c_nmass,logr_f_leaf_cn,logr_f_soil_n, logr_f_LMA,logr_f_anpp,old_logr_f_BNPP,old_logr_f_Nuptake) 
vcmax_all2_f_logr <- subset(vcmax_all2_f_logr,exp!="itatinga_fk")
vcmax_all2_cf_logr <- merge_multi_logr(logr_cf_vcmax,logr_cf_jmax,logr_cf_narea,logr_cf_nmass,logr_cf_leaf_cn,logr_cf_soil_n,logr_cf_LMA,logr_cf_anpp,old_logr_cf_BNPP,old_logr_cf_Nuptake) 
# this two are equal and could be obtained completely! = cf -f 
vcmax_all2_f_logr[,2:ncol(vcmax_all2_f_logr)] <- vcmax_all2_cf_logr[,2:ncol(vcmax_all2_cf_logr)] - vcmax_all2_f_logr[,2:ncol(vcmax_all2_f_logr)]
#now calculating sensitivity coefficient 
vcmax_all2_f_logr$exp # duke2,euroface4_pa,pe,pn, new_zealand_face
#first - show them
subset(logr_cf_vcmax,exp=="duke2_cf")$co2_e-subset(logr_cf_vcmax,exp=="duke2_cf")$co2_a
subset(logr_cf_vcmax,exp=="euroface4_pa_cf")$co2_e-subset(logr_cf_vcmax,exp=="euroface4_pa_cf")$co2_a
subset(logr_cf_vcmax,exp=="euroface4_pe_cf")$co2_e-subset(logr_cf_vcmax,exp=="euroface4_pe_cf")$co2_a
subset(logr_cf_vcmax,exp=="euroface4_pn_cf")$co2_e-subset(logr_cf_vcmax,exp=="euroface4_pn_cf")$co2_a
subset(logr_cf_vcmax,exp=="new_zealand_face_cf")$co2_e-subset(logr_cf_vcmax,exp=="new_zealand_face_cf")$co2_a
subset(logr_cf_vcmax,exp=="setres_cf")$co2_e-subset(logr_cf_vcmax,exp=="setres_cf")$co2_a
aco2list <- c(200,182,182,182,111,NA)
#now, divide them to calculate final sensitivity coefficient!
vcmax_all2_f_logr[,2:ncol(vcmax_all2_f_logr)] <- vcmax_all2_f_logr[,2:ncol(vcmax_all2_f_logr)] / log(aco2list)
vcmax_all2_f_logr$condition <- "co2+high_N"

vcmax_all2_final <- dplyr::bind_rows(vcmax_all2_c,vcmax_all2_cf,vcmax_all2_cw,vcmax_all2_f_logr) 
vcmax_all2_final$Nup_BNPP <- vcmax_all2_final$Nuptake/vcmax_all2_final$BNPP
vcmax_all2_final$ANPP_BNPP <- vcmax_all2_final$anpp/vcmax_all2_final$BNPP
dim(vcmax_all2_final)

ecm_csv <- read.csv("/Users/yunpeng/data/gcme/kevin/orig_vcmax/new_ecm_types.csv")
ecm_csv2 <- ecm_csv[,c("exp","ecm_type","ecosystem")]
ecm_csv2$rep <- duplicated(ecm_csv2$exp)
ecm_csv2 <- subset(ecm_csv2,rep==FALSE)
ecm_csv2 <- ecm_csv2[,c("exp","ecm_type","ecosystem")]
vcmax_all2_final_ecm <- merge(vcmax_all2_final,ecm_csv2,by=c("exp"),all.x=TRUE)

vcmax_all2_final_ecm$ecosystem[vcmax_all2_final_ecm$ecosystem=="heathland"] <- "grassland"
vcmax_all2_final_ecm$ecosystem[vcmax_all2_final_ecm$ecosystem=="temperate_forest"] <- "forest"
vcmax_all2_final_ecm$ecosystem[vcmax_all2_final_ecm$ecosystem=="shrubland"] <- "forest"

#one outlier has problem? - must something wrong with magnitude
vcmax_all2_final_ecm$anpp[vcmax_all2_final_ecm$anpp< -5 ] <- NA

obj <- subset(vcmax_all2_final_ecm,condition!="co2+T")
#obj <- subset(vcmax_all2_final_ecm,condition=="co2")

#condition, ecm
p <- list()
for(i in c(1:9,11,12)){
  p[[i]] <- ggplot(obj,aes_string(x=names(obj)[i+2], y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=condition,shape=ecm_type),size=3)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+theme_classic()}
# geom_text(aes(label=exp),hjust=1, vjust=0,check_overlap = T)
plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[11]],p[[12]],nrow=4,label_size = 15)+ theme(plot.background=element_rect(fill="white", color="white"))
#ggsave(paste("~/data/output_gcme/colin/multi_new_vc_combination.jpg",sep=""),width = 15, height = 15)

#figure 2 - Jmax, Nmass, LMA, ANPP, BNPP, Nup/BNPP
plot_grid(p[[1]],p[[3]],p[[6]],p[[7]],p[[8]],p[[11]],nrow=2,label_size = 15)+ theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/final_fig2_vcmax.jpg",sep=""),width = 12, height = 6)


#condition, ecosystem
p <- list()
for(i in c(1:9,11,12)){
  p[[i]] <- ggplot(obj,aes_string(x=names(obj)[i+2], y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=ecosystem,shape=condition),size=3)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+theme_classic()}
plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[11]],p[[12]],nrow=4,label_size = 15)+ theme(plot.background=element_rect(fill="white", color="white"))

#for jmax
jmax_all2_c <- merge_multi(logr_c_jmax,logr_c_vcmax,logr_c_narea,logr_c_nmass,logr_c_leaf_cn,logr_c_soil_n, logr_c_LMA,logr_c_anpp,old_logr_c_BNPP,old_logr_c_Nuptake) 
jmax_all2_c$condition <- "co2"
jmax_all2_cf <- merge_multi(logr_cf_jmax,logr_cf_vcmax,logr_cf_narea,logr_cf_nmass,logr_cf_leaf_cn,logr_cf_soil_n,logr_cf_LMA,logr_cf_anpp,old_logr_cf_BNPP,old_logr_cf_Nuptake) 
jmax_all2_cf$condition <- "co2+low_N"
#one plot from cw_soil_n has weired varice - deleting them now
jmax_all2_cw <- merge_multi(logr_cw_jmax,logr_cw_vcmax,logr_cw_narea,logr_cw_nmass,logr_cw_leaf_cn,logr_cw_soil_n[c(1:3,8:10),], logr_cw_LMA,logr_cw_anpp,old_logr_cw_BNPP,old_logr_cw_Nuptake) 
jmax_all2_cw$condition <- "co2+T"

#now, looking at additional plots - but just logr
jmax_all2_f_logr <- merge_multi_logr(logr_f_jmax,logr_f_vcmax,logr_f_narea,logr_c_nmass,logr_f_leaf_cn,logr_f_soil_n, logr_f_LMA,logr_f_anpp,old_logr_f_BNPP,old_logr_f_Nuptake) 
jmax_all2_cf_logr <- merge_multi_logr(logr_cf_jmax,logr_cf_vcmax,logr_cf_narea,logr_cf_nmass,logr_cf_leaf_cn,logr_cf_soil_n,logr_cf_LMA,logr_cf_anpp,old_logr_cf_BNPP,old_logr_cf_Nuptake) 
# this two are equal and could be obtained completely! = cf -f 
jmax_all2_f_logr[,2:ncol(jmax_all2_f_logr)] <- jmax_all2_cf_logr[,2:ncol(jmax_all2_cf_logr)] - jmax_all2_f_logr[,2:ncol(jmax_all2_f_logr)]
#now calculating sensitivity coefficient 
jmax_all2_f_logr$exp # duke2,euroface4_pa,pe,pn, new_zealand_face
#first - show them
subset(logr_cf_jmax,exp=="duke2_cf")$co2_e-subset(logr_cf_jmax,exp=="duke2_cf")$co2_a
subset(logr_cf_jmax,exp=="euroface4_pa_cf")$co2_e-subset(logr_cf_jmax,exp=="euroface4_pa_cf")$co2_a
subset(logr_cf_jmax,exp=="euroface4_pe_cf")$co2_e-subset(logr_cf_jmax,exp=="euroface4_pe_cf")$co2_a
subset(logr_cf_jmax,exp=="euroface4_pn_cf")$co2_e-subset(logr_cf_jmax,exp=="euroface4_pn_cf")$co2_a
subset(logr_cf_jmax,exp=="setres_cf")$co2_e-subset(logr_cf_jmax,exp=="setres_cf")$co2_a
bco2list <- c(200,182,182,182,NA)
#now, divide them to calculate final sensitivity coefficient!
jmax_all2_f_logr[,2:ncol(jmax_all2_f_logr)] <- jmax_all2_f_logr[,2:ncol(jmax_all2_f_logr)] / log(bco2list)
jmax_all2_f_logr$condition <- "co2+high_N"

jmax_all2_final <- dplyr::bind_rows(jmax_all2_c,jmax_all2_cf,jmax_all2_cw,jmax_all2_f_logr) 
jmax_all2_final$Nup_BNPP <- jmax_all2_final$Nuptake/jmax_all2_final$BNPP
jmax_all2_final$ANPP_BNPP <- jmax_all2_final$anpp/jmax_all2_final$BNPP

jmax_all2_final_ecm <- merge(jmax_all2_final,ecm_csv2,by=c("exp"),all.x=TRUE)
jmax_all2_final_ecm$ecosystem[jmax_all2_final_ecm$ecosystem=="heathland"] <- "grassland"
jmax_all2_final_ecm$ecosystem[jmax_all2_final_ecm$ecosystem=="temperate_forest"] <- "forest"
jmax_all2_final_ecm$ecosystem[jmax_all2_final_ecm$ecosystem=="shrubland"] <- "forest"

#one outlier has problem? - must something wrong with magnitude
jmax_all2_final_ecm$anpp[jmax_all2_final_ecm$anpp< -5 ] <- NA

#obj <- jmax_all2_final_ecm
obj <- subset(jmax_all2_final_ecm,condition!="co2+T")

#condition, ecm
p <- list()
for(i in c(1:9,11,12)){
  p[[i]] <- ggplot(obj,aes_string(x=names(obj)[i+2], y="jmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=condition,shape=ecm_type),size=3)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+theme_classic()}
# geom_text(aes(label=exp),hjust=1, vjust=0,check_overlap = T)
#plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[11]],p[[12]],nrow=4,label_size = 15)+ theme(plot.background=element_rect(fill="white", color="white"))
#ggsave(paste("~/data/output_gcme/colin/multi_new_j_combination.jpg",sep=""),width = 15, height = 15)

#figure 2_jmax - Narea, Nmass, LMA, ANPP, BNPP, Nup/BNPP
plot_grid(p[[2]],p[[3]],p[[6]],p[[7]],p[[8]],p[[11]],nrow=2,label_size = 25)+ theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/final_fig2_jmax.jpg",sep=""),width = 12, height = 6)


#now, for jmax/vcmax
jmax_all2_final_ecm$jmax_vcmax <- jmax_all2_final_ecm$jmax/jmax_all2_final_ecm$vcmax
obj <- jmax_all2_final_ecm
#obj <- subset(jmax_all2_final_ecm,condition!="co2+T")

for(i in c(1:10,12,13)){
  p[[i]] <- ggplot(obj,aes_string(x=names(obj)[i+1], y="jmax_vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=condition,shape=ecm_type),size=3)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+theme_classic()}
# geom_text(aes(label=exp),hjust=1, vjust=0,check_overlap = T)
plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],p[[10]],p[[12]],p[[13]],nrow=4,label_size = 15)+ theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/multi_new_j_v_combination.jpg",sep=""),width = 15, height = 15)

#combine to one results
agg_meta_plots <- function(df,type_name,logr,log_var){
  df$logr <- df$middle
  df$logr_var <- ((df$middle-df$ymin)/1.96)^2
  df$logr_var[df$logr_var==0] <- NA
  
  explist <- unique(eval(parse(text=paste("df$", type_name, sep = ""))))
  mylist <- list() #create an empty list
  df$logr <- eval(parse(text=paste("df$", logr, sep = "")))
  df$logr_var <- eval(parse(text=paste("df$", log_var, sep = "")))
  for (i in 1:length(explist)){
    out_meta <- metafor::rma.mv( logr, logr_var, method = "REML", control = list(stepadj=0.3), 
                                 data= subset(df,eval(parse(text=paste(type_name)))==explist[i]))
    out_meta_quarter <- metafor::rma.mv( logr, logr_var, method = "REML",level=50,control = list(stepadj=0.3),
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

agg_meta_plots_ratio <- function(df,logr){
  final_data <- eval(parse(text=paste("df$", logr, sep = "")))
  df_box <- tibble(
      type_name=logr, middle = mean(final_data,na.rm=TRUE), ymin   = quantile(final_data,probs=c(.025,.975),na.rm=TRUE)[1], ymax   = quantile(final_data,probs=c(.025,.975),na.rm=TRUE)[2],
      ymin_quarter=quantile(final_data,probs=c(.25,.75),na.rm=TRUE)[1], ymax_quarter = quantile(final_data,probs=c(.25,.75),na.rm=TRUE)[2],
      no_plots = sum(is.na(jv_ratio$jmax_vcmax)==FALSE))
  return(df_box)
}

all_logr_c_vcmax <- agg_meta_sen_coef(logr_c_vcmax);all_logr_c_vcmax$response <- "vcmax"
all_logr_c_jmax <- agg_meta_sen_coef(logr_c_jmax);all_logr_c_jmax$response <- "jmax"
all_logr_c_LMA <- agg_meta_sen_coef(logr_c_LMA);all_logr_c_LMA$response <- "LMA"
all_logr_c_narea <- agg_meta_sen_coef(logr_c_narea);all_logr_c_narea$response <- "narea"
all_logr_c_nmass <- agg_meta_sen_coef(logr_c_nmass);all_logr_c_nmass$response <- "nmass"
all_logr_c_leaf_cn <- agg_meta_sen_coef(logr_c_leaf_cn);all_logr_c_leaf_cn$response <- "leafcn"
all_logr_c_anpp <- agg_meta_sen_coef(logr_c_anpp);all_logr_c_anpp$response <- "anpp"
all_logr_c_lai <- agg_meta_sen_coef(logr_c_lai);all_logr_c_lai$response <- "LAI"
all_old_logr_c_BNPP <- agg_meta_sen_coef(old_logr_c_BNPP);all_old_logr_c_BNPP$response <- "bnpp"
all_old_logr_c_Nuptake <- agg_meta_sen_coef(old_logr_c_Nuptake);all_old_logr_c_Nuptake$response <- "Nuptake"
all_old_logr_c_Asat <- agg_meta_sen_coef(old_logr_c_Asat);all_old_logr_c_Asat$response <- "Asat"
#all_old_logr_c_Amax <- agg_meta_sen_coef(old_logr_c_Amax);all_old_logr_c_Amax$response <- "Amax"

jv_ratio <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),
                            list(all_logr_c_vcmax[,c("exp","middle")],all_logr_c_jmax[,c("exp","middle")]))
jv_ratio$jmax_vcmax <- jv_ratio$middle.y/jv_ratio$middle.x
all_logr_c_jvratio <- agg_meta_plots_ratio(jv_ratio,"jmax_vcmax")

anpp_bnpp_ratio <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),
                  list(all_logr_c_anpp[,c("exp","middle")],all_old_logr_c_BNPP[,c("exp","middle")]))
anpp_bnpp_ratio$anpp_bnpp <- anpp_bnpp_ratio$middle.x/anpp_bnpp_ratio$middle.y
#one outlier - exclude
anpp_bnpp_ratio <- subset(anpp_bnpp_ratio,anpp_bnpp<10 & anpp_bnpp > -10)
all_logr_c_anpp_bnpp_ratio <- agg_meta_plots_ratio(anpp_bnpp_ratio,"anpp_bnpp")

nuptake_bnpp_ratio <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),
                         list(all_old_logr_c_Nuptake[,c("exp","middle")],all_old_logr_c_BNPP[,c("exp","middle")]))
nuptake_bnpp_ratio$nuptake_bnpp <- nuptake_bnpp_ratio$middle.x/nuptake_bnpp_ratio$middle.y
all_logr_c_nuptake_bnpp_ratio <- agg_meta_plots_ratio(nuptake_bnpp_ratio,"nuptake_bnpp")

all_variables <-  rbind(all_old_logr_c_Nuptake,all_old_logr_c_BNPP,all_logr_c_anpp,all_logr_c_lai,all_logr_c_leaf_cn,all_logr_c_LMA,all_logr_c_narea,all_logr_c_nmass,all_old_logr_c_Asat,all_logr_c_jmax,all_logr_c_vcmax)

overall_response <- agg_meta_plots(all_variables,"response","logr","logr_var")

overall_response_final <- dplyr::bind_rows(all_logr_c_jvratio,all_logr_c_anpp_bnpp_ratio,all_logr_c_nuptake_bnpp_ratio,overall_response)

overall_response_final$type_name <- factor(overall_response_final$type_name, levels = overall_response_final$type_name)

overall_response_final %>%
  ggplot( aes(x=type_name, y=middle)) + 
  geom_crossbar(aes(x=type_name, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5) +
  geom_crossbar(aes(x=type_name, y=middle, ymin=ymin_quarter, ymax=ymax_quarter),width = 0.5) +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-2,2)+
  labs(x="", y="CO2 - Sensitivity coefficient") + theme_classic()+coord_flip() +theme(axis.text=element_text(size=12))
ggsave(paste("~/data/output_gcme/colin/final_fig1.jpg",sep=""),width = 10, height = 5)

overall_response_final2 <- overall_response_final
overall_response_final2$ymin[overall_response_final2$type_name=="jmax_vcmax"] <- NA;overall_response_final2$ymax[overall_response_final2$type_name=="jmax_vcmax"] <- NA;overall_response_final2$ymin_quarter[overall_response_final2$type_name=="jmax_vcmax"] <- NA;overall_response_final2$ymax_quarter[overall_response_final2$type_name=="jmax_vcmax"] <- NA
overall_response_final2$ymin[overall_response_final2$type_name=="anpp_bnpp"] <- NA;overall_response_final2$ymax[overall_response_final2$type_name=="anpp_bnpp"] <- NA;overall_response_final2$ymin_quarter[overall_response_final2$type_name=="anpp_bnpp"] <- NA;overall_response_final2$ymax_quarter[overall_response_final2$type_name=="anpp_bnpp"] <- NA
overall_response_final2$ymin[overall_response_final2$type_name=="nuptake_bnpp"] <- NA;overall_response_final2$ymax[overall_response_final2$type_name=="nuptake_bnpp"] <- NA;overall_response_final2$ymin_quarter[overall_response_final2$type_name=="nuptake_bnpp"] <- NA;overall_response_final2$ymax_quarter[overall_response_final2$type_name=="nuptake_bnpp"] <- NA

overall_response_final2 %>%
  ggplot( aes(x=type_name, y=middle)) + 
  geom_crossbar(aes(x=type_name, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5) +
  geom_crossbar(aes(x=type_name, y=middle, ymin=ymin_quarter, ymax=ymax_quarter),width = 0.5) +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-0.6,0.6)+
  labs(x="", y="CO2 - Sensitivity coefficient") + theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))
ggsave(paste("~/data/output_gcme/colin/final_fig1_update.jpg",sep=""),width = 10, height = 5)


#now, combined with prediction data from gcme_vcmax
#the prediction of vcmax only includes c, w and cw - but having same value in cf
prediction <- read.csv("/Users/yunpeng/data/gcme/kevin/forcing/pred_vcmax.csv")
pred_vcmax <- subset(prediction,response=="vcmax")
vcmax_plotmean <- aggregate(pred_vcmax,by=list(pred_vcmax$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","pred_vcmax25_coef")]
names(vcmax_plotmean) <- c("exp","pred_vcmax"); vcmax_plotmean <- na.omit(vcmax_plotmean)

#
final_prediction <- Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),list(vcmax_all2_final_ecm,vcmax_plotmean))

a1 <- ggplot(subset(final_prediction,condition=="co2"),aes(x=vcmax, y=pred_vcmax)) + 
  geom_point(aes(color=ecm_type,shape=ecosystem),size=3)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_text(aes(label=exp),hjust=1, vjust=0,check_overlap = T)+labs(x="Observation vcmax", y="Prediction vcmax")+
  geom_smooth(color="black",method="lm",se=T)+theme_classic()

logr_c_vcmax2 <- Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),list(logr_c_vcmax,final_prediction[,c("exp","ecm_type","ecosystem","pred_vcmax")]))
#remove setres_cf as no prediction available
logr_c_vcmax2 <- subset(logr_c_vcmax2,exp!="setres_cf") # but already not included

vcmax_ecm <- agg_meta_sen_coef(subset(logr_c_vcmax2,ecm_type=="ECM")); vcmax_ecm$response <- "ECM"
vcmax_am <- agg_meta_sen_coef(subset(logr_c_vcmax2,ecm_type=="AM")); vcmax_am$response <- "AM"
vcmax_grassland <- agg_meta_sen_coef(subset(logr_c_vcmax2,ecosystem.y=="grassland")); vcmax_grassland$response <- "grassland"
vcmax_forest <- agg_meta_sen_coef(subset(logr_c_vcmax2,ecosystem.y=="forest")); vcmax_forest$response <- "forest"
vcmax_cropland <- agg_meta_sen_coef(subset(logr_c_vcmax2,ecosystem.y=="cropland")); vcmax_cropland$response <- "cropland"
vcmax_alltypes <- dplyr::bind_rows(vcmax_ecm,vcmax_am,vcmax_grassland,vcmax_cropland,vcmax_forest)

final_vcmax_prediction <- agg_meta_plots(vcmax_alltypes,"response","logr","logr_var")
final_vcmax_prediction$prediction <- NA
final_vcmax_prediction$prediction[final_vcmax_prediction$type_name=="ECM"] <- mean(final_prediction$pred_vcmax[final_prediction$ecm_type=="ECM"],na.rm=TRUE)
final_vcmax_prediction$prediction[final_vcmax_prediction$type_name=="AM"] <- mean(final_prediction$pred_vcmax[final_prediction$ecm_type=="AM"],na.rm=TRUE)
final_vcmax_prediction$prediction[final_vcmax_prediction$type_name=="grassland"] <- mean(final_prediction$pred_vcmax[final_prediction$ecosystem=="grassland"],na.rm=TRUE)
final_vcmax_prediction$prediction[final_vcmax_prediction$type_name=="forest"] <- mean(final_prediction$pred_vcmax[final_prediction$ecosystem=="forest"],na.rm=TRUE)
final_vcmax_prediction$prediction[final_vcmax_prediction$type_name=="cropland"] <- mean(final_prediction$pred_vcmax[final_prediction$ecosystem=="cropland"],na.rm=TRUE)

final_vcmax_prediction$type_name <- factor(final_vcmax_prediction$type_name, levels = final_vcmax_prediction$type_name)

a2 <- final_vcmax_prediction %>%
  ggplot( aes(x=type_name, y=middle)) + 
  geom_crossbar(aes(x=type_name, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5) +
  geom_boxplot(aes(x=type_name, y=prediction),color="red",size=1)+ 
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-0.6,0.6)+
  labs(x="", y="vcmax") + theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))


#jmax
pred_jmax <- subset(prediction,response=="jmax")
jmax_plotmean <- aggregate(pred_jmax,by=list(pred_jmax$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","pred_jmax25_coef")]
names(jmax_plotmean) <- c("exp","pred_jmax"); jmax_plotmean <- na.omit(jmax_plotmean)

#
final_prediction2 <- Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),list(jmax_all2_final_ecm,jmax_plotmean))
final_prediction2$ecosystem[final_prediction2$ecosystem=="rice_paddies"] <- "cropland"
a3 <- ggplot(subset(final_prediction2,condition=="co2"),aes(x=jmax, y=pred_jmax)) + 
  geom_point(aes(color=ecm_type,shape=ecosystem),size=3)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_text(aes(label=exp),hjust=1, vjust=0,check_overlap = T)+labs(x="Observation jmax", y="Prediction jmax")+
  geom_smooth(color="black",method="lm",se=T)+theme_classic()

logr_c_jmax2 <- Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),list(logr_c_jmax,final_prediction2[,c("exp","ecm_type","ecosystem","pred_jmax")]))
#remove setres_cf as no prediction available
logr_c_jmax2 <- subset(logr_c_jmax2,exp!="setres_cf") # but already not included

jmax_ecm <- agg_meta_sen_coef(subset(logr_c_jmax2,ecm_type=="ECM")); jmax_ecm$response <- "ECM"
jmax_am <- agg_meta_sen_coef(subset(logr_c_jmax2,ecm_type=="AM")); jmax_am$response <- "AM"
jmax_grassland <- agg_meta_sen_coef(subset(logr_c_jmax2,ecosystem.y=="grassland")); jmax_grassland$response <- "grassland"
jmax_forest <- agg_meta_sen_coef(subset(logr_c_jmax2,ecosystem.y=="forest")); jmax_forest$response <- "forest"
jmax_cropland <- agg_meta_sen_coef(subset(logr_c_jmax2,ecosystem.y=="cropland")); jmax_cropland$response <- "cropland"
jmax_alltypes <- dplyr::bind_rows(jmax_ecm,jmax_am,jmax_grassland,jmax_cropland,jmax_forest)

final_jmax_prediction <- agg_meta_plots(jmax_alltypes,"response","logr","logr_var")
final_jmax_prediction$prediction <- NA
final_jmax_prediction$prediction[final_jmax_prediction$type_name=="ECM"] <- mean(final_prediction2$pred_jmax[final_prediction2$ecm_type=="ECM"],na.rm=TRUE)
final_jmax_prediction$prediction[final_jmax_prediction$type_name=="AM"] <- mean(final_prediction2$pred_jmax[final_prediction2$ecm_type=="AM"],na.rm=TRUE)
final_jmax_prediction$prediction[final_jmax_prediction$type_name=="grassland"] <- mean(final_prediction2$pred_jmax[final_prediction2$ecosystem=="grassland"],na.rm=TRUE)
final_jmax_prediction$prediction[final_jmax_prediction$type_name=="forest"] <- mean(final_prediction2$pred_jmax[final_prediction2$ecosystem=="forest"],na.rm=TRUE)
final_jmax_prediction$prediction[final_jmax_prediction$type_name=="cropland"] <- mean(final_prediction2$pred_jmax[final_prediction2$ecosystem=="cropland"],na.rm=TRUE)

final_jmax_prediction$type_name <- factor(final_jmax_prediction$type_name, levels = final_jmax_prediction$type_name)

a4 <- final_jmax_prediction %>%
  ggplot( aes(x=type_name, y=middle)) + 
  geom_crossbar(aes(x=type_name, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5) +
  geom_boxplot(aes(x=type_name, y=prediction),color="red",size=1)+ 
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-0.6,0.6)+
  labs(x="", y="jmax") + theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

#jmax/vcmax
final_prediction3 <- Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),
                       list(subset(final_prediction2,condition=="co2"),
                            subset(final_prediction,condition=="co2")[,c("exp","pred_vcmax")]))
final_prediction3$pred_jmax_vcmax <- final_prediction3$pred_jmax/final_prediction3$pred_vcmax

a5 <- ggplot(final_prediction3,aes(x=jmax_vcmax, y=pred_jmax_vcmax)) + 
  geom_point(aes(color=ecm_type,shape=ecosystem),size=3)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_text(aes(label=exp),hjust=1, vjust=0,check_overlap = T)+labs(x="Observation jmax/vcmax", y="Prediction jmax/vcmax")+
  geom_smooth(color="black",method="lm",se=T)+theme_classic()
a5

a6 <- final_prediction3 %>%
  ggplot( aes(x=ecm_type, y=jmax_vcmax)) + 
  geom_boxplot(aes(x=ecm_type, y=jmax_vcmax),size=1)+ 
  geom_boxplot(aes(x=ecm_type, y=pred_jmax_vcmax),color="red",size=1)+ 
  geom_boxplot(aes(x=ecosystem, y=jmax_vcmax),size=1)+ 
  geom_boxplot(aes(x=ecosystem, y=pred_jmax_vcmax),color="red",size=1)+ 
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-3,3)+
  labs(x="", y="jmax/vcmax") + theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))
a6

plot_grid(a1,a2,a3,a4,a5,a6,nrow=3,label_size = 15)+ theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/final_fig3.jpg",sep=""),width = 15, height = 20)
