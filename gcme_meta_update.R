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

#convert lai_max to lai
kevin_othervars$response[kevin_othervars$response=="lai_max"] <- "lai"

sitename <- kevin_othervars%>% group_by(site)  %>% summarise(number = n())
varname <- kevin_othervars%>% group_by(response)  %>% summarise(number = n())

#leaf_n needs division - but leaf_p doesn't need (as it is all pmass)
nmass <- subset(kevin_othervars,response=="leaf_n") %>%         
  filter(Unit %in% c("g","g_100g","g_g","g_kg","mg_g","mg_kg","mmol_g","ug_mg","umol_g"))
nmass$response <- "nmass"
logr_c_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="c")));
logr_w_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="w")));
logr_f_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="f")));
logr_d_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="d")));
logr_cw_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="cw")));
logr_cf_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="cf")));
logr_cd_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="cd")));
logr_df_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="df")));
logr_dw_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="dw")))

narea <- subset(kevin_othervars,response=="leaf_n") %>%         
  filter(Unit %in% c("g_m2","g_m3","g_pot","mg_cm2","mg_m2","mmol_m2","ug_cm2","ug_cm3","umol_m2"))
narea$response <- "narea"

logr_c_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="c")));
logr_w_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="w")));
logr_f_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="f")));
logr_d_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="d")));
logr_cw_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="cw")));
logr_cf_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="cf")));
logr_cd_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="cd")));
logr_df_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="df")));
logr_dw_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="dw")))

#remove one plot from response ratio of root_shoot_ratio in advance - otherwise it will be NA
kevin_othervars$response[kevin_othervars$response=="root_shoot_ratio" &
                        kevin_othervars$exp=="maricopaface_cotton91_c"] <- NA

#other vars - all created now
for (i in 1:nrow(varname)) {
  tryCatch({
    varname1 <- varname$response[i]
    df_c <- subset(kevin_othervars,treatment=="c" & response==varname1)
    assign(paste("logr_c_", varname1,sep=""), as_tibble(response_ratio_v2(df_c)))
    
    df_f <- subset(kevin_othervars,treatment=="f"& response==varname1)
    assign(paste("logr_f_", varname1,sep=""), as_tibble(response_ratio_v2(df_f)))
    
    df_w <- subset(kevin_othervars,treatment=="w"& response==varname1)
    assign(paste("logr_w_", varname1,sep=""), as_tibble(response_ratio_v2(df_w)))
    
    df_d <- subset(kevin_othervars,treatment=="d"& response==varname1)
    assign(paste("logr_d_", varname1,sep=""), as_tibble(response_ratio_v2(df_d)))
    
    df_cf <- subset(kevin_othervars,treatment=="cf" & response==varname1)
    assign(paste("logr_cf_", varname1,sep=""), as_tibble(response_ratio_v2(df_cf)))
    
    df_cw <- subset(kevin_othervars,treatment=="cw" &response==varname1)
    assign(paste("logr_cw_", varname1,sep=""), as_tibble(response_ratio_v2(df_cw)))
    
    df_fw <- subset(kevin_othervars,treatment=="fw" & response==varname1)
    assign(paste("logr_fw_", varname1,sep=""), as_tibble(response_ratio_v2(df_fw)))
    
    df_cd <- subset(kevin_othervars,treatment=="cd" & response==varname1)
    assign(paste("logr_cd_", varname1,sep=""), as_tibble(response_ratio_v2(df_cd)))
    
    df_df <- subset(kevin_othervars,treatment=="df" & response==varname1)
    assign(paste("logr_df_", varname1,sep=""), as_tibble(response_ratio_v2(df_df)))
    
    df_dw <- subset(kevin_othervars,treatment=="dw" & response==varname1)
    assign(paste("logr_dw_", varname1,sep=""), as_tibble(response_ratio_v2(df_dw)))
    
    df_cfw <- subset(kevin_othervars,treatment=="cfw" & response==varname1)
    assign(paste("logr_cfw_", varname1,sep=""), as_tibble(response_ratio_v2(df_cfw)))
  }, error=function(e){})} 

#anpp filtering with unit
logr_c_anpp <-  logr_c_anpp %>% filter(Unit %in% c("g_m2d","g_m2gs","g_m2y","gc_m2y","kg_ha_y","kg_m2y","kgc_m2y","mg_ha_y","t_ha_y"))
logr_f_anpp <-  logr_f_anpp %>% filter(Unit %in% c("g_m2d","g_m2gs","g_m2y","gc_m2y","kg_ha_y","kg_m2y","kgc_m2y","mg_ha_y","t_ha_y"))
logr_d_anpp <-  logr_d_anpp %>% filter(Unit %in% c("g_m2d","g_m2gs","g_m2y","gc_m2y","kg_ha_y","kg_m2y","kgc_m2y","mg_ha_y","t_ha_y"))
logr_w_anpp <-  logr_w_anpp %>% filter(Unit %in% c("g_m2d","g_m2gs","g_m2y","gc_m2y","kg_ha_y","kg_m2y","kgc_m2y","mg_ha_y","t_ha_y"))
logr_cw_anpp <-  logr_cw_anpp %>% filter(Unit %in% c("g_m2d","g_m2gs","g_m2y","gc_m2y","kg_ha_y","kg_m2y","kgc_m2y","mg_ha_y","t_ha_y"))
logr_cf_anpp <-  logr_cf_anpp %>% filter(Unit %in% c("g_m2d","g_m2gs","g_m2y","gc_m2y","kg_ha_y","kg_m2y","kgc_m2y","mg_ha_y","t_ha_y"))
logr_cd_anpp <-  logr_cd_anpp %>% filter(Unit %in% c("g_m2d","g_m2gs","g_m2y","gc_m2y","kg_ha_y","kg_m2y","kgc_m2y","mg_ha_y","t_ha_y"))
logr_fw_anpp <-  logr_fw_anpp %>% filter(Unit %in% c("g_m2d","g_m2gs","g_m2y","gc_m2y","kg_ha_y","kg_m2y","kgc_m2y","mg_ha_y","t_ha_y"))
logr_df_anpp <-  logr_df_anpp %>% filter(Unit %in% c("g_m2d","g_m2gs","g_m2y","gc_m2y","kg_ha_y","kg_m2y","kgc_m2y","mg_ha_y","t_ha_y"))
logr_dw_anpp <-  logr_dw_anpp %>% filter(Unit %in% c("g_m2d","g_m2gs","g_m2y","gc_m2y","kg_ha_y","kg_m2y","kgc_m2y","mg_ha_y","t_ha_y"))
logr_cfw_anpp <-  logr_cfw_anpp %>% filter(Unit %in% c("g_m2d","g_m2gs","g_m2y","gc_m2y","kg_ha_y","kg_m2y","kgc_m2y","mg_ha_y","t_ha_y"))

#how about combing soil N and soil total N?
#logr_c_soil_n_combined <-  rbind(logr_c_soil_n,logr_c_soil_total_n)
#logr_c_soil_n_combined <-  logr_c_soil_total_n
#logr_c_soil_n_combined$response <- "soil_N"

#something wrong with soil N - just convert
logr_cw_soil_n$logr_var[(logr_cw_soil_n$exp)=="ua2007_cw"] <- NA


#LMA
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
sla <- subset(kevin_LMA,response=="sla")
lma <- subset(kevin_LMA,response=="lma")
#convert SLA to LMA format 
sla$ambient <- 1/sla$ambient; sla$elevated <- 1/sla$elevated; sla$ambient_Sd <- 1/sla$ambient_Sd
sla$elevated_Sd <- 1/sla$elevated_Sd;sla$ambient_Se <- 1/sla$ambient_Se;sla$elevated_Se <- 1/sla$elevated_Se
LMA <- dplyr::bind_rows(lma,sla)
LMA$response <- "LMA"
logr_c_LMA <- as_tibble(response_ratio_v2(subset(LMA,treatment=="c")));
logr_w_LMA <-  as_tibble(response_ratio_v2(subset(LMA,treatment=="w")));
logr_d_LMA <-  as_tibble(response_ratio_v2(subset(LMA,treatment=="d")));
logr_f_LMA <-  as_tibble(response_ratio_v2(subset(LMA,treatment=="f")));
logr_cw_LMA <-  as_tibble(response_ratio_v2(subset(LMA,treatment=="cw")));
logr_cd_LMA <-  as_tibble(response_ratio_v2(subset(LMA,treatment=="cd")));
logr_cf_LMA <-  as_tibble(response_ratio_v2(subset(LMA,treatment=="cf")));
logr_df_LMA <-  as_tibble(response_ratio_v2(subset(LMA,treatment=="df")));
logr_dw_LMA <-  as_tibble(response_ratio_v2(subset(LMA,treatment=="dw")))

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
#kevin$warmQ_e2 # either 0, 1, or 2.2

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

logr_c_vcmax <- as_tibble(response_ratio_v2(subset(kevin_z,treatment=="c"&response=="vcmax")));
logr_w_vcmax <-  as_tibble(response_ratio_v2(subset(kevin_z,treatment=="w"&response=="vcmax")));
logr_f_vcmax <-  as_tibble(response_ratio_v2(subset(kevin_z,treatment=="f"&response=="vcmax")));
logr_d_vcmax <- as_tibble(response_ratio_v2(subset(kevin_z,treatment=="d"&response=="vcmax")));
logr_cw_vcmax <-  as_tibble(response_ratio_v2(subset(kevin_z,treatment=="cw"&response=="vcmax")));
logr_cf_vcmax <-  as_tibble(response_ratio_v2(subset(kevin_z,treatment=="cf"&response=="vcmax")));
logr_cd_vcmax <-  as_tibble(response_ratio_v2(subset(kevin_z,treatment=="cd"&response=="vcmax")));
logr_df_vcmax <-  as_tibble(response_ratio_v2(subset(kevin_z,treatment=="df"&response=="vcmax")));
logr_dw_vcmax <-  as_tibble(response_ratio_v2(subset(kevin_z,treatment=="dw"&response=="vcmax")))

logr_c_jmax <- as_tibble(response_ratio_v2(subset(kevin_z,treatment=="c"&response=="jmax")));
logr_w_jmax <-  as_tibble(response_ratio_v2(subset(kevin_z,treatment=="w"&response=="jmax")));
logr_f_jmax <-  as_tibble(response_ratio_v2(subset(kevin_z,treatment=="f"&response=="jmax")));
logr_d_jmax <-  as_tibble(response_ratio_v2(subset(kevin_z,treatment=="d"&response=="jmax")));
logr_cw_jmax <-  as_tibble(response_ratio_v2(subset(kevin_z,treatment=="cw"&response=="jmax")));
logr_cf_jmax <-  as_tibble(response_ratio_v2(subset(kevin_z,treatment=="cf"&response=="jmax")));
logr_cd_jmax <-  as_tibble(response_ratio_v2(subset(kevin_z,treatment=="cd"&response=="jmax")));
logr_df_jmax <-  as_tibble(response_ratio_v2(subset(kevin_z,treatment=="df"&response=="jmax")));
logr_dw_jmax <-  as_tibble(response_ratio_v2(subset(kevin_z,treatment=="dw"&response=="jmax")))


#now, newly adding bnpp, Nuptake, Asat and npp
df <- read_csv("~/data/gcme/data_received_190325/NewData_wide_CORRECTED2.csv") %>%
  mutate( ambient_Sd  = as.numeric(ambient_Sd),  ambient_Se  = as.numeric(ambient_Se), 
          elevated_Sd = as.numeric(elevated_Sd), elevated_Se = as.numeric(elevated_Se) )

#for anpp, bnpp, Asat, c13
df$response <- NA

#correct unit so that it is not filtered
df$Unit[is.na(df$Unit)==TRUE] <- "NA_unit"

#remove unit g/m2 and g/tree, as not possible in anpp, bnpp, Asat and c13!
df <- subset(df,Unit!="g/m2" & Unit!="g/tree")

#bnpp
#not used: dead_fine_root_biomass, DOC_exudation, exudation, fine_root_density
df$response[df$Data_type=="BNPP"] <- "BNPP"; # all data is correct unit
df$response[df$Data_type=="fine_root_biomass" & df$Unit=="g_/_m2_year"] <- "BNPP";
df$response[df$Data_type=="fine_root_production"& df$Unit!="cm"& df$Unit!="gC/m2"& df$Unit!="km/m3"& df$Unit!="m/m2"] <- "BNPP";
df$response[df$Data_type=="fine_root_turnover"& df$Unit!="km/m3"& df$Unit!="m/m2"] <- "BNPP"

#npp
df$response[df$Data_type=="NPP"] <- "NPP" # all unit is correct

#gpp
df$response[df$Data_type=="GPP"] <- "GPP" # all unit is correct

#wue - take care about unit!
df$response[df$Data_type=="EWUE"|df$Data_type=="instantaneous_water_use_efficiency"|df$Data_type=="leaf-level_instantaneous_water_use_efficiency"|
               df$Data_type=="water_use_efficiency"|df$Data_type=="EWUE"] <- "WUE"

df$response[df$Data_type=="Ci"] <- "ci"

df$response[df$Data_type=="Ci/Ca"] <- "ci_ca"

#anpp - just help to fill 3 plots from popface
df$response[df$Data_type=="ANPP"] <- "ANPP"
#Asat
df$response[df$Data_type=="Asat"] <- "Asat";df$response[df$Data_type=="Amax"] <- "Amax";df$response[df$Data_type=="Anet"] <- "Anet"
#c13
df$response[df$Data_type=="d13C_isotope_ratio_hardwood_leaflitter"] <- "c13";df$response[df$Data_type=="d13C_isotope_ratio_leaflitter"] <- "c13";
df$response[df$Data_type=="leaf_d13C"] <- "c13";df$response[df$Data_type=="leaf_litter_d13C"] <- "c13";df$response[df$Data_type=="litter_stable_13C_isotope_signature"] <- "c13"
#N uptake
df$response[df$Data_type=="N_uptake"] <- "Nuptake";df$response[df$Data_type=="NH4+_uptake"] <- "Nuptake";df$response[df$Data_type=="NO3-_uptake"] <- "Nuptake";
df$response[df$Data_type=="gross_N_mineralization"] <- "Nuptake"

#root-shoot ratio (only defined in below plot - otherwise will be failed in response ratio calculation)
df$response[df$Data_type=="root-shoot_ratio" &(df$exp_nam=="POPFACE_pa"|df$exp_nam=="POPFACE_pe"|df$exp_nam=="POPFACE_pn"|df$exp_nam=="EUROFACE4_pa"|df$exp_nam=="EUROFACE4_pe"|df$exp_nam=="EUROFACE4_pn")] <- "root_shoot"

#gs
df$response[df$Data_type=="gs"|df$Data_type=="stomatal_conductance_(gs)"|df$Data_type=="stomatal_conductance"|df$Data_type=="canopy_conductance"] <- "gs"


#some unit is wrong - but given we had much less plots, and given that N uptake should always be flux?
#subset(df,Data_type=="N_uptake")%>% group_by(Unit)  %>% summarise(number = n()) 

#litter mass
#df$response[df$Data_type=="litter_prod"] <- "litter";df$response[df$Data_type=="litter_mass"] <- "litter"

df_only <- subset(df,is.na(response)==FALSE)

#newly added soil N - copied from beni
selvars <- c("mineral_soil_N", "soil_mineral_N", "soil_NH4-N", "soil_NO3-N","soil_NO3-N_", "soil_solution_mineral_N", "soil_solution_NH4+", "soil_solution_NO3-", "soil_NO3-_", "soil_NO3-", "soil_NH4+", "soil_NH4+-N", "resin_N", "resin_NH4+", "resin_NO3-", "rhizosphere_NO3-N", "rhizosphere_NH4-N")
#selvars2 <- c("mineral_soil_N", "soil_mineral_N", "soil_NH4-N", "soil_NO3-N","soil_NO3-N_",  "soil_NO3-_", "soil_NO3-", "soil_NH4+", "soil_NH4+-N", "resin_N", "resin_NH4+", "resin_NO3-", "rhizosphere_NO3-N", "rhizosphere_NH4-N")

df_only <- df_only %>% 
  bind_rows(
    df %>% 
      filter(Data_type %in% selvars) %>% 
      mutate(response = "ninorg"))

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

df_only[grep("riceface", df_only$exp),]$exp <- gsub("\\,", "", df_only[grep("riceface", df_only$exp),]$exp)
#rice face has more lma, narea,nmass data - needs updated



#3. divide into dataframe
for (i in 1:nrow(varname)) {
  tryCatch({
    varname1 <- varname$response[i]
    df_c <- subset(df_only,treatment=="c" & response==varname1)
    assign(paste("old_logr_c_", varname1,sep=""), as_tibble(response_ratio_v2(df_c)))
    
    df_f <- subset(df_only,treatment=="f"& response==varname1)
    assign(paste("old_logr_f_", varname1,sep=""), as_tibble(response_ratio_v2(df_f)))
    
    df_w <- subset(df_only,treatment=="w"& response==varname1)
    assign(paste("old_logr_w_", varname1,sep=""), as_tibble(response_ratio_v2(df_w)))

    df_d <- subset(df_only,treatment=="d"& response==varname1)
    assign(paste("old_logr_d_", varname1,sep=""), as_tibble(response_ratio_v2(df_d)))
    
    df_cf <- subset(df_only,treatment=="cf" & response==varname1)
    assign(paste("old_logr_cf_", varname1,sep=""), as_tibble(response_ratio_v2(df_cf)))
    
    df_cw <- subset(df_only,treatment=="cw" &response==varname1)
    assign(paste("old_logr_cw_", varname1,sep=""), as_tibble(response_ratio_v2(df_cw)))

    df_cd <- subset(df_only,treatment=="cd" &response==varname1)
    assign(paste("old_logr_cd_", varname1,sep=""), as_tibble(response_ratio_v2(df_cd)))
        
    df_fw <- subset(df_only,treatment=="fw" & response==varname1)
    assign(paste("old_logr_fw_", varname1,sep=""), as_tibble(response_ratio_v2(df_fw)))
    
    df_df <- subset(df_only,treatment=="df" &response==varname1)
    assign(paste("old_logr_df_", varname1,sep=""), as_tibble(response_ratio_v2(df_df)))
    
    df_dw <- subset(df_only,treatment=="dw" & response==varname1)
    assign(paste("old_logr_dw_", varname1,sep=""), as_tibble(response_ratio_v2(df_dw)))
    
    df_cfw <- subset(df_only,treatment=="cfw" & response==varname1)
    assign(paste("old_logr_cfw_", varname1,sep=""), as_tibble(response_ratio_v2(df_cfw)))
  }, error=function(e){})}
# never used old_logr_c_ANPP - since it adds more data from logr_c_anpp!

#check unit
Asat_unit <- old_logr_c_Asat %>% group_by(exp,Unit)  %>% summarise(number = n())
gs_unit <- old_logr_c_gs %>% group_by(exp,Unit)  %>% summarise(number = n())
Asat_gs <- merge(Asat_unit,gs_unit,by=c("exp"),all.x=TRUE)

#correct unit and calculate Ca - Asat/gs
#remove some strange unit that cannot be combined -->now all in umol/m2/s
old_logr_c_Asat$ambient[old_logr_c_Asat$exp=="st_face_ld_c" & (old_logr_c_Asat$Unit=="mmol/g*s" | old_logr_c_Asat$Unit=="mol/g*s")] <- NA
old_logr_c_Asat$elevated[old_logr_c_Asat$exp=="st_face_ld_c" & (old_logr_c_Asat$Unit=="mmol/g*s" | old_logr_c_Asat$Unit=="mol/g*s")] <- NA

old_logr_c_Asat$ambient[old_logr_c_Asat$exp=="st_face_pu_c" & (old_logr_c_Asat$Unit=="mmol/g*s" | old_logr_c_Asat$Unit=="mol/g*s")] <- NA
old_logr_c_Asat$elevated[old_logr_c_Asat$exp=="st_face_pu_c" & (old_logr_c_Asat$Unit=="mmol/g*s" | old_logr_c_Asat$Unit=="mol/g*s")] <- NA
old_logr_c_Asat <- subset(old_logr_c_Asat,is.na(ambient)==FALSE)

#gs: making unit within plot consistent (from mmol to mol)
old_logr_c_gs$ambient[old_logr_c_gs$exp=="euroface4_pe_c" & old_logr_c_gs$Unit=="mmol/m2/s"] <- old_logr_c_gs$ambient[old_logr_c_gs$exp=="euroface4_pe_c" & old_logr_c_gs$Unit=="mmol/m2/s"]/1000
old_logr_c_gs$ambient[old_logr_c_gs$exp=="facts_ii_face3_pt_c" & old_logr_c_gs$Unit=="mmol/m2/s"] <- old_logr_c_gs$ambient[old_logr_c_gs$exp=="facts_ii_face3_pt_c" & old_logr_c_gs$Unit=="mmol/m2/s"]/1000
old_logr_c_gs$ambient[old_logr_c_gs$exp=="facts_ii_face4_bp_c" & old_logr_c_gs$Unit=="mmol/m2/s"] <- old_logr_c_gs$ambient[old_logr_c_gs$exp=="facts_ii_face4_bp_c" & old_logr_c_gs$Unit=="mmol/m2/s"]/1000
old_logr_c_gs$ambient[old_logr_c_gs$exp=="nevada_desert_face_c" & old_logr_c_gs$Unit=="mmol/m2/sec"] <- old_logr_c_gs$ambient[old_logr_c_gs$exp=="nevada_desert_face_c" & old_logr_c_gs$Unit=="mmol/m2/sec"]/1000
old_logr_c_gs$ambient[old_logr_c_gs$exp=="ornerp_liqui2_c" & old_logr_c_gs$Unit=="mmol/m2*s"] <- old_logr_c_gs$ambient[old_logr_c_gs$exp=="ornerp_liqui2_c" & old_logr_c_gs$Unit=="mmol/m2*s"]/1000

#gs: all converting from from H2O to CO2
old_logr_c_gs$ambient[old_logr_c_gs$exp=="brandbjerg_c"] <- old_logr_c_gs$ambient[old_logr_c_gs$exp=="brandbjerg_c"]*44/18
old_logr_c_gs$ambient[old_logr_c_gs$exp=="riceface_china_33n_120e_ar_1_cv"] <- old_logr_c_gs$ambient[old_logr_c_gs$exp=="riceface_china_33n_120e_ar_1_cv"]*44/18
old_logr_c_gs$ambient[old_logr_c_gs$exp=="riceface_china_33n_120e_ar_1_c"] <- old_logr_c_gs$ambient[old_logr_c_gs$exp=="riceface_china_33n_120e_ar_1_c"]*44/18

old_logr_c_gs$elevated[old_logr_c_gs$exp=="euroface4_pe_c" & old_logr_c_gs$Unit=="mmol/m2/s"] <- old_logr_c_gs$elevated[old_logr_c_gs$exp=="euroface4_pe_c" & old_logr_c_gs$Unit=="mmol/m2/s"]/1000
old_logr_c_gs$elevated[old_logr_c_gs$exp=="facts_ii_face3_pt_c" & old_logr_c_gs$Unit=="mmol/m2/s"] <- old_logr_c_gs$elevated[old_logr_c_gs$exp=="facts_ii_face3_pt_c" & old_logr_c_gs$Unit=="mmol/m2/s"]/1000
old_logr_c_gs$elevated[old_logr_c_gs$exp=="facts_ii_face4_bp_c" & old_logr_c_gs$Unit=="mmol/m2/s"] <- old_logr_c_gs$elevated[old_logr_c_gs$exp=="facts_ii_face4_bp_c" & old_logr_c_gs$Unit=="mmol/m2/s"]/1000
old_logr_c_gs$elevated[old_logr_c_gs$exp=="nevada_desert_face_c" & old_logr_c_gs$Unit=="mmol/m2/sec"] <- old_logr_c_gs$elevated[old_logr_c_gs$exp=="nevada_desert_face_c" & old_logr_c_gs$Unit=="mmol/m2/sec"]/1000
old_logr_c_gs$elevated[old_logr_c_gs$exp=="ornerp_liqui2_c" & old_logr_c_gs$Unit=="mmol/m2*s"] <- old_logr_c_gs$elevated[old_logr_c_gs$exp=="ornerp_liqui2_c" & old_logr_c_gs$Unit=="mmol/m2*s"]/1000

old_logr_c_gs$elevated[old_logr_c_gs$exp=="brandbjerg_c"] <- old_logr_c_gs$elevated[old_logr_c_gs$exp=="brandbjerg_c"]*44/18
old_logr_c_gs$elevated[old_logr_c_gs$exp=="riceface_china_33n_120e_ar_1_cv"] <- old_logr_c_gs$elevated[old_logr_c_gs$exp=="riceface_china_33n_120e_ar_1_cv"]*44/18
old_logr_c_gs$elevated[old_logr_c_gs$exp=="riceface_china_33n_120e_ar_1_c"] <- old_logr_c_gs$elevated[old_logr_c_gs$exp=="riceface_china_33n_120e_ar_1_c"]*44/18

#gs data must be converting from mmol to mol
old_logr_c_gs$Unit[old_logr_c_gs$exp=="st_face_ld_c"]
old_logr_c_gs$ambient[old_logr_c_gs$exp=="st_face_ld_c"] <- old_logr_c_gs$ambient[old_logr_c_gs$exp=="st_face_ld_c"]/1000
old_logr_c_gs$ambient[old_logr_c_gs$exp=="st_face_pu_c"] <- old_logr_c_gs$ambient[old_logr_c_gs$exp=="st_face_pu_c"]/1000

old_logr_c_gs$elevated[old_logr_c_gs$exp=="st_face_ld_c"] <- old_logr_c_gs$elevated[old_logr_c_gs$exp=="st_face_ld_c"]/1000
old_logr_c_gs$elevated[old_logr_c_gs$exp=="st_face_pu_c"] <- old_logr_c_gs$elevated[old_logr_c_gs$exp=="st_face_pu_c"]/1000

#hist(old_logr_c_gs$ambient[old_logr_c_gs$exp=="nevada_desert_face_c"])
#already check distribution of all plots --> good

Asat_mean <- aggregate(old_logr_c_Asat,by=list(old_logr_c_Asat$exp,old_logr_c_Asat$Unit), FUN=mean, na.rm=TRUE)[,c("Group.1","co2_a","co2_e","ambient","elevated")]
gs_mean <- aggregate(old_logr_c_gs,by=list(old_logr_c_gs$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","co2_a","co2_e","ambient","elevated")]
Asat_gs_mean <- merge(Asat_mean,gs_mean,by=c("Group.1"),all.x=TRUE)
Asat_gs_mean$ambient_ratio <- Asat_gs_mean$ambient.x/Asat_gs_mean$ambient.y
Asat_gs_mean$elevated_ratio <- Asat_gs_mean$elevated.x/Asat_gs_mean$elevated.y
#since n_plot = 1, magnitude must be wrong but we don't know!
Asat_gs_mean$ambient_ratio[Asat_gs_mean$Group.1=="riceface_china_33n_120e_ar_1_cv"] <- NA
Asat_gs_mean$ambient_ratio[Asat_gs_mean$Group.1=="riceface_china_33n_120e_ar_1_c"] <- NA
Asat_gs_mean$elevated_ratio[Asat_gs_mean$Group.1=="riceface_china_33n_120e_ar_1_cv"] <- NA
Asat_gs_mean$elevated_ratio[Asat_gs_mean$Group.1=="riceface_china_33n_120e_ar_1_c"] <- NA
summary(Asat_gs_mean)

Asat_gs_mean$ambient_ci <- Asat_gs_mean$co2_a.x - Asat_gs_mean$ambient_ratio 
Asat_gs_mean$elevated_ci <- Asat_gs_mean$co2_e.x - Asat_gs_mean$elevated_ratio 

Asat_gs_mean$ci <- log(Asat_gs_mean$elevated_ci/Asat_gs_mean$ambient_ci)/log(Asat_gs_mean$co2_e.x/Asat_gs_mean$co2_a.x) 
Asat_gs_mean_final <- subset(Asat_gs_mean,is.na(ci)==FALSE)[,c("Group.1","ci","co2_a.x","co2_e.x","ambient_ci","elevated_ci")]
names(Asat_gs_mean_final) <- c("exp","ci","co2_a","co2_e","ambient_ci","elevated_ci")
Asat_gs_mean_final$data <- "Asat/gs"

#then ci, ca
old_logr_c_ci_mean <- aggregate(old_logr_c_ci,by=list(old_logr_c_ci$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","co2_a","co2_e","ambient","elevated")]
old_logr_c_ci_mean$ci <- log(old_logr_c_ci_mean$elevated/old_logr_c_ci_mean$ambient)/log(old_logr_c_ci_mean$co2_e/old_logr_c_ci_mean$co2_a) 
old_logr_c_ci_mean <- old_logr_c_ci_mean[,c("Group.1","ci","co2_a","co2_e","ambient","elevated")]
names(old_logr_c_ci_mean) <- c("exp","ci","co2_a","co2_e","ambient_ci","elevated_ci")
old_logr_c_ci_mean$data <- "ci"

old_logr_c_ci_ca_mean  <- aggregate(old_logr_c_ci_ca,by=list(old_logr_c_ci_ca$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","co2_a","co2_e","ambient","elevated")]
old_logr_c_ci_ca_mean$ci <- log((old_logr_c_ci_ca_mean$elevated*old_logr_c_ci_ca_mean$co2_e)/(old_logr_c_ci_ca_mean$ambient*old_logr_c_ci_ca_mean$co2_a))/log(old_logr_c_ci_ca_mean$co2_e/old_logr_c_ci_ca_mean$co2_a) 
old_logr_c_ci_ca_mean$ambient <- old_logr_c_ci_ca_mean$ambient*old_logr_c_ci_ca_mean$co2_a
old_logr_c_ci_ca_mean$elevated <- old_logr_c_ci_ca_mean$elevated*old_logr_c_ci_ca_mean$co2_e

old_logr_c_ci_ca_mean <- old_logr_c_ci_ca_mean[,c("Group.1","ci","co2_a","co2_e","ambient","elevated")]
names(old_logr_c_ci_ca_mean) <- c("exp","ci","co2_a","co2_e","ambient_ci","elevated_ci")
old_logr_c_ci_ca_mean$data <- "ci"
old_logr_c_ci_ca_mean <- subset(old_logr_c_ci_ca_mean,exp=="swissface_trifolium2_c") #remove another since it was already covered

#finally, c13

old_logr_c_c13$delta_ambient <- ((-8.4+1000)/(old_logr_c_c13$ambient+1000)-1)*1000 # -8.4: For the same year the mean atmospheric d13CCO2 reported from Mauna Loa is -8.4‰ (Keeling et al., 2001).
old_logr_c_c13$delta_elevated <- ((-8.4+1000)/(old_logr_c_c13$elevated+1000)-1)*1000 

old_logr_c_c13$ci_ambient <- old_logr_c_c13$co2_a * (old_logr_c_c13$delta_ambient - 4.4)/(27-4.4)
old_logr_c_c13$ci_elevated <- old_logr_c_c13$co2_e * (old_logr_c_c13$delta_elevated - 4.4)/(27-4.4)

old_logr_c_c13_mean <- aggregate(old_logr_c_c13,by=list(old_logr_c_c13$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","co2_a","co2_e","ci_ambient","ci_elevated")]
old_logr_c_c13_mean$ci <- log(old_logr_c_c13_mean$ci_elevated/old_logr_c_c13_mean$ci_ambient)/log(old_logr_c_c13_mean$co2_e/old_logr_c_c13_mean$co2_a) 
old_logr_c_c13_mean <- old_logr_c_c13_mean[,c("Group.1","ci","co2_a","co2_e","ci_ambient","ci_elevated")]
names(old_logr_c_c13_mean) <- c("exp","ci","co2_a","co2_e","ambient_ci","elevated_ci")
old_logr_c_c13_mean$data <- "c13"

ci_final <- dplyr::bind_rows(Asat_gs_mean_final,old_logr_c_ci_mean,old_logr_c_ci_ca_mean,old_logr_c_c13_mean)
#mark 4 replcated plots
ci_final$data[ci_final$exp=="soyfacesoy1_c" & ci_final$data=="Asat/gs"] <- "Asat/gs_replicated"
ci_final$data[ci_final$exp=="soyfacesoy2_c"& ci_final$data=="Asat/gs"] <- "Asat/gs_replicated"
ci_final$data[ci_final$exp=="swissface_trifolium2_c"& ci_final$data=="Asat/gs"] <- "Asat/gs_replicated"
ci_final$data[ci_final$exp=="facts_ii_face3_pt_c"& ci_final$data=="Asat/gs"] <- "Asat/gs_replicated"
ci_final$ci_increase <- ci_final$elevated_ci-ci_final$ambient_ci

ci_final_removal <- subset(ci_final,data!="Asat/gs_replicated")
#now, revising figures
#this function creates to combine cf plots (high-N) into c-only plot
# log (cf/a) - log(f/a) = log(cf/f)
# log (cw/a) - log(w/a) = log(cw/w)
# log (cd/a) - log(d/a) = log(cd/d)
#below is cfwd
combine_co2 <- function(logr_c_var,logr_f_var,logr_w_var,logr_d_var,logr_cf_var,logr_cw_var,logr_cd_var,name){
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
  new_df$condition <- "(co2 + Nfer)/Nfer"
    
  #co2 + warming / warming
  all_logr_w_lma <- agg_meta(logr_w_var)[,c("exp","middle")]
  all_logr_cw_lma <- agg_meta(logr_cw_var)[,c("exp","middle")]
  all_logr_cw_lma$exp_old <- all_logr_cw_lma$exp 
  all_logr_cw_lma$exp <-  paste(substr(all_logr_cw_lma$exp,1,nchar(all_logr_cw_lma$exp)-2),"w",sep="")
  all_logr_w_lma$exp;all_logr_cw_lma$exp
  all_logr_w_lma2 <- all_logr_w_lma %>% filter(exp %in%all_logr_cw_lma$exp)
  new_df2 <- merge(all_logr_w_lma2,all_logr_cw_lma,by=c("exp"),all.x=TRUE)
  for (i in 1:nrow(new_df2)) {
    new_df2$co2[i] <- subset(logr_cw_var,exp==new_df2$exp_old[i])$co2_e[1]-
      subset(logr_cw_var,exp==new_df2$exp_old[i])$co2_a[1]
  } 
  new_df2$middle <- (new_df2$middle.y - new_df2$middle.x)/log(new_df2$co2)
  new_df2$exp <- new_df2$exp_old
  new_df2$condition <- "(co2 + warming)/warming"

  #co2 + drought / drought
  all_logr_d_lma <- agg_meta(logr_d_var)[,c("exp","middle")]
  all_logr_cd_lma <- agg_meta(logr_cd_var)[,c("exp","middle")]
  all_logr_cd_lma$exp_old <- all_logr_cd_lma$exp 
  all_logr_cd_lma$exp <-  paste(substr(all_logr_cd_lma$exp,1,nchar(all_logr_cd_lma$exp)-2),"d",sep="")
  all_logr_d_lma$exp;all_logr_cd_lma$exp
  all_logr_d_lma2 <- all_logr_d_lma %>% filter(exp %in%all_logr_cd_lma$exp)
  new_df3 <- merge(all_logr_d_lma2,all_logr_cd_lma,by=c("exp"),all.x=TRUE)
  for (i in 1:nrow(new_df3)) {
    new_df3$co2[i] <- subset(logr_cd_var,exp==new_df3$exp_old[i])$co2_e[1]-
      subset(logr_cd_var,exp==new_df3$exp_old[i])$co2_a[1]
  } 
  new_df3$middle <- (new_df3$middle.y - new_df3$middle.x)/log(new_df3$co2) # log (cd/a) - log(d/a) = log(cd/d)
  new_df3$exp <- new_df3$exp_old
  new_df3$condition <- "(co2 + drought)/drought"
  
  lma_plot <- dplyr::bind_rows(all_logr_c_lma,new_df[,c("exp","middle","condition")],new_df2[,c("exp","middle","condition")],new_df3[,c("exp","middle","condition")]) 
  names(lma_plot) <- c("exp",name,"condition")
  return(lma_plot)
}

#below is cfd
combine_co2_cfd <- function(logr_c_var,logr_f_var,logr_w_var,logr_d_var,logr_cf_var,logr_cw_var,logr_cd_var,name){
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
  new_df$condition <- "(co2 + Nfer)/Nfer"
  
  #co2 + drought / drought
  all_logr_d_lma <- agg_meta(logr_d_var)[,c("exp","middle")]
  all_logr_cd_lma <- agg_meta(logr_cd_var)[,c("exp","middle")]
  all_logr_cd_lma$exp_old <- all_logr_cd_lma$exp 
  all_logr_cd_lma$exp <-  paste(substr(all_logr_cd_lma$exp,1,nchar(all_logr_cd_lma$exp)-2),"d",sep="")
  all_logr_d_lma$exp;all_logr_cd_lma$exp
  all_logr_d_lma2 <- all_logr_d_lma %>% filter(exp %in%all_logr_cd_lma$exp)
  new_df3 <- merge(all_logr_d_lma2,all_logr_cd_lma,by=c("exp"),all.x=TRUE)
  for (i in 1:nrow(new_df3)) {
    new_df3$co2[i] <- subset(logr_cd_var,exp==new_df3$exp_old[i])$co2_e[1]-
      subset(logr_cd_var,exp==new_df3$exp_old[i])$co2_a[1]
  } 
  new_df3$middle <- (new_df3$middle.y - new_df3$middle.x)/log(new_df3$co2) # log (cd/a) - log(d/a) = log(cd/d)
  new_df3$exp <- new_df3$exp_old
  new_df3$condition <- "(co2 + drought)/drought"
  
  lma_plot <- dplyr::bind_rows(all_logr_c_lma,new_df[,c("exp","middle","condition")],new_df3[,c("exp","middle","condition")]) 
  names(lma_plot) <- c("exp",name,"condition")
  return(lma_plot)
}

#below is cfw
combine_co2_cfw <- function(logr_c_var,logr_f_var,logr_w_var,logr_d_var,logr_cf_var,logr_cw_var,logr_cd_var,name){
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
  new_df$condition <- "(co2 + Nfer)/Nfer"
  #co2 + warming / warming
  all_logr_w_lma <- agg_meta(logr_w_var)[,c("exp","middle")]
  all_logr_cw_lma <- agg_meta(logr_cw_var)[,c("exp","middle")]
  all_logr_cw_lma$exp_old <- all_logr_cw_lma$exp 
  all_logr_cw_lma$exp <-  paste(substr(all_logr_cw_lma$exp,1,nchar(all_logr_cw_lma$exp)-2),"w",sep="")
  all_logr_w_lma$exp;all_logr_cw_lma$exp
  all_logr_w_lma2 <- all_logr_w_lma %>% filter(exp %in%all_logr_cw_lma$exp)
  new_df2 <- merge(all_logr_w_lma2,all_logr_cw_lma,by=c("exp"),all.x=TRUE)
    for (i in 1:nrow(new_df2)) {
      new_df2$co2[i] <- subset(logr_cw_var,exp==new_df2$exp_old[i])$co2_e[1]-
        subset(logr_cw_var,exp==new_df2$exp_old[i])$co2_a[1]
    } 
    new_df2$middle <- (new_df2$middle.y - new_df2$middle.x)/log(new_df2$co2)
    new_df2$exp <- new_df2$exp_old
    new_df2$condition <- "(co2 + warming)/warming"
  
  lma_plot <- dplyr::bind_rows(all_logr_c_lma,new_df[,c("exp","middle","condition")],new_df2[,c("exp","middle","condition")]) 
  names(lma_plot) <- c("exp",name,"condition")
  return(lma_plot)
}

#below is cf
combine_co2_cf <- function(logr_c_var,logr_f_var,logr_w_var,logr_d_var,logr_cf_var,logr_cw_var,logr_cd_var,name){
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
  new_df$condition <- "(co2 + Nfer)/Nfer"
  
  lma_plot <- dplyr::bind_rows(all_logr_c_lma,new_df[,c("exp","middle","condition")]) 
  names(lma_plot) <- c("exp",name,"condition")
  return(lma_plot)
}

combine_co2_c <- function(logr_c_var,logr_f_var,logr_w_var,logr_d_var,logr_cf_var,logr_cw_var,logr_cd_var,name){
  all_logr_c_lma <- agg_meta_sen_coef(logr_c_var)[,c("exp","middle")]
  all_logr_c_lma$condition <- "co2"
  lma_plot <- all_logr_c_lma
  names(lma_plot) <- c("exp",name,"condition")
  return(lma_plot)
}

lma_plot <- combine_co2_c(logr_c_LMA,logr_f_LMA,logr_w_LMA,logr_d_LMA,logr_cf_LMA,logr_cw_LMA,logr_cd_LMA,"LMA")
vcmax_plot <- combine_co2(logr_c_vcmax,logr_f_vcmax,logr_w_vcmax,logr_d_vcmax,logr_cf_vcmax,logr_cw_vcmax,logr_cd_vcmax,"vcmax")
jmax_plot <- combine_co2(logr_c_jmax,logr_f_jmax,logr_w_jmax,logr_d_jmax,logr_cf_jmax,logr_cw_jmax,logr_cd_jmax,"jmax")
narea_plot <- combine_co2(logr_c_narea,logr_f_narea,logr_w_narea,logr_d_narea,logr_cf_narea,logr_cw_narea,logr_cd_narea,"narea")
nmass_plot <- combine_co2(logr_c_nmass,logr_f_nmass,logr_w_nmass,logr_d_nmass,logr_cf_nmass,logr_cw_nmass,logr_cd_nmass,"nmass")
leaf_cn_plot <- combine_co2(logr_c_leaf_cn,logr_f_leaf_cn,logr_w_leaf_cn,logr_d_leaf_cn,logr_cf_leaf_cn,logr_cw_leaf_cn,logr_cd_leaf_cn,"leaf_cn")
rm(logr_cw_anpp) #it is nothing available here but affects function
anpp_plot <- combine_co2_cfd(logr_c_anpp,logr_f_anpp,logr_w_anpp,logr_d_anpp,logr_cf_anpp,logr_cw_anpp,logr_cd_anpp,"anpp")
lai_plot <- combine_co2(logr_c_lai,logr_f_lai,logr_w_lai,logr_d_lai,logr_cf_lai,logr_cw_lai,logr_cd_lai,"lai")
bnpp_plot <- combine_co2_c(old_logr_c_BNPP,old_logr_f_BNPP,old_logr_w_BNPP,old_logr_d_BNPP,old_logr_cf_BNPP,old_logr_cw_BNPP,old_logr_cd_BNPP,"bnpp")
Nuptake_plot <- combine_co2(old_logr_c_Nuptake,old_logr_f_Nuptake,old_logr_w_Nuptake,old_logr_d_Nuptake,old_logr_cf_Nuptake,old_logr_cw_Nuptake,old_logr_cd_Nuptake,"Nuptake")
npp_plot <- combine_co2_cf(old_logr_c_NPP,old_logr_f_NPP,old_logr_w_NPP,old_logr_d_NPP,old_logr_cf_NPP,old_logr_cw_NPP,old_logr_cd_NPP,"npp")
soilN_plot <- combine_co2(logr_c_soil_n,logr_f_soil_n,logr_w_soil_n,logr_d_soil_n,logr_cf_soil_n,logr_cw_soil_n,logr_cd_soil_n,"soilN")
soil_total_N_plot <- combine_co2_cfw(logr_c_soil_total_n,logr_f_soil_total_n,logr_w_soil_total_n,logr_d_soil_total_n,logr_cf_soil_total_n,logr_cw_soil_total_n,logr_cd_soil_total_n,"soil_totalN")
#root_shoot_plot <- combine_co2_c(logr_c_root_shoot_ratio,logr_f_root_shoot_ratio,logr_w_root_shoot_ratio,logr_d_root_shoot_ratio,logr_cf_root_shoot_ratio,logr_cw_root_shoot_ratio,logr_cd_root_shoot_ratio,"root_shoot") 
old_root_shoot_plot <- combine_co2_c(old_logr_c_root_shoot,old_logr_f_root_shoot,old_logr_w_root_shoot,old_logr_d_root_shoot,old_logr_cf_root_shoot,old_logr_cw_root_shoot,old_logr_cd_root_shoot,"root_shoot") 
old_anpp_plot <- combine_co2_cf(old_logr_c_ANPP,old_logr_f_ANPP,old_logr_w_ANPP,old_logr_d_ANPP,old_logr_cf_ANPP,old_logr_cw_ANPP,old_logr_cd_ANPP,"old_anpp")
#some are 0 which is wrong
old_logr_c_ninorg$logr[old_logr_c_ninorg$logr==0] <- NA
old_ninorg_plot <- combine_co2_cfw(old_logr_c_ninorg,old_logr_f_ninorg,old_logr_w_ninorg,old_logr_d_ninorg,old_logr_cf_ninorg,old_logr_cw_ninorg,old_logr_cd_ninorg,"old_ninorg")

#gpp
gpp_plot <- combine_co2_c(old_logr_c_GPP,old_logr_f_GPP,old_logr_w_GPP,old_logr_d_GPP,old_logr_cf_GPP,old_logr_cw_GPP,old_logr_cd_GPP,"gpp")

#Asat, Amax, c13 and gs
Asat_plot <- combine_co2(old_logr_c_Asat,old_logr_f_Asat,old_logr_w_Asat,old_logr_d_Asat,old_logr_cf_Asat,old_logr_cw_Asat,old_logr_cd_Asat,"Asat")
Amax_plot <- combine_co2_c(old_logr_c_Amax,old_logr_f_Amax,old_logr_w_Amax,old_logr_d_Amax,old_logr_cf_Amax,old_logr_cw_Amax,old_logr_cd_Amax,"Amax")
c13_plot <- combine_co2_c(old_logr_c_c13,old_logr_f_c13,old_logr_w_c13,old_logr_d_c13,old_logr_cf_c13,old_logr_cw_c13,old_logr_cd_c13,"c13")
gs_plot <- combine_co2(old_logr_c_gs,old_logr_f_gs,old_logr_w_gs,old_logr_d_gs,old_logr_cf_gs,old_logr_cw_gs,old_logr_cd_gs,"gs")

#wue, ci, c1
wue_plot <- combine_co2_c(old_logr_c_WUE,old_logr_f_WUE,old_logr_w_WUE,old_logr_d_WUE,old_logr_cf_WUE,old_logr_cw_WUE,old_logr_cd_WUE,"wue")
ci_plot <- combine_co2_c(old_logr_c_ci,old_logr_f_ci,old_logr_w_ci,old_logr_d_ci,old_logr_cf_ci,old_logr_cw_ci,old_logr_cd_ci,"ci")
ci_ca_plot <- combine_co2_c(old_logr_c_ci_ca,old_logr_f_ci_ca,old_logr_w_ci_ca,old_logr_d_ci_ca,old_logr_cf_ci_ca,old_logr_cw_ci_ca,old_logr_cd_ci_ca,"ci_ca")

vcmax_plot_nfer <- combine_co2_c(logr_cf_vcmax,logr_f_vcmax,logr_w_vcmax,logr_d_vcmax,logr_cf_vcmax,logr_cw_vcmax,logr_cd_vcmax,"vcmax")
#check N fertilization's simple effect
vcmax_plot %>% filter(exp %in% c("duke2_c","euroface4_pa_c","euroface4_pe_c","euroface4_pn_c","new_zealand_face_c"))
vcmax_plot_nfer

#only additionally fill popface's 3 plots
new_anpp <- subset(old_anpp_plot,exp=="popface_pa_c"|exp=="popface_pe_c"|exp=="popface_pn_c")
names(new_anpp) <- c("exp","anpp","condition")
anpp_plot <- rbind(anpp_plot,new_anpp)

vcmax_main <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp","condition"),all.x=TRUE),
                    list(vcmax_plot,jmax_plot,anpp_plot,lma_plot,narea_plot,nmass_plot,
                         leaf_cn_plot,lai_plot,bnpp_plot,Nuptake_plot,npp_plot,soilN_plot,soil_total_N_plot,old_root_shoot_plot))

anpp_compare <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp","condition"),all.x=TRUE),
                    list(anpp_plot,old_anpp_plot))
#needs double check
vcmax_main$soilN[vcmax_main$soilN==0] <- NA
vcmax_main$soilN[vcmax_main$soilN< -4] <- NA

vcmax_main$combined_soilN <- vcmax_main$soilN
vcmax_main$combined_soilN[is.na(vcmax_main$soil_totalN)==FALSE] <- vcmax_main$soil_totalN[is.na(vcmax_main$soil_totalN)==FALSE]
vcmax_main$nuptake_bnpp <-vcmax_main$Nuptake - vcmax_main$bnpp
vcmax_main$bnpp_anpp <-vcmax_main$bnpp - vcmax_main$anpp

p <- list()
for(i in c(1:(ncol(vcmax_main)-4))){
  p[[i]] <- ggplot(vcmax_main,aes_string(x=names(vcmax_main)[i+4],
                                         y=names(vcmax_main)[3])) +
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=condition),size=3)+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))}

#bnpp/anpp and root/shoot
aa1 <- p[[15]]+geom_text(aes(label=exp,hjust=0, vjust=-1,check_overlap = F))+theme(legend.position = "none")
aa2 <- p[[12]]+geom_text(aes(label=exp,hjust=0, vjust=-1,check_overlap = F))+theme(legend.position = "none")
plot_grid(aa1,aa2,nrow=1,label_size = 15)+
  theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/egu_update_2figs.jpg",sep=""),width = 20, height = 10)

plot_grid(p[[4]],p[[3]],p[[1]],p[[8]],p[[14]],p[[7]],p[[6]],p[[13]],p[[9]],nrow=3,label_size = 15)+
  theme(plot.background=element_rect(fill="white", color="white"))

ggsave(paste("~/data/output_gcme/colin/egu_update_v.jpg",sep=""),width = 20, height = 10)


jmax_main <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp","condition"),all.x=TRUE),
                    list(jmax_plot,vcmax_plot,anpp_plot,lma_plot,narea_plot,nmass_plot,
                         leaf_cn_plot,lai_plot,bnpp_plot,Nuptake_plot,npp_plot,soilN_plot,soil_total_N_plot))
jmax_main$soilN[jmax_main$soilN==0] <- NA
jmax_main$soilN[jmax_main$soilN< -4] <- NA

jmax_main$combined_soilN <- jmax_main$soilN
jmax_main$combined_soilN[is.na(jmax_main$soil_totalN)==FALSE] <- jmax_main$soil_totalN[is.na(jmax_main$soil_totalN)==FALSE]
jmax_main$nuptake_bnpp <-jmax_main$Nuptake - jmax_main$bnpp

p <- list()
for(i in c(1:(ncol(jmax_main)-4))){
  p[[i]] <- ggplot(jmax_main,aes_string(x=names(jmax_main)[i+4],
                                         y=names(jmax_main)[3])) +
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=condition),size=3)+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))}

plot_grid(p[[4]],p[[3]],p[[1]],p[[8]],p[[13]],p[[7]],p[[6]],p[[12]],p[[9]],nrow=3,label_size = 15)+
  theme(plot.background=element_rect(fill="white", color="white"))

ggsave(paste("~/data/output_gcme/colin/egu_update_j.jpg",sep=""),width = 20, height = 10)

#now, lai ~anpp

lai_main <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp","condition"),all.x=TRUE),
                    list(lai_plot,vcmax_plot,jmax_plot,anpp_plot,lma_plot,narea_plot,nmass_plot,
                         leaf_cn_plot,bnpp_plot,Nuptake_plot,npp_plot,soilN_plot,soil_total_N_plot))
lai_main$soilN[lai_main$soilN==0] <- NA
lai_main$soilN[lai_main$soilN< -4] <- NA

lai_main$combined_soilN <- lai_main$soilN
lai_main$combined_soilN[is.na(lai_main$soil_totalN)==FALSE] <- lai_main$soil_totalN[is.na(lai_main$soil_totalN)==FALSE]
lai_main$nuptake_bnpp <-lai_main$Nuptake - lai_main$bnpp

p <- list()
for(i in c(1:(ncol(lai_main)-4))){
  p[[i]] <- ggplot(lai_main,aes_string(x=names(lai_main)[i+3],
                                         y=names(lai_main)[3])) +
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=condition),size=3)+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))}

plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],
          p[[6]],p[[7]],p[[8]],p[[9]],
          p[[10]],p[[13]],
          nrow=4,label_size = 15)+
  theme(plot.background=element_rect(fill="white", color="white"))

ggsave(paste("~/data/output_gcme/colin/egu_update_lai.jpg",sep=""),width = 20, height = 15)

lai1<- ggplot(vcmax_main,aes_string(x="vcmax",
                           y="anpp")) +
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=lai),size=3)+
  geom_text(aes(label=format(round(lai, 3), nsmall = 3)),hjust=1, vjust=0,check_overlap = T)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="black",method="lm",se=F)

lai2 <- ggplot(vcmax_main,aes_string(x="vcmax",
                                    y="bnpp")) +
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=lai),size=3)+
  geom_text(aes(label=format(round(lai, 3), nsmall = 3)),hjust=1, vjust=0,check_overlap = T)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="black",method="lm",se=F)

lai3<- ggplot(vcmax_main,aes_string(x="vcmax",
                             y="nmass")) +
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=lai),size=3)+
  geom_text(aes(label=format(round(lai, 3), nsmall = 3)),hjust=1, vjust=0,check_overlap = F)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="black",method="lm",se=F)

plot_grid(lai1,lai2,lai3,nrow=1,label_size = 15)+
  theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/egu_update_lai2.jpg",sep=""),width = 15, height = 5)

#now, soil inorganic n ~ vcmax, jmax, nmass
ninorg_main <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp","condition"),all.x=TRUE),
                  list(old_ninorg_plot,vcmax_plot,jmax_plot,nmass_plot))
#remove biocon
#ninorg_main$old_ninorg[ninorg_main$old_ninorg< - 0.4] <- NA
p <- list()
for(i in c(1:3)){
  p[[i]] <- ggplot(ninorg_main,aes_string(y=names(ninorg_main)[i+3],
                                       x=names(ninorg_main)[3])) +
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=condition),size=3)+
    geom_text(aes(label=exp,hjust=1, vjust=0,check_overlap = T))+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+labs(x="Inorganic soil N")+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))}
  
plot_grid(p[[1]],p[[2]],p[[3]],nrow=1,label_size = 15)+
  theme(plot.background=element_rect(fill="white", color="white"))

ggsave(paste("~/data/output_gcme/colin/egu_update_soil.jpg",sep=""),width = 20, height = 5)

#newly combine from combination
#have a look at soil data
avail_N <- subset(df_only,response=="ninorg") %>% filter(exp %in%vcmax_main$exp)
avail_N2 <- avail_N[,c("exp","ambient","elevated","co2_a","co2_e","Data_type","Unit")]
avail_N2$variable <- NA
avail_N2$variable[avail_N2$Data_type=="mineral_soil_N"|avail_N2$Data_type=="soil_mineral_N"|avail_N2$Data_type=="resin_N"] <- "mineral_N"
avail_N2$variable[avail_N2$Data_type=="soil_NH4"|avail_N2$Data_type=="soil_NH4+"|avail_N2$Data_type=="soil_NH4+-N"|avail_N2$Data_type=="soil_NH4-N"] <- "NH4"
avail_N2$variable[avail_N2$Data_type=="soil_NO3-N"|avail_N2$Data_type=="soil_NO3-N_"|avail_N2$Data_type=="soil_NO3-_"|avail_N2$Data_type=="soil_NO3-"] <- "NO3"
#now already exclude solution N, and |avail_N2$Data_type=="resin_NH4+"|avail_N2$Data_type=="rhizosphere_NH4-N"

avail_N2 <- subset(avail_N2,is.na(variable)==FALSE)
avail_N3 <- aggregate(avail_N2,by=list(avail_N2$exp,avail_N2$variable), FUN=mean, na.rm=TRUE)
avail_N3$exp <- avail_N3$Group.1;avail_N3$variable <- avail_N3$Group.2; avail_N3 <- avail_N3[,c("exp","ambient","elevated","co2_a","co2_e","variable")]
avail_option <- subset(avail_N3,exp=="biocon_c"|exp=="brandbjerg_c"|exp=="duke2_c"|exp=="giface_c")
#perfectly divided by 4 + 4 --> combine them
avail_option
avail_option[1,2] <- avail_option[1,2] +avail_option[5,2];avail_option[1,3] <- avail_option[1,3] +avail_option[5,3];
avail_option[2,2] <- avail_option[2,2] +avail_option[6,2];avail_option[2,3] <- avail_option[2,3] +avail_option[6,3];
avail_option[3,2] <- avail_option[3,2] +avail_option[7,2];avail_option[3,3] <- avail_option[3,3] +avail_option[7,3];
avail_option[4,2] <- avail_option[4,2] +avail_option[8,2];avail_option[4,3] <- avail_option[4,3] +avail_option[8,3]
avail_option_final <- avail_option[1:4,]
avail_option_final$variable <- "NH4+NO3"
soil_mineral <- rbind(avail_N3[1:7,],avail_option_final)
soil_mineral$mineral_soilN <- log(soil_mineral$elevated/soil_mineral$ambient)/log(soil_mineral$co2_e/soil_mineral$co2_a)

soil_NH4 <- subset(avail_N3,variable=="NH4")
soil_NH4$NH4 <- log(soil_NH4$elevated/soil_NH4$ambient)/log(soil_NH4$co2_e/soil_NH4$co2_a)

soil_NO3 <- subset(avail_N3,variable=="NO3")
soil_NO3$NO3 <- log(soil_NO3$elevated/soil_NO3$ambient)/log(soil_NO3$co2_e/soil_NO3$co2_a)

new_soil_vcmax <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),
                     list(vcmax_main,soil_mineral[c("exp","variable","mineral_soilN")],soil_NH4[,c("exp","NH4")],soil_NO3[,c("exp","NO3")]))

ggplot(new_soil_vcmax,aes_string(y="vcmax",x="mineral_soilN")) +
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=variable),size=3)+
  geom_text(aes(label=exp,hjust=1, vjust=0,check_overlap = T))+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="black",method="lm",se=F)+labs(x="Inorganic soil N")+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))
ggsave(paste("~/data/output_gcme/colin/egu_update_soil_final.jpg",sep=""),width = 10, height = 5)


#now, combined with prediction data from gcme_vcmax
#the prediction of vcmax only includes c, w and cw - but having same value in cf
prediction <- read.csv("/Users/yunpeng/data/gcme/kevin/forcing/pred_vcmax.csv")
pred_vcmax <- subset(prediction,response=="vcmax")
vcmax_plotmean <- aggregate(pred_vcmax,by=list(pred_vcmax$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","pred_vcmax25_coef")]
names(vcmax_plotmean) <- c("exp","pred_vcmax"); vcmax_plotmean <- na.omit(vcmax_plotmean)

#the prediction for c, cf and cd is correct. But cw is not correct
#merged with pft and ecm types
ecm_csv <- read.csv("/Users/yunpeng/data/gcme/kevin/orig_vcmax/new_ecm_types_updated.csv")
vcmax_obs <- merge(vcmax_main,ecm_csv,by=c("exp"),all.x=TRUE)
vcmax_obs$ecosystem[vcmax_obs$ecosystem=="heathland"] <- "grassland"
vcmax_obs$ecosystem[vcmax_obs$ecosystem=="temperate_forest"] <- "forest"
vcmax_obs$ecosystem[vcmax_obs$ecosystem=="shrubland"] <- "forest"

vcmax_obs_pred <- Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),
                         list(vcmax_obs,vcmax_plotmean))

#correct 3 co2+warming effects prediction = (logr_cw - logr_w)/log(co2 increase). 
#This fits observation!
# riceface_japan_l_2008_3938_14057_cw's warming effect is 0 - needs revisitsed predction
# what is warmQ_e1, warmQ_e2, warmQ_e3?
#vcmax_obs_pred$pred_vcmax[vcmax_obs_pred$exp=="brandbjerg_cw"] <- (subset(pred_vcmax,exp=="brandbjerg_cw")$pred_vcmax25_logr[1]-subset(pred_vcmax,exp=="brandbjerg_w")$pred_vcmax25_logr[1])/log(subset(pred_vcmax,exp=="brandbjerg_cw")$co2_e[1]-subset(pred_vcmax,exp=="brandbjerg_cw")$co2_a[1])
#vcmax_plotmean$pred_vcmax[vcmax_obs_pred$exp=="riceface_japan_l_2007_3938_14057_cw"] <- (subset(pred_vcmax,exp=="riceface_japan_l_2007_3938_14057_cw")$pred_vcmax25_logr[1]-subset(pred_vcmax,exp=="riceface_japan_l_2007_3938_14057_w")$pred_vcmax25_logr[1])/log(subset(pred_vcmax,exp=="riceface_japan_l_2007_3938_14057_cw")$co2_e[1]-subset(pred_vcmax,exp=="riceface_japan_l_2007_3938_14057_cw")$co2_a[1])
#vcmax_plotmean$pred_vcmax[vcmax_obs_pred$exp=="riceface_japan_l_2008_3938_14057_cw"] <- (subset(pred_vcmax,exp=="riceface_japan_l_2008_3938_14057_cw")$pred_vcmax25_logr[1]-subset(pred_vcmax,exp=="riceface_japan_l_2008_3938_14057_w")$pred_vcmax25_logr[1])/log(subset(pred_vcmax,exp=="riceface_japan_l_2008_3938_14057_cw")$co2_e[1]-subset(pred_vcmax,exp=="riceface_japan_l_2008_3938_14057_cw")$co2_a[1])

#but only include co2-effect only now
vcmax_obs_pred_co2 <- subset(vcmax_obs_pred,condition=="co2")
vcmax_obs_pred_co2$response <- "all"
dim(vcmax_obs_pred_co2)
logr_c_vcmax2 <- Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),
                        list(logr_c_vcmax,vcmax_obs_pred_co2[,c("exp","ecm_type","ecosystem","pred_vcmax")]))
#remove setres_cf as no prediction available
#logr_c_vcmax2 <- subset(logr_c_vcmax2,exp!="setres_cf") # but already not included

vcmax_grassland <- agg_meta_sen_coef(subset(logr_c_vcmax2,ecosystem.y=="grassland")); vcmax_grassland$response <- "grassland"
vcmax_forest <- agg_meta_sen_coef(subset(logr_c_vcmax2,ecosystem.y=="forest")); vcmax_forest$response <- "forest"
vcmax_cropland <- agg_meta_sen_coef(subset(logr_c_vcmax2,ecosystem.y=="cropland")); vcmax_cropland$response <- "cropland"
vcmax_alltypes <- dplyr::bind_rows(vcmax_grassland,vcmax_cropland,vcmax_forest)
#remove a few sites with no available logr_var??? - but not applied since it not affects sen coef
final_vcmax <- agg_meta_plots(vcmax_alltypes,"response","logr","logr_var")
#created a unified one and help to combine to preediction column then 
vcmax_alltypes2 <- vcmax_alltypes; vcmax_alltypes2$response <- "alltypes"

final_observation <- agg_meta_plots(vcmax_alltypes2,"response","logr","logr_var")
final_observation$type_name <- "all"
#aggregate by mean is more closed to prediction???
#final_observation$middle<- mean(vcmax_alltypes2$middle)
#final_observation$ymin<- mean(vcmax_alltypes2$middle)
#final_observation$ymax<- mean(vcmax_alltypes2$middle)

pred_plotmean <- mean(vcmax_obs_pred_co2$pred_vcmax)
final_prediction <- tibble(type_name="all",middle=pred_plotmean,ymin=pred_plotmean,ymax=pred_plotmean)

a1 <- final_vcmax %>%
  ggplot( aes(x=type_name, y=middle)) + 
  geom_crossbar(aes(x=type_name, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5) +
  geom_crossbar(data=final_prediction,aes(x=type_name, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="red") +
  geom_crossbar(data=final_observation,aes(x=type_name, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="black") +
  geom_point(data=vcmax_alltypes,aes(x=response, y=middle, size= 1.96/(middle-ymin)), alpha = 0.6, width = 0.5) +#geom_point(data=subset(vcmax_alltypes,is.na(middle_scaled)==TRUE),aes(x=response, y=middle,size=5,color="no_SE_info"),width = 0.5) +
  geom_point(data=vcmax_obs_pred_co2,aes(x=response, y=pred_vcmax), alpha = 0.6, width = 0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(x="", y="Sensitivity coefficient of vcmax",size=expression(paste("Standard Error"^{-1}))) +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))
a1

#jmax
pred_jmax <- subset(prediction,response=="jmax")
jmax_plotmean <- aggregate(pred_jmax,by=list(pred_jmax$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","pred_jmax25_coef")]
names(jmax_plotmean) <- c("exp","pred_jmax"); jmax_plotmean <- na.omit(jmax_plotmean)

#the prediction for c, cf and cd is correct. But cw is not correct
#merged with pft and ecm types
jmax_obs <- merge(jmax_main,ecm_csv,by=c("exp"),all.x=TRUE)
jmax_obs$ecosystem[jmax_obs$ecosystem=="heathland"] <- "grassland"
jmax_obs$ecosystem[jmax_obs$ecosystem=="temperate_forest"] <- "forest"
jmax_obs$ecosystem[jmax_obs$ecosystem=="shrubland"] <- "forest"

jmax_obs_pred <- Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),
                        list(jmax_obs,jmax_plotmean))


#but only include co2-effect only now
jmax_obs_pred_co2 <- subset(jmax_obs_pred,condition=="co2")
jmax_obs_pred_co2$response <- "all"
dim(jmax_obs_pred_co2)
logr_c_jmax2 <- Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),
                       list(logr_c_jmax,jmax_obs_pred_co2[,c("exp","ecm_type","ecosystem","pred_jmax")]))
#remove setres_cf as no prediction available
#logr_c_vcmax2 <- subset(logr_c_vcmax2,exp!="setres_cf") # but already not included

jmax_grassland <- agg_meta_sen_coef(subset(logr_c_jmax2,ecosystem.y=="grassland")); jmax_grassland$response <- "grassland"
jmax_forest <- agg_meta_sen_coef(subset(logr_c_jmax2,ecosystem.y=="forest")); jmax_forest$response <- "forest"
jmax_cropland <- agg_meta_sen_coef(subset(logr_c_jmax2,ecosystem.y=="cropland")); jmax_cropland$response <- "cropland"
jmax_alltypes <- dplyr::bind_rows(jmax_grassland,jmax_cropland,jmax_forest)
#remove a few sites with no available logr_var??? - but not applied since it not affects sen coef
final_jmax <- agg_meta_plots(jmax_alltypes,"response","logr","logr_var")
#created a unified one and help to combine to preediction column then 
jmax_alltypes2 <- jmax_alltypes; jmax_alltypes2$response <- "alltypes"

#for obs
final_observation_j <- agg_meta_plots(jmax_alltypes2,"response","logr","logr_var")
final_observation_j$type_name <- "all"

pred_plotmean_j <- mean(jmax_obs_pred_co2$pred_jmax)
final_prediction_j <- tibble(type_name="all",middle=pred_plotmean_j,ymin=pred_plotmean_j,ymax=pred_plotmean_j)
a2 <- final_jmax %>%
  ggplot( aes(x=type_name, y=middle)) + 
  geom_crossbar(aes(x=type_name, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5) +
  geom_crossbar(data=final_prediction_j,aes(x=type_name, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="red") +
  geom_crossbar(data=final_observation_j,aes(x=type_name, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="black") +
  geom_point(data=jmax_alltypes,aes(x=response, y=middle, size= 1.96/(middle-ymin)), alpha = 0.6, width = 0.5) +
  geom_point(data=jmax_obs_pred_co2,aes(x=response, y=pred_jmax), alpha = 0.6, width = 0.5,color="red") + #geom_point(data=subset(jmax_alltypes,is.na(middle_scaled)==TRUE),aes(x=response, y=middle,size=5,color="no_SE_info"),width = 0.5) +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(x="", y="Sensitivity coefficient of jmax",size=expression(paste("Standard Error"^{-1}))) +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))
a2

#finally, merged vcmax and jmax
final_vj <- merge(final_vcmax,final_jmax, by = c("type_name"),all.x=TRUE)
final_prediction_vj <- merge(final_prediction,final_prediction_j, by = c("type_name"),all.x=TRUE)
final_observation_vj <- merge(final_observation,final_observation_j, by = c("type_name"),all.x=TRUE)
vj_alltypes <- merge(vcmax_alltypes,jmax_alltypes, by = c("exp"),all.x=TRUE)
vj_obs_pred_co2 <- merge(vcmax_obs_pred_co2,jmax_obs_pred_co2, by = c("exp"),all.x=TRUE)

#another choice? aggregate directly?
final_observation_vj$middle<- mean(vj_alltypes$middle.y-vj_alltypes$middle.x,na.rm=TRUE)
final_observation_vj$ymin<- mean(vj_alltypes$middle.y-vj_alltypes$middle.x,na.rm=TRUE)
final_observation_vj$ymax<- mean(vj_alltypes$middle.y-vj_alltypes$middle.x,na.rm=TRUE)

a3 <- final_vj %>%
  ggplot( aes(x=type_name, y=middle.y-middle.x)) + 
  #geom_crossbar(aes(x=type_name, y=middle.y-middle.x, ymin=middle.y-middle.x, ymax=middle.y-middle.x), alpha = 0.6, width = 0.5) +
  geom_crossbar(data=final_prediction_vj,aes(x=type_name, y=middle.y-middle.x, ymin=middle.y-middle.x, ymax=middle.y-middle.x), alpha = 0.6, width = 0.5,color="red") +
  geom_crossbar(data=final_observation_vj,aes(x=type_name, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="black") +
  geom_point(data=vj_alltypes,aes(x=response.x, y=middle.y-middle.x), alpha = 0.6, width = 0.5) +
  geom_point(data=vj_obs_pred_co2,aes(x=response.x, y=pred_jmax-pred_vcmax), alpha = 0.6, width = 0.5,color="red") +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(x="", y="Sensitivity coefficient of jmax/vcmax") +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))
a3
#something wrong with cropland

plot_grid(a1,a2,a3,nrow=1,label_size = 15)+theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/egu_update_gcme.jpg",sep=""),width = 15, height = 5)

# now, newly adding smith's data

#new file
#GlycineCE -->ecocells_w (NOT ENTERED)
#OakOTC --> mi_c
#ORNL -->ornerp_liqui_c
#PineOTC --> Christchurch_pr_c
#RichmondGH --> australiacotton_c
# UIAPine --> ua_otc_c
#PHACE --> phace_c

#old file
#AspenFACE -->Rhine-aspenFACE_c (important!)
#BilyKriz --> Bily_Kriz_c
#Headley (remove Quercus	rubra and only keep) --> Headley_qp_c
#viesalm --> Vielsalm_c
  
smith_all <- read.csv("/Users/yunpeng/data/smith_keenan_gcb/gcb_co2/pred_vcmax.csv")

smith_all_removal <- subset(smith_all,exp_nam!="NZFACE" &exp_nam!="BioCON" & exp_nam!="ChinaRiceFACE" & exp_nam!="DukeFACE" & exp_nam!="EUROPOPFACE" & exp_nam!="NevadaFACE" & exp_nam!="SwissFACE")
smith_all_simple <- smith_all_removal[,c("SiteID","pft","sen_coef_v","sen_coef_j","pred_vcmax25_coef","pred_jmax25_coef","lon","lat","z")]

smith_all_plotmean <- aggregate(smith_all_simple,by=list(smith_all_simple$SiteID,smith_all_simple$pft), FUN=mean, na.rm=TRUE)[,c("Group.1","Group.2","sen_coef_v","sen_coef_j","pred_vcmax25_coef","pred_jmax25_coef","lon","lat","z")]
names(smith_all_plotmean) <- c("exp","ecosystem","vcmax","jmax","pred_vcmax","pred_jmax","lon","lat","z")
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
smith_combined <- smith_all_plotmean[,c("exp","ecosystem","vcmax","jmax")]

smith_all_plotmean$obs_jv <- smith_all_plotmean$jmax -  smith_all_plotmean$vcmax
smith_all_plotmean$pred_jv <- smith_all_plotmean$pred_jmax -  smith_all_plotmean$pred_vcmax
smith_all_plotmean$response <- "all"

final_vc <- tibble(ecosystem="all",obs_vcmax=mean(smith_all_plotmean$vcmax,na.rm=TRUE),
                   obs_jmax=mean(smith_all_plotmean$jmax,na.rm=TRUE),
                   pred_vcmax=mean(smith_all_plotmean$pred_vcmax,na.rm=TRUE),
                   pred_jmax=mean(smith_all_plotmean$pred_jmax,na.rm=TRUE))
#vcmax
b1 <- smith_all_plotmean %>%
  ggplot( aes(x=ecosystem, y=vcmax)) +
  geom_boxplot()+
  geom_point()+
  geom_boxplot(aes(x=response, y=vcmax))+
  geom_point(aes(x=response, y=vcmax))+
  geom_boxplot(aes(x=response, y=pred_vcmax),color="red")+
  geom_point(aes(x=response, y=pred_vcmax),color="red")+
  geom_hline( yintercept=0.0, size=0.5 )+ylim(-1,1)+
  labs(x="", y="Smith data") +theme_classic()+
  coord_flip() 
b1

#jmax
b2 <- smith_all_plotmean %>%
  ggplot( aes(x=ecosystem, y=jmax)) +
  geom_boxplot()+
  geom_point()+
  geom_boxplot(aes(x=response, y=jmax))+
  geom_point(aes(x=response, y=jmax))+
  geom_boxplot(aes(x=response, y=pred_jmax),color="red")+
  geom_point(aes(x=response, y=pred_jmax),color="red")+
  geom_hline( yintercept=0.0, size=0.5 )+ylim(-1,1)+
  labs(x="", y="Smith data") +theme_classic()+
  coord_flip() 
b2

#jmax/vcmax
b3 <- smith_all_plotmean %>%
  ggplot( aes(x=ecosystem, y=obs_jv)) +
  geom_boxplot()+
  geom_point()+
  geom_boxplot(aes(x=response, y=obs_jv))+
  geom_point(aes(x=response, y=obs_jv))+
  geom_boxplot(aes(x=response, y=pred_jv),color="red")+
  geom_point(aes(x=response, y=pred_jv),color="red")+
  geom_hline( yintercept=0.0, size=0.5 )+ylim(-1,1)+
  labs(x="", y="Smith data") +theme_classic()+
  coord_flip() 
b3


plot_grid(a1,a2,a3,b1,b2,b3,nrow=2,label_size = 15)+theme(plot.background=element_rect(fill="white", color="white"))

ggsave(paste("~/data/output_gcme/colin/egu_update_overall.jpg",sep=""),width = 15, height = 10)

#finally, combination of vcmax~jmax
vcmax_obs_pred_co2$ref <- "GCME"
smith_all_plotmean$ref <- "smith"
vcmax_final <- dplyr::bind_rows(vcmax_obs_pred_co2,smith_all_plotmean)

c1 <- ggplot(vcmax_obs_pred_co2,aes_string(y="vcmax", x="jmax")) +
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=ecosystem,shape=ref),size=3)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="black",method="lm",se=F)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))

c2 <- ggplot(smith_all_plotmean,aes_string(y="vcmax", x="jmax")) +
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=ecosystem,shape=ref),size=3)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="black",method="lm",se=F)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))

c3 <- ggplot(vcmax_final,aes_string(y="vcmax", x="jmax")) +
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=ecosystem,shape=ref),size=3)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="black",method="lm",se=F)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))

plot_grid(c1,c2,c3,nrow=1,label_size = 15)+theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/egu_update_overall2.jpg",sep=""),width = 15, height = 5)

#alternative method -simply using mean()
simple_mean <- function(logr_c_vcmax,logr_f_vcmax,logr_cf_vcmax,name){
  logr_c_vcmax$sen_coef <- logr_c_vcmax$logr/log(logr_c_vcmax$co2_e/logr_c_vcmax$co2_a)
  df_mean <- aggregate(logr_c_vcmax,by=list(logr_c_vcmax$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","sen_coef")]
  df_mean$condition <- "co2"
  logr_cf_vcmax$sen_coef <- logr_cf_vcmax$logr/log(logr_cf_vcmax$co2_e/logr_cf_vcmax$co2_a)
  df_mean2 <- aggregate(logr_cf_vcmax,by=list(logr_cf_vcmax$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","sen_coef","co2_a","co2_e")]
  df_mean3 <- aggregate(logr_f_vcmax,by=list(logr_f_vcmax$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","logr")]
  df_mean3$Group.1 <- paste(substr(df_mean3$Group.1,1,nchar(df_mean3$Group.1)-1),"cf",sep="")
  
  df_mean_additional <-Reduce(function(x,y) merge(x = x, y = y, by = c("Group.1"),all.x=TRUE),list(df_mean2,df_mean3))
  df_mean_additional$sen_coef <- df_mean_additional$sen_coef - df_mean_additional$logr/(log(df_mean_additional$co2_e-df_mean_additional$co2_a))
  df_mean_additional <- df_mean_additional[,c("Group.1","sen_coef")]
  df_mean_additional$condition <- "(co2 + Nfer)/Nfer"
  df_mean_final <- dplyr::bind_rows(df_mean,df_mean_additional)
  names(df_mean_final) <- c("exp",name,"condition")
  return(df_mean_final)
}

simple_mean_c <- function(logr_c_vcmax,name){
  logr_c_vcmax$sen_coef <- logr_c_vcmax$logr/log(logr_c_vcmax$co2_e/logr_c_vcmax$co2_a)
  df_mean <- aggregate(logr_c_vcmax,by=list(logr_c_vcmax$exp), FUN=mean, na.rm=TRUE)[,c("Group.1","sen_coef")]
  names(df_mean) <- c("exp",name)
  df_mean$condition<- "co2"
  return(df_mean)
}


vcmax_mean <- simple_mean(logr_c_vcmax,logr_f_vcmax,logr_cf_vcmax,"vcmax")
jmax_mean <-simple_mean(logr_c_jmax,logr_f_jmax,logr_cf_jmax,"jmax")
nmass_mean <-simple_mean(logr_c_nmass,logr_f_nmass,logr_cf_nmass,"nmass")
anpp_mean <-simple_mean(logr_c_anpp,logr_f_anpp,logr_cf_anpp,"anpp")
bnpp_mean <-simple_mean_c(old_logr_c_BNPP,"bnpp")
lai_mean <-simple_mean(logr_c_lai,logr_f_lai,logr_cf_lai,"lai")
soil_inorg_mean <-simple_mean(old_logr_c_ninorg,old_logr_f_ninorg,old_logr_cf_ninorg,"soil_inorg_N")
Asat_mean <-simple_mean(old_logr_c_Asat,old_logr_f_Asat,old_logr_cf_Asat,"Asat")
c13_mean <-simple_mean_c(old_logr_c_c13,"c13")


final_mean <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp","condition"),all.x=TRUE),
                    list(vcmax_mean,jmax_mean,nmass_mean,anpp_mean,bnpp_mean,lai_mean,soil_inorg_mean))

p <- list()
for(i in c(1:6)){
  p[[i]] <- ggplot(subset(final_mean),aes_string(x=names(final_mean)[i+3],
                                         y="vcmax")) +
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=condition),size=3)+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+labs(y="vcmax")+
    #geom_text(aes(label=exp),hjust=1, vjust=0,check_overlap = T)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))}

plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],nrow=2,label_size = 15)+
  theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/egu_update_alternative.jpg",sep=""),width = 15, height = 10)

#and update of soil N
logr_c_nmass
old_logr_c_ninorg$sen_coef <- old_logr_c_ninorg$logr/log(old_logr_c_ninorg$co2_e/old_logr_c_ninorg$co2_a)

mean_soilN <- aggregate(old_logr_c_ninorg,by=list(old_logr_c_ninorg$exp,old_logr_c_ninorg$Data_type), FUN=mean, na.rm=TRUE)[,c("Group.1","Group.2","sen_coef")]
names(mean_soilN) <- c("exp","variable","soil_inorg_N")
mean_soilN$condition<- "co2"
ninorg_main_v2 <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp","condition"),all.x=TRUE),
                        list(mean_soilN,vcmax_plot,jmax_plot,nmass_plot))
#remove biocon
#ninorg_main$old_ninorg[ninorg_main$old_ninorg< - 0.4] <- NA
p <- list()
for(i in c(1:3)){
  p[[i]] <- ggplot(ninorg_main_v2,aes_string(y=names(ninorg_main_v2)[i+4],
                                             x="soil_inorg_N")) +
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=variable),size=3)+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+labs(x="Inorganic soil N")+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))}

plot_grid(p[[1]],p[[2]],p[[3]],nrow=1,label_size = 15)
ggsave(paste("~/data/output_gcme/colin/egu_update_alternative_soilN.jpg",sep=""),width = 15, height = 5)


#combined GCME + Smith
smith_combined$condition <- "co2_smith"
smith_combined$exp <- tolower(smith_combined$exp)
names(smith_combined) <- c("exp","ecosystem_smith","vcmax_smith","jmax_smith","condition_smith")
smith_vcmax <- smith_combined[,c("exp","vcmax_smith","condition_smith")]
smith_jmax <- smith_combined[,c("exp","jmax_smith","ecosystem_smith")]
#vcmax_plot2 <- rbind(vcmax_plot,smith_vcmax)
#jmax_plot2 <- rbind(jmax_plot,smith_jmax)

vcmax_all <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp","condition"),all=TRUE),
                    list(vcmax_plot,jmax_plot,Asat_plot,anpp_plot,lma_plot,narea_plot,nmass_plot,
                         leaf_cn_plot,lai_plot,bnpp_plot,Nuptake_plot,npp_plot,soilN_plot,soil_total_N_plot,old_root_shoot_plot,old_ninorg_plot,gpp_plot))
vcmax_all$rep <- duplicated(vcmax_all$exp)
vcmax_all <- subset(vcmax_all,rep=="FALSE")
#vcmax_all_ci <- merge(vcmax_all,ci_final_removal,by=c("exp"),all.x=TRUE)
vcmax_all_ci <- merge(vcmax_all,ci_final,by=c("exp"),all.x=TRUE)

vcmax_all_ci_smith <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),list(vcmax_all_ci,smith_vcmax,smith_jmax))
vcmax_all_ci_smith$vcmax[is.na(vcmax_all_ci_smith$vcmax_smith)==FALSE] <- vcmax_all_ci_smith$vcmax_smith[is.na(vcmax_all_ci_smith$vcmax_smith)==FALSE]
vcmax_all_ci_smith$jmax[is.na(vcmax_all_ci_smith$jmax_smith)==FALSE] <- vcmax_all_ci_smith$jmax_smith[is.na(vcmax_all_ci_smith$jmax_smith)==FALSE]
vcmax_all_ci_smith_soil <- merge(vcmax_all_ci_smith,new_soil_vcmax[,c("exp","mineral_soilN")],by=c("exp"),all.x=TRUE)

vcmax_all_ci$x_increase <- vcmax_all_ci$elevated_ci/vcmax_all_ci$co2_e - vcmax_all_ci$ambient_ci/vcmax_all_ci$co2_a

ggplot(subset(vcmax_all_ci,data!="c13" & data!="ci"),aes_string(y="vcmax",x="x_increase")) +
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(size=3)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="black",method="lm",se=F)+
  geom_text(aes(label=substr(exp, 1, 11)),hjust=-0.1, vjust=0,check_overlap = F)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))

ggplot(subset(vcmax_all_ci, data=="ci"),aes_string(y="vcmax",x="x_increase")) +
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(size=3)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="black",method="lm",se=F)+
  geom_text(aes(label=substr(exp, 1, 11)),hjust=0, vjust=0,check_overlap = F)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))

ggplot(vcmax_all_ci,aes_string(y="vcmax",x="old_ninorg")) +
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(size=3)+xlab("Inorganic N")+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="black",method="lm",se=F)+
  geom_text(aes(label=substr(exp, 1, 11)),hjust=-0.1, vjust=0,check_overlap = F)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))

for(i in c(c(4,9,6,12,11,32))){
  p[[i]] <- ggplot(subset(vcmax_all_ci_smith_soil,condition=="co2"|condition=="(co2 + Nfer)/Nfer"),aes_string(x=names(vcmax_all_ci_smith_soil)[i],
                                             y="vcmax")) +
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=condition),size=3)+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))}
plot_grid(p[[4]],p[[9]],p[[6]],p[[12]],p[[11]],p[[32]],nrow=2,label_size = 15)
ggsave(paste("~/data/output_gcme/colin/egu_update_finalvcmax.jpg",sep=""),width = 25, height = 15)

#check lacked points

vcmax_all_ci_smith_plot <- merge(vcmax_all_ci_smith,unique(df_only[,c("exp","exp_nam")]),unique(kevin_othervars[,c("exp","site")]),by=c("exp"),all.x=TRUE)

vcmax_all_ci_smith_plot <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp"),all.x=TRUE),
                       list(vcmax_all_ci_smith,unique(df_only[,c("exp","exp_nam")]),unique(kevin_othervars[,c("exp","site")])))
vcmax_all_ci_smith_plot$exp_nam <- tolower(vcmax_all_ci_smith_plot$exp_nam)
vcmax_all_ci_smith_plot$site[is.na(vcmax_all_ci_smith_plot$site)==TRUE] <- vcmax_all_ci_smith_plot$exp_nam[is.na(vcmax_all_ci_smith_plot$site)==TRUE]
vcmax_all_ci_smith_plot$site[is.na(vcmax_all_ci_smith_plot$site)==TRUE] <- vcmax_all_ci_smith_plot$exp[is.na(vcmax_all_ci_smith_plot$site)==TRUE]
vcmax_all_ci_smith_plot[grep("riceface_japan", vcmax_all_ci_smith_plot$site),]$site <- "riceface_japan"
vcmax_all_ci_smith_plot[grep("riceface_china", vcmax_all_ci_smith_plot$site),]$site <- "riceface_china"

#vcmax, jmax, anpp, bnpp, npp, nmass, lai, soil N 

missing_vcmax <- unique(subset(vcmax_all_ci_smith_plot,condition=="co2"&is.na(vcmax)==TRUE & is.na(anpp)==FALSE)[,c("condition","site")])
missing_anpp <- unique(subset(vcmax_all_ci_smith_plot,condition=="co2"&is.na(vcmax)==FALSE & is.na(anpp)==TRUE)[,c("condition","site")])
missing_bnpp <- unique(subset(vcmax_all_ci_smith_plot,condition=="co2"&is.na(vcmax)==FALSE & is.na(bnpp)==TRUE)[,c("condition","site")])
missing_npp <- unique(subset(vcmax_all_ci_smith_plot,condition=="co2"&is.na(vcmax)==FALSE & is.na(npp)==TRUE)[,c("condition","site")])
missing_lai <- unique(subset(vcmax_all_ci_smith_plot,condition=="co2"&is.na(vcmax)==FALSE & is.na(lai)==TRUE)[,c("condition","site")])
missing_nmass <- unique(subset(vcmax_all_ci_smith_plot,condition=="co2"&is.na(vcmax)==FALSE & is.na(nmass)==TRUE)[,c("condition","site")])
missing_old_ninorg <- unique(subset(vcmax_all_ci_smith_plot,condition=="co2"&is.na(vcmax)==FALSE & is.na(old_ninorg)==TRUE)[,c("condition","site")])
missing_ci <- unique(subset(vcmax_all_ci_smith_plot,condition=="co2"&is.na(vcmax)==FALSE & is.na(ci)==TRUE)[,c("condition","site")])
names(missing_vcmax) <- c("vcmax","site");names(missing_anpp) <- c("anpp","site");names(missing_bnpp) <- c("bnpp","site");
names(missing_npp) <- c("npp","site");names(missing_lai) <- c("lai","site");names(missing_old_ninorg) <- c("old_ninorg","site");
names(missing_ci) <- c("ci","site");names(missing_nmass) <- c("nmass","site")
missing_together <-Reduce(function(x,y) merge(x = x, y = y, by = c("site"),all=TRUE),
                                 list(missing_vcmax,missing_anpp,missing_bnpp,missing_npp,missing_lai,missing_old_ninorg,missing_nmass,missing_ci))
missing_together_final <- data.frame(lapply(missing_together, gsub, pattern = "co2", replacement = "missing"))

#now, warming
# a look
aa <- kevin_z %>% group_by(treatment,exp)  %>% summarise(number = n())

#w, cw, dw
combine_warming <- function(logr_c_var,logr_w_var,logr_d_var,logr_cw_var,logr_dw_var,name){
  all_logr_w_lma <- agg_meta(logr_w_var)[,c("exp","middle")]
  all_logr_w_lma$condition <- "warming"
  
  #cw / c
  all_logr_c_lma <- agg_meta(logr_c_var)[,c("exp","middle")]
  all_logr_cw_lma <- agg_meta(logr_cw_var)[,c("exp","middle")]
  all_logr_cw_lma$exp_old <- all_logr_cw_lma$exp 
  all_logr_cw_lma$exp <-  paste(substr(all_logr_cw_lma$exp,1,nchar(all_logr_cw_lma$exp)-2),"c",sep="")
  all_logr_c_lma$exp;all_logr_cw_lma$exp
  all_logr_c_lma2 <- all_logr_c_lma %>% filter(exp %in%all_logr_cw_lma$exp)
  new_df <- merge(all_logr_c_lma2,all_logr_cw_lma,by=c("exp"),all.x=TRUE)
  new_df$middle <- (new_df$middle.y - new_df$middle.x)
  new_df$exp <- new_df$exp_old
  new_df$condition <- "(warming + co2)/co2"
  
  #dw / d
  all_logr_d_lma <- agg_meta(logr_d_var)[,c("exp","middle")]
  all_logr_dw_lma <- agg_meta(logr_dw_var)[,c("exp","middle")]
  all_logr_dw_lma$exp_old <- all_logr_dw_lma$exp 
  all_logr_dw_lma$exp <-  paste(substr(all_logr_dw_lma$exp,1,nchar(all_logr_dw_lma$exp)-2),"d",sep="")
  all_logr_d_lma$exp;all_logr_dw_lma$exp
  all_logr_d_lma2 <- all_logr_d_lma %>% filter(exp %in%all_logr_dw_lma$exp)
  new_df2 <- merge(all_logr_d_lma2,all_logr_dw_lma,by=c("exp"),all.x=TRUE)
  new_df2$middle <- (new_df2$middle.y - new_df2$middle.x)
  new_df2$exp <- new_df2$exp_old
  new_df2$condition <- "(warming + drought)/drought"
    
  lma_plot <- dplyr::bind_rows(all_logr_w_lma,new_df[,c("exp","middle","condition")],new_df2[,c("exp","middle","condition")]) 
  names(lma_plot) <- c("exp",name,"condition")
  return(lma_plot)
}

combine_warming_wd <- function(logr_c_var,logr_w_var,logr_d_var,logr_cw_var,logr_dw_var,name){
  all_logr_w_lma <- agg_meta(logr_w_var)[,c("exp","middle")]
  all_logr_w_lma$condition <- "warming"
  
  #dw / d
  all_logr_d_lma <- agg_meta(logr_d_var)[,c("exp","middle")]
  all_logr_dw_lma <- agg_meta(logr_dw_var)[,c("exp","middle")]
  all_logr_dw_lma$exp_old <- all_logr_dw_lma$exp 
  all_logr_dw_lma$exp <-  paste(substr(all_logr_dw_lma$exp,1,nchar(all_logr_dw_lma$exp)-2),"d",sep="")
  all_logr_d_lma$exp;all_logr_dw_lma$exp
  all_logr_d_lma2 <- all_logr_d_lma %>% filter(exp %in%all_logr_dw_lma$exp)
  new_df2 <- merge(all_logr_d_lma2,all_logr_dw_lma,by=c("exp"),all.x=TRUE)
  new_df2$middle <- (new_df2$middle.y - new_df2$middle.x)
  new_df2$exp <- new_df2$exp_old
  new_df2$condition <- "(warming + drought)/drought"
  
  lma_plot <- dplyr::bind_rows(all_logr_w_lma,new_df2[,c("exp","middle","condition")]) 
  names(lma_plot) <- c("exp",name,"condition")
  return(lma_plot)
}

combine_warming_wc <- function(logr_c_var,logr_w_var,logr_d_var,logr_cw_var,logr_dw_var,name){
  all_logr_w_lma <- agg_meta(logr_w_var)[,c("exp","middle")]
  all_logr_w_lma$condition <- "warming"
  
  #cw / c
  all_logr_c_lma <- agg_meta(logr_c_var)[,c("exp","middle")]
  all_logr_cw_lma <- agg_meta(logr_cw_var)[,c("exp","middle")]
  all_logr_cw_lma$exp_old <- all_logr_cw_lma$exp 
  all_logr_cw_lma$exp <-  paste(substr(all_logr_cw_lma$exp,1,nchar(all_logr_cw_lma$exp)-2),"c",sep="")
  all_logr_c_lma$exp;all_logr_cw_lma$exp
  all_logr_c_lma2 <- all_logr_c_lma %>% filter(exp %in%all_logr_cw_lma$exp)
  new_df <- merge(all_logr_c_lma2,all_logr_cw_lma,by=c("exp"),all.x=TRUE)
  new_df$middle <- (new_df$middle.y - new_df$middle.x)
  new_df$exp <- new_df$exp_old
  new_df$condition <- "(warming + co2)/co2"
  
  
  lma_plot <- dplyr::bind_rows(all_logr_w_lma,new_df[,c("exp","middle","condition")]) 
  names(lma_plot) <- c("exp",name,"condition")
  return(lma_plot)
}

combine_warming_w <- function(logr_c_var,logr_w_var,logr_d_var,logr_cw_var,logr_dw_var,name){
  all_logr_w_lma <- agg_meta(logr_w_var)[,c("exp","middle")]
  all_logr_w_lma$condition <- "warming"
  
  lma_plot <- all_logr_w_lma
  names(lma_plot) <- c("exp",name,"condition")
  return(lma_plot)
}

vcmax_plot_w <- combine_warming(logr_c_vcmax,logr_w_vcmax,logr_d_vcmax,logr_cw_vcmax,logr_dw_vcmax,"vcmax")
jmax_plot_w <- combine_warming(logr_c_jmax,logr_w_jmax,logr_d_jmax,logr_cw_jmax,logr_dw_jmax,"jmax")
lma_plot_w <- combine_warming(logr_c_LMA,logr_w_LMA,logr_d_LMA,logr_cw_LMA,logr_dw_LMA,"LMA")
narea_plot_w <- combine_warming(logr_c_narea,logr_w_narea,logr_d_narea,logr_cw_narea,logr_dw_narea,"narea")
nmass_plot_w <- combine_warming(logr_c_nmass,logr_w_nmass,logr_d_nmass,logr_cw_nmass,logr_dw_nmass,"nmass")
leaf_cn_plot_w <- combine_warming(logr_c_leaf_cn,logr_w_leaf_cn,logr_d_leaf_cn,logr_cw_leaf_cn,logr_dw_leaf_cn,"leaf_cn")
anpp_plot_w <- combine_warming_wd(logr_c_anpp,logr_w_anpp,logr_d_anpp,logr_cw_anpp,logr_dw_anpp,"anpp")
lai_plot_w <- combine_warming_wc(logr_c_lai,logr_w_lai,logr_d_lai,logr_cw_lai,logr_dw_lai,"lai")
bnpp_plot_w <- combine_warming_w(old_logr_c_BNPP,old_logr_w_BNPP,old_logr_d_BNPP,old_logr_cw_BNPP,old_logr_dw_BNPP,"bnpp")
Nuptake_plot_w <- combine_warming_w(old_logr_c_Nuptake,old_logr_w_Nuptake,old_logr_d_Nuptake,old_logr_cw_Nuptake,old_logr_dw_Nuptake,"Nuptake")
npp_plot_w <- combine_warming_w(old_logr_c_NPP,old_logr_w_NPP,old_logr_d_NPP,old_logr_cw_NPP,old_logr_dw_NPP,"npp")
soilN_plot_w <- combine_warming(logr_c_soil_n,logr_w_soil_n,logr_d_soil_n,logr_cw_soil_n,logr_dw_soil_n,"soilN")
soil_total_N_plot_w <- combine_warming_w(logr_c_soil_total_n,logr_w_soil_total_n,logr_d_soil_total_n,logr_cw_soil_total_n,logr_dw_soil_total_n,"soil_totalN")

vcmax_warming <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp","condition"),all.x=TRUE),
                    list(vcmax_plot_w,jmax_plot_w,anpp_plot_w,lma_plot_w,narea_plot_w,nmass_plot_w,
                         leaf_cn_plot_w,lai_plot_w,bnpp_plot_w,Nuptake_plot_w,npp_plot_w,soilN_plot_w,soil_total_N_plot_w))

for(i in c(1:(ncol(vcmax_warming)-4))){
  p[[i]] <- ggplot(vcmax_warming,aes_string(x=names(vcmax_warming)[i+4],
                                         y="vcmax")) +
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=condition),size=3)+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+geom_text(aes(label=exp),hjust=1, vjust=0,check_overlap = T)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))}

plot_grid(p[[4]],p[[3]],p[[2]],p[[5]],p[[8]],p[[10]],nrow=2,label_size = 15)+
  theme(plot.background=element_rect(fill="white", color="white"))

ggsave(paste("~/data/output_gcme/colin/egu_update_v_warming.jpg",sep=""),width = 20, height = 10)

for(i in c(1:(ncol(vcmax_warming)-4))){
  p[[i]] <- ggplot(vcmax_warming,aes_string(x=names(vcmax_warming)[i+4],
                                            y="jmax")) +
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=condition),size=3)+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+geom_text(aes(label=exp),hjust=1, vjust=0,check_overlap = T)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))}

plot_grid(p[[4]],p[[3]],p[[2]],p[[5]],p[[8]],p[[10]],nrow=2,label_size = 15)+
  theme(plot.background=element_rect(fill="white", color="white"))

ggsave(paste("~/data/output_gcme/colin/egu_update_j_warming.jpg",sep=""),width = 20, height = 10)


#now, fertilization
#f, cf, df

combine_fertilization <- function(logr_c_var,logr_f_var,logr_d_var,logr_cf_var,logr_df_var,name){
  all_logr_f_lma <- agg_meta(logr_f_var)[,c("exp","middle")]
  all_logr_f_lma$condition <- "fertilization"
  
  #cf / c
  all_logr_c_lma <- agg_meta(logr_c_var)[,c("exp","middle")]
  all_logr_cf_lma <- agg_meta(logr_cf_var)[,c("exp","middle")]
  all_logr_cf_lma$exp_old <- all_logr_cf_lma$exp 
  all_logr_cf_lma$exp <-  paste(substr(all_logr_cf_lma$exp,1,nchar(all_logr_cf_lma$exp)-2),"c",sep="")
  all_logr_c_lma$exp;all_logr_cf_lma$exp
  all_logr_c_lma2 <- all_logr_c_lma %>% filter(exp %in%all_logr_cf_lma$exp)
  new_df <- merge(all_logr_c_lma2,all_logr_cf_lma,by=c("exp"),all.x=TRUE)
  new_df$middle <- (new_df$middle.y - new_df$middle.x)
  new_df$exp <- new_df$exp_old
  new_df$condition <- "(fertilization + co2)/co2"
  
  #df / d
  all_logr_d_lma <- agg_meta(logr_d_var)[,c("exp","middle")]
  all_logr_df_lma <- agg_meta(logr_df_var)[,c("exp","middle")]
  all_logr_df_lma$exp_old <- all_logr_df_lma$exp 
  all_logr_df_lma$exp <-  paste(substr(all_logr_df_lma$exp,1,nchar(all_logr_df_lma$exp)-2),"d",sep="")
  all_logr_d_lma$exp;all_logr_df_lma$exp
  all_logr_d_lma2 <- all_logr_d_lma %>% filter(exp %in%all_logr_df_lma$exp)
  new_df2 <- merge(all_logr_d_lma2,all_logr_df_lma,by=c("exp"),all.x=TRUE)
  new_df2$middle <- (new_df2$middle.y - new_df2$middle.x)
  new_df2$exp <- new_df2$exp_old
  new_df2$condition <- "(fertilization + drought)/drought"
  
  lma_plot <- dplyr::bind_rows(all_logr_f_lma,new_df[,c("exp","middle","condition")],new_df2[,c("exp","middle","condition")]) 
  names(lma_plot) <- c("exp",name,"condition")
  return(lma_plot)
}
#above is not possible for either vcmax and jmax (not combination effect of drought + fer), so ignore them
combine_fertilization_fc <- function(logr_c_var,logr_f_var,logr_d_var,logr_cf_var,logr_df_var,name){
  all_logr_f_lma <- agg_meta(logr_f_var)[,c("exp","middle")]
  all_logr_f_lma$condition <- "fertilization"
  
  #cf / c
  all_logr_c_lma <- agg_meta(logr_c_var)[,c("exp","middle")]
  all_logr_cf_lma <- agg_meta(logr_cf_var)[,c("exp","middle")]
  all_logr_cf_lma$exp_old <- all_logr_cf_lma$exp 
  all_logr_cf_lma$exp <-  paste(substr(all_logr_cf_lma$exp,1,nchar(all_logr_cf_lma$exp)-2),"c",sep="")
  all_logr_c_lma$exp;all_logr_cf_lma$exp
  all_logr_c_lma2 <- all_logr_c_lma %>% filter(exp %in%all_logr_cf_lma$exp)
  new_df <- merge(all_logr_c_lma2,all_logr_cf_lma,by=c("exp"),all.x=TRUE)
  new_df$middle <- (new_df$middle.y - new_df$middle.x)
  new_df$exp <- new_df$exp_old
  new_df$condition <- "(fertilization + co2)/co2"

  lma_plot <- dplyr::bind_rows(all_logr_f_lma,new_df[,c("exp","middle","condition")]) 
  names(lma_plot) <- c("exp",name,"condition")
  return(lma_plot)
}

combine_fertilization_f <- function(logr_c_var,logr_f_var,logr_d_var,logr_cf_var,logr_df_var,name){
  all_logr_f_lma <- agg_meta(logr_f_var)[,c("exp","middle")]
  all_logr_f_lma$condition <- "fertilization"

  lma_plot <- all_logr_f_lma
  names(lma_plot) <- c("exp",name,"condition")
  return(lma_plot)
}

vcmax_plot_f <- combine_fertilization_fc(logr_c_vcmax,logr_f_vcmax,logr_d_vcmax,logr_cf_vcmax,logr_df_vcmax,"vcmax")
jmax_plot_f <- combine_fertilization_fc(logr_c_jmax,logr_f_jmax,logr_d_jmax,logr_cf_jmax,logr_df_jmax,"jmax")
lma_plot_f <- combine_fertilization_fc(logr_c_LMA,logr_f_LMA,logr_d_LMA,logr_cf_LMA,logr_df_LMA,"LMA")
narea_plot_f <- combine_fertilization_fc(logr_c_narea,logr_f_narea,logr_d_narea,logr_cf_narea,logr_df_narea,"narea")
nmass_plot_f <- combine_fertilization_fc(logr_c_nmass,logr_f_nmass,logr_d_nmass,logr_cf_nmass,logr_df_nmass,"nmass")
leaf_cn_plot_f <- combine_fertilization_fc(logr_c_leaf_cn,logr_f_leaf_cn,logr_d_leaf_cn,logr_cf_leaf_cn,logr_df_leaf_cn,"leaf_cn")
anpp_plot_f <- combine_fertilization_fc(logr_c_anpp,logr_f_anpp,logr_d_anpp,logr_cf_anpp,logr_df_anpp,"anpp")
lai_plot_f <- combine_fertilization_fc(logr_c_lai,logr_f_lai,logr_d_lai,logr_cf_lai,logr_df_lai,"lai")
bnpp_plot_f <- combine_fertilization_f(old_logr_c_BNPP,old_logr_f_BNPP,old_logr_d_BNPP,old_logr_cf_BNPP,old_logr_df_BNPP,"bnpp")
Nuptake_plot_f <- combine_fertilization_fc(old_logr_c_Nuptake,old_logr_f_Nuptake,old_logr_d_Nuptake,old_logr_cf_Nuptake,old_logr_df_Nuptake,"Nuptake")
npp_plot_f <- combine_fertilization_fc(old_logr_c_NPP,old_logr_f_NPP,old_logr_d_NPP,old_logr_cf_NPP,old_logr_df_NPP,"npp")
soilN_plot_f <- combine_fertilization_fc(logr_c_soil_n,logr_f_soil_n,logr_d_soil_n,logr_cf_soil_n,logr_df_soil_n,"soilN")
soil_total_N_plot_f <- combine_fertilization_fc(logr_c_soil_total_n,logr_f_soil_total_n,logr_d_soil_total_n,logr_cf_soil_total_n,logr_df_soil_total_n,"soil_totalN")

vcmax_fertilization <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp","condition"),all.x=TRUE),
                       list(vcmax_plot_f,jmax_plot_f,anpp_plot_f,lma_plot_f,narea_plot_f,nmass_plot_f,
                            leaf_cn_plot_f,lai_plot_f,bnpp_plot_f,Nuptake_plot_f,npp_plot_f,soilN_plot_f,soil_total_N_plot_f))

for(i in c(1:(ncol(vcmax_fertilization)-4))){
  p[[i]] <- ggplot(vcmax_fertilization,aes_string(x=names(vcmax_fertilization)[i+4],
                                            y="vcmax")) +
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=condition),size=3)+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))}

plot_grid(p[[4]],p[[3]],p[[2]],p[[1]],p[[6]],p[[11]],nrow=2,label_size = 15)+
  theme(plot.background=element_rect(fill="white", color="white"))

ggsave(paste("~/data/output_gcme/colin/egu_update_v_fer.jpg",sep=""),width = 20, height = 10)


for(i in c(1:(ncol(vcmax_fertilization)-4))){
  p[[i]] <- ggplot(vcmax_fertilization,aes_string(x=names(vcmax_fertilization)[i+4],
                                                  y="jmax")) +
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=condition),size=3)+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))}

plot_grid(p[[4]],p[[3]],p[[2]],p[[1]],p[[6]],p[[11]],nrow=2,label_size = 15)+
  theme(plot.background=element_rect(fill="white", color="white"))

ggsave(paste("~/data/output_gcme/colin/egu_update_j_fer.jpg",sep=""),width = 20, height = 10)


