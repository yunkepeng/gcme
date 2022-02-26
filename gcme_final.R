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
  #some are just all NA for original data's standard deviation - needs division -> in this way just calculate normal mean, lower and upper band
  for (i in 1:length(explist)){ # if over half of samples doens't have variance, or N =1 then we should not rely on this anyways.
    if (sum(is.na(subset(df,exp==explist[i])$logr_var),na.rm=TRUE)/length(is.na(subset(df,exp==explist[i])$logr_var)) >=threshold |
        length(subset(df,exp==explist[i])$logr_var) ==1){
      mean_value <- mean(subset(df,exp==explist[i])$logr,na.rm=TRUE)
      n_samples <- length(subset(df,exp==explist[i])$logr_var)
      df_box <- tibble(
        exp=explist[i], middle = mean_value,
        ymin   = NA,
        ymax   = NA,
        ymin_quarter = NA,
        ymax_quarter   = NA,
        variance_info ="No",no=n_samples,logr_var=NA)
    } else {
      #where year is random factor
      out_meta <-  subset(df,exp==explist[i]) %>% 
        mutate( id = 1:nrow(subset(df,exp==explist[i])) ) %>% 
        metafor::rma.mv( logr, logr_var, method = "REML", control = list(stepadj=0.3),  random = ~ 1 | exp_nam / id, data = .)
      
      out_meta_quarter <- subset(df,exp==explist[i]) %>% 
        mutate( id = 1:nrow(subset(df,exp==explist[i])) ) %>% 
        metafor::rma.mv( logr, logr_var, method = "REML",level=50, control = list(stepadj=0.3),  random = ~ 1 | exp_nam / id, data = .)
      
      n_samples <- length(subset(df,exp==explist[i])$logr_var)
      
      df_box <- tibble(
        exp=explist[i], middle = out_meta$b[1,1], ymin   = out_meta$ci.lb, ymax   = out_meta$ci.ub,
        ymin_quarter=out_meta_quarter$ci.lb, ymax_quarter = out_meta_quarter$ci.ub,
        variance_info ="Yes",no=n_samples,logr_var = (out_meta$se * sqrt(out_meta$k))^2)}
    mylist[[i]] <- df_box}
  output <- do.call("rbind",mylist)
  return(output)
}

#for response ratio
agg_meta <- function(df){
  threshold <- 0.5
  explist <- unique(df$exp)
  mylist <- list() #create an empty list
  #some are just all NA for original data's standard deviation - needs division -> in this way just calculate normal mean, lower and upper band
  for (i in 1:length(explist)){ # if over half of samples doens't have variance, or N =1 then we should not rely on this anyways.
    if (sum(is.na(subset(df,exp==explist[i])$logr_var),na.rm=TRUE)/length(is.na(subset(df,exp==explist[i])$logr_var)) >=threshold |
        length(subset(df,exp==explist[i])$logr_var) ==1){
      mean_value <- mean(subset(df,exp==explist[i])$logr,na.rm=TRUE)
      n_samples <- length(subset(df,exp==explist[i])$logr_var)
      df_box <- tibble(
        exp=explist[i], middle = mean_value,
        ymin   = NA,
        ymax   = NA,
        ymin_quarter = NA,
        ymax_quarter   = NA,
        variance_info ="No",no=n_samples,logr_var = NA)
    } else {
      #where year is random factor
      out_meta <-  subset(df,exp==explist[i]) %>% 
        mutate( id = 1:nrow(subset(df,exp==explist[i])) ) %>% 
        metafor::rma.mv( logr, logr_var, method = "REML", control = list(stepadj=0.3),  random = ~ 1 | exp_nam / id, data = .)
      
      out_meta_quarter <- subset(df,exp==explist[i]) %>% 
        mutate( id = 1:nrow(subset(df,exp==explist[i])) ) %>% 
        metafor::rma.mv( logr, logr_var, method = "REML",level=50, control = list(stepadj=0.3),  random = ~ 1 | exp_nam / id, data = .)
      
      n_samples <- length(subset(df,exp==explist[i])$logr_var)
      
      df_box <- tibble(
        exp=explist[i], middle = out_meta$b[1,1], ymin   = out_meta$ci.lb, ymax   = out_meta$ci.ub,
        ymin_quarter=out_meta_quarter$ci.lb, ymax_quarter = out_meta_quarter$ci.ub,
        variance_info ="Yes",no=n_samples,logr_var = (out_meta$se * sqrt(out_meta$k))^2)}
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
logr_f_LMA <-  as_tibble(response_ratio_v2(subset(LMA,treatment=="f")));
logr_cf_LMA <-  as_tibble(response_ratio_v2(subset(LMA,treatment=="cf")))

#anpp, bnpp, nmass, LAI, soil N
kevin_othervars <- read.csv("/Users/yunpeng/data/gcme/kevin_20220222/MESI_2022.csv")
kevin_othervars <- rename(kevin_othervars, c(ambient = x_c, elevated=x_t, ambient_Sd=sd_c, elevated_Sd=sd_t,ambient_Se=se_c,elevated_Se=se_t,n_plots=rep_c,
                                             z=elevation, co2_a=c_c, co2_e=c_t, nfertQ_a = n_c, nfertQ_e = n_t, pfertQ_a = p_c, pfertQ_e = p_t,kfertQ_a = k_c, kfertQ_e = k_t,
                                             warmQ_e1 = w_t1, warmQ_e2 = w_t2, warmQ_e3 = w_t3, Unit=x_units))
kevin_othervars$ambient <-as.numeric(kevin_othervars$ambient)
kevin_othervars$elevated <-as.numeric(kevin_othervars$elevated)
kevin_othervars$ambient_Sd  <-as.numeric(kevin_othervars$ambient_Sd)
kevin_othervars$elevated_Sd  <-as.numeric(kevin_othervars$elevated_Sd)
kevin_othervars$ambient_Se <- as.numeric(kevin_othervars$ambient_Se)
kevin_othervars$elevated_Se <- as.numeric(kevin_othervars$elevated_Se)
kevin_othervars$n_plots  <-as.numeric(kevin_othervars$n_plots)
kevin_othervars$z <- as.numeric(kevin_othervars$z)

#convert lai_max to lai
kevin_othervars$response[kevin_othervars$response=="lai_max"] <- "lai"

#check numbers of sites and variables 
sitename <- kevin_othervars%>% group_by(site)  %>% summarise(number = n())
varname <- kevin_othervars%>% group_by(response)  %>% summarise(number = n())

#leaf_n needs division - but leaf_p doesn't need (as it is all pmass)
nmass <- subset(kevin_othervars,response=="leaf_n") %>%         
  filter(Unit %in% c("g","g_100g","g_g","g_kg","mg_g","mg_kg","mmol_g","ug_mg","umol_g"))
nmass$response <- "nmass"
logr_c_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="c")));
logr_f_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="f")));
logr_cf_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="cf")))

narea <- subset(kevin_othervars,response=="leaf_n") %>%         
  filter(Unit %in% c("g_m2","g_m3","g_pot","mg_cm2","mg_m2","mmol_m2","ug_cm2","ug_cm3","umol_m2"))
narea$response <- "narea"

logr_c_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="c")));
logr_f_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="f")));
logr_cf_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="cf")))

#other vars - all created now
for (i in 1:nrow(varname)) {
  tryCatch({
    varname1 <- varname$response[i]
    df_c <- subset(kevin_othervars,treatment=="c" & response==varname1)
    assign(paste("logr_c_", varname1,sep=""), as_tibble(response_ratio_v2(df_c)))
    
    df_f <- subset(kevin_othervars,treatment=="f"& response==varname1)
    assign(paste("logr_f_", varname1,sep=""), as_tibble(response_ratio_v2(df_f)))
    
    df_cf <- subset(kevin_othervars,treatment=="cf" & response==varname1)
    assign(paste("logr_cf_", varname1,sep=""), as_tibble(response_ratio_v2(df_cf)))

  }, error=function(e){})} 


#check soil data
kevin_othervars[grep("soil", kevin_othervars$response),]

aa <- kevin_othervars[grep("nh4", kevin_othervars$response),]%>% group_by(response,site,Unit)  %>% summarise(number = n())
