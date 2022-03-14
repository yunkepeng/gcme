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
library(Hmisc)

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

combine_co2_c <- function(logr_c_var,logr_f_var,logr_cf_var,name){
  all_logr_c_lma <- agg_meta_sen_coef(logr_c_var)[,c("exp","middle")]
  all_logr_c_lma$condition <- "co2"
  lma_plot <- all_logr_c_lma
  names(lma_plot) <- c("exp",name,"condition")
  return(lma_plot)
}


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
kevin_othervars$exp_nam <- kevin_othervars$site


#filter plots
vcmax_plot <- read.csv("/Users/yunpeng/data/gcme/kevin_20220222/vcmax_plotmean.csv")$exp
jmax_plot <- read.csv("/Users/yunpeng/data/gcme/kevin_20220222/jmax_plotmean.csv")$exp
photo_plot <- unique(c(vcmax_plot,jmax_plot,c("duke2_f","euroface4_pa_cf","euroface4_pe_cf","euroface4_pn_cf","biforface_c")))
kevin_othervars_cf <- kevin_othervars %>% filter(exp %in% photo_plot)
#convert plot properly 
kevin_othervars_cf$exp_nam[is.na(kevin_othervars_cf$exp_nam)==TRUE] <- substr(subset(kevin_othervars_cf,is.na(exp_nam)==TRUE)$exp,1,nchar(subset(kevin_othervars_cf,is.na(exp_nam)==TRUE)$exp)-2)


#have a look
aa <- kevin_othervars_cf%>% group_by(response,Unit,exp)  %>% summarise(number = n())

#leaf_n needs division - but leaf_p doesn't need (as it is all pmass)
nmass <- subset(kevin_othervars_cf,response=="leaf_n") %>%         
  filter(Unit %in% c("g","g_100g","g_g","g_kg","mg_g","mg_kg","mmol_g","ug_mg","umol_g"))
nmass$response <- "nmass"
logr_c_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="c")));
logr_f_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="f")));
logr_cf_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="cf")))

narea <- subset(kevin_othervars_cf,response=="leaf_n") %>%         
  filter(Unit %in% c("g_m2","g_m3","g_pot","mg_cm2","mg_m2","mmol_m2","ug_cm2","ug_cm3","umol_m2"))
narea$response <- "narea"

logr_c_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="c")));
logr_f_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="f")));
logr_cf_narea <- as_tibble(response_ratio_v2(subset(narea,treatment=="cf")))

#anpp and anpp_grain are both anpp.
kevin_othervars_cf$response[kevin_othervars_cf$response=="anpp_grain"] <- "anpp"
unique(subset(kevin_othervars_cf,response=="anpp")$Unit)

#bnpp - make sure data between them are not mixed.
kevin_othervars_cf$output[kevin_othervars_cf$response=="coarse_root_production"|
                            kevin_othervars_cf$response=="fine_root_turnover"|
                            kevin_othervars_cf$response=="root_production"|
                            (kevin_othervars_cf$response=="fine_root_biomass" & kevin_othervars_cf$Unit=="g_m2y")|
                            (kevin_othervars_cf$response=="fine_root_production" & (kevin_othervars_cf$Unit=="gc_m2yr"|kevin_othervars_cf$Unit=="mm_cm2d"|kevin_othervars_cf$Unit=="g_m2y"| kevin_othervars_cf$Unit=="gc_m2y"|kevin_othervars_cf$Unit=="gdw_m2y"))] <- "bnpp"

bnpp_dataset <- subset(kevin_othervars_cf,output=="bnpp"&(response=="root_production"|exp=="christchurch_pr_c"|exp=="nevada_desert_face_c"|exp=="popface_pa_c"|exp=="popface_pe_c"|exp=="popface_pn_c"| (exp=="facts_ii_face3_pt_c"&response=="fine_root_production")|exp=="biforface_c"))

unique(bnpp_dataset[,c("exp","response","Unit")])

#soil mineral N in dry-mass: 7+1
#7 plots: nh4 + no3
soil_nh4 <- subset(kevin_othervars_cf,response=="soil_nh4-n") %>% group_by(exp,Unit) %>% summarise(co2_a = mean(co2_a), co2_e = mean(co2_e), ambient = mean(ambient), elevated = mean(elevated))
soil_no3 <- subset(kevin_othervars_cf,response=="soil_no3-n") %>% group_by(exp,Unit) %>% summarise(co2_a = mean(co2_a), co2_e = mean(co2_e),ambient = mean(ambient), elevated = mean(elevated))
soil_nh4no3 <- na.omit(merge(soil_nh4,soil_no3,by=c("exp","Unit"),all.x=TRUE))
soil_nh4no3$logr <- log((soil_nh4no3$elevated.x+soil_nh4no3$elevated.y)/(soil_nh4no3$ambient.x+soil_nh4no3$ambient.y))
soil_nh4no3$soil_mineral_N <- log((soil_nh4no3$elevated.x+soil_nh4no3$elevated.y)/(soil_nh4no3$ambient.x+soil_nh4no3$ambient.y))/log(soil_nh4no3$co2_e.x/soil_nh4no3$co2_a.x)

duke2_cf_soil <- (soil_nh4no3$logr[soil_nh4no3$exp=="duke2_cf"]-soil_nh4no3$logr[soil_nh4no3$exp=="duke2_f"])/log(soil_nh4no3$co2_e.x[soil_nh4no3$exp=="duke2_cf"]/soil_nh4no3$co2_a.x[soil_nh4no3$exp=="duke2_cf"])
soil_mineral_1 <- subset(soil_nh4no3,exp!="duke2_cf" & exp!="duke2_f")[,c("exp","soil_mineral_N")]
soil_mineral_2 <- tibble(exp="duke2_cf",soil_mineral_N =duke2_cf_soil)

#1 plot: soil_in
soil_mineral_3 <- tibble(exp="facts_ii_face3_pt_c",soil_mineral_N =agg_meta_sen_coef(response_ratio_v2(subset(kevin_othervars_cf,response=="soil_in")))$middle)
soil_mineral_dry <- rbind(soil_mineral_1,soil_mineral_2,soil_mineral_3)
soil_mineral_dry$type <- "dry"

#soil mineral N in solution: 2+1
solution_mineral <-  agg_meta_sen_coef(response_ratio_v2(subset(kevin_othervars_cf,response=="soil_solution_mineral_n")))
soil_mineral_4 <- tibble(exp=solution_mineral$exp,soil_mineral_N=solution_mineral$middle)

soil_solution_nh4 <- subset(kevin_othervars_cf,response=="soil_solution_nh4") %>% group_by(exp,Unit) %>% summarise(co2_a = mean(co2_a), co2_e = mean(co2_e), ambient = mean(ambient), elevated = mean(elevated))
soil_solution_no3 <- subset(kevin_othervars_cf,response=="soil_solution_no3") %>% group_by(exp,Unit) %>% summarise(co2_a = mean(co2_a), co2_e = mean(co2_e),ambient = mean(ambient), elevated = mean(elevated))
soil_solution_nh4no3 <- merge(soil_solution_nh4,soil_solution_no3,by=c("exp","Unit"),all.x=TRUE)
soil_solution_nh4no3$soil_mineral_N <- log((soil_solution_nh4no3$elevated.x+soil_solution_nh4no3$elevated.y)/(soil_solution_nh4no3$ambient.x+soil_solution_nh4no3$ambient.y))/log(soil_solution_nh4no3$co2_e.x/soil_solution_nh4no3$co2_a.x)
soil_mineral_5 <- tibble(exp=soil_solution_nh4no3$exp,soil_mineral_N=soil_solution_nh4no3$soil_mineral_N)

soil_mineral_wet <- rbind(soil_mineral_4,soil_mineral_5)
soil_mineral_wet$type <- "wet"

soil_mineral <- rbind(soil_mineral_dry,soil_mineral_wet)
soil_mineral_plotmean <- soil_mineral%>% group_by(exp,type)  %>% summarise(soil_mineral_N = mean(soil_mineral_N))

#others are all the same
subset(kevin_othervars_cf,response=="lai")
subset(kevin_othervars_cf,response=="gpp") 
subset(kevin_othervars_cf,response=="soil_total_n")
subset(kevin_othervars_cf,response=="soil_cn_min_layer")
subset(kevin_othervars_cf,response=="agb") # Unit needs identified
subset(kevin_othervars_cf,response=="agb_n") # Unit needs identified
subset(kevin_othervars_cf,response=="asat") 
subset(kevin_othervars_cf,response=="anet")
subset(kevin_othervars_cf,response=="bgb") # Unit needs identified
subset(kevin_othervars_cf,response=="bgb_coarse") # Unit needs identified
subset(kevin_othervars_cf,response=="bgb_n") # Unit needs identified
subset(kevin_othervars_cf,response=="fine_root_biomass" & Unit !="g_m2y") # Unit needs identified
subset(kevin_othervars_cf,response=="fine_root_n") # Unit needs identified
subset(kevin_othervars_cf,response=="leaf_biomass") 
subset(kevin_othervars_cf,response=="leaf_c") 
subset(kevin_othervars_cf,response=="leaf_cn") 
subset(kevin_othervars_cf,response=="leaf_litter_c") 
subset(kevin_othervars_cf,response=="leaf_litter_cn") 
subset(kevin_othervars_cf,response=="leaf_litter_n") 
subset(kevin_othervars_cf,response=="leaf_nue")
subset(kevin_othervars_cf,response=="litter_biomass")
subset(kevin_othervars_cf,response=="litterfall")
subset(kevin_othervars_cf,response=="mb")
subset(kevin_othervars_cf,response=="mbc")
subset(kevin_othervars_cf,response=="mbn")
subset(kevin_othervars_cf,response=="r_leaf")
subset(kevin_othervars_cf,response=="r_root")
subset(kevin_othervars_cf,response=="r_soil")
subset(kevin_othervars_cf,response=="root_n_uptake")
subset(kevin_othervars_cf,response=="root_shoot_ratio")
subset(kevin_othervars_cf,response=="soc")
subset(kevin_othervars_cf,response=="soil_c")
subset(kevin_othervars_cf,response=="soil_total_c") # soil_total_c and soil_c?
subset(kevin_othervars_cf,response=="total_biomass")
subset(kevin_othervars_cf,response=="total_biomass_n")
subset(kevin_othervars_cf,response=="wood_n")
subset(kevin_othervars_cf,response=="soil_potential_net_n_mineralization")


#check numbers of sites and variables 
varname <- kevin_othervars_cf%>% group_by(response)  %>% summarise(number = n())

#other vars - all created now
for (i in 1:nrow(varname)) {
  tryCatch({
    varname1 <- varname$response[i]
    df_c <- subset(kevin_othervars_cf,treatment=="c" & response==varname1)
    assign(paste("logr_c_", varname1,sep=""), as_tibble(response_ratio_v2(df_c)))
    
    df_f <- subset(kevin_othervars_cf,treatment=="f"& response==varname1)
    assign(paste("logr_f_", varname1,sep=""), as_tibble(response_ratio_v2(df_f)))
    
    df_cf <- subset(kevin_othervars_cf,treatment=="cf" & response==varname1)
    assign(paste("logr_cf_", varname1,sep=""), as_tibble(response_ratio_v2(df_cf)))
  }, error=function(e){})} 

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
kevin_LMA <- kevin_LMA %>% filter(exp %in% photo_plot)


#THIS dangerous! However, since we don't know year of measurement, we could only assume it as average
kevin_LMA$Year[is.na(kevin_LMA$Year)==TRUE] <- 2
sla <- subset(kevin_LMA,response=="sla")
lma <- subset(kevin_LMA,response=="lma")
#convert SLA to LMA format 
sla$ambient <- 1/sla$ambient; sla$elevated <- 1/sla$elevated; sla$ambient_Sd <- 1/sla$ambient_Sd
sla$elevated_Sd <- 1/sla$elevated_Sd;sla$ambient_Se <- 1/sla$ambient_Se;sla$elevated_Se <- 1/sla$elevated_Se
LMA <- dplyr::bind_rows(lma,sla)
LMA$response <- "LMA"
LMA2 <- merge(LMA,unique(kevin_othervars_cf[,c("exp","exp_nam")]),by=c("exp"),all.x=TRUE)
logr_c_LMA <- as_tibble(response_ratio_v2(subset(LMA2,treatment=="c")));
logr_f_LMA <-  as_tibble(response_ratio_v2(subset(LMA2,treatment=="f")));
logr_cf_LMA <-  as_tibble(response_ratio_v2(subset(LMA2,treatment=="cf")))


#######final aggregation 
lma_plot <- combine_co2_cf(logr_c_LMA,logr_f_LMA,logr_cf_LMA,"LMA")
narea_plot <- combine_co2_cf(logr_c_narea,logr_f_narea,logr_cf_narea,"narea")
nmass_plot <- combine_co2_cf(logr_c_nmass,logr_f_nmass,logr_cf_nmass,"nmass")

logr_c_bnpp <- response_ratio_v2(bnpp_dataset)
bnpp_plot <- combine_co2_c(logr_c_bnpp,logr_c_bnpp,logr_c_bnpp,"bnpp")
bnpp_plot$bnpp[bnpp_plot$exp=="duke_c"] <- log(subset(bnpp_dataset,exp=="duke_c")$elevated/subset(bnpp_dataset,exp=="duke_c")$ambient)/log(563/363)
bnpp_plot_final <- merge(bnpp_plot,unique(bnpp_dataset[,c("exp","response")]),by=c("exp"),all.x=TRUE)

soil_mineral_plotmean$condition <- "co2"
soil_mineral_plotmean$condition[soil_mineral_plotmean$exp=="duke2_cf"] <- "fertilization"
soil_mineral_plotmean <- soil_mineral_plotmean[,c("exp","condition","soil_mineral_N","type")]
soil_mineral_plotmean_dry <- subset(soil_mineral_plotmean,type=="dry")

anpp_plot <- combine_co2_cf(logr_c_anpp,logr_f_anpp,logr_cf_anpp,"anpp")
lai_plot <- combine_co2_cf(logr_c_lai,logr_f_lai,logr_cf_lai,"lai")
soil_total_n_plot <- combine_co2_cf(logr_c_soil_total_n,logr_f_soil_total_n,logr_cf_soil_total_n,"soil_total_n")
soil_cn_min_layer_plot <- combine_co2_c(logr_c_soil_cn_min_layer,logr_f_soil_cn_min_layer,logr_cf_soil_cn_min_layer,"soil_cn_min_layer")
agb_plot <- combine_co2_cf(logr_c_agb,logr_f_agb,logr_cf_agb,"agb")
agb_n_plot <- combine_co2_c(logr_c_agb_n,logr_f_agb_n,logr_cf_agb_n,"agb_n")
asat_plot <- combine_co2_c(logr_c_asat,logr_f_asat,logr_cf_asat,"asat")
anet_plot <- combine_co2_c(logr_c_anet,logr_f_anet,logr_cf_anet,"anet")
bgb_plot <- combine_co2_c(logr_c_bgb,logr_f_bgb,logr_cf_bgb,"bgb")
bgb_coarse_plot <- combine_co2_cf(logr_c_bgb_coarse,logr_f_bgb_coarse,logr_cf_bgb_coarse,"bgb_coarse")
bgb_n_plot <- combine_co2_cf(logr_c_bgb_n,logr_f_bgb_n,logr_cf_bgb_n,"bgb_n")
fine_root_biomass_plot <- combine_co2_cf(logr_c_fine_root_biomass,logr_f_fine_root_biomass,logr_cf_fine_root_biomass,"fine_root_biomass")
fine_root_n_plot <- combine_co2_c(logr_c_fine_root_n,logr_f_fine_root_n,logr_cf_fine_root_n,"fine_root_n")
leaf_biomass_plot <- combine_co2_c(logr_c_leaf_biomass,logr_f_leaf_biomass,logr_cf_leaf_biomass,"leaf_biomass")
leaf_c_plot <- combine_co2_c(logr_c_leaf_c,logr_f_leaf_c,logr_cf_leaf_c,"leaf_c")
leaf_cn_plot <- combine_co2_c(logr_c_leaf_cn,logr_f_leaf_cn,logr_cf_leaf_cn,"leaf_cn")
leaf_litter_c_plot <- combine_co2_c(logr_c_leaf_litter_c,logr_f_leaf_litter_c,logr_cf_leaf_litter_c,"leaf_litter_c")
leaf_litter_cn_plot <- combine_co2_c(logr_c_leaf_litter_cn,logr_f_leaf_litter_cn,logr_cf_leaf_litter_cn,"leaf_litter_cn")
leaf_litter_n_plot <- combine_co2_c(logr_c_leaf_litter_n,logr_f_leaf_litter_n,logr_cf_leaf_litter_n,"leaf_litter_n")
leaf_nue_plot <- combine_co2_c(logr_c_leaf_nue,logr_f_leaf_nue,logr_cf_leaf_nue,"leaf_nue")
litter_biomass_plot <- combine_co2_c(logr_c_litter_biomass,logr_f_litter_biomass,logr_cf_litter_biomass,"litter_biomass")
litterfall_plot <- combine_co2_c(logr_c_litterfall,logr_f_litterfall,logr_cf_litterfall,"litterfall")
mb_plot <- combine_co2_cf(logr_c_mb,logr_f_mb,logr_cf_mb,"mb")
mbc_plot <- combine_co2_c(logr_c_mbc,logr_f_mbc,logr_cf_mbc,"mbc")
mbn_plot <- combine_co2_c(logr_c_mbn,logr_f_mbn,logr_cf_mbn,"mbn")
r_leaf_plot <- combine_co2_c(logr_c_r_leaf,logr_f_r_leaf,logr_cf_r_leaf,"r_leaf")
r_root_plot <- combine_co2_c(logr_c_r_root,logr_f_r_root,logr_cf_r_root,"r_root")
r_soil_plot <- combine_co2_cf(logr_c_r_soil,logr_f_r_soil,logr_cf_r_soil,"r_soil")
root_n_uptake_plot <- combine_co2_c(logr_c_root_n_uptake,logr_f_root_n_uptake,logr_cf_root_n_uptake,"root_n_uptake")
root_shoot_ratio_plot <- combine_co2_c(logr_c_root_shoot_ratio,logr_f_root_shoot_ratio,logr_cf_root_shoot_ratio,"root_shoot_ratio")
soc_plot <- combine_co2_cf(logr_c_soc,logr_f_soc,logr_cf_soc,"soc")
soil_c_plot <- combine_co2_c(logr_c_soil_c,logr_f_soil_c,logr_cf_soil_c,"soil_c")
soil_total_c_plot <- combine_co2_c(logr_c_soil_total_c,logr_f_soil_total_c,logr_cf_soil_total_c,"soil_total_c")
total_biomass_plot <- combine_co2_c(logr_c_total_biomass,logr_f_total_biomass,logr_cf_total_biomass,"total_biomass")
total_biomass_n_plot <- combine_co2_c(logr_c_total_biomass_n,logr_f_total_biomass_n,logr_cf_total_biomass_n,"total_biomass_n")
wood_n_plot <- combine_co2_c(logr_c_wood_n,logr_f_wood_n,logr_cf_wood_n,"wood_n")
soil_potential_net_n_mineralization_plot <- combine_co2_c(logr_c_soil_potential_net_n_mineralization,logr_f_soil_potential_net_n_mineralization,logr_cf_soil_potential_net_n_mineralization,"soil_potential_net_n_mineralization")

vcmax_all_plot <- read.csv("/Users/yunpeng/data/gcme/kevin_20220222/vcmax_plotmean.csv")[,c("ecosystem","exp","middle","condition","variance_info")]
vcmax_all_plot$condition[is.na(vcmax_all_plot$condition)==TRUE] <- "co2"
names(vcmax_all_plot) <- c("ecosystem","exp","vcmax","condition","variance_info")
jmax_all_plot <- read.csv("/Users/yunpeng/data/gcme/kevin_20220222/jmax_plotmean.csv")[,c("exp","middle","condition")]
jmax_all_plot$condition[is.na(jmax_all_plot$condition)==TRUE] <- "co2"
names(jmax_all_plot) <- c("exp","jmax","condition")

final_mean <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp","condition"),all.x=TRUE),
                    list(vcmax_all_plot,jmax_all_plot,lma_plot,narea_plot,nmass_plot,bnpp_plot_final,soil_mineral_plotmean_dry,
                         anpp_plot,lai_plot,soil_total_n_plot,soil_cn_min_layer_plot,agb_plot,agb_n_plot,asat_plot,anet_plot,
                         bgb_plot,bgb_coarse_plot,bgb_n_plot,fine_root_biomass_plot,fine_root_n_plot,leaf_biomass_plot,leaf_c_plot,leaf_cn_plot,
                         leaf_litter_c_plot,leaf_litter_cn_plot,leaf_litter_n_plot,leaf_nue_plot,litter_biomass_plot,litterfall_plot,
                         mb_plot,mbc_plot,mbn_plot,r_leaf_plot,r_root_plot,r_soil_plot,root_n_uptake_plot,root_shoot_ratio_plot,soc_plot,soil_c_plot,
                         soil_total_c_plot,total_biomass_plot,total_biomass_n_plot,wood_n_plot,soil_potential_net_n_mineralization_plot))

#variables interested...
c("soil_nh4-n","soil_no3-n","soil_in")
c("soil_solution_mineral_n","soil_solution_nh4","soil_solution_no3")


soil_min <- kevin_othervars_cf %>%
  filter(response %in%c("soil_nh4-n","soil_no3-n","soil_in"))
soil_min_c <- subset(soil_min,treatment=="c")
new_soil_min <- agg_meta_sen_coef(response_ratio_v2(soil_min_c))
test<- (merge(final_mean,new_soil_min[,c("exp","middle")],by=c("exp"),all.x=TRUE))

ggplot(test,aes_string(x="middle", y="vcmax")) +
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=condition),size=3)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="black",method="lm",se=F)+labs(x="mineral soil N")+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))

#further process

#1. Remove a negative Asat (since it is Asat at constant co2)
final_mean$asat[final_mean$exp=="new_zealand_face_c"] <- NA
#2. remove high Tleaf measurement's vcmax (actually, two papers: darbah_et_al_2010a and darbah_et_al_2010b)
final_mean$vcmax[final_mean$exp=="facts_ii_face3_pt_c"] <- agg_meta_sen_coef(response_ratio_v2(subset(logr_c_vcmax,exp=="facts_ii_face3_pt_c" & citation!="darbah_et_al_2010b" & citation!="darbah_et_al_2010a")))$middle
final_mean$vcmax[final_mean$exp=="facts_ii_face4_bp_c"] <- agg_meta_sen_coef(response_ratio_v2(subset(logr_c_vcmax,exp=="facts_ii_face4_bp_c" & citation!="darbah_et_al_2010b")))$middle

#3. anpp modification
# when removing g/m2 - it becomes better
anpp_new <- subset(kevin_othervars_cf,response=="anpp"&Unit!="t_ha"&Unit!="g_m2"&Unit!="gc_m2"&Unit!="mg"&Unit!="g_plant")
unique(anpp_new[,c("exp","Unit")]);unique(anpp_new[,c("Unit")])

anpp_new_c <- response_ratio_v2(subset(anpp_new,treatment=="c"))
#anpp_new_f <- response_ratio_v2(subset(anpp_new,treatment=="f"))
anpp_new_cf <- response_ratio_v2(subset(anpp_new,treatment=="cf"))

anpp_new_vcmax <- merge(vcmax_all_plot,combine_co2_c(anpp_new_c,anpp_new_f,anpp_new_cf,"anpp_new"),by = c("exp","condition"),all.x=TRUE)

#Add additional grassland anpp - summing up samples within year and include
unique(logr_c_anpp$Unit)
aa <- subset(logr_c_anpp,Unit=="g_m2"|Unit=="gc_m2"|Unit=="t_ha")
unique(aa[,c("exp","start_year","sampling_year")])
# 5 additional plots available as the others are forest --> they are just one sample, so use it directly
add_anpp <- final_mean[,c("exp","anpp")] %>% filter(exp %in%c("riceface_japan_ko_2013_3558_13960_c","riceface_japan_a_2003_3938_14057_c","riceface_japan_a_2004_3938_14057_c","soyfacesoy1_c","soyfacesoy2_c"))
#only 3 samples - include them
anpp_new_vcmax$anpp_new[anpp_new_vcmax$exp=="riceface_japan_ko_2013_3558_13960_c"] <- add_anpp$anpp[add_anpp$exp=="riceface_japan_ko_2013_3558_13960_c"] 
anpp_new_vcmax$anpp_new[anpp_new_vcmax$exp=="soyfacesoy1_c"] <- add_anpp$anpp[add_anpp$exp=="soyfacesoy1_c"] 
anpp_new_vcmax$anpp_new[anpp_new_vcmax$exp=="soyfacesoy2_c"] <- add_anpp$anpp[add_anpp$exp=="soyfacesoy2_c"] 

#popface is g/m2 but forest, we don't know why so only just remove it!
final_mean <- merge(final_mean,anpp_new_vcmax[,c("exp","anpp_new")],by = c("exp"),all.x=TRUE)
final_mean$anpp <- final_mean$anpp_new
final_mean <- subset( final_mean, select = -anpp_new )


#add two plots from Cesar -->after adding them --becoming weaker!
#Cesar_anpp <- read.csv("/Users/yunpeng/data/gcme/cesar/ANPP.csv")
#phace <- log(112.6493/111.8533)/log(600/384)
#Aspen <- (log(669.6334/517.9303)/log(550/360) + log(621.7131/429.4247)/log(550/360))/2

#further look

#final_mean$anpp[final_mean$exp=="phace_c"] <- phace
#final_mean$anpp[final_mean$exp=="rhine-aspenface_c"] <- Aspen

csvfile <- paste("/Users/yunpeng/data/gcme/kevin/final_mean.csv")
write.csv(final_mean, csvfile, row.names = TRUE)


final_data <- final_mean[,c(4,6,7:10,12,14:50)]

res2 <- rcorr(as.matrix(final_data))
p_value <- as.data.frame(res2$P[,c("vcmax","jmax")])

vcmax_sig <- subset(p_value[order(p_value[,c("vcmax")] ),],vcmax<= 0.3)
jmax_sig <- subset(p_value[order(p_value[,c("jmax")] ),],jmax<= 0.3)

vcmax_select <- rownames(vcmax_sig)
jmax_select <- rownames(jmax_sig)

p <- list()
s <- list()
for(i in c(1:length(vcmax_select))){
  p[[i]] <- ggplot(final_mean,aes_string(x=vcmax_select[i], y="vcmax")) +
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=condition),size=3)+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+labs(y="vcmax")+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))}


#add mineral soil of vcmax
s[[1]] <- ggplot(final_mean,aes_string(x="soil_mineral_N", y="vcmax")) +
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=condition),size=3)+ #or type
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  labs(y="vcmax")+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))



plot_grid(p[[1]],p[[2]],p[[11]],p[[4]],p[[3]],p[[5]],p[[10]],p[[14]],s[[1]],p[[22]],p[[23]],p[[19]],
          nrow=4,label_size = 15)+
  theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/final_fig2_vcmax_simple.jpg",sep=""),width = 20, height = 20)


#jmax
b <- list()
#add mineral soil of jmax
s[[2]] <- ggplot(final_mean,aes_string(x="soil_mineral_N", y="jmax")) +
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=condition),size=3)+ #or type
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))

for(i in c(1:length(jmax_select))){
  b[[i]] <- ggplot(final_mean,aes_string(x=jmax_select[i], y="jmax")) +
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=condition),size=3)+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+labs(y="jmax")+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))}


plot_grid(b[[2]],b[[8]],b[[10]],b[[4]],b[[3]],b[[9]],
          nrow=2,label_size = 15)+theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/final_fig2_jmax.jpg",sep=""),width = 20, height = 20)


plot_grid(p[[1]],p[[2]],p[[4]],p[[3]],p[[5]],b[[2]],b[[4]],b[[3]],b[[9]],
          nrow=3,label_size = 15)+theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/fig2_final.jpg",sep=""),width = 20, height = 15)

plot_grid(p[[22]],p[[23]],s[[1]],b[[10]],b[[8]],s[[2]],
          nrow=2,label_size = 15)+theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/figs1_final.jpg",sep=""),width = 20, height = 10)



#new part - warming, fertilization and light
subset(kevin_othervars,treatment=="w")
subset(kevin_othervars,treatment=="s")
subset(kevin_othervars,treatment=="f")
kevin_othervars_wsf <- subset(kevin_othervars,treatment=="w"|treatment=="s"|treatment=="f")

#check numbers of sites and variables 
varname2 <- kevin_othervars_wsf%>% group_by(response)  %>% summarise(number = n())

#other vars - all created now
for (i in 1:nrow(varname2)) {
  tryCatch({
    varname2a <- varname2$response[i]
    df_w <- subset(kevin_othervars_wsf,treatment=="w" & response==varname2a)
    assign(paste("logr_w_", varname2a,sep=""), as_tibble(response_ratio_v2(df_w)))
    
    df_s <- subset(kevin_othervars_wsf,treatment=="s"& response==varname2a)
    assign(paste("logr_s_", varname2a,sep=""), as_tibble(response_ratio_v2(df_s)))
    
    df_f <- subset(kevin_othervars_wsf,treatment=="f" & response==varname2a)
    assign(paste("logr_f_", varname2a,sep=""), as_tibble(response_ratio_v2(df_f)))
  }, error=function(e){})} 


response_ratio_v2(logr_f_vcmax)[,c("logr","logr_var","exp")] %>%
  ggplot(aes(x=exp, y=logr)) +
  geom_point(aes(size=1/logr_var)) +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(x="", y="fertilization on vcmax",size=expression(paste("Variance"^{-1})))+
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

response_ratio_v2(logr_s_npp)[,c("logr","logr_var","exp")] %>%
  ggplot(aes(x=exp, y=logr)) +
  geom_point(aes(size=1/logr_var)) +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(x="", y="Shading on NPP",size=expression(paste("Variance"^{-1})))+
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=20))
  
#check
ggplot(final_mean,aes_string(x="jmax", y="vcmax")) +
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=condition),size=3)+
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  geom_smooth(color="black",method="lm",se=F)+labs(y="vcmax")+
  geom_text(aes(label=substr(exp, 1, 11)),hjust=0, vjust=0,check_overlap = F)+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"))

logr_c_vcmax[grep("trifolium", logr_c_vcmax$dominant_species),]$exp
final_mean
#swiss_face_trifolium
#aspen