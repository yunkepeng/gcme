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
vcmax_plot <- read.csv("/Users/yunpeng/data/gcme/kevin_20220222/vcmax_plotmean_v2.csv")$exp
jmax_plot <- read.csv("/Users/yunpeng/data/gcme/kevin_20220222/jmax_plotmean_v2.csv")$exp
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
                            (kevin_othervars_cf$response=="fine_root_production")] <- "bnpp"

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

vcmax_all_plot <- read.csv("/Users/yunpeng/data/gcme/kevin_20220222/vcmax_plotmean_v2.csv")[,c("ecosystem","exp","middle","condition","variance_info")]
vcmax_all_plot$condition[is.na(vcmax_all_plot$condition)==TRUE] <- "co2"
names(vcmax_all_plot) <- c("ecosystem","exp","vcmax","condition","variance_info")
jmax_all_plot <- read.csv("/Users/yunpeng/data/gcme/kevin_20220222/jmax_plotmean_v2.csv")[,c("exp","middle","condition")]
jmax_all_plot$condition[is.na(jmax_all_plot$condition)==TRUE] <- "co2"
names(jmax_all_plot) <- c("exp","jmax","condition")

final_mean <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp","condition"),all.x=TRUE),
                    list(vcmax_all_plot,jmax_all_plot,lma_plot,narea_plot,nmass_plot,bnpp_plot_final,soil_mineral_plotmean_dry,
                         anpp_plot,lai_plot,soil_total_n_plot,soil_cn_min_layer_plot,agb_plot,agb_n_plot,asat_plot,anet_plot,
                         bgb_plot,bgb_coarse_plot,bgb_n_plot,fine_root_biomass_plot,fine_root_n_plot,leaf_biomass_plot,leaf_c_plot,leaf_cn_plot,
                         leaf_litter_c_plot,leaf_litter_cn_plot,leaf_litter_n_plot,leaf_nue_plot,litter_biomass_plot,litterfall_plot,
                         mb_plot,mbc_plot,mbn_plot,r_leaf_plot,r_root_plot,r_soil_plot,root_n_uptake_plot,root_shoot_ratio_plot,soc_plot,soil_c_plot,
                         soil_total_c_plot,total_biomass_plot,total_biomass_n_plot,wood_n_plot,soil_potential_net_n_mineralization_plot))

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

#add popface's soil 
old_data <- read_csv("~/data/gcme/data_received_190325/NewData_wide_CORRECTED2.csv") %>%
  mutate( ambient_Sd  = as.numeric(ambient_Sd),  ambient_Se  = as.numeric(ambient_Se), 
          elevated_Sd = as.numeric(elevated_Sd), elevated_Se = as.numeric(elevated_Se),
          co2_a  = as.numeric(co2_a),  co2_e  = as.numeric(co2_e), 
          ambient  = as.numeric(ambient),  elevated  = as.numeric(elevated))

old_data$exp <- tolower(old_data$exp_nam)
popface <-subset(old_data, (exp_nam=="POPFACE_pa"|exp_nam=="POPFACE_pe"|exp_nam=="POPFACE_pn")&Data_type=="soil_mineral_N")
new_popface <- as.data.frame(agg_meta_sen_coef(response_ratio_v2(popface))[,c("exp","middle")])
names(new_popface) <- c("exp","soil_mineral_N")
final_mean$soil_mineral_N[final_mean$exp=="popface_pa_c"] <- new_popface$soil_mineral_N[new_popface$exp=="popface_pa"]
final_mean$soil_mineral_N[final_mean$exp=="popface_pe_c"] <- new_popface$soil_mineral_N[new_popface$exp=="popface_pe"]
final_mean$soil_mineral_N[final_mean$exp=="popface_pn_c"] <- new_popface$soil_mineral_N[new_popface$exp=="popface_pn"]

#finally, merged with ECM and ecosystem
#merge
#include ecm from gcme
ecm_csv <- na.omit(read.csv("/Users/yunpeng/data/gcme/kevin/orig_vcmax/new_ecm_types_updated.csv")[,c("exp","ecm_type")])

#include ecm from smith (see line 637 from code gcme_final_vcmas)
smith <- read.csv("/Users/yunpeng/data/smith_keenan_gcb/gcb_co2/pred_vcmax.csv")
smith_all_plotmean <- subset(smith,exp_nam!="NZFACE" &exp_nam!="BioCON" & exp_nam!="ChinaRiceFACE" & exp_nam!="DukeFACE" & exp_nam!="EUROPOPFACE" & exp_nam!="NevadaFACE" & exp_nam!="SwissFACE")
smith_all_plotmean$ecm_type <- smith_all_plotmean$Nac
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
#three are replicated - modify them
smith_all_plotmean$ecm_type[smith_all_plotmean$exp=="grassotc"] <- "Nfix" # 2 N-fixing vs. 4 AM
smith_all_plotmean$ecm_type[smith_all_plotmean$exp=="rhine-aspenface_c"] <- "ECM" # two ECM > 1 AM
smith_all_plotmean$ecm_type[smith_all_plotmean$exp=="macchia"] <- "AM" # two AM > 1 ECM
coord_smith <- unique(smith_all_plotmean[,c("exp","lon","lat","z")])
ecm_smith <- unique(smith_all_plotmean[,c("exp","ecm_type")])

ecm_all <- unique(rbind(ecm_smith,ecm_csv))

#merge to combine them
final_mean <- merge(final_mean,ecm_all,by=c("exp"),all.x=TRUE)

#some others are Nfix, change manually
unique(smith$exp_nam[smith$Nac=="Nfix"])
final_mean$ecm_type[final_mean$exp=="biocon_c"]<- "Nfix"
final_mean$ecm_type[final_mean$exp=="swissface_trifolium2_c"]<- "Nfix"
final_mean$ecm_type[final_mean$exp=="new_zealand_face_c"] <- "Nfix"

#make two divisions - ecosystem
final_mean$ecosystem_level <- final_mean$ecosystem
final_mean$ecosystem_level[final_mean$ecm_type=="Nfix"] <- "N-fixing"
final_mean$ecosystem_level[final_mean$condition=="Fertilization"] <- "High-N"

#make two divisions - ecm
final_mean$ecm_level <- final_mean$ecm_type
final_mean$ecm_level[final_mean$ecm_type=="Nfix"] <- "N-fixing"
final_mean$ecm_level[final_mean$condition=="Fertilization"] <- "High-N"

#corret one outlier of bgb_coarse
final_mean$bgb_coarse[final_mean$bgb_coarse< -5] <- NA
final_mean$bgb[final_mean$bgb< -5] <- NA

#add two plots from Cesar -->after adding them --becoming weaker!
#Cesar_anpp <- read.csv("/Users/yunpeng/data/gcme/cesar/ANPP.csv")
phace <- log(112.6493/111.8533)/log(600/384)
Aspen <- (log(669.6334/517.9303)/log(550/360) + log(621.7131/429.4247)/log(550/360))/2

#further look
final_mean$anpp[final_mean$exp=="phace_c"] <- phace
final_mean$anpp[final_mean$exp=="rhine-aspenface_c"] <- Aspen

csvfile <- paste("/Users/yunpeng/data/gcme/kevin/final_mean_v2.csv")
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
final_mean$ecosystem_level[final_mean$ecosystem_level=="forest"] <- "Others"
final_mean$ecosystem_level[final_mean$ecosystem_level=="grassland"] <- "Others"
final_mean$ecosystem_level[final_mean$ecosystem_level=="cropland"] <- "Others"



for(i in c(1:length(vcmax_select))){
  p[[i]] <- ggplot(final_mean,aes_string(x=vcmax_select[i], y="vcmax")) +
    geom_hline(yintercept=0)+geom_vline(xintercept=0)+
    geom_point(aes(color=ecosystem_level),size=3)+
    stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
    geom_smooth(color="black",method="lm",se=F)+labs(y=~paste(V[cmax]))+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"),legend.position="none")}

#add mineral soil of vcmax
s[[1]] <- ggplot(final_mean,aes_string(x="soil_mineral_N", y="vcmax")) +
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=ecosystem_level),size=3)+ #or type
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))+
  labs(y="vcmax")+
  theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"),legend.position="none")

#plot_grid(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[12]],s[[1]],nrow=3,label_size = 15)+theme(plot.background=element_rect(fill="white", color="white"))
#ggsave(paste("~/data/output_gcme/colin/fig2_vcmax_v2.jpg",sep=""),width = 15, height = 15)

#test
names(final_mean)
ggplot(final_mean,aes_string(x="soil_mineral_N", y="vcmax")) +
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  geom_point(aes(color=ecosystem_level),size=3)+ #or type
  stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")))

#for egu
a1 <- ggplot(final_mean,aes_string(x="jmax", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=ecosystem_level),size=3)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+geom_smooth(color="black",method="lm",se=F)+theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(V[cmax]))+labs(x=~paste(J[cmax]))
a2 <- ggplot(final_mean,aes_string(x="nmass", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=ecosystem_level),size=3)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+geom_smooth(color="black",method="lm",se=F)+theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(V[cmax]))+labs(x=~paste(N[mass]))
a3 <- ggplot(final_mean,aes_string(x="LMA", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=ecosystem_level),size=3)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+geom_smooth(color="black",method="lm",se=F)+theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(V[cmax]))+labs(x=~paste(LMA))
a4 <- ggplot(final_mean,aes_string(x="anpp", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=ecosystem_level),size=3)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+geom_smooth(color="black",method="lm",se=F)+theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(V[cmax]))+labs(x=~paste(ANPP))
a5 <- ggplot(final_mean,aes_string(x="bnpp", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=ecosystem_level),size=3)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+geom_smooth(color="black",method="lm",se=F)+theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(V[cmax]))+labs(x=~paste(BNPP))
a6 <- ggplot(final_mean,aes_string(x="root_shoot_ratio", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=ecosystem_level),size=3)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+geom_smooth(color="black",method="lm",se=F)+theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(V[cmax]))+labs(x=~paste(Root/Shoot))
a7 <- ggplot(final_mean,aes_string(x="litter_biomass", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=ecosystem_level),size=3)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+geom_smooth(color="black",method="lm",se=F)+theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(V[cmax]))+labs(x=~paste("litter_biomass"))
a8 <- ggplot(final_mean,aes_string(x="soil_mineral_N", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=ecosystem_level),size=3)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+geom_smooth(color="black",method="lm",se=F)+theme(axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"),legend.position="none")+
  labs(y=~paste(V[cmax]))+labs(x=~paste("Soil inorganic N"))

final1_legend <- ggplot(final_mean,aes_string(x="soil_mineral_N", y="vcmax")) +geom_hline(yintercept=0)+geom_vline(xintercept=0)+geom_point(aes(color=ecosystem_level),size=3)+stat_cor(aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")),size=7)+geom_smooth(color="black",method="lm",se=F)+theme(text = element_text(size=30),axis.text=element_text(size=25),axis.title=element_text(size=25,face="bold"))+
  labs(y=~paste(V[cmax]))+labs(x=~paste("Soil inorganic N"))+scale_colour_discrete(" ")

legend_info <- as_ggplot(get_legend(final1_legend))

plot_grid(a1,a2,a3,a4,a5,a6,a7,a8,legend_info,nrow=3,label_size = 15)+
  theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/fig2_MS.jpg",sep=""),width = 15, height = 15)

#site distribution
kevin <- unique(read.csv("/Users/yunpeng/data/gcme/kevin/orig_vcmax/JunkePeng_11252021.csv")[,c("lon","lat")])
smith <- unique(read.csv("/Users/yunpeng/data/smith_keenan_gcb/gcb_co2/pred_vcmax.csv")[,c("lon","lat")])
vcmax_sites <- rbind(kevin,smith)

newmap <- getMap(resolution = "low")
sp::plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)

points(vcmax_sites$lon,vcmax_sites$lat, col="green", pch=16,cex=2)

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
    theme(axis.text=element_text(size=20),axis.title=element_text(size=20,face="bold"),legend.position="none")}


#plot_grid(b[[2]],b[[3]],b[[5]],b[[6]],b[[7]],b[[8]],b[[9]],s[[2]], nrow=3,label_size = 15)+theme(plot.background=element_rect(fill="white", color="white"))
#ggsave(paste("~/data/output_gcme/colin/final_fig2_jmax_v2.jpg",sep=""),width = 20, height = 20)

#new part - warming, fertilization and light
kevin_othervars_wsf <- subset(kevin_othervars, is.na(exp)==FALSE)

###pre-processing: 
#anpp and anpp_grain are both anpp.
kevin_othervars_wsf$response[kevin_othervars_wsf$response=="anpp_grain"] <- "anpp"
unique(subset(kevin_othervars_wsf,response=="anpp")$Unit)

#some tiny change: leaf_n needs division - but leaf_p doesn't need (as it is all pmass)
nmass_wsf <- subset(kevin_othervars_wsf,response=="leaf_n") %>% filter(Unit %in% c("g","g_100g","g_g","g_kg","mg_g","mg_kg","mmol_g","ug_mg","umol_g"))
nmass_wsf$response <- "nmass"
logr_w_nmass <- as_tibble(response_ratio_v2(subset(nmass_wsf,treatment=="w")));
all_logr_c_nmass <- as_tibble(response_ratio_v2(subset(nmass_wsf,treatment=="c")))
logr_d_nmass <- as_tibble(response_ratio_v2(subset(nmass_wsf,treatment=="d")));
logr_f_nmass <- as_tibble(response_ratio_v2(subset(nmass_wsf,treatment=="f")))


#bnpp - make sure data between them are not mixed.
kevin_othervars_wsf$output[kevin_othervars_wsf$response=="coarse_root_production"|
                             kevin_othervars_wsf$response=="fine_root_turnover"|
                             (kevin_othervars_wsf$response=="root_production")|
                            (kevin_othervars_wsf$response=="fine_root_biomass" & kevin_othervars_wsf$Unit=="g_m2y")|
                            (kevin_othervars_wsf$response=="fine_root_production")] <- "bnpp" # this g/m2 should be in cropland

bnpp_dataset <- subset(kevin_othervars_wsf,output=="bnpp")
unique(bnpp_dataset[,c("response","Unit")])

bnpp_dataset$response <- "bnpp"

all_logr_c_bnpp <- as_tibble(response_ratio_v2(subset(bnpp_dataset,treatment=="c")));
logr_w_bnpp <- as_tibble(response_ratio_v2(subset(bnpp_dataset,treatment=="w")))
logr_f_bnpp <- as_tibble(response_ratio_v2(subset(bnpp_dataset,treatment=="f")))
logr_d_bnpp <- as_tibble(response_ratio_v2(subset(bnpp_dataset,treatment=="d")))

#create a column to combine exp and exp_nam
exp_name <- unique(kevin_othervars[,c("exp","exp_nam")])

#other vars - all created now
unique(subset(kevin_othervars_wsf,response=="vcmax")$treatment)

logr_f_vcmax <- subset(kevin_othervars_wsf,treatment=="f" & response=="vcmax")
logr_d_vcmax <- subset(kevin_othervars_wsf,treatment=="d" & response=="vcmax")
logr_w_vcmax <- subset(kevin_othervars_wsf,treatment=="w" & response=="vcmax")

logr_w_jmax <- subset(kevin_othervars_wsf,treatment=="w" & response=="jmax")
logr_f_jmax <- subset(kevin_othervars_wsf,treatment=="f" & response=="jmax")
logr_d_jmax <- subset(kevin_othervars_wsf,treatment=="d" & response=="jmax")

logr_f_anpp <- subset(kevin_othervars_wsf,treatment=="f" & response=="anpp")
logr_d_anpp <- subset(kevin_othervars_wsf,treatment=="d" & response=="anpp")
logr_w_anpp <- subset(kevin_othervars_wsf,treatment=="w" & response=="anpp")

logr_f_root_shoot_ratio <- subset(kevin_othervars_wsf,treatment=="f" & response=="root_shoot_ratio")
logr_d_root_shoot_ratio <- subset(kevin_othervars_wsf,treatment=="d" & response=="root_shoot_ratio")
logr_w_root_shoot_ratio <- subset(kevin_othervars_wsf,treatment=="w" & response=="root_shoot_ratio")
logr_c_root_shoot_ratio <- subset(kevin_othervars_wsf,treatment=="c" & response=="root_shoot_ratio")

#create a function, with a input of list of variable, and output a potential figure

#create function
output_point <- function(list_of_variable){
  test <- list_of_variable
  point_list <- list() 
  boxplot_list <- list() 
  for (i in 1:length(test)){ 
    
    data1 <- assign(test[i], get(test[i]))
    a1 <- agg_meta(response_ratio_v2(data1))
    a2 <- merge(a1[,c("middle","logr_var","exp","no")],exp_name,by=c("exp"),all.x=TRUE)

    a2$type_name <- test[i]
    
    output2 <- tibble(test[i],middle=mean(a1$middle,na.rm=TRUE))
    
    point_list[[i]] <- a2
    boxplot_list[[i]] <- output2
    }
  point_list_final <- do.call("rbind",point_list)
  boxplot_list_final <- do.call("rbind",boxplot_list)
  return(point_list_final)
}

vcmax_point <- output_point(c("logr_c_vcmax","logr_f_vcmax","logr_d_vcmax","logr_w_vcmax"))

vc_points <- vcmax_point %>% ggplot( aes(x=type_name, y=middle)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(y="vcmax",x=" ") + theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

jmax_point <- output_point(c("logr_c_jmax","logr_f_jmax","logr_d_jmax","logr_w_jmax"))

j_points <- jmax_point %>% ggplot( aes(x=type_name, y=middle)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(y="jmax",x=" ") + theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))


#now, additionally, shading
walker <- read.csv("/Users/yunpeng/data/leaf_traits/Walker/LEAF_PHOTOSYNTHESIS_TRAITS_1224/data/Leaf_Photosynthesis_Traits.csv")

#vcmax
light_vcmax <- subset(walker,Treatment_light!="not applicable") %>% 
  group_by(lon,lat,species,Treatment_light)  %>% summarise(mean = mean(Vcmax))

light1 <- subset(light_vcmax,Treatment_light=="shade"|Treatment_light=="sun")
light1a <- merge(subset(light1,Treatment_light=="shade"),subset(light1,Treatment_light=="sun"),
      by=c("lon","lat","species"),all.x=TRUE)
light1a$middle <- log(light1a$mean.y/light1a$mean.x)
light1a$type_name <- "shade_to_sun"

light2 <- subset(light_vcmax,Treatment_light=="high"|Treatment_light=="low")
light2a <- merge(subset(light2,Treatment_light=="low"),subset(light2,Treatment_light=="high"),
                 by=c("lon","lat","species"),all.x=TRUE)
light2a$middle <- log(light2a$mean.y/light2a$mean.x)
light2a$type_name <- "low_to_high_light"

#light3 <- subset(light_vcmax,Treatment_light=="lower"|Treatment_light=="middle"|Treatment_light=="upper")
#light3a <-Reduce(function(x,y) merge(x = x, y = y,by=c("lon","lat","species"),all.x=TRUE),
#                    list(subset(light3,Treatment_light=="lower"),subset(light3,Treatment_light=="middle"),subset(light3,Treatment_light=="upper")))
#light3a$middle <- log(light3a$mean/light3a$mean.x) # divide upper to lower - not considering middle
#light3a$type_name <- "lower_to_upper_canopy"

light_vcmax_points <- dplyr::bind_rows(light1a[,c("lon","lat","species","middle","type_name")],
                 light2a[,c("lon","lat","species","middle","type_name")]) 

light_vcmax_fig <- light_vcmax_points %>% ggplot( aes(x=type_name, y=middle)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(y="light effect on vcmax",x=" ") + theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

#jmax
light_jmax <- subset(walker,Treatment_light!="not applicable" & Jmax>0) %>% 
  group_by(lon,lat,species,Treatment_light)  %>% summarise(mean = mean(Jmax))

light1 <- subset(light_jmax,Treatment_light=="shade"|Treatment_light=="sun")
light1a <- merge(subset(light1,Treatment_light=="shade"),subset(light1,Treatment_light=="sun"),
                 by=c("lon","lat","species"),all.x=TRUE)
light1a$middle <- log(light1a$mean.y/light1a$mean.x)
light1a$type_name <- "shade_to_sun"

light2 <- subset(light_jmax,Treatment_light=="high"|Treatment_light=="low")
light2a <- merge(subset(light2,Treatment_light=="low"),subset(light2,Treatment_light=="high"),
                 by=c("lon","lat","species"),all.x=TRUE)
light2a$middle <- log(light2a$mean.y/light2a$mean.x)
light2a$type_name <- "low_to_high_light"

#light3 <- subset(light_jmax,Treatment_light=="lower"|Treatment_light=="middle"|Treatment_light=="upper")
#light3a <-Reduce(function(x,y) merge(x = x, y = y,by=c("lon","lat","species"),all.x=TRUE),
#                 list(subset(light3,Treatment_light=="lower"),subset(light3,Treatment_light=="middle"),subset(light3,Treatment_light=="upper")))
#light3a$middle <- log(light3a$mean/light3a$mean.x) # divide upper to lower - not considering middle
#light3a$type_name <- "lower_to_upper_canopy"

light_jmax_points <- dplyr::bind_rows(light1a[,c("lon","lat","species","middle","type_name")],
                                light2a[,c("lon","lat","species","middle","type_name")]) 
#when it is 0, it means both are -10000, removed it
light_jmax_fig <- light_jmax_points %>% ggplot( aes(x=type_name, y=middle)) + geom_boxplot()+
  geom_point(alpha = 0.6, width = 0.5) +geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(y="light effect on jmax",x=" ") + theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

#combine
warming <- read.csv("/Users/yunpeng/data/gcme/kevin_20220222/vcmaxjmax_warming_plotmean_v2.csv")
warming$treatment <- "warming"

vcmax_dry <- subset(vcmax_point,type_name=="logr_d_vcmax")
jmax_dry <- subset(jmax_point,type_name=="logr_d_jmax")
vj_dry<- (merge(vcmax_dry[,c("exp","middle")],jmax_dry[,c("exp","middle")],
              by=c("exp"),all.x=TRUE,all.y=TRUE))
names(vj_dry) <- c("exp","obs_vcmax","obs_jmax")
vj_dry$obs_jmax_vcmax <- vj_dry$obs_jmax - vj_dry$obs_vcmax
vj_dry$treatment <- "drying";vj_dry$pred_vcmax <- NA ;vj_dry$pred_jmax <- NA ;vj_dry$pred_jmax_vcmax <- NA 
    
vcmax_fer <- subset(vcmax_point,type_name=="logr_f_vcmax");vcmax_fer$treatment <- "fertilization"
jmax_fer <- subset(jmax_point,type_name=="logr_f_jmax");jmax_fer$treatment <- "fertilization"
vj_fer<- (merge(vcmax_fer[,c("exp","middle")],jmax_fer[,c("exp","middle")],
                by=c("exp"),all.x=TRUE,all.y=TRUE))
names(vj_fer) <- c("exp","obs_vcmax","obs_jmax")
vj_fer$obs_jmax_vcmax <- vj_fer$obs_jmax - vj_fer$obs_vcmax
vj_fer$treatment <- "fertilization";vj_fer$pred_vcmax <- NA ;vj_fer$pred_jmax <- NA ;vj_fer$pred_jmax_vcmax <- NA 

light_vcmax_points$exp <- paste(light_vcmax_points$lon,light_vcmax_points$lat,light_vcmax_points$species,sep="_")
light_jmax_points$exp <- paste(light_jmax_points$lon,light_jmax_points$lat,light_jmax_points$species,sep="_")
vj_light <- (merge(light_vcmax_points[,c("exp","middle","type_name")],light_jmax_points[,c("exp","middle","type_name")],
                by=c("exp","type_name"),all.x=TRUE,all.y=TRUE))
names(vj_light) <- c("exp","type_name","obs_vcmax","obs_jmax")

vj_light$obs_jmax_vcmax <- vj_light$obs_jmax - vj_light$obs_vcmax
vj_light$treatment <- "light";vj_light$pred_vcmax <- NA ;vj_light$pred_jmax <- NA ;vj_light$pred_jmax_vcmax <- NA 

ymin_max <- function(vj_all){
  vj_all$ymin_vc_obs <- quantile(vj_all$obs_vcmax, 0.25,na.rm=T)
  vj_all$ymax_vc_obs <- quantile(vj_all$obs_vcmax, 0.75,na.rm=T)
  vj_all$mean_vc_obs <- median(vj_all$obs_vcmax,na.rm=T)
  
  vj_all$ymin_j_obs <- quantile(vj_all$obs_jmax, 0.25,na.rm=T)
  vj_all$ymax_j_obs <- quantile(vj_all$obs_jmax, 0.75,na.rm=T)
  vj_all$mean_j_obs <- median(vj_all$obs_jmax,na.rm=T)
  
  vj_all$ymin_jv_obs <- quantile(vj_all$obs_jmax_vcmax, 0.25,na.rm=T)
  vj_all$ymax_jv_obs <- quantile(vj_all$obs_jmax_vcmax, 0.75,na.rm=T)
  vj_all$mean_jv_obs <- median(vj_all$obs_jmax_vcmax,na.rm=T)
  
  vj_all$ymin_vc_pred <- quantile(vj_all$pred_vcmax, 0.25,na.rm=T)
  vj_all$ymax_vc_pred <- quantile(vj_all$pred_vcmax, 0.75,na.rm=T)
  vj_all$mean_vc_pred <- median(vj_all$pred_vcmax,na.rm=T)
  
  vj_all$ymin_j_pred <- quantile(vj_all$pred_jmax, 0.25,na.rm=T)
  vj_all$ymax_j_pred <- quantile(vj_all$pred_jmax, 0.75,na.rm=T)
  vj_all$mean_j_pred <- median(vj_all$pred_jmax,na.rm=T)
  
  vj_all$ymin_jv_pred <- quantile(vj_all$pred_jmax_vcmax, 0.25,na.rm=T)
  vj_all$ymax_jv_pred <- quantile(vj_all$pred_jmax_vcmax, 0.75,na.rm=T)
  vj_all$mean_jv_pred <- median(vj_all$pred_jmax_vcmax,na.rm=T)
  
  return(vj_all[1,c("treatment","ymin_vc_obs","ymax_vc_obs","mean_vc_obs",
                    "ymin_j_obs","ymax_j_obs","mean_j_obs",
                    "ymin_jv_obs","ymax_jv_obs","mean_jv_obs",
                    "ymin_vc_pred","ymax_vc_pred","mean_vc_pred",
                    "ymin_j_pred","ymax_j_pred","mean_j_pred",
                    "ymin_jv_pred","ymax_jv_pred","mean_jv_pred")])
}

box_warming <- ymin_max(warming)
box_dry <- ymin_max(vj_dry)
box_fer <- ymin_max(vj_fer)
box_light <- ymin_max(vj_light)
box_type <- dplyr::bind_rows(box_warming,box_dry,box_fer,box_light)

vj_all <- dplyr::bind_rows(warming,vj_dry,vj_fer,vj_light)

others_v <- vj_all %>%
  ggplot( aes(x=treatment, y=obs_vcmax)) +
  geom_boxplot(width = 0.5)+
  geom_crossbar(data=box_type,aes(x=treatment,y=mean_vc_obs, ymin=ymin_vc_obs, ymax=ymax_vc_obs), alpha = 0.6, width = 0.5,color="black") +
  geom_crossbar(data=box_type,aes(x=treatment,y=mean_vc_pred, ymin=ymin_vc_pred, ymax=ymax_vc_pred), alpha = 0.6, width = 0.5,color="red") +
  #geom_point(aes(color=source),alpha = 0.6, width = 0.5,size=2) +
  geom_point(alpha = 0.6, width = 0.5,size=2) +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(x="", y="Vcmax response") +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

others_j <- vj_all %>%
  ggplot( aes(x=treatment, y=obs_jmax)) +
  geom_boxplot(width = 0.5)+
  geom_crossbar(data=box_type,aes(x=treatment,y=mean_j_obs, ymin=ymin_j_obs, ymax=ymax_j_obs), alpha = 0.6, width = 0.5,color="black") +
  geom_crossbar(data=box_type,aes(x=treatment,y=mean_j_pred, ymin=ymin_j_pred, ymax=ymax_j_pred), alpha = 0.6, width = 0.5,color="red") +
  #geom_point(aes(color=source),alpha = 0.6, width = 0.5,size=2) +
  geom_point(alpha = 0.6, width = 0.5,size=2) +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(x="", y="Jmax response") +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

others_jv <- vj_all %>%
  ggplot( aes(x=treatment, y=obs_jmax_vcmax)) +
  geom_boxplot(width = 0.5)+
  geom_crossbar(data=box_type,aes(x=treatment,y=mean_jv_obs, ymin=ymin_jv_obs, ymax=ymax_jv_obs), alpha = 0.6, width = 0.5,color="black") +
  geom_crossbar(data=box_type,aes(x=treatment,y=mean_jv_pred, ymin=ymin_jv_pred, ymax=ymax_jv_pred), alpha = 0.6, width = 0.5,color="red") +
  #geom_point(aes(color=source),alpha = 0.6, width = 0.5,size=2) +
  geom_point(alpha = 0.6, width = 0.5,size=2) +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(x="", y="Jmax/Vcmax response") +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

plot_grid(others_v,others_j,others_jv,nrow=1,label_size = 15)+theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/final_fig1b_alternative_v2.jpg",sep=""),width = 15, height = 5)


#N-fertilized under CO2 experiment
#vcmax
Nfer_vcmax <- subset(walker,Treatment_N!="not applicable" &Treatment_CO2!="not applicable" & Vcmax>0) %>% 
  group_by(lon,lat,species,Treatment_N,Treatment_CO2)  %>% summarise(mean = mean(Vcmax))

highN_vcmax <- merge(subset(Nfer_vcmax,Treatment_N=="high"&Treatment_CO2=="amb"),
                            subset(Nfer_vcmax,Treatment_N=="high"&Treatment_CO2=="ele"),
                 by=c("lon","lat","species","Treatment_N"),all.x=TRUE)

highN_vcmax$middle <- log(highN_vcmax$mean.y/highN_vcmax$mean.x)

highN_vcmax$type_name <- "CO2 response at high N"
highN_vcmax$site_species <- paste(highN_vcmax$lon,highN_vcmax$lat,highN_vcmax$species,sep="_")

lowN_vcmax <- merge(subset(Nfer_vcmax,Treatment_N=="low"&Treatment_CO2=="amb"),
                     subset(Nfer_vcmax,Treatment_N=="low"&Treatment_CO2=="ele"),
                     by=c("lon","lat","species","Treatment_N"),all.x=TRUE)

lowN_vcmax$middle <- log(lowN_vcmax$mean.y/lowN_vcmax$mean.x)

lowN_vcmax$type_name <- "CO2 response at low N"
lowN_vcmax$site_species <- paste(lowN_vcmax$lon,lowN_vcmax$lat,lowN_vcmax$species,sep="_")

nfer_gcme <- final_mean %>% filter(exp %in%c("duke2_cf","euroface4_pa_cf","euroface4_pe_cf","euroface4_pn_cf","new_zealand_face_cf",
                                             "duke2_c","euroface4_pa_c","euroface4_pe_c","euroface4_pn_c","new_zealand_face_c"))
nfer_gcme2 <- nfer_gcme[,c("vcmax","condition","exp")]
names(nfer_gcme2) <- c("middle","type_name","site_species")
nfer_gcme2$type_name[nfer_gcme2$type_name=="co2"] <- "CO2 response without N"
nfer_gcme2$type_name[nfer_gcme2$type_name=="Fertilization"] <- "CO2 response at high N"
nfer_gcme2$site_species <- c("duke2","duke2","euroface-pa","euroface-pa",
                             "euroface-pe","euroface-pe","euroface-pn","euroface-pn",
                             "new_zealand","new_zealand")

vcmax_co2nfer_points <- dplyr::bind_rows(lowN_vcmax[,c("site_species","middle","type_name")],
                                       highN_vcmax[,c("site_species","middle","type_name")],
                                       nfer_gcme2) 

vcmax_co2nfer_fig <- vcmax_co2nfer_points %>% ggplot( aes(x=site_species, y=middle)) + 
  geom_point(aes(color=type_name),size=3) +geom_hline( yintercept=0.0, size=0.5)+
  labs(y="CO2 effect on vcmax at Fertilized points",x=" ") + theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

a1 <- subset(vcmax_co2nfer_points,type_name=="CO2 response at high N")[,c("middle","site_species")]
a2 <- subset(vcmax_co2nfer_points,type_name=="CO2 response at low N")[,c("middle","site_species")]
a3 <- subset(vcmax_co2nfer_points,type_name=="CO2 response without N")[,c("middle","site_species")]

final_nfer_c <- Reduce(function(x,y) merge(x = x, y = y,by=c("site_species"),all.x=TRUE),list(a1,a2,a3))
names(final_nfer_c) <- c("exp","high_N","low_N","without_N")
final_nfer_c$final_N_high <- final_nfer_c$high_N
final_nfer_c$final_N_low <- final_nfer_c$low_N
final_nfer_c$final_N_low[is.na(final_nfer_c$without_N)==F] <- final_nfer_c$without_N[is.na(final_nfer_c$without_N)==F]

#now, jmax
Nfer_jmax <- subset(walker,Treatment_N!="not applicable" &Treatment_CO2!="not applicable" & Jmax>0) %>% 
  group_by(lon,lat,species,Treatment_N,Treatment_CO2)  %>% summarise(mean = mean(Jmax))

highN_jmax <- merge(subset(Nfer_jmax,Treatment_N=="high"&Treatment_CO2=="amb"),
                     subset(Nfer_jmax,Treatment_N=="high"&Treatment_CO2=="ele"),
                     by=c("lon","lat","species","Treatment_N"),all.x=TRUE)

highN_jmax$middle <- log(highN_jmax$mean.y/highN_jmax$mean.x)

highN_jmax$type_name <- "CO2 response at high N"
highN_jmax$site_species <- paste(highN_jmax$lon,highN_jmax$lat,highN_jmax$species,sep="_")

lowN_jmax <- merge(subset(Nfer_jmax,Treatment_N=="low"&Treatment_CO2=="amb"),
                    subset(Nfer_jmax,Treatment_N=="low"&Treatment_CO2=="ele"),
                    by=c("lon","lat","species","Treatment_N"),all.x=TRUE)

lowN_jmax$middle <- log(lowN_jmax$mean.y/lowN_jmax$mean.x)

lowN_jmax$type_name <- "CO2 response at low N"
lowN_jmax$site_species <- paste(lowN_jmax$lon,lowN_jmax$lat,lowN_jmax$species,sep="_")

nfer_gcme <- final_mean %>% filter(exp %in%c("duke2_cf","euroface4_pa_cf","euroface4_pe_cf","euroface4_pn_cf",
                                             "duke2_c","euroface4_pa_c","euroface4_pe_c","euroface4_pn_c"))

nfer_gcme2 <- nfer_gcme[,c("jmax","condition","exp")]
names(nfer_gcme2) <- c("middle","type_name","site_species")
nfer_gcme2$type_name[nfer_gcme2$type_name=="co2"] <- "CO2 response without N"
nfer_gcme2$type_name[nfer_gcme2$type_name=="Fertilization"] <- "CO2 response at high N"
nfer_gcme2$site_species <- c("duke2","duke2","euroface-pa","euroface-pa",
                             "euroface-pe","euroface-pe","euroface-pn","euroface-pn")

jmax_co2nfer_points <- dplyr::bind_rows(lowN_jmax[,c("site_species","middle","type_name")],
                                         highN_jmax[,c("site_species","middle","type_name")],
                                         nfer_gcme2) 

b1 <- subset(jmax_co2nfer_points,type_name=="CO2 response at high N")[,c("middle","site_species")]
b2 <- subset(jmax_co2nfer_points,type_name=="CO2 response at low N")[,c("middle","site_species")]
b3 <- subset(jmax_co2nfer_points,type_name=="CO2 response without N")[,c("middle","site_species")]

final_nfer_j <- Reduce(function(x,y) merge(x = x, y = y,by=c("site_species"),all.x=TRUE),list(b1,b2,b3))
names(final_nfer_j) <- c("exp","high_N","low_N","without_N")
final_nfer_j$final_N_high <- final_nfer_j$high_N
final_nfer_j$final_N_low <- final_nfer_j$low_N
final_nfer_j$final_N_low[is.na(final_nfer_j$without_N)==F] <- final_nfer_j$without_N[is.na(final_nfer_j$without_N)==F]

final_nfer_vj <- merge(final_nfer_c[,c("exp","final_N_high","final_N_low")],
                    final_nfer_j[,c("exp","final_N_high","final_N_low")],
                    by=c("exp"),all.x=TRUE,all.y=TRUE)
names(final_nfer_vj) <- c("exp","vcmax_high","vcmax_low","jmax_high","jmax_low")

final_nfer_vj$jv_high <- final_nfer_vj$jmax_high-final_nfer_vj$vcmax_high
final_nfer_vj$jv_low <- final_nfer_vj$jmax_low -final_nfer_vj$vcmax_low
final_nfer_vj$treatment <- "Nfer+co2"

box1 <- tibble(treatment = "Nfer+co2",middle=median(final_nfer_vj$vcmax_high),ymin=quantile(final_nfer_vj$vcmax_high, 0.25),
                            ymax=quantile(final_nfer_vj$vcmax_high, 0.75))
box2 <- tibble(treatment = "Nfer+co2",middle=median(final_nfer_vj$vcmax_low),ymin=quantile(final_nfer_vj$vcmax_low, 0.25),
               ymax=quantile(final_nfer_vj$vcmax_low, 0.75))

f1 <- final_nfer_vj %>%
  ggplot( aes(x=treatment, y=vcmax_high)) +
  geom_boxplot(width = 0.5,color="red")+
  geom_boxplot(aes(x=treatment, y=vcmax_low),width = 0.5,color="blue")+
  geom_crossbar(data=box1,aes(x=treatment,y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="red") +
  geom_crossbar(data=box2,aes(x=treatment,y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="blue") +
  geom_point(alpha = 0.6, width = 0.5,size=2,color="red") +
  geom_point(aes(x=treatment, y=vcmax_low),alpha = 0.6, width = 0.5,size=2,color="blue")  +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(x="", y="Vcmax response") +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))


box3 <- tibble(treatment = "Nfer+co2",middle=median(final_nfer_vj$jmax_high,na.rm = T),ymin=quantile(final_nfer_vj$jmax_high, 0.25,na.rm = T),
               ymax=quantile(final_nfer_vj$jmax_high, 0.75,na.rm = T))
box4 <- tibble(treatment = "Nfer+co2",middle=median(final_nfer_vj$jmax_low,na.rm = T),ymin=quantile(final_nfer_vj$jmax_low, 0.25,na.rm = T),
               ymax=quantile(final_nfer_vj$jmax_low, 0.75,na.rm = T))

f2 <- final_nfer_vj %>%
  ggplot( aes(x=treatment, y=jmax_high)) +
  geom_boxplot(width = 0.5,color="red")+
  geom_boxplot(aes(x=treatment, y=jmax_low),width = 0.5,color="blue")+
  geom_crossbar(data=box3,aes(x=treatment,y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="red") +
  geom_crossbar(data=box4,aes(x=treatment,y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="blue") +
  geom_point(alpha = 0.6, width = 0.5,size=2,color="red") +
  geom_point(aes(x=treatment, y=jmax_low),alpha = 0.6, width = 0.5,size=2,color="blue")  +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(x="", y="Jmax response") +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

box5 <- tibble(treatment = "Nfer+co2",middle=median(final_nfer_vj$jv_high,na.rm = T),ymin=quantile(final_nfer_vj$jv_high, 0.25,na.rm = T),
               ymax=quantile(final_nfer_vj$jv_high, 0.75,na.rm = T))
box6 <- tibble(treatment = "Nfer+co2",middle=median(final_nfer_vj$jv_low,na.rm = T),ymin=quantile(final_nfer_vj$jv_low, 0.25,na.rm = T),
               ymax=quantile(final_nfer_vj$jv_low, 0.75,na.rm = T))

f3 <- final_nfer_vj %>%
  ggplot( aes(x=treatment, y=jv_high)) +
  geom_boxplot(width = 0.5,color="red")+
  geom_boxplot(aes(x=treatment, y=jv_low),width = 0.5,color="blue")+
  geom_crossbar(data=box5,aes(x=treatment,y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="red") +
  geom_crossbar(data=box6,aes(x=treatment,y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5,color="blue") +
  geom_point(alpha = 0.6, width = 0.5,size=2,color="red") +
  geom_point(aes(x=treatment, y=jv_low),alpha = 0.6, width = 0.5,size=2,color="blue")  +
  geom_hline( yintercept=0.0, size=0.5)+ ylim(-1,1)+
  labs(x="", y="Jmax/Vcmax response") +
  theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

plot_grid(f1,f2,f3,nrow=1,label_size = 15)+theme(plot.background=element_rect(fill="white", color="white"))
ggsave(paste("~/data/output_gcme/colin/final_fig1c_alternative_v2.jpg",sep=""),width = 15, height = 5)

#check final_mean for PCA
final2 <- final_mean[,c("exp","lai","vcmax","jmax","narea","LMA","nmass","bnpp","anpp","root_shoot_ratio")]
final2$count_na <- rowSums(is.na(final2))
final3 <- subset(final2,count_na<6)
final3 <- unique(final3[order(final3$count_na),])
head(final3)
summary(final3)

#include lon/lat - remove one repoeated data
coord_all <- dplyr::bind_rows(coord_smith,unique((read.csv("/Users/yunpeng/data/gcme/kevin/forcing/pred_vcmax_v2.csv")[,c("lon","lat","z","exp")])))
subset(coord_all,exp=="new_zealand_face_c")
coord_all$z[coord_all$z==16 & coord_all$exp=="new_zealand_face_c"] <- NA
coord_all <- na.omit(coord_all)

final4 <- merge(final3,coord_all,by=c("exp"),all.x=TRUE)

unique(final4[,c("lon","lat","z")])

#check
aa <- subset(kevin_othervars,exp=="eucface_c")%>% group_by(response,dominant_species,sampling_year,sampling_date)  %>% summarise(number = n())


###eucface - LAI
#from https://onlinelibrary.wiley.com/doi/full/10.1111/gcb.13151?casa_token=6CKcWQ_OHHwAAAAA%3AEsLJPJXb45rz2WxE807-NvACiQmFkELScHJiV_eaRUEPd0psT7co5ZnJp8Mo7CKaPFt4H6dkKe8XqZFXJw
eucface_lai_df <- read.csv("/Users/yunpeng/data/Duursma_gcb/EucFACE_DUURSMA_GCB_LEAFAREAINDEX/data/FACE_RA_P0037_GAPFRACLAI_20121026-20150225_L2.csv")
#150 as given in paper, also consistent with our df
subset(logr_c_vcmax,exp=="eucface_c")$co2_e[1];subset(logr_c_vcmax,exp=="eucface_c")$co2_a[1]
eucface_lai <- log(mean(subset(eucface_lai_df,treatment=="elevated")$LAI,na.rm=TRUE)/mean(subset(eucface_lai_df,treatment=="ambient")$LAI,na.rm=TRUE))/log(540/394)

###biocon_c - LAI
#LMA, everything is LMA (g/cm2) now, though shown as cm2/g - converting to g/m2
#lma: g/m2
lma_a <-10000*mean(subset(logr_c_LMA,exp=="biocon_c"&Unit=="cm2/g")$ambient)
lma_e <-10000*mean(subset(logr_c_LMA,exp=="biocon_c"&Unit=="cm2/g")$elevated,na.rm=TRUE)
check <- subset(kevin_othervars,exp=="biocon_c")%>% group_by(response,Unit)  %>% summarise(number = n())
#agb: g/m2
agb_a <- mean(subset(kevin_othervars,exp=="biocon_c" & response=="agb" & Unit=="g_m2")$ambient)
agb_e <- mean(subset(kevin_othervars,exp=="biocon_c" & response=="agb" & Unit=="g_m2")$elevated)
agb_a/lma_a;agb_e/lma_e # looks ok
subset(logr_c_vcmax,exp=="biocon_c")$co2_a;subset(logr_c_vcmax,exp=="biocon_c")$co2_e
biocon_lai <- log((agb_e/lma_e)/(agb_a/lma_a))/log(570/367)

#??? biocon_c - root/shoot
check

###brandbjerg_c
#LAI

#anpp


#df_only <- read_csv("~/data/gcme/data_received_190325/NewData_wide_CORRECTED2.csv") %>%
#  mutate( ambient_Sd  = as.numeric(ambient_Sd),  ambient_Se  = as.numeric(ambient_Se), 
#          elevated_Sd = as.numeric(elevated_Sd), elevated_Se = as.numeric(elevated_Se) )
#df_only$exp <- tolower(df_only$prev_name)
#df_only$response <- df_only$Data_type
#2. convert to lower case
#check <- subset(df_only,exp=="biocon_c")%>% group_by(response)  %>% summarise(number = n())

