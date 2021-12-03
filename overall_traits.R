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

#a test of function to unify all important steps: samples -> per plot -> all plots
rm(list=ls())
#1. response_ratio includes metafor::escalc() to calculate response ratio 
response_ratio <- function(df,variable_name){
  #-----------------------------------------------------------------------
  # Input: 
  # name: df, variable_name,..
  # df should include: c("exp_nam","ambient","ambient_Sd","ambient_Se",
  # "elevated","elevated_Se","elevated_Se","elevated_Sd","n_plots")
  #-----------------------------------------------------------------------
  df_sub <- df %>%         
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
  
  df_sub$id <- 1:nrow(df_sub)
  df_sub$my_varnam <- variable_name
  return(df_sub)
  #-----------------------------------------------------------------------
  # Output: df_sub
  #-----------------------------------------------------------------------
}

#2. agg_plot includes MAd::agg to aggregate samples into plot mean
#NOT used
#since agg_meta was done - not used this because the aggregation was not useful after discussing with beni
agg_plot <- function(df_sub,variable_name){
  #-----------------------------------------------------------------------
  # Input: 
  # name: df_sub that output from response_ratio, variable_name.
  #-----------------------------------------------------------------------
  df_agg <- df_sub %>% 
    filter(!is.na(logr_var) & !is.na(logr)) %>% 
    dplyr::select(-id) %>% # this is necessary because the dumb function agg() uses my 'id' otherwise
    mutate( id = paste(exp_nam, my_varnam, sep="_XXX_")) %>% 
    MAd::agg( id = id, es = logr, var = logr_var, n.1 = n_plots, n.2 = n_plots, cor = 1.0, method = "BHHR", data = . ) %>% 
    as_tibble() %>% 
    mutate( id = str_split(id, "_XXX_") ) %>% 
    mutate( exp_nam = purrr::map_chr(id, 1),
            my_varnam = purrr::map_chr(id, 2) ) %>% 
    dplyr::select(exp_nam, my_varnam, es, var) %>% 
    
    ## add number of plots column and my_varnam
    left_join( df_sub %>% 
                 group_by( exp_nam, my_varnam ) %>%
                 summarise( n_plots = sum(n_plots) ) %>% 
                 dplyr::select( exp_nam, my_varnam, n_plots ),
               by = c("exp_nam", "my_varnam") ) %>% 
    rename( logr = es, logr_var = var ) %>% 
    mutate( logr_se = sqrt(logr_var)/sqrt(n_plots) )
  return(df_agg)
  #-----------------------------------------------------------------------
  # Output: df_agg 
  #-----------------------------------------------------------------------
}

#3. an overall boxplot of vcmax25-co2: based on each plot - needs further work
#create a function that input original df and plot-mean df and output a metabox
#change to no year
agg_meta_noyear <- function(df){
  explist <- unique(df$exp_nam)
  mylist <- list() #create an empty list
  #some are just all NA for original data's standard deviation - needs division -> in this way just calculate normal mean, lower and upper band
  for (i in 1:length(explist)){
    if (all(is.na(subset(df,exp_nam==explist[i])$logr_var)) == TRUE){
      uncertainty <- qt(.975,9)*sd(as.vector(subset(df,exp_nam==explist[i])$logr))/sqrt(length(as.vector(subset(df,exp_nam==explist[i])$logr)))
      mean_value <- mean(subset(df,exp_nam==explist[i])$logr,na.rm=TRUE)
      df_box <- tibble(
        exp_nam=explist[i], middle = mean_value,
        ymin   = mean_value-uncertainty,
        ymax   = mean_value+uncertainty,
        middle_scaled = NA, ymin_scaled   = NA, ymax_scaled   = NA,variance_info ="No")
    } else {# define if uncertainty is 0, when N = 1
      if (length(subset(df,exp_nam==explist[i])$logr_var) ==1){
        uncertainty <- 0
        mean_value <- mean(subset(df,exp_nam==explist[i])$logr,na.rm=TRUE)
        df_box <- tibble(
          exp_nam=explist[i], middle = mean_value,
          ymin   = mean_value-uncertainty,
          ymax   = mean_value+uncertainty,
          middle_scaled = NA, ymin_scaled   = NA, ymax_scaled   = NA,variance_info ="No")
      } else 
        #where year is random factor
        out_meta <- metafor::rma.mv( logr, logr_var, method = "REML",
                                      control = list(stepadj=0.3),data = subset(df,exp_nam==explist[i]))
      out_meta_scaled <- predict( out_meta, transf=exp )
      df_box <- tibble(
        exp_nam=explist[i], middle = out_meta$b[1,1], ymin   = out_meta$ci.lb, ymax   = out_meta$ci.ub,
        middle_scaled = out_meta_scaled$pred, ymin_scaled   = out_meta_scaled$ci.lb, ymax_scaled   = out_meta_scaled$ci.ub,variance_info ="Yes")
    }
    mylist[[i]] <- df_box}
  output <- do.call("rbind",mylist)
  return(output)
}
agg_meta_plots <- function(df,type_name,logr,log_var){
  explist <- unique(eval(parse(text=paste("df$", type_name, sep = ""))))
  mylist <- list() #create an empty list
  df$logr <- eval(parse(text=paste("df$", logr, sep = "")))
  df$logr_var <- eval(parse(text=paste("df$", log_var, sep = "")))
  for (i in 1:length(explist)){
    out_meta <- metafor::rma.mv( logr, logr_var, method = "REML", control = list(stepadj=0.3),
                                 data= subset(df,eval(parse(text=paste(type_name)))==explist[i]))
    df_box <- tibble(
      type_name=explist[i], middle = out_meta$b[1,1], ymin   = out_meta$ci.lb, ymax   = out_meta$ci.ub)
    mylist[[i]] <- df_box
  }
  output <- do.call("rbind",mylist)
  return(output)
}

#above is just a test - now starts

df <- read_csv("~/data/gcme/data_received_190325/NewData_wide_CORRECTED2.csv") %>%
  mutate( ambient_Sd  = as.numeric(ambient_Sd),  ambient_Se  = as.numeric(ambient_Se), 
          elevated_Sd = as.numeric(elevated_Sd), elevated_Se = as.numeric(elevated_Se) )
df <- subset(df,treatment=="f")
dim(df)
#now, shifting to paralleled measurements of LMA, Narea, Nmass, soil N, ANPP, NPP, biomass
LMA <- subset(df,Data_type=="leaf_mass_per_area"|Data_type=="leaf_mass_per_unit_leaf_area"|Data_type=="LMA")
SLA <- subset(df, Data_type=="specific_leaf_area"|Data_type=="SLA")
#LMA_a <- agg_meta_noyear(response_ratio(LMA,"LMA"))
LMA_all <- agg_plot(response_ratio(LMA,"LMA"),"LMA")

#convert SLA to LMA format 
#SLA$ambient <- 1/SLA$ambient; SLA$ambient_Se <- 1/SLA$ambient_Se; SLA$ambient_Sd <- 1/SLA$ambient_Sd
#SLA$elevated <- 1/SLA$elevated;SLA$elevated_Se <- 1/SLA$elevated_Se;SLA$elevated_Sd <- 1/SLA$elevated_Sd;
#SLA_a <- agg_meta_noyear(response_ratio(SLA,"SLA"))
SLA_all <- agg_plot(response_ratio(SLA,"SLA"),"SLA")

#LMA_all <- rbind(LMA_a,SLA_a)
#LMA_final <- LMA_all[,c("exp_nam","middle")]
#LMA_final <- LMA_all[,c("exp_nam","logr")]
#names(LMA_final) <- c("exp_nam","logLMA")

#leaf C/N
leafCN <- subset(df,Data_type=="leaf_C:N")
#leafCN_all <- agg_meta_noyear(response_ratio(leafCN,"leafCN"))
#leafCN_final <- leafCN_all[,c("exp_nam","middle")]
leafCN_all <- agg_plot(response_ratio(leafCN,"leafCN"),"leafCN")
leafCN_final <- leafCN_all[,c("exp_nam","logr")]
names(leafCN_final) <- c("exp_nam","logleafCN")

#Narea,Nmass
leaf_N <- subset(df,Data_type=="leaf_N")
leafNmass <-  subset(leaf_N,Unit=="%"|Unit=="g/g"|Unit=="g/kg"|Unit=="mg_/_g"|Unit=="mg/g"|Unit=="g/mg"|
                       Unit=="mg/g_"|Unit=="mg/kg"|Unit=="mol/g"|Unit=="mmol/g")
leafNarea <-  subset(leaf_N,Unit!="%" & Unit!="g/g" & Unit!="g/kg" & Unit!="mg_/_g" & Unit!="mg/g" & Unit!="g/mg" &
                       Unit!="mg/g_" & Unit!="mg/kg" & Unit!="mol/g" & Unit!="mmol/g")
#Nmass_all <- agg_meta_noyear(response_ratio(leafNmass,"Nmass"))
#Narea_all <- agg_meta_noyear(response_ratio(leafNarea,"Narea"))
#Nmass_final <- Nmass_all[,c("exp_nam","middle")]
#Narea_final <- Narea_all[,c("exp_nam","middle")]

Nmass_all <- agg_plot(response_ratio(leafNmass,"Nmass"),"Nmass")
Nmass_all_final <- Nmass_all[,c("exp_nam","logr")]

Narea_all <- agg_plot(response_ratio(leafNarea,"Narea"),"Narea")
Narea_all_final <- Narea_all[,c("exp_nam","logr")]

#names(Nmass_final) <- c("exp_nam","logNmass")
#names(Narea_final) <- c("exp_nam","logNarea")


#ANPP
ANPP_all <- subset(df,Data_type=="ANPP")
#ANPP_all <- agg_meta_noyear(response_ratio(ANPP_all,"ANPP"))
#ANPP_final <- ANPP_all[,c("exp_nam","middle")]
ANPP_all <- agg_plot(response_ratio(ANPP_all,"ANPP"),"ANPP")
ANPP_all_final <- ANPP_all[,c("exp_nam","logr")]
names(ANPP_all_final) <- c("exp_nam","logANPP")

#soil N
#SoilN_all <- agg_meta_noyear(response_ratio(subset(df,Data_type=="soil_N"),"soil_N"))
#soilN_final <- SoilN_all[,c("exp_nam","middle")]
soilN_all <- agg_plot(response_ratio(subset(df,Data_type=="soil_N"),"soil_N"),"soil_N")
soilN_all_final <- soilN_all[,c("exp_nam","logr")]
names(soilN_all_final) <- c("exp_nam","logsoilN")

soil_mineral_N_all <- agg_plot(response_ratio(subset(df,Data_type=="soil_mineral_N"),"soil_mineral_N"),"soil_mineral_N")
soil_mineral_N_all_final <- soil_mineral_N_all[,c("exp_nam","logr")]
names(soil_mineral_N_all_final) <- c("exp_nam","soil_mineral_N")

soil_NH4_N_all <- agg_plot(response_ratio(subset(df,Data_type=="soil_NH4-N"),"soil_NH4"),"soil_NH4")
soil_NH4_N_all_final <- soil_NH4_N_all[,c("exp_nam","logr")]
names(soil_NH4_N_all_final) <- c("exp_nam","soil_NH4")


soil_NO3_N_all <- agg_plot(response_ratio(subset(df,Data_type=="soil_NO3-N"),"soil_NO3-N"),"soil_NO3-N")
soil_NO3_N_all_final <- soil_NO3_N_all[,c("exp_nam","logr")]
names(soil_NO3_N_all_final) <- c("exp_nam","soil_NO3_N")

soil_CN_all <- agg_plot(response_ratio(subset(df,Data_type=="mineral_soil_C:N"),"mineral_soil_C:N"),"mineral_soil_C:N")
soil_CN_all_final <- soil_CN_all[,c("exp_nam","logr")]
names(soil_CN_all_final) <- c("exp_nam","soil_NO3_N")

#LAI
LAI_all <- subset(df,Data_type=="LAI"|Data_type=="max_LAI"|Data_type=="maximum_LAI")
#LAI_all <- agg_meta_noyear(response_ratio(LAI_all,"LAI"))
#LAI_all_final <- LAI_all[,c("exp_nam","middle")]
LAI_all <- agg_plot(response_ratio(LAI_all,"LAI"),"LAI")
LAI_all_final <- LAI_all[,c("exp_nam","logr")]
names(ANPP_all_final) <- c("exp_nam","logLAI")

#additionally having c13, n15
#check c13
unique((df[grep("13", df$Data_type),])$Data_type)
C13_all <- subset(df,Data_type=="d13C_isotope_ratio_hardwood_leaflitter"|
                    Data_type=="d13C_isotope_ratio_leaflitter"|
                    Data_type=="leaf_d13C"|Data_type=="leaf_litter_d13C"|
                    Data_type=="litter_stable_13C_isotope_signature")
C13_all_positive <- subset(df,Data_type=="leaf_?13C")
#C13_all <- agg_meta_noyear(response_ratio(C13_all,"c13"))
C13_all <- agg_plot(response_ratio(C13_all,"C13"),"C13")
C13_all_final <- C13_all[,c("exp_nam","logr")]
names(C13_all_final) <- c("exp_nam","logC13")

unique((df[grep("15", df$Data_type),])$Data_type)
N15_all <- subset(df,Data_type=="leaf_d15N")
#N15_all <- agg_meta_noyear(response_ratio(N15_all,"N15"))
N15_all <- agg_plot(response_ratio(N15_all,"N15"),"N15")
N15_all_final <- N15_all[,c("exp_nam","logr")]
names(N15_all_final) <- c("exp_nam","logN15")

#
vcmax_all <- subset(df,Data_type=="maximal_Rubisco_activity_(Vcmax)"|Data_type=="Vcmax"|
                      Data_type=="Vcmax,_maximum_carboxylation_rate_of_Rubisco")
vcmax_all <- agg_plot(response_ratio(vcmax_all,"vcmax"),"vcmax")
vcmax_all_final <- vcmax_all[,c("exp_nam","logr")]
names(vcmax_all_final) <- c("exp_nam","vcmax")

Jmax_all <- subset(df,Data_type=="Jmax_"|Data_type=="Jmax"|
                     Data_type=="Jmax,_maximum_rate_of_electron_transport_for_regeneration_of_RuBP")
Jmax_all <- agg_plot(response_ratio(Jmax_all,"jmax"),"jmax")
Jmax_all_final <- Jmax_all[,c("exp_nam","logr")]
names(Jmax_all_final) <- c("exp_nam","jmax")

leafP_all <- subset(df,Data_type=="leaf_P")
leafP_all <- agg_plot(response_ratio(leafP_all,"leaf_P"),"leaf_P")

NPP_all <- subset(df,Data_type=="NPP")
NPP_all <- agg_plot(response_ratio(NPP_all,"NPP"),"NPP")

NPPleaf_all <- subset(df,Data_type=="ANPP_foliage")
NPPleaf_all <- agg_plot(response_ratio(NPPleaf_all,"ANPP_foliage"),"ANPP_foliage")

BNPP_all <- subset(df,Data_type=="BNPP")
BNPP_all <- agg_plot(response_ratio(BNPP_all,"BNPP"),"BNPP")

soil_P_all <- subset(df,Data_type=="soil_P")
soil_P_all <- agg_plot(response_ratio(soil_P_all,"soil_P"),"soil_P")


#now, combine
a0 <- agg_meta_plots(SLA_all,"my_varnam","logr","logr_var")
a1 <- agg_meta_plots(LMA_all,"my_varnam","logr","logr_var")
a2 <- agg_meta_plots(leafCN_all,"my_varnam","logr","logr_var")
a3 <- agg_meta_plots(Narea_all,"my_varnam","logr","logr_var")
a4 <- agg_meta_plots(ANPP_all,"my_varnam","logr","logr_var")
a5 <- agg_meta_plots(C13_all,"my_varnam","logr","logr_var")
a6 <- agg_meta_plots(N15_all,"my_varnam","logr","logr_var")
a7 <- agg_meta_plots(vcmax_all,"my_varnam","logr","logr_var")
a8 <- agg_meta_plots(Jmax_all,"my_varnam","logr","logr_var")
a9 <- agg_meta_plots(soilN_all,"my_varnam","logr","logr_var")
a10 <- agg_meta_plots(soil_mineral_N_all,"my_varnam","logr","logr_var")
a11 <- agg_meta_plots(soil_NH4_N_all,"my_varnam","logr","logr_var")
a12 <- agg_meta_plots(soil_NO3_N_all,"my_varnam","logr","logr_var")
a13 <- agg_meta_plots(soil_CN_all,"my_varnam","logr","logr_var")
a14 <- agg_meta_plots(soil_mineral_N_all,"my_varnam","logr","logr_var")
a15 <- agg_meta_plots(leafP_all,"my_varnam","logr","logr_var")
a16 <- agg_meta_plots(NPP_all,"my_varnam","logr","logr_var")
a17 <- agg_meta_plots(NPPleaf_all_all,"my_varnam","logr","logr_var")
a18 <- agg_meta_plots(BNPP_all,"my_varnam","logr","logr_var")
a19 <- agg_meta_plots(soil_P_all,"my_varnam","logr","logr_var")
a20 <- agg_meta_plots(LAI_all,"my_varnam","logr","logr_var")

#a_final <- rbind(a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20)
#co2: a_final <- rbind(a0,a1,a2,a3,a4,a5,a7,a8,a9,a10,a11,a12,a13,a14,a16,a18,a20)
#warming: a_final <- rbind(a2,a3,a4,a7,a8,a9,a11,a12,a16,a18,a20)
#fertilization: a_final <- rbind(a0,a1,a3,a4,a7,a8,a9,a10,a11,a12,a13,a14,a15,a16,a18,a19,a20)

a_final %>%
  ggplot( aes(x=type_name, y=middle)) +
  geom_crossbar(aes(x=type_name, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5 ) +
  geom_hline( yintercept=0.0, size=0.5 )+
  labs(x="", y="N fertilization effect", size=expression(paste("Error"^{-1}))) +
  coord_flip() 
ggsave(paste("~/data/output_gcme/colin/overall_fertilization.jpg",sep=""),width = 10, height = 5)
