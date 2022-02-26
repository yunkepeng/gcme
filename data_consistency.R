#soil mineral N, LAI, measured Ci, nmass (narea), LMA, vcmax, jmax, anpp, bnpp, npp
#compare mean_value, n_plot first -> if in trouble --> compare units
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
  dplyr::filter(Unit %in% c("g","g_100g","g_g","g_kg","mg_g","mg_kg","mmol_g","ug_mg","umol_g"))
nmass$response <- "nmass"
logr_c_nmass <- as_tibble(response_ratio_v2(subset(nmass,treatment=="c")))

#other vars - all created now
for (i in 1:nrow(varname)) {
  tryCatch({
    varname1 <- varname$response[i]
    df_c <- subset(kevin_othervars,treatment=="c" & response==varname1)
    assign(paste("logr_c_", varname1,sep=""), as_tibble(response_ratio_v2(df_c)))
    }, error=function(e){})} 

#anpp filtering with unit
logr_c_anpp <-  logr_c_anpp %>% filter(Unit %in% c("g_m2d","g_m2gs","g_m2y","gc_m2y","kg_ha_y","kg_m2y","kgc_m2y","mg_ha_y","t_ha_y"))

#something wrong with soil N - just convert
logr_cw_soil_n$logr_var[(logr_cw_soil_n$exp)=="ua2007_cw"] <- NA

#final, vcmax, jmax and jmax/vcmax
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
logr_c_jmax <- as_tibble(response_ratio_v2(subset(kevin_z,treatment=="c"&response=="jmax")))

#old data processing below -> see overall_traits.R

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

#anpp - just help to fill 3 plots from popface
df$response[df$Data_type=="ANPP"] <- "ANPP"

#anpp - just help to fill 3 plots from popface
df$response[df$Data_type=="LAI"] <- "LAI"
df$response[df$Data_type=="max_LAI"] <- "LAI"
df$response[df$Data_type=="maximum_LAI"] <- "LAI"

df_only <- subset(df,is.na(response)==FALSE)

#newly added soil N - copied from beni
selvars <- c("mineral_soil_N", "soil_mineral_N", "soil_NH4-N", "soil_NO3-N","soil_NO3-N_", "soil_solution_mineral_N", "soil_solution_NH4+", "soil_solution_NO3-", "soil_NO3-_", "soil_NO3-", "soil_NH4+", "soil_NH4+-N", "resin_N", "resin_NH4+", "resin_NO3-", "rhizosphere_NO3-N", "rhizosphere_NH4-N")

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
  }, error=function(e){})}
# never used old_logr_c_ANPP - since it adds more data from logr_c_anpp!


#check 5 variables: anpp, bnpp, nmass, LAI and soil inorganic N
new_old_anpp <- merge(logr_c_anpp%>% group_by(exp)%>%summarise(number = n(),mean_run =mean(logr,na.rm=TRUE)),
      old_logr_c_ANPP%>% group_by(exp)%>%summarise(number = n(),mean_run =mean(logr,na.rm=TRUE)),
      by=c("exp"),all.y=TRUE)
names(new_old_anpp) <- c("exp","new_dataset_no","new_dataset_mean","old_dataset_no","old_dataset_mean")
new_old_anpp

new_old_bnpp <- merge(logr_c_bnpp%>% group_by(exp)%>%summarise(number = n(),mean_run =mean(logr,na.rm=TRUE)),
                      old_logr_c_BNPP%>% group_by(exp)%>%summarise(number = n(),mean_run =mean(logr,na.rm=TRUE)),
                      by=c("exp"),all.y=TRUE)
names(new_old_bnpp) <- c("exp","new_dataset_no","new_dataset_mean","old_dataset_no","old_dataset_mean")
new_old_bnpp

new_old_lai <- merge(logr_c_lai%>% group_by(exp)%>%summarise(number = n(),mean_run =mean(logr,na.rm=TRUE)),
                      old_logr_c_LAI%>% group_by(exp)%>%summarise(number = n(),mean_run =mean(logr,na.rm=TRUE)),
                      by=c("exp"),all.y=TRUE)
names(new_old_lai) <- c("exp","new_dataset_no","new_dataset_mean","old_dataset_no","old_dataset_mean")
new_old_lai

#check vcmax
#final, vcmax, jmax and jmax/vcmax
kevin <- read.csv("/Users/yunpeng/data/gcme/kevin_20220222/MESI_2022.csv")
kevin <- subset(kevin,response=="vcmax")
kevin$x_c <- as.numeric(kevin$x_c)
aa <- kevin%>% group_by(site)  %>% summarise(number = n(),mean_run =mean(x_c,na.rm=TRUE))

kevin_old <- read.csv("/Users/yunpeng/data/gcme/kevin/orig_vcmax/JunkePeng_11252021.csv")
kevin_old <- subset(kevin_old,response=="vcmax")

bb <- kevin_old%>% group_by(site)  %>% summarise(number = n(),mean_run =mean(x_c,na.rm=TRUE))
check <-merge(aa,bb,,by=c("site"),all.x=TRUE)s