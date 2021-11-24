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

#a test of function to unify all important steps: samples -> per plot -> all plots
final_co2 <- read.csv("~/data/example_df.csv")
#final_co2<- final_co2[,c("exp_nam","ambient","ambient_Sd","ambient_Se","elevated","elevated_Se","elevated_Se","elevated_Sd","n_plots")]

#1. response_ratio includes metafor::escalc() to calculate response ratio 
response_ratio <- function(df,variable_name){
  #-----------------------------------------------------------------------
  # Input: 
  # name: df, variable_name,..
  # df should include: c("exp_nam","ambient","ambient_Sd","ambient_Se",
  # "elevated","elevated_Se","elevated_Se","elevated_Sd","n_plots")
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
  df_c_sub$my_varnam <- variable_name
  return(df_c_sub)
  #-----------------------------------------------------------------------
  # Output: df_c_sub
  #-----------------------------------------------------------------------
}
aa <- response_ratio(final_co2,"vcmax25")

#2. agg_plot includes MAd::agg to aggregate samples into plot mean
#NOT used
#since agg_meta was done - not used this because the aggregation was not useful after discussing with beni
agg_plot <- function(df_c_sub,variable_name){
  #-----------------------------------------------------------------------
  # Input: 
  # name: df_c_sub that output from response_ratio, variable_name.
  #-----------------------------------------------------------------------
  df_c_agg <- df_c_sub %>% 
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
    left_join( df_c_sub %>% 
                 group_by( exp_nam, my_varnam ) %>%
                 summarise( n_plots = sum(n_plots) ) %>% 
                 dplyr::select( exp_nam, my_varnam, n_plots ),
               by = c("exp_nam", "my_varnam") ) %>% 
    rename( logr = es, logr_var = var ) %>% 
    mutate( logr_se = sqrt(logr_var)/sqrt(n_plots) )
  return(df_c_agg)
  #-----------------------------------------------------------------------
  # Output: df_c_agg 
  #-----------------------------------------------------------------------
}
bb <- agg_plot(aa,"vcmax25")

#3. an overall boxplot of vcmax25-co2: based on each plot - needs further work
#create a function that input original df and plot-mean df and output a metabox

agg_meta <- function(df){
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
                                     random = ~ 1 | Year, slab = Year, control = list(stepadj=0.3),data = subset(df,exp_nam==explist[i]))
      out_meta_scaled <- predict( out_meta, transf=exp )
      df_box <- tibble(
        exp_nam=explist[i], middle = out_meta$b[1,1], ymin   = out_meta$ci.lb, ymax   = out_meta$ci.ub,
        middle_scaled = out_meta_scaled$pred, ymin_scaled   = out_meta_scaled$ci.lb, ymax_scaled   = out_meta_scaled$ci.ub,variance_info ="Yes")
    }
  mylist[[i]] <- df_box}
  output <- do.call("rbind",mylist)
  return(output)
}


#do meta-analysis for all each bplot
df_metabox <- agg_meta(aa)
summary(df_metabox$middle)
aa %>%
  ggplot( aes(x=exp_nam, y=logr)) +
  geom_jitter( color = rgb(0,0,0,0.3), aes( size = 1/logr_se ), position = position_jitter(w = 0.2, h = 0) ) +
  geom_boxplot(alpha=0.2,outlier.shape = NA,color="red") + # to compare
  geom_crossbar( data = df_metabox, aes(x=exp_nam, y=middle, ymin=ymin, ymax=ymax), fill = "grey80", alpha = 0.6, width = 0.5 ) +
  geom_hline( yintercept=0.0, size=0.5 ) +
  labs(x="", y="Log Response Ratio", size=expression(paste("Error"^{-1}))) +
  coord_flip() +
  ylim(-1,1) 
# big difference!


#above is just a test - now starts
df <- read_csv("~/data/gcme/data_received_190325/NewData_wide_CORRECTED2.csv") %>%
  mutate( ambient_Sd  = as.numeric(ambient_Sd),  ambient_Se  = as.numeric(ambient_Se), 
          elevated_Sd = as.numeric(elevated_Sd), elevated_Se = as.numeric(elevated_Se) )

# save experiments names
df_experiments <- df %>% dplyr::select(exp_nam, prev_name) %>% distinct()

## Determine all experiments that have more than 1 year data
list_exp_gt1yr <- df %>% filter(!is.na(Year)) %>% group_by(exp_nam) %>% summarise(nyears=max(Year)) %>% filter(nyears>1) %>% dplyr::select(exp_nam) %>% unlist() %>% unname()

#include coordinates and experimental details
coord <- read.csv("~/data/gcme/data_received_190325/table_var_exp_names.csv")[,c("Site.Name","System.Type","Latitude","Lat","Longitude","Lon","Elevation","AM","EcM","ErM")]
for (i in 1:nrow(coord)){
  if (coord$Lat[i] == "S"){
    coord$Latitude[i] <- -as.numeric(coord$Latitude[i])
  } else {
    coord$Latitude[i] <- as.numeric(coord$Latitude[i])
  }
  if (coord$Lon[i] == "W"){
    coord$Longitude[i] <- -as.numeric(coord$Longitude[i])
  } else {
    coord$Longitude[i] <- as.numeric(coord$Longitude[i])
  }
}
coord <- coord[,c("Site.Name","System.Type","Latitude","Longitude","Elevation","AM","EcM","ErM")]
names(coord) <- c("exp_nam","pft","lat","lon","z","AM","EcM","ErM")
coord$rep <- duplicated(coord$exp_nam)
coord <- subset(coord,rep=="FALSE")# remove 11 duplicate site

coord2 <- coord[,c("exp_nam","lat","lon","z")]
csvfile <- paste("~/data/gcme/data_received_190325/table_var_exp_names_coordinates.csv")
write.csv(coord2, csvfile, row.names = TRUE)

df_experiments_coord <- merge(df_experiments,coord,by=c("exp_nam"),all.x=TRUE)

experimental <- read.csv("~/data/gcme/data_received_190325/table_var_exp_names_experiments.csv")[,c("Site.Name","prev_name","Start_Year","End_Year","Cquantity_Control","Cquantity_Treatment",
                                                                                                    "Fumigation_type","Nutrients","Type","N.Quantity.control","N.Quantity.treatment",
                                                                                                    "N.Quantity.unit","P.Quantity.control","P.Quantity.treatment","P.Quantity.unit")]

names(experimental) <- c("exp_nam","prev_name","Start_Year","End_Year","Cquantity_Control","Cquantity_Treatment",
                         "Fumigation_type","Nutrients","Type","N.Quantity.control","N.Quantity.treatment",
                         "N.Quantity.unit","P.Quantity.control","P.Quantity.treatment","P.Quantity.unit")

df_experiments_coord_exp <- merge(df_experiments_coord,experimental,by=c("exp_nam","prev_name"),all.x=TRUE)

check_na <- subset(df_experiments_coord_exp,is.na(lon)==TRUE) %>% dplyr::select(exp_nam) %>% distinct()
dim(check_na)
check_yes <- subset(df_experiments_coord_exp,is.na(lon)==FALSE) %>% dplyr::select(exp_nam) %>% distinct()
dim(check_yes)
#output to site coordinates csv
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/coord_df.csv")
write.csv(df_experiments_coord_exp, csvfile, row.names = TRUE)

#for leaf ci/ca, leaf N, vcmax, jmax, c13
Asat_leafN_c13_all <- subset(df,Data_type=="Ci"|Data_type=="Ci/Ca"|
                               Data_type=="leaf_N"|Data_type=="Asat"|Data_type=="Amax"|
                               Data_type=="CO2_assimilation_rate_(Asat)"|Data_type=="CO2_assimilation_rate_(Amax)"|
                               Data_type=="maximal_Rubisco_activity_(Vcmax)"|Data_type=="Vcmax"|
                               Data_type=="Vcmax,_maximum_carboxylation_rate_of_Rubisco"|Data_type=="Vcmax25"|
                               Data_type=="Jmax"|Data_type=="Jmax25"|Data_type=="Jmax_"|
                               Data_type=="Jmax,_maximum_rate_of_electron_transport_for_regeneration_of_RuBP"|
                               Data_type=="leaf_[C]"|Data_type=="leaf_?13C"|Data_type=="leaf_d13C")

Asat_all2 <- df[grep("photo", df$Data_type),] # for photosynthesis capacity (Asat)

Asat_all3 <- df[grep("water_use", df$Data_type),] # for ci/ca

Asat_leafN_final <- rbind(Asat_leafN_c13_all,Asat_all2,Asat_all3)

leaf_all_coord <- merge(Asat_leafN_final,df_experiments_coord_exp,by=c("exp_nam","prev_name"),all.x=TRUE)
dim(Asat_leafN_final)
dim(leaf_all_coord)

coord_Asat <- leaf_all_coord %>% dplyr::select(exp_nam,lat,lon,z) %>% distinct()
dim(coord_Asat)
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/coord_Asat_leafN.csv")
write.csv(coord_Asat, csvfile, row.names = TRUE)
#105 sites overall!
#then, after looking at sources from..., this was then filled by 
#/Users/yunpeng/data/gcme/yunke_coord/coord_Asat_leafN_filled.csv
#ingest elevation
leaf_sites <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/coord_Asat_leafN_filled.csv")
leaf_sites <- subset(leaf_sites,is.na(lon)==FALSE)# remove an AUS site without coordinates
dim(leaf_sites)
devtools::load_all("~/yunkepeng/gcme/pmodel/ingestr/")
aaa <- aggregate(leaf_sites,by=list(leaf_sites$lon,leaf_sites$lat), FUN=mean, na.rm=TRUE)
aaa <- aaa[,c("lon","lat")]
aaa$sitename <- paste("c",1:length(aaa$lon),sep="")

df_etopo <- ingest(aaa,source = "etopo1",dir = "~/data/etopo/" )
aaa$elv <- as.numeric(as.data.frame(df_etopo$data))
aaa$elv[aaa$elv<0] <- 0 #two points are negative, change to 0
aaa <- aaa[,c("lon","lat","elv")]
leaf_sites_z <- merge(leaf_sites,aaa,by=c("lon","lat"),all.x=TRUE)
plot(leaf_sites_z$z~leaf_sites_z$elv) # looks ok - now interploate original elevation value with etopo elevation
leaf_sites_z$z[is.na(leaf_sites_z$z)==TRUE] <- leaf_sites_z$elv[is.na(leaf_sites_z$z)==TRUE] 
leaf_sites_z_final <- leaf_sites_z[,c("exp_nam","lon","lat","z")]

df_leafsites <- merge(df,leaf_sites_z_final,by=c("exp_nam"),all.x=TRUE)
dim(df_leafsites)

#some years were missing, interpolating them primarily basing on max measurement year from the same plot and same variable
# secondarily basing on max measurement year from the same plot
NA_years <- as.data.frame(subset(df_leafsites,is.na(Sampling_Year)==TRUE)%>% group_by(exp_nam,Data_type)%>% summarise(number = n()))

for (i in 1:nrow(NA_years)){
  exp_name <- NA_years[i,1]
  variable_name <- NA_years[i,2]
  if (is.infinite(max(df_leafsites$Sampling_Year[is.na(df_leafsites$Sampling_Year)==FALSE & df_leafsites$exp_nam==exp_name & df_leafsites$Data_type==variable_name]))==FALSE){
    df_leafsites$Sampling_Year[is.na(df_leafsites$Sampling_Year)==TRUE & df_leafsites$exp_nam==exp_name & df_leafsites$Data_type==variable_name] <-max(df_leafsites$Sampling_Year[is.na(df_leafsites$Sampling_Year)==FALSE & df_leafsites$exp_nam==exp_name & df_leafsites$Data_type==variable_name])
  } else {
    df_leafsites$Sampling_Year[is.na(df_leafsites$Sampling_Year)==TRUE & df_leafsites$exp_nam==exp_name & df_leafsites$Data_type==variable_name] <-max(df_leafsites$Sampling_Year[is.na(df_leafsites$Sampling_Year)==FALSE & df_leafsites$exp_nam==exp_name])
  }
}

#only a few missing - no worries!

#now, start working with vcmax data

vcmax_all <- subset(df_leafsites,Data_type=="maximal_Rubisco_activity_(Vcmax)"|Data_type=="Vcmax"|
                      Data_type=="Vcmax,_maximum_carboxylation_rate_of_Rubisco")
vcmax25_all <- subset(df_leafsites,Data_type=="Vcmax25") # all unit is ok

Jmax_all <- subset(df_leafsites,Data_type=="Jmax_"|Data_type=="Jmax"|
                     Data_type=="Jmax,_maximum_rate_of_electron_transport_for_regeneration_of_RuBP")

Jmax25_all <- subset(df_leafsites,Data_type=="Jmax25" & Unit!="g/m2y") # remove one 'outlier' unit

Asat_all <- df_leafsites[grep("Asat", df_leafsites$Data_type),]

Amax_all <- df_leafsites[grep("Amax", df_leafsites$Data_type),]

#now, convert unit separately
#this three value is consistent, but its value should be umol/m2/s. Let's use this as final data of Asat
summary(subset(Asat_all,Unit=="?mol_CO2_m?2_s?1")$ambient)
summary(subset(Asat_all,Unit=="mol_CO2_/_m2_s")$ambient)
summary(subset(Asat_all,Unit=="mol/m2*s")$ambient)

Asat_all_final <- subset(Asat_all,Unit=="?mol_CO2_m?2_s?1"|
                           Unit=="mol_CO2_/_m2_s"|
                           Unit=="mol/m2*s")

Amax_all_final <- subset(Asat_all,Unit=="?mol_CO2_m?2_s?1"|
                           Unit=="mol_/_m2_s"|
                           Unit=="mol/m2*s")

vcmax_all %>% group_by(Unit)  %>% summarise(mean = mean(ambient,na.rm=TRUE))
summary(subset(vcmax_all,Unit=="mol/m2_x_s")$ambient)
summary(subset(vcmax_all,is.na(Unit)==TRUE)$ambient)

#I think all vcmax is umol/m2/s, in no doubt! So does Jmax
vcmax_all_final <- vcmax_all
vcmax25_all_final <- vcmax25_all

Jmax_all_final <- Jmax_all
Jmax25_all_final <- Jmax25_all

#one point shows start and sampling year was 1998 and 1997 separately, wrong. Convert to same year
Asat_all_final$Sampling_Year[Asat_all_final$Sampling_Year<Asat_all_final$Start_Year] <- 1998
Amax_all_final$Sampling_Year[Amax_all_final$Sampling_Year<Amax_all_final$Start_Year] <- 1998
Asat_all_final$Year <- Asat_all_final$Sampling_Year-Asat_all_final$Start_Year
Amax_all_final$Year <- Amax_all_final$Sampling_Year-Amax_all_final$Start_Year

vcmax_all_final$Start_Year[is.na(vcmax_all_final$Start_Year)==TRUE] <- 1995 # one place was missing - interpolate as equal to sampling_year
vcmax_all_final$Year <- vcmax_all_final$Sampling_Year-vcmax_all_final$Start_Year
Jmax_all_final$Year <- Jmax_all_final$Sampling_Year-Jmax_all_final$Start_Year
vcmax25_all_final$Year <- vcmax25_all_final$Sampling_Year-vcmax25_all_final$Start_Year

Jmax25_all_final$Start_Year[is.na(Jmax25_all_final$Start_Year)==TRUE]  <- 1995 # one place was missing - interpolate as equal to sampling_year
Jmax25_all_final$Year <- Jmax25_all_final$Sampling_Year-Jmax25_all_final$Start_Year

#now, all done! start to do one-point method - see gwr_for_climates_Asat.R
csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/Asat_all_df.csv")
write.csv(Asat_all_final, csvfile, row.names = TRUE)

csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/Amax_all_df.csv")
write.csv(Amax_all_final, csvfile, row.names = TRUE)

csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/vcmax_all_df.csv")
write.csv(vcmax_all_final, csvfile, row.names = TRUE)

csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/jmax_all_df.csv")
write.csv(Jmax_all_final, csvfile, row.names = TRUE)

csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/vcmax25_all_df.csv")
write.csv(vcmax25_all_final, csvfile, row.names = TRUE)

csvfile <- paste("/Users/yunpeng/data/gcme/yunke_coord/jmax25_all_df.csv")
write.csv(Jmax25_all_final, csvfile, row.names = TRUE)

#now, after gwr_for_climates_Asat.R by collecting all climate data for vcmax, jmax, Asat, Amax...
#basing on proposal by below priority
#1. vcmax25
#2. vcmax and then converting to vcmax25 by Tg
#3. Asat and then converting to vcmax by one-point method, and convert to vcmax25 by Tg
#so does jmax..

#1.vcmax25
summary(vcmax25_all_final)
vcmax25_all_final%>% group_by(exp_nam,treatment)  %>% summarise(number = n())
#Brandbjerg, c, d, w --> 7 types

#2.vcmax
summary(vcmax_all_final)
vcmax_all_final%>% group_by(exp_nam,treatment)  %>% summarise(number = n())
#Brandbjerg, c, d, w --> 7 types but more samples --> so still use vcmax rather than vcmax25 in Brandbjerg
#now, include Tg
tmn_df <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/climates/vcmax_tmn.csv")
tmx_df <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/climates/vcmax_tmx.csv")
lat <- tmn_df$lat
tmx <- tmx_df[,c(8:ncol(tmx_df))]
tmn <- tmn_df[,c(8:ncol(tmn_df))]
s1 <- (-23.1+ -17.3)/2
s2 <- (-17.3 + -8)/2
s3 <- (-8 + 4.1)/2
s4 <- (4.1 + 14.8)/2
s5 <- (14.8 + 21.9)/2
s6 <- (21.9 + 23.2)/2
s7 <- (23.2 + 18.3)/2
s8 <- (18.3 + 8.6)/2
s9 <-  (8.6 + -2.8)/2
s10 <- (-2.8 + -14.1)/2
s11 <- (-14.1 + -21.6)/2
s12 <- (-21.6 + -23.1)/2

s <- c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)
xx <- data.frame(matrix(nrow=nrow(tmx), ncol=ncol(tmx)))
output_Tg <- data.frame(matrix(nrow=nrow(tmx), ncol=ncol(tmx)))
#xx = acos(h), h = hour angle of the sun
for (a in 1:12){ 
  month_no <- seq(from = 1, to = ncol(tmx), by = 12)+a-1
  xx[1:nrow(tmx),month_no] <- -tan(pi*lat/180)*tan(s[a]*pi/180)
}

#check each part of Tg formula
part1 <- (0.5+((1-xx^2)^(0.5))/(2*acos(xx)))[,1:12]
part2 <- (0.5-((1-xx^2)^(0.5))/(2*acos(xx)))[,1:12]
summary(part1)
summary(part2)

#the percentage of tmx was dominated overall
Tg_site <- tmx*(0.5+((1-xx^2)^(0.5))/(2*acos(xx)))+ tmn*(0.5-((1-xx^2)^(0.5))/(2*acos(xx)))
Tg_site[Tg_site =="NaN"] <- NA
Tg_site[Tg_site < 0] <- NA

tmn_df$Tg <- rowMeans(Tg_site,na.rm=TRUE)
vcmax_Tg <- tmn_df[,c("lon","lat","z","year_start","year_end","Tg")]

#keep the design same as forcing
vcmax_all_final$year_start <- vcmax_all_final$Start_Year
vcmax_all_final$year_end <- vcmax_all_final$Sampling_Year
vcmax_all_final$year_start[vcmax_all_final$year_end>2016] <- 2007
vcmax_all_final$year_end[vcmax_all_final$year_end>2016] <- 2016

vcmax_all_final_Tg <- merge(vcmax_all_final,vcmax_Tg,by=c("lon","lat","z","year_start","year_end"),all.x=TRUE)
summary(vcmax_all_final_Tg)

#now, calculate response ratio directly
vcmax_all_final_Tg$vcmax25_a <- vcmax_all_final_Tg$ambient*exp((65330/8.314)*((1/(vcmax_all_final_Tg$Tg+273.15))-(1/298.15)))
vcmax_all_final_Tg$vcmax25_e <- vcmax_all_final_Tg$elevated*exp((65330/8.314)*((1/(vcmax_all_final_Tg$Tg+273.15))-(1/298.15)))
vcmax_all_final_Tg$response_ratio <- log(vcmax_all_final_Tg$vcmax25_e/vcmax_all_final_Tg$vcmax25_a)

#use alternative way (kattge and Knorr, 2007) - assume Tcgrowth = Tleaf = Tg
tcgrowth <- vcmax_all_final_Tg$Tg + 273.15
tkleaf <- vcmax_all_final_Tg$Tg + 273.15
dent <-  668.39-1.07*tcgrowth
fva <- exp( 71513 * (tkleaf - 298.15) / (298.15 * 8.3145 * tkleaf) )
fvb <- (1.0 + exp( (298.15 * dent - 200000)/(8.3145 * 298.15) ) ) / (1.0 + exp( (tkleaf * dent - 200000)/(8.3145 * tkleaf) ) )
fv  <- fva * fvb
vcmax_all_final_Tg$vcmax25_a_inst <- vcmax_all_final_Tg$ambient/fv
vcmax_all_final_Tg$vcmax25_e_inst <- vcmax_all_final_Tg$elevated/fv
vcmax_all_final_Tg$response_ratio_inst <- log(vcmax_all_final_Tg$vcmax25_e_inst/vcmax_all_final_Tg$vcmax25_a_inst)
vcmax_all_final_Tg$vcmax_a <- vcmax_all_final_Tg$ambient
vcmax_all_final_Tg$vcmax_e <- vcmax_all_final_Tg$elevated
#nothing changed here...just update the formula..
hist(vcmax_all_final_Tg$response_ratio)


#now, final combination
vcmax_all_final_Tg$method <- "vcmax"

#finally, check unique basing on exp_nam, prev_name,ambient,elevated - removing them
vcmax25_final2 <- vcmax_all_final_Tg[!duplicated(vcmax_all_final_Tg[c("exp_nam","prev_name","ambient","elevated")]),]

#now, firstly look at co2 compulsory effect
vcmax25_final2%>% group_by(treatment)  %>% summarise(number = n())
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="c"] <- "co2"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="cd"] <- "co2 + drought"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="cf"] <- "co2 + fertilization"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="ci"] <- "co2 + irrigation"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="cw"] <- "co2+ 1-degree warming"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="cw3"] <- "co2+ 3-degree warming"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="cwd"] <- "co2 + warming + drought"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="d"] <- "drought"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="f"] <- "fertilization"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="i"] <- "irrigation"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="w"] <- "1-degree warming"
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="w2"] <- " 2-degree warming "
vcmax25_final2$treatment_label[vcmax25_final2$treatment=="wd"] <- "warming + drought"

#now, jmax!!
#1.jmax25
summary(Jmax25_all_final)
Jmax25_all_final%>% group_by(exp_nam,treatment)  %>% summarise(number = n())
#3 sites, 12 types

#2.jmax
summary(Jmax_all_final)
Jmax_all_final%>% group_by(exp_nam,treatment)  %>% summarise(number = n())
#lots of samples than jmax25 - so use jmax at the end
#now, include Tg
tmn_df <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/climates/jmax_tmn.csv")
tmx_df <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/climates/jmax_tmx.csv")
lat <- tmn_df$lat
tmx <- tmx_df[,c(8:ncol(tmx_df))]
tmn <- tmn_df[,c(8:ncol(tmn_df))]
s1 <- (-23.1+ -17.3)/2
s2 <- (-17.3 + -8)/2
s3 <- (-8 + 4.1)/2
s4 <- (4.1 + 14.8)/2
s5 <- (14.8 + 21.9)/2
s6 <- (21.9 + 23.2)/2
s7 <- (23.2 + 18.3)/2
s8 <- (18.3 + 8.6)/2
s9 <-  (8.6 + -2.8)/2
s10 <- (-2.8 + -14.1)/2
s11 <- (-14.1 + -21.6)/2
s12 <- (-21.6 + -23.1)/2

s <- c(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12)
xx <- data.frame(matrix(nrow=nrow(tmx), ncol=ncol(tmx)))
output_Tg <- data.frame(matrix(nrow=nrow(tmx), ncol=ncol(tmx)))
#xx = acos(h), h = hour angle of the sun
for (a in 1:12){ 
  month_no <- seq(from = 1, to = ncol(tmx), by = 12)+a-1
  xx[1:nrow(tmx),month_no] <- -tan(pi*lat/180)*tan(s[a]*pi/180)
}

#check each part of Tg formula
part1 <- (0.5+((1-xx^2)^(0.5))/(2*acos(xx)))[,1:12]
part2 <- (0.5-((1-xx^2)^(0.5))/(2*acos(xx)))[,1:12]
summary(part1)
summary(part2)

#the percentage of tmx was dominated overall
Tg_site <- tmx*(0.5+((1-xx^2)^(0.5))/(2*acos(xx)))+ tmn*(0.5-((1-xx^2)^(0.5))/(2*acos(xx)))
Tg_site[Tg_site =="NaN"] <- NA
Tg_site[Tg_site < 0] <- NA

tmn_df$Tg <- rowMeans(Tg_site,na.rm=TRUE)
Jmax_Tg <- tmn_df[,c("lon","lat","z","year_start","year_end","Tg")]

#keep the design same as forcing
Jmax_all_final$year_start <- Jmax_all_final$Start_Year
Jmax_all_final$year_end <- Jmax_all_final$Sampling_Year
Jmax_all_final$year_start[Jmax_all_final$year_end>2016] <- 2007
Jmax_all_final$year_end[Jmax_all_final$year_end>2016] <- 2016

Jmax_all_final_Tg <- merge(Jmax_all_final,Jmax_Tg,by=c("lon","lat","z","year_start","year_end"),all.x=TRUE)
summary(Jmax_all_final_Tg)

#now, calculate response ratio directly
Jmax_all_final_Tg$Jmax25_a <- Jmax_all_final_Tg$ambient*exp((43540/8.314)*((1/(Jmax_all_final_Tg$Tg+273.15))-(1/298.15)))
Jmax_all_final_Tg$Jmax25_e <- Jmax_all_final_Tg$elevated*exp((43540/8.314)*((1/(Jmax_all_final_Tg$Tg+273.15))-(1/298.15)))
Jmax_all_final_Tg$response_ratio <- log(Jmax_all_final_Tg$Jmax25_e/Jmax_all_final_Tg$Jmax25_a)

#use alternative way (kattge and Knorr, 2007) - assume Tcgrowth = Tleaf = Tg
tcgrowth <- Jmax_all_final_Tg$Tg + 273.15
tkleaf <- Jmax_all_final_Tg$Tg + 273.15
dent <-  659.7-0.7*tcgrowth
fva <- exp( 49884 * (tkleaf - 298.15) / (298.15 * 8.3145 * tkleaf) )
fvb <- (1.0 + exp( (298.15 * dent - 200000)/(8.3145 * 298.15) ) ) / (1.0 + exp( (tkleaf * dent - 200000)/(8.3145 * tkleaf) ) )
fv  <- fva * fvb
Jmax_all_final_Tg$jmax25_a_inst <- Jmax_all_final_Tg$ambient/fv
Jmax_all_final_Tg$jmax25_e_inst <- Jmax_all_final_Tg$elevated/fv
Jmax_all_final_Tg$response_ratio_inst <- log(Jmax_all_final_Tg$jmax25_e_inst/Jmax_all_final_Tg$jmax25_a_inst)

Jmax_all_final_Tg$Jmax_a <- Jmax_all_final_Tg$ambient
Jmax_all_final_Tg$Jmax_e <- Jmax_all_final_Tg$elevated
Jmax_all_final_Tg$method <- "jmax"

Jmax25_final2 <- Jmax_all_final_Tg

#now, firstly look at co2 compulsory effect
Jmax25_final2%>% group_by(treatment)  %>% summarise(number = n())
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="c"] <- "co2"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="cd"] <- "co2 + drought"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="cf"] <- "co2 + fertilization"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="ci"] <- "co2 + irrigation"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="cw"] <- "co2+ 1-degree warming"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="cw3"] <- "co2+ 3-degree warming"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="cwd"] <- "co2 + warming + drought"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="d"] <- "drought"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="f"] <- "fertilization"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="i"] <- "irrigation"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="w"] <- "1-degree warming"
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="w2"] <- " 2-degree warming "
Jmax25_final2$treatment_label[Jmax25_final2$treatment=="wd"] <- "warming + drought"

#now, start forcing
#forcing R file was available at ~/yunkepeng/gcme/vcmax_forcing.R
#forcing output was available in /Users/yunpeng/data/gcme/vcmax_forcing/
#siteinfo was available in /Users/yunpeng/data/gcme/yunke_coord/forcing_vcmax.csv
summary(vcmax25_final2)
summary(Jmax25_final2)
#c,cw,cw3,w,w2 were firstly selected to test its response in pmodel
vcmax25_warmingco2 <- subset(vcmax25_final2,treatment=="c"|treatment=="cw"|treatment=="cw3"|treatment=="w"|treatment=="w2")
Jmax25_warmingco2 <- subset(Jmax25_final2,treatment=="c"|treatment=="cw"|treatment=="cw3"|treatment=="w"|treatment=="w2")

subset(vcmax25_final2,treatment=="c"|treatment=="cw"|treatment=="cw3") %>% group_by(co2_a,co2_e)  %>% summarise(number = n())
subset(Jmax25_warmingco2,treatment=="c"|treatment=="cw"|treatment=="cw3") %>% group_by(co2_a,co2_e)  %>% summarise(number = n())

subset(vcmax25_final2,treatment=="w"|treatment=="w2"|treatment=="cw"|treatment=="cw3") %>% group_by(warm_a,warm_e,warmQ_e)  %>% summarise(number = n())
subset(Jmax25_warmingco2,treatment=="w"|treatment=="w2"|treatment=="cw"|treatment=="cw3") %>% group_by(warm_a,warm_e,warmQ_e)  %>% summarise(number = n())
#using warm_e, about how many temperature increased

#combine with siteinfo
siteinfo <- read.csv("/Users/yunpeng/data/gcme/yunke_coord/forcing_vcmax.csv")
siteinfo <- siteinfo[,c("lon","lat","elv","year_start","year_end","sitename")]
names(siteinfo) <- c("lon","lat","z","year_start","year_end","sitename")
vcmax25_warmingco2_siteinfo <- merge(vcmax25_warmingco2,siteinfo,by=c("lon","lat","z","year_start","year_end"),all.x=TRUE)
Jmax25_warmingco2_siteinfo <- merge(Jmax25_warmingco2,siteinfo,by=c("lon","lat","z","year_start","year_end"),all.x=TRUE)

#now, start p-model
#devtools::load_all("~/yunkepeng/gcme/pmodel/ingestr/")

#default inst_vcmax
devtools::load_all("/Users/yunpeng/yunkepeng/gcme/pmodel/rsofun/")
#changing to use adjusted parameters for inst_vcmax from Smith and Keenan 2020 GCB - doesn't change anything!
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

forcing_path <- "/Users/yunpeng/data/gcme/vcmax_forcing/"
vcmax25_warmingco2_siteinfo$elv <- vcmax25_warmingco2_siteinfo$z


for (i in 1:nrow(vcmax25_warmingco2_siteinfo)){
  siteinfo_site <- vcmax25_warmingco2_siteinfo[i,c("sitename","lon","lat","elv","year_start","year_end")]
  sitename <- siteinfo_site$sitename
  
  df1 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for ambient
  df2 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for elevated
  
  #Q1 - add warming directly
  df2$temp <- df2$temp+vcmax25_warmingco2_siteinfo$warm_e[i]
  
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
  
  mean_vcmax25_ambient <- mean(modlist1$vcmax25)*1000000
  mean_vcmax25_elevated <- mean(modlist2$vcmax25)*1000000
  mean_vcmax_ambient <- mean(modlist1$vcmax)*1000000
  mean_vcmax_elevated <- mean(modlist2$vcmax)*1000000
  
  #mean_vcmax25_ambient <- max(modlist1$vcmax25)*1000000
  #mean_vcmax25_elevated <- max(modlist2$vcmax25)*1000000
  
  vcmax25_warmingco2_siteinfo$mean_vcmax25_ambient[i] <- mean_vcmax25_ambient
  vcmax25_warmingco2_siteinfo$mean_vcmax25_elevated[i] <- mean_vcmax25_elevated
  vcmax25_warmingco2_siteinfo$mean_vcmax_ambient[i] <- mean_vcmax_ambient
  vcmax25_warmingco2_siteinfo$mean_vcmax_elevated[i] <- mean_vcmax_elevated
  vcmax25_warmingco2_siteinfo$pred_response_ratio[i] <- log(vcmax25_warmingco2_siteinfo$mean_vcmax25_elevated[i]/vcmax25_warmingco2_siteinfo$mean_vcmax25_ambient[i])
  vcmax25_warmingco2_siteinfo$PPFD[i] <- mean(df1$ppfd,na.rm=TRUE)*1000000 # in umol/m2/s
  vcmax25_warmingco2_siteinfo$vpd[i] <- mean(df1$vpd,na.rm=TRUE)/1000 #in kPa
  vcmax25_warmingco2_siteinfo$temp[i] <- mean(df1$temp,na.rm=TRUE) #in degree celcius
}
summary(vcmax25_warmingco2_siteinfo$pred_response_ratio)

# do the same for jmax
Jmax25_warmingco2_siteinfo$elv <- Jmax25_warmingco2_siteinfo$z

for (i in 1:nrow(Jmax25_warmingco2_siteinfo)){
  siteinfo_site <- Jmax25_warmingco2_siteinfo[i,c("sitename","lon","lat","elv","year_start","year_end")]
  sitename <- siteinfo_site$sitename
  
  df1 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for ambient
  df2 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for elevated
  
  #Q1 - add warming directly
  df2$temp <- df2$temp+Jmax25_warmingco2_siteinfo$warm_e[i]
  
  #Q2 - if having co2 then replacing it directly - if not then using co2 at measurement year
  if (is.na(Jmax25_warmingco2_siteinfo$co2_a[i])==TRUE|is.na(Jmax25_warmingco2_siteinfo$co2_e[i])==TRUE) { 
    df1$co2 <- df1$co2
    df2$co2 <- df2$co2
  } else {
    df1$co2 <- Jmax25_warmingco2_siteinfo$co2_a[i] # set to ambient 
    df2$co2 <- Jmax25_warmingco2_siteinfo$co2_e[i] # set to elevated 
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
  
  mean_jmax25_ambient <- mean(modlist1$jmax25)*1000000
  mean_jmax25_elevated <- mean(modlist2$jmax25)*1000000
  mean_jmax_ambient <- mean(modlist1$jmax)*1000000
  mean_jmax_elevated <- mean(modlist2$jmax)*1000000
  #mean_jmax25_ambient <- max(modlist1$jmax25)*1000000
  #mean_jmax25_elevated <- max(modlist2$jmax25)*1000000
  
  Jmax25_warmingco2_siteinfo$mean_jmax25_ambient[i] <- mean_jmax25_ambient
  Jmax25_warmingco2_siteinfo$mean_jmax25_elevated[i] <- mean_jmax25_elevated
  Jmax25_warmingco2_siteinfo$mean_jmax_ambient[i] <- mean_jmax_ambient
  Jmax25_warmingco2_siteinfo$mean_jmax_elevated[i] <- mean_jmax_elevated
  Jmax25_warmingco2_siteinfo$pred_response_ratio[i] <- log(Jmax25_warmingco2_siteinfo$mean_jmax25_elevated[i]/Jmax25_warmingco2_siteinfo$mean_jmax25_ambient[i])
  Jmax25_warmingco2_siteinfo$PPFD[i] <- mean(df1$ppfd,na.rm=TRUE)*1000000 # in umol/m2/s
  Jmax25_warmingco2_siteinfo$vpd[i] <- mean(df1$vpd,na.rm=TRUE)/1000 #in kPa
  Jmax25_warmingco2_siteinfo$temp[i] <- mean(df1$temp,na.rm=TRUE) #in degree celcius
}

#now, manually merged with ECM/AM plots
naming_expnam <- vcmax25_warmingco2_siteinfo %>% group_by(exp_nam) %>% summarise(number=n()) 
naming_expnam
#csvfile <- paste("~/data/gcme/data_received_190325/co2_plot.csv")
#write.csv(naming_expnam, csvfile, row.names = TRUE)
#newly update sitename, so that it can be merged with ECM/AM info properly
naming_expnam_filled <- read.csv("/Users/yunpeng/data/gcme/data_received_190325/co2_plot_filled.csv")
naming_expnam_filled <- naming_expnam_filled[,c("exp_nam_update","exp_nam")]
names(naming_expnam_filled) <- c("Site.Name","exp_nam")
coord_org <- read.csv("~/data/gcme/data_received_190325/table_var_exp_names.csv")
df_experiments_coord <- merge(coord_org,naming_expnam_filled,by=c("Site.Name"),all.x=TRUE)
df_experiments_coord_co2 <- subset(df_experiments_coord,is.na(exp_nam)==FALSE)
dim(df_experiments_coord_co2)
#now, divide them in advance
df_experiments_coord_co2$ecm_type <- "NA"
df_experiments_coord_co2$ecm_type[df_experiments_coord_co2$EcM==TRUE] <- "ECM"
df_experiments_coord_co2$ecm_type[df_experiments_coord_co2$AM==TRUE] <- "AM"
df_experiments_coord_co2$ecm_type[df_experiments_coord_co2$ErM==TRUE] <- "ErM"
df_experiments_coord_co2$ecm_type[df_experiments_coord_co2$Moss==TRUE] <- "Moss"

df_experiments_coord_co2$pft_info <- "NA"
df_experiments_coord_co2$pft_info[df_experiments_coord_co2$Fern==TRUE] <-"fern"
df_experiments_coord_co2$pft_info[df_experiments_coord_co2$C3.grass==TRUE] <-"grassland"
df_experiments_coord_co2$pft_info[df_experiments_coord_co2$C4.grass==TRUE] <-"grassland"
df_experiments_coord_co2$pft_info[df_experiments_coord_co2$C3.herb==TRUE] <-"grassland"
df_experiments_coord_co2$pft_info[df_experiments_coord_co2$C4.herb==TRUE] <-"grassland"
df_experiments_coord_co2$pft_info[df_experiments_coord_co2$pft_info=="NA"] <- "forest"

df_experiments_coord_co2$c3c4 <- "NA"
df_experiments_coord_co2$c3c4[df_experiments_coord_co2$C4.grass==TRUE] <-"c4"
df_experiments_coord_co2$c3c4[df_experiments_coord_co2$C4.herb==TRUE] <-"c4"
df_experiments_coord_co2$c3c4[df_experiments_coord_co2$c3c4=="NA"] <- "c3"

#some additionals from Cesar's paper - page 14 - SI- Table 1 - https://static-content.springer.com/esm/art%3A10.1038%2Fs41558-019-0545-2/MediaObjects/41558_2019_545_MOESM1_ESM.pdf
#now, fully checked
df_experiments_coord_co2$ecm_type[df_experiments_coord_co2$exp_nam=="FACTS_II_FACE3_pt"] <-"ECM"
df_experiments_coord_co2$ecm_type[df_experiments_coord_co2$exp_nam=="FACTS_II_FACE4_bp"] <-"ECM"
df_experiments_coord_co2$ecm_type[df_experiments_coord_co2$exp_nam=="EucFACE"] <-"ECM"
df_experiments_coord_co2$ecm_type[df_experiments_coord_co2$exp_nam=="EUROFACE4_pn"] <-"ECM"
df_experiments_coord_co2$ecm_type[df_experiments_coord_co2$exp_nam=="EUROFACE4_pe"] <-"ECM"
df_experiments_coord_co2$ecm_type[df_experiments_coord_co2$exp_nam=="EUROFACE4_pa"] <-"ECM"
df_experiments_coord_co2$ecm_type[df_experiments_coord_co2$exp_nam=="POPFACE_pa"] <-"ECM"
df_experiments_coord_co2$ecm_type[df_experiments_coord_co2$exp_nam=="POPFACE_pe"] <-"ECM"
df_experiments_coord_co2$ecm_type[df_experiments_coord_co2$exp_nam=="POPFACE_pn"] <-"ECM"
df_experiments_coord_co2$ecm_type[df_experiments_coord_co2$exp_nam=="Nevada_Desert_FACE"] <-"AM"
df_experiments_coord_co2$ecm_type[df_experiments_coord_co2$exp_nam=="Brandbjerg"] <-"AM"



#final merge
vcmax25_warmingco2_siteinfo_ecm <- merge(vcmax25_warmingco2_siteinfo,df_experiments_coord_co2,by=c("exp_nam"),all.x=TRUE)
#manually add rest of points basing on Terrer et al. 2016 New Phytol- fig.1
vcmax25_warmingco2_siteinfo_ecm$ecm_type[vcmax25_warmingco2_siteinfo_ecm$exp_nam=="BioCON"] <-"AM"
vcmax25_warmingco2_siteinfo_ecm$pft_info[vcmax25_warmingco2_siteinfo_ecm$exp_nam=="BioCON"] <-"grassland"
vcmax25_warmingco2_siteinfo_ecm$c3c4[vcmax25_warmingco2_siteinfo_ecm$exp_nam=="BioCON"] <-"c3"

# do the same for jmax
Jmax25_warmingco2_siteinfo_ecm <- merge(Jmax25_warmingco2_siteinfo,df_experiments_coord_co2,by=c("exp_nam"),all.x=TRUE)
#manually add rest of points basing on Terrer et al. 2016 New Phytol- fig.1
Jmax25_warmingco2_siteinfo_ecm$ecm_type[Jmax25_warmingco2_siteinfo_ecm$exp_nam=="BioCON"] <-"AM"
Jmax25_warmingco2_siteinfo_ecm$pft_info[Jmax25_warmingco2_siteinfo_ecm$exp_nam=="BioCON"] <-"grassland"
Jmax25_warmingco2_siteinfo_ecm$c3c4[Jmax25_warmingco2_siteinfo_ecm$exp_nam=="BioCON"] <-"c3"

#now, it is time to apply 3 functions!
#firstly, only working on vcmax25
dim(vcmax25_warmingco2_siteinfo_ecm)
dim(Jmax25_warmingco2_siteinfo_ecm)
#not yet divided into co2 and warming...
vc25_data <- response_ratio(vcmax25_warmingco2_siteinfo_ecm,"vcmax25")
co2_vcmax25 <- subset(vc25_data,treatment=="c")
#within same plot and treatment -> now agg plot (basing on logr, logr_var and Year as the only random factor)
vc25_data_sitemean <- agg_meta(co2_vcmax25)
#now, output figure
#treatment divided by c, cw, w, w2, cf, f
##three three below should be the same --> so the response ratio should just be directly calculated as from original data!
#summary(co2_vcmax25$vcmax25_e_inst/co2_vcmax25$vcmax25_a_inst)
#summary(co2_vcmax25$vcmax25_e/co2_vcmax25$vcmax25_a)
#summary(co2_vcmax25$elevated/co2_vcmax25$ambient)

#check ecm/am type is all correct after check cesar's si table
vcmax25_warmingco2_siteinfo_ecm
aa <- vcmax25_warmingco2_siteinfo_ecm %>% group_by(exp_nam,ecm_type) %>% summarise(number=n())
#Something should be changed: AM: Brandbjerg
#Something should be changed: ECM: POPFACE_pa, POPFACE_pe, POPFACE_pn

#fig.1 co2 effect on vcmax25
#not added:   geom_jitter(alpha=0.2, aes(size = logr_se), position = position_jitter(w = 0.2, h = 0) ) +
co2_vcmax25 %>%
  ggplot( aes(x=exp_nam, y=logr)) +
  geom_jitter(alpha=0.2, aes(size = 1/logr_se)) +
  geom_crossbar( data = vc25_data_sitemean, aes(x=exp_nam, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5 ) +
  geom_boxplot(aes(x=exp_nam, y=pred_response_ratio,color=ecm_type),size=1)+ # pft_info,c3c4
  geom_hline( yintercept=0.0, size=0.5 )+
  labs(x="", y="CO2 effects on vcmax25", size=expression(paste("Error"^{-1}))) +
  coord_flip() 
ggsave(paste("~/data/output_gcme/colin/co2_1.jpg",sep=""))

#now, jmax25
Jmax25_data <- response_ratio(Jmax25_warmingco2_siteinfo_ecm,"jmax25")
co2_Jmax25 <- subset(Jmax25_data,treatment=="c")
Jmax25_data_sitemean <- agg_meta(co2_Jmax25)
#now, output figure
#treatment divided by c, cw, w, w2, cf, f

#fig.2 co2 effect on vcmax25
#not added:   geom_jitter(alpha=0.2, aes(size = logr_se), position = position_jitter(w = 0.2, h = 0) ) +
co2_Jmax25 %>%
  ggplot( aes(x=exp_nam, y=logr)) +
  geom_jitter(alpha=0.2, aes(size = 1/logr_se)) +
  geom_crossbar( data = Jmax25_data_sitemean, aes(x=exp_nam, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5 ) +
  geom_boxplot(aes(x=exp_nam, y=pred_response_ratio,color=ecm_type),size=1)+ # pft_info,c3c4
  geom_hline( yintercept=0.0, size=0.5 )+
  labs(x="", y="CO2 effects on Jmax25", size=expression(paste("Error"^{-1}))) +
  coord_flip() 
ggsave(paste("~/data/output_gcme/colin/co2_2.jpg",sep=""))

#fig.3 measured or predicted response ratio of vc25 vs. jmax25 at plot-level
#here sensitivity coefficient doesn't account for uncertainty when 
vcmax25_pred <- aggregate(co2_vcmax25,by=list(co2_vcmax25$exp_nam,co2_vcmax25$pft_info,co2_vcmax25$c3c4,co2_vcmax25$ecm_type), FUN=mean, na.rm=TRUE)
vcmax25_pred_part <- vcmax25_pred[,c("Group.1","Group.2","Group.3","Group.4",
                                     "vcmax25_a_inst","vcmax25_e_inst",
                                     "mean_vcmax25_ambient","mean_vcmax25_elevated","pred_response_ratio",
                                     "co2_a","co2_e","Tg","PPFD","temp","vpd","n_plots")]
names(vcmax25_pred_part) <- c("exp_nam","pft_final","c3c4","ecm_type",
                              "vcmax25_a","vcmax25_e",
                              "pred_vcmax25_a","pred_vcmax25_e","logr_pred_vcmax25",
                              "co2_a","co2_e","Tg","PPFD","temp","vpd","n_plots")

Jmax25_pred <- aggregate(co2_Jmax25,by=list(co2_Jmax25$exp_nam,co2_Jmax25$pft_info,co2_Jmax25$c3c4,co2_Jmax25$ecm_type), FUN=mean, na.rm=TRUE)
Jmax25_pred_part <- Jmax25_pred[,c("Group.1","Group.2","Group.3","Group.4",
                                     "jmax25_a_inst","jmax25_e_inst",
                                     "mean_jmax25_ambient","mean_jmax25_elevated","pred_response_ratio",
                                     "co2_a","co2_e","Tg","PPFD","temp","vpd")]
names(Jmax25_pred_part) <- c("exp_nam","pft_final_j","c3c4_j","ecm_type_j",
                              "jmax25_a","jmax25_e",
                              "pred_jmax25_a","pred_jmax25_e","logr_pred_jmax25",
                              "co2_a_j","co2_e_j","Tg_j","PPFD_j","temp_j","vpd_j")

vc25_data_sitemean_obs <- vc25_data_sitemean[,c("exp_nam","middle","ymin","ymax")]
names(vc25_data_sitemean_obs) <- c("exp_nam","logr_obs_vcmax25","ymin_logr_obs_vcmax25","ymax_logr_obs_vcmax25")
Jmax25_data_sitemean_obs <- Jmax25_data_sitemean[,c("exp_nam","middle","ymin","ymax")]
names(Jmax25_data_sitemean_obs) <- c("exp_nam","logr_obs_jmax25","ymin_logr_obs_jmax25","ymax_logr_obs_jmax25")

photo_final <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp_nam"),all.x=TRUE),
                 list(vcmax25_pred_part,Jmax25_pred_part,vc25_data_sitemean_obs,Jmax25_data_sitemean_obs))

ggplot(photo_final) +
  geom_point(aes(x=logr_obs_jmax25, y=logr_obs_vcmax25,shape=pft_final,color=ecm_type,size=2))+
  geom_smooth(aes(x=logr_obs_jmax25, y=logr_obs_vcmax25),color="black",method="lm")+
  geom_smooth(aes(x=logr_pred_jmax25, y=logr_pred_vcmax25),color="red",method="lm")+
  labs(x = "response ratio jmax25",y="response ratio vcmax25")  +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))#+geom_text(aes(logr_obs_jmax25, y=logr_obs_vcmax25,label=exp_nam))
ggsave(paste("~/data/output_gcme/colin/co2_3.jpg",sep=""))

#fig.3 sensitivity cofficient - basing on logr_obs_vcmax25
#just check - one plot only has one co2_a (co2_e) set!
dim(co2_vcmax25 %>% group_by(exp_nam,co2_a,co2_e) %>% summarise(number=n()))
photo_final$coef_r_vcmax25 <- photo_final$logr_obs_vcmax25/(log(photo_final$co2_e/photo_final$co2_a))
photo_final$coef_r_jmax25 <- photo_final$logr_obs_jmax25/(log(photo_final$co2_e_j/photo_final$co2_a_j))

photo_final$coef_r_vcmax25_ymax <- photo_final$ymax_logr_obs_vcmax25/(log(photo_final$co2_e/photo_final$co2_a))
photo_final$coef_r_jmax25_ymax <- photo_final$ymax_logr_obs_jmax25/(log(photo_final$co2_e_j/photo_final$co2_a_j))
photo_final$coef_r_vcmax25_ymin <- photo_final$ymin_logr_obs_vcmax25/(log(photo_final$co2_e/photo_final$co2_a))
photo_final$coef_r_jmax25_ymin <- photo_final$ymin_logr_obs_jmax25/(log(photo_final$co2_e_j/photo_final$co2_a_j))

coef_r_vcmax25 <- photo_final$coef_r_vcmax25
coef_r_jmax25 <- photo_final$coef_r_jmax25

logppfd <- log(photo_final$PPFD)
tg <- photo_final$Tg
logvpd <- log(photo_final$vpd)

t1 <- lm(coef_r_vcmax25~logppfd+tg+logvpd)
t2 <- lm(coef_r_jmax25~logppfd+tg+logvpd)
summary(t1)
summary(t2)

photo_final$coef_r_pred_vcmax25 <- log(photo_final$pred_vcmax25_e/photo_final$pred_vcmax25_a)/(log(photo_final$co2_e/photo_final$co2_a))
photo_final$coef_r_pred_jmax25 <- log(photo_final$pred_jmax25_e/photo_final$pred_jmax25_a)/(log(photo_final$co2_e_j/photo_final$co2_a_j))

#direct comparison for response ratio and sensitivity coefficient of vcmax25 
ggplot(photo_final) +
  geom_point(aes(x=logr_obs_vcmax25, y=logr_pred_vcmax25,shape=pft_final,color=ecm_type,size=2))+
  geom_smooth(aes(x=logr_obs_vcmax25, y=logr_pred_vcmax25),color="black",method="lm")+
  labs(x = "logr vcmax25[obs]",y="logr vcmax25[pred]")  +theme_classic()+xlim(-0.2,0.2)+ylim(-0.2,0.2)+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))#+geom_text(aes(logr_obs_jmax25, y=logr_obs_vcmax25,label=exp_nam))

ggplot(photo_final) +
  geom_point(aes(x=coef_r_vcmax25, y=coef_r_pred_vcmax25,shape=pft_final,color=ecm_type,size=2))+
  geom_smooth(aes(x=coef_r_vcmax25, y=coef_r_pred_vcmax25),color="black",method="lm")+
  labs(x = "coef_r vcmax25[obs]",y="coef_r vcmax25[pred]")  +theme_classic()+xlim(-0.2,0.2)+ylim(-0.2,0.2)+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))#+geom_text(aes(logr_obs_jmax25, y=logr_obs_vcmax25,label=exp_nam))

#compare coef_r_vcmax25 in observation and prediction
# first - compute stanard deviation and then variance
photo_final$logr_var_vcmax25 <- ((photo_final$ymax_logr_obs_vcmax25-photo_final$logr_obs_vcmax25)/1.96)^2
photo_final$logr_var_jmax25 <- ((photo_final$ymax_logr_obs_jmax25-photo_final$logr_obs_jmax25)/1.96)^2
#create a function to divide them into specific cateory
#input:df, data_type, logr and logr_var
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

meta_box_vc25_ecm <- agg_meta_plots(photo_final,"ecm_type","logr_obs_vcmax25","logr_var_vcmax25")
meta_box_j25_ecm <- agg_meta_plots(photo_final,"ecm_type","logr_obs_jmax25","logr_var_jmax25")

photo_final %>%
  ggplot( aes(x=ecm_type, y=logr_obs_vcmax25)) +
  geom_jitter(alpha=0.2, aes(size = 1/(sqrt(logr_var_vcmax25)/sqrt(n_plots)))) + # this is inverse to standard error
  geom_crossbar( data = meta_box_vc25_ecm, aes(x=type_name, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5 ) +
  geom_boxplot(aes(x=ecm_type, y=logr_pred_vcmax25),color="red",size=1,outlier.shape = NA)+ # pft_info,c3c4
  geom_hline( yintercept=0.0, size=0.5 )+
  labs(x="", y="response ratio of vcmax25, ECM types", size=expression(paste("Error"^{-1}))) +
  coord_flip() 
ggsave(paste("~/data/output_gcme/colin/co2_4a.jpg",sep=""))

photo_final %>%
  ggplot(aes(x=ecm_type, y=logr_obs_jmax25)) +
  geom_jitter(alpha=0.2, aes(size = 1/(sqrt(logr_var_jmax25)/sqrt(n_plots)))) + # this is inverse to standard error
  geom_crossbar( data = meta_box_j25_ecm, aes(x=type_name, y=middle, ymin=ymin, ymax=ymax), alpha = 0.6, width = 0.5 ) +
  geom_boxplot(aes(x=ecm_type, y=logr_pred_jmax25),color="red",size=1,outlier.shape = NA)+ # pft_info,c3c4
  geom_hline( yintercept=0.0, size=0.5 )+
  labs(x="", y="response ratio of jmax25, ECM types", size=expression(paste("Error"^{-1}))) +
  coord_flip() 
ggsave(paste("~/data/output_gcme/colin/co2_4b.jpg",sep=""))

#co2 on jmax25/vcmax25
photo_final$jv_obs <- log((photo_final$jmax25_e/photo_final$vcmax25_e)/(photo_final$jmax25_a/photo_final$vcmax25_a))
photo_final$jv_pred <- log((photo_final$pred_jmax25_e/photo_final$pred_vcmax25_e)/(photo_final$pred_jmax25_a/photo_final$pred_vcmax25_a))

ggplot(photo_final) +
  geom_point(aes(x=jv_obs, y=jv_pred,shape=pft_final,color=ecm_type,size=2))+
  geom_smooth(aes(x=jv_obs, y=jv_pred),color="black",method="lm")+
  geom_vline(xintercept = 0)+
  labs(x = "logr_J/V[obs]",y="logr_J/V[pred]")  +theme_classic()+
  theme(axis.text=element_text(size=20),axis.title =element_text(size=20))#+geom_text(aes(logr_obs_jmax25, y=logr_obs_vcmax25,label=exp_nam))
ggsave(paste("~/data/output_gcme/colin/co2_5.jpg",sep=""))

photo_final %>%
  ggplot( aes(x=ecm_type, y=jv_obs)) + #ecm_type, c3c4, pft_info
  geom_boxplot(outlier.shape = NA) +
  geom_boxplot(aes(x=ecm_type, y=jv_pred),color="red",outlier.shape = NA) +
  labs(x="", y="Sensitivity coefficient of jmax25/vcmax25 in response to co2") +
  geom_hline( yintercept=0.0, size=0.5)+
  coord_flip() 
ggsave(paste("~/data/output_gcme/colin/co2_4c.jpg",sep=""))

#model bias
hist(co2_vcmax25$vcmax25_a)

#second part - warming
temp_vcmax25 <- subset(vc25_data,treatment=="w")
temp_jmax25 <- subset(Jmax25_data,treatment=="w")

tempco2_vcmax25 <- subset(vc25_data,treatment=="cw")
tempco2_jmax25 <- subset(Jmax25_data,treatment=="cw")

tempppco2_vcmax25 <- subset(vc25_data,treatment=="cw3")
tempppco2_jmax25 <- subset(Jmax25_data,treatment=="cw")

temp_vcmax25 %>%
  ggplot( aes(x=exp_nam, y=logr)) + #ecm_type, c3c4, pft_info
  geom_boxplot(outlier.shape = NA) +
  geom_jitter( color = rgb(0,0,0,0.3), aes( size = 1/logr_se ), position = position_jitter(w = 0.2, h = 0) ) +
  geom_boxplot(aes(x=exp_nam, y=pred_response_ratio),color="red",outlier.shape = NA)+
  geom_hline( yintercept=0.0, size=0.5 )+scale_fill_manual(values=c("yellow","purple","brown","green"))+
  labs(x="", y="1-degree warming on vcmax25", size=expression(paste("Standard Error"))) +
  coord_flip() 
ggsave(paste("~/data/output_gcme/colin/temp_1.jpg",sep=""))

temp_jmax25 %>%
  ggplot( aes(x=exp_nam, y=logr)) + #ecm_type, c3c4, pft_info
  geom_boxplot(outlier.shape = NA) +
  geom_jitter( color = rgb(0,0,0,0.3), aes( size = 1/logr_se ), position = position_jitter(w = 0.2, h = 0) ) +
  geom_boxplot(aes(x=exp_nam, y=pred_response_ratio),color="red",outlier.shape = NA)+
  geom_hline( yintercept=0.0, size=0.5 )+scale_fill_manual(values=c("yellow","purple","brown","green"))+
  labs(x="", y="1-degree warming on jmax25", size=expression(paste("Standard Error"))) +
  coord_flip() 
ggsave(paste("~/data/output_gcme/colin/temp_2.jpg",sep=""))

tempco2_vcmax25 %>%
  ggplot( aes(x=exp_nam, y=logr)) + #ecm_type, c3c4, pft_info
  geom_boxplot(outlier.shape = NA) +
  geom_jitter( color = rgb(0,0,0,0.3), aes( size = 1/logr_se ), position = position_jitter(w = 0.2, h = 0) ) +
  geom_boxplot(aes(x=exp_nam, y=pred_response_ratio),color="red",outlier.shape = NA)+
  geom_hline( yintercept=0.0, size=0.5 )+scale_fill_manual(values=c("yellow","purple","brown","green"))+
  labs(x="", y="1-degree warming + co2 on vcmax25", size=expression(paste("Standard Error"))) +
  coord_flip() 
ggsave(paste("~/data/output_gcme/colin/temp_1a.jpg",sep=""))

tempco2_jmax25 %>%
  ggplot( aes(x=exp_nam, y=logr)) + #ecm_type, c3c4, pft_info
  geom_boxplot(outlier.shape = NA) +
  geom_jitter( color = rgb(0,0,0,0.3), aes( size = 1/logr_se ), position = position_jitter(w = 0.2, h = 0) ) +
  geom_boxplot(aes(x=exp_nam, y=pred_response_ratio),color="red",outlier.shape = NA)+
  geom_hline( yintercept=0.0, size=0.5 )+scale_fill_manual(values=c("yellow","purple","brown","green"))+
  labs(x="", y="1-degree warming +co2 on jmax25", size=expression(paste("Standard Error"))) +
  coord_flip() 
ggsave(paste("~/data/output_gcme/colin/temp_2a.jpg",sep=""))

# n fertilization..will do later
vcmax25_fertilization <- subset(vcmax25_final2,treatment=="f"|treatment=="cf")
Jmax25_fertilization <- subset(Jmax25_final2,treatment=="f"|treatment=="cf")

#final merge vcmax25 and jmax25
vcmax25_fertilization_ecm <- merge(vcmax25_fertilization,df_experiments_coord_co2,by=c("exp_nam"),all.x=TRUE)
vcmax25_fertilization_ecm%>% group_by(AM)  %>% summarise(number = n())
vcmax25_fertilization_ecm%>% group_by(EcM)  %>% summarise(number = n())
vcmax25_fertilization_ecm$ecm_type <- "NA"
vcmax25_fertilization_ecm$ecm_type[vcmax25_fertilization_ecm$EcM==TRUE] <- "ECM"
vcmax25_fertilization_ecm$ecm_type[vcmax25_fertilization_ecm$AM==TRUE] <- "AM"
vcmax25_fertilization_ecm$ecm_type[vcmax25_fertilization_ecm$ErM==TRUE] <- "ErM"
vcmax25_fertilization_ecm$ecm_type[vcmax25_fertilization_ecm$Moss==TRUE] <- "Moss"
#also defines pft and c3/c4
vcmax25_fertilization_ecm$pft_info <- "NA"
vcmax25_fertilization_ecm$pft_info[vcmax25_fertilization_ecm$Fern==TRUE] <-"fern"
vcmax25_fertilization_ecm$pft_info[vcmax25_fertilization_ecm$C3.grass==TRUE] <-"grassland"
vcmax25_fertilization_ecm$pft_info[vcmax25_fertilization_ecm$C4.grass==TRUE] <-"grassland"
vcmax25_fertilization_ecm$pft_info[vcmax25_fertilization_ecm$C3.herb==TRUE] <-"grassland"
vcmax25_fertilization_ecm$pft_info[vcmax25_fertilization_ecm$C4.herb==TRUE] <-"grassland"
vcmax25_fertilization_ecm$pft_info[vcmax25_fertilization_ecm$pft_info=="NA"] <- "forest"
vcmax25_fertilization_ecm$c3c4 <- "NA"
vcmax25_fertilization_ecm$c3c4[vcmax25_fertilization_ecm$C4.grass==TRUE] <-"c4"
vcmax25_fertilization_ecm$c3c4[vcmax25_fertilization_ecm$C4.herb==TRUE] <-"c4"
vcmax25_fertilization_ecm$c3c4[vcmax25_fertilization_ecm$c3c4=="NA"] <- "c3"

jmax25_fertilization_ecm <- merge(Jmax25_fertilization,df_experiments_coord_co2,by=c("exp_nam"),all.x=TRUE)
jmax25_fertilization_ecm%>% group_by(AM)  %>% summarise(number = n())
jmax25_fertilization_ecm%>% group_by(EcM)  %>% summarise(number = n())
jmax25_fertilization_ecm$ecm_type <- "NA"
jmax25_fertilization_ecm$ecm_type[jmax25_fertilization_ecm$EcM==TRUE] <- "ECM"
jmax25_fertilization_ecm$ecm_type[jmax25_fertilization_ecm$AM==TRUE] <- "AM"
jmax25_fertilization_ecm$ecm_type[jmax25_fertilization_ecm$ErM==TRUE] <- "ErM"
jmax25_fertilization_ecm$ecm_type[jmax25_fertilization_ecm$Moss==TRUE] <- "Moss"
#also defines pft and c3/c4
jmax25_fertilization_ecm$pft_info <- "NA"
jmax25_fertilization_ecm$pft_info[jmax25_fertilization_ecm$Fern==TRUE] <-"fern"
jmax25_fertilization_ecm$pft_info[jmax25_fertilization_ecm$C3.grass==TRUE] <-"grassland"
jmax25_fertilization_ecm$pft_info[jmax25_fertilization_ecm$C4.grass==TRUE] <-"grassland"
jmax25_fertilization_ecm$pft_info[jmax25_fertilization_ecm$C3.herb==TRUE] <-"grassland"
jmax25_fertilization_ecm$pft_info[jmax25_fertilization_ecm$C4.herb==TRUE] <-"grassland"
jmax25_fertilization_ecm$pft_info[jmax25_fertilization_ecm$pft_info=="NA"] <- "forest"
jmax25_fertilization_ecm$c3c4 <- "NA"
jmax25_fertilization_ecm$c3c4[jmax25_fertilization_ecm$C4.grass==TRUE] <-"c4"
jmax25_fertilization_ecm$c3c4[jmax25_fertilization_ecm$C4.herb==TRUE] <-"c4"
jmax25_fertilization_ecm$c3c4[jmax25_fertilization_ecm$c3c4=="NA"] <- "c3"

vcmax25_fertilization_ecm$factors #remove cfi

jmax25_fertilization_ecm$factors

vc25_f_data <- response_ratio(vcmax25_fertilization_ecm,"vcmax25")
vc25_f_data_sitemean <- agg_plot(vc25_f_data,"vcmax25")

j25_f_data <- response_ratio(jmax25_fertilization_ecm,"jmax25")
j25_f_data_sitemean <- agg_plot(vc25_data,"vcmax25")

#needs to divide f and cf of treatment --> then f further into d/f magnitudes
subset(vc25_f_data,treatment=="f") %>%
  ggplot( aes(x=exp_nam, y=logr,fill=ecm_type)) + #ecm_type, c3c4, pft_info
  geom_jitter(alpha=0.2, aes(size = logr_se), position = position_jitter(w = 0.2, h = 0) ) +
  geom_boxplot(outlier.shape = NA,alpha=0.3) +
  geom_hline( yintercept=0.0, size=0.5 )+scale_fill_manual(values=c("yellow","purple","brown","green"))+
  labs(x="", y="N fertilization effects on vcmax25", size=expression(paste("Standard Error"))) +
  coord_flip() 
ggsave(paste("~/data/output_gcme/colin/Nfer_1.jpg",sep=""))

ggplot(data=subset(vc25_f_data,treatment=="f"),aes(x=exp_nam, y=logr)) + 
  geom_boxplot(outlier.shape = NA,alpha=0.5)+
  geom_boxplot(data=subset(vc25_f_data,treatment=="cf"),color="red",outlier.shape = NA,alpha=0.5)+
  geom_hline( yintercept=0.0, size=0.5 )+scale_fill_manual(values=c("yellow","purple","brown","green"))+
  labs(x="", y="N (+eco2) fertilization effects on vcmax25") +
  coord_flip() 
ggsave(paste("~/data/output_gcme/colin/Nfer_1a.jpg",sep=""))

vc25_f_data$exp_nam_years <- paste(vc25_f_data$exp_nam,sep = "_",vc25_f_data$Year)
subset(vc25_f_data,treatment=="f") %>%
  ggplot( aes(x=exp_nam_years, y=logr,fill=ecm_type)) + #ecm_type, c3c4, pft_info
  geom_jitter(alpha=0.2, aes(size = logr_se), position = position_jitter(w = 0.2, h = 0) ) +
  geom_boxplot(outlier.shape = NA,alpha=0.3) +
  geom_hline( yintercept=0.0, size=0.5 )+scale_fill_manual(values=c("yellow","purple","brown","green"))+
  labs(x="", y="N fertilization effects on vcmax25", size=expression(paste("Standard Error"))) +
  coord_flip() 
ggsave(paste("~/data/output_gcme/colin/Nfer_1b.jpg",sep=""))


subset(j25_f_data,treatment=="f") %>%
  ggplot( aes(x=exp_nam, y=logr,fill=ecm_type)) + #ecm_type, c3c4, pft_info
  geom_jitter(alpha=0.2, aes(size = logr_se), position = position_jitter(w = 0.2, h = 0) ) +
  geom_boxplot(outlier.shape = NA,alpha=0.3) +
  geom_hline( yintercept=0.0, size=0.5 )+scale_fill_manual(values=c("yellow","purple","brown","green"))+
  labs(x="", y="N fertilization effects on jmax25", size=expression(paste("Standard Error"))) +
  coord_flip() 
ggsave(paste("~/data/output_gcme/colin/Nfer_2.jpg",sep=""))

ggplot(data=subset(j25_f_data,treatment=="f"),aes(x=exp_nam, y=logr)) + 
  geom_boxplot(outlier.shape = NA,alpha=0.5)+
  geom_boxplot(data=subset(j25_f_data,treatment=="cf"),color="red",outlier.shape = NA,alpha=0.5)+
  geom_hline( yintercept=0.0, size=0.5 )+scale_fill_manual(values=c("yellow","purple","brown","green"))+
  labs(x="", y="N (+eco2) fertilization effects on jmax25") +
  coord_flip() 
ggsave(paste("~/data/output_gcme/colin/Nfer_2a.jpg",sep=""))

vc25_f_data$nfertQ_e
vc25_f_data%>% group_by(exp_nam,nfertQ_e)  %>% summarise(number = n())# all reasonable unit


#merged with lma, narea, lai, anpp, check site specific info (low soil availbility)? check original reference?
LMA_SLA <- subset(df,Data_type=="leaf_mass_per_area"|Data_type=="leaf_mass_per_unit_leaf_area"|Data_type=="LMA"|
                    Data_type=="specific_leaf_area"|Data_type=="SLA")
LMA_SLA%>% group_by(Unit)  %>% summarise(number = n())# all reasonable unit
LMA_SLA$LMA_response_ratio <- log(LMA_SLA$elevated/LMA_SLA$ambient)
#because they are inverse!
LMA_SLA$LMA_response_ratio[LMA_SLA$Data_type=="specific_leaf_area"] <- log(LMA_SLA$ambient[LMA_SLA$Data_type=="specific_leaf_area"]/LMA_SLA$elevated[LMA_SLA$Data_type=="specific_leaf_area"])
LMA_SLA$LMA_response_ratio[LMA_SLA$Data_type=="SLA"] <- log(LMA_SLA$ambient[LMA_SLA$Data_type=="SLA"]/LMA_SLA$elevated[LMA_SLA$Data_type=="SLA"])

LMA_SLA_final <-LMA_SLA[,c("exp_nam","treatment","LMA_response_ratio")]
LMA_SLA_final <- aggregate(LMA_SLA_final,by=list(LMA_SLA_final$exp_nam,LMA_SLA_final$treatment),
                           FUN=mean, na.rm=TRUE)
LMA_SLA_final <-LMA_SLA_final[,c("Group.1","Group.2","LMA_response_ratio")]
names(LMA_SLA_final) <- c("exp_nam","treatment","LMA_response_ratio")

#how about leaf N
leaf_N <- subset(df,Data_type=="leaf_N")
leaf_N$leafN_response_ratio <- log(leaf_N$elevated/leaf_N$ambient)

leafNmass <-  subset(leaf_N,Unit=="%"|Unit=="g/g"|Unit=="g/kg"|Unit=="mg_/_g"|Unit=="mg/g"|Unit=="g/mg"|
                       Unit=="mg/g_"|Unit=="mg/kg"|Unit=="mol/g"|Unit=="mmol/g")
leafNarea <-  subset(leaf_N,Unit!="%" & Unit!="g/g" & Unit!="g/kg" & Unit!="mg_/_g" & Unit!="mg/g" & Unit!="g/mg" &
                       Unit!="mg/g_" & Unit!="mg/kg" & Unit!="mol/g" & Unit!="mmol/g")

leafNmass_final <- aggregate(leafNmass,by=list(leafNmass$exp_nam,leafNmass$treatment),
                           FUN=mean, na.rm=TRUE)
leafNmass_final <-leafNmass_final[,c("Group.1","Group.2","leafN_response_ratio")]
names(leafNmass_final) <- c("exp_nam","treatment","Nmass_response_ratio")

leafNarea_final <- aggregate(leafNarea,by=list(leafNarea$exp_nam,leafNarea$treatment),
                             FUN=mean, na.rm=TRUE)
leafNarea_final <-leafNarea_final[,c("Group.1","Group.2","leafN_response_ratio")]
names(leafNarea_final) <- c("exp_nam","treatment","Narea_response_ratio")

vc25_f_data_traits <-Reduce(function(x,y) merge(x = x, y = y, by = c("exp_nam","treatment"),all.x=TRUE),
                 list(vc25_f_data,LMA_SLA_final,leafNarea_final,leafNmass_final))

subset(vc25_f_data_traits,treatment=="f") %>%
  ggplot( aes(x=exp_nam, y=logr,fill=ecm_type)) + #ecm_type, c3c4, pft_info
  geom_jitter(alpha=0.2, aes(size = logr_se), position = position_jitter(w = 0.2, h = 0) ) +
  geom_boxplot(outlier.shape = NA,alpha=0.3) +
  geom_boxplot(aes(x=exp_nam, y=LMA_response_ratio),color="green",outlier.shape = NA,alpha=1)+
  geom_boxplot(aes(x=exp_nam, y=Nmass_response_ratio),color="blue",outlier.shape = NA,alpha=1)+
  geom_hline( yintercept=0.0, size=0.5 )+scale_fill_manual(values=c("yellow","purple","brown","green"))+
  labs(x="", y="N fertilization effects on vcmax25 - green:LMA, blue:Nmass", size=expression(paste("Standard Error"))) +
  coord_flip() 
ggsave(paste("~/data/output_gcme/colin/Nfer_3.jpg",sep=""))

#now, start soil N collection
summary(vcmax25_final2)
summary(Jmax25_final2)
all_expnam <- aggregate(vcmax25_final2,by=list(vcmax25_final2$exp_nam), FUN=mean, na.rm=TRUE)$Group.1
df_part <- df %>% filter(exp_nam %in% all_expnam)
df_part_soil <- df_part[grep("soil", df_part$Data_type),] # all soil data - just see
#check how many variables
df_part_soil%>% group_by(Data_type,exp_nam) %>% summarise(number=n()) %>%
  group_by(Data_type) %>% summarise(number=n()) %>% arrange(desc(number), .by_group = TRUE)
#