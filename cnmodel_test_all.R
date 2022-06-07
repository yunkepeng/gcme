library(rsofun)
library(dplyr)
pars <- list(
  kphio                 = 0.04607080,
  soilm_par_a           = 2.75687824,
  soilm_par_b           = 1.68140444,
  tau_acclim_tempstress = 7.35259044,
  par_shape_tempstress  = 0.09863961,
  f_nretain             = 0.500000,
  fpc_tree_max          = 0.950000,
  growtheff             = 0.600000,
  r_root                = 0.913000,
  r_sapw                = 0.044000,
  exurate               = 0.003000,
  k_decay_leaf_base     = 1.00000,
  k_decay_leaf_width    = 2.00000,
  k_decay_root          = 1.00000,
  k_decay_labl          = 0.00000,
  k_decay_sapw          = 1.00000,
  r_cton_root           = 37.0000,
  r_cton_wood           = 100.000,
  ncw_min               = 0.056,
  r_n_cw_v              = 0.4,
  r_ctostructn_leaf     = 35.0000,
  kbeer                 = 0.500000,
  gddbase               = 5.0,
  ramp                  = 0.0,
  phentype              = 2.0,
  perc_k1               = 5.0,        
  thdiff_wp             = 0.2,          
  thdiff_whc15          = 0.8,
  thdiff_fc             = 0.4,          
  forg                  = 0.01,
  wbwp                  = 0.029,  
  por                   = 0.421,    
  fsand                 = 0.82,      
  fclay                 = 0.06,      
  fsilt                 = 0.12,  
  kA                    = 107,     
  kalb_sw               = 0.17,    
  kalb_vis              = 0.03,    
  kb                    = 0.20,    
  kc                    = 0.25,    
  kCw                   = 1.05,    
  kd                    = 0.50,    
  ke                    = 0.0167,  
  keps                  = 23.44,   
  kWm                   = 220.0,   
  kw                    = 0.26,    
  komega                = 283.0,
  maxmeltrate           = 3.0,
  klitt_af10            = 1.2,
  klitt_as10            = 0.35,
  klitt_bg10            = 0.35,
  kexu10                = 50.0,
  ksoil_fs10            = 0.021,
  ksoil_sl10            = 7.0e-04,
  ntoc_crit1            = 0.45,
  ntoc_crit2            = 0.76,
  cton_microb           = 10.0,
  cton_soil             = 9.77,
  fastfrac              = 0.985,
  eff_nup               = 0.600000,
  minimumcostfix        = 1.000000,
  fixoptimum            = 25.15000,
  a_param_fix           = -3.62000,
  b_param_fix           = 0.270000,
  maxnitr               = 0.1,
  non                   = 0.01,
  n2on                  = 0.0005,
  kn                    = 83.0,
  kdoc                  = 17.0,
  docmax                = 1.0,
  dnitr2n2o             = 0.01,
  beta                  = 146.000000,
  rd_to_vcmax           = 0.01400000,
  tau_acclim            = 10,
  tmppar                = 9999)

plot_level<- read.csv("~/data/gcme/kevin/forcing/forcing_info_allnew.csv")
forcing_path <- "/Users/yunpeng/data/gcme/kevin/forcing/climate_large/"

#problem sites
wrong_sites <- c(9,11,12,14,18,19,22,25,28,29,32,33,34,37,39,40,41,42,45,48,49,50,62,65,66,67,
  70,77,78,79,80,81,82,83,85,86,89,92,93)
seq2 <- c(1:93)
correct_sites <- seq2[!seq2 %in% wrong_sites]

#
plot_level$vcmax25 <- NA
plot_level$lai <- NA
plot_level$npp <- NA
plot_level$narea <- NA
plot_level$lma <- NA

#this plots not work
for (i in wrong_sites){
  print(i)
  siteinfo_site <- plot_level[i,c("sitename","lon","lat","elv","year_start","year_end")]
  sitename <- siteinfo_site$sitename
  
  df1 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for ambient
  df1$date <- as.Date(df1$date);df1$doy <- 1:length(df1$date)
  forcing <- rsofun::p_model_drivers
  library(dplyr)
  forcing$forcing[[1]] <- as_tibble(df1[,c("date","temp","prec","vpd","ppfd","patm","ccov_int","ccov","snow","rain","fapar","co2","doy","tmin","tmax")])
  forcing$sitename <- df1$sitename[1]
  forcing$site_info[[1]]$lon <- df1$lon[1]
  forcing$site_info[[1]]$lat <- df1$lat[1]
  forcing$site_info[[1]]$elv <-  df1$elv[1]
  forcing$site_info[[1]]$date_start <- df1$date[1]
  forcing$site_info[[1]]$date_end <- df1$date[length(df1$date)]
  forcing$params_siml[[1]]$firstyeartrend <- 1991
  forcing$params_siml[[1]]$nyeartrend <- 20
  tmp <- forcing %>% mutate(forcing = purrr::map(forcing, ~mutate(., fharv = 0.0, dno3 = 0.1,dnh4 = 0.1)))
  tmp$params_siml[[1]]$spinupyears <- 1500
  tmp$params_siml[[1]]$recycle <- 1
  modlist1 <- runread_pmodel_f(tmp,par = pars)#ambient
  modlist1 <- as.data.frame(modlist1$data[[1]])
  plot_level$vcmax25[i] <- max(modlist1$vcmax25,na.rm=TRUE)*1000000
  plot_level$lai[i] <- max(modlist1$lai,na.rm=TRUE)
  plot_level$npp[i] <- max(modlist1$npp,na.rm=TRUE)
  plot_level$narea[i] <- max(modlist1$narea,na.rm=TRUE)
  plot_level$lma[i] <- max(modlist1$lma,na.rm=TRUE)
}

#this plots work well
for (i in correct_sites){
  print(i)
  siteinfo_site <- plot_level[i,c("sitename","lon","lat","elv","year_start","year_end")]
  sitename <- siteinfo_site$sitename
  
  df1 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for ambient
  df1$date <- as.Date(df1$date);df1$doy <- 1:length(df1$date)
  forcing <- rsofun::p_model_drivers
  library(dplyr)
  forcing$forcing[[1]] <- as_tibble(df1[,c("date","temp","prec","vpd","ppfd","patm","ccov_int","ccov","snow","rain","fapar","co2","doy","tmin","tmax")])
  forcing$sitename <- df1$sitename[1]
  forcing$site_info[[1]]$lon <- df1$lon[1]
  forcing$site_info[[1]]$lat <- df1$lat[1]
  forcing$site_info[[1]]$elv <-  df1$elv[1]
  forcing$site_info[[1]]$date_start <- df1$date[1]
  forcing$site_info[[1]]$date_end <- df1$date[length(df1$date)]
  forcing$params_siml[[1]]$firstyeartrend <- 1991
  forcing$params_siml[[1]]$nyeartrend <- 20
  tmp <- forcing %>% mutate(forcing = purrr::map(forcing, ~mutate(., fharv = 0.0, dno3 = 0.1,dnh4 = 0.1)))
  tmp$params_siml[[1]]$spinupyears <- 1500
  tmp$params_siml[[1]]$recycle <- 1
  modlist1 <- runread_pmodel_f(tmp,par = pars)#ambient
  modlist1 <- as.data.frame(modlist1$data[[1]])
  plot_level$vcmax25[i] <- max(modlist1$vcmax25,na.rm=TRUE)*1000000
  plot_level$lai[i] <- max(modlist1$lai,na.rm=TRUE)
  plot_level$npp[i] <- max(modlist1$npp,na.rm=TRUE)
  plot_level$narea[i] <- max(modlist1$narea,na.rm=TRUE)
  plot_level$lma[i] <- max(modlist1$lma,na.rm=TRUE)
}


#distribution of plots
library(maps)
library(rworldmap)

newmap <- getMap(resolution = "low")
sp::plot(newmap, xlim = c(-180, 180), ylim = c(-75, 75), asp = 1)
points(plot_level[correct_sites,c("lon")],plot_level[correct_sites,c("lat")], col="green", pch=16,cex=2)
points(plot_level[wrong_sites,c("lon")],plot_level[wrong_sites,c("lat")], col="red", pch=16,cex=2)
