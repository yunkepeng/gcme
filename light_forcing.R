walker <- read.csv("/Users/yunpeng/data/leaf_traits/Walker/LEAF_PHOTOSYNTHESIS_TRAITS_1224/data/Leaf_Photosynthesis_Traits.csv")

#vcmax
light_vcmax <- subset(walker,Treatment_light!="not applicable") %>% 
  group_by(lon,lat,z,Year,species,Treatment_light)  %>% summarise(mean = mean(Vcmax))

#prepare forcing
aaa <- as.data.frame(unique(light_vcmax[,c("lon","lat","z")]))

aaa$year_start <- 1991
aaa$year_end <- 2010

aaa$elv <- aaa$z
aaa$sitename <- paste("light",1:nrow(aaa),sep="")

siteinfo_final <- aaa[,c("sitename","lon","lat","elv","year_start","year_end")]

devtools::load_all("~/yunkepeng/gcme/pmodel/ingestr/")
devtools::load_all("~/yunkepeng/gcme/pmodel/rsofun/")

for (i in 1:(nrow(siteinfo_final))){
  df_watch <- ingest(
    siteinfo  = siteinfo_final[i,],
    source    = "watch_wfdei",
    getvars   = c("temp", "prec", "ppfd", "vpd", "patm"), 
    dir       = "/Volumes/My Passport/data/watch_wfdei/",
    settings  = list(correct_bias = "worldclim", dir_bias = "/Volumes/My Passport/data/worldclim/"))
  
  df_cru <- ingest(
    siteinfo  = siteinfo_final[i,],
    source    = "cru",
    getvars   = "ccov",
    dir       = "/Volumes/My Passport/data/cru/ts_4.01/")
  
  df_co2 <- ingest(
    siteinfo  = siteinfo_final[i,],
    source  = "co2_mlo",
    verbose = FALSE)
  
  df_co2_final <- as.data.frame(df_co2$data)
  
  df_co2_final2 <- df_co2_final[!(format(df_co2_final$date,"%m") == "02" & format(df_co2_final$date, "%d") == "29"), , drop = FALSE] # make columns consistent
  
  co2 <- df_co2_final2$co2
  
  ddf_meteo <- as_tibble(cbind(as.data.frame(df_watch$data),as.data.frame(df_cru$data)[,c("ccov_int","ccov")],co2))
  
  ddf_meteo$fapar <- 1
  ddf_meteo$tmax <- ddf_meteo$temp
  ddf_meteo$tmin <- ddf_meteo$temp
  
  ddf_meteo$sitename <- siteinfo_final[i,c("sitename")]
  ddf_meteo$lon <- siteinfo_final[i,c("lon")]
  ddf_meteo$lat <-  siteinfo_final[i,c("lat")]
  ddf_meteo$elv <-  siteinfo_final[i,c("elv")]
  ddf_meteo$year_start <-  siteinfo_final[i,c("year_start")]
  ddf_meteo$year_end <-  siteinfo_final[i,c("year_end")]
  
  csvfile <- paste("~/data/leaf_traits/Walker/forcing/climate/",siteinfo_final$sitename[i],".csv",sep = "")
  write.csv(ddf_meteo, csvfile, row.names = TRUE)
  print(i)    
}

#work on rsofun
siteinfo_final

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

forcing_path <- "~/data/leaf_traits/Walker/forcing/climate/"

for (i in 1:nrow(siteinfo_final)){
  siteinfo_final_one <- siteinfo_final[i,c("sitename","lon","lat","elv","year_start","year_end")]
  sitename <- siteinfo_final_one$sitename
  
  df1 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for ambient
  df2 <- read.csv(paste(forcing_path,sitename,".csv",sep="")) # this is for elevated
  
  #Q1 - add warming directly
  df2$ppfd <- df1$ppfd*1.5
  
  #Q2 - if having co2 then replacing it directly - if not then using co2 at measurement year
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
      firstyeartrend = siteinfo_final_one$year_start,
      nyeartrend = siteinfo_final_one$year_end-siteinfo_final_one$year_start+1), 
    siteinfo_final_one, 
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
      firstyeartrend = siteinfo_final_one$year_start,
      nyeartrend = siteinfo_final_one$year_end-siteinfo_final_one$year_start+1), 
    siteinfo_final_one, 
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
  
  siteinfo_final$mean_vcmax25_ambient[i] <- mean_vcmax25_ambient
  siteinfo_final$mean_vcmax25_elevated[i] <- mean_vcmax25_elevated
  siteinfo_final$mean_vcmax_ambient[i] <- mean_vcmax_ambient
  siteinfo_final$mean_vcmax_elevated[i] <- mean_vcmax_elevated
  siteinfo_final$pred_vcmax_logr[i] <- log(siteinfo_final$mean_vcmax_elevated[i]/siteinfo_final$mean_vcmax_ambient[i])
  siteinfo_final$pred_vcmax25_logr[i] <- log(siteinfo_final$mean_vcmax25_elevated[i]/siteinfo_final$mean_vcmax25_ambient[i])

  siteinfo_final$mean_jmax25_ambient[i] <- mean_jmax25_ambient
  siteinfo_final$mean_jmax25_elevated[i] <- mean_jmax25_elevated
  siteinfo_final$mean_jmax_ambient[i] <- mean_jmax_ambient
  siteinfo_final$mean_jmax_elevated[i] <- mean_jmax_elevated
  siteinfo_final$pred_jmax_logr[i] <- log(siteinfo_final$mean_jmax_elevated[i]/siteinfo_final$mean_jmax_ambient[i])
  siteinfo_final$pred_jmax25_logr[i] <- log(siteinfo_final$mean_jmax25_elevated[i]/siteinfo_final$mean_jmax25_ambient[i])

  siteinfo_final$PPFD[i] <- mean(df1$ppfd,na.rm=TRUE)*1000000 # in umol/m2/s
  siteinfo_final$vpd[i] <- mean(df1$vpd,na.rm=TRUE)/1000 #in kPa
  siteinfo_final$temp[i] <- mean(df1$temp,na.rm=TRUE) #in degree celcius
}
summary(siteinfo_final)

light_walker <- siteinfo_final

walker <- read.csv("/Users/yunpeng/data/leaf_traits/Walker/LEAF_PHOTOSYNTHESIS_TRAITS_1224/data/Leaf_Photosynthesis_Traits.csv")

#vcmax
light_vcmax <- subset(walker,Treatment_light!="not applicable") %>% 
  group_by(lon,lat,species,Treatment_light)  %>% summarise(mean = mean(Vcmax))

subset(walker,Treatment_light!="not applicable") %>% 
  group_by(lon,lat,species,Treatment_light,Irradiance)  %>% summarise(mean = mean(Vcmax))

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

light_vcmax_points <- merge(light_vcmax_points,light_walker, by=c("lon","lat"),all.x=TRUE)

l1 <- light_vcmax_points %>% ggplot( aes(x=type_name, y=middle)) + geom_boxplot()+
  geom_boxplot(aes(x=type_name, y=pred_vcmax25_logr),color="red")+
  geom_point(alpha = 0.6, width = 0.5) +geom_hline( yintercept=0.0, size=0.5)+ ylim(-0.5,1)+
  labs(y="light effect on vcmax",x=" ") + theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))



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

light_jmax_points <- dplyr::bind_rows(light1a[,c("lon","lat","species","middle","type_name")],
                                      light2a[,c("lon","lat","species","middle","type_name")]) 

light_jmax_points <- merge(light_jmax_points,light_walker, by=c("lon","lat"),all.x=TRUE)

#when it is 0, it means both are -10000, removed it
l2 <- light_jmax_points %>% ggplot( aes(x=type_name, y=middle)) + geom_boxplot()+
  geom_boxplot(aes(x=type_name, y=pred_jmax25_logr),color="red")+
  geom_point(alpha = 0.6, width = 0.5) +geom_hline( yintercept=0.0, size=0.5)+ ylim(-0.5,1)+
  labs(y="light effect on jmax",x=" ") + theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

#ratio
jvratio <- merge(light_jmax_points[,c("lon","lat","species","middle","pred_jmax25_logr")],light_vcmax_points,by=c("lon","lat","species"),all.x=TRUE)
jvratio$jv_obs <- jvratio$middle.x-jvratio$middle.y
jvratio$jv_pred <- jvratio$pred_jmax25_logr.x-jvratio$pred_vcmax25_logr

l3 <- jvratio %>% ggplot( aes(x=type_name, y=jv_obs)) + geom_boxplot()+
  geom_boxplot(aes(x=type_name, y=jv_pred),color="red")+
  geom_point(alpha = 0.6, width = 0.5) +geom_hline( yintercept=0.0, size=0.5)+ ylim(-0.5,1)+
  labs(y="light effect on jmax/vcmax",x=" ") + theme_classic()+coord_flip()+theme(axis.text=element_text(size=12))

plot_grid(l1,l2,l3,nrow=1,label_size = 15)+
  theme(plot.background=element_rect(fill="white", color="white"))

ggsave(paste("~/data/output_gcme/colin/light1.jpg",sep=""),width = 15, height = 5)
