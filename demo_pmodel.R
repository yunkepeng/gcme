sitename <- "a1"
lon <- -1.32
lat <- 51.46
elv <- 124
year_start <- 2008
year_end <- 2012

siteinfo2 <- data.frame(sitename,lon,lat,elv,year_start,year_end)

devtools::load_all("~/yunkepeng/gcme/pmodel/ingestr/")
devtools::load_all("~/yunkepeng/gcme/pmodel/rsofun/")

df_watch <- ingest(
  siteinfo  = siteinfo2,
  source    = "watch_wfdei",
  getvars   = c("temp", "prec", "ppfd", "vpd", "patm"), 
  dir       = "/Volumes/My Passport/data/watch_wfdei/",
  settings  = list(correct_bias = "worldclim", dir_bias = "/Volumes/My Passport/data/worldclim/"))

df_cru <- ingest(
  siteinfo  = siteinfo2,
  source    = "cru",
  getvars   = "ccov",
  dir       = "/Volumes/My Passport/data/cru/ts_4.01/")

df_co2 <- ingest(
  siteinfo  = siteinfo2,
  source  = "co2_mlo",
  verbose = FALSE)

df_co2_final <- as.data.frame(df_co2$data)

df_co2_final2 <- df_co2_final[!(format(df_co2_final$date,"%m") == "02" & format(df_co2_final$date, "%d") == "29"), , drop = FALSE] # make columns consistent

co2 <- df_co2_final2$co2

ddf_meteo <- as_tibble(cbind(as.data.frame(df_watch$data),as.data.frame(df_cru$data)[,c("ccov_int","ccov")],co2))

ddf_meteo$fapar <- 1
ddf_meteo$tmax <- ddf_meteo$temp
ddf_meteo$tmin <- ddf_meteo$temp
head(ddf_meteo)
dim(ddf_meteo)
#then output ddf_meteo as csv


#now, we can work on predicting gpp
df_soiltexture <- bind_rows(
  top    = tibble(layer = "top",    fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1),
  bottom = tibble(layer = "bottom", fsand = 0.4, fclay = 0.3, forg = 0.1, fgravel = 0.1))
params_modl <- list(
  kphio           = 0.09423773,
  soilm_par_a     = 0.33349283,
  soilm_par_b     = 1.45602286,
  tau_acclim_tempstress = 10,
  par_shape_tempstress  = 0.0)
forcing <- ddf_meteo
head(forcing)


siteinfo2 <- data.frame(sitename,lon,lat,elv,year_start,year_end)

modlist <- run_pmodel_f_bysite( 
  sitename="a1", 
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
    firstyeartrend = 2008,
    nyeartrend = 1), 
  siteinfo = siteinfo2, 
  forcing, 
  df_soiltexture, 
  params_modl = params_modl, 
  makecheck = TRUE)
mean(modlist$jmax25)*1000000
