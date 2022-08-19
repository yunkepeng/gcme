#forcing.R

input all site-level data to download forcing climate, that will be further used in rsofun. The forcing data was downloaded in separate file. Site-name info was also well saved so that it can be used in rsofun.


1.  GCME dataset from Kevin: for co2 effect and warming on Vcmax and Jmax. Warming (soil warming) was not used in this paper finally.
	- All forcing file: ~/data/gcme/kevin/forcing/climate/
	- Sitename info: ~/data/gcme/kevin/forcing/forcing_info.csv

2.  Smith and Keenan GCB dataset sent from Nick: for co2 effect on Vcmax and Jmax. 
    - All forcing file: ~/data/smith_keenan_gcb/gcb_co2/forcing/climate/
    - Sitename info: ~/data/smith_keenan_gcb/gcb_co2/forcing/forcing_info.csv

3.  Kumarathunge et al. publicaly data: for air warming effect on vcmax25 and jmax25
	- All forcing file: ~/data/Kumarathunge_2020_newphy/forcing_warming/
	- Sitename info: ~/data/Kumarathunge_2020_newphy/forcing_info.csv

4. Walker 2014 Light treatment data: for light effect on vcmax and jmax
	- All forcing file: ~/data/leaf_traits/Walker/forcing_light/
	-  Sitename info: ~/data/leaf_traits/Walker/forcing_info.csv


Then for inputting all forcing file at each site, to predict Sensitivity coefficients of vcmax and jmax and their responses to co2, warming and light at each site-level.

5. Output this p-model prediction csv in /Users/yunpeng/data/gcme/MS_data/prediction.csv
