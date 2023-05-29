#The input data to run all figures of chapter 2 MS is available in Euler:

("~/data/gcme/MS_data/plot_data.csv")


exp: experimental name

lon: longitude

lat: latitude

condition: light = light treatment (divided into low vs. high light, or shading vs. sun; see typename), co2 = co2 treatment, warming  = warming treatment, Fertilization = Fertilization effect (but still only evalualtes co2 effect). High and low N = high vs. low fertilization effect (but still only evalualtes co2 effect)

	Note: For fertilization effect: There are 5 paralleled sites that measures eCO2-Vcmax effect at Nfertilizations vs. no fertilizations effect. There are another 7 paralleled sites that measures eCO2-Vcmax effect at high N vs. low N. See 'type name' variable.



All below values are sensitivity coefficient, but definitions different between treamtents: 

For condition == co2, sensitivity coefficients = log(vcmax-ele/vcmax-amb)/log(co2-elv/co2-amb)

For condition == light, sensitivity coefficients = log(vcmax-ele/vcmax-amb)

For condition == warming, sensitivity coefficients = log(vcmax25-ele/vcmax25-amb)/ (T-elv - Tamb)

For condition == Fertilization, sensitivity coefficients = log (Vcmax[CO2 + N] / Vcmax[N]) / log (CO2[ele] / CO2 [amb]) 

Below are variables full name - all dimensionless - since it is sensitivity coefficient

vcmax: maximum rate of carbyxolation capacity

jmax: maximum rate of electron transport

jmax_vcmax: jmax/vcmax

pred_vcmax: predicted maximum rate of carbyxolation capacity

pred_jmax: predicted maximum rate of electron transport
 
pred_jmax_vcmax: predicted jmax/vcmax

ecosystem: forest, grassland, cropland

ref: reference

LMA: leaf-mass-per-area

narea: leaf nitrogen per area

nmass: leaf mass per area

bnpp: belowground net primary production

root_shoot_ratio: root/shoot ratio

soil_mineral_N: soil mineral N

anpp: abovegruond net primary production

lai: leaf area index

ecm_type: to divide categories as N-fixing or non N-fixing species.
