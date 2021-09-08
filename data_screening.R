#the processing code of gcme was deried from: /Users/yunpeng/yunkepeng/cnreview/dataanalysis.Rmd

library(readr)
library(dplyr)
library(metafor)  
library(ggplot2)
library(stringr)

df <- read_csv("~/data/gcme/data_received_190325/NewData_wide_CORRECTED2.csv") %>%
  
  # something is wrong in columns ambient_Sd, ambient_Se, elevated...
  mutate( ambient_Sd  = as.numeric(ambient_Sd),  ambient_Se  = as.numeric(ambient_Se), 
          elevated_Sd = as.numeric(elevated_Sd), elevated_Se = as.numeric(elevated_Se) )

# XXX testing
# filter(exp_nam == "SwissFACE_lolium2" & Data_type == "aboveground_biomass")
# filter(exp_nam == "SwissFACE_trifolium2" & Data_type == "aboveground_biomass")

# save experiments names
df_experiments <- df %>% select(exp_nam, prev_name) %>% distinct()

### Duplicate rows
#Some rows are simply duplicates. This is dangerous.
keyvars <- c("exp_nam", "factors", "treatment", "Data_type", "Unit", "Sampling_date", "Start_Year", "Year", "n_plots")
valuevars <-  c("ambient", "ambient_Se", "ambient_Sd", "elevated", "elevated_Se", "elevated_Sd")

df %>% 
  filter(exp_nam == "SwissFACE_trifolium2" & Data_type == "aboveground_biomass") %>% 
  filter(treatment=="cf") %>% 
  select(one_of(keyvars, valuevars)) %>% 
  print()

#This may be remediated automatically by using only columns that are distinct w.r.t. `valuevars` and `keyvars` (basically all columns in original data, except `ALIAS`).
print(nrow(df))
df <- df %>% 
  distinct_at(vars(one_of(keyvars, valuevars))) %>%
  mutate( id=1:n() )   # create new ID key (before this was 'ALIAS')
print(nrow(df))

df %>% 
  filter(exp_nam == "SwissFACE_trifolium2" & Data_type == "aboveground_biomass") %>% 
  filter(treatment=="cf") %>% 
  select(one_of(keyvars, valuevars)) %>% 
  print()

#Still problems?
df %>% 
  filter(exp_nam == "SwissFACE_lolium2" & Data_type == "belowground_biomass") %>% 
  filter(treatment=="cf", Year==6) %>% 
  select(one_of(keyvars, valuevars)) %>% 
  print()

# Distinct ambients
df %>% 
  filter(exp_nam == "SwissFACE_lolium2" & Data_type == "Rh" ) %>% 
  print()

#How to perform a conversion of our current wide format to a long format? Let's do it for one example. An example for how the wide format is implemented in the current dataset:
df_wide <- df %>% 
  filter(exp_nam=="SwissFACE_lolium2", Data_type=="aboveground_biomass") %>%
  select(exp_nam, factors, treatment, Data_type) %>% 
  distinct() %>%
  mutate( id=1:nrow(.),
          ambient=rep("A", 3), elevated=c("B", "C", "D"),
          ambient_Sd=rep("sdA", 3), elevated_Sd=c("sdB", "sdC", "sdD"),
          ambient_Se=rep("seA", 3), elevated_Se=c("seB", "seC", "seD")) %>%
  print()
#To convert this to a long format, first get all available factors and create a vector. The available factors are:
df %>% 
  filter(!(factors %in% c("x","x99"))) %>% 
  select(factors) %>% 
  unique() %>% 
  unlist() %>% 
  unname() %>% 
  print()
#Create a vector with all single factors
factors_avl <- df %>%
  filter(!(factors %in% c("x","x99"))) %>%      # Something weird with this
  select(factors) %>% 
  unique() %>% 
  mutate( factors_sep = strsplit(factors, NULL) ) %>% 
  select(factors_sep) %>% 
  unlist() %>% 
  unname() %>% 
  unique() %>% 
  print()
#Convert from wide to long, adding columns for all factors for which we have ambient and elevated levels' data.
source("~/yunkepeng/cnreview/R/wide_to_long_gcme.R")
# run it
keyvars <- c("id", "exp_nam", "factors", "treatment", "Data_type")
valuevars <-  c("ambient", "ambient_Se", "ambient_Sd", "elevated", "elevated_Se", "elevated_Sd")
df_long_example <- df_wide %>%
  select(one_of(keyvars, valuevars)) %>%
  wide_to_long_gcme( keyvars ) %>%
  print()

### Converting the entire table
keyvars <- c("id", "exp_nam", "factors", "treatment", "level", "Data_type", "Unit", "Start_Year", "Year",  "n_plots")
valuevars <-  c("ambient", "ambient_Se", "ambient_Sd", "elevated", "elevated_Se", "elevated_Sd")
df_long <- df %>% 
  select(one_of(keyvars, valuevars)) %>%
  # filter(exp_nam=="BioCON" & Data_type=="leaf_N" & Year==4) %>%  # xxx debug
  wide_to_long_gcme(keyvars) %>% 
  write_csv(path="~/data/gcme/data_received_190325/NewData_LONG.csv")

#Thus, we have a long-format data table, that has *less* than twice as many rows as the original table, because the "absolute ambient" is no longer repeated:
print(nrow(df)*2)
nrow(df_long)

#### Test
check_cf <- function(df_wide, df_long, verbose=FALSE){
  
  ## check cf experiments
  df_cf <- df_wide %>% filter(treatment %in% c("f","c","cf","fc"))
  
  out <- tibble()
  
  all_exp <- df_cf$exp_nam %>% unique
  
  for (iexp in all_exp){
    
    df_sub <- df_cf %>% filter(exp_nam==iexp)
    
    all_var <- df_sub$Data_type %>% unique
    
    for (ivar in all_var){
      
      df_subsub <- df_sub %>% filter(Data_type==ivar)
      
      all_var_unit <- df_subsub$Unit %>% unique
      
      for (iunit in all_var_unit){
        
        df_subsubsub <- df_subsub %>% filter(Unit==iunit)
        
        # elevated-f treatment
        correct_f <- df_subsubsub %>% 
          filter(treatment=="f") %>% 
          select(exp_nam, Data_type, Unit, Year, mean=elevated, se=elevated_Se, sd=elevated_Sd) %>%
          all_equal( ., df_long %>% filter(f & !c & exp_nam==iexp & Data_type==ivar & Unit==iunit) %>% select(exp_nam, Data_type, Unit, Year, mean, se, sd) )
        
        # elevated-c treatment
        correct_c <- df_subsubsub %>% 
          filter(treatment=="c") %>% 
          select(exp_nam, Data_type, Unit, Year, mean=elevated, se=elevated_Se, sd=elevated_Sd) %>%
          all_equal( ., df_long %>% filter(!f & c & exp_nam==iexp & Data_type==ivar & Unit==iunit) %>% select(exp_nam, Data_type, Unit, Year, mean, se, sd) )
        
        # elevated-cf treatment
        correct_cf <- df_subsubsub %>% 
          filter(treatment=="cf") %>% 
          select(exp_nam, Data_type, Unit, Year, mean=elevated, se=elevated_Se, sd=elevated_Sd) %>%
          all_equal( ., df_long %>% filter(f & c & exp_nam==iexp & Data_type==ivar & Unit==iunit) %>% select(exp_nam, Data_type, Unit, Year, mean, se, sd) )
        
        # (absolute) control: WARNING: DUPLICATED ROWS
        correct_0 <- df_subsubsub %>% 
          distinct(exp_nam, Data_type, Unit, Year, ambient, ambient_Se, ambient_Sd, .keep_all=TRUE) %>% 
          select(exp_nam, Data_type, Unit, Year, mean=ambient, se=ambient_Se, sd=ambient_Sd) %>%
          all_equal( ., df_long %>% filter(!f & !c & exp_nam==iexp & Data_type==ivar & Unit==iunit) %>% select(exp_nam, Data_type, Unit, Year, mean, se, sd) )
        
        if (verbose){
          if (!identical(correct_0, TRUE)){
            rlang::warn(correct_0)
            print("WIDE:")
            df_subsubsub %>%
              distinct(exp_nam, Data_type, Unit, Year, ambient, ambient_Se, ambient_Sd, .keep_all=TRUE) %>%
              select(exp_nam, Data_type, Unit, Year, mean=ambient, se=ambient_Se, sd=ambient_Sd) %>%
              print()
            print("LONG:")
            df_long %>%
              filter(!f & !c & exp_nam==iexp & Data_type==ivar & Unit==iunit) %>%
              select(exp_nam, Data_type, Unit, Year, mean, se, sd) %>%
              print()
            print("---------")
          }
        }
        
        out <- bind_rows(out, c(factors="cf", exp_nam=iexp, Data_type=ivar, Unit=iunit, correct_0=correct_0, correct_c=correct_c, correct_f=correct_f, correct_cf=correct_cf))
        
      }
    }
  }
  
  return(out)
}

#There are still lots of problems with identifying the ambients (see column `correct_0` in `out`). See for example:
## note row 18 is missing in wide
out <- check_cf( filter(df, exp_nam=="BioCON", Data_type=="leaf_N", Year==4), df_long, verbose = TRUE)

#- Make sure the wide-to-long conversion works fine. That is, `out` should have `TRUE` in all columns. 
#Let's implement this, by converting our example back to a CO2-wide format:
long_to_wide_gcme <- function(df_long, keyvar){

  joinvars <- names(df_long)[-which(names(df_long) %in% c("mean", "sd", "se", keyvar))]
  factors_all <- df_long %>%
    filter(!(factors %in% c("x","x99"))) %>%      # Something weird with this
    select(factors) %>% 
    unique() %>% 
    mutate( factors_sep = strsplit(factors, NULL) ) %>% 
    select(factors_sep) %>% 
    unlist() %>% 
    unname() %>% 
    unique()
  factorvars <- factors_all[-which(factors_all==keyvar)]
    
  # take all data where factor 'c' is TRUE
  df_wide <- dplyr::filter(df_long, eval(parse(text=keyvar)) ) %>% 
    
    # call this 'elevated'
    dplyr::rename(elevated=mean, elevated_sd=sd, elevated_se=se) %>% 
    
    # remove column 'c', is no longer used
    dplyr::select(-keyvar) %>% 
    
    # merge this with the corresponding row where all other factors are the same, while 'c' is FALSE
    left_join( 
      
      # take all data where factor 'c' is FALSE
      dplyr::filter(df_long, !(eval(parse(text=keyvar)))) %>% 

        # call this 'ambient'
        dplyr::rename(ambient=mean, ambient_sd=sd, ambient_se=se) %>% 

        # remove column 'c', is no longer used
        dplyr::select(-keyvar),
      
      ## merge by all other columns
      by = joinvars ) %>% 
    
    # order columns
    select( c(joinvars, ambient, elevated, ambient_se, elevated_se, ambient_sd, elevated_sd)  )

  return(df_wide)
}

df_long_example %>%
  long_to_wide_gcme("c") %>%
  print()

#Let's apply this to the full data frame.
keyvars <- c("ALIAS", "exp_nam", "factors", "Data_type", "Unit", "Sampling_date", "n_plots")
valuevars <-  c("mean", "sd", "se")
df_wide_c <- df_long %>% 
  # select(one_of(keyvars, valuevars, factors_avl)) %>%
  # select(-prev_name, -TT_Nut_Detail) %>% 
  long_to_wide_gcme("c") %>% 
  write_csv(path="~/data/gcme/data_received_190325/created/NewData_WIDE_c.csv")


df_wide_f <- df_long %>% 
  # select(one_of(keyvars, valuevars, factors_avl)) %>%
  # select(-prev_name, -TT_Nut_Detail) %>% 
  long_to_wide_gcme("f") %>% 
  write_csv(path="~/data/gcme/data_received_190325/created/NewData_WIDE_f.csv")

#data selection
## Determine all experiments that have more than 1 year data
list_exp_gt1yr <- df %>% 
  filter(!is.na(Year)) %>% 
  group_by(exp_nam) %>% 
  summarise(nyears=max(Year)) %>% 
  filter(nyears>1) %>% 
  select(exp_nam) %>% 
  unlist() %>% 
  unname()

df_c <- df %>%
  
  # ## Take this info from experiments below
  # select(-Fumigation_type, -Vegetation_type) %>% 
  
  ## Add prev_name back
  left_join( df_experiments, by = "exp_nam") %>% 
  
  ## filter experiments with only manipulated CO2
  ## (no other factors manipulated, strong reduction of data)
  filter(treatment=="c") %>% 
  
  ## More than 1 year data
  filter(exp_nam %in% list_exp_gt1yr) %>% 
  
  ## Combine with experiments meta info
  left_join( 
    read_csv("~/data/gcme/data_received_190325/table_var_exp_names_experiments.csv") %>% 
      select(prev_name, Fumigation_type=my_fumigation_type, Vegetation_type),
    by = c("prev_name")
  ) %>% 
  
  ## Filter only Fumigation_type OTC or FACE
  filter( Fumigation_type %in% c("OTC", "FACE") ) %>%
  
  {.}


#now, leaf n 
df_leafn <- read_csv("~/data/gcme/data_received_190325/table_var_exp_names_leafn.csv") %>% 
  select(prev_name=expnam, varnam, Data_type, Source_Reference,ambient,elevated) %>%
  distinct() %>% 
  filter(varnam %in% c("nmass_leaf", "narea_leaf", "rcton_leaf") ) %>% 
  
  ## add 'exp_nam' (new experiment name) from data table
  left_join( df_experiments, by="prev_name" ) %>% 
  
  ## filter only experiments selected by general filters above
  filter( exp_nam %in% unique(df_c$exp_nam) ) %>% 
  
  ## add experiments meta information from experiments table
  left_join( 
    read_csv("~/data/gcme/data_received_190325/table_var_exp_names_experiments.csv") %>% 
      select(prev_name, Fumigation_type=my_fumigation_type, Cquantity_Control, Cquantity_Treatment, Vegetation_type, Plants_specs, Start_Year, End_Year, Remarks),
    by = c("prev_name")
  )
