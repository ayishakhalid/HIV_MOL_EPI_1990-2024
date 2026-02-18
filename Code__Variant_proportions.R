###############################################################################
### HIV MOL EPI 1990-2024 - Data analysis - Variant proportions
### Ayisha Khalid
### R Version 4.4.1

### Created: 14 August 2025
###############################################################################


###############################################################################
### Set up ####################################################################

## Libraries
library(tidyverse)
library(readxl)
library(writexl)

options(scipen=999) 



###############################################################################
### Data ##############################################################

database_time <- read_xlsx("~/Desktop/HIV Analysis/Data/HIV_DB_Time_Full_24Sep2025.xlsx", guess_max = 11000)
plhiv <- read_xlsx("~/Desktop/HIV Analysis/Data/UNAIDS_PLHIV.xlsx")

## Group Taiwan and Hong Kong into China
database_time %>% filter(str_detect(sample_collection_sites_country_1, "Taiwan|Hong Kong|China")) %>% group_by(sample_collection_sites_country_1) %>% count()
database_time <- database_time %>% mutate(sample_collection_sites_country_1 = case_when(
  str_detect(sample_collection_sites_country_1, "Taiwan|Hong Kong") ~ "China", 
  TRUE ~ sample_collection_sites_country_1))
database_time %>% filter(str_detect(sample_collection_sites_country_1, "Taiwan|Hong Kong|China")) %>% group_by(sample_collection_sites_country_1) %>% count()



################################################################################
### Distribution of variants ###################################################

## Data 
database_variants <- database_time %>% rename_with(~ gsub("hiv1_group_m_", "hiv1_", .x, fixed = TRUE), starts_with("hiv1"))
database_variants <- database_variants %>% 
  rename(unspec_crfs = unspecified_crf_number, 
         unspec_recombs = unspecified_recombinants, 
         urfs = undefined_urf_number,
         hiv1_a_sum = sum_hiv1_a,
         hiv1_f_sum = sum_hiv1_f,
         other_crfs = sum_other_crfs) 

database_variants <- database_variants %>% relocate(c(other_crfs, urfs, unspec_recombs), .after = crf07_bc)
database_variants <- database_variants %>% relocate(c(hiv1_a_sum, hiv1_f_sum), .before = hiv1_g)
database_variants <- database_variants %>% relocate(c(hiv1_b, hiv1_c, hiv1_d), .after = hiv1_a_sum)
colnames(database_variants)


## Pivot 
database_variants_wider_full <- database_variants %>% 
  pivot_longer(cols = hiv1_a:unspec_recombs,  
               names_to = "variant", 
               values_to = "variant_n", 
               values_drop_na = FALSE)
database_variants_wider <- database_variants_wider_full %>% 
  select(region, sample_collection_sites_country_1, 
         year_period, 
         variant, variant_n) %>%
  rename(country = sample_collection_sites_country_1)
unique(database_variants_wider$variant)


# Variants to keep
variants_keep <- c("hiv1_a_sum", 
                   "hiv1_b", "hiv1_c", "hiv1_d", 
                   "hiv1_f_sum", 
                   "hiv1_g", "hiv1_h", "hiv1_j", "hiv1_k", "hiv1_l", 
                   "crf01_ae", "crf02_ag", "crf07_bc",
                   "other_crfs", 
                   "urfs", "unspec_recombs")


database_variants_wider_key <- database_variants_wider %>%  ## Keeping key variables
  filter(variant %in% variants_keep) %>%
  mutate(variant = factor(variant, levels = variants_keep))



### Country-level ##############################################################

## Variant proportions
variants_country <- database_variants_wider_key %>%
  group_by(region, country, year_period, variant) %>% 
  summarise(variant_n = sum(variant_n, na.rm = TRUE))
variants_country <- variants_country %>%
  group_by(region, country, year_period) %>%
  mutate(total_n = sum(variant_n)) %>%
  mutate(variant_pct = (variant_n/total_n)*100)

## Variance
variants_country <- variants_country %>%
  mutate(
    p_i = variant_n/total_n, # same as variant_pct
    var_country = p_i * (1 - p_i) / total_n) 


## Number of PLHIV
plhiv_mean_country <- plhiv %>% 
  group_by(year_period, region, country) %>% 
  summarise(plhiv_mean_country = mean(plhiv))

variants_plhiv_country <- left_join(variants_country, plhiv_mean_country) 
variants_plhiv_country <- variants_plhiv_country %>%
  mutate(plhiv_variant_country = (variant_pct/100) * plhiv_mean_country) 



### Region-level ##############################################################

## Variant proportions
variants_region <- variants_plhiv_country %>%
  group_by(year_period, region, variant) %>%
  mutate(plhiv_variant_region = sum(plhiv_variant_country, na.rm = TRUE)) 

variants_region <- variants_region %>% 
  group_by(year_period, region) %>%
  mutate(plhiv_region_data = sum(plhiv_variant_country, na.rm = TRUE )) %>% ungroup() 

variants_region <- variants_region %>%
  group_by(year_period, region, variant) %>%
  summarise(
    plhiv_variant_region = unique(plhiv_variant_region), #retain 
    plhiv_region_data = unique(plhiv_region_data), #retain
    variant_pct_region = (plhiv_variant_region / plhiv_region_data)*100,
    
    ## Confidence intervals 
    var_reg = sum( var_country * ( (plhiv_mean_country / unique(plhiv_region_data))^2 ), na.rm = TRUE ),
    se_reg = sqrt(var_reg),
    ci_low = pmax(0, (plhiv_variant_region / plhiv_region_data) - 1.96 * se_reg),
    ci_high = pmin(1, (plhiv_variant_region / plhiv_region_data) + 1.96 * se_reg),
    ci_low_pct_reg  = 100 * ci_low,
    ci_high_pct_reg = 100 * ci_high,
    .groups = "drop")


## Number of PLHIV
plhiv_summeans_region <- plhiv_mean_country %>% 
  group_by(year_period, region) %>%
  summarise(plhiv_region_unaids = sum(plhiv_mean_country))
variants_plhiv_region <- left_join(variants_region, plhiv_summeans_region)  

variants_plhiv_region <- variants_plhiv_region %>%
  mutate(plhiv_variant_region_adj = (variant_pct_region/100)*plhiv_region_unaids)



### Global ##############################################################

## Variant proportions
variants_global <- variants_plhiv_region %>%
  group_by(year_period, variant) %>%
  summarise(
    
    ## Number of PLHIV
    plhiv_variant_global = sum(plhiv_variant_region_adj, na.rm = TRUE),
    plhiv_global = sum(plhiv_region_unaids, na.rm = TRUE),
    
    variant_pct_global = (plhiv_variant_global / plhiv_global) * 100,
    
    # Confidence intervals 
    var_global = sum( var_reg * ( (plhiv_region_unaids / unique(plhiv_global))^2 ), na.rm = TRUE ),
    se_global = sqrt(var_global),
    ci_low = pmax(0, (plhiv_variant_global / plhiv_global) - 1.96 * se_global),
    ci_high = pmin(1, (plhiv_variant_global / plhiv_global) + 1.96 * se_global),
    ci_low_pct_global  = 100 * ci_low,
    ci_high_pct_global = 100 * ci_high,
    .groups = "drop")


