## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ">#"
)
library(tidywater)
library(tidyr)
library(dplyr)
library(ggplot2)
library(furrr)
library(purrr)
# Uncomment the following line for parallel processing.
# plan(multisession)

## ----setup, warning=FALSE-----------------------------------------------------
# Use define_water to prepare for tidywater analysis
no_alum_water <- define_water(ph = 8.3, temp = 18, alk = 150)

# Dose 30 mg/L of alum
alum_30 <- no_alum_water %>%
  chemdose_ph(alum = 30) %>%
  solvedose_ph(target_ph = 8, chemical = "naoh")

alum_30 # Caustic dose required to raise pH to 8 when 30 mg/L of alum is added

# Dose 20 mg/L of alum
alum_20 <- no_alum_water %>%
  chemdose_ph(alum = 20) %>%
  solvedose_ph(target_ph = 8, chemical = "naoh")

alum_20 # Caustic dose required to raise pH to 8 when 20 mg/L of alum is added

## ----warning=FALSE------------------------------------------------------------
# Set a range of alum doses

alum_doses <- tibble(alum = seq(10, 100, 10))

# Use tidywater's built-in synthetic data, water_df, for this example
raw_water <- water_df %>%
  define_water_chain() %>%
  balance_ions_chain() %>%
  # Join alum doses to create several dosing scenarios.
  cross_join(alum_doses)

## ----warning=FALSE------------------------------------------------------------
# 1. Use an existing column in your data frame to dose a chemical.
#    Here, we use the alum column as the dosed chemical.
dose_column_water <- raw_water %>%
  chemdose_ph_chain(input_water = "balanced_water") %>% # The function recognizes the 'alum' column as the chemical dose
  pluck_water(input_water = "dosed_chem_water", parameter = "ph") %>%
  select(-c(defined_water, balanced_water))

head(dose_column_water)

# 2. Dose a chemical in the function. Rename the alum column so it doesn't get used in the function
dose_argument_water <- raw_water %>%
  rename(coagulant = alum) %>%
  chemdose_ph_chain(input_water = "balanced_water", alum = 30) %>%
  pluck_water(input_water = "dosed_chem_water", parameter = "ph") %>%
  select(-c(defined_water, balanced_water))

head(dose_argument_water)

## ----warning=FALSE------------------------------------------------------------
# 1. Use existing columns in your dataframe to set a target pH and the chemicals to dose
raise_ph <- tibble(
  chemical = c("naoh", "mgoh2"),
  target_ph = c(8, 8)
)
solve_column <- raw_water %>%
  chemdose_ph_chain(input_water = "balanced_water") %>%
  cross_join(raise_ph) %>%
  solvedose_ph_once(input_water = "dosed_chem_water") %>%
  select(-c(defined_water:dosed_chem_water))

head(solve_column)

# 2. Set the target pH and chemical needed to raise the pH inside the function
solve_argument <- raw_water %>%
  chemdose_ph_chain(input_water = "balanced_water") %>%
  solvedose_ph_once(input_water = "dosed_chem_water", chemical = "naoh", target_ph = 8) %>%
  select(-c(defined_water:dosed_chem_water))

head(solve_argument)

## ----warning=FALSE------------------------------------------------------------
dosed_caustic_water <- raw_water %>%
  chemdose_ph_chain(input_water = "balanced_water", output_water = "alum_dosed") %>%
  solvedose_ph_once(input_water = "alum_dosed", chemical = "naoh", target_ph = 8) %>%
  rename(
    naoh = dose_required,
    coagulant = alum
  ) %>% # rename alum column so it doesn't get dosed twice
  chemdose_ph_chain(input_water = "alum_dosed", output_water = "caustic_dosed") %>%
  pluck_water(input_water = "caustic_dosed", "ph") %>%
  select(-c(defined_water:chemical))

head(dosed_caustic_water)

## ----warning=FALSE------------------------------------------------------------
enhanced_coag_water <- raw_water %>%
  mutate(alum = 30) %>%
  chemdose_ph_chain(input_water = "balanced_water", output_water = "alum_dosed", hcl = 10) %>%
  pluck_water("alum_dosed", "ph") %>%
  solvedose_ph_once(input_water = "alum_dosed", target_ph = 8, chemical = "naoh", output_column = "naoh") %>%
  select(-c(alum, hcl)) %>% # remove chemical columns so they don't get dosed again in the next line.
  chemdose_ph_chain(input_water = "alum_dosed", output_water = "ph_adjusted") %>%
  select(-c(defined_water:alum_dosed, target_ph, chemical))

head(enhanced_coag_water)

