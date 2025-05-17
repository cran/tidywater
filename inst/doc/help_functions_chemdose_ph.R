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

alum_doses <- tibble(alum_dose = seq(20, 60, 10))

# Use tidywater's built-in synthetic data, water_df, for this example
raw_water <- water_df %>%
  slice_head(n = 2) %>%
  define_water_chain("raw") %>%
  # Join alum doses to create several dosing scenarios.
  cross_join(alum_doses)

## ----warning=FALSE------------------------------------------------------------
dose_water <- raw_water %>%
  mutate(hcl = 10) %>%
  chemdose_ph_chain(input_water = "raw", alum = alum_dose) %>%
  pluck_water(input_water = c("raw", "dosed_chem_water"), parameter = "ph") %>%
  select(-c(raw, dosed_chem_water))

head(dose_water)

dose_water <- raw_water %>%
  chemdose_ph_chain(input_water = "raw", alum = alum_dose, hcl = 5) %>%
  pluck_water(input_water = c("raw", "dosed_chem_water"), parameter = "ph") %>%
  select(-c(raw, dosed_chem_water))

head(dose_water)

## ----warning=FALSE------------------------------------------------------------
solve_ph <- raw_water %>%
  chemdose_ph_chain("raw", alum = alum_dose) %>%
  mutate(target_ph = 8) %>%
  solvedose_ph_once(input_water = "dosed_chem_water", chemical = c("naoh", "mgoh2")) %>%
  select(-c(raw, dosed_chem_water))

head(solve_ph)

## ----warning=FALSE------------------------------------------------------------
dosed_caustic_water <- raw_water %>%
  chemdose_ph_chain(input_water = "raw", output_water = "alum_dosed", alum = alum_dose) %>%
  solvedose_ph_once(input_water = "alum_dosed", chemical = "naoh", target_ph = 8) %>%
  chemdose_ph_chain(input_water = "alum_dosed", output_water = "caustic_dosed", naoh = dose_required) %>%
  pluck_water(input_water = "caustic_dosed", "ph")

head(dosed_caustic_water)

