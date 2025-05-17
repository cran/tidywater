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
# plan(multisession)

## ----warning=FALSE, echo=TRUE-------------------------------------------------
# Read in data from Wells A and B
raw_wells_water <- tibble(
  Well = c("A", "B"),
  ph = c(8, 9),
  alk = c(100, 150),
  temp = c(18, 19),
  ca = c(5, 10),
  cond = c(500, 900),
  tds = c(300, 500),
  na = c(100, 200),
  k = c(0, 20),
  cl = c(0, 30),
  so4 = c(0, 0)
) %>%
  define_water_chain() %>%
  balance_ions_chain()

raw_wells_water

## ----fig.width=7--------------------------------------------------------------
# Ion plot before balance_ions_chain was applied
raw_wells_water$defined_water[[1]] %>%
  plot_ions()
# Plot of balanced ions
raw_wells_water$balanced_water[[1]] %>%
  plot_ions()

## ----warning=FALSE------------------------------------------------------------
# Blend "vertically": blends the data in well A's row with that of well B's.
# The pluck function from the purrr package is useful for indexing a water class column
### First, index the water column using the name or number of the column (ie "balanced_water" or 3 (column number))
### Next, index the row

blended_wells_water <- blend_waters(
  waters = c(
    pluck(raw_wells_water, "balanced_water", 1),
    pluck(raw_wells_water, 3, 2)
  ),
  ratios = c(.5, .5)
)
# outputs a water class object.
blended_wells_water

## ----warning=FALSE------------------------------------------------------------
# Assume wells can contribute 2.5 MGD each
groundwater <- tibble(Wells_flow = c(0, 2.5, 5))
# Blending scenarios and the resulting source water ratios
scenarios <- tibble(
  surface_flow = seq(2, 20, 2),
  River_flow = c(seq(2, 10, 2), rep(10, 5)),
  Lake_flow = c(rep(0, 5), seq(2, 10, 2)),
) %>%
  mutate(group = row_number()) %>%
  cross_join(groundwater) %>%
  mutate(
    total_flow = River_flow + Lake_flow + Wells_flow,
    River_ratio = River_flow / total_flow,
    Lake_ratio = Lake_flow / total_flow,
    Wells_ratio = Wells_flow / total_flow
  )

## ----warning=FALSE------------------------------------------------------------
Wells_water <- tibble(wells = c(blended_wells_water))

River_water <- tibble(
  ph = 7, temp = 20, alk = 200, tds = 950, cond = 1400,
  tot_hard = 300, na = 100, cl = 150, so4 = 200
) %>%
  define_water_chain() %>%
  balance_ions_chain(output_water = "river") %>%
  select(-defined_water)

Lake_water <- tibble(
  ph = 7.5, temp = 19, alk = 180, tds = 900, cond = 1000,
  tot_hard = 350, ca_hard = 250, na = 100, cl = 100, so4 = 150
) %>%
  define_water_chain() %>%
  balance_ions_chain(output_water = "lake") %>%
  select(-defined_water)

## ----warning=FALSE------------------------------------------------------------
blend_water <- scenarios %>%
  cross_join(Wells_water) %>%
  cross_join(River_water) %>%
  cross_join(Lake_water) %>%
  blend_waters_chain(
    waters = c("wells", "river", "lake"),
    ratios = c("Wells_ratio", "River_ratio", "Lake_ratio")
  )

## ----fig.width= 7-------------------------------------------------------------
plotting_data <- blend_water %>%
  pluck_water(input_water = "blended_water", "tot_hard")

# Plot the results!
ggplot(plotting_data, aes(x = total_flow, y = blended_water_tot_hard, color = as.character(Wells_flow))) +
  geom_point() +
  labs(
    y = "Hardness (mg/L as CaCO3)", color = "Well Flow (MGD)",
    x = "Total Plant Flow (MGD)"
  )

## ----warning=FALSE------------------------------------------------------------
# For most operating systems, especially Windows, use this at the beginning of your script
# We recommend removing the `workers` argument to use your computer's full power.
plan(multisession, workers = 2)

# rest of script

# At the end of the script, here's an option to explicitly close the multisession processing
plan(sequential)

