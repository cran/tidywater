## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(tidywater)

## ----echo=TRUE----------------------------------------------------------------
empty_water <- define_water()

## ----echo=TRUE----------------------------------------------------------------
print(empty_water)

## ----warning=FALSE------------------------------------------------------------
my_water <- define_water(ph = 7.5, alk = 100, temp = 20, na = 50, ca = 50)
my_water
my_water@na
my_water@hco3

## -----------------------------------------------------------------------------
summarize_wq(my_water)

## ----echo=TRUE, fig.width=7---------------------------------------------------
plot_ions(my_water)

## ----, warning=FALSE, fig.width=7---------------------------------------------
balanced_water <- my_water %>% balance_ions()
plot_ions(balanced_water)

## -----------------------------------------------------------------------------
my_water@cl # We did not input any chloride in the original water

balanced_water@cl # The balanced water now contains chloride

## -----------------------------------------------------------------------------
convert_units(value = balanced_water@cl, formula = "cl", startunit = "M", endunit = "mg/L")

## -----------------------------------------------------------------------------
# The ionic strength slot was NA in the original water because we did not
# provide enough information to calculate it
my_water@is

# Input TDS or conductivity
new_water1 <- define_water(ph = 7.5, alk = 100, temp = 20, na = 50, ca = 50, tds = 100)
new_water2 <- define_water(ph = 7.5, alk = 100, temp = 20, na = 50, ca = 50, cond = 200)
# Input more known ions
new_water3 <- define_water(ph = 7.5, alk = 100, temp = 20, na = 50, ca = 50, so4 = 100)

new_water1@is
new_water2@is
new_water3@is

## ----warning=FALSE------------------------------------------------------------
# Calculate hardness or calcium hardness
hard_water <- define_water(8, 20, 100, tot_hard = 150)

# total hardness in mg/L CaCO3
hard_water@tot_hard

# calcium hardness
convert_units(value = hard_water@ca, formula = "ca", startunit = "M", endunit = "mg/L CaCO3")

# magnesium hardness
convert_units(value = hard_water@mg, formula = "mg", startunit = "M", endunit = "mg/L CaCO3")

## ----warning=FALSE------------------------------------------------------------
# Calculate TOC and DOC
toc_water <- define_water(8, 20, 100, toc = 3)
toc_water@toc # mg/L
toc_water@doc # mg/L

doc_water <- define_water(8, 20, 100, doc = 1.3)
doc_water@toc # mg/L
doc_water@doc # mg/L

