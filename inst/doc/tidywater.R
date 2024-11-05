## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(tidywater)

## -----------------------------------------------------------------------------
mywater <- define_water(
  ph = 7, temp = 15, alk = 100, tot_hard = 100, na = 100, cl = 80,
  cond = 100,
  toc = 3, uv254 = .02, br = 50
)

## -----------------------------------------------------------------------------
dosed_water <- chemdose_ph(mywater, hcl = 5, alum = 20)
mywater@ph
dosed_water@ph

## -----------------------------------------------------------------------------
coag_water <- chemdose_toc(dosed_water, alum = 20)

dosed_water@doc
coag_water@doc

## -----------------------------------------------------------------------------
caustic_req <- solvedose_ph(coag_water, target_ph = 8.6, chemical = "naoh")

fin_water <- chemdose_ph(coag_water, naoh = caustic_req)

## -----------------------------------------------------------------------------
dist_water <- chemdose_ph(fin_water, naocl = 4) %>%
  chemdose_dbp(cl2 = 4, time = 24, treatment = "coag")

summarize_wq(dist_water, "dbps")

