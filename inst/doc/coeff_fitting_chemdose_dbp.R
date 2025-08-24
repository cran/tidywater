## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(tidywater)
library(ggplot2)
library(stats)
library(dplyr)

## -----------------------------------------------------------------------------
influent <- define_water(ph = 7.5, temp = 20, alk = 50, toc = 4, uv254 = .2, br = 50) %>%
  chemdose_dbp(cl2 = 2, time = 8, treatment = "coag", cl_type = "chlorine", location = "plant")

## ----initial-plot, echo=TRUE, message=FALSE, warning=FALSE--------------------
# Given sample data
dbp_data <- data.frame(
  ph = c(7, 7.04, 7.5, 7.9, 7.2, 7.25),
  temp = rep(20, 6),
  alk = c(50, 66, 75, 80, 55, 100),
  toc = c(4, 8, 4, 10, 5, 6),
  uv254 = c(0.1, 0.2, 0.2, 0.1, 0.05, 0.2),
  br = rep(50, 6),
  fin_chcl3 = c(62, 138.5, 71, 206, 83.5, 100)
) %>%
  define_water_df("input")

# Predicted chcl3 concentrations using default coefficients
plot_unfit_dbps <- dbp_data %>%
  chemdose_dbp_df("input", cl2 = 2, time = 8, treatment = "raw", cl_type = "chlorine", location = "plant") %>%
  pluck_water(input_waters = c("disinfected"), parameter = c("chcl3"))

ggplot(plot_unfit_dbps, aes(x = fin_chcl3, y = disinfected_chcl3)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Observed vs predicted chloroform",
    x = "Observed chloroform",
    y = "Default Predicted chloroform"
  ) +
  expand_limits(x = 0, y = 0)

## ----chloroform-plot, echo=TRUE, message=FALSE, warning=FALSE-----------------
# Model we're trying to fit
model_fn <- function(params, data) {
  custom_coeff <- data.frame(
    ID = "chcl3", "A" = params[1], "a" = params[2], "b" = params[3], "c" = params[4], "d" = params[5], "e" = params[6],
    "f" = params[7], ph_const = NA
  )
  temp_data <- data %>%
    chemdose_dbp_df("input",
      cl2 = 2, time = 8, treatment = "raw", cl_type = "chlorine", location = "plant",
      coeff = custom_coeff
    ) %>%
    pluck_water(input_waters = c("disinfected"), parameter = c("chcl3")) %>%
    mutate(diff = (disinfected_chcl3 - fin_chcl3)^2)
  sum(temp_data$diff)
}

# optimize the coefficients

test <- optim(par = c(6.237e-2, 1.617, -0.094, -0.175, 0.607, 1.403, 0.306), fn = model_fn, data = dbp_data, control = list(maxit = 100))
coeffs <- data.frame(
  ID = "chcl3",
  A = test$par[1],
  a = test$par[2],
  b = test$par[3],
  c = test$par[4],
  d = test$par[5],
  e = test$par[6],
  f = test$par[7],
  ph_const = NA
)

dbp_data <- dbp_data %>%
  chemdose_dbp_df(input_water = "input", output_water = "coeff_disinfected", cl2 = 2, time = 8, coeff = coeffs) %>%
  pluck_water(input_waters = c("coeff_disinfected"), parameter = c("chcl3"))

ggplot() +
  geom_point(data = plot_unfit_dbps, aes(x = fin_chcl3, y = disinfected_chcl3, color = "Default coeffs")) +
  geom_point(data = dbp_data, aes(x = fin_chcl3, y = coeff_disinfected_chcl3, color = "Custom coeffs")) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") +
  labs(
    title = "Observed vs predicted chloroform",
    x = "Observed chloroform",
    y = "Predicted chloroform"
  ) +
  scale_color_manual(
    name = "Prediction type",
    breaks = c("Default coeffs", "Custom coeffs"),
    values = c("black", "blue")
  ) +
  expand_limits(x = 0, y = 0)

