# Chlorine/Chloramine Decay Modeling functions
# These functions predict chlorine residual concentration given reaction time

#' @title Calculate chlorine decay
#'
#' @description calculates the decay of chlorine or chloramine based on the U.S. EPA's
#' Water Treatment Plant Model (U.S. EPA, 2001).
#'
#' @details Required arguments include an object of class "water" created by \code{\link{define_water}},
#' applied chlorine/chloramine dose, type, reaction time, and treatment applied (options include "raw" for
#' no treatment, or "coag" for coagulated water). The function also requires additional water quality
#' parameters defined in \code{define_water} including TOC and UV254. The output is a new "water" class
#' with the calculated total chlorine value stored in the 'free_chlorine' or 'combined_chlorine' slot,
#' depending on what type of chlorine is dosed. When modeling residual concentrations
#' through a unit process, the U.S. EPA Water Treatment Plant Model applies a correction factor based on the
#' influent and effluent residual concentrations (see U.S. EPA (2001) equation 5-118) that may need to be
#' applied manually by the user based on the output.
#'
#' @source U.S. EPA (2001)
#' @source See references list at: \url{https://github.com/BrownandCaldwell-Public/tidywater/wiki/References}
#'
#' @param water Source water object of class "water" created by \code{\link{define_water}}
#' @param cl2_dose Applied chlorine or chloramine dose (mg/L as cl2). Model results are valid for doses between 0.995 and 41.7 mg/L for raw water,
#' and for doses between 1.11 and 24.7 mg/L for coagulated water.
#' @param time Reaction time (hours). Chlorine decay model results are valid for reaction times between 0.25 and 120 hours.Chloramine decay model
#' does not have specified boundary conditions.
#' @param treatment Type of treatment applied to the water. Options include "raw" for no treatment (default), "coag" for
#' water that has been coagulated or softened.
#' @param cl_type Type of chlorination applied, either "chlorine" (default) or "chloramine".
#' @examples
#' example_cl2 <- suppressWarnings(define_water(8, 20, 66, toc = 4, uv254 = 0.2)) %>%
#'   chemdose_chlordecay(cl2_dose = 2, time = 8)
#' @export
#' @returns An updated disinfectant residual in the free_chlorine or combined chlorine water slot in units of M.
#' Use \code{\link{convert_units}} to convert to mg/L.
#'
chemdose_chlordecay <- function(water, cl2_dose, time, treatment = "raw", cl_type = "chlorine") {
  validate_water(water, c("toc", "uv254"))

  toc <- water@toc
  uv254 <- water@uv254

  # Handle missing arguments with warnings (not all parameters are needed for all models).
  if (missing(cl2_dose)) {
    stop("Missing value for chlorine dose. Please check the function inputs required to calculate chlorine/chloramine decay.")
  }

  if (missing(time)) {
    stop("Missing value for reaction time. Please check the function inputs required to calculate chlorine/chloramine decay.")
  }

  if (!(cl_type %in% c("chlorine", "chloramine"))) {
    stop("cl_type should be 'chlorine' or 'chloramine'. Please check the spelling for cl_type to calculate chlorine/chloramine decay.")
  }

  # chlorine decay model
  if (cl_type == "chlorine") {
    if (!(treatment %in% c("raw", "coag"))) {
      stop("The treatment type should be 'raw' or 'coag'. Please check the spelling for treatment.")
    }

    # toc warnings
    if (treatment == "raw" & (toc < 1.2 | toc > 16)) {
      warning("TOC is outside the model bounds of 1.2 <= TOC <= 16 mg/L for raw water.")
    }

    if (treatment == "coag" & (toc < 1.0 | toc > 11.1)) {
      warning("TOC is outside the model bounds of 1.0 <= TOC <= 11.1 mg/L for coagulated water.")
    }

    # uv254 warnings
    if (treatment == "raw" & (uv254 < 0.010 | uv254 > 0.730)) {
      warning("UV254 is outside the model bounds of 0.010 <= UV254 <= 0.730 cm-1 for raw water.")
    }

    if (treatment == "coag" & (uv254 < 0.012 | uv254 > 0.250)) {
      warning("UV254 is outside the model bounds of 0.012 <= UV254 <= 0.250 cm-1 for coagulated water.")
    }

    # cl2_dose warnings
    if (treatment == "raw" & (cl2_dose < 0.995 | cl2_dose > 41.7)) {
      warning("Chlorine dose is outside the model bounds of 0.995 <= cl2_dose <= 41.7 mg/L for raw water.")
    }

    if (treatment == "coag" & (cl2_dose < 1.11 | cl2_dose > 24.7)) {
      warning("Chlorine dose is outside the model bounds of 1.11 <= cl2_dose <= 24.7 mg/L for coagulated water.")
    }

    # time warning
    if (time < 0.25 | time > 120) {
      warning("For chlorine decay estimate, reaction time is outside the model bounds of 0.25 <= time <= 120 hours.")
    }

    # get coefficients from defined clcoeffs table
    if (treatment == "raw") {
      coeffs <- subset(tidywater::cl2coeffs, treatment == "chlorine_raw")
    } else if (treatment == "coag") {
      coeffs <- subset(tidywater::cl2coeffs, treatment == "chlorine_coag")
    }

    # define function for chlorine decay
    # U.S. EPA (2001) equation 5-113 (raw) and equation 5-117 (coag)
    solve_decay <- function(ct, a, b, cl2_dose, uv254, time, c, toc) {
      a * cl2_dose * log(cl2_dose / ct) - b * (cl2_dose / uv254)^c * toc * time + cl2_dose - ct
    }

    # chloramine decay model
  } else if (cl_type == "chloramine") {
    # define function for chloramine decay
    # U.S. EPA (2001) equation 5-120
    solve_decay <- function(ct, a, b, cl2_dose, uv254, time, c, toc) {
      a * cl2_dose * log(cl2_dose / ct) - b * uv254 * time + cl2_dose - ct
    }

    coeffs <- subset(tidywater::cl2coeffs, treatment == "chloramine")
  }

  # if dose is 0, do not run uniroot function
  if (cl2_dose == 0) {
    ct <- 0
  } else {
    root_ct <- stats::uniroot(solve_decay,
      interval = c(0, cl2_dose),
      a = coeffs$a,
      b = coeffs$b,
      c = coeffs$c,
      cl2_dose = cl2_dose,
      uv254 = uv254,
      toc = toc,
      time = time,
      tol = 1e-14
    )

    ct <- root_ct$root
  }

  # Convert final result to molar
  if (cl_type == "chlorine") {
    water@free_chlorine <- convert_units(ct, "cl2", "mg/L", "M")
  } else if (cl_type == "chloramine") {
    water@combined_chlorine <- convert_units(ct, "cl2", "mg/L", "M")
  }


  return(water)
}



#' Apply `chemdose_chlordecay`function within a data frame and output a data frame
#'
#' This function allows \code{\link{chemdose_chlordecay}} to be added to a piped data frame.
#' Its output is a data frame containing columns for free_chlorine or combined_chlorine (depending on chlorine type).
#'
#' The data input comes from a `water` class column, as initialized in \code{\link{define_water_chain}}.
#'
#' If the input data frame has a chlorine dose column (cl2) or time column (time), the function will use those columns. Note:
#' The function can only take cl2 and time inputs as EITHER a column or as function arguments, not both.
#'
#' tidywater functions cannot be added after this function because they require a `water` class input.
#'
#' For large datasets, using `fn_once` or `fn_chain` may take many minutes to run. These types of functions use the furrr package
#' for the option to use parallel processing and speed things up. To initialize parallel processing, use
#' `plan(multisession)` or `plan(multicore)` (depending on your operating system) prior to your piped code with the
#' `fn_once` or `fn_chain` functions. Note, parallel processing is best used when your code block takes more than a minute to run,
#' shorter run times will not benefit from parallel processing.
#'
#' @param df a data frame containing a water class column, which has already been computed using
#' \code{\link{define_water_once}}. The df may include a column named for the applied chlorine dose (cl2),
#' and a column for time in hours.
#' @param input_water name of the column of water class data to be used as the input for this function. Default is "defined_water".
#' @param cl2_dose Applied chlorine or chloramine dose (mg/L as cl2). Model results are valid for doses between 0.995 and 41.7 mg/L for raw water,
#' and for doses between 1.11 and 24.7 mg/L for coagulated water.
#' @param time Reaction time (hours). Chlorine decay model results are valid for reaction times between 0.25 and 120 hours. Chloramine decay model
#' does not have specified boundary conditions.
#' @param treatment Type of treatment applied to the water. Options include "raw" for no treatment (default), "coag" for
#' water that has been coagulated or softened.
#' @param cl_type Type of chlorination applied, either "chlorine" (default) or "chloramine".
#'
#' @seealso \code{\link{chemdose_chlordecay}}
#'
#' @examples
#'
#' library(purrr)
#' library(furrr)
#' library(tidyr)
#' library(dplyr)
#'
#' example_df <- water_df %>%
#'   mutate(br = 50) %>%
#'   define_water_chain() %>%
#'   balance_ions_chain() %>%
#'   chemdose_chlordecay_once(input_water = "balanced_water", cl2_dose = 4, time = 8)
#'
#' example_df <- water_df %>%
#'   mutate(br = 50) %>%
#'   define_water_chain() %>%
#'   balance_ions_chain() %>%
#'   mutate(
#'     cl2_dose = seq(2, 24, 2),
#'     time = 30
#'   ) %>%
#'   chemdose_chlordecay_once(input_water = "balanced_water")
#'
#' example_df <- water_df %>%
#'   mutate(br = 80) %>%
#'   define_water_chain() %>%
#'   balance_ions_chain() %>%
#'   mutate(time = 8) %>%
#'   chemdose_chlordecay_once(
#'     input_water = "balanced_water", cl2_dose = 6, treatment = "coag",
#'     cl_type = "chloramine"
#'   )
#' \donttest{
#' # Initialize parallel processing
#' plan(multisession, workers = 2) # Remove the workers argument to use all available compute
#' example_df <- water_df %>%
#'   mutate(br = 50) %>%
#'   define_water_chain() %>%
#'   balance_ions_chain() %>%
#'   chemdose_chlordecay_once(input_water = "balanced_water", cl2_dose = 4, time = 8)
#'
#' # Optional: explicitly close multisession processing
#' plan(sequential)
#' }
#'
#' @import dplyr
#' @importFrom tidyr unnest
#' @export
#'
#' @returns A data frame with updated chlorine residuals.

chemdose_chlordecay_once <- function(df, input_water = "defined_water", cl2_dose = 0, time = 0,
                                     treatment = "raw", cl_type = "chlorine") {
  temp_cl2 <- chlor <- NULL # Quiet RCMD check global variable note
  output <- df %>%
    chemdose_chlordecay_chain(
      input_water = input_water, output_water = "temp_cl2",
      cl2_dose, time, treatment, cl_type
    ) %>%
    mutate(chlor = furrr::future_map(temp_cl2, convert_water)) %>%
    unnest(chlor) %>%
    select(-temp_cl2)
}

#' Apply `chemdose_chlordecay` within a data frame and output a column of `water` class to be chained to other tidywater functions
#'
#' This function allows \code{\link{chemdose_chlordecay}} to be added to a piped data frame.
#' Its output is a `water` class, and can therefore be used with "downstream" tidywater functions.
#' free_chlorine or combined_chlorine slots will be updated depending on chlorine type.
#'
#' The data input comes from a `water` class column, as initialized in \code{\link{define_water_chain}}.
#'
#' If the input data frame has a chlorine dose column (cl2_dose) or time column (time), the function will use those columns. Note:
#' The function can only take cl2_dose and time inputs as EITHER a column or as function arguments, not both.
#'
#'  For large datasets, using `fn_once` or `fn_chain` may take many minutes to run. These types of functions use the furrr package
#'  for the option to use parallel processing and speed things up. To initialize parallel processing, use
#'  `plan(multisession)` or `plan(multicore)` (depending on your operating system) prior to your piped code with the
#'  `fn_once` or `fn_chain` functions. Note, parallel processing is best used when your code block takes more than a minute to run,
#'  shorter run times will not benefit from parallel processing.
#'
#' @param df a data frame containing a water class column, which has already been computed using
#' \code{\link{define_water_chain}}. The df may include a column named for the applied chlorine dose (cl2_dose),
#' and a column for time in hours.
#' @param input_water name of the column of water class data to be used as the input for this function. Default is "defined_water".
#' @param output_water name of the output column storing updated parameters with the class, water. Default is "disinfected_water".
#' @param cl2_dose Applied chlorine or chloramine dose (mg/L as cl2). Model results are valid for doses between 0.995 and 41.7 mg/L for raw water,
#' and for doses between 1.11 and 24.7 mg/L for coagulated water.
#' @param time Reaction time (hours). Chlorine decay model results are valid for reaction times between 0.25 and 120 hours. Chloramine decay model
#' does not have specified boundary conditions.
#' @param treatment Type of treatment applied to the water. Options include "raw" for no treatment (default), "coag" for
#' water that has been coagulated or softened.
#' @param cl_type Type of chlorination applied, either "chlorine" (default) or "chloramine".
#'
#' @seealso \code{\link{chemdose_chlordecay}}
#'
#' @examples
#'
#' library(purrr)
#' library(furrr)
#' library(tidyr)
#' library(dplyr)
#'
#' example_df <- water_df %>%
#'   mutate(br = 50) %>%
#'   define_water_chain() %>%
#'   balance_ions_chain() %>%
#'   chemdose_chlordecay_chain(input_water = "balanced_water", cl2_dose = 4, time = 8)
#'
#' example_df <- water_df %>%
#'   mutate(br = 50) %>%
#'   define_water_chain() %>%
#'   balance_ions_chain() %>%
#'   mutate(
#'     cl2_dose = seq(2, 24, 2),
#'     time = 30
#'   ) %>%
#'   chemdose_chlordecay_chain(input_water = "balanced_water")
#'
#' example_df <- water_df %>%
#'   mutate(br = 80) %>%
#'   define_water_chain() %>%
#'   balance_ions_chain() %>%
#'   mutate(time = 8) %>%
#'   chemdose_chlordecay_chain(
#'     input_water = "balanced_water", cl2_dose = 6, treatment = "coag",
#'     cl_type = "chloramine"
#'   )
#'
#' \donttest{
#' # Initialize parallel processing
#' plan(multisession, workers = 2) # Remove the workers argument to use all available compute
#' example_df <- water_df %>%
#'   mutate(br = 50) %>%
#'   define_water_chain() %>%
#'   balance_ions_chain() %>%
#'   chemdose_chlordecay_chain(input_water = "balanced_water", cl2_dose = 4, time = 8)
#'
#' # Optional: explicitly close multisession processing
#' plan(sequential)
#' }
#'
#' @import dplyr
#' @export
#'
#' @returns A data frame containing a water class column with updated chlorine residuals.

chemdose_chlordecay_chain <- function(df, input_water = "defined_water", output_water = "disinfected_water",
                                      cl2_dose = 0, time = 0, treatment = "raw", cl_type = "chlorine") {
  ID <- NULL # Quiet RCMD check global variable note

  arguments <- construct_helper(
    df, list("cl2_dose" = cl2_dose, "time" = time),
    list("treatment" = treatment, "cl_type" = cl_type)
  )

  output <- df %>%
    subset(select = !names(df) %in% c("cl2_dose", "time", "treatment", "cl_type")) %>%
    mutate(
      ID = row_number()
    ) %>%
    left_join(arguments, by = "ID") %>%
    select(-ID) %>%
    mutate(!!output_water := furrr::future_pmap(
      list(
        water = !!as.name(input_water),
        cl2_dose = cl2_dose,
        time = time,
        treatment = treatment,
        cl_type = cl_type
      ),
      chemdose_chlordecay
    ))
}
