% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/solvedose_ph_alk.R
\name{solvedose_ph_once}
\alias{solvedose_ph_once}
\title{Apply `solvedose_ph` to a dataframe and create a new column with numeric dose}
\usage{
solvedose_ph_once(
  df,
  input_water = "defined_water",
  output_column = "dose_required",
  target_ph = NULL,
  chemical = NULL
)
}
\arguments{
\item{df}{a data frame containing a water class column, which has already been computed using
\code{\link{define_water_chain}}. The df may include a column with names for each of the chemicals being dosed.}

\item{input_water}{name of the column of water class data to be used as the input. Default is "defined_water".}

\item{output_column}{name of the output column storing doses in mg/L. Default is "dose_required".}

\item{target_ph}{set a goal for pH using the function argument or a data frame column}

\item{chemical}{select the chemical to be used to reach the desired pH using function argument or data frame column}
}
\value{
A data frame containing the original data frame and columns for target pH, chemical dosed, and required chemical dose.
}
\description{
This function allows \code{\link{solvedose_ph}} to be added to a piped data frame.
Its output is a chemical dose in mg/L.
}
\details{
The data input comes from a `water` class column, initialized in \code{\link{define_water}} or \code{\link{balance_ions}}.

If the input data frame has column(s) named "target_ph" or "chemical", the function will use the column(s)
as function argument(s). If these columns aren't present, specify "target_ph" or "chemical" as function arguments.
The chemical names must match the chemical names as displayed in \code{\link{solvedose_ph}}.
To see which chemicals can be dosed, see \code{\link{solvedose_ph}}.

 For large datasets, using `fn_once` or `fn_chain` may take many minutes to run. These types of functions use the furrr package
 for the option to use parallel processing and speed things up. To initialize parallel processing, use
 `plan(multisession)` or `plan(multicore)` (depending on your operating system) prior to your piped code with the
 `fn_once` or `fn_chain` functions. Note, parallel processing is best used when your code block takes more than a minute to run,
 shorter run times will not benefit from parallel processing.
}
\examples{

library(purrr)
library(furrr)
library(tidyr)
library(dplyr)

example_df <- water_df \%>\%
  define_water_chain() \%>\%
  mutate(
    target_ph = 10,
    chemical = rep(c("naoh", "mgoh2"), 6)
  ) \%>\%
  solvedose_ph_once(input_water = "defined_water")

example_df <- water_df \%>\%
  define_water_chain() \%>\%
  solvedose_ph_once(input_water = "defined_water", target_ph = 8.8, chemical = "naoh")


example_df <- water_df \%>\%
  define_water_chain() \%>\%
  mutate(target_ph = seq(9, 10.1, .1)) \%>\%
  solvedose_ph_once(chemical = "naoh")

\donttest{
# Initialize parallel processing
plan(multisession, workers = 2) # Remove the workers argument to use all available compute
example_df <- water_df \%>\%
  define_water_chain() \%>\%
  mutate(target_ph = seq(9, 10.1, .1)) \%>\%
  solvedose_ph_once(chemical = "naoh")

# Optional: explicitly close multisession processing
plan(sequential)
}

}
\seealso{
\code{\link{solvedose_ph}}
}
