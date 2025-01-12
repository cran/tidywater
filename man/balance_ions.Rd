% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/balance_ions.R
\name{balance_ions}
\alias{balance_ions}
\title{Add Na, K, Cl, or SO4 to balance overall charge in a water}
\usage{
balance_ions(water)
}
\arguments{
\item{water}{Water created with define_water, which may have some ions set to 0 when unknown}
}
\value{
A water class object with updated ions to balance water charge.
}
\description{
This function takes a water defined by \code{\link{define_water}} and balances charge.
}
\details{
If more cations are needed, sodium will be added, unless a number for sodium is already provided and potassium is 0, then it will add potassium. Similarly,
anions are added using chloride, unless sulfate is 0. If calcium and magnesium are not specified when defining a water with
\code{\link{define_water}}, they will default to 0 and not be changed by this function.  This function is purely mathematical.
User should always check the outputs to make sure values are reasonable for the input source water.
}
\examples{
water_defined <- define_water(7, 20, 50, 100, 80, 10, 10, 10, 10, tot_po4 = 1) \%>\%
  balance_ions()

}
