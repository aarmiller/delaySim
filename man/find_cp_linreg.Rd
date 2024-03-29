% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cp_functions.R
\name{find_cp_linreg}
\alias{find_cp_linreg}
\title{Find the change point in count data using linear regression models}
\usage{
find_cp_linreg(
  data,
  var_name = "n_miss_visits",
  method = "lm",
  eval_criteria = "AIC",
  return_miss_only = FALSE,
  specify_cp = NULL,
  week_period = FALSE
)
}
\arguments{
\item{data}{A dataset of visit counts}

\item{var_name}{The name of the count variable to find the change-point for}

\item{method}{The method used to fit curves before and after the changepoint. Options include "lm",
"lm_quad", "lm_cube", "quad", "cube", "exp", "spline"}

\item{eval_criteria}{The evaluation criteria used to find change points}

\item{return_miss_only}{Logical argument to only return the tibbles of miss visit counts}

\item{specify_cp}{Set a specific change point you want to use instead of searching for optimal change point. Enter a postive integer value
repersenting the days before the index on which you you want to specify the change point. (e.g. 100 would be 100 days before the index)}

\item{week_period}{Logical to incorporate a "day of the week" effect into
the linear model. Note this is only sensible for one-day period aggregation.}
}
\description{
Find the change point in count data using linear regression models
}
\examples{
cp_result_original <- final_time_map \%>\%
prep_cp_data(event_name = "any_ssd", start_day = 1, by_days = 1) \%>\%
find_cp_linreg(var_name="n_miss_visits", method="lm_cube")

}
