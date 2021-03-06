% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/volta.R
\name{plot_peak2_data}
\alias{plot_peak2_data}
\title{Plot Peak 2 CV vs. PMT voltages & find minimum "knee" voltage}
\usage{
plot_peak2_data(
  x,
  report_dir,
  image_dir = report_dir,
  trans_fun = log10,
  replace_with_na_condition = ~is.infinite(.x) | is.nan(.x),
  remove_empty_cols = TRUE,
  save_png = FALSE,
  png... = list(),
  plot_series... = list(),
  create_smooth_variables... = list(),
  plot_derivative = FALSE,
  points... = list(),
  debug = FALSE,
  xlsx_expression = NULL
)
}
\arguments{
\item{x}{List object generated by \link{get_peak2_data}}

\item{report_dir}{[Coming soon]}

\item{image_dir}{...}

\item{trans_fun}{...}

\item{replace_with_na_condition}{...}

\item{remove_empty_cols}{...}

\item{save_png}{...}

\item{plot_series...}{...}

\item{create_smooth_variables...}{...}

\item{plot_derivative}{...}

\item{points...}{...}

\item{debug}{...}

\item{xlsx_expression}{...}
}
\value{

}
\description{
Creates a plot of Peak 2 CV vs. PMT voltages for each channel for each instrument & finds minimum "knee" voltage for optimal resolution sensitivity. Optionally generates a summary report as an Excel workbook.
}
\examples{
\dontrun{
## There will totally be code here soon.
}
}
