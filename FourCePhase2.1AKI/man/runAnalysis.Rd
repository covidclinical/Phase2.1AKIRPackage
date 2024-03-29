% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runAnalysis.R
\name{runAnalysis}
\alias{runAnalysis}
\title{Runs the analytic workflow for the AKI project}
\usage{
runAnalysis(
  is_obfuscated = TRUE,
  factor_cutoff = 5,
  ckd_cutoff = 2.25,
  restrict_models = FALSE,
  docker = TRUE,
  input = "/4ceData/Input",
  siteid_nodocker = "",
  skip_qc = FALSE,
  offline = FALSE,
  allow_no_prior_cr = FALSE,
  use_rrt_surrogate = TRUE,
  print_rrt_surrogate = FALSE,
  modified_severity = FALSE,
  debug_on = FALSE,
  date_cutoff = "2020-09-10",
  lab_date_cutoff = "2021-09-10",
  custom_output = FALSE,
  custom_output_dir = "/4ceData/Output"
)
}
\description{
Runs the analytic workflow for the AKI project
}
\keyword{4CE}
