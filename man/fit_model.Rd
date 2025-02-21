% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_model.R
\name{fit_model}
\alias{fit_model}
\title{fit_model}
\usage{
fit_model(
  X,
  Y,
  allowed_dependence,
  reg_method = "linear",
  tau = NULL,
  method = "quantile",
  penalty = NULL,
  path_tracking = NULL
)
}
\arguments{
\item{X}{matrix or data.frame of n observation of the decision variable of size d.}

\item{Y}{matrix or data.frame of n observation of the p objectifs.}

\item{reg_method}{only "linear" available yet.}

\item{tau}{vector of size NCOL(Y) containing the risk in ]0,1[ of each objectif.}

\item{method}{either "quantile" or "expected"}

\item{penalty}{NULL by default}

\item{path_tracking}{path where to write the step of the running function.}
}
\description{
fit_model
}
