% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/criteria.R
\name{bestseparated}
\alias{bestseparated}
\title{Criteria for measuring the quality of clustering solution.}
\usage{
bestseparated(mean, sigma, groups, ...)
}
\arguments{
\item{mean}{A matrix, containing the mean vectors for the clusters.}

\item{sigma}{An array, containing the covariance matrices for the clusters.}

\item{groups}{The number of clusters.}

\item{...}{Not used.}
}
\value{
A number.
}
\description{
Summarise a model-based clustering solution by a single number.
}

