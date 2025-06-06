% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ffpo_2d.R
\name{ffpo_2d}
\alias{ffpo_2d}
\title{Defining partially observed bidimensional functional data terms in VDPO formulae}
\usage{
ffpo_2d(X, miss_points, missing_points, nbasis = rep(15, 4), bdeg = rep(3, 4))
}
\arguments{
\item{X}{partially observed bidimensional functional covariate \code{matrix}.}

\item{miss_points, missing_points}{\code{list} of missing observation points.
See 'Details' for more information about the difference in structure between both.}

\item{nbasis}{number of basis to be used.}

\item{bdeg}{degree of the basis to be used.}
}
\value{
The function is interpreted in the formula of a \code{VDPO} model.
\code{list} containing the following elements:
\itemize{
\item \code{B_ffpo2d} design matrix.
\item \code{Phi_ffpo2d} bidimensional B-spline basis used for the functional coefficient.
\item \code{M_ffpo2d} the \code{missing_points} used as input in the function.
\item \code{nbasis} number of the basis used.
}
}
\description{
Auxiliary function used to define \code{ffpo_2d} terms within \code{VDPO} model
formulae.
}
\details{
The difference between
miss_points and missing_points is the format in which the data is
presented.

\code{miss_points} is a \code{list} of \code{list}s where each inner list corresponds
to the observation points in the y-axis and contains the observation points
of the missing values for the x-axis. \code{miss_points} acts as a guide for
identifying and addressing missing observations in functional data and is used
for properly calculating the inner product matrix.

\code{missing_points} is a \code{list} where each element is a \code{matrix}
containing the missing observations points.
}
