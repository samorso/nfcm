% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihood.R
\name{splines.control}
\alias{splines.control}
\title{Control for spline}
\usage{
splines.control(
  df = NULL,
  knots = seq(0.2, 0.8, by = 0.2),
  degree = 3L,
  intercept = FALSE,
  Boundary.knots = c(0, 1),
  periodic = FALSE,
  derivs = 0L,
  integral = FALSE
)
}
\arguments{
\item{df}{degree of freedom (can be \code{NULL});}

\item{knots}{intenal breakpoints that define the splines;}

\item{degree}{the degree of the piecewise polynomial;}

\item{intercept}{boolean;}

\item{Boundary.knots}{boundary points at which to anchor the splines;}

\item{periodic}{boolean;}

\item{derivs}{the order of the derivative the spline;}

\item{integral}{boolean.}
}
\description{
Control for the spline basis, currently from
the \code{\link{splines2}} package.
}
\seealso{
\code{\link[splines2]{bSpline}}, \code{\link[splines2]{mSpline}},
\code{\link[splines2]{cSpline}}, \code{\link[splines2]{iSpline}}
}
