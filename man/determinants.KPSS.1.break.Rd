% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils.Variables.R
\name{determinants.KPSS.1.break}
\alias{determinants.KPSS.1.break}
\title{Construct determinant variables for \link{KPSS.1.break}}
\usage{
determinants.KPSS.1.break(model, n.obs, break.point)
}
\arguments{
\item{model}{A scalar equal to
\itemize{
\item 1: Model with trend, break in const,
\item 2: Model with const and trend, break in const,
\item 3: Model with const and trend, break in trend,
\item 4: Model with const and trend, break in const and trend.
}}

\item{n.obs}{Number of observations.}

\item{break.point}{Break point.}
}
\value{
Matrix of determinant variables.
}
\description{
Construct determinant variables for \link{KPSS.1.break}
}
\details{
Procedure to compute deterministic terms
for KPSS with 1 structural break.
}
\keyword{internal}
