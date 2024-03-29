% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SSR.matrix.R
\name{SSR.recursive}
\alias{SSR.recursive}
\title{Calculate SSR recursively}
\usage{
SSR.recursive(y, x, beg, end, width = 2)
}
\arguments{
\item{y}{A dependent variable.}

\item{x}{Explanatory variables.}

\item{beg, end}{The start and the end of SSR calculating period.}

\item{width}{Minimum spacing between the breaks.}
}
\value{
The vector of calculated recursive SSR.
}
\description{
Calculate SSR recursively
}
\references{
Brown, R. L., J. Durbin, and J. M. Evans.
“Techniques for Testing the Constancy of Regression Relationships over Time.”
Journal of the Royal Statistical Society.
Series B (Methodological) 37, no. 2 (1975): 149–92.
}
\keyword{internal}
