% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NW.loocv.R
\name{NW.loocv}
\alias{NW.loocv}
\title{NW.loocv - LOO-CV for h in Nadaraya–Watson kernel regression.}
\usage{
NW.loocv(y, x, kernel = "unif")
}
\arguments{
\item{y}{LHS dependent variable.}

\item{x}{RHS explanation variable.}

\item{kernel}{Needed kernel, currently only `unif` and `gauss`.}
}
\value{
A list of arguments as well as the estimated bandwidth `h`.
}
\description{
NW.loocv - LOO-CV for h in Nadaraya–Watson kernel regression.
}
\details{
See Harvey, Leybourne and Zu (2021) for more information.
}