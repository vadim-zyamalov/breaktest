% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils.Estimates.R
\name{NW.loocv}
\alias{NW.loocv}
\title{LOO-CV for h in Nadaraya–Watson kernel regression.}
\usage{
NW.loocv(y, x, kernel = "unif")
}
\arguments{
\item{y}{A dependent variable.}

\item{x}{An explanatory variable.}

\item{kernel}{Needed kernel, currently only \code{unif} and \code{gauss}:
\itemize{
\item \code{unif}: \eqn{K(x) = \left\{\begin{array}{ll}
1 & \frac{|x - x_i|}{h} \leq 1 \\
0 & \textrm{otherwize}
\end{array}\right.}
\item \code{gauss}: \eqn{\Phi(\frac{x - x_i}{h})}
}}
}
\value{
A list of arguments as well as the estimated bandwidth \code{h}.
}
\description{
LOO-CV for h in Nadaraya–Watson kernel regression.
}
\references{
Harvey, David I., S. Leybourne, Stephen J., and Yang Zu.
“Nonparametric Estimation of the Variance Function
in an Explosive Autoregression Model.”
School of Economics. University of Nottingham, 2022.
}
\keyword{internal}
