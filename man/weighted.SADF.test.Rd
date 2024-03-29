% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weighted.SADF.test.R
\name{weighted.SADF.test}
\alias{weighted.SADF.test}
\alias{weighted.GSADF.test}
\title{Weighted supremum ADF test}
\usage{
weighted.SADF.test(
  y,
  trim = 0.01 + 1.8/sqrt(length(y)),
  const = TRUE,
  alpha = 0.05,
  iter = 4 * 200,
  urs = TRUE,
  seed = round(10^4 * sd(y))
)

weighted.GSADF.test(
  y,
  trim = 0.01 + 1.8/sqrt(length(y)),
  const = TRUE,
  alpha = 0.05,
  iter = 4 * 200,
  urs = TRUE,
  seed = round(10^4 * sd(y))
)
}
\arguments{
\item{y}{A time series of interest.}

\item{trim}{The trimming parameter to find the lower and upper bounds of
possible break dates.}

\item{const}{Whether the constant needs to be included.}

\item{alpha}{A significance level of interest.}

\item{iter}{Nnumber of iterations.}

\item{urs}{Use \verb{union of rejections} strategy.}

\item{seed}{A seed parameter for the random number generator.}
}
\value{
An object of type \code{sadf}. It's a list of:
\itemize{
\item \code{y},
\item \code{trim},
\item \code{const},
\item \code{alpha},
\item \code{iter},
\item \code{urs},
\item \code{seed},
\item \code{sigma.sq}: the estimated variance,
\item \code{BZ.values}: a series of BZ-statistic,
\item \code{supBZ.value}: the maximum of \code{supBZ.values},
\item \code{supBZ.bootstsrap.values}: bootstrapped supremum BZ values,
\item \code{supBZ.cr.value}: supremum BZ \eqn{\alpha} critical value,
\item \code{p.value},
\item \code{is.explosive}: 1 if \code{supBZ.value} is greater than \code{supBZ.cr.value}.
}

if \code{urs} is \code{TRUE} the following items are also included:
\itemize{
\item vector of \eqn{t}-values,
\item the value of the SADF test statistic,
\item \code{SADF.bootstrap.values}: bootstrapped SADF values,
\item \code{U.value}: union test statistic value,
\item \code{U.bootstrap.values}: bootstrapped series of \code{U.value},
\item \code{U.cr.value}: critical value of \code{U.value}.
}
}
\description{
Weighted supremum ADF test
}
\details{
Refactored original code by Kurozumi et al.
}
\references{
Harvey, David I., Stephen J. Leybourne, and Yang Zu.
“Testing Explosive Bubbles with Time-Varying Volatility.”
Econometric Reviews 38, no. 10 (November 26, 2019): 1131–51.
https://doi.org/10.1080/07474938.2018.1536099.

Kurozumi, Eiji, Anton Skrobotov, and Alexey Tsarev.
“Time-Transformed Test for Bubbles under Non-Stationary Volatility.”
Journal of Financial Econometrics, April 23, 2022.
https://doi.org/10.1093/jjfinec/nbac004.
}
