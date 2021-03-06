% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weibullfns.R
\name{weib_icdf}
\alias{weib_icdf}
\title{Inverse CDF for the Weibull distribution}
\usage{
weib_icdf(u, theta, gamma)
}
\arguments{
\item{u}{Numerical value(s) to be converted to Weibull variable(s)}

\item{theta}{Scale parameter \eqn{\theta}}

\item{gamma}{Shape parameter \eqn{\gamma}}
}
\value{
Output is a value or vector of values
from the inverse CDF of the Weibull distribution.
}
\description{
\code{weib_icdf} returns a value from the inverse CDF of the
Weibull distribution.
}
\details{
This function uses the Weibull density of the form
\deqn{f(t)=\theta t^(\gamma - 1)exp(-\theta/\gamma t^(\gamma))}
to get the inverse CDF
\deqn{F^(-1)(u)=(-\gamma/\theta log(1-u))^(1/\gamma).} It can be
implemented directly and is also called by the functions
\code{\link{weib_memsim}} and \code{\link{weib_cdfsim}}.
}
\examples{
simdta <- weib_icdf(u = runif(10), theta = 0.05, gamma = 2)

}
