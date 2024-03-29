% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlts.R
\name{portman.Q}
\alias{portman.Q}
\title{Ljung-Box test for whiteness in a time series.}
\usage{
portman.Q(x, K)
}
\arguments{
\item{x}{A time series (vector without missing values).}

\item{K}{the maximum lag of the ACF to be used in the test.}
}
\value{
A vector is returned consisting of the asymptotic chi-square value,
the associated d.f. and asymptotic p.val for the test of whiteness.
}
\description{
portman.Q uses the cumulative ACF to test for whiteness of a time series.
}
\details{
This is the Ljung-Box version of the the Portmanteau test for whiteness
(Tong 1990). It may in particular be useful to test for whiteness in the
residuals from time series models.
}
\examples{

   data(plodia)

   portman.Q(sqrt(plodia), K = 10) 

   fit <- ar(sqrt(plodia)) 
   portman.Q(na.omit(fit$resid), K = 10) 
}
\references{
Tong, H. (1990) Non-linear time series : a dynamical system
approach. Clarendon Press, Oxford.
}
\keyword{ts}
