% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlts.R
\name{specar.ci}
\alias{specar.ci}
\title{Confidence interval for the ar-spectrum and the dominant period.}
\usage{
specar.ci(x, order, resamp = 500, nfreq = 100, echo = TRUE)
}
\arguments{
\item{x}{A time series without missing values.}

\item{order}{a scalar representing the order to be considered. If
\code{"aic"} the order is selected automatically using the AIC criterion.}

\item{resamp}{the number of resamples of the ar-coefficients from the
covariance matrix.}

\item{nfreq}{the number of points at which to save the value for the power
spectrum (and confidence envelope).}

\item{echo}{If \code{TRUE}, a counter for each nrun shows the progress.}
}
\value{
An object of class "specar.ci" is returned consisting of the
following components: \item{order}{the ar-order.} \item{spectrum$freq}{the
spectral frequencies.} \item{spectrum$spec}{the estimated power-spectrum of
the data.} \item{resamp$spectrum}{gives the quantile summary for the
resampling distribution of the spectral powers.} \item{resamp$maxfreq}{the
full vector of output for the resampled max.frequencies.}
}
\description{
A function to estimate a "confidence interval" for the power spectrum and in
particular a confidence interval for the dominant period. The function uses
resampling of the autoregressive parameters to attain the estimate.
}
\details{
A "confidence interval" for the periodogram is obtained by resampling the
ar-coefficients using the variance-covariance matrix from the ar.mle object.

If a zero'th order process is chosen by using the AIC criterion, a first
order process will be used.

If the dynamics is highly nonlinear, the parametric estimate of the power
spectrum may be inappropriate.
}
\examples{

   data(plodia)


    fit <- specar.ci(sqrt(plodia), order=3, resamp=10) 

    \dontrun{plot(fit, period=FALSE)}

    summary(fit)
}
\seealso{
\code{\link{plot.specar.ci}} \code{\link{summary.specar.ci}}
}
\keyword{ts}
