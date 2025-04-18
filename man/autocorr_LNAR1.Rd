% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux_LNAR1.r
\name{autocorr_LNAR1}
\alias{autocorr_LNAR1}
\title{Autcorrelation function of a LNAR1 process}
\usage{
autocorr_LNAR1(k, lmean, lsd, phi)
}
\arguments{
\item{k}{the lag}

\item{lmean}{log mean}

\item{lsd}{log standard deviation}

\item{phi}{correlation parameter of the latent AR(1)}
}
\value{
the autocorrelation at \code{k}
}
\description{
Autcorrelation function of a LNAR1 process
}
