% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dpgrowmm.R
\name{mmIplusDpPost}
\alias{mmI}
\alias{mmIplusDpPost}
\title{Bayesian mixed effects model with a DP prior on by-subject effects and zero mean independent Gaussian priors on multiple membership effects}
\usage{
mmIplusDpPost(y, X, Z, Wcase, Wsubject, subjects, niter, nburn, nthin,
  strength.mm, shapealph, ratebeta)
}
\arguments{
\item{y}{An \emph{N x 1} response (of subject-measure cases)}

\item{X}{Fixed effects design matrix}

\item{Z}{Random effects design matrix.  Assumed grouped by \code{subjects}}

\item{Wcase}{An \emph{N x 1} multiple membership weight matrix to map supplemental random effects}

\item{Wsubject}{An \emph{P.aff x S} multiple membership weight matrix with rows equal to number of unique affected subjects}

\item{subjects}{An \emph{N x 1} set of subject identifiers}

\item{niter}{The number of MCMC iterations}

\item{nburn}{The number of MCMC burn-in iterations to discard}

\item{nthin}{The step increment of MCMC samples to return}

\item{strength.mm}{The shape and rate parameters for the \eqn{\Gamma} prior on the CAR precision parameter, \eqn{\tau_{\gamma}}}

\item{shapealph}{The shape parameter for the \eqn{\Gamma} prior on the DP concentration parameter.
The rate parameter is set of \code{1}.}

\item{ratebeta}{The rate parameter for the \eqn{\Gamma} prior on the DP concentration parameter. Default value is \code{1}.}
}
\value{
res A list object containing MCMC runs for all model parameters.
}
\description{
An internal function to \code{\link{dpgrowmm}}
}
\note{
Intended as an internal function for \code{\link{dpgrowmm}}
}
\author{
Terrance Savitsky \email{tds151@gmail.com}
}
\seealso{
\code{\link{dpgrowmm}}
}

