\name{ddpPost}
\alias{ddpPost}
\title{Run a Bayesian mixed effects model for by-subject random effects with DDP prior}
\usage{
ddpPost(y, X, Z, subject, dosemat, numt, typet, Omega, omegaplus, n.iter,
  n.burn, n.thin, shapealph, ratebeta, M.init)
}
\arguments{
  \item{y}{An \emph{N x 1} response (of subject-measure
  cases)}

  \item{X}{Fixed effects design matrix}

  \item{Z}{Random effects design matrix.  Assumed grouped
  by \code{subjects}}

  \item{subject}{An \emph{N x 1} set of subject
  identifiers}

  \item{dosemat}{An \emph{P x T} Anova or Multiple
  Membership design matrix linking treatment dosages to
  subjects where \emph{T} is the total number dosages
  across all treatments + 1 for an intercept.  This
  formulation assumes there is a hold-out dose for each
  treatment.  e.g. the null dosage.}

  \item{numt}{A numeric vector of length equal to the
  number of treatments that contains the number of dosages
  for each treatment.}

  \item{typet}{A numeric vector of length equal to the
  number of treatments that contains the base distribution
  for each treatment.  \code{1 = "car"}, \code{2 = "mvn"},
  \code{3 = "ind"}}

  \item{Omega}{A list object of length equal to the number
  of treatments with \code{"car"} selected for base
  distribution.  Each entry is an \code{numt[m] x numt[m]}
  numeric CAR adjacency matrix for the dosages of treatment
  \code{m}.}

  \item{omegaplus}{A list object of length equal to the
  number of treatments under \code{"car"} containing
  numeric vectors that are rowSums of the corresponding
  matrix element in \code{Omega}.}

  \item{n.iter}{The number of MCMC iterations}

  \item{n.burn}{The number of MCMC burn-in iterations to
  discard}

  \item{n.thin}{The step increment of MCMC samples to
  return}

  \item{shapealph}{The shape parameter for the \eqn{\Gamma}
  prior on the DP concentration parameter.}

  \item{ratebeta}{The rate parameter for the \eqn{\Gamma}
  prior on the DP concentration parameter.}

  \item{M.init}{Initial MCMC chain scalar value for number
  of by-subject clusters.  If excluded defaults to
  \code{length(unique(subject))}.}
}
\value{
res A list object containing MCMC runs for all model
parameters.
}
\description{
An internal function to \code{\link{ddpgrow}}
}
\note{
Intended as an internal function for \code{\link{ddpgrow}}
}
\author{
Terrance Savitsky \email{tds151@gmail.com}
}
\seealso{
\code{\link{dpgrow}}
}

