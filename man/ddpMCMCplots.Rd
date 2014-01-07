\name{ddpMCMCplots}
\alias{ddpMCMCplots}
\alias{ddpplots}
\title{generate plots of posterior samples under \code{ddpgrow} model}
\usage{
ddpMCMCplots(subjecti.u, labt = NULL, typet, numt, theta.summary, lambda.mean,
  pmvn.mean = NULL, taucar.summary = NULL, alphacar.summary = NULL,
  Taucar = NULL, Alphacar = NULL, tauind.summary = NULL, Tauind = NULL,
  M, Taue, Deviance)
}
\arguments{
  \item{subjecti.u}{A vector of length \code{P}, number of
  unique subjects, containing unique set of user input
  values for \code{subject}.}

  \item{labt}{An vector object (of the same length as
  \code{typet}) providing user names for each treatment.
  The names are used in plot objects.}

  \item{typet}{A numeric vector of length equal to the
  number of treatments that specifies prior option for each
  treatment.  Options must be one of: 1 = "car", 2 = "mvn",
  3 = "ind".}

  \item{numt}{A numeric vector with same length as
  \code{typet} where each entry counts the number of doses
  for that treatment.}

  \item{theta.summary}{A list object of \code{q} elements,
  each containing an \code{P x 3} matrix of
  c(2.5\%,50\%,97.5\%) quantile summaries for each subject
  of the applicable subject random effect parameter.
  \code{P} = number of subjects, \code{q} = number of
  random effect parameters, per subject.}

  \item{lambda.mean}{A \code{q x q} matrix of mean values
  of the polyomial order covariance matrix, \code{Lambda},
  returned from \code{ddpgrow}.}

  \item{pmvn.mean}{A list object of length equal to the
  number of treatments with \code{"mvn" \%in\% typetreat}.
  Each list object contains an \code{numt[m] x numt[m]}
  matrix of mean elements of the "mvn" treatments
  covariance matrices, \code{Pmvn}, returned from
  \code{ddpgrow}.}

  \item{taucar.summary}{A \code{numcar x 3} numeric matrix
  of c(2.5\%,50\%,97.5\%) quantile summaries for the scale
  parameter for \code{numcar} treatments where \code{"car"
  \%in\% typetreat}.}

  \item{alphacar.summary}{A \code{numcar x 3} numeric
  matrix of c(2.5\%,50\%,97.5\%) quantile summaries for the
  strength parameter for \code{numcar} treatments where
  \code{"car" \%in\% typetreat}.}

  \item{Taucar}{An \code{iter.keep x numcar} matrix of
  posterior samples capturing the CAR precision parameter
  for each treatment where \code{"car" \%in\% typetreat}.}

  \item{Alphacar}{An \code{iter.keep x numcar} matrix of
  posterior samples capturing the CAR strength parameter
  for each treatment where \code{"car" \%in\% typetreat}.}

  \item{tauind.summary}{A list object of length equal to
  \code{numcar}, the number of treatment with \code{"in"
  \%in\% typetreat}.  Each list element contains an
  \code{numt[m] x 3} matrix of c(2.5\%,50\%,97.5\%)
  quantile summaries for the dosage scale parameters
  associated to that treatment \code{m}.}

  \item{Tauind}{A list object of length \code{numind}, the
  number treatments where \code{"car" \%in\% typetreat}
  with each element holding an \code{nkeep x numt[m]}
  matrix of sampled by-dose precision parameters for
  treatment \code{m}.}

  \item{M}{The \code{iter.keep x 1} matrix of posterior
  samples for the parameter capturing the number of
  clusters formed under the DP prior on the client
  effects.}

  \item{Taue}{\code{iter.keep x 1} matrix of posterior
  samples capturing the precision parameter for the model
  error term.}

  \item{Deviance}{\code{iter.keep x 1} matrix of posterior
  samples for the model deviance.}
}
\value{
A list of plot objects of class \code{ggplot2} including:
\item{p.theta}{stacked plots of b0,...,b(q-1) - vertical
lines for each client span 2.5% - 97.5% values with mean
noted.} \item{p.M}{MCMC trace plot of M, number of
clusters.} \item{p.acar}{MCMC trace plots of alphacar, the
CAR strength parameter. Plot is faceted for more than one
set of treatments under \code{"car"}.} \item{p.tcar}{MCMC
trace plots of alphacar, the CAR precision parameter. Plot
is faceted for more than one set of treatments under
\code{"car"}.} \item{p.taue}{MCMC trace plots of tau.e.}
\item{p.dev}{MCMC trace plots of deviance.}
\item{p.lam}{Heatmap (tiled) plot of posterior mean
covariance, Lambda, of polynomial orders of random effects,
Delta_i.} \item{p.mvn}{Heatmap plots of posterior mean
covariance, Pmvn, of treatment dosages under \code{"mvn"}
base distribution.  Plots are faceted by treatment.}
\item{p.aband}{95\% credible bands for alphacar CAR
strength parameters associated to treatments under
\code{"car"} base distribution.} \item{p.tband}{95\%
credible bands for taucar CAR precision parameters
associated to treatments under \code{"car"} base
distribution.} \item{p.iband}{95\% credible bands for
tauind dose precision parameters for \code{"ind"} base
distribution, faceted by treatment (if more than one).}
}
\description{
Constructs plots for subject effects, theta, with credible
intervals as well as selected trace plots. Returns a list
of objects of class \code{ggplot}.
}
\note{
Intended as an internal function for \code{\link{ddpgrow}}
}
\author{
Terrance Savitsky \email{tds151@gmail.com}
}
\seealso{
\code{\link{ddpgrow}}, \code{\link{dpgrowmm}},
\code{\link{dpgrow}}, \code{\link{dpgrowmult}}
}

