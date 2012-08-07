\name{ddp_quantiles}
\alias{ddp_quantiles}
\title{Produce quantile summaries of model posterior samples}
\usage{
  ddp_quantiles(model.output, dosemat, Nfixed, Nrandom,
    Nsubject, typet)
}
\arguments{
  \item{model.output}{A list vector of objects returned by
  MCMC sampling functions.  e.g. mmCplusDpPost for
  \code{option = "mmcar"}.}

  \item{dosemat}{An \code{P x (T+1)} \code{matrix} object
  that maps \code{subjects} to treatment dosages.  The
  first column should be an intercept column (filled with
  1's).}

  \item{Nfixed}{Number of total fixed effects, both
  time-based and nuisance.}

  \item{Nrandom}{Number of total random effects, both
  time-based and nuisance, all grouped by subject.}

  \item{Nsubject}{Number of unique subjects (on which
  repeated measures are observed).}

  \item{typet}{A numeric vector of length equal to the
  number of treatments that contains the base distribution
  for each treatment.  \code{1 = "car"}, \code{2 = "mvn"},
  \code{3 = "ind"}.}
}
\value{
  A list object containing quantile summaries for all
  sampled model parameters.  \item{deviance.summary}{vector
  of length 3 summarizing quantiles for model deviance.}
  \item{beta.summary}{\code{Nfixed x 3} quantile summaries
  of model fixed effects.} \item{alpha.summary}{quantile
  summary of model global intercept parameter.}
  \item{bmat.summary}{list object of length \code{Nrandom},
  each cell containing a \code{P x 3} matrix of by-subject
  parameter quantile summaries.}
  \item{tauu.summary}{quantile summary for prior precision
  parameter employed for multiple membership random
  effects. An \code{nty x 3} matrix in the case of
  \code{nty} multiple membership effect terms.}
  \item{taue.summary}{quantile summary for model error
  precision parameter.} \item{taub}{\code{Nrandom x 3}
  quantile summaries for subject effect precision
  parameters.} \item{u.summary}{\code{S x 3} quantile
  summaries for multiple membership random effect
  parameters. A list of such matrices in the case of
  \code{nty} multiple membership effect terms.}
  \item{mm.summary}{\code{P.aff x 3} quantile summaries
  derived from multiplying the affected subject weight
  matrix by the multiple membership random effects.}
  \item{Dbar}{Model fit statistics.} \item{pD}{Model fit
  statistics.} \item{pV}{Model fit statistics.}
  \item{DIC}{Model fit statistics.} \item{lpml}{Model fit
  statistics.}
}
\description{
  Inputs MCMC samples for model parameters and constructs
  c(2.5\%,50\%,97.5\%) quantile summaries.
}
\author{
  Terrance Savitsky \email{tds151@gmail.com}
}
\seealso{
  \code{\link{dpgrowmm}}, \code{\link{dpgrow}}
}

