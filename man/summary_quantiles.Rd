\name{summary_quantiles}
\alias{quantiles}
\alias{summary_quantiles}
\title{Produce quantile summaries of model posterior samples}
\usage{
  summary_quantiles(model.output, Nfixed, Nrandom,
    Nsubject, Nsubj.aff = NULL)
}
\arguments{
  \item{model.output}{A list vector of objects returned by
  MCMC sampling functions.  e.g. mmCplusDpPost for option =
  "mmcar".}

  \item{Nfixed}{Number of total fixed effects, both
  time-based and nuisance.}

  \item{Nrandom}{Number of total random effects, both
  time-based and nuisance, all grouped by subject.}

  \item{Nsubject}{Number of unique subjects (on which
  repeated measures are observed).}

  \item{Nsubj.aff}{Number of subjects, \code{P.aff},
  receiving multiple membership effects}
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

