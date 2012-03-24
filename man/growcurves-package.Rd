\docType{package}
\name{growcurves-package}
\alias{growcurves}
\alias{growcurves-package}
\alias{package-growcurves}
\title{Bayesian semiparametric growth curve models with employment of multiple membership random effects.}
\description{
  \tabular{ll}{ Package: \tab growcurves\cr Type: \tab
  Package\cr Version: \tab 0.1.1\cr Date: \tab
  2012-04-11\cr License: \tab GPL (>= 2) \cr LazyLoad: \tab
  yes\cr }
}
\details{
  Employs a Dirichlet Process (DP) prior on the set of
  by-subject random effect parameters under repeated waves
  of measurements to allow the number of random effect
  parameters specified per subject, \code{q}, to be equal
  to the number of measurement waves, \code{T}.  Random
  effects are grouped by subject and all \code{q}
  parameters receive the DP prior.  Additional sets of
  random effects under different (than subject) grouping
  factors are accommodated through multiple membership
  constructions, each of which maps a set of effects under
  a given grouping back to affected subjects.

  CORE FUNCTIONS

  \code{\link{dpgrow}} performs Bayesian mixed effects
  estimation for data characterized by repeated subject
  measures (typically in time).  The user inputs a
  \code{subject} identifier vector, a vector of \code{time}
  measurements, and a \code{trt} vector for treatment/group
  assignments.  Fixed and random effects are then
  automatically generated and both subject and treatment
  level growth curves are constructed.

  \code{\link{dpgrowmm}} is very similar to \code{dpgrow},
  but it adds an additional set of random effects which
  aren't grouped by subject that may be used to inject
  treatment dosage or attendance patterns that is mapped
  back to clients via a multiple membership weight matrix.

  \code{\link{dpgrowmult}} is very similar to
  \code{dpgrowmm}, but it allows more than one set of
  multiple membership effects. Each multiple membership
  effects term (block) may apply to any sub-set of subjects
  through specification of the weight matrix and
  identification of affected subjects for that term.

  \code{\link{growplot}} uses model outputs from
  \code{dpgrow}, \code{dpgrowmm} and \code{dpgrowmult} to
  provide by-subject growth curves in two forms: 1. Growth
  curves aggregated under specified groupings; 2. Indvidual
  growth curves plotted along with data for selected (or
  random subset of) subjects.

  \code{\link{trtplot}} uses model outputs from
  \code{dpgrow}, \code{dpgrowmm} and \code{dpgrowmult} to
  plot a distribution for fixed mean effects difference
  between two selected treatments across a range of chosen
  models for a one or more chosen time points. Outputs
  include a set of boxplots for each time point that span
  95% credible interval.

  \code{\link{effectsplot}} uses model outputs from
  \code{dpgrow}, \code{dpgrowmm} and \code{dpgrowmult} to
  overlay plots of multiple membership effects for a given
  term under use of different prior formulations and/or
  from distinctly formulated models (e.g. with varying
  numbers of multiple membership terms).

  SIMULATED DATA SETS

  There are 3 simulated data sets available in order to
  allow exploration of the engine and associated accessor
  functions.

  \code{\link{datsim}} Simulated dataset with two treatment
  arms (treatment and control) composed from a model with a
  Dirichlet process (DP) prior on the set of client effects
  and a single MM term under a \code{"mmcar"} formulation.
  Structured to express similar properties as the case
  example in both Savitsky and Paddock (2012) references,
  below.

  \code{\link{datsimcov}} Of similar structre to
  \code{simdat}, only the data generating model now
  additionally employs 2 nuisance fixed effects.

  \code{\link{datsimmult}} Simulated data under 4 treatment
  arms generated from a model with now 4 multiple
  membership terms.  The terms are generated under
  \code{c("mmi","mmcar","mmi","mmcar")} prior formulations.

  BENCHMARK DATA SETS \code{\link{datbrghtterms}} Data
  derived from BRIGHT study reviewed in reference and
  includes BDI-II depressive symptom scores for client
  experimental units.  Associated data objects are included
  to facilitate runs under engine functions \code{dpgrow},
  \code{dpgrowmm} and \code{dpgrowmult}.
}
\examples{
\dontrun{
## extract simulated dataset
library(growcurves)
data(datsim)
## attach(datsim)
## run dpgrow mixed effects model, returning object of class "dpgrow"
shape.dp			= 4
res				= dpgrow(y = datsim$y, subject = datsim$subject, trt = datsim$trt, time = datsim$time,
						n.random = 3, n.fix_degree = 2, n.iter = 10000, n.burn = 2000,
						n.thin = 10, shape.dp = shape.dp, option = "dp")
## Each plot is a "ggplot2" object saved in a list to plot.results
plot.results			= plot(res) ## includes subject and treatment growth curves
## Extract credible interval summary (2.5\%, mean, 97.5\%) values for each sampled parameter set,
## including deviance.  Includes fit statistics, such as Dbar, DIC, pD, lpml.  Note: DIC is
## the DIC3 of Celeaux et. al. (2006) for option = "dp".  Finally, the constructed fixed and random effects
## matrices, X and Z, are returned with growth curve covariates appended to user submitted nuisance covariates.
summary.results		= summary(res)
## View the summary results in the console
print(summary.results)
## Collect posterior sampled values over the (n.iter - n.burn) retained iterations for each
## sampled parameter.
samples.posterior		= samples(res)
residuals			= resid(res) ## model residuals (y - fit)
## Model with DP on clients effects, but now INCLUDE session random effects
## in a multiple membership construction communicated with the N x S matrix, W.subj.aff.
## Returns object, res.mm, of class "dpgrowmm".
shape.dp			= 3
strength.mm			= 0.001
res.mm			= dpgrowmm(y = datsim$y, subject = datsim$subject, trt = datsim$trt, time = datsim$time, n.random = 3,
						Omega = datsim$Omega, group = datsim$group, subj.aff = datsim$subj.aff,
						W.subj.aff = datsim$W.subj.aff, n.iter = 10000, n.burn = 2000, n.thin = 10,
						strength.mm = strength.mm, shape.dp = shape.dp,
						option = "mmcar")
plot.results			= plot(res.mm)
}
}
\author{
  Terrance Savitsky \email{tds151@gmail.com} Susan Paddock
  \email{paddock@rand.org}
}
\references{
  S. M. Paddock and T. D. Savitsky (2012) Bayesian
  Hierarchical Semiparametric Modeling of Longitudinal
  Post-treatment Outcomes from Open-enrollment Therapy
  Groups, Invited re-submission to: JRSS Series A
  (Statistics in Society).

  T. D. Savitsky and S. M. Paddock (2012) Visual Sufficient
  Statistics for Repeated Measures data with growcurves for
  R, submitted to: Journal of Statistical Software.
}
\keyword{package}

