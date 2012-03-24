\name{dpgrowmult}
\alias{dpgrowmult}
\alias{dpgrowmult.default}
\title{Bayesian semiparametric growth curve models under employment of more-than-\code{1} multiple membership random effects (block) term.}
\usage{
  dpgrowmult(y, subject, trt, time, n.random, n.fix_degree,
    formula, random.only, data, Omega, group, subj.aff,
    W.subj.aff, n.iter, n.burn, n.thin, strength.mm,
    shape.dp, plot.out, option, ulabs)
}
\arguments{
  \item{y}{A univariate continuous response, specified as
  an \emph{N x 1} matrix or vector, where \code{N} captures
  the number of subject-time cases (repeated subject
  measures).  Data may reflect unequal number of measures
  per subject.  Missing occasions are left out as no
  \code{NA} values are allowed.}

  \item{subject}{The objects on which repeated measures are
  conducted that serves as the random effects grouping
  factor.  Input as an \emph{N x 1} matrix or vector of
  subject-measure cases in either integer or character
  formt; e.g. \code{(1,1,1,2,2,3,3,3,...,n,n,n)}, where
  \code{n} is the total number of subjects.}

  \item{trt}{An integer or character vector of length
  \code{N} (number of cases) indicating treatment group
  assignments for each case.  May also be input as length
  \code{P} vector, where \code{P} is the number of unique
  subjects, indicating subject group assignment.  Multiple
  treatment groups are allowed and if the vector is entered
  as numeric, e.g. \code{(0,1,2,3,..)}, the lowest numbered
  group is taken as baseline (captured by global fixed
  effects).  If entered in character format, the first
  treatment entry is taken as baseline.  If the are no
  treatment (vs. control) groups, then this vector may be
  excluded (set to NULL).}

  \item{time}{A univariate vector of length \code{N},
  capturing the time points associated to each by-subject
  measure.  Mav leave blank if only one time point (no
  repeated measures).}

  \item{n.random}{The desired number of subject random
  effect terms, \code{q}.  Since a DP prior is used on
  client effects, may be set equal to the number of
  measurement waves, \code{T}.  The \code{y, trt, time}
  vectors will together be used to create both fixed and
  random effect design matrices.  The random effects matrix
  will be of the the form, \code{(1, time, ... ,
  time^(n.random - 1))} (grouped, by \code{subject}).  This
  formulation is a growth curve model that allows
  assessment of by-treatment effects and by-client growth
  curves.}

  \item{n.fix_degree}{The desired polynomial order in time
  to use for generating time-based fix effects.  The fixed
  effects matrix will be constructed as, \code{(time, ...,
  time^(n.fix_degree), trt_1,time*trt_1, ...
  ,time^(n.fix_degree)*trt_l, trt_L,...,
  time^(n.fix_degree)*trt_L)}.  This formulation is a
  growth curve model that allows assessment of by-treatment
  effects and by-client growth curves.  If
  \code{is.null(n.fix_degree) | n.fix_degree == 0 &
  is.null(trt)} time-by-treatment fixed effects and growth
  curves are not generated.}

  \item{formula}{Nuisance fixed and random effects may be
  entered in \code{formula} with the following format,
  \code{y ~ x_1 + x_2*x_3 | z_1*z_2 } as an object of class
  \code{formula}.  The bar, \code{|}, separates fixed and
  random effects.  If it is only desired to enter either
  fixed or random effects, but not both then the \code{|}
  may be omitted.  Note: the nuisance random effects are
  assumed to be grouped by subject.  The fixed and random
  effects values may change with each repeated measure;
  however, within subject growth curves will keep constant
  \code{z} and \code{x} values between measurement waves.
  It is possible to bypass the growth curve construction by
  leaving \code{y, trt, time, n.random, n.fix_degree} blank
  and entering only \code{formula}, instead.  The model
  output plots, will, however exclude growth curves in that
  event.  If a formula is input (which requires response,
  \code{y}) then the separate entry of \code{y} may be
  omitted.  If the parameter \code{y} is input, it will be
  over-written by that from \code{formula}.}

  \item{random.only}{A Boolean variable indicating whether
  the input formula contains random (for fixed) effects in
  the case that only one set are entered.  If excluded and
  \code{formula} is entered without a \code{|},
  \code{random.only} defaults to \code{FALSE}.}

  \item{data}{a \code{data.frame} containing the variables
  with names as specified in \code{formula}, including the
  response, \code{y}.}

  \item{Omega}{A list object of length equal to the number
  of multiple membership (MM) effect terms chosen with the
  \code{"mmcar"} option.  List element \code{i} contains an
  \emph{S[i] x S[i]} numeric matrix to encode the CAR
  adjacency matrix, where \code{S[i]} is the number of
  effects mapped to subjects for list component \code{i}.
  This input is required only under \code{option =
  "mmcar"}.}

  \item{group}{A list object of length equal to the number
  of MM terms chosen with prior formulation options
  \code{"mmcar"} or \code{"mmigrp"}.  List element \code{i}
  contains a numeric or character vector of length
  \code{S[i]}, providing group identifiers for each of
  \code{S[i]} effects in term \code{i}. (e.g.
  \code{(1,1,1,2,2,...)}.  If there is only a single group
  for term [i], this element should be loaded with an
  \code{S[i] x 1} vector of a single value.}

  \item{subj.aff}{A list object of length equal to the
  number of total MM terms.  List element \code{i} contains
  a \code{Paff[i] x 1} vector subset of \code{subject}
  composed with unique subject identifiers that are linked
  to the effects in term \code{i}; e.g. one or more
  treatment cohorts.  \code{P.aff[i]} is the length of the
  unique subjects linked to the effects in MM term
  \code{i}.  If all subjects are to receive the mapping of
  multiple membership effects then \code{Paff[i]} should
  contain a list of the unique subjects in \code{subject}.}

  \item{W.subj.aff}{A list object of length equal to the
  number of MM terms.  List element \code{i} contains a
  \emph{P.aff[i] x S[i]} numeric matrix that maps a set of
  random effects to affected subjects
  (\code{subj.aff[[i]]}).  It is assumed that the row order
  is the same as the order of \code{subj.aff[[i]]}.  If
  \code{W.subj.aff[[i]]} is a multiple membership weight
  matrix, then the rows will sum to 1, though this is not
  required.  The rows of \code{W.subj.aff} may
  alternatively be formulated with indicators for whether
  each of \code{S} treatment dosages are linked to a given
  \code{subject}.}

  \item{n.iter}{Total number of MCMC iterations.}

  \item{n.burn}{Number of MCMC iterations to discard.
  \code{dpgrow} will return \code{(n.iter - n.burn)}
  posterior samples.}

  \item{n.thin}{Gap between successive sampling iterations
  to save.}

  \item{strength.mm}{Sets both the shape and rate parameter
  for a \code{tau_{gamma} ~ G(strength.mm,strength.mm)}
  prior on the precision parameter of either a CAR
  (\emph{gamma ~ CAR(tau_gamma)}) or independent
  (\emph{gamma ~ N(0,tau_gamma^(-1)I_S}) prior on the set
  of \code{S} multiple membership effects. Defaults to
  \code{strength.mm = 0.01}.}

  \item{shape.dp}{Shape parameter under a \emph{c ~
  G(shape.dp, 1)} prior on the concentration parameter of
  the DP (prior on the set of random effects parameters,
  \emph{b_1, ..., b_n ~ DP(c,G_0)} where \code{n} is the
  total number of subjects.}

  \item{plot.out}{A boolean variable indicating whether
  user wants to return plots with output results.  Defaults
  to \code{TRUE}.}

  \item{option}{A character vector of length equal to the
  total number of multiple membership terms that supplies
  the prior formulation choice for each term.  The elements
  of \code{option} are confined to choose from among
  \code{c("mmcar","mmi","mmigrp","mmdp")}.  Any element of
  this choice set may be selected multiple times as
  desired.  For example, to add 3 multiple membership terms
  with effects in the first term under a DP prior, the
  second term under a CAR prior and the third also under a
  CAR prior, the entry would be, \code{option =
  c("mmdp","mmcar","mmcar")}.  The corresponding list
  entries should conform this choice for \code{option}.
  The order of sampled effect values returned conforms to
  this order of input in \code{option} (and the
  corresponding \code{subj.aff} and \code{W.subj.aff}).}

  \item{ulabs}{A vector of the same length as \code{option}
  containing desired labels for each multiple membership
  term.  These label values are employed in returned plot
  objects.  If left blank, \code{ulabs} is set to a
  sequential number vector starting at \code{1}.}
}
\value{
  S3 \code{dpgrowmult} object, for which many methods are
  available to return and view results.  Generic functions
  applied to an object, \code{res} of class \code{dpgrow},
  includes: \item{summary(res)}{returns \code{call}, the
  function call made to \code{dpgrowmult} and
  \code{summary.results}, which contains a list of objects
  that include \emph{95\%} credible intervals for each set
  of sampled parameters, specified as (\code{2.5\%}, mean,
  \emph{97.5\%}, including fixed and random effects.  Also
  contains model fit statistics, including \code{DIC} (and
  associated \code{Dbar}, \code{Dhat}, \code{pD},
  \code{pV}), as well as the log pseudo marginal likelihood
  (LPML), a leave-one-out fit statistic.  Note that
  \code{DIC} is constructed as \code{DIC3} (see Celeaux et.
  al. 2006), where the conditional likehihood evaluated at
  the posterior mode is replaced by the marginal predictive
  density.  Lastly, the random and fixed effects design
  matrices, \code{X, Z}, are returned that include both the
  user input nuisance covariates appended to the time and
  treatment-based covariates constructed by
  \code{dpgrowmult}.} \item{print(summary(res))}{prints
  contents of summary to console.} \item{plot(res)}{returns
  results plots, including the set of subject random
  effects values and credible intervals, a sample of
  by-subject growth curves, mean growth curves split by
  each treatment and control, as well as selected trace
  plots for number of clusters and for precision parameters
  for the likehilood and random effects.  Lastly, a trace
  plot for the deviance statistic is also included.}
  \item{samples(res)}{contains (\code{n.iter - n.burn})
  posterior sampling iterations for every model parameter,
  including fixed and random effects.}
  \item{resid(res)}{contains the model residuals.}
}
\description{
  Employs a Dirichlet Process (DP) prior on the set of
  by-subject random effect parameters under repeated waves
  of measurements to allow the number of random effect
  parameters specified per subject, \code{q}, to be equal
  to the number of measurement waves, \code{T}.  Random
  effects are grouped by subject and all \code{q}
  parameters receive the DP prior.  Additional sets of
  possibly more than \code{1} multiple membership effect
  terms are included, each with a separate weight/design
  matrix that maps the effects back to clients.  A variety
  of prior formulations are available for the effects in
  each multiple membership term.
}
\note{
  The intended focus for this package are data where both
  number of subjects and number of repeated measures are
  limited.  A DP prior is placed on the by-subject random
  effects to borrow strength across subjects for each
  estimation of each subject's growth curve.  The
  imposition of the DP prior also allows the resulting
  posterior distributions over the subject random effects
  to be non-Gaussian.  The \code{dpgrowmult} function
  generalizes \code{dpgrowmm} by allowing more than one
  multiple membership effects term.
}
\examples{
\dontrun{
## extract simulated dataset
library(growcurves)
data(datsimmult)
## Model with DP on clients effects, but now INCLUDE session random effects
## in a multiple membership construction communicated with the N x S matrix, W.subj.aff.
## Returns object, res.mm, of class "dpgrowmm".
shape.dp			= 3
res.mult			= dpgrowmult(y = datsimmult$y, subject = datsimmult$subject, trt = datsimmult$trt, time = datsimmult$time,
						n.random = datsimmult$n.random, Omega = datsimmult$Omega, group = datsimmult$group,
						subj.aff = datsimmult$subj.aff, W.subj.aff = datsimmult$W.subj.aff, n.iter = 10000,
						n.burn = 2000, n.thin = 10, shape.dp = shape.dp, option = c("mmi","mmcar"))
plot.results			= plot(res.mult) ## ggplot2 plot objects, including growth curves
summary.results		= summary(res.mult) ## parameter credible intervals and fit statistics
samples.posterior		= samples(res.mult) ## posterior sampled values
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
  Groups, invited re-submission to: JRSS Series A
  (Statistics in Society).

  T. D. Savitsky and S. M. Paddock (2012) Visual Sufficient
  Statistics for Repeated Measures data with growcurves for
  R, submitted to: Journal of Statistical Software.
}
\seealso{
  \code{\link{dpgrowmm}}
}
\keyword{model}

