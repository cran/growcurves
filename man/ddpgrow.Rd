\name{ddpgrow}
\alias{ddpgrow}
\alias{ddpgrow.default}
\title{Bayesian semiparametric growth curve models.}
\usage{
  ddpgrow(y, subject, trt, time, n.random, n.fix_degree,
    formula, random.only, data, dosemat, numdose,
    typetreat, labt, Omega, n.iter, n.burn, n.thin,
    shape.dp, rate.dp, M.init, plot.out)
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

  \item{trt}{An integer or character matrix/vector of
  length \code{N} (number of cases) indicating treatment
  arm assignments for each case.  May also be input as
  length \code{n} vector, where \code{n} is the number of
  unique subjects, indicating subject arm assignment.
  Multiple treatment arms are allowed and if the vector is
  entered as numeric, e.g. \code{(0,1,2,3,..)}, the lowest
  numbered arm is taken as baseline (captured by global
  fixed effects).  If entered in character format, the
  first treatment entry is taken as baseline.  If the are
  no treatment (vs. control) arm, then this vector is
  composed of a single value for all entries.}

  \item{time}{A univariate vector of length \code{N},
  capturing the time points associated to each by-subject
  measure.}

  \item{n.random}{The desired number of subject random
  effect terms, \code{q}.  May be set equal to the number
  of measurement waves, \code{T}.  The \code{y, trt, time}
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
  time^(n.fix_degree)*trt_L)}.}

  \item{formula}{Nuisance fixed and random effects may be
  entered in \code{formula} with the following format,
  \code{y ~ x_1 + x_2*x_3 | z_1*z_2 }.  The bar, \code{|},
  separates fixed and random effects.  If it is only
  desired to enter either fixed or random effects, but not
  both then the \code{|} may be omitted.  Note: the
  nuisance random effects are assumed to be grouped by
  subject.  The fixed and random effects valules may change
  with each repeated measure; however, within subject
  growth curves will keep constant \code{z} and \code{x}
  values between measurement waves.  It is possible to
  bypass the growth curve construction by leaving \code{y,
  trt, time, n.random} blank and entering only
  \code{formula}, instead.  The model output plots, will,
  however exclude growth curves in that event.  If a
  formula is input and a response, \code{y}, is included,
  then the parameter input \code{y} may be omitted.  If the
  \code{y} input is included, it will be over-written by
  that from \code{formula}.}

  \item{random.only}{A Boolean variable indicating whether
  the input formula contains random (for fixed) effects in
  the case that only one set are entered.  If excluded and
  \code{formula} is entered without a \code{|},
  \code{random.only} defaults to \code{FALSE}.}

  \item{data}{A \code{data.frame} containing the variables
  named in \code{formula}.}

  \item{dosemat}{An \code{P x (T+1)} \code{matrix} object
  that maps \code{subjects} to treatment dosages.  The
  first column should be an intercept column (filled with
  1's).}

  \item{numdose}{A vector object containing the number of
  dosages for each treatment.  So the length should be the
  same as \code{typetreat}.}

  \item{typetreat}{A vector object specifying the prior
  formulation for each treatment.  The choices for prior
  formulations are \code{c("car","mvn","ind")}.}

  \item{labt}{An optional vector object (of the same length
  as \code{typetreat}) providing user names for each
  treatment.  The names are used in plot objects.  If
  \code{NULL}, then the numerical order of treatment
  entries are used.}

  \item{Omega}{A list object of length equal to the number
  of treatments chosen with the \code{"car" \%in\%
  typetreat}.  List element \code{m} contains an
  \emph{numdose[m] x numdose[m]} numeric matrix to encode
  the CAR adjacency matrix, where \code{numdose[m]} is the
  number of dosages receiving the multivariate CAR prior.
  This input is required only under \code{"car" \%in\%
  typetreat}.}

  \item{n.iter}{Total number of MCMC iterations.}

  \item{n.burn}{Number of MCMC iterations to discard.
  \code{ddpgrow} will return \code{(n.iter - n.burn)}
  posterior samples.}

  \item{n.thin}{Gap between successive sampling iterations
  to save.}

  \item{shape.dp}{Shape parameter under a \emph{c ~
  G(shape.dp, rate.dp)} prior on the concentration
  parameter of the DP (prior on the set of random effects
  parameters, \emph{b_1, ..., b_n ~ DP(c,G_0)} where
  \code{n} is the total number of subjects.}

  \item{rate.dp}{Rate parameter under a \emph{c ~
  G(shape.dp, rate.dp)} prior on the concentration
  parameter of the DP.}

  \item{M.init}{Scalar value capturing number of initial
  subject clusters to kick-off MCMC chain.}

  \item{plot.out}{A boolean variable indicating whether
  user wants to return plots with output results.  Defaults
  to \code{TRUE}.}
}
\value{
  S3 \code{dpgrow} object, for which many methods are
  available to return and view results.  Generic functions
  applied to an object, \code{res} of class \code{dpgrow},
  includes: \item{summary(res)}{ returns \code{call}, the
  function call made to \code{dpgrow} and
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
  treatment-based covariates constructed by \code{dpgrow}.}
  \item{print(summary(res))}{ prints contents of summary to
  console.} \item{plot(res)}{ returns results plots,
  including the set of subject random effects values and
  credible intervals, a sample of by-subject growth curves,
  mean growth curves split by each treatment and control,
  as well as selected trace plots for number of clusters
  and for precision parameters for the likehilood and
  random effects.  Lastly, a trace plot for the deviance
  statistic is also included.} \item{samples(res)}{
  contains (\code{n.iter - n.burn}) posterior sampling
  iterations for every model parameter, including fixed and
  random effects.} \item{resid(res)}{ contains the model
  residuals.}
}
\description{
  Employs an anova Dependent Dirichlet Process (DDP) prior
  on the set of by-subject random effect parameters with
  dependence indexed by multiple membership effects under
  repeated waves of measurements to allow the number of
  random effect parameters specified per subject, \code{q},
  to be equal to the number of measurement waves, \code{T}.
  Random effects are grouped by subject and all \code{q}
  parameters receive the DP prior.  The resulting joint
  marginal distribution over the data is a DP mixture.
}
\note{
  The intended focus for this package are data where both
  number of subjects and number of repeated measures are
  limited.  This function places a DDP prior on the set of
  subject effects.  This means that the unnknown (random)
  prior on subject effects is indexed by the subject dosage
  patterns across one or more treatments.
}
\examples{
\dontrun{
## extract simulated dataset
library(growcurves)
data(datddpsim)
## attach(datddpsim)
## run dpgrow mixed effects model; returns object of class "ddpgrow"
shape.dp			= 4
res				= ddpgrow(y = dat$y, subject = dat$subject, trt = dat$trt, time = dat$time,
						typetreat = c("mvn","car","ind","car"), numdose = dat$numdose,
						labt = dat$labt, dosemat = dat$dosemat, Omega = dat$Omega, n.random = dat$n.random,
						n.fix_degree = 2, n.iter = 10000, n.burn = 2000,
						n.thin = 10, shape.dp = 1)
plot.results			= plot(res) ## ggplot2 plot objects, including growth curves
summary.results		= summary(res) ## parameter credible intervals and fit statistics
samples.posterior		= samples(res) ## posterior sampled values
}
}
\author{
  Terrance Savitsky \email{tds151@gmail.com} Susan Paddock
  \email{paddock@rand.org}
}
\references{
  T. D. Savitsky and S. M. Paddock (2011) Bayesian
  Hierarchical Semiparametric Modeling of Longitudinal
  Post-treatment Outcomes from Open-enrollment Therapy
  Groups, submitted to: JRSS Series A (Statistics in
  Society).

  T. D. Savitsky and S. M. Paddock (2011) Visual Sufficient
  Statistics for Repeated Measures data with growcurves for
  R, submitted to: Journal of Statistical Software.
}
\seealso{
  \code{\link{dpgrowmult}}, \code{\link{dpgrowmm}},
  \code{\link{dpgrow}}
}
\keyword{model}

