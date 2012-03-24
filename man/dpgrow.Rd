\name{dpgrow}
\alias{dpgrow}
\alias{dpgrow.default}
\title{Bayesian semiparametric growth curve models.}
\usage{
  dpgrow(y, subject, trt, time, n.random, n.fix_degree,
    formula, random.only, data, n.iter, n.burn, n.thin,
    shape.dp, plot.out, option)
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
  group assignments for each case.  May also be input as
  length \code{P} vector, where \code{P} is the number of
  unique subjects, indicating subject group assignment.
  Multiple treatment groups are allowed and if the vector
  is entered as numeric, e.g. \code{(0,1,2,3,..)}, the
  lowest numbered group is taken as baseline (captured by
  global fixed effects).  If entered in character format,
  the first treatment entry is taken as baseline.  If the
  are no treatment (vs. control) groups, then this input
  may be excluded (set to NULL).}

  \item{time}{A univariate vector of length \code{N},
  capturing the time points associated to each by-subject
  measure.  Mav leave blank if only one time point (no
  repeated measures).}

  \item{n.random}{The desired number of subject random
  effect terms, \code{q}.  Under \code{option = "dp"} may
  be set equal to the number of measurement waves,
  \code{T}.  The \code{y, trt, time} vectors will together
  be used to create both fixed and random effect design
  matrices.  The random effects matrix will be of the the
  form, \code{(1, time, ... , time^(n.random - 1))}
  (grouped, by \code{subject}).  This formulation is a
  growth curve model that allows assessment of by-treatment
  effects and by-client growth curves.}

  \item{n.fix_degree}{The desired polynomial order in time
  to use for generating time-based fix effects.  The fixed
  effects matrix will be constructed as, \code{(time, ...,
  time^(n.fix_degree), trt_1,time*trt_1, ...
  ,time^(n.fix_degree)*trt_l, trt_L,...,
  time^(n.fix_degree)*trt_L)}.  If
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

  \item{n.iter}{Total number of MCMC iterations.}

  \item{n.burn}{Number of MCMC iterations to discard.
  \code{dpgrow} will return \code{(n.iter - n.burn)}
  posterior samples.}

  \item{n.thin}{Gap between successive sampling iterations
  to save.}

  \item{shape.dp}{Shape parameter under a \emph{c ~
  G(shape.dp, 1)} prior on the concentration parameter of
  the DP (prior on the set of random effects parameters,
  \emph{b_1, ..., b_n ~ DP(c,G_0)} where \code{n} is the
  total number of subjects.}

  \item{plot.out}{A boolean variable indicating whether
  user wants to return plots with output results.  Defaults
  to \code{TRUE}.}

  \item{option}{Modeling option, of which there are two: 1.
  \code{dp} places a DP prior on the set of subject random
  effects; 2. \code{lgm} places the usual independent
  Gaussian priors on the set of random effects.}
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
  (LPML), a leave-one-out fit statistic.  Note that for
  \code{option = "dp"}, \code{DIC} is constructed as
  \code{DIC3} (see Celeaux et. al. 2006), where the
  conditional likehihood evaluated at the posterior mode is
  replaced by the marginal predictive density.  Lastly, the
  random and fixed effects design matrices, \code{X, Z},
  are returned that include both the user input nuisance
  covariates appended to the time and treatment-based
  covariates constructed by \code{dpgrow}.}
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
  Employs a Dirichlet Process (DP) prior on the set of
  by-subject random effect parameters under repeated waves
  of measurements to allow the number of random effect
  parameters specified per subject, \code{q}, to be equal
  to the number of measurement waves, \code{T}.  Random
  effects are grouped by subject and all \code{q}
  parameters receive the DP prior.  The resulting joint
  marginal distribution over the data is a DP mixture.
}
\note{
  The intended focus for this package are data where both
  number of subjects and number of repeated measures are
  limited.  A DP prior is placed on the by-subject random
  effects to borrow strength across subjects for each
  estimation of each subject's growth curve.  The
  imposition of the DP prior also allows the resulting
  posterior distributions over the subject random effects
  to be non-Gaussian.  The \code{dpgrow} function is very
  similar to \code{dpgrowmm}; only the latter includes a
  separate set of random effects not grouped by subject
  (e.g. for treatment dosages allocated to subjects) mapped
  back to subject-time cases through a multiple membership
  design matrix.  The \code{dpgrowmult} function
  generalizes \code{dpgrowmm} by allowing more than one
  multiple membership effects term.  See Savitsky and
  Paddock (2011) for detailed model constructions.
}
\examples{
\dontrun{
## extract simulated dataset
library(growcurves)
data(datsim)
## attach(datsim)
## run dpgrow mixed effects model; returns object of class "dpgrow"
shape.dp			= 4
res				= dpgrow(y = datsim$y, subject = datsim$subject, trt = datsim$trt, time = datsim$time,
						n.random = datsim$n.random, n.fix_degree = 2, n.iter = 10000, n.burn = 2000,
						n.thin = 10, shape.dp = shape.dp, option = "dp")
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
  S. M. Paddock and T. D. Savitsky (2012) Bayesian
  Hierarchical Semiparametric Modeling of Longitudinal
  Post-treatment Outcomes from Open-enrollment Therapy
  Groups, submitted to: JRSS Series A (Statistics in
  Society).

  T. D. Savitsky and S. M. Paddock (2011) Visual Sufficient
  Statistics for Repeated Measures data with growcurves for
  R, submitted to: Journal of Statistical Software.
}
\seealso{
  \code{\link{dpgrowmm}}
}
\keyword{model}

