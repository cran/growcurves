\name{ddpEffectsplot}
\alias{ddpEffectsplot}
\title{Plot by-subject and by-treatment posterior mean values for dossage random effects}
\usage{
ddpEffectsplot(object, subjects.plot = NULL, n.plot = 3, trts.plot = NULL,
  x.axis.label = NULL, smoother = TRUE, re.order = TRUE,
  cred.intervals = TRUE, map.group = NULL, n.dose.plot = 5,
  orderto = NULL)
}
\arguments{
  \item{object}{A \code{ddpgrow} object.}

  \item{subjects.plot}{A vector of subjects for performing
  plots that is composed of some subset of the
  \code{subject} vector input for modeling.  If left blank,
  a random subset is chosen from \code{subject}.}

  \item{n.plot}{An optional scalar input for number of
  randomly generated subjects to plot (if
  \code{subjects.plot} is left blank).}

  \item{trts.plot}{A vector of focus treatments to use for
  plotting.}

  \item{x.axis.label}{An optional scalar character entry to
  label the treatment(s) dosages}

  \item{smoother}{A scalar boolean input indicating whether
  to co-plot a smoother line with point values.}

  \item{re.order}{A scalar boolean input indicating whether
  to sort the plots of effects in order of increasing
  value.}

  \item{cred.intervals}{A boolean scalar indicating whether
  the by-subject effects plots should include credible
  intervals.}

  \item{map.group}{A \code{matrix} or \code{data.frame}
  object containing a grouping of subjects that will be
  used to produce an additional set of effect plots that
  aggregate subjects by the grouping structure.  The first
  column containing subject identifiers for all subjects
  modeled in \code{object}.  The second column contains the
  desired desired group identifiers that may be of type
  character or numeric.}

  \item{n.dose.plot}{Optional numeric input for number of
  randomly chosen doses for which to plot effects growth
  curves.}

  \item{orderto}{A numeric vector of length equal to the
  total number of dosages across all treatments that
  conveys an order to be used by-dosage growth curve plots
  within cluster and treatment.}
}
\value{
A list object containing a faceted set of heat maps (one
per subject), a faceted set of effect point plots, and the
associated \code{data.frame} objct.  \item{dat.se}{A
\code{data.frame} object used to generate the trt-subject
faceted plots for effect means Fields are titled,
\code{c("order","dose","trt","subject","effects")}.}
\item{dat.ci}{A \code{data.frame} object used to generate
the trt-subject faceted plots for effect credible intervals
Fields are titled,
\code{c("order","dose","trt","subject","quantile","effects")}.}
\item{dat.clust}{A \code{data.frame} object used to
generate the trt-group faceted plots for effect means
Fields are titled,
\code{c("order","dose","trt","cluster","effects")}.}
\item{dat.clust.ci}{A \code{data.frame} object used to
generate the trt-group faceted plots for effect credible
intervals Fields are titled,
\code{c("order","dose","trt","cluster","quantile","effects")}.}
\item{dat.gc}{A \code{data.frame} object used to generate
the by-dose growth curves (for multivariate polynomial
effects) Fields are titled,
\code{c("fit","time","cluster","trt","dose")}.}
\item{p.hm}{A \code{ggplot2} object of heat maps for mean
random effect values, faceted by trt and subject
combinations.} \item{p.pp}{A \code{ggplot2} object of point
plots for mean random effect values or credible intervals,
faceted by trt and subject combinations.} \item{pc.m}{A
\code{ggplot2} object of point plots for mean random effect
values, faceted by trt and group combinations.}
\item{pc.ci}{A \code{ggplot2} object of point plots for
credible intervals of random effect values, faceted by trt
and cluster combinations.} \item{pc.gc}{A \code{ggplot2}
growth curve plots for each dose where the does effects are
multivariate polynomial.}
}
\description{
Each \code{ddpgrow} object contains posterior mean
estimates for the \code{q} x \code{T} matrix of by-subject
random effects that is extracted from the \code{ddpgrow}
object that is input to \code{ddpEffectsplot}. This
function produces a \code{q} x \code{T.m} heat map plot of
posterior mean effect values for the dosages in treatment
\code{m} faceted on a set of chosen subjects.  The
resulting plot produces a heatmap for each trt-subject
combination.  Both a ggplot2 plot object and a
\code{data.frame} object are returned.
}
\author{
Terrance Savitsky \email{tds151@gmail.com}
}
\seealso{
\code{\link{ddpgrow}}, \code{\link{dpgrow}},
\code{\link{dpgrowmm}}, \code{\link{dpgrowmult}}
}

