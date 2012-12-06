\name{effectsplot}
\alias{effectsplot}
\title{Plot comparison of Effect parameters of a Multiple membership (MM) term under varied prior formulations}
\usage{
  effectsplot(objects, mm.terms = NULL, prior.labs = NULL,
    center = TRUE, axis.labs = NULL, group.plot = NULL,
    smoother = FALSE, order = FALSE, orderto = NULL,
    label.mm = NULL)
}
\arguments{
  \item{objects}{A list object where each entry is an
  object with class within any of either \code{dpgrowmm} or
  \code{dgrowmult}.  returned from a model run.  The list
  length is equal to the number of models on which a
  comparison of multiple membership term effects is
  desired.}

  \item{mm.terms}{An optional character vector of MM term
  names or numeric vector of MM term positions targeted for
  plotting.  The entries a term name for each
  \code{dpgrowmult} object and group name for each
  \code{dpgrowmm} object.  If a character vector is used,
  names must correspond to those used in package engine
  functions that generated each object.  Length is equal to
  number of \code{objects}.}

  \item{prior.labs}{An optional character vector providing
  names for the prior formulations under comparison.
  Length is equal to number of objects.}

  \item{center}{A scalar boolean input to indicate whether
  effect values in all objects should be centered.
  Defaults to \code{TRUE}.}

  \item{axis.labs}{An optional vector input of length
  \code{2} corresponding to \code{x.axis} and \code{y.axis}
  labels, respectively.}

  \item{group.plot}{An optional numeric vector that groups
  or re-groups the MM effects in each term selected in
  \code{mm.terms}.}

  \item{smoother}{A scalar boolean input indicating whether
  to co-plot a smoother line with point values.}

  \item{order}{A scalar boolean input indicating whether to
  sort the plots of effects in order of increasing (mean)
  value.}

  \item{orderto}{A list object of length equal to number of
  MM terms (for \code{dpgrowmult} objects), or groups (for
  \code{dpgrowmm} objects), each holding a vector of length
  \code{S[n]}, the number of MM effects for object \code{n}
  that conveys a desired order for the MM effects under
  option \code{order = TRUE} as an alternative to ordering
  by increasing mean value.}

  \item{label.mm}{An optional list object, each containing
  a vector of labels for the MM random effects in each term
  (\code{dpgrowmult}) or each group (\code{dpgrowmm}).
  Used for plotting all terms or groups for input objects.
  The length of \code{label.mm} is equal to the number of
  terms or groups per object.}
}
\value{
  A list object containing a plot objects and associated
  plot data.frames \item{dat.term}{A \code{data.frame}
  object used to generate a plot comparing a specific term
  or group from an MM object with field titles,
  \code{c("session","group","low","mean","high","prior")},
  where \code{session} denotes the MM effect labels,
  \code{group} the sub-group labels within the focus MM
  terms, and \code{prior}, the model identifiers.}
  \item{p.term}{A \code{ggplot2} object of effect plots for
  chosen object terms.} \item{p.all}{A \code{ggplot2} plots
  effects, by term and/or compatible group, for each
  object.  Only renders if the number of effects in terms
  and/or groups are compatible between objects.}
  \item{dat.all}{A \code{data.frame} object used to
  generate \code{p.all} with field titles,
  \code{c("prior","session","block","low","mean","high")},
  where \code{block} denotes labels for MM terms or groups
  included.}
}
\description{
  Produces an overlaid plot of multiple membership effect
  values under chosen prior formulations.  The effect
  values for a focus (or all) multiple membership term(s)
  or group(s) may be compared under different model
  formulations (one vs. multiple terms, or multivariate vs.
  univariate) for differing prior formulations (e.g.
  "mmcar" vs. "mmi").
}
\author{
  Terrance Savitsky \email{tds151@gmail.com}
}
\seealso{
  \code{\link{dpgrowmm}}, \code{\link{dpgrowmult}},
  \code{\link{ddpgrow}}
}

