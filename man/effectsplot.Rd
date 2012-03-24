\name{effectsplot}
\alias{effectsplot}
\title{Plot comparison of Effect parameters of a Multiple membership (MM) term under varied prior formulations}
\usage{
  effectsplot(objects, mm.terms = NULL, prior.labs = NULL,
    center = TRUE, axis.labs = NULL, group.plot = NULL,
    smoother = FALSE)
}
\arguments{
  \item{objects}{A list object where each entry is either a
  \code{dpgrowmm} or \code{dgrowmult} object returned from
  a model run.  The list length is equal to the number of
  models on which a comparison of multiple membership term
  effects is desired.}

  \item{mm.terms}{An optional character vector of MM term
  names or numeric vector of MM term positions targeted for
  plotting.  The entries are term names for
  \code{dpgrowmult} objects and group names for
  \code{dpgrowmm} objects.  If a character vector is used,
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
  chosen object terms.} \item{p.all}{A \code{ggplot2} only
  rendered where all objects are of class \code{dpgrowmult}
  and have common number of terms and effects.  Plot
  compares effects across objects.} \item{dat.all}{A
  \code{data.frame} object used to generate \code{p.all}
  with field titles,
  \code{c("prior","session","block","low","mean","high")},
  where \code{block} denotes labels for MM terms or groups
  included.}
}
\description{
  Produces an overlaid plot of multiple membership effect
  values under chosen prior formulations.  The effect
  values for a focus (or all) multiple membership term(s)
  may be compared under different model formulations (one
  vs. multiple terms) and for differing prior formulations
  (e.g. "mmcar" vs. "mmdp").
}
\author{
  Terrance Savitsky \email{tds151@gmail.com}
}
\seealso{
  \code{\link{dpgrowmm}}, \code{\link{dpgrowmult}},
}

