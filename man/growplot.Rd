\name{growplot}
\alias{growplot}
\title{Plot by-subject and by-group growth curves}
\usage{
  growplot(object, compare.objects = NULL,
    subjects.plot = NULL, groups.plot = NULL,
    subjects.subset = NULL, subjects.random = TRUE,
    x.lab = "Time", y.lab = "Fit",
    main.label = "Main_Model", title.lab = NULL)
}
\arguments{
  \item{object}{A \code{dpgrow}, \code{dgrowmm} or
  \code{dpgrowmult} object obtained from running the
  appropriate modeling function.}

  \item{compare.objects}{An optional list of \code{dpgrow},
  \code{dgrowmm} or \code{dpgrowmm} objects to employ as
  comparison models for by-subject plotting}

  \item{subjects.plot}{A vector of subjects to use for
  by-group plot that is composed of some subset of
  \code{subject}.  Defaults to all subjects modeled in
  \code{object}.}

  \item{groups.plot}{A vector with associated group
  identifiers for \code{subjects.plot} if the grouping is
  different from \code{trt} used to run the \code{dpgrow}
  or \code{dpgrowmm} model.  The entered grouping does not
  have to relate to that used for modeling.  Defaults to
  use of treatment groups modeled in \code{object} (if not
  entered).}

  \item{subjects.subset}{A vector of a subset of
  \code{subjects.plot} to use for by-subject plotting for
  readability.  The full \code{subjects.plot} set is used
  for by group plotting. If left blank, the full
  \code{subjects.plot} vector is used for by-subject
  plotting.}

  \item{subjects.random}{A boolean scalar.  If \code{TRUE}
  a random subset (of 10) is selected from
  \code{subjects.plot} for by-subject plotting.  Leave
  blank if enter \code{subjects.subset}.  Defaults to
  \code{TRUE}.}

  \item{x.lab}{Optional title for x-axis.  Default =
  "Time"}

  \item{y.lab}{Optional title for y-axis.  Default = "Fit"}

  \item{main.label}{An optional character model label for
  \code{object} to use in by-subject plots for comparison
  with models from \code{compare.objects}}

  \item{title.lab}{option plot title. A vector of 2
  character entries is allowed.  The first entry is the
  title for the group aggregated plot.  The second entry is
  the title for the plot of selected subject growth curves.
  If \code{title.lab} contains a single entry, then it is
  used as the title for both plots.  Otherwise, defaults to
  \code{NULL}.}
}
\value{
  A list object containing a plot objects, data.frame
  object from which it is constructed, and data.frame with
  actual data values co-plotted.  \item{dat.gc}{A
  \code{data.frame} object used to generate the within
  subject predicted growth curves for \code{object}.  Used
  for aggregated plot of growth curves by group returned in
  \code{p.gctrt} Fields are titled,
  \code{c("fit","time","subject","trt")}.} \item{dat.igc}{A
  \code{data.frame} object containing within subject
  predicted growth curves under models in \code{object} and
  \code{compare.objects}.  Used for by-subject growth
  curves plot returned in \code{p.gcsub}.}
  \item{dat.data}{A \code{data.frame} object containing the
  actual data observations for plotted subjects.  Field
  titles are the same as for \code{dat.gc}.}
  \item{p.gctrt}{A \code{ggplot2} object of subjects
  aggregated by group.} \item{p.gcsub}{A \code{ggplot2}
  object of subjects.}
}
\description{
  Produces a by-subject plot of predicted growth curves
  with associated data values and an aggregated by-group
  plot of growth curves along with a smoother line for each
  group based on user input.  Facilates inference for
  different growth curve patterns based on subsets of
  subjects and various subject groupings.
}
\author{
  Terrance Savitsky \email{tds151@gmail.com}
}
\seealso{
  \code{\link{dpgrowmm}}, \code{\link{dpgrow}},
  \code{\link{dpgrowmult}}, \code{\link{growthCurve}}
}

