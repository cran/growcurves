\name{getmf}
\alias{get.mf}
\alias{getmf}
\title{Produce fixed and random effects design matrices from single formula input}
\usage{
  getmf(formula, random.only, data, na.action)
}
\arguments{
  \item{formula}{A formula of format \code{y ~ x_1 +
  x_2*x_3 | z_1*z_2 } where \code{|} separates fixed (to
  the left of \code{|}) and random effects.}

  \item{random.only}{A boolean scalar used in the case that
  either fixed or random effects are entered in
  \code{formula}, but not both, which case the \code{|} is
  not entered (e.g. \code{y ~ x_1 + x_2*x_3 }.  Then, if
  \code{random.only == TRUE} the variables on the
  right-hand side are interpreted to be random effects;
  otherwise fixed for use in \code{\link{dpgrow}} and
  \code{\link{dpgrowmm}}.}

  \item{data}{Associated data.frame containing names
  variables in \code{formula}}

  \item{na.action}{Should be left blank for use in
  \code{\link{dpgrow}} and \code{\link{dpgrowmm}}, where is
  automatically set to \code{na.fail}.}
}
\value{
  res A list object containing \code{list(y = y, x = x, z =
  z, m = m, mf = mf)}.
}
\description{
  An internal function to \code{\link{dpgrow}} and
  \code{\link{dpgrowmm}}
}
\note{
  Intended as an internal function for \code{\link{dpgrow}}
  and \code{\link{dpgrowmm}}
}
\author{
  Terrance Savitsky \email{tds151@gmail.com}
}
\seealso{
  \code{\link{dpgrow}} \code{\link{dpgrowmm}}
}

