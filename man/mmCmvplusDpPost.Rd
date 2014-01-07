\name{mmCmvplusDpPost}
\alias{mmCmv}
\alias{mmCmvplusDpPost}
\title{Bayesian mixed effects model with a DP prior on by-subject effects and CAR prior on a multivariate set of multiple membership effects}
\usage{
mmCmvplusDpPost(y, X, Z, H, Wcase, Wsubject, Omega, omegaplus, groups, subjects,
  niter, nburn, nthin, strength.mm, corsess, shapealph, ratebeta, typemm)
}
\arguments{
  \item{y}{An \emph{N x 1} response (of subject-measure
  cases)}

  \item{X}{Fixed effects design matrix}

  \item{Z}{Random effects design matrix.  Assumed grouped
  by \code{subjects}}

  \item{H}{Multivariate MM effects design matrix.}

  \item{Wcase}{An \emph{N x 1} multiple membership weight
  matrix to map supplemental random effects}

  \item{Wsubject}{An \emph{P.aff x S} multiple membership
  weight matrix with rows equal to number of unique
  affected subjects}

  \item{Omega}{An \emph{S x S} unnormalized adjacency
  matrix with entries equal to 1 where two effects
  communicate and 0, otherwise.  Diagonal elements are
  zero}

  \item{omegaplus}{\emph{S x 1} vector of row sums of
  \code{Omega}}

  \item{groups}{\emph{S x 1} vector of group identifiers
  for each effect.  Effects within each group communicate.
  Effects don't communicate across groups.  Not used under
  "mmi" prior (though input is required).}

  \item{subjects}{An \emph{N x 1} set of subject
  identifiers}

  \item{niter}{The number of MCMC iterations}

  \item{nburn}{The number of MCMC burn-in iterations to
  discard}

  \item{nthin}{The step increment of MCMC samples to
  return}

  \item{strength.mm}{The shape and rate parameters for the
  \eqn{\Gamma} prior on the CAR precision parameter,
  \eqn{\tau_{\gamma}}.}

  \item{corsess}{A single value to set the prior
  correlations among the multivariate \code{q = ncol(H)}
  orders for the MM effects.  where \eqn{\tau_{\gamma}} is
  replaced by the \code{q x q}, \eqn{\Lambda}.}

  \item{shapealph}{The shape parameter for the \eqn{\Gamma}
  prior on the DP concentration parameter.  The rate
  parameter is set of \code{1}.}

  \item{ratebeta}{The rate parameter for the \eqn{\Gamma}
  prior on the DP concentration parameter. Default value is
  \code{1}.}

  \item{typemm}{An indicator the prior formulation
  specified for the multivariate MM effects term.  Set
  \code{typemm = 0} for \code{"mmi"} and \code{typemm = 1}
  for \code{"mmcar"}.}
}
\value{
res A list object containing MCMC runs for all model
parameters.
}
\description{
An internal function to \code{\link{dpgrowmm}}
}
\note{
Intended as an internal function for \code{\link{dpgrowmm}}
}
\author{
Terrance Savitsky \email{tds151@gmail.com}
}
\seealso{
\code{\link{dpgrow}}
}

