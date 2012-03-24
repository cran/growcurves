\name{relabel}
\alias{relabel}
\title{Relabel user vector input to sequential numerical}
\usage{
  relabel(label.input, start)
}
\arguments{
  \item{label.input}{A vector in numerical or character
  format that the user desires to convert to sequential
  numeric.}

  \item{start}{An integer representing the starting value
  of the sequential sequence for the new label vector.}
}
\value{
  A list object containing the new sequential vector,
  \code{label.new}, and a vector of unique input values,
  \code{labeli.u}.
}
\description{
  Takes user input for vector in either character or
  numerical format and converts it to sequential numeric
  format. For use in functions requiring sequential
  numerical format.  Returns new sequential numerical
  vector and vector of unique values inputted by user.  The
  latter is used to label plot variables returned to the
  user.
}
\author{
  Terrance Savitsky \email{tds151@gmail.com}
}
\seealso{
  \code{\link{dpgrowmm}}, \code{\link{dpgrow}}
}

