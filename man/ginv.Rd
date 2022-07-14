% file MASS/man/ginv.Rd
% copyright (C) 1994-2022 W. N. Venables and B. D. Ripley
%
\name{ginv}
\alias{ginv}
\title{
Generalized Inverse of a Matrix
}
\description{
Calculates the Moore-Penrose generalized inverse of a matrix
\code{X}.
}
\usage{
ginv(X, tol = sqrt(.Machine$double.eps))
}
\arguments{
\item{X}{
Matrix for which the Moore-Penrose inverse is required.
}
\item{tol}{
A relative tolerance to detect zero singular values.
}}
\value{
A MP generalized inverse matrix for \code{X}.
}
\references{
  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S-PLUS.} Fourth Edition. Springer.
}
\seealso{
\code{\link{solve}}, \code{\link{svd}}, \code{\link{eigen}}
}
\keyword{algebra}
