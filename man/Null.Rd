% file MASS/man/Null.Rd
% copyright (C) 1994-2014 W. N. Venables and B. D. Ripley
%
\name{Null}
\alias{Null}
\title{
  Null Spaces of Matrices
}
\description{
  Given a matrix, \code{M}, find a matrix \code{N} giving a basis for the
  (left) null space.  That is \code{crossprod(N, M) = t(N) \%*\% M}
  is an all-zero matrix and \code{N} has the maximum number of linearly
  independent columns.
}
\usage{
Null(M)
}
\arguments{
  \item{M}{
    Input matrix.  A vector is coerced to a 1-column matrix.
  }
}
\details{
  For a basis for the (right) null space
  \eqn{\{x : Mx = 0\}}{{x : Mx = 0}},
  use \code{Null(t(M))}.
}
\value{
  The matrix \code{N} with the basis for the (left) null space, or a
  matrix with zero columns if the matrix \code{M} is square and of
  maximal rank.
}
\references{
  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.} Fourth edition.  Springer.
}
\seealso{
  % avoid link to Matrix
  \code{\link[base]{qr}}, \code{\link{qr.Q}}.
}
\examples{
# The function is currently defined as
function(M)
{
    tmp <- qr(M)
    set <- if(tmp$rank == 0L) seq_len(ncol(M)) else  -seq_len(tmp$rank)
    qr.Q(tmp, complete = TRUE)[, set, drop = FALSE]
}
}
\keyword{algebra}
