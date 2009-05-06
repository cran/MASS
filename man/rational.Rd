% file MASS/man/rational.Rd
% copyright (C) 1994-9 W. N. Venables and B. D. Ripley
%
\name{rational}
\alias{rational}
\alias{.rat}
\title{
Rational Approximation
}
\description{
Find rational approximations to the components of a real numeric
object using a standard continued fraction method.
}
\usage{
rational(x, cycles = 10, max.denominator = 2000, \dots)
}
\arguments{
\item{x}{
Any object of mode numeric. Missing values are now allowed.
}
\item{cycles}{
The maximum number of steps to be used in the continued fraction
approximation process.
}
\item{max.denominator}{
An early termination criterion.  If any partial denominator
exceeds \code{max.denominator} the continued fraction stops at that point.
}
\item{\dots}{
arguments passed to or from other methods.
}}
\value{
A numeric object with the same attributes as \code{x} but with entries
rational approximations to the values.  This effectively rounds
relative to the size of the object and replaces very small
entries by zero.
}
\details{
Each component is first expanded in a continued fraction of the
form


\code{x = floor(x) + 1/(p1 + 1/(p2 + \dots)))}


where \code{p1}, \code{p2}, \dots are positive integers, terminating either
at \code{cycles} terms or when a \code{pj > max.denominator}.  The
continued fraction is then re-arranged to retrieve the numerator
and denominator as integers and the ratio returned as the value.
}
\seealso{
\code{\link{fractions}}
}
\examples{
X <- matrix(runif(25), 5, 5)
solve(X, X/5)
##             [,1]        [,2]       [,3]        [,4]        [,5]
## [1,]  2.0000e-01  3.7199e-17 1.2214e-16  5.7887e-17 -8.7841e-17
## [2,] -1.1473e-16  2.0000e-01 7.0955e-17  2.0300e-17 -1.0566e-16
## [3,]  2.7975e-16  1.3653e-17 2.0000e-01 -1.3397e-16  1.5577e-16
## [4,] -2.9196e-16  2.0412e-17 1.5618e-16  2.0000e-01 -2.1921e-16
## [5,] -3.6476e-17 -3.6430e-17 3.6432e-17  4.7690e-17  2.0000e-01

## rational(solve(X, X/5))
##      [,1] [,2] [,3] [,4] [,5]
## [1,]  0.2  0.0  0.0  0.0  0.0
## [2,]  0.0  0.2  0.0  0.0  0.0
## [3,]  0.0  0.0  0.2  0.0  0.0
## [4,]  0.0  0.0  0.0  0.2  0.0
## [5,]  0.0  0.0  0.0  0.0  0.2
}
\keyword{math}
