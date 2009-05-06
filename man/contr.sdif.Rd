\name{contr.sdif}
\alias{contr.sdif}
\title{
Successive Differences contrast coding
}
\description{
A coding for unordered factors based on successive differences.
}
\usage{
contr.sdif(n, contrasts = TRUE)
}
\arguments{
\item{n}{
The number of levels required.
}
\item{contrasts}{
Should there be \code{n - 1} columns orthogonal to the mean (the default)
or \code{n} columns spanning the space.
}}
\value{
If \code{contrasts is} \code{TRUE}, a matrix with \code{n} rows and \code{n - 1} columns,
and the \code{n} by \code{n} identity matrix if \code{contrasts} is \code{FALSE}.
}
\details{
The contrast coefficients are chosen so that the coded coefficients
in a one-way layout are the differences between the means of the
second and first levels, the third and second levels, and so on.
}
\references{
Venables, W. N. and Ripley, B. D. (2002)
\emph{Modern Applied Statistics with S.}
Fourth Edition, Springer.
}
\seealso{
\code{\link{contr.treatment}}, \code{\link{contr.sum}}, \code{\link{contr.helmert}}
}
\examples{
contr.sdif(6)
}
\keyword{models}
