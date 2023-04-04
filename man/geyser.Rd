% file MASS/man/geyser.Rd
% copyright (C) 1994-2014 W. N. Venables and B. D. Ripley
%
\name{geyser}
\alias{geyser}
\usage{
geyser
}
\title{Old Faithful Geyser Data}
\description{
  A version of the eruptions data from the \sQuote{Old Faithful} geyser
  in Yellowstone National  Park,  Wyoming. This version comes from
  Azzalini and Bowman (1990) and is of continuous measurement from August
  1 to August 15, 1985.

  Some nocturnal duration measurements were coded as 2, 3 or 4 minutes,
  having originally been described as \sQuote{short}, \sQuote{medium}
  or \sQuote{long}.
}
\format{A data frame with 299 observations on 2 variables.
  \tabular{lll}{
    \code{duration}  \tab numeric  \tab Eruption time in mins \cr
    \code{waiting}   \tab numeric  \tab Waiting time for this eruption \cr
  }
}
\note{
  The \code{waiting} time was incorrectly described as the time to the
  next eruption in the original files, and corrected for \pkg{MASS}
  version 7.3-30.
}
\seealso{
  \code{\link{faithful}}.

  CRAN package \pkg{sm}.
}
\references{
  Azzalini, A. and Bowman, A. W. (1990) A look at some
  data on the Old Faithful geyser.  \emph{Applied Statistics}
  \bold{39}, 357--365.

  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.} Fourth edition.  Springer.

}
\keyword{datasets}
