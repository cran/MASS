% file MASS/man/newcomb.Rd
% copyright (C) 1994-2020 W. N. Venables and B. D. Ripley
%
\name{newcomb}
\alias{newcomb}
\title{
  Newcomb's Measurements of the Passage Time of Light
}
\description{
  A numeric vector giving the \sQuote{Third Series} of measurements of the
  passage time of light recorded by Newcomb in 1882.  The given values
  divided by 1000 plus 24.8 give the time in millionths of a second for
  light to traverse a known distance. The \sQuote{true} value is now
  considered to be 33.02.

  The dataset is given in the order in Staudte and Sheather.  Stigler
  (1977, Table 5) gives the dataset as
\preformatted{
    28 26 33 24 34 -44 27 16 40 -2 29 22 24 21 25 30 23 29 31 19
    24 20 36 32 36 28 25 21 28 29 37 25 28 26 30 32 36 26 30 22
    36 23 27 27 28 27 31 27 26 33 26 32 32 24 39 28 24 25 32 25
    29 27 28 29 16 23
}
  However, order is not relevant to its use as an example of robust
  estimation.  (Thanks to Anthony Unwin for bringing this difference to our
  attention.)
}
\usage{
newcomb
}
\source{
  S. M. Stigler (1973)
  Simon Newcomb, Percy Daniell, and the history of robust estimation
  1885--1920.
  \emph{Journal of the American Statistical Association} \bold{68}, 872--879.

  S. M. Stigler (1977)
  Do robust estimators work with \emph{real} data?
  \emph{Annals of Statistics}, \bold{5}, 1055--1098.

  R. G. Staudte and S. J. Sheather (1990)
  \emph{Robust Estimation and Testing.} Wiley.
}
\keyword{datasets}
