\name{glmmPQL}
\alias{glmmPQL}
\title{
  Fit Generalized Linear Mixed Models via PQL
}
\description{
  Fit a GLMM model with multivariate normal random effects, using
  Penalized Quasi-Likelihood.
}
\usage{
glmmPQL(fixed, random, family, data, correlation, weights,
        control, niter = 10, verbose = TRUE, \dots)
}
\arguments{
\item{fixed}{
a two-sided linear formula giving fixed-effects part of the model.
}
\item{random}{
a formula or list of formulae describing the random effects.
}
\item{family}{
a GLM family.
}
\item{data}{
an optional data frame, list or environment used as the first place to find
variables in the formulae, \code{weights} and if present in
\code{\dots}, \code{subset}.
}
\item{correlation}{
an optional correlation structure.
}
\item{weights}{
optional case weights as in \code{glm}.
}
\item{control}{
an optional argument to be passed to \code{lme}.
}
\item{niter}{
maximum number of iterations.
}
\item{verbose}{
logical: print out record of iterations?
}
\item{\dots}{
Further arguments for \code{lme}.
}}
\value{
  A object of class \code{c("glmmPQL", "lme")}: see \code{\link[nlme]{lmeObject}}.
}
\details{
  \code{glmmPQL} works by repeated calls to \code{\link[nlme]{lme}}, so
  namespace \CRANpkg{nlme} will be loaded at first use.  (Before 2015 it
  used to attach \code{nlme} but nowadays only loads the namespace.)

  Unlike \code{lme}, \code{\link{offset}} terms are allowed in
  \code{fixed} -- this is done by pre- and post-processing the calls to
  \code{lme}.

  Note that the returned object inherits from class \code{"lme"} and
  that most generics will use the method for that class.  As from
  version 3.1-158, the fitted values have any offset included, as do
  the results of calling \code{\link{predict}}.
}
\references{
  Schall, R. (1991) Estimation in generalized linear models with
  random effects.
  \emph{Biometrika}
  \bold{78}, 719--727.

  Breslow, N. E. and Clayton, D. G. (1993) Approximate inference in
  generalized linear mixed models.
  \emph{Journal of the American Statistical Association}
  \bold{88}, 9--25.

  Wolfinger, R. and O'Connell, M. (1993) Generalized linear mixed models: a
  pseudo-likelihood approach.
  \emph{Journal of Statistical Computation and Simulation}
  \bold{48}, 233--243.

  Venables, W. N. and Ripley, B. D. (2002)
  \emph{Modern Applied Statistics with S.} Fourth edition.  Springer.
}
\seealso{
  \code{\link[nlme]{lme}}
}
\examples{
summary(glmmPQL(y ~ trt + I(week > 2), random = ~ 1 | ID,
                family = binomial, data = bacteria))

## an example of an offset: the coefficient of 'week' changes by one.
summary(glmmPQL(y ~ trt + week, random = ~ 1 | ID,
               family = binomial, data = bacteria))
summary(glmmPQL(y ~ trt + week + offset(week), random = ~ 1 | ID,
                family = binomial, data = bacteria))
}
\keyword{models}
