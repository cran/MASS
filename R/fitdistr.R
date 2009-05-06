# file MASS/R/fitdistr.R
# copyright (C) 2002-2007 W. N. Venables and B. D. Ripley
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
fitdistr <- function(x, densfun, start, ...)
{
    myfn <- function(parm, ...) -sum(log(dens(parm, ...)))
    mylogfn <- function(parm, ...) -sum(dens(parm, ..., log = TRUE))
    mydt <- function(x, m, s, df, log) dt((x-m)/s, df, log = TRUE) - log(s)

    Call <- match.call(expand.dots=TRUE)
    if(missing(start)) start <- NULL
    dots <- names(list(...))
    dots <- dots[!is.element(dots, c("upper", "lower"))]
    if(missing(x) || length(x) == 0L || mode(x) != "numeric")
        stop("'x' must be a non-empty numeric vector")
    if(missing(densfun) || !(is.function(densfun) || is.character(densfun)))
        stop("'densfun' must be supplied as a function or name")
    n <- length(x)
    if(is.character(densfun)) {
        distname <- tolower(densfun)
        densfun <-
            switch(distname,
                   "beta" = dbeta,
                   "cauchy" = dcauchy,
                   "chi-squared" = dchisq,
                   "exponential" = dexp,
                   "f" = df,
                   "gamma" = dgamma,
                   "geometric" = dgeom,
                   "log-normal" = dlnorm,
                   "lognormal" = dlnorm,
                   "logistic" = dlogis,
                   "negative binomial" = dnbinom,
                   "normal" = dnorm,
                   "poisson" = dpois,
                   "t" = mydt,
                   "weibull" = dweibull,
                   NULL)
        if(is.null(densfun)) stop("unsupported distribution")
        if(distname %in% c("lognormal",  "log-normal")) {
            if(!is.null(start))
                stop("supplying pars for the log-Normal is not supported")
            if(any(x <= 0))
                stop("need positive values to fit a log-Normal")
            lx <- log(x)
            sd0 <- sqrt((n-1)/n)*sd(lx)
            mx <- mean(lx)
            estimate <- c(mx, sd0)
            sds <- c(sd0/sqrt(n), sd0/sqrt(2*n))
            names(estimate) <- names(sds) <- c("meanlog", "sdlog")
            return(structure(list(estimate = estimate, sd = sds, n = n,
				  loglik = sum(dlnorm(x, mx, sd0, log=TRUE))),
                             class = "fitdistr"))
        }
        if(distname == "normal") {
            if(!is.null(start))
                stop("supplying pars for the Normal is not supported")
            sd0 <- sqrt((n-1)/n)*sd(x)
            mx <- mean(x)
            estimate <- c(mx, sd0)
            sds <- c(sd0/sqrt(n), sd0/sqrt(2*n))
            names(estimate) <- names(sds) <- c("mean", "sd")
            return(structure(list(estimate = estimate, sd = sds, n = n,
				  loglik = sum(dnorm(x, mx, sd0, log=TRUE))),
                             class = "fitdistr"))
        }
        if(distname == "poisson") {
            if(!is.null(start))
                stop("supplying pars for the Poisson is not supported")
            estimate <- mean(x)
            sds <- sqrt(estimate/n)
            names(estimate) <- names(sds) <- "lambda"
            return(structure(list(estimate = estimate, sd = sds, n = n,
				  loglik = sum(dpois(x, estimate, log=TRUE))),
                             class = "fitdistr"))
        }
        if(distname == "exponential") {
            if(!is.null(start))
                stop("supplying pars for the exponential is not supported")
            estimate <- 1/mean(x)
            sds <- estimate/sqrt(n)
            names(estimate) <- names(sds) <- "rate"
            return(structure(list(estimate = estimate, sd = sds, n = n,
				  loglik = sum(dexp(x, estimate, log=TRUE))),
                             class = "fitdistr"))
        }
        if(distname == "geometric") {
            if(!is.null(start))
                stop("supplying pars for the geometric is not supported")
            estimate <- 1/(1 + mean(x))
            sds <- estimate * sqrt((1-estimate)/n)
            names(estimate) <- names(sds) <- "prob"
            return(structure(list(estimate = estimate, sd = sds, n = n,
				  loglik = sum(dexp(x, estimate, log=TRUE))),
                             class = "fitdistr"))
        }
        if(distname == "weibull" && is.null(start)) {
            ## log-Weibull is Gumbel, so start from that
            m <- mean(log(x)); v <- var(log(x))
            shape <- 1.2/sqrt(v); scale <- exp(m + 0.572/shape)
            start <- list(shape = shape, scale = scale)
            start <- start[!is.element(names(start), dots)]
        }
        if(distname == "gamma" && is.null(start)) {
            m <- mean(x); v <- var(x)
            start <- list(shape = m^2/v, rate = m/v)
            start <- start[!is.element(names(start), dots)]
        }
        if(distname == "negative binomial" && is.null(start)) {
            m <- mean(x); v <- var(x)
            size <- if(v > m) m^2/(v - m) else 100
            start <- list(size = size, mu = m)
            start <- start[!is.element(names(start), dots)]
        }
        if(is.element(distname, c("cauchy", "logistic")) && is.null(start)) {
            start <- list(location = median(x), scale = IQR(x)/2)
            start <- start[!is.element(names(start), dots)]
        }
        if(distname == "t" && is.null(start)) {
            start <- list(m = median(x), s = IQR(x)/2, df = 10)
            start <- start[!is.element(names(start), dots)]
        }
    }
    if(is.null(start) || !is.list(start))
        stop("'start' must be a named list")
    nm <- names(start)
    ## reorder arguments to densfun
    f <- formals(densfun)
    args <- names(f)
    m <- match(nm, args)
    if(any(is.na(m)))
        stop("'start' specifies names which are not arguments to 'densfun'")
    formals(densfun) <- c(f[c(1, m)], f[-c(1, m)])
    dens <- function(parm, x, ...) densfun(x, parm, ...)
    if((l <- length(nm)) > 1L)
        body(dens) <-
            parse(text = paste("densfun(x,",
                  paste("parm[", 1L:l, "]", collapse = ", "),
                  ", ...)"))
    Call[[1L]] <- as.name("optim")
    Call$densfun <- Call$start <- NULL
    Call$x <- x # want local variable as eval in this frame
    Call$par <- start
    Call$fn <- if("log" %in% args) mylogfn else myfn
    Call$hessian <- TRUE
    if(is.null(Call$method)) {
        if(any(c("lower", "upper") %in% names(Call))) Call$method <- "L-BFGS-B"
        else if (length(start) > 1L) Call$method <- "BFGS"
        else Call$method <- "Nelder-Mead"
    }
    res <- eval.parent(Call)
    if(res$convergence > 0L) stop("optimization failed")
    sds <- sqrt(diag(solve(res$hessian)))
    structure(list(estimate = res$par, sd = sds,
                   loglik = - res$value, n = n), class = "fitdistr")
}

logLik.fitdistr <- function(object, REML = FALSE, ...)
{
    if (REML) stop("only 'REML = FALSE' is implemented")
    val <- object$loglik
    attr(val, "nobs") <- object$n
    attr(val, "df") <- length(object$estimate)
    class(val) <- "logLik"
    val
}


print.fitdistr <-
    function(x, digits = getOption("digits"), ...)
{
    ans <- format(rbind(x$estimate, x$sd), digits=digits)
    ans[1L, ] <- sapply(ans[1L, ], function(x) paste("", x))
    ans[2L, ] <- sapply(ans[2L, ], function(x) paste("(", x, ")", sep=""))
    ## only used for digits
    dn <- dimnames(ans)
    dn[[1L]] <- rep("", 2L)
    dn[[2L]] <- paste(substring("      ", 1L, (nchar(ans[2L,]) - nchar(dn[[2L]])) %/% 2), dn[[2L]])
    dn[[2L]] <- paste(dn[[2L]], substring("      ", 1L, (nchar(ans[2L,]) - nchar(dn[[2L]])) %/% 2))
    dimnames(ans) <- dn
    print(ans, quote = FALSE)
    invisible(x)
}

coef.fitdistr <- function(object, ...) object$estimate

