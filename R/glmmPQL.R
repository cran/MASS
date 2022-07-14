# file MASS/R/glmmPQL.R
# copyright (C) 2002-2022 W. N. Venables and B. D. Ripley
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
glmmPQL <- function(fixed, random, family, data, correlation, weights,
                    control, niter = 10, verbose = TRUE, ...)
{
    if(!requireNamespace("nlme", quietly = TRUE))
        stop("package 'nlme' is essential")
    ## family
    if(is.character(family)) family <- get(family)
    if(is.function(family)) family <- family()
    if(is.null(family$family)) {
	print(family)
	stop("'family' not recognized")
    }
    m <- mcall <- Call <- match.call()
    nm <- names(m)[-1L]
    keep <- is.element(nm, c("weights", "data", "subset", "na.action"))
    for(i in nm[!keep]) m[[i]] <- NULL
    allvars <-
        if (is.list(random))
            c(all.vars(fixed), names(random),
              unlist(lapply(random, function(x) all.vars(formula(x)))))
        else c(all.vars(fixed), all.vars(random))
    ## allvars does not contain offset term.
    Terms <- if(missing(data)) terms(fixed) else terms(fixed, data = data)
    offt <- attr(Terms, "offset")
    have_offset <- length(offt) >= 1L
    if(have_offset) {
        offvars <- as.character(attr(Terms, "variables"))[offt + 1L]
        allvars <- c(allvars, offvars)
    }
    ## add variables in correlation argument, if any
    if (!missing(correlation) && !is.null(attr(correlation,"formula")))
        allvars <- c(allvars, all.vars(attr(correlation,"formula")))
    ## substitute back actual formula (rather than a variable name)
    Call$fixed <- eval(fixed); Call$random <- eval(random)
    m$formula <- as.formula(paste("~", paste(allvars, collapse="+")))
    environment(m$formula) <- environment(fixed)
    m$drop.unused.levels <- TRUE
    m[[1L]] <- quote(stats::model.frame)
    mf <- eval.parent(m)
    off <- model.offset(mf)
    if(is.null(off)) off <- 0
    ## as from 2022-04, nlme refuses fixed with offsets,
    ## rather than ignoring them.
    fixed2 <-
        if(have_offset) {
            tf <- drop.terms(Terms, offt, keep.response = TRUE)
            reformulate(attr(tf,"term.labels"), response = fixed[[2L]],
                        intercept = attr(tf, "intercept"),
                        env = environment(fixed))
        } else fixed
    wts <-  model.weights(mf)
    if(is.null(wts)) wts <- rep(1, nrow(mf))
    mf$wts <- wts
    fit0 <- glm(formula = fixed, family = family, data = mf, weights = wts, ...)
    w <- fit0$prior.weights
    eta <- fit0$linear.predictors
    zz <- eta + fit0$residuals - off
    wz <- fit0$weights
    fam <- family

    nm <- names(mcall)[-1L]
    keep <- is.element(nm, c("fixed", "random", "data", "subset",
                             "na.action", "control"))
    for(i in nm[!keep]) mcall[[i]] <- NULL
    fixed2[[2L]] <- quote(zz)
    mcall[["fixed"]] <- fixed2
    mcall[[1L]] <- quote(nlme::lme)
    mcall$random <- random
    mcall$method <- "ML"
    if(!missing(correlation))
        mcall$correlation <- correlation
    mcall$weights <- quote(nlme::varFixed(~invwt))
    mf$zz <- zz
    mf$invwt <- 1/wz
    mcall$data <- mf

    for(i in seq_len(niter)) {
        if(verbose) message(gettextf("iteration %d", i), domain = NA)
        fit <- eval(mcall)
        etaold <- eta
        ## update zz and invwt
        eta <- fitted(fit) + off
        if(sum((eta-etaold)^2) < 1e-6*sum(eta^2)) break;
        mu <- fam$linkinv(eta)
        mu.eta.val <- fam$mu.eta(eta)
        mf$zz <- eta + (fit0$y - mu)/mu.eta.val - off
        wz <- w * mu.eta.val^2 / fam$variance(mu)
        mf$invwt <- 1/wz
        mcall$data <- mf
    }
    attributes(fit$logLik) <- NULL # needed for some versions of nlme
    fit$call <- Call
    fit$family <- family
    fit$logLik <- as.numeric(NA)
    if(have_offset) {
        ## fit$fitted is a matrix
        fit$fitted <- fit$fitted + off # will recycle
        fit$have_offset <- TRUE
        fit$fixed2 <- fixed2
        fit$offvars <- offvars
    }
    oldClass(fit) <- c("glmmPQL", oldClass(fit))
    fit
}

predict.glmmPQL <-
  function(object, newdata = NULL, type = c("link", "response"),
           level = Q, na.action = na.pass, ...)
{
    type <- match.arg(type)
    Q <- object$dims$Q # default for level
    if(missing(newdata)) {
        pred <- fitted(object, level = level)
        pred <- switch(type,
                       link = pred,
                       response = object$family$linkinv(pred))
        if(!is.null(object$na.action))
            pred <- napredict(object$na.action, pred)
    } else {
        have_offset <- isTRUE(object$have_offset)
        class(object) <- class(object)[-1L] # call lme method
        if(have_offset) {
            ## remove offset(s) from call
            mcall <- object$call
            mcall[["fixed"]] <- object$fixed2
            object$call <- mcall
            pred <- predict(object, newdata, level = level,
                            na.action = na.action)
            ## Now add back offset(s)
            offvars <- sub("offset[(]", "", object$offvars)
            offvars <- sub("[)]", "", offvars)
            for (n in offvars) pred <- pred + eval(as.name(n), newdata)
        } else
            pred <- predict(object, newdata, level = level,
                            na.action = na.action)
        switch(type,
               response = {pred <- object$family$linkinv(pred)},
               link =)
    }
    pred
}

logLik.glmmPQL <- function (object, ...)
{
    p <- object$dims$ncol[object$dims$Q + 1]
    N <- object$dims$N
    val <- as.numeric(NA)
    attr(val, "nall") <- N
    attr(val, "nobs") <- N
    attr(val, "df") <- p + length(coef(object[["modelStruct"]])) + 1
    class(val) <- "logLik"
    val
}

anova.glmmPQL <- function (object, ...)
    stop("'anova' is not available for PQL fits")

## Suggested by Sebastian Meyer 2022-03

formula.glmmPQL <- function(x, ...) eval(x$call$fixed)
