##generic
checkConv <- function(mod, ...) {
    UseMethod("checkConv", mod)
}



##default
checkConv.default <- function(mod, ...) {
    stop("\nFunction not yet defined for this object class\n")
}



##betareg
checkConv.betareg <- function(mod, ...) {
    if(mod$optim$convergence == 0) {
        conv <- TRUE
    } else {conv <- FALSE}

    msg <- mod$optim$message
    out <- list(converged = conv, message = msg)
    class(out) <- "checkConv"
    return(out)
}



##clm
checkConv.clm <- function(mod, ...) {
    if(mod$convergence$code == 0) {
        conv <- TRUE
    } else {conv <- FALSE}

    msg <- mod$convergence$alg.message
    out <- list(converged = conv, message = msg)
    class(out) <- "checkConv"
    return(out)
}



##clmm
checkConv.clmm <- function(mod, ...) {
    if(mod$optRes$convergence == 0) {
        conv <- TRUE
    } else {conv <- FALSE}

    msg <- mod$optRes$message
    out <- list(converged = conv, message = msg)
    class(out) <- "checkConv"
    return(out)
}



##glm
checkConv.glm <- function(mod, ...) {
    if(mod$converged) {
        conv <- TRUE
    } else {conv <- FALSE}

    msg <- NULL ##object does not include a message from IWLS algorithm
    out <- list(converged = conv, message = msg)
    class(out) <- "checkConv"
    return(out)
}



##glmmTMB
checkConv.glmmTMB <- function(mod, ...) {
    if(mod$fit$convergence == 0) {
        conv <- TRUE
    } else {conv <- FALSE}

    msg <- mod$fit$message ##object does not include a message from IWLS algorithm
    out <- list(converged = conv, message = msg)
    class(out) <- "checkConv"
    return(out)
}



##hurdle
checkConv.hurdle <- function(mod, ...) {
    if(mod$converged) {
        conv <- TRUE
    } else {conv <- FALSE}

    msg <- NULL ##object does not include a message from IWLS algorithm
    out <- list(converged = conv, message = msg)
    class(out) <- "checkConv"
    return(out)
}



##lavaan
checkConv.lavaan <- function(mod, ...) {
    if(mod@Fit@converged) {
        conv <- TRUE
    } else {conv <- FALSE}

    msg <- NULL ##object does not include a message from IWLS algorithm
    out <- list(converged = conv, message = msg)
    class(out) <- "checkConv"
    return(out)
}



##maxlikefit
checkConv.maxlikeFit <- function(mod, ...) {
    if(mod$optim$convergence == 0) {
        conv <- TRUE
    } else {conv <- FALSE}

    msg <- mod$optim$message
    out <- list(converged = conv, message = msg)
    class(out) <- "checkConv"
    return(out)
}



##merMod
checkConv.merMod <- function(mod, ...) {
    if(mod@optinfo$conv$opt == 0) {
        conv <- TRUE
    } else {conv <- FALSE}

    msg <- NULL ##object does not include a message from IWLS algorithm
    out <- list(converged = conv, message = msg)
    class(out) <- "checkConv"
    return(out)
}



##lmerModLmerTest
checkConv.lmerModLmerTest <- function(mod, ...) {
    if(mod@optinfo$conv$opt == 0) {
        conv <- TRUE
    } else {conv <- FALSE}

    msg <- NULL ##object does not include a message from IWLS algorithm
    out <- list(converged = conv, message = msg)
    class(out) <- "checkConv"
    return(out)
}



##multinom
checkConv.multinom <- function(mod, ...) {
    if(mod$convergence == 0) {
        conv <- TRUE
    } else {conv <- FALSE}

    msg <- NULL ##object does not include a message from IWLS algorithm
    out <- list(converged = conv, message = msg)
    class(out) <- "checkConv"
    return(out)
}



##nls
checkConv.nls <- function(mod, ...) {
    if(mod$convInfo$isConv) {
        conv <- TRUE
    } else {conv <- FALSE}

    msg <- mod$convInfo$stopMessage
    out <- list(converged = conv, message = msg)
    class(out) <- "checkConv"
    return(out)
}



##polr
checkConv.polr <- function(mod, ...) {
    if(mod$convergence == 0) {
        conv <- TRUE
    } else {conv <- FALSE}

    msg <- NULL ##object does not include a message from IWLS algorithm
    out <- list(converged = conv, message = msg)
    class(out) <- "checkConv"
    return(out)
}



##unmarked
checkConv.unmarkedFit <- function(mod, ...) {
    if(mod@opt$convergence == 0) {
        conv <- TRUE
    } else {conv <- FALSE}

    msg <- mod@opt$message
    out <- list(converged = conv, message = msg)
    class(out) <- "checkConv"
    return(out)
}



##zeroinfl
checkConv.zeroinfl <- function(mod, ...) {
    if(mod$converged) {
        conv <- TRUE
    } else {conv <- FALSE}

    msg <- NULL ##object does not include a message from IWLS algorithm
    out <- list(converged = conv, message = msg)
    class(out) <- "checkConv"
    return(out)
}



print.checkConv <- function(x, ...) {
  cat("\nConverged: ", x$converged, "\n")
  if(!is.null(x$message)) {
      cat("(", x$message, ")", "\n", sep = "")
  }
}
