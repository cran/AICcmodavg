##Goodness-of-fit test based on the chi-square
##chi-square
Nmix.chisq <- function(mod, ...) {
  UseMethod("Nmix.chisq", mod)
}


Nmix.chisq.default <- function(mod, ...) {
  stop("\nFunction not yet defined for this object class\n")
} 


##PCount
Nmix.chisq.unmarkedFitPCount <- function(mod, ...) {

  ##extract original data from model object
  obs <- mod@data@y
  ##extract fitted values
  fits <- fitted(mod)
  ##check if sites were removed from analysis
  sr <- mod@sitesRemoved
  if(length(sr) > 0) {
    obs <- obs[-sr, ]
    fits <- fits[-sr, ]
  }
  
  ##add NA's where fitted values are NA
  #obs[is.na(fits)] <- NA
  ##compute chi-square
  chi.sq <- sum((obs - fits)^2/fits, na.rm = TRUE) #added argument na.rm = TRUE when NA's occur
  result <- list(chi.square = chi.sq, model.type = class(mod)[1])
  class(result) <- "Nmix.chisq"
  return(result)
}


##PCO
Nmix.chisq.unmarkedFitPCO <- function(mod, ...) {

  ##extract original data from model object
  obs <- mod@data@y
  ##extract fitted values
  fits <- fitted(mod)
  ##check if sites were removed from analysis
  sr <- mod@sitesRemoved
  if(length(sr) > 0) {
    obs <- obs[-sr, ]
    fits <- fits[-sr, ]
  }
  
  ##add NA's where fitted values are NA
  #obs[is.na(fits)] <- NA
  ##compute chi-square
  chi.sq <- sum((obs - fits)^2/fits, na.rm = TRUE) #added argument na.rm = TRUE when NA's occur
  result <- list(chi.square = chi.sq, model.type = class(mod)[1])
  class(result) <- "Nmix.chisq"
  return(result)
}


##unmarkedFitMPois
Nmix.chisq.unmarkedFitMPois <- function(mod, ...) {

  ##extract original data from model object
  obs <- mod@data@y
  ##extract fitted values
  fits <- fitted(mod)
  ##check if sites were removed from analysis
  sr <- mod@sitesRemoved
  if(length(sr) > 0) {
    obs <- obs[-sr, ]
    fits <- fits[-sr, ]
  }
  
  ##add NA's where fitted values are NA
  #obs[is.na(fits)] <- NA
  ##compute chi-square
  chi.sq <- sum((obs - fits)^2/fits, na.rm = TRUE) #added argument na.rm = TRUE when NA's occur
  result <- list(chi.square = chi.sq, model.type = class(mod)[1])
  class(result) <- "Nmix.chisq"
  return(result)
}


##unmarkedFitDS
Nmix.chisq.unmarkedFitDS <- function(mod, ...) {

  ##extract original data from model object
  obs <- mod@data@y
  ##extract fitted values
  fits <- fitted(mod)
  ##check if sites were removed from analysis
  sr <- mod@sitesRemoved
  if(length(sr) > 0) {
    obs <- obs[-sr, ]
    fits <- fits[-sr, ]
  }
  
  ##add NA's where fitted values are NA
  #obs[is.na(fits)] <- NA
  ##compute chi-square
  chi.sq <- sum((obs - fits)^2/fits, na.rm = TRUE) #added argument na.rm = TRUE when NA's occur
  result <- list(chi.square = chi.sq, model.type = class(mod)[1])
  class(result) <- "Nmix.chisq"
  return(result)
}


##unmarkedFitGDS
Nmix.chisq.unmarkedFitGDS <- function(mod, ...) {

  ##extract original data from model object
  obs <- mod@data@y
  ##extract fitted values
  fits <- fitted(mod)
  ##check if sites were removed from analysis
  sr <- mod@sitesRemoved
  if(length(sr) > 0) {
    obs <- obs[-sr, ]
    fits <- fits[-sr, ]
  }
  
  ##add NA's where fitted values are NA
  #obs[is.na(fits)] <- NA
  ##compute chi-square
  chi.sq <- sum((obs - fits)^2/fits, na.rm = TRUE) #added argument na.rm = TRUE when NA's occur
  result <- list(chi.square = chi.sq, model.type = class(mod)[1])
  class(result) <- "Nmix.chisq"
  return(result)
}


##unmarkedFitGPC
Nmix.chisq.unmarkedFitGPC <- function(mod, ...) {

  ##extract original data from model object
  obs <- mod@data@y
  ##extract fitted values
  fits <- fitted(mod)
  ##check if sites were removed from analysis
  sr <- mod@sitesRemoved
  if(length(sr) > 0) {
    obs <- obs[-sr, ]
    fits <- fits[-sr, ]
  }
  
  ##add NA's where fitted values are NA
  #obs[is.na(fits)] <- NA
  ##compute chi-square
  chi.sq <- sum((obs - fits)^2/fits, na.rm = TRUE) #added argument na.rm = TRUE when NA's occur
  result <- list(chi.square = chi.sq, model.type = class(mod)[1])
  class(result) <- "Nmix.chisq"
  return(result)
}


##unmarkedFitGMM
Nmix.chisq.unmarkedFitGMM <- function(mod, ...) {

  ##extract original data from model object
  obs <- mod@data@y
  ##extract fitted values
  fits <- fitted(mod)
  ##check if sites were removed from analysis
  sr <- mod@sitesRemoved
  if(length(sr) > 0) {
    obs <- obs[-sr, ]
    fits <- fits[-sr, ]
  }
  
  ##add NA's where fitted values are NA
  #obs[is.na(fits)] <- NA
  ##compute chi-square
  chi.sq <- sum((obs - fits)^2/fits, na.rm = TRUE) #added argument na.rm = TRUE when NA's occur
  result <- list(chi.square = chi.sq, model.type = class(mod)[1])
  class(result) <- "Nmix.chisq"
  return(result)
}



##multmixOpen
Nmix.chisq.unmarkedFitMMO <- function(mod, ...) {

  ##extract original data from model object
  obs <- mod@data@y
  ##extract fitted values
  fits <- fitted(mod)
  ##check if sites were removed from analysis
  sr <- mod@sitesRemoved
  if(length(sr) > 0) {
    obs <- obs[-sr, ]
    fits <- fits[-sr, ]
  }
  
  ##add NA's where fitted values are NA
  #obs[is.na(fits)] <- NA
  ##compute chi-square
  chi.sq <- sum((obs - fits)^2/fits, na.rm = TRUE) #added argument na.rm = TRUE when NA's occur
  result <- list(chi.square = chi.sq, model.type = class(mod)[1])
  class(result) <- "Nmix.chisq"
  return(result)
}



##distsampOpen
Nmix.chisq.unmarkedFitDSO <- function(mod, ...) {

  ##extract original data from model object
  obs <- mod@data@y
  ##extract fitted values
  fits <- fitted(mod)
  ##check if sites were removed from analysis
  sr <- mod@sitesRemoved
  if(length(sr) > 0) {
    obs <- obs[-sr, ]
    fits <- fits[-sr, ]
  }
  
  ##add NA's where fitted values are NA
  #obs[is.na(fits)] <- NA
  ##compute chi-square
  chi.sq <- sum((obs - fits)^2/fits, na.rm = TRUE) #added argument na.rm = TRUE when NA's occur
  result <- list(chi.square = chi.sq, model.type = class(mod)[1])
  class(result) <- "Nmix.chisq"
  return(result)
}



##generic
Nmix.gof.test <- function(mod, nsim = 5, plot.hist = TRUE,
                          report = NULL, parallel = TRUE, ncores,
                          cex.axis = 1, cex.lab = 1, cex.main = 1,
                          lwd = 1,...){
  UseMethod("Nmix.gof.test", mod)
}


Nmix.gof.test.default <- function(mod, nsim = 5, plot.hist = TRUE,
                                  report = NULL, parallel = TRUE, ncores,
                                  cex.axis = 1, cex.lab = 1, cex.main = 1,
                                  lwd = 1,...){
  stop("\nFunction not yet defined for this object class\n")
}


##PCount
Nmix.gof.test.unmarkedFitPCount <- function(mod, nsim = 5, plot.hist = TRUE,
                                            report = NULL, parallel = TRUE, ncores,
                                            cex.axis = 1, cex.lab = 1, cex.main = 1,
                                            lwd = 1, ...){
  ##more bootstrap samples are recommended (e.g., 1000, 5000, or 10 000)

  ##extract model type
  model.type <- Nmix.chisq(mod)$model.type

  ##if NULL, don't print test statistic at each iteration
  if(is.null(report)) {
      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, parallel = parallel, ncores = ncores)
  } else {

      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, report = report, parallel = parallel, ncores = ncores)

  }
  
  ##determine significance
  p.value <- sum(out@t.star >= out@t0)/nsim
  if(p.value == 0) {
      p.display <- paste("<", round(1/nsim, digits = 4))
  } else {
      p.display = paste("=", round(p.value, digits = 4))
  }

  ##create plot
    if(plot.hist) {

        hist(out@t.star,
             main = as.expression(substitute("Bootstrapped "*chi^2*" fit statistic ("*nsim*" samples)",
                                             list(nsim = nsim))),
             xlim = range(c(out@t.star, out@t0)),
             xlab = paste("Simulated statistic ", "(observed = ", round(out@t0, digits = 2), ")", sep = ""),
             cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        title(main = bquote(paste(italic(P), .(p.display))), line = 0.5,
              cex.main = cex.main)
        abline(v = out@t0, lty = "dashed", col = "red", lwd = lwd)
    }
  
  ##estimate c-hat
  c.hat.est <- out@t0/mean(out@t.star)

  ##assemble result
  gof.out <- list(model.type = model.type, chi.square = out@t0, t.star = out@t.star,
                  p.value = p.value, c.hat.est = c.hat.est, nsim = nsim)
  class(gof.out) <- "Nmix.chisq"
  return(gof.out)  
}


##PCO
Nmix.gof.test.unmarkedFitPCO <- function(mod, nsim = 5, plot.hist = TRUE,
                                         report = NULL, parallel = TRUE, ncores,
                                         cex.axis = 1, cex.lab = 1, cex.main = 1,
                                         lwd = 1, ...){
  ##more bootstrap samples are recommended (e.g., 1000, 5000, or 10 000)

  ##extract model type
  model.type <- Nmix.chisq(mod)$model.type

  ##if NULL, don't print test statistic at each iteration
  if(is.null(report)) {
      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, parallel = parallel, ncores = ncores)
  } else {

      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, report = report, parallel = parallel,
                     ncores = ncores)
      
  }
  
  ##determine significance
  p.value <- sum(out@t.star >= out@t0)/nsim
  if(p.value == 0) {
      p.display <- paste("<", round(1/nsim, digits = 4))
  } else {
      p.display = paste("=", round(p.value, digits = 4))
  }

  ##create plot
    if(plot.hist) {
        hist(out@t.star,
             main = as.expression(substitute("Bootstrapped "*chi^2*" fit statistic ("*nsim*" samples)",
                                             list(nsim = nsim))),
             xlim = range(c(out@t.star, out@t0)),
             xlab = paste("Simulated statistic ", "(observed = ", round(out@t0, digits = 2), ")", sep = ""),
             cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        title(main = bquote(paste(italic(P), .(p.display))), line = 0.5,
              cex.main = cex.main)
        abline(v = out@t0, lty = "dashed", col = "red", lwd = lwd)
    }
  
  ##estimate c-hat
  c.hat.est <- out@t0/mean(out@t.star)

  ##assemble result
  gof.out <- list(model.type = model.type, chi.square = out@t0, t.star = out@t.star,
                  p.value = p.value, c.hat.est = c.hat.est, nsim = nsim)
  class(gof.out) <- "Nmix.chisq"
  return(gof.out)  
}


##DS
Nmix.gof.test.unmarkedFitDS <- function(mod, nsim = 5, plot.hist = TRUE,
                                        report = NULL, parallel = TRUE, ncores,
                                        cex.axis = 1, cex.lab = 1, cex.main = 1,
                                        lwd = 1, ...){
  ##more bootstrap samples are recommended (e.g., 1000, 5000, or 10 000)

  ##extract model type
  model.type <- Nmix.chisq(mod)$model.type

  ##if NULL, don't print test statistic at each iteration
  if(is.null(report)) {
      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, parallel = parallel, ncores = ncores)
  } else {

      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, report = report, parallel = parallel,
                     ncores = ncores)

  }

  ##determine significance
  p.value <- sum(out@t.star >= out@t0)/nsim
  if(p.value == 0) {
      p.display <- paste("<", round(1/nsim, digits = 4))
  } else {
      p.display = paste("=", round(p.value, digits = 4))
  }

  ##create plot
    if(plot.hist) {

        hist(out@t.star,
             main = as.expression(substitute("Bootstrapped "*chi^2*" fit statistic ("*nsim*" samples)",
                                             list(nsim = nsim))),
             xlim = range(c(out@t.star, out@t0)),
             xlab = paste("Simulated statistic ", "(observed = ", round(out@t0, digits = 2), ")", sep = ""),
             cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        title(main = bquote(paste(italic(P), .(p.display))), line = 0.5,
              cex.main = cex.main)
        abline(v = out@t0, lty = "dashed", col = "red", lwd = lwd)
    }
  
  ##estimate c-hat
  c.hat.est <- out@t0/mean(out@t.star)

  ##assemble result
  gof.out <- list(model.type = model.type, chi.square = out@t0, t.star = out@t.star,
                  p.value = p.value, c.hat.est = c.hat.est, nsim = nsim)
  class(gof.out) <- "Nmix.chisq"
  return(gof.out)  
}


##GDS
Nmix.gof.test.unmarkedFitGDS <- function(mod, nsim = 5, plot.hist = TRUE,
                                         report = NULL, parallel = TRUE, ncores,
                                         cex.axis = 1, cex.lab = 1, cex.main = 1,
                                         lwd = 1, ...){
  ##more bootstrap samples are recommended (e.g., 1000, 5000, or 10 000)

  ##extract model type
  model.type <- Nmix.chisq(mod)$model.type

  ##if NULL, don't print test statistic at each iteration
  if(is.null(report)) {
      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, parallel = parallel, ncores = ncores)
  } else {

      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, report = report, parallel = parallel,
                     ncores = ncores)

  }
  
  ##determine significance
  p.value <- sum(out@t.star >= out@t0)/nsim
  if(p.value == 0) {
      p.display <- paste("<", round(1/nsim, digits = 4))
  } else {
      p.display = paste("=", round(p.value, digits = 4))
  }

  ##create plot
    if(plot.hist) {
        
        hist(out@t.star,
             main = as.expression(substitute("Bootstrapped "*chi^2*" fit statistic ("*nsim*" samples)",
                                             list(nsim = nsim))),
             xlim = range(c(out@t.star, out@t0)),
             xlab = paste("Simulated statistic ", "(observed = ", round(out@t0, digits = 2), ")", sep = ""),
             cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        title(main = bquote(paste(italic(P), .(p.display))), line = 0.5,
              cex.main = cex.main)
        abline(v = out@t0, lty = "dashed", col = "red", lwd = lwd)
    }
  
  ##estimate c-hat
  c.hat.est <- out@t0/mean(out@t.star)

  ##assemble result
  gof.out <- list(model.type = model.type, chi.square = out@t0, t.star = out@t.star,
                  p.value = p.value, c.hat.est = c.hat.est, nsim = nsim)
  class(gof.out) <- "Nmix.chisq"
  return(gof.out)  
}


##GMM
Nmix.gof.test.unmarkedFitGMM <- function(mod, nsim = 5, plot.hist = TRUE,
                                         report = NULL, parallel = TRUE, ncores,
                                         cex.axis = 1, cex.lab = 1, cex.main = 1,
                                         lwd = 1, ...){
  ##more bootstrap samples are recommended (e.g., 1000, 5000, or 10 000)

  ##extract model type
  model.type <- Nmix.chisq(mod)$model.type

  ##if NULL, don't print test statistic at each iteration
  if(is.null(report)) {
      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, parallel = parallel, ncores = ncores)
  } else {

      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, report = report, parallel = parallel,
                     ncores = ncores)

  }
  
  ##determine significance
  p.value <- sum(out@t.star >= out@t0)/nsim
  if(p.value == 0) {
      p.display <- paste("<", round(1/nsim, digits = 4))
  } else {
      p.display = paste("=", round(p.value, digits = 4))
  }

  ##create plot
    if(plot.hist) {
        
        hist(out@t.star,
             main = as.expression(substitute("Bootstrapped "*chi^2*" fit statistic ("*nsim*" samples)",
                                             list(nsim = nsim))),
             xlim = range(c(out@t.star, out@t0)),
             xlab = paste("Simulated statistic ", "(observed = ", round(out@t0, digits = 2), ")", sep = ""),
             cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        title(main = bquote(paste(italic(P), .(p.display))), line = 0.5,
              cex.main = cex.main)
        abline(v = out@t0, lty = "dashed", col = "red", lwd = lwd)
    }
  
  ##estimate c-hat
  c.hat.est <- out@t0/mean(out@t.star)

  ##assemble result
  gof.out <- list(model.type = model.type, chi.square = out@t0, t.star = out@t.star,
                  p.value = p.value, c.hat.est = c.hat.est, nsim = nsim)
  class(gof.out) <- "Nmix.chisq"
  return(gof.out)  
}


##GPC
Nmix.gof.test.unmarkedFitGPC <- function(mod, nsim = 5, plot.hist = TRUE,
                                         report = NULL, parallel = TRUE, ncores,
                                         cex.axis = 1, cex.lab = 1, cex.main = 1,
                                         lwd = 1,...){
  ##more bootstrap samples are recommended (e.g., 1000, 5000, or 10 000)

  ##extract model type
  model.type <- Nmix.chisq(mod)$model.type

  ##if NULL, don't print test statistic at each iteration
  if(is.null(report)) {
      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, parallel = parallel, ncores = ncores)
  } else {

      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, report = report, parallel = parallel,
                     ncores = ncores)

  }
  
  ##determine significance
  p.value <- sum(out@t.star >= out@t0)/nsim
  if(p.value == 0) {
      p.display <- paste("<", round(1/nsim, digits = 4))
  } else {
      p.display = paste("=", round(p.value, digits = 4))
  }

  ##create plot
    if(plot.hist) {
        
        hist(out@t.star,
             main = as.expression(substitute("Bootstrapped "*chi^2*" fit statistic ("*nsim*" samples)",
                                             list(nsim = nsim))),
             xlim = range(c(out@t.star, out@t0)),
             xlab = paste("Simulated statistic ", "(observed = ", round(out@t0, digits = 2), ")", sep = ""),
             cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        title(main = bquote(paste(italic(P), .(p.display))), line = 0.5,
              cex.main = cex.main)
        abline(v = out@t0, lty = "dashed", col = "red", lwd = lwd)
    }
  
  ##estimate c-hat
  c.hat.est <- out@t0/mean(out@t.star)

  ##assemble result
  gof.out <- list(model.type = model.type, chi.square = out@t0, t.star = out@t.star,
                  p.value = p.value, c.hat.est = c.hat.est, nsim = nsim)
  class(gof.out) <- "Nmix.chisq"
  return(gof.out)  
}


##MPois
Nmix.gof.test.unmarkedFitMPois <- function(mod, nsim = 5, plot.hist = TRUE,
                                           report = NULL, parallel = TRUE, ncores,
                                           cex.axis = 1, cex.lab = 1, cex.main = 1,
                                           lwd = 1, ...){
  ##more bootstrap samples are recommended (e.g., 1000, 5000, or 10 000)

  ##extract model type
  model.type <- Nmix.chisq(mod)$model.type

  ##if NULL, don't print test statistic at each iteration
  if(is.null(report)) {
      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, parallel = parallel, ncores = ncores)
  } else {

      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, report = report, parallel = parallel,
                     ncores = ncores)

  }
  
  ##determine significance
  p.value <- sum(out@t.star >= out@t0)/nsim
  if(p.value == 0) {
      p.display <- paste("<", round(1/nsim, digits = 4))
  } else {
      p.display = paste("=", round(p.value, digits = 4))
  }

  ##create plot
    if(plot.hist) {
              
        hist(out@t.star,
             main = as.expression(substitute("Bootstrapped "*chi^2*" fit statistic ("*nsim*" samples)",
                                             list(nsim = nsim))),
             xlim = range(c(out@t.star, out@t0)),
             xlab = paste("Simulated statistic ", "(observed = ", round(out@t0, digits = 2), ")", sep = ""),
             cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        title(main = bquote(paste(italic(P), .(p.display))), line = 0.5,
              cex.main = cex.main)
        abline(v = out@t0, lty = "dashed", col = "red", lwd = lwd)
    }
  
  ##estimate c-hat
  c.hat.est <- out@t0/mean(out@t.star)

  ##assemble result
  gof.out <- list(model.type = model.type, chi.square = out@t0, t.star = out@t.star,
                  p.value = p.value, c.hat.est = c.hat.est, nsim = nsim)
  class(gof.out) <- "Nmix.chisq"
  return(gof.out)  
}



##multmixOpen
Nmix.gof.test.unmarkedFitMMO <- function(mod, nsim = 5, plot.hist = TRUE,
                                         report = NULL, parallel = TRUE, ncores,
                                         cex.axis = 1, cex.lab = 1, cex.main = 1,
                                         lwd = 1, ...){
  ##more bootstrap samples are recommended (e.g., 1000, 5000, or 10 000)

  ##extract model type
  model.type <- Nmix.chisq(mod)$model.type

  ##if NULL, don't print test statistic at each iteration
  if(is.null(report)) {
      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, parallel = parallel, ncores = ncores)
  } else {

      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, report = report, parallel = parallel,
                     ncores = ncores)
      
  }
  
  ##determine significance
  p.value <- sum(out@t.star >= out@t0)/nsim
  if(p.value == 0) {
      p.display <- paste("<", round(1/nsim, digits = 4))
  } else {
      p.display = paste("=", round(p.value, digits = 4))
  }

  ##create plot
    if(plot.hist) {

        hist(out@t.star,
             main = as.expression(substitute("Bootstrapped "*chi^2*" fit statistic ("*nsim*" samples)",
                                             list(nsim = nsim))),
             xlim = range(c(out@t.star, out@t0)),
             xlab = paste("Simulated statistic ", "(observed = ", round(out@t0, digits = 2), ")", sep = ""),
             cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        title(main = bquote(paste(italic(P), .(p.display))), line = 0.5,
              cex.main = cex.main)
        abline(v = out@t0, lty = "dashed", col = "red", lwd = lwd)
    }
  
  ##estimate c-hat
  c.hat.est <- out@t0/mean(out@t.star)

  ##assemble result
  gof.out <- list(model.type = model.type, chi.square = out@t0, t.star = out@t.star,
                  p.value = p.value, c.hat.est = c.hat.est, nsim = nsim)
  class(gof.out) <- "Nmix.chisq"
  return(gof.out)  
}



##distsampOpen
Nmix.gof.test.unmarkedFitDSO <- function(mod, nsim = 5, plot.hist = TRUE,
                                         report = NULL, parallel = TRUE, ncores,
                                         cex.axis = 1, cex.lab = 1, cex.main = 1,
                                         lwd = 1, ...){
  ##more bootstrap samples are recommended (e.g., 1000, 5000, or 10 000)

  ##extract model type
  model.type <- Nmix.chisq(mod)$model.type

  ##if NULL, don't print test statistic at each iteration
  if(is.null(report)) {
      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, parallel = parallel, ncores = ncores)
  } else {

      ##compute GOF P-value
      out <- parboot(mod, statistic = function(i) Nmix.chisq(i)$chi.square,
                     nsim = nsim, report = report, parallel = parallel,
                     ncores = ncores)
      
  }
  
  ##determine significance
  p.value <- sum(out@t.star >= out@t0)/nsim
  if(p.value == 0) {
      p.display <- paste("<", round(1/nsim, digits = 4))
  } else {
      p.display = paste("=", round(p.value, digits = 4))
  }

  ##create plot
    if(plot.hist) {
  
        hist(out@t.star,
             main = as.expression(substitute("Bootstrapped "*chi^2*" fit statistic ("*nsim*" samples)",
                                             list(nsim = nsim))),
             xlim = range(c(out@t.star, out@t0)),
             xlab = paste("Simulated statistic ", "(observed = ", round(out@t0, digits = 2), ")", sep = ""),
             cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main)
        title(main = bquote(paste(italic(P), .(p.display))), line = 0.5,
               cex.main = cex.main)
        abline(v = out@t0, lty = "dashed", col = "red", lwd = lwd)
    }
  
  ##estimate c-hat
  c.hat.est <- out@t0/mean(out@t.star)

  ##assemble result
  gof.out <- list(model.type = model.type, chi.square = out@t0, t.star = out@t.star,
                  p.value = p.value, c.hat.est = c.hat.est, nsim = nsim)
  class(gof.out) <- "Nmix.chisq"
  return(gof.out)  
}



##print method
print.Nmix.chisq <- function(x, digits.vals = 2, digits.chisq = 4, ...) {
  cat("\nChi-square goodness-of-fit for N-mixture model of \'", x$model.type, "\' class\n", sep = "")
  cat("\nObserved chi-square statistic =", round(x$chi.square, digits = digits.chisq), "\n")
  if(length(x) > 2){
    cat("Number of bootstrap samples =", x$nsim)
    cat("\nP-value =", x$p.value)
    cat("\n\nQuantiles of bootstrapped statistics:\n")
    print(quantile(x$t.star), digits = digits.vals)
    cat("\nEstimate of c-hat =", round(x$c.hat.est, digits = digits.vals), "\n")
  }
  cat("\n")
}
