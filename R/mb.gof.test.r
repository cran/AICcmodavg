##MacKenzie and Bailey goodness of fit test for single season occupancy models
mb.chisq <- function(mod, print.table = TRUE) {

##for single-season occupancy models
  if(!identical(class(mod)[1], "unmarkedFitOccu")) {
    stop("\nThis function is only appropriate for single-season occupancy models\n")
  }
  
##step 1:
##extract detection histories
y.raw <- getData(mod)@y

##if some rows are all NA and sites are discarded, adjust sample size accordingly
N.raw <- nrow(y.raw)
#if(all NA) {N - number of rows with all NA}
##identify sites without data
na.raw <- apply(X = y.raw, MARGIN = 1, FUN = function(i) all(is.na(i)))
  
##remove sites without data
y.data <- y.raw[!na.raw, ]

N <- N.raw - sum(na.raw)

#N is required for computations in the end
T <- ncol(y.data)

det.hist <- apply(X = y.data, MARGIN = 1, FUN = function(i) paste(i, collapse = ""))

##compute predicted values of occupancy
preds.psi <- predict(mod, type = "state")$Predicted



##compute predicted values of p
preds.p <- matrix(data = predict(mod, type = "det")$Predicted,
                  ncol = T, byrow = TRUE)

##assemble in data.frame
out.hist <- data.frame(det.hist, preds.psi)

##identify unique histories
un.hist <- unique(det.hist)
n.un.hist <- length(un.hist)

##identify if missing values occur
na.vals <- length(grep(pattern = "NA", x = un.hist)) > 0

if(na.vals) {

  ##identify each history with NA
  id.na <- grep(pattern = "NA", x = un.hist)
  id.det.hist.na <- grep(pattern = "NA", x = det.hist)
  
  ##cohorts with NA
  cohort.na <- sort(un.hist[id.na])
  n.cohort.na <- length(cohort.na)
  ##determine cohorts that will be grouped together (same missing value)
  unique.na <- gsub(pattern = "NA", replacement = "N", x = cohort.na)
  ##determine which visit has missing value
  na.visits <- sapply(strsplit(x = unique.na, split = ""), FUN = function(i) paste(ifelse(i == "N", 1, 0), collapse = ""))
  ##add cohort labels for histories
  names(cohort.na) <- na.visits  
  ##number of histories in each cohort
  n.hist.missing.cohorts <- table(na.visits)
  ##number of missing cohorts
  n.missing.cohorts <- length(n.hist.missing.cohorts)
  out.hist.na <- out.hist[id.det.hist.na, ]
  out.hist.na$det.hist <- droplevels(out.hist.na$det.hist)

  ##groupings in out.hist.na
  just.na <- sapply(X = out.hist.na$det.hist, FUN = function(i) gsub(pattern = "1", replacement = "0", x = i))
  out.hist.na$coh <- sapply(X = just.na, FUN = function(i) gsub(pattern = "NA", replacement = "1", x = i))                                           
  ##number of sites in each missing cohort
  freqs.missing.cohorts <- table(out.hist.na$coh)
    
  ##number of sites with each history
  na.freqs <- table(det.hist[id.det.hist.na]) 
  preds.p.na <- preds.p[id.det.hist.na, ]
  
  ##cohorts without NA
  cohort.not.na <- sort(un.hist[-id.na])
  out.hist.not.na <- out.hist[-id.det.hist.na, ]
  out.hist.not.na$det.hist <- droplevels(out.hist.not.na$det.hist)
  n.cohort.not.na <- length(cohort.not.na)
  n.sites.not.na <- length(det.hist) - length(id.det.hist.na)
  preds.p.not.na <- preds.p[-id.det.hist.na, ]

  
} else {
  cohort.not.na <- sort(un.hist)
  out.hist.not.na <- out.hist
  preds.p.not.na <- preds.p
  n.cohort.not.na <- length(cohort.not.na)
  n.sites.not.na <- length(det.hist)
}
    
##for each missing data cohort, determine number of sites for each
  
  ##iterate over each site for each unique history
if(n.cohort.not.na > 0) {  ##expected frequencies for non-missing data
  exp.freqs <- rep(NA, n.cohort.not.na)
  names(exp.freqs) <- cohort.not.na ########################################SORT ENCOUNTER HISTORIES CHECK THAT ORDER IS IDENTICAL TO OBSERVED FREQS
  
  ##iterate over detection histories
  for (i in 1:n.cohort.not.na) {
    eq.solved <- rep(NA, n.sites.not.na)
    select.hist <- cohort.not.na[i]
    ##strip all values
    strip.hist <- unlist(strsplit(select.hist, split = ""))
    
    ##translate each visit in probability statement
    hist.mat <- matrix(NA, nrow = n.sites.not.na, ncol = T)

    ##iterate over sites
    for(j in 1:n.sites.not.na) {

      hist.mat[j, ] <- ifelse(strip.hist == "1", preds.p.not.na[j, ],
                              ifelse(strip.hist == "0", 1 - preds.p.not.na[j, ],
                                     0))
      
      ##combine into equation
      combo.p <- paste(hist.mat[j, ], collapse = "*")
      ##for history without detection
      if(sum(as.numeric(strip.hist)) == 0) {
        combo.first <- paste(c(out.hist.not.na[j, "preds.psi"], combo.p), collapse = "*")
        combo.psi.p <- paste((1 - out.hist.not.na[j, "preds.psi"]), "+", combo.first)
      } else {
        combo.psi.p <- paste(c(out.hist.not.na[j, "preds.psi"], combo.p), collapse = "*")       
      }
      eq.solved[j] <- eval(parse(text = as.expression(combo.psi.p)))
    }
    
    exp.freqs[i] <- sum(eq.solved, na.rm = TRUE)
  }


  ##for each detection history, compute observed frequencies
  freqs <- table(out.hist.not.na$det.hist)
  
  out.freqs <- matrix(NA, nrow = n.cohort.not.na, ncol = 4)
  colnames(out.freqs) <- c("Cohort", "Observed", "Expected", "Chi-square")
  rownames(out.freqs) <- names(freqs)
  ##cohort
  out.freqs[, 1] <- 0
  ##observed
  out.freqs[, 2] <- freqs
  ##expected
  out.freqs[, 3] <- exp.freqs
  ##chi-square
  out.freqs[, 4] <- ((out.freqs[, "Observed"] - out.freqs[, "Expected"])^2)/out.freqs[, "Expected"]
}

##if missing values
if(na.vals) {
##create list to store the chisquare for each cohort
  missing.cohorts <- list( )
  ##check if preds.p.na has only 1 row and change to matrix
  if(!is.matrix(preds.p.na)) {preds.p.na <- matrix(data = preds.p.na, nrow = 1)}
  
  for(m in 1:n.missing.cohorts) {
    ##select cohort
    select.cohort <- out.hist.na[which(out.hist.na$coh == names(freqs.missing.cohorts)[m]), ]
    select.preds.p.na <- preds.p.na[which(out.hist.na$coh == names(freqs.missing.cohorts)[m]), ]
    ##replace NA's with 1 to remove from likelihood
    if(!is.matrix(select.preds.p.na)) {select.preds.p.na <- matrix(data = select.preds.p.na, nrow = 1)}
    select.preds.p.na[, gregexpr(pattern = "N", text = gsub(pattern = "NA", replacement = "N", x = select.cohort$det.hist[1]))[[1]]] <- 1
    n.total.sites <- nrow(select.cohort)
    freqs.na <- table(droplevels(select.cohort$det.hist))
    cohort.na.un <- sort(unique(select.cohort$det.hist))
    n.hist.na <- length(freqs.na)
    exp.na <- rep(NA, n.hist.na)
    names(exp.na) <- cohort.na.un  

          for(i in 1:n.hist.na) {
            ##number of sites in given history
            n.sites.hist <- freqs.na[i] ##this should be number of sites for each history
            eq.solved <- rep(NA, n.total.sites)
            ##replace NA's with N
            select.hist <- gsub(pattern = "NA", replacement = "N", x = cohort.na.un[i])
            ##strip all values
            strip.hist <- unlist(strsplit(select.hist, split = ""))
            ##translate each visit in probability statement
            hist.mat <- matrix(NA, nrow = n.total.sites, ncol = T)
            ##iterate over sites
            for(j in 1:n.total.sites) {
              hist.mat[j, ] <- ifelse(strip.hist == "1", select.preds.p.na[j, ],
                                      ifelse(strip.hist == "0", 1 - select.preds.p.na[j, ], 1))
              ##replace NA by 1 (missing visit is removed from likelihood)
###################################################
###for missing value, remove occasion
###################################################      
              ##combine into equation
              combo.p <- paste(hist.mat[j, ], collapse = "*")
              ##for history without detection
              if(sum(as.numeric(gsub(pattern = "N", replacement = "0", x = strip.hist))) == 0) {
                combo.first <- paste(c(select.cohort[j, "preds.psi"], combo.p), collapse = "*")
                combo.psi.p <- paste((1 - select.cohort[j, "preds.psi"]), "+", combo.first)
              } else {
                combo.psi.p <- paste(c(select.cohort[j, "preds.psi"], combo.p), collapse = "*")
              }  

              eq.solved[j] <- eval(parse(text = as.expression(combo.psi.p)))
            }
            
            exp.na[i] <- sum(eq.solved, na.rm = TRUE)
          }

    ##compute chisq for missing data cohorts
    
    ##for each detection history, compute observed frequencies
    out.freqs.na <- matrix(NA, nrow = n.hist.na, ncol = 4)
    colnames(out.freqs.na) <- c("Cohort", "Observed", "Expected", "Chi-square")
    rownames(out.freqs.na) <- cohort.na.un
    ##cohort
    out.freqs.na[, 1] <- m
    ##observed
    out.freqs.na[, 2] <- freqs.na
    ##expected
    out.freqs.na[, 3] <- exp.na
    ##chi-square
    out.freqs.na[, 4] <- ((out.freqs.na[, "Observed"] - out.freqs.na[, "Expected"])^2)/out.freqs.na[, "Expected"]
    
    missing.cohorts[[m]] <- list(out.freqs.na = out.freqs.na)
  }
  
  
}
##test statistic is chi-square for all possible detection histories
##for observed detection histories, chisq = sum((obs - exp)^2)/exp = Y

##to avoid computing all possible detection histories, it is possible to obtain it by subtraction:
#for unobserved detection histories, chisq = sum((obs - exp)^2)/exp = sum(0 - exp)^2/exp = sum(exp) = X
##X = N.sites - sum(exp values from observed detection histories)
##Thus, test statistic = Y + X = Y + N - sum(exp values from observed detection histories)

##compute partial chi-square for observed detection histories (Y)
#chisq.obs.det <- sum(((out.freqs[, "observed"] - out.freqs[, "expected"])^2)/out.freqs[, "expected"])

##compute partial chi-square for unobserved detection histories (X)
if(na.vals) {
  chisq.missing <- do.call("rbind", lapply(missing.cohorts, FUN = function(i) i$out.freqs.na))
  if(n.cohort.not.na > 0) {
    chisq.unobs.det <- N - sum(out.freqs[, "Expected"]) - sum(chisq.missing[, "Expected"])
    chisq.table <- rbind(out.freqs, chisq.missing)
  } else {
    chisq.unobs.det <- N - sum(chisq.missing[, "Expected"])
    chisq.table <- chisq.missing
  }  
} else {
  chisq.unobs.det <- N - sum(out.freqs[, "Expected"])
  chisq.na <- 0
  chisq.table <- out.freqs
}

##test statistic (Y + X = Y - sum(exp observed detection histories)
chisq <- sum(chisq.table[, "Chi-square"]) + chisq.unobs.det

if(print.table) {
  out <- list(chisq.table = chisq.table, chi.square = chisq)
} else {
  out <- list(chi.square = chisq)
}

class(out) <- "mb.chisq"
return(out)
}



mb.gof.test <- function(mod, nsim = 5, plot.hist = TRUE){#more bootstrap samples are recommended (e.g., 1000, 5000, or 10 000)

  ##extract table from fitted model
  mod.table <- mb.chisq(mod)

  ##compute GOF P-value
  out <- parboot(mod, statistic = function(i) mb.chisq(i)$chi.square,
                 nsim = nsim)
  ##determine significance
  p.value <- sum(out@t.star >= out@t0)/nsim
  if(p.value == 0) {
    p.display <- paste("<", 1/nsim)
  } else {
    p.display = paste("=", round(p.value, digits = 4))
  }

  ##create plot
  if(plot.hist) {
  hist(out@t.star, main = paste("Bootstrapped MacKenzie and Bailey fit statistic (", nsim, " samples)", sep = ""),
       xlim = range(c(out@t.star, out@t0)), xlab = paste("Simulated statistic ", "(observed = ", round(out@t0, digits = 2), ")", sep = ""))
  title(main = bquote(paste(italic(P), " ", .(p.display))), line = 0.5)
  abline(v = out@t0, lty = "dashed", col = "red")
}
  ##estimate c-hat
  c.hat.est <- out@t0/mean(out@t.star)


##assemble result
gof.out <- list(chisq.table = mod.table$chisq.table, chi.square = mod.table$chi.square,
                t.star = out@t.star, p.value = p.value, c.hat.est = c.hat.est, nsim = nsim)
  class(gof.out) <- "mb.chisq"
  return(gof.out)  
}



print.mb.chisq <- function(x, digits.vals = 2, digits.chisq = 4, ...) {
  cat("\nMacKenzie and Bailey goodness-of-fit for single-season occupancy model\n")
  cat("\nPearson chi-square table:\n\n")
  print(round(x$chisq.table, digits = digits.vals))
  cat("\nChi-square statistic =", round(x$chi.square, digits = digits.chisq), "\n")
  if(length(x) > 2){
    cat("Number of bootstrap samples =", x$nsim)
    cat("\nP-value =", x$p.value)
    cat("\n\nQuantiles of bootstrapped statistics:\n")
    print(quantile(x$t.star), digits = digits.vals)
    cat("\nEstimate of c-hat =", round(x$c.hat.est, digits = digits.vals), "\n\n")
  }
}
