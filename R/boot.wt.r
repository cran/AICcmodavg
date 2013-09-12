boot.wt <- function(cand.set, modnames, sort = TRUE, c.hat = 1,
                  second.ord = TRUE, nobs = NULL, nsim = 100){
  ##nsim Burnham and Anderson 2002, suggest 1000 or more iterations

  ##extract data from first model
  ##aictab( ) already checks that data sets are identical for all models
  data <- eval(cand.set[[1]]$call$data, environment(formula(cand.set[[1]])))

  ##check for model class
  ##OK for lm, glm, multinom, polr, rlm
  results <- NULL
  known <- 0
  ##extract classes
  mod.class <- unlist(lapply(X = cand.set, FUN = class))
  ##check if all are identical
  check.class <- unique(mod.class)

  ##determine if lm or glm
  if(identical(check.class, "lm") || identical(check.class, c("glm", "lm"))) {
    known <- 1
  }   

  ##determine if multinom
  if(identical(sort(check.class), c("multinom", "nnet"))) {
    known <- 1
  }   

  ##determine if polr
  if(identical(check.class, "polr")) {
    known <- 1
  }   
      
  ##determine if rlm
  if(identical(check.class, c("rlm", "lm")))  {
    known <- 1
  }
  
  ##warn if class is neither lm, glm, multinom, polr
  if(sum(known) < 1) {stop("\nFunction not yet defined for this object class\n")}

  
  ##create vector to store top model
  top <- character(nsim)
  
  for(i in 1:nsim) {
    ##resample data
    new.data <- data[sample(x = rownames(data), replace = TRUE), ]

    ##run models based on contents of list
    ##extract formulas
    ##model.forms <- lapply(X = cand.set, FUN = formula)

    ##to get complete call
    ##model.calls <- lapply(X = cand.set, FUN = getCall)
    
    ##store results
    results <- lapply(X = cand.set, FUN = function(j) update(j, data = new.data))

    ##compute AIC table
    out <- aictab(cand.set = results, modnames = modnames, sort = TRUE, second.ord = second.ord,
                  c.hat = c.hat, nobs = nobs)

    ##determine highest-ranked model
    top[i] <- as.character(out$Modnames[1])
  }

  ##compute selection frequencies for each model
  rel.freqs <- table(top)/nsim

  
  ##check whether all models appear in table
  if(length(cand.set) != length(rel.freqs)) {
    ##assign 0 freqs for models never appearing at first rank
    all.freqs <- rep(0, times = length(modnames))
    names(all.freqs) <- modnames
    ##iterate over observed models
    for (k in 1:length(rel.freqs)){
      all.freqs[names(rel.freqs)[k]] <- rel.freqs[k]
    }
    rel.freqs <- all.freqs
  }
    

  ##original model selection
  orig.aic <- aictab(cand.set = cand.set, modnames = modnames, second.ord = second.ord,
                     sort = FALSE, c.hat = c.hat, nobs = nobs)

  ##rename column for relative frequencies
  names(orig.aic)[7] <- "PiWt"
  
  ##add column for pi_i
  orig.aic$PiWt<- rel.freqs

  ##reorder if required
  if (sort) {
    orig.aic <- orig.aic[order(orig.aic[, 3]), ]
    #orig.aic$Cum.PiWt <- cumsum(orig.aic[, 7])
  }
  
  class(orig.aic) <- c("aictab", "data.frame")
  return(orig.aic)
}
