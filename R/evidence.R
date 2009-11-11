evidence <-
  function(aic.table, model.high="top", model.low) {
    if(!identical(class(aic.table)[1], "aictab")) {stop("The input object must be of class 'aic.tab'")}

    if(identical(model.high, "top")) {
      mod1.AICwt <- aic.table[1,6]
    } else {mod1.AICwt <- aic.table[which(aic.table$Modnames==paste(model.high)),6]}

    mod2.AICwt <- aic.table[which(aic.table$Modnames==paste(model.low)),6]
    ev.ratio <- mod1.AICwt/mod2.AICwt
    ev.ratio.list <- list("Model.high" = paste(model.high), "Model.low" = paste(model.low), "Ev.ratio" = ev.ratio)
    class(ev.ratio.list) <- c("evidence", "list")
    return(ev.ratio.list)
  }

print.evidence <- function(x, digits = 2, ...) {
  cat("\nEvidence ratio between models '", x$Model.high,"' and '", x$Model.low, "':\n")
  cat(round(x$Ev.ratio, digits=digits), "\n\n")
}
