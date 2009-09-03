print.aictab <-
function(x, digits=4, ...) {
  cat("\nModel selection based on", colnames(x)[3], ":\n")
  if (ncol(x) > 7) {cat("(c-hat estimate = ", x$c_hat[1], ")\n")}
  cat("\n")
  nice.tab <- cbind(x[,2], x[,3], x[,4], x[,6], x[,7])
  colnames(nice.tab) <- colnames(x)[c(2,3,4,6,7)]
  rownames(nice.tab) <- x[,1]
  print(round(nice.tab, digits=digits)) #select rounding off with digits argument
  cat("\n")
}

