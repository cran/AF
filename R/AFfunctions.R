############## Common functions ##############
globalVariables(".SD")
is.binary <- function(v) {

  x <- unique(v)
  if((length(x) - sum(is.na(x)) == 2L) & (max(x) == 1 & min(x) == 0)) TRUE
  else FALSE

}

# The function "aggr" aggregate the data table x by clusterid
aggr <- function(x, clusters){

  temp <- data.table(x)
  temp <- as.matrix(temp[, j = lapply(.SD, sum), by = clusters])[, -1]

}

############## Summary and print functions ##############
#' @export
print.AF<-function(x, ...){
  if(!x$n.cluster == 0) {
    Std.Error <- "Robust SE"
    se <- "cluster-robust standard error"
  }
  else {
    Std.Error <- "Std.Error"
    se <- "standard error"
  }
  cat("\nEstimated attributable fraction (AF) and", se, ":", "\n")
  cat("\n")
  table.est <- cbind(x$AF.est, sqrt(x$AF.var))
  colnames(table.est) <- c("AF", Std.Error)
  r <- rep("", , length(x$AF.est))
  rownames(table.est) <- c(r)
  modelcall <- as.character(x$objectcall[1])
  if(modelcall == "coxph") {
    table.est <- cbind(x$times, table.est)
    colnames(table.est) <- c("Time", "AF", Std.Error)
    print.default(table.est)
  }
  else {
    print.default(table.est)
  }
}

CI.AF <- function(AF, Std.Error, confidence.level, CI.transform){
  if(CI.transform == "untransformed"){
    lower <- AF - abs(qnorm((1 - confidence.level) / 2)) * Std.Error
    upper <- AF + abs(qnorm((1 - confidence.level) / 2)) * Std.Error
  }
  if(CI.transform == "log"){
    lower <- AF * exp( - abs(qnorm((1 - confidence.level) / 2)) * Std.Error / AF)
    upper <- AF * exp(abs(qnorm((1 - confidence.level) / 2)) * Std.Error / AF)
  }
  if(CI.transform == "logit"){
    logit <- function(x) log(x / (1 - x))
    lower <- exp(logit(AF) - abs(qnorm((1 - confidence.level) / 2)) * Std.Error / (AF * (1 - AF))) / (1 + exp(logit(AF) - abs(qnorm((1 - confidence.level) / 2)) * Std.Error / (AF * (1 - AF))))
    upper <- exp(logit(AF) + abs(qnorm((1 - confidence.level) / 2)) * Std.Error / (AF * (1 - AF))) / (1 + exp(logit(AF) + abs(qnorm((1 - confidence.level) / 2)) * Std.Error / (AF * (1 - AF))))
  }
  CI.AF <- cbind(lower, upper)
  return(CI.AF)
}

#' @title Summary function for objects of class "\code{AF}".
#' @description Gives a summary of the AF estimate(s) including z-value, p-value and confidence interval(s).
#' @param object an object of class \code{AF} from \code{\link{AFglm}}, \code{\link{AFcoxph}} or \code{\link{AFclogit}} functions.
#' @param confidence.level user-specified confidence level for the confidence intervals. If not specified it defaults to 95 percent. Should be specified in decimals such as 0.95 for 95 percent.
#' @param CI.transform user-specified transformation of the Wald confidence interval(s). Options are \code{untransformed}, \code{log} and \code{logit}. If not specified untransformed will be calculated.
#' @param digits maximum number of digits.
#' @param ... further arguments to be passed to the summary function. See \code{\link[base]{summary}}.
#' @author Elisabeth Dahlqwist, Arvid \enc{Sjölander}{Sjolander}
#' @importFrom graphics legend lines plot.default
#' @export
summary.AF <- function(object, digits = max(3L, getOption("digits") - 3L),
                       confidence.level, CI.transform, ...){
  if(missing(confidence.level)) confidence.level <- 0.95
  if(missing(CI.transform)) CI.transform <- "untransformed"
  se <- sqrt(object$AF.var)
  zvalue <- object$AF.est / sqrt(object$AF.var)
  pvalue <- 2 * pnorm( - abs(zvalue))
  confidence.interval <- CI.AF(AF = object$AF.est, Std.Error = se,
                            confidence.level = confidence.level,
                            CI.transform = CI.transform)
  colnames(confidence.interval) <- c("Lower limit", "Upper limit")

  if(!object$n.cluster == 0) Std.Error <- "Robust SE"
      else Std.Error <- "Std.Error"
  AF <- cbind(object$AF.est, se, zvalue, pvalue)
  colnames(AF) <- c("AF estimate", Std.Error, "z value", "Pr(>|z|)")

  modelcall <- as.character(object$objectcall[1])
  if(modelcall == "glm") method = "Logistic regression"
  if(modelcall == "coxph") method = "Cox Proportional Hazards model"
  if(modelcall == "gee") method = "Conditional logistic regression"
  if(modelcall == "clogit") method = "Conditional logistic regression"
  if(modelcall == "parfrailty") method = "Weibull gamma-frailty model"

  if(modelcall == "coxph"  |  modelcall == "parfrailty"){
    ans <- list(AF = AF, times = object$times,
                CI.transform = CI.transform, confidence.level = confidence.level,
                confidence.interval = confidence.interval, n.obs = object$n,
                n.cases = object$n.cases, n.cluster = object$n.cluster,
                modelcall = modelcall, objectcall = object$objectcall, method = method, formula = object$formula,
                exposure = object$exposure, outcome = object$outcome, object = object,
                sandwich = object$sandwich, Std.Error = se, times = object$times)
  }
  else{
    ans <- list(AF = AF, times = object$times,
                CI.transform = CI.transform, confidence.level = confidence.level,
                confidence.interval = confidence.interval, n.obs = object$n,
                n.cases = object$n.cases, n.cluster = object$n.cluster,
                modelcall = modelcall, objectcall = object$objectcall, method = method, formula = object$formula,
                exposure = object$exposure, outcome = object$outcome, object = object,
                sandwich = object$sandwich, Std.Error = se)
  }
  class(ans) <- "summary.AF"
  return(ans)
}

#' @export
print.summary.AF <- function(x, digits = max(3L, getOption("digits") - 3L),
                             ...){
  cat("Call: ", "\n")
  print.default(x$objectcall)
  if(!x$n.cluster == 0) Std.Error <- "Robust SE"
  else Std.Error <- "Std.Error"
  if(x$CI.transform == "log") x$CI.transform <- "log transformed"
  if(x$CI.transform == "logit") x$CI.transform <- "logit transformed"
  level <- x$confidence.level * 100
  CI.text <- paste0(as.character(level),"%")
  cat("\nEstimated attributable fraction (AF) and", x$CI.transform, CI.text,  "Wald CI:", "\n")
  cat("\n")
  table.est <- cbind(x$AF, x$confidence.interval)
  colnames(table.est) <- c("AF", Std.Error, "z value", "Pr(>|z|)", "Lower limit", "Upper limit")
  r <- rep("", , nrow(x$AF))
  rownames(table.est) <- c(r)
  modelcall <- as.character(x$call[1])
  if(x$modelcall == "coxph"  |  x$modelcall == "parfrailty"){
    table.est <- cbind(x$times, table.est)
    colnames(table.est) <- c("Time", "AF", Std.Error, "z value", "Pr(>|z|)",
                             "Lower limit", "Upper limit")
    print.default(table.est)
  }
  else {
    print.default(table.est)
  }
  cat("\nExposure", ":", x$exposure, "\n")

  if(x$modelcall == "coxph" | x$modelcall == "parfrailty") outcome <- "Event   "
  else outcome <- "Outcome "
  #cat("\n")
  cat(outcome, ":", x$outcome, "\n")

  cat("\n")
  table.nr <- cbind(x$n.obs, x$n.cases)
  rownames(table.nr) <- c("")
  if(x$modelcall == "coxph" | x$modelcall == "parfrailty") number <- "Events"
  else number <- "Cases"
  colnames(table.nr) <- c("Observations", number)
  if (x$n.cluster == 0) print.default(table.nr)
  else{
    table.nr.cluster <- cbind(table.nr, x$n.cluster)
    colnames(table.nr.cluster) <- c("Observations", number, "Clusters")
    print.default(table.nr.cluster)
  }
  cat("\nMethod for confounder adjustment: ", x$method, "\n")
  formula <- as.character(x$formula)
  cat("\nFormula: ", formula[2], formula[1], formula[3], "\n")

  return(table.est)
}

#' @title Plot function for objects of class "\code{AF}" from the function \code{AFcoxph} or \code{AF.ch}.
#' @description Creates a simple scatterplot for the AF function with time sequence (specified by the user as \code{times} in the \code{\link{AFcoxph}} function) on the x-axis and the AF function estimate on the y-axis.
#' @param x an object of class \code{AF} from the \code{\link{AFcoxph}} function.
#' @param CI if TRUE confidence intervals are estimated and ploted in the graph.
#' @param confidence.level user-specified confidence level for the confidence intervals. If not specified it defaults to 95 percent. Should be specified in decimals such as 0.95 for 95 percent.
#' @param CI.transform user-specified transformation of the Wald confidence interval(s). Options are \code{untransformed}, \code{log} and \code{logit}. If not specified untransformed will be calculated.
#' @param xlab label on the x-axis. If not specified the label \emph{"Time"} will be displayed.
#' @param main main title of the plot. If not specified the lable \emph{"Estimate of the attributable fraction function"} will be displayed.
#' @param ylim limits on the y-axis of the plot. If not specified the minimum value of the lower bound of the confidence interval will be used as the minimal value and the maximum value of the upper bound of the confidence interval will be used as the maximum of y-axis of the plot.
#' @param ... further arguments to be passed to the plot function. See \code{\link[graphics]{plot}}.
#' @author Elisabeth Dahlqwist, Arvid \enc{Sjölander}{Sjolander}
#' @importFrom graphics legend lines plot.default
#' @export
plot.AF <- function(x, CI = TRUE, confidence.level,
                    CI.transform, xlab, main, ylim, ...){
  modelcall <- as.character(x$objectcall[1])
  if(modelcall != "coxph" & modelcall != "parfrailty")
    stop("Plot function is only available for the attributable fraction function. That is objects from the AFcoxph or AFparfrailty functions", call. = FALSE)
  if(missing(confidence.level)) confidence.level <- 0.95
  if(missing(CI.transform)) CI.transform <- "untransformed"
  if(missing(xlab)) xlab <- "Time"
  if(missing(main)) main <- ""
  if(CI == TRUE){
    confidence.interval <- CI.AF(AF = x$AF.est, Std.Error = sqrt(x$AF.var),
                                 confidence.level = confidence.level,
                                 CI.transform = CI.transform)
  if(missing(ylim)) ylim <- c(min(confidence.interval), max(confidence.interval))
  plot.default(x$times, x$AF.est, main = main,
               ylab = "Attributable fraction function" , xlab = xlab, ylim = ylim, pch = 19,
               lty = 1, type = "o", ...)
    lines( x$times, confidence.interval[, 2], lty = 2)
    lines( x$times, confidence.interval[, 1], lty = 2)
    level <- confidence.level * 100
    CI <- paste0(as.character(level),"% Conf. Interval")
    if(CI.transform == "log") transform <- "(log transformed)"
    if(CI.transform == "logit") transform <- "(logit transformed)"
    if(CI.transform == "untransformed") transform <- ""
    legend("topright", legend = c("AF estimate", CI, transform), pch = c(19, NA, NA), lty = c(1, 2, 0), bty = "n")
  }
}
