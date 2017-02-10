############## AF function for matched and unmatched case-control #####################
#' @title Attributable fraction for mached and non-matched case-control sampling designs. NOTE! Deprecated function. Use \code{\link[AF]{AFglm}} (for unmatched case-control studies) or \code{\link[AF]{AFclogit}} (for matched case-control studies).
#' @description \code{AF.cc} estimates the model-based adjusted attributable fraction for data from matched and non-matched case-control sampling designs.
#' @param formula an object of class "\code{formula}" (or one that can be coerced to that class): a symbolic description of the model used for confounder adjustment. The exposure and confounders should be specified as independent (right-hand side) variables. The outcome should be specified as dependent (left-hand side) variable. The formula is used to object a logistic regression by \code{\link[stats]{glm}} for non-matched case-control and conditional logistic regression by \code{\link[drgee]{gee}} (in package \code{\link[drgee]{drgee}}) for matched case-control.
#' @param data an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in \code{data}, the variables are taken from environment (\code{formula}), typically the environment from which the function is called.
#' @param exposure the name of the exposure variable as a string. The exposure must be binary (0/1) where unexposed is coded as 0.
#' @param matched a logical that specifies if the sampling design is matched (TRUE) or non-matched (FALSE) case-control. Default setting is non-matched (\code{matched = FALSE}).
#' @param clusterid the name of the cluster identifier variable as a string, if data are clustered (e.g. matched).
#' @return \item{AF.est}{estimated attributable fraction.}
#' @return \item{AF.var}{estimated variance of \code{AF.est}. The variance is obtained by combining the delta methods with the sandwich formula.}
#' @return \item{log.or}{a vector of the estimated log odds ratio for every individual. \code{log.or} contains the estimated coefficient for the exposure variable \code{X} for every level of the confounder \code{Z} as specified by the user in the formula. If the model to be estimated is
#'  \deqn{logit\{Pr(Y=1|X,Z)\} = \alpha+\beta{X}+\gamma{Z}}{logit {Pr(Y=1|X,Z)} = \alpha + \beta X + \gamma Z}
#'   then \code{log.or} is the estimate of \eqn{\beta}.
#'   If the model to be estimated is
#'   \deqn{logit\{Pr(Y=1|X,Z)\}=\alpha+\beta{X}+\gamma{Z}+\psi{XZ}}{logit{Pr(Y=1|X,Z)} = \alpha + \beta X +\gamma Z +\psi XZ}
#'   then \code{log.odds} is the estimate of
#'    \eqn{\beta + \psi{Z}}{\beta + \psi Z}.}
#' @return \item{object}{the fitted model. Fitted using logistic regression, \code{\link{glm}}, for non-matched case-control and conditional logistic regression, \code{\link[drgee]{gee}}, for matched case-control.}
#' @details \code{Af.cc} estimates the attributable fraction for a binary outcome \code{Y}
#' under the hypothetical scenario where a binary exposure \code{X} is eliminated from the population.
#' The estimate is adjusted for confounders \code{Z} by logistic regression for unmatched case-control (\code{\link[stats]{glm}}) and conditional logistic regression for matched case-control (\code{\link[drgee]{gee}}).
#' The estimation assumes that the outcome is rare so that the risk ratio can be approximated by the odds ratio, for details see Bruzzi et. al.
#' Let the AF be defined as
#' \deqn{AF = 1 - \frac{Pr(Y_0=1)}{Pr(Y = 1)}}{AF = 1 - Pr(Y0 = 1) / Pr(Y = 1)}
#' where \eqn{Pr(Y_0=1)}{Pr(Y0 = 1)} denotes the counterfactual probability of the outcome if
#' the exposure would have been eliminated from the population. If \code{Z} is sufficient for confounding control then the probability \eqn{Pr(Y_0=1)}{Pr(Y0 = 1)} can be expressed as
#' \deqn{Pr(Y_0=1)=E_Z\{Pr(Y=1\mid{X}=0,Z)\}.}{Pr(Y0=1) = E_z{Pr(Y = 1 | X = 0, Z)}.}
#' Using Bayes' theorem this implies that the AF can be expressed as
#' \deqn{AF = 1-\frac{E_Z\{Pr(Y=1\mid X=0,Z)\}}{Pr(Y=1)}=1-E_Z\{RR^{-X}(Z)\mid{Y = 1}\}}{
#' AF = 1 - E_z{Pr( Y = 1 | X = 0, Z)} / Pr(Y = 1) = 1 - E_z{RR^{-X} (Z) | Y = 1}}
#' where \eqn{RR(Z)} is the risk ratio \deqn{\frac{Pr(Y=1\mid{X=1,Z})}{Pr(Y=1\mid{X=0,Z})}.}{Pr(Y = 1 | X = 1,Z)/Pr(Y=1 | X = 0, Z).}
#' Moreover, the risk ratio can be approximated by the odds ratio if the outcome is rare. Thus,
#' \deqn{ AF \approx 1 - E_Z\{OR^{-X}(Z)\mid{Y = 1}\}.}{AF is approximately 1 - E_z{OR^{-X}(Z) | Y = 1}.}
#' The odds ratio is estimated by logistic regression or conditional logistic regression.
#' If \code{clusterid} is supplied, then a clustered sandwich formula is used in all variance calculations.
#' @author Elisabeth Dahlqwist, Arvid \enc{Sjölander}{Sjolander}
#' @seealso The new and more general version of the function: \code{\link[AF]{AFglm}} for non-matched and \code{\link[AF]{AFclogit}} for matched case-control sampling designs. \code{\link[stats]{glm}} and \code{\link[drgee]{gee}} used for fitting the logistic regression model (for non-matched case-control) and the conditional logistic regression model (for matched case-control).
#' @references Bruzzi, P., Green, S. B., Byar, D., Brinton, L. A., and Schairer, C. (1985). Estimating the population attributable risk for multiple risk factors using case-control data. \emph{American Journal of Epidemiology} \bold{122}, 904-914.
#' @examples
#' expit <- function(x) 1 / (1 + exp( - x))
#' NN <- 1000000
#' n <- 500
#'
#' # Example 1: non matched case-control
#' # Simulate a sample from a non matched case-control sampling design
#' # Make the outcome a rare event by setting the intercept to -6
#' intercept <- -6
#' Z <- rnorm(n = NN)
#' X <- rbinom(n = NN, size = 1, prob = expit(Z))
#' Y <- rbinom(n = NN, size = 1, prob = expit(intercept + X + Z))
#' population <- data.frame(Z, X, Y)
#' Case <- which(population$Y == 1)
#' Control <- which(population$Y == 0)
#' # Sample cases and controls from the population
#' case <- sample(Case, n)
#' control <- sample(Control, n)
#' data <- population[c(case, control), ]
#'
#' # Estimation of the attributable fraction
#' AF.cc_est <- AF.cc(formula = Y ~ X + Z + X * Z, data = data, exposure = "X")
#' summary(AF.cc_est)
#'
#' # Example 2: matched case-control
#' # Duplicate observations in order to create a matched data sample
#' # Create an unobserved confounder U common for each pair of individuals
#' U  <- rnorm(n = NN)
#' Z1 <- rnorm(n = NN)
#' Z2 <- rnorm(n = NN)
#' X1 <- rbinom(n = NN, size = 1, prob = expit(U + Z1))
#' X2 <- rbinom(n = NN, size = 1, prob = expit(U + Z2))
#' Y1 <- rbinom(n = NN, size = 1, prob = expit(intercept + U + Z1 + X1))
#' Y2 <- rbinom(n = NN, size = 1, prob = expit(intercept + U + Z2 + X2))
#' # Select discordant pairs
#' discordant <- which(Y1!=Y2)
#' id <- rep(1:n, 2)
#' # Sample from discordant pairs
#' incl <- sample(x = discordant, size = n, replace = TRUE)
#' data <- data.frame(id = id, Y = c(Y1[incl], Y2[incl]), X = c(X1[incl], X2[incl]),
#'                    Z = c(Z1[incl], Z2[incl]))
#'
#' # Estimation of the attributable fraction
#' AF.cc_match <- AF.cc(formula = Y ~ X + Z + X * Z, data = data,
#'                          exposure = "X", clusterid = "id", matched = TRUE)
#' summary(AF.cc_match)
#' @import drgee
#' @export
AF.cc<-function(formula, data, exposure, clusterid,
                matched = FALSE){
  warning("NOTE! Deprecated function. Use AFglm (for unmatched case-control studies) or AFclogit (for matched case-control studies).", call = FALSE)
  call <- match.call()
  mm <- match(c("formula", "data", "exposure", "clusterid", "matched"), names(call), 0L)
  #### Preparation of dataset ####
  ## Delete rows with missing on variables in the model ##
  rownames(data) <- 1:nrow(data)
  m <- model.matrix(object = formula, data = data)
  complete <- as.numeric(rownames(m))
  data <- data[complete, ]
  outcome <- as.character(terms(formula)[[2]])
  if(matched == TRUE){
    ni.vals <- ave(as.vector(data[, outcome]), data[, clusterid], FUN = function(y) {
      length(unique(y[which(!is.na(y))]))
    })
    compl.rows <- (ni.vals > 1)
    data <- data[compl.rows, ]
  }
  ## Checks ##
  if(is.binary(data[, outcome]) == FALSE)
    stop("Only binary outcome (0/1) is accepted.", call. = FALSE)
  if(is.binary(data[, exposure]) == FALSE)
    stop("Only binary exposure (0/1) is accepted.", call. = FALSE)
  if(max(all.vars(formula[[3]]) == exposure) == 0)
    stop("The exposure variable is not included in the formula.", call. = FALSE)
  if(missing(clusterid)) n.cluster <- 0
  else n.cluster <- length(unique(data[, clusterid]))
  #### Methods for non-matched or matched sampling designs ####
  n <- nrow(data)
  n.cases <- sum(data[, outcome])
  if (!missing(clusterid))
   data <- data[order(data[, clusterid]), ]
  data0 <- data
  data0[, exposure] <- 0
  #### Estimate model ####
  if(matched == FALSE)
    object <- glm(formula = formula, family = binomial, data = data)
  if(matched == TRUE)
    object <- gee(formula, link = "logit", data, cond = TRUE, clusterid = clusterid)
  npar <- length(object$coef)
  ## Design matrices ##
  if(matched == FALSE){
    design <- as.matrix(model.matrix(object = delete.response(terms(object)), data = data))
    design0 <- as.matrix(model.matrix(object = delete.response(terms(object)), data = data0))
  }
  if(matched == TRUE){
     design <- as.matrix(model.matrix(object = formula, data = data)[, - 1])
    design0 <- as.matrix(model.matrix(object = formula, data = data0)[, - 1])
  }
  ## Create linear predictors to estimate the log odds ratio ##
  diff.design <- design0 - design
  linearpredictor <- design  %*% coef(object)
  linearpredictor0 <- design0 %*% coef(object)
  #log odds ratio#
  log.or <- linearpredictor - linearpredictor0
  ## Estimate approximate AF ##
  AF.est   <- 1 - sum(data[, outcome] * exp( - log.or)) / sum(data[, outcome])
  #### Meat: score equations ####
  ## Score equation 1 ## individual estimating equations of the estimate of AF
  score.AF <- data[, outcome] * (exp( - log.or) - AF.est)
  ## Score equation 2 ## individual estimating equations from conditional logistic reg.
  if(matched == FALSE)
    pred.diff <- data[, outcome] - predict(object, newdata = data, type = "response")
  if(matched == TRUE)
    pred.diff <- object$res
  score.beta <- design * pred.diff
  score.equations <- cbind(score.AF, score.beta)
  if (!missing(clusterid))
    score.equations <- aggregate(score.equations, list(data[, clusterid]), sum)[, - 1]
  meat <- var(score.equations, na.rm=TRUE)
  #### Bread: hessian of score equations ####
  ## Hessian of score equation 1 ##
  #### Estimating variance using Sandwich estimator ####
  hessian.AF1 <- - data[, outcome]

  hessian.AF2 <- (design0 - design) * as.vector(data[, outcome] * exp( - log.or))
  if (!missing(clusterid)){
    if(length(all.vars(formula[[3]]))>1){
      hessian.AF <- cbind(mean(aggregate(hessian.AF1, list(data[, clusterid]), sum)[, - 1], na.rm=TRUE)
                          , t(colMeans(aggregate(hessian.AF2
                                                 , list(data[, clusterid]), sum)[, - 1], na.rm = TRUE)))
    }
    if(length(all.vars(formula[[3]]))==1){
      hessian.AF <- cbind(mean(aggregate(hessian.AF1, list(data[, clusterid]), sum)[, - 1], na.rm=TRUE)
                          , t(mean(aggregate(hessian.AF2
                                             , list(data[, clusterid]), sum)[, - 1], na.rm = TRUE)))
    }
  }
  else
    hessian.AF <- cbind(mean(hessian.AF1), t(colMeans(hessian.AF2, na.rm = TRUE)))
  hessian.beta <- cbind(matrix(rep(0, npar), nrow = npar, ncol = 1), - solve(vcov(object = object)) / n)
  ### Bread ###
  bread <- rbind(hessian.AF, hessian.beta)
  #### Sandwich ####
  if (!missing(clusterid))
    sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) * n.cluster/ n ^ 2 ) [1:2, 1:2]
  else
    sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) / n) [1:2, 1:2]
  AF.var <- sandwich[1, 1]

  clusterid <- data[, clusterid]
  #### Output ####
  out <- c(list(hessian.beta = hessian.beta, hessian.AF= hessian.AF,clusterid = clusterid, score.equations= score.equations, hessian.beta = hessian.beta, bread = bread, meat = meat, AF.est = AF.est, AF.var = AF.var, log.or = log.or,
                objectcall = object$call, call = call, exposure = exposure, outcome = outcome, object = object,
                sandwich = sandwich, formula = formula,
                n = n, n.cases = n.cases, n.cluster = n.cluster))
  class(out) <- "AF"
  return(out)
}
