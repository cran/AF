############## AF function for a glm object #####################
#' @title Attributable fraction estimation based on a logistic regression model from a \code{glm} object (commonly used for cross-sectional or case-control sampling designs).
#' @description \code{AFglm} estimates the model-based adjusted attributable fraction for data from a logistic regression model in the form of a \code{\link[stats]{glm}} object. This model is commonly used for data from a cross-sectional or non-matched case-control sampling design.
#' @param object a fitted logistic regression model object of class "\code{\link[stats]{glm}}".
#' @param data an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment (\code{formula}), typically the environment from which the function is called.
#' @param exposure the name of the exposure variable as a string. The exposure must be binary (0/1) where unexposed is coded as 0.
#' @param clusterid the name of the cluster identifier variable as a string, if data are clustered. Cluster robust standard errors will be calculated.
#' @param case.control can be set to \code{TRUE} if the data is from a non-matched case control study. By default \code{case.control} is set to \code{FALSE} which is used for cross-sectional sampling designs.
#' @return \item{AF.est}{estimated attributable fraction.}
#' @return \item{AF.var}{estimated variance of \code{AF.est}. The variance is obtained by combining the delta method with the sandwich formula.}
#' @return \item{P.est}{estimated factual proportion of cases; \eqn{Pr(Y=1)}. Returned by default when \code{case.control = FALSE}.}
#' @return \item{P.var}{estimated variance of \code{P.est}. The variance is obtained by the sandwich formula. Returned by default when \code{case.control = FALSE}.}
#' @return \item{P0.est}{estimated counterfactual proportion of cases if exposure would be eliminated; \eqn{Pr(Y_0=1)}{Pr(Y0=1)}. Returned by default when \code{case.control = FALSE}.}
#' @return \item{P0.var}{estimated variance of \code{P0.est}. The variance is obtained by the sandwich formula. Returned by default when \code{case.control = FALSE}.}
#' @return \item{log.or}{a vector of the estimated log odds ratio for every individual. \code{log.or} contains the estimated coefficient for the exposure variable \code{X} for every level of the confounder \code{Z} as specified by the user in the formula. If the model to be estimated is
#'  \deqn{logit\{Pr(Y=1|X,Z)\} = \alpha+\beta{X}+\gamma{Z}}{logit {Pr(Y=1|X,Z)} = \alpha + \beta X + \gamma Z}
#'   then \code{log.or} is the estimate of \eqn{\beta}.
#'   If the model to be estimated is
#'   \deqn{logit\{Pr(Y=1|X,Z)\}=\alpha+\beta{X}+\gamma{Z}+\psi{XZ}}{logit{Pr(Y=1|X,Z)} = \alpha + \beta X +\gamma Z +\psi XZ}
#'   then \code{log.odds} is the estimate of
#'    \eqn{\beta + \psi{Z}}{\beta + \psi Z}. Only returned if argument \code{case.control} is set to \code{TRUE}.}
#' @details \code{AFglm} estimates the attributable fraction for a binary outcome \code{Y}
#' under the hypothetical scenario where a binary exposure \code{X} is eliminated from the population.
#' The estimate is adjusted for confounders \code{Z} by logistic regression using the (\code{\link[stats]{glm}}) function.
#' The estimation strategy is different for cross-sectional and case-control sampling designs even if the underlying logististic regression model is the same.
#' For cross-sectional sampling designs the AF can be defined as
#' \deqn{AF=1-\frac{Pr(Y_0=1)}{Pr(Y=1)}}{AF = 1 - Pr(Y0 = 1) / Pr(Y = 1)}
#' where \eqn{Pr(Y_0=1)}{Pr(Y0 = 1)} denotes the counterfactual probability of the outcome if
#' the exposure would have been eliminated from the population and \eqn{Pr(Y = 1)} denotes the factual probability of the outcome.
#' If \code{Z} is sufficient for confounding control, then \eqn{Pr(Y_0=1)}{Pr(Y0 = 1)} can be expressed as
#'  \eqn{E_Z\{Pr(Y=1\mid{X=0,Z})\}.}{E_z{Pr(Y = 1  |X = 0,Z)}.}
#' The function uses logistic regression to estimate \eqn{Pr(Y=1\mid{X=0,Z})}{Pr(Y=1|X=0,Z)}, and the marginal sample distribution of \code{Z}
#' to approximate the outer expectation (\enc{Sjölander}{Sjolander} and Vansteelandt, 2012).
#' For case-control sampling designs the outcome prevalence is fixed by sampling design and absolute probabilities (\code{P.est} and \code{P0.est}) can not be estimated.
#' Instead adjusted log odds ratios (\code{log.or}) are estimated for each individual.
#' This is done by setting \code{case.control} to \code{TRUE}. It is then assumed that the outcome is rare so that the risk ratio can be approximated by the odds ratio.
#' For case-control sampling designs the AF be defined as (Bruzzi et. al)
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
#' If \code{clusterid} is supplied, then a clustered sandwich formula is used in all variance calculations.
#' @author Elisabeth Dahlqwist, Arvid \enc{Sjölander}{Sjolander}
#' @seealso \code{\link[stats]{glm}} used for fitting the logistic regression model. For conditional logistic regression (commonly for data from a matched case-control sampling design) see \code{\link[AF]{AFclogit}}.
#' @references Bruzzi, P., Green, S. B., Byar, D., Brinton, L. A., and Schairer, C. (1985). Estimating the population attributable risk for multiple risk factors using case-control data. \emph{American Journal of Epidemiology} \bold{122}, 904-914.
#' @references Greenland, S. and Drescher, K. (1993). Maximum Likelihood Estimation of the Attributable Fraction from logistic Models. \emph{Biometrics} \bold{49}, 865-872.
#' @references \enc{Sjölander}{Sjolander}, A. and Vansteelandt, S. (2011). Doubly robust estimation of attributable fractions. \emph{Biostatistics} \bold{12}, 112-121.
#' @examples
#' # Simulate a cross-sectional sample
#'
#' expit <- function(x) 1 / (1 + exp( - x))
#' n <- 1000
#' Z <- rnorm(n = n)
#' X <- rbinom(n = n, size = 1, prob = expit(Z))
#' Y <- rbinom(n = n, size = 1, prob = expit(Z + X))
#'
#' # Example 1: non clustered data from a cross-sectional sampling design
#' data <- data.frame(Y, X, Z)
#'
#' # Fit a glm object
#' fit <- glm(formula = Y ~ X + Z + X * Z, family = binomial, data = data)
#'
#' # Estimate the attributable fraction from the fitted logistic regression
#' AFglm_est <- AFglm(object = fit, data = data, exposure = "X")
#' summary(AFglm_est)
#'
#' # Example 2: clustered data from a cross-sectional sampling design
#' # Duplicate observations in order to create clustered data
#' id <- rep(1:n, 2)
#' data <- data.frame(id = id, Y = c(Y, Y), X = c(X, X), Z = c(Z, Z))
#'
#' # Fit a glm object
#' fit <- glm(formula = Y ~ X + Z + X * Z, family = binomial, data = data)
#'
#' # Estimate the attributable fraction from the fitted logistic regression
#' AFglm_clust <- AFglm(object = fit, data = data,
#'                          exposure = "X", clusterid = "id")
#' summary(AFglm_clust)
#'
#'
#' # Example 3: non matched case-control
#' # Simulate a sample from a non matched case-control sampling design
#' # Make the outcome a rare event by setting the intercept to -6
#'
#' expit <- function(x) 1 / (1 + exp( - x))
#' NN <- 1000000
#' n <- 500
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
#' # Fit a glm object
#' fit <- glm(formula = Y ~ X + Z + X * Z, family = binomial, data = data)
#'
#' # Estimate the attributable fraction from the fitted logistic regression
#' AFglm_est_cc <- AFglm(object = fit, data = data, exposure = "X", case.control = TRUE)
#' summary(AFglm_est_cc)
#' @importFrom stats aggregate as.formula ave binomial coef delete.response family model.matrix pnorm predict qnorm residuals stepfun terms var vcov
#' @import data.table
#' @export
AFglm <- function(object, data, exposure, clusterid, case.control = FALSE){
  call <- match.call()
  # Warning if the object is not a glm object
  if(!(as.character(object$call[1]) == "glm"))
    stop("The object is not a glm object", call. = FALSE)
  # Warning if the object is not a logistic regression
  if(!(object$family[1] == "binomial" & object$family[2] == "logit"))
    stop("The object is not a logistic regression", call. = FALSE)
  #### Preparation of dataset ####
  formula <- object$formula
  #data <- object$data
  npar <- length(object$coef)

  ## Delete rows with missing on variables in the model ##
  rownames(data) <- 1:nrow(data)
  m <- model.matrix(object = formula, data = data)
  complete <- as.numeric(rownames(m))
  data <- data[complete, ]
  outcome <- as.character(terms(formula)[[2]])
  n <- nrow(data)
  n.cases <- sum(data[, outcome])
  clusters <- data[, clusterid]

  if(missing(clusterid)) n.cluster <- 0
  else {
    n.cluster <- length(unique(data[, clusterid]))
  }

  ## Checks ##
  if(!class(exposure) == "character") 
    stop("Exposure must be a string.", call. = FALSE)
  if(!is.binary(data[, exposure]))
    stop("Only binary exposure (0/1) is accepted.", call. = FALSE)
  if(max(all.vars(formula[[3]]) == exposure) == 0)
    stop("The exposure variable is not included in the formula.", call. = FALSE)

  # Create dataset data0 for counterfactual X = 0
  data0 <- data
  data0[, exposure] <- 0

  ## Design matrices ##
  design <- model.matrix(object = delete.response(terms(object)), data = data)
  design0 <- model.matrix(object = delete.response(terms(object)), data = data0)

  #### Meat: score equations ####
  ## If sampling design is case-control ##
  if (case.control == TRUE){
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
    pred.diff <- data[, outcome] - predict(object, newdata = data, type = "response")
    score.beta <- design * pred.diff
    score.equations <- cbind(score.AF, score.beta)
    if (!missing(clusterid)){
      score.equations <- score.equations
      score.equations <- aggr(score.equations, clusters = clusters)
    }
    meat <- var(score.equations, na.rm=TRUE)
    #### Bread: hessian of score equations ####
    ## Hessian of score equation 1 ##
    #### Estimating variance using Sandwich estimator ####
    hessian.AF1 <- - data[, outcome]
    hessian.AF2 <- (design0 - design) * as.vector(data[, outcome] * exp( - log.or))
    hessian.AF <- cbind(mean(hessian.AF1), t(colMeans(hessian.AF2, na.rm = TRUE)))
    hessian.beta <- cbind(matrix(rep(0, npar), nrow = npar, ncol = 1), - solve(vcov(object = object)) / n)
    ### Bread ###
    bread <- rbind(hessian.AF, hessian.beta)
    #### Sandwich ####
    if (!missing(clusterid))
      sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) * n.cluster / n^2 ) [1:2, 1:2]
    else
      sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) / n) [1:2, 1:2]
    AF.var <- sandwich[1, 1]
    #### Output ####
    out <- c(list(AF.est = AF.est, AF.var = AF.var, log.or = log.or,
                  objectcall = object$call, call = call, exposure = exposure, outcome = outcome, object = object,
                  sandwich = sandwich, formula = formula,
                  n = n, n.cases = n.cases, n.cluster = n.cluster))
  }
  ## If sampling design is cross-sectional ##
  else {
    ## Score equation 1 ##
    score.P <- data[, outcome]
    pred.Y  <- predict(object, newdata = data, type = "response")
    ## Score equation 2 ##
    score.P0 <- predict(object, newdata = data0, type = "response")
    ## Score equation 3 ##
    score.beta <- design * (score.P - pred.Y)
    ### Meat ###
    score.equations <- cbind(score.P, score.P0, score.beta)
    if (!missing(clusterid)){
      score.equations <- score.equations
      score.equations <- aggr(score.equations, clusters = clusters)
    }
    meat <- var(score.equations, na.rm = TRUE)
    #### Bread: hessian of score equations ####
    ## Hessian of score equation 1 ##
    hessian.P <- matrix(c(- 1, 0, rep(0,npar)), nrow = 1, ncol = 2 + npar)
    ## Hessian of score equation 2 ##
    g <- family(object)$mu.eta
    dmu.deta <- g(predict(object = object, newdata = data0))
    deta.dbeta <- design0
    dmu.dbeta <- dmu.deta * deta.dbeta
    hessian.P0 <- matrix(c(0, - 1, colMeans(dmu.dbeta)), nrow = 1, ncol = 2 + npar)
    ## Hessian of score equation 3 ##
    hessian.beta <- cbind(matrix(rep(0, npar * 2), nrow = npar, ncol = 2)
                          , - solve(vcov(object = object)) / n)
    ### Bread ###
    bread <- rbind(hessian.P, hessian.P0, hessian.beta)
    #### Sandwich ####
    if (!missing(clusterid))
      sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) * n.cluster / n^2 ) [1:2, 1:2]
    else
      sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) / n) [1:2, 1:2]
    #### Point estimate of AF ####
    P.est  <- mean(score.P, na.rm = TRUE)
    P0.est <- mean(score.P0, na.rm = TRUE)
    AF.est <- 1 - P0.est / P.est
    ## Delta method for variance estimate ##
    gradient <- as.matrix(c(P0.est / P.est ^ 2, - 1 / P.est), nrow = 2, ncol = 1)
    AF.var <- t(gradient) %*% sandwich %*% gradient
    P.var <- sandwich[1, 1]
    P0.var <- sandwich[2, 2]

    objectcall <- object$call
    #### Output ####
    out <- c(list(AF.est = AF.est, AF.var = AF.var, P.est = P.est, P0.est = P0.est, P.var = P.var,
                  P0.var = P0.var, objectcall = objectcall, call = call, exposure = exposure, outcome = outcome,
                  object = object, sandwich = sandwich, gradient = gradient, formula = formula,
                  n = n, n.cases = n.cases, n.cluster = n.cluster))
  }
  class(out) <- "AF"
  return(out)
}
