############## AF function for cross-sectional sampling design #####################
#' @title Attributable fraction for cross-sectional sampling designs. NOTE! Deprecated function. Use \code{\link[AF]{AFglm}}.
#' @description \code{AF.cs} estimates the model-based adjusted attributable fraction for data from cross-sectional sampling designs.
#' @param formula an object of class "\code{\link{formula}}" (or one that can be coerced to that class): a symbolic description of the model used for adjusting for confounders. The exposure and confounders should be specified as independent (right-hand side) variables. The outcome should be specified as dependent (left-hand side) variable. The formula is used to object a logistic regression by \code{\link[stats]{glm}}.
#' @param data an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment (\code{formula}), typically the environment from which the function is called.
#' @param exposure the name of the exposure variable as a string. The exposure must be binary (0/1) where unexposed is coded as 0.
#' @param clusterid the name of the cluster identifier variable as a string, if data are clustered.
#' @return \item{AF.est}{estimated attributable fraction.}
#' @return \item{AF.var}{estimated variance of \code{AF.est}. The variance is obtained by combining the delta method with the sandwich formula.}
#' @return \item{P.est}{estimated factual proportion of cases; \eqn{Pr(Y=1)}.}
#' @return \item{P.var}{estimated variance of \code{P.est}. The variance is obtained by the sandwich formula.}
#' @return \item{P0.est}{estimated counterfactual proportion of cases if exposure would be eliminated; \eqn{Pr(Y_0=1)}{Pr(Y0=1)}.}
#' @return \item{P0.var}{estimated variance of \code{P0.est}. The variance is obtained by the sandwich formula.}
#' @return \item{object}{the fitted model. Fitted using logistic regression, \code{\link{glm}}.}
#' @details \code{Af.cs} estimates the attributable fraction for a binary outcome \code{Y}
#' under the hypothetical scenario where a binary exposure \code{X} is eliminated from the population.
#' The estimate is adjusted for confounders \code{Z} by logistic regression (\code{\link{glm}}).
#' Let the AF be defined as
#' \deqn{AF=1-\frac{Pr(Y_0=1)}{Pr(Y=1)}}{AF = 1 - Pr(Y0 = 1) / Pr(Y = 1)}
#' where \eqn{Pr(Y_0=1)}{Pr(Y0 = 1)} denotes the counterfactual probability of the outcome if
#' the exposure would have been eliminated from the population and \eqn{Pr(Y = 1)} denotes the factual probability of the outcome.
#' If \code{Z} is sufficient for confounding control, then \eqn{Pr(Y_0=1)}{Pr(Y0 = 1)} can be expressed as
#'  \eqn{E_Z\{Pr(Y=1\mid{X=0,Z})\}.}{E_z{Pr(Y = 1  |X = 0,Z)}.}
#' The function uses logistic regression to estimate \eqn{Pr(Y=1\mid{X=0,Z})}{Pr(Y=1|X=0,Z)}, and the marginal sample distribution of \code{Z}
#' to approximate the outer expectation (\enc{Sjölander}{Sjolander} and Vansteelandt, 2012).
#' If \code{clusterid} is supplied, then a clustered sandwich formula is used in all variance calculations.
#' @author Elisabeth Dahlqwist, Arvid \enc{Sjölander}{Sjolander}
#' @seealso The new and more general version of the function: \code{\link[AF]{AFglm}}.
#' @references Greenland, S. and Drescher, K. (1993). Maximum Likelihood Estimation of the Attributable Fraction from logistic Models. \emph{Biometrics} \bold{49}, 865-872.
#' @references \enc{Sjölander}{Sjolander}, A. and Vansteelandt, S. (2011). Doubly robust estimation of attributable fractions. \emph{Biostatistics} \bold{12}, 112-121.
#' @examples
#' # Simulate a cross-sectional sample
#' expit <- function(x) 1 / (1 + exp( - x))
#' n <- 1000
#' Z <- rnorm(n = n)
#' X <- rbinom(n = n, size = 1, prob = expit(Z))
#' Y <- rbinom(n = n, size = 1, prob = expit(Z + X))
#'
#' # Example 1: non clustered data from a cross-sectional sampling design
#' data <- data.frame(Y, X, Z)
#'
#' # Estimation of the attributable fraction
#' AF.cs_est <- AF.cs(formula = Y ~ X + Z + X * Z, data = data, exposure = "X")
#' summary(AF.cs_est)
#'
#' # Example 2: clustered data from a cross-sectional sampling design
#' # Duplicate observations in order to create clustered data
#' id <- rep(1:n, 2)
#' data <- data.frame(id = id, Y = c(Y, Y), X = c(X, X), Z = c(Z, Z))
#'
#' # Estimation of the attributable fraction
#' AF.cs_clust <- AF.cs(formula = Y ~ X + Z + X * Z, data = data,
#'                          exposure = "X", clusterid = "id")
#' summary(AF.cs_clust)
#' @importFrom stats aggregate ave binomial coef delete.response family glm model.matrix pnorm predict qnorm residuals stepfun terms var vcov
#' @export
AF.cs<- function(formula, data, exposure, clusterid){
  warning("NOTE! Deprecated function. Use AFglm.")
  call <- match.call()
  mm <- match(c("formula", "data", "exposure", "clusterid"), names(call), 0L)
  #### Preparation of dataset ####
  ## Delete rows with missing on variables in the model ##
  rownames(data) <- 1:nrow(data)
  m <- model.matrix(object = formula, data = data)
  complete <- as.numeric(rownames(m))
  data <- data[complete, ]
  outcome <- as.character(terms(formula)[[2]])
  n <- nrow(data)
  n.cases <- sum(data[, outcome])
  if(missing(clusterid)) n.cluster <- 0
  else {
    n.cluster <- length(unique(data[, clusterid]))
    data <- data[order(data[, clusterid]), ]
  }
  ## Checks ##
  if(!is.binary(data[, outcome]))
    stop("Only binary outcome (0/1) is accepted.", call. = FALSE)
  if(!is.binary(data[, exposure]))
    stop("Only binary exposure (0/1) is accepted.", call. = FALSE)
  if(max(all.vars(formula[[3]]) == exposure) == 0)
    stop("The exposure variable is not included in the formula.", call. = FALSE)
  ## Counterfactual dataset ##
  data0 <- data
  data0[, exposure] <- 0
  #### Estimate model ####
  object <- glm(formula = formula, family = binomial, data = data)
  npar <- length(object$coef)
  ## Design matrices ##
  design <- model.matrix(object = delete.response(terms(object)), data = data)
  design0 <- model.matrix(object = delete.response(terms(object)), data = data0)
  #### Meat: score equations ####
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
    score.equations <- aggregate(score.equations, list(data[, clusterid]), sum)[, - 1]
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
  #### Output ####
  out <- c(list(AF.est = AF.est, AF.var = AF.var, P.est = P.est, P0.est = P0.est, P.var = P.var,
                P0.var = P0.var, objectcall = object$call, call = call, exposure = exposure, outcome = outcome,
                object = object, sandwich = sandwich, gradient = gradient, formula = formula,
                n = n, n.cases = n.cases, n.cluster = n.cluster))
  class(out) <- "AF"
  return(out)
}
