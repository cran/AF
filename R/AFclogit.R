############## AF function for a clogit object #####################
#' @title Attributable fraction estimation based on a conditional logistic regression model as a \code{clogit} object (commonly used for matched case-control sampling designs).
#' @description \code{AFclogit} estimates the model-based adjusted attributable fraction from a conditional logistic regression model in form of a \code{\link[survival]{clogit}} object. This model is model is commonly used for data from matched case-control sampling designs.
#' @param object a fitted conditional logistic regression model object of class "\code{\link[survival]{clogit}}".
#' @param data an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in \code{data}, the variables are taken from environment (\code{formula}), typically the environment from which the function is called.
#' @param exposure the name of the exposure variable as a string. The exposure must be binary (0/1) where unexposed is coded as 0.
#' @param clusterid the name of the cluster identifier variable as a string. Because conditional logistic regression is only used for clustered data, this argument must be supplied.
#' @return \item{AF.est}{estimated attributable fraction.}
#' @return \item{AF.var}{estimated variance of \code{AF.est}. The variance is obtained by combining the delta methods with the sandwich formula.}
#' @return \item{log.or}{a vector of the estimated log odds ratio for every individual. \code{log.or} contains the estimated coefficient for the exposure variable \code{X} for every level of the confounder \code{Z} as specified by the user in the formula. If the model to be estimated is
#'  \deqn{logit\{Pr(Y=1|X,Z)\} = \alpha+\beta{X}+\gamma{Z}}{logit {Pr(Y=1|X,Z)} = \alpha + \beta X + \gamma Z}
#'   then \code{log.or} is the estimate of \eqn{\beta}.
#'   If the model to be estimated is
#'   \deqn{logit\{Pr(Y=1|X,Z)\}=\alpha+\beta{X}+\gamma{Z}+\psi{XZ}}{logit{Pr(Y=1|X,Z)} = \alpha + \beta X +\gamma Z +\psi XZ}
#'   then \code{log.odds} is the estimate of
#'    \eqn{\beta + \psi{Z}}{\beta + \psi Z}.}
#' @details \code{AFclogit} estimates the attributable fraction for a binary outcome \code{Y}
#' under the hypothetical scenario where a binary exposure \code{X} is eliminated from the population.
#' The estimate is adjusted for confounders \code{Z} by conditional logistic regression.
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
#' The odds ratio is estimated by conditional logistic regression.
#' The function \code{\link[drgee]{gee}} in the \code{drgee} package is used to get the score contributions for each cluster and the hessian.
#' A clustered sandwich formula is used in the variance calculation.
#' @author Elisabeth Dahlqwist, Arvid \enc{Sjölander}{Sjolander}
#' @seealso \code{\link[survival]{clogit}} used for fitting the conditional logistic regression model for matched case-control designs. For non-matched case-control designs see \code{\link[AF]{AFglm}}.
#' @references Bruzzi, P., Green, S. B., Byar, D., Brinton, L. A., and Schairer, C. (1985). Estimating the population attributable risk for multiple risk factors using case-control data. \emph{American Journal of Epidemiology} \bold{122}, 904-914.
#' @examples
#' expit <- function(x) 1 / (1 + exp( - x))
#' NN <- 1000000
#' n <- 500
#'
#' # Example 1: matched case-control
#' # Duplicate observations in order to create a matched data sample
#' # Create an unobserved confounder U common for each pair of individuals
#' intercept <- -6
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
#' # Fit a clogit object
#' fit <- clogit(formula = Y ~ X + Z + X * Z + strata(id), data = data)
#'
#' # Estimate the attributable fraction from the fitted conditional logistic regression
#' AFclogit_est <- AFclogit(fit, data, exposure = "X", clusterid="id")
#' summary(AFclogit_est)
#' @import survival drgee data.table
#' @export
AFclogit<-function(object, data, exposure, clusterid){
  call <- match.call()
  # Warning if the object is not a clogit object
  objectcall <- object$userCall[1]
  if(!(class(object)[1])=="clogit")
    stop("The object is not a clogit object", call. = FALSE)
  if(missing(clusterid))
    stop("Argument 'clusterid' must be provided by the user", call. = FALSE)
  #### Preparation of variables ####
  formula <- object$formula
  npar <- length(object$coef)

  ## Delete rows with missing on variables in the model ##
  rownames(data) <- 1:nrow(data)
  m <- model.matrix(object = formula, data = data)
  complete <- as.numeric(rownames(m))
  data <- data[complete, ]
  outcome <- as.character(terms(formula)[[2]])[3]
  variables <- attr(object$coefficients, "names")
  ## Create a formula which can be used to create a design matrix
  formula.model <- as.formula(paste(outcome, "~", paste(variables, collapse=" + ")))
  ni.vals <- ave(as.vector(data[, outcome]), data[, clusterid], FUN = function(y) {
    length(unique(y[which(!is.na(y))]))
  })
  compl.rows <- (ni.vals > 1)
  data <- data[compl.rows, ]

  ## Checks ##
  if(is.binary(data[, outcome]) == FALSE)
    stop("Only binary outcome (0/1) is accepted.", call. = FALSE)
  if(is.binary(data[, exposure]) == FALSE)
    stop("Only binary exposure (0/1) is accepted.", call. = FALSE)
  if(max(all.vars(formula[[3]]) == exposure) == 0)
    stop("The exposure variable is not included in the formula.", call. = FALSE)

  #### Methods for non-matched or matched sampling designs ####
  n <- nrow(data)
  n.cases <- sum(data[, outcome])
  n.cluster <- length(unique(data[, clusterid]))
  data <- data[order(data[, clusterid]), ]

  # Create dataset data0 for counterfactual X = 0s
  data0 <- data
  data0[, exposure] <- 0
  clusters <- data[, clusterid]

  ## Design matrices ##
  design <- model.matrix(object = formula.model, data = data)[, - 1, drop = FALSE]
  design0 <- model.matrix(object = formula.model, data = data0)[, - 1, drop = FALSE]

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
  pred.diff <- getScoreResidualsFromClogit(fit = object,
                                           y = data[, outcome],
                                           x = design,
                                           id = clusters)

  if(missing(pred.diff)) warning("Use the latest version of package 'drgee'", call. = FALSE)
  score.beta <- pred.diff$U
  score.equations <- cbind(score.AF, score.beta)
  score.equations <- aggr(x = score.equations, clusters = clusters)
  meat <- var(score.equations, na.rm=TRUE)
  #### Bread: hessian of score equations ####
  ### Hessian of score equation 1 ##
  #### Estimating variance using Sandwich estimator ####
  ### Aggregate data ###
  hessian.AF1 <- - data[, outcome]
  hessian.AF1 <- aggr(x = hessian.AF1, clusters = clusters)
  hessian.AF2 <- cbind(as.matrix((design0 - design) * as.vector(data[, outcome] * exp( - log.or))))
  hessian.AF2 <- aggr(x = hessian.AF2, clusters = clusters)
  hessian.AF <-  cbind(mean(hessian.AF1), t(colMeans(hessian.AF2)))
  hessian.beta <- cbind(matrix(rep(0, npar), nrow = npar, ncol = 1), pred.diff$dU.sum / n)
  ### Bread ###
  bread <- rbind(hessian.AF, hessian.beta)

  #### Sandwich ####
  sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) * n.cluster/ n ^ 2 )

  AF.var <- sandwich[1, 1]

  #### Output ####
  out <- c(list(hessian.beta = hessian.beta, hessian.AF = hessian.AF, clusterid = clusterid,
                score.equations = score.equations, hessian.beta = hessian.beta, bread = bread,
                meat = meat, AF.est = AF.est, AF.var = AF.var, log.or = log.or,
                objectcall = objectcall, call = call, exposure = exposure, outcome = outcome,
                object = object, sandwich = sandwich, formula = formula, n = n, n.cases = n.cases,
                n.cluster = n.cluster))

  class(out) <- "AF"
  return(out)
}
