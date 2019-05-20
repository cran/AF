############## AF function for a parfrailty object #####################
#' @title Attributable fraction function based on a Weibull gamma-frailty model as a \code{\link[stdReg]{parfrailty}} object (commonly used for cohort sampling family designs with time-to-event outcomes).
#' @description \code{AFparfrailty} estimates the model-based adjusted attributable fraction function from a shared Weibull gamma-frailty model in form of a  \code{\link[stdReg]{parfrailty}} object. This model is commonly used for data from cohort sampling familty designs with time-to-event outcomes.
#' @param object a fitted Weibull gamma-parfrailty object of class "\code{\link[stdReg]{parfrailty}}".
#' @param data an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment (\code{formula}), typically the environment from which the function is called.
#' @param exposure the name of the exposure variable as a string. The exposure must be binary (0/1) where unexposed is coded as 0.
#' @param times a scalar or vector of time points specified by the user for which the attributable fraction function is estimated. If not specified the observed death times will be used.
#' @param clusterid the name of the cluster identifier variable as a string, if data are clustered.
#' @return \item{AF.est}{estimated attributable fraction function for every time point specified by \code{times}.}
#' @return \item{AF.var}{estimated variance of \code{AF.est}. The variance is obtained by combining the delta methods with the sandwich formula.}
#' @return \item{S.est}{estimated factual survival function; \eqn{S(t)}.}
#' @return \item{S.var}{estimated variance of \code{S.est}. The variance is obtained by the sandwich formula.}
#' @return \item{S0.est}{estimated counterfactual survival function if exposure would be eliminated; \eqn{S_0(t)}{S0(t)}.}
#' @return \item{S0.var}{estimated variance of \code{S0.est}. The variance is obtained by the sandwich formula.}
#' @details \code{AFparfrailty} estimates the attributable fraction for a time-to-event outcome
#' under the hypothetical scenario where a binary exposure \code{X} is eliminated from the population.
#' The estimate is adjusted for confounders \code{Z} by the shared frailty model (\code{\link[stdReg]{parfrailty}}).
#' The baseline hazard is assumed to follow a Weibull distribution and the unobserved shared frailty effects \code{U} are assumed to be gamma distributed.
#' Let the AF function be defined as
#' \deqn{AF=1-\frac{\{1-S_0(t)\}}{\{1-S(t)\}}}{AF = 1 - {1 - S0(t)} / {1 - S(t)}}
#' where \eqn{S_0(t)}{S0(t)} denotes the counterfactual survival function for the event if
#' the exposure would have been eliminated from the population at baseline and \eqn{S(t)} denotes the factual survival function.
#' If \code{Z} and \code{U} are sufficient for confounding control, then \eqn{S_0(t)}{S0(t)} can be expressed as \eqn{E_Z\{S(t\mid{X=0,Z })\}}{E_z{S(t|X=0,Z)}}.
#' The function uses a fitted Weibull gamma-frailty model to estimate \eqn{S(t\mid{X=0,Z})}{S(t|X=0,Z)}, and the marginal sample distribution of \code{Z}
#' to approximate the outer expectation. A clustered sandwich formula is used in all variance calculations.
#' @author Elisabeth Dahlqwist, Arvid \enc{Sjölander}{Sjolander}
#' @seealso \code{\link[stdReg]{parfrailty}} used for fitting the Weibull gamma-frailty and \code{\link[stdReg]{stdParfrailty}} used for standardization of a \code{parfrailty} object.
#' @examples
#'# Example 1: clustered data with frailty U
#' expit <- function(x) 1 / (1 + exp( - x))
#' n <- 100
#' m <- 2
#' alpha <- 1.5
#' eta <- 1
#' phi <- 0.5
#' beta <- 1
#' id <- rep(1:n,each=m)
#' U <- rep(rgamma(n, shape = 1 / phi, scale = phi), each = m)
#' Z <- rnorm(n * m)
#' X <- rbinom(n * m, size = 1, prob = expit(Z))
#' # Reparametrize scale as in rweibull function
#' weibull.scale <- alpha / (U * exp(beta * X)) ^ (1 / eta)
#' t <- rweibull(n * m, shape = eta, scale = weibull.scale)
#'
#' # Right censoring
#' c <- runif(n * m, 0, 10)
#' delta <- as.numeric(t < c)
#' t <- pmin(t, c)
#'
#' data <- data.frame(t, delta, X, Z, id)
#'
#' # Fit a parfrailty object
#' library(stdReg)
#' fit <- parfrailty(formula = Surv(t, delta) ~ X + Z + X * Z, data = data, clusterid = "id")
#' summary(fit)
#'
#' # Estimate the attributable fraction from the fitted frailty model
#'
#' time <- c(seq(from = 0.2, to = 1, by = 0.2))
#'
#' AFparfrailty_est <- AFparfrailty(object = fit, data = data, exposure = "X",
#'                                   times = time, clusterid = "id")
#' summary(AFparfrailty_est)
#' plot(AFparfrailty_est, CI = TRUE, ylim=c(0.1,0.7))
#' @import survival data.table stdReg
#' @importFrom stats model.extract model.frame
#' @export
AFparfrailty <- function(object, data, exposure, times, clusterid){
  call <- match.call()
  formula <- object$formula
  npar <- length(object$est)

  ## Delete rows with missing on variables in the model ##
  rownames(data) <- 1:nrow(data)
  m <- model.frame(formula, data = data)
  complete <- as.numeric(rownames(m))
  data <- data[complete, ]

  ## Find names of outcome
  rr <- rownames(attr(terms(formula), "factors"))[1]
  temp <- gregexpr(", ", rr)[[1]]
  if(length(temp == 1)){
    outcome <- substr(rr, temp[1] + 2, nchar(rr) - 1)
  }
  if(length(temp) == 2){
    outcome <- substr(rr, temp[2] + 2, nchar(rr) - 1)
  }

  ## Define end variable and event variable
  Y <- model.extract(frame = m, "response")
  if(ncol(Y) == 2){
    endvar <- Y[, 1]
    eventvar <- Y[, 2]
  }
  if(ncol(Y) == 3){
    endvar <- Y[, 2]
    eventvar <- Y[, 3]
  }

  ## Defining parameters and variables##
  logalpha <- object$est[1]
  alpha <- exp(logalpha)
  logeta <- object$est[2]
  eta <- exp(logeta)
  logphi <- object$est[3]
  phi <- exp(logphi)
  beta <- object$est[(3 + 1):npar]

  # Assign value to t if missing
  if(missing(times)){
      times <- endvar[eventvar == 1]
  }

  times <- sort(times)

  n <- nrow(data)
  n.cases <- sum(eventvar)
  n.cluster <- object$ncluster

  ## Counterfactual dataset ##
  data0 <- data
  data0[, exposure] <- 0

  ## Design matrices ##
  design <- model.matrix(object = formula, data = data)[, -1, drop=FALSE]
  design0 <- model.matrix(object = formula, data = data0)[, -1, drop=FALSE]
  clusters <-  data[, clusterid]

  ### Estimate the survival functions ###
  ########### order of beta has to be the same as order of design matrix, not fixed
  predX <- design %*% beta
  pred0 <- design0 %*% beta

  ## One point and variance estimate for each time t in times ##
  S.est <- vector(length = length(times))
  S0.est <- vector(length = length(times))
  AF.var <- vector(length = length(times))
  S.var <- vector(length = length(times))
  S0.var <- vector(length = length(times))

  # Loop over all t in times
  for (i in 1:length(times)){

    t <- times[i]
    H0t <- (t / alpha) ^ eta

    ### Survival functions
    temp <- 1 + phi * H0t * exp(predX)
    temp0 <- 1 + phi * H0t * exp(pred0)
    surv <- temp ^ ( - 1 / phi)
    surv0 <- temp0 ^ ( - 1 / phi)
    S.est[i] <- mean(surv, na.rm = TRUE)
    S0.est[i] <- mean(surv0, na.rm = TRUE)

    ## Score functions
    sres <- surv - S.est[i]
    sres0 <- surv0 - S0.est[i]
    Scores.S <- cbind(sres, sres0)
    Scores.S <- aggr(x = Scores.S, clusters = clusters)
    coefres <- object$score
    res <- cbind(Scores.S, coefres)

    meat <- var(res, na.rm = TRUE)

    ### Hessian for the factual survival function
    dS.dlogalpha <- sum(H0t * eta * exp(predX) / temp ^ (1 / phi + 1)) / n.cluster
    dS.dlogeta <- sum(-H0t * exp(predX) * log(t / alpha) * eta / temp ^ (1 / phi + 1)) / n.cluster
    dS.dlogphi <- sum(log(temp) / (phi * temp ^ (1 / phi)) - H0t * exp(predX) / temp ^ (1 / phi + 1)) / n.cluster
    dS.dbeta <- colSums(-H0t * as.vector(exp(predX)) * design / as.vector(temp) ^ (1 / phi + 1)) / n.cluster

    ### Hessian for the counterfactual survival function
    dS0.dlogalpha <- sum(H0t * eta * exp(pred0) / temp0 ^ (1 / phi + 1)) / n.cluster
    dS0.dlogeta <- sum( - H0t * exp(pred0) * log(t / alpha) * eta / temp0 ^ (1 / phi + 1)) /n.cluster
    dS0.dlogphi <- sum(log(temp0) / (phi * temp0 ^ (1 / phi)) - H0t * exp(pred0) / temp0 ^ (1 / phi + 1)) / n.cluster
    dS0.dbeta <- colSums(-H0t * as.vector(exp(pred0)) * design0 / as.vector(temp0) ^ (1 / phi + 1)) / n.cluster

    #Note: the term n/n.cluster is because SI.logalpha, SI.logeta, SI.logphi,
    #and SI.beta are clustered, which they are not in stdCoxph
    S.hessian <- cbind(-diag(2) * n / n.cluster, rbind(dS.dlogalpha, dS0.dlogalpha), rbind(dS.dlogeta, dS0.dlogeta),
                rbind(dS.dlogphi, dS0.dlogphi), rbind(dS.dbeta, dS0.dbeta))

    par.hessian <- cbind(matrix(0, nrow = npar, ncol = 2), -solve(object$vcov) / n.cluster)
    bread <- rbind(S.hessian, par.hessian)
    sandwich <- (solve(bread) %*% meat %*% t(solve(bread)) / n.cluster)[1:2, 1:2]

    #### Estimate of variance using the delta method ####
    gradient <- as.matrix(c( - (1 - S0.est[i]) / (1 - S.est[i]) ^ 2, 1 / (1 - S.est[i])), nrow = 2, ncol = 1)
    AF.var[i] <- t(gradient) %*% sandwich %*% gradient
    S.var[i] <- sandwich[1, 1]
    S0.var[i] <- sandwich[2, 2]
    }

  ### The AF function estimate ###
  AF.est <- 1 - (1 - S0.est) / (1 - S.est)

  out <- c(list(AF.est = AF.est, AF.var = AF.var, S.est = S.est,
                S0.est = S0.est, S.var = S.var, S0.var = S0.var,
                objectcall = object$call, call = call, exposure = exposure, outcome = outcome, object = object,
                sandwich = sandwich, gradient = gradient, formula = formula,
                n = n, n.cases = n.cases, n.cluster = n.cluster,  times = times))
  class(out) <- "AF"
  return(out)

}


