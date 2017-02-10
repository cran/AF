############## AF function for a coxph object #####################
#' @title Attributable fraction function based on a Cox Proportional Hazard regression model as a \code{coxph} object (commonly used for cohort sampling designs with time-to-event outcomes).
#' @description \code{AFcoxph} estimates the model-based adjusted attributable fraction function from a Cox Proportional Hazard regression model in form of a  \code{\link[survival]{coxph}} object. This model is commonly used for data from cohort sampling designs with time-to-event outcomes.
#' @param object a fitted Cox Proportional Hazard regression model object of class "\code{\link[survival]{coxph}}". Method for handling ties must be breslow since this is assumed in the calculation of the standard errors. No special terms such as \code{cluster}, \code{strata} and \code{tt} is allowed in the formula for the fitted object.
#' @param data an optional data frame, list or environment (or object coercible by \code{as.data.frame} to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment (\code{formula}), typically the environment from which the function is called.
#' @param exposure the name of the exposure variable as a string. The exposure must be binary (0/1) where unexposed is coded as 0.
#' @param times a scalar or vector of time points specified by the user for which the attributable fraction function is estimated. If not specified the observed event times will be used.
#' @param clusterid the name of the cluster identifier variable as a string, if data are clustered. Cluster robust standard errors will be calculated.
#' @return \item{AF.est}{estimated attributable fraction function for every time point specified by \code{times}.}
#' @return \item{AF.var}{estimated variance of \code{AF.est}. The variance is obtained by combining the delta methods with the sandwich formula.}
#' @return \item{S.est}{estimated factual survival function; \eqn{S(t)}.}
#' @return \item{S.var}{estimated variance of \code{S.est}. The variance is obtained by the sandwich formula.}
#' @return \item{S0.est}{estimated counterfactual survival function if exposure would be eliminated; \eqn{S_0(t)}{S0(t)}.}
#' @return \item{S0.var}{estimated variance of \code{S0.est}. The variance is obtained by the sandwich formula.}
#' @details \code{AFcoxph} estimates the attributable fraction for a time-to-event outcome
#' under the hypothetical scenario where a binary exposure \code{X} is eliminated from the population. The estimate is adjusted for confounders \code{Z}
#' by the Cox proportional hazards model (\code{\link[survival]{coxph}}). Let the AF function be defined as
#' \deqn{AF=1-\frac{\{1-S_0(t)\}}{\{1-S(t)\}}}{AF = 1 - {1 - S0(t)} / {1 - S(t)}}
#' where \eqn{S_0(t)}{S0(t)} denotes the counterfactual survival function for the event if
#' the exposure would have been eliminated from the population at baseline and \eqn{S(t)} denotes the factual survival function.
#' If \code{Z} is sufficient for confounding control, then \eqn{S_0(t)}{S0(t)} can be expressed as \eqn{E_Z\{S(t\mid{X=0,Z })\}}{E_z{S(t|X=0,Z)}}.
#' The function uses a fitted Cox proportional hazards regression to estimate \eqn{S(t\mid{X=0,Z})}{S(t|X=0,Z)}, and the marginal sample distribution of \code{Z}
#' to approximate the outer expectation (\enc{Sjölander}{Sjolander} and Vansteelandt, 2014).  If \code{clusterid} is supplied, then a clustered sandwich formula is used in all variance calculations.
#' @author Elisabeth Dahlqwist, Arvid \enc{Sjölander}{Sjolander}
#' @seealso \code{\link[survival]{coxph}} and \code{\link[survival]{Surv}} used for fitting the Cox proportional hazards model.
#' @references Chen, L., Lin, D. Y., and Zeng, D. (2010). Attributable fraction functions for censored event times. \emph{Biometrika} \bold{97}, 713-726.
#' @references \enc{Sjölander}{Sjolander}, A. and Vansteelandt, S. (2014). Doubly robust estimation of attributable fractions in survival analysis. \emph{Statistical Methods in Medical Research}. doi: 10.1177/0962280214564003.
#' @examples
#' # Simulate a sample from a cohort sampling design with time-to-event outcome
#' expit <- function(x) 1 / (1 + exp( - x))
#' n <- 500
#' time <- c(seq(from = 0.2, to = 1, by = 0.2))
#' Z <- rnorm(n = n)
#' X <- rbinom(n = n, size = 1, prob = expit(Z))
#' Tim <- rexp(n = n, rate = exp(X + Z))
#' C <- rexp(n = n, rate = exp(X + Z))
#' Tobs <- pmin(Tim, C)
#' D <- as.numeric(Tobs < C)
#' #Ties created by rounding
#' Tobs <- round(Tobs, digits = 2)
#'
#' # Example 1: non clustered data from a cohort sampling design with time-to-event outcomes
#' data <- data.frame(Tobs, D, X,  Z)
#'
#' # Fit a Cox PH regression model
#' fit <- coxph(formula = Surv(Tobs, D) ~ X + Z + X * Z, data = data, ties="breslow")
#'
#' # Estimate the attributable fraction from the fitted Cox PH regression model
#' AFcoxph_est <- AFcoxph(fit, data=data, exposure ="X", times = time)
#' summary(AFcoxph_est)
#'
#' # Example 2: clustered data from a cohort sampling design with time-to-event outcomes
#' # Duplicate observations in order to create clustered data
#' id <- rep(1:n, 2)
#' data <- data.frame(Tobs = c(Tobs, Tobs), D = c(D, D), X = c(X, X), Z = c(Z, Z), id = id)
#'
#' # Fit a Cox PH regression model
#' fit <- coxph(formula = Surv(Tobs, D) ~ X + Z + X * Z, data = data, ties="breslow")
#'
#' # Estimate the attributable fraction from the fitted Cox PH regression model
#' AFcoxph_clust <- AFcoxph(object = fit, data = data,
#'                          exposure = "X", times = time, clusterid = "id")
#' summary(AFcoxph_clust)
#' plot(AFcoxph_clust, CI = TRUE)
#' 
#' # Estimate the attributable fraction from the fitted Cox PH regression model, time unspecified
#' AFcoxph_clust_no_time <- AFcoxph(object = fit, data = data,
#'                          exposure = "X", clusterid = "id")
#' summary(AFcoxph_clust_no_time)
#' plot(AFcoxph_clust, CI = TRUE)
#' @import survival data.table
#' @export
AFcoxph <- function(object, data, exposure, times, clusterid){
  call <- match.call()
  #### Preparation of dataset ####
  formula <- object$formula
  vars <- as.character(attr(terms(formula),"variables"))[-1]
  npar <- length(object$coef)
  # Warning if the object is not a glm object
  if(!(as.character(object$call[1]) == "coxph") | !is.null(object$userCall))
    stop("The object is not a coxph object", call. = FALSE)
  # Warning if specials are in the object formula
  specials <- pmatch(c("strata(","cluster(","tt("), attr(terms(formula),"variables"))
  if(any(!is.na(specials)))
    stop("No special terms are allowed in the formula")

  ## Delete rows with missing on variables in the model ##
  rownames(data) <- 1:nrow(data)
  m <- model.matrix(object = formula, data = data)
  complete <- as.numeric(rownames(m))
  data <- data[complete, ]
  
  ## Define object.detail
  object.detail <- coxph.detail(object = object)

  ## If times is missing ##
  if(missing(times))
    times <- object.detail$time
  ## Checks ##
  if(!object$method=="breslow")
    stop("Only breslow method for handling ties is allowed.", call. = FALSE)
  if(!is.binary(data[, exposure]))
    stop("Only binary exposure (0/1) is accepted.", call. = FALSE)
  if(max(all.vars(formula[[3]]) == exposure) == 0)
    stop("The exposure variable is not included in the formula.", call. = FALSE)
  if(missing(clusterid)) n.cluster <- 0
  else n.cluster <- length(unique(data[, clusterid]))

  ## Find names of end variable and event variable
  rr <- rownames(attr(terms(formula), "factors"))[1]
  temp <- gregexpr(", ", rr)[[1]]
  if(length(temp == 1)){
    endvar <- substr(rr, 6, temp[1] - 1)
    eventvar <- substr(rr, temp[1] + 2, nchar(rr) - 1)
  }
  if(length(temp) == 2){
    endvar <- substr(rr, temp[1] + 2, temp[2] - 1)
    eventvar <- substr(rr, temp[2] + 2, nchar(rr) - 1)
  }
  n <- nrow(data)
  n.cases <- sum(data[, eventvar])
  clusters <- data[, clusterid]
  npar <- length(object$coef)

  # Sort on "end-variable"
  ord <- order(data[, endvar])
  data <- data[ord, ]

  # Create dataset data0 for counterfactual X = 0
  data0 <- data
  data0[, exposure] <- 0

  ## Design matrices ##
  design <- as.matrix(model.matrix(object = delete.response(terms(object)), data = data)[, -1])
  design0 <- as.matrix(model.matrix(object = delete.response(terms(object)), data = data0)[, -1])

  ### Estimate the survival functions ###
  ## Hazard increment ##
  dH0 <- object.detail$hazard
  H0 <- cumsum(dH0)
  ## Baseline hazard function ##
  H0step <- stepfun(object.detail$time, c(0, H0))
  H0res <- rep(0, n)
  dH0.untied <- rep(dH0, object.detail$nevent) / rep(object.detail$nevent, object.detail$nevent)
  H0res[data[, eventvar] == 1] <- dH0.untied * n #handle ties

  ## Predict based on the Cox PH model ##
  epred <- predict(object = object, newdata = data, type = "risk")
  epred0 <- predict(object = object, newdata = data0, type = "risk")

  ### Meat ###
  ## Score equation 4 ## for the Cox PH model (made outside of loop)
  score.beta <- as.matrix(residuals(object = object, type = "score"))
  score.beta <- score.beta[ord, ]
  ## Weighted mean of the variable at event for all at risk at that time ##
  E <- matrix(0, nrow = n, ncol = npar)
  means <- as.matrix(object.detail$means)
  means <- means[rep(1:nrow(means), object.detail$nevent), ] #handle ties
  E[data[, eventvar] == 1, ] <- means
  ## One point and variance estimate for each time t in times ##
  S.est <- vector(length = length(times))
  S0.est <- vector(length = length(times))
  AF.var <- vector(length = length(times))
  S.var <- vector(length = length(times))
  S0.var <- vector(length = length(times))

  # Loop over all t in times
  for (i in 1:length(times)){
    t <- times[i]
    #### Meat: score equations ####
    ## Score equation 1 ## for the factual survival function
    score.S <- exp( - H0step(t) * epred)
    ## Score equation 2 ## for the counterfactual survival function
    score.S0 <- exp( - H0step(t) * epred0)
    ## Score equation 3 ##  for the breslow estimator
    score.H0 <- H0res * (data[, endvar] <= t)
    ## Score equation 4 ## for the Cox PH model (made outside of loop)

    ### Meat ###
    score.equations <- cbind(score.S, score.S0, score.H0, score.beta)
    if (!missing(clusterid)){
      score.equations <- score.equations
      score.equations <- aggr(score.equations, clusters = clusters)
    }

    meat <- var(score.equations, na.rm = TRUE)

    #### Bread: hessian of score equations ####
    ## Hessian of score equation 1 ##
    hessian.S <- c(-1, 0, mean(epred * score.S), colMeans(design * H0step(t) * epred * score.S))
    ## Hessian of score equation 2 ##
    hessian.S0 <- c(0, -1, mean(epred0 * score.S0), colMeans(design0 * H0step(t) * epred0 * score.S0))
    ## Hessian of score equation 3 ##
    hessian.H0 <- c(rep(0,2), - 1, - colMeans(E * score.H0, na.rm = TRUE))
    ## Hessian of score equation 4 ##
    hessian.beta <- cbind(matrix(0, nrow = npar, ncol = 3), - solve(vcov(object = object)) / n)
    ### Bread ###
    bread<-rbind(hessian.S, hessian.S0, hessian.H0, hessian.beta)
    ### Sandwich ###
    if (!missing(clusterid))
      sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) * n.cluster/ n^2 ) [1:2, 1:2]
    else
      sandwich <- (solve (bread) %*% meat %*% t(solve (bread)) / n) [1:2, 1:2]
    #### For point estimate ####
    S.est[i] <- mean(x = score.S, na.rm = TRUE)
    S0.est[i] <- mean(x = score.S0, na.rm = TRUE)
    #### Estimate of variance using the delta method ####
    gradient <- as.matrix(c( - (1 - S0.est[i]) / (1 - S.est[i]) ^ 2, 1 / (1 - S.est[i]))
                          , nrow = 2, ncol = 1)
    AF.var[i] <- t(gradient) %*% sandwich %*% gradient
    S.var[i] <- sandwich[1, 1]
    S0.var[i] <- sandwich[2, 2]
  }

  ### The AF function estimate ###
  AF.est <- 1 - (1 - S0.est) / (1 - S.est)

  #### Output ####
  #func <- AF.cc
  out <- c(list(AF.est = AF.est, AF.var = AF.var, S.est = S.est,
                S0.est = S0.est, S.var = S.var, S0.var = S0.var,
                objectcall = object$call, call = call, exposure = exposure, outcome = eventvar,
                object = object, sandwich = sandwich, gradient = gradient, formula = formula,
                n = n, n.cases = n.cases, n.cluster = n.cluster,  times = times))
  class(out) <- "AF"
  return(out)
}
