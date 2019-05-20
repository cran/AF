AFgest <- function(object, data){
  call <- match.call()
  link <- object$input$link
  inputcall <- object$input$estmethod
  objectcall <- object$call
  # Warning if the object is not a gest object
  if(!(as.character(inputcall) == "g"))
    stop("The object is not the G-estimator", call. = FALSE)
  # Warning if the object is not a gest object with log or logit link
  if(!(link == "log" | link == "logit"))
    stop("The link is not log or logit", call. = FALSE)
  #######################
  # Preparation of object to fit format
  psi <- object$est
  X <- object$input$X
  Y <- object$input$Y
  #######################
  fitY <- object$input$fitY.LZX
  fitZ <- object$fitZ.L
  formulaZ <- fitZ$formula
  Z <- as.character(formulaZ[2])
  formula <- object$input$formula
  #### remove NA based on additional model
  temp <- model.matrix(object = formula, data = data)
  if(formula == ~1) rownames(temp) <- rownames(data)
  ## Delete rows with missing on variables in the model ##
  designpsi <- expand(temp, rownames(data))
  npsi <- ncol(designpsi)
  nZ <- length(coef(fitZ))
  n <- nrow(designpsi)
  n.cases <- sum(data[ , Y], na.rm = TRUE)
  if(object$converged == TRUE) {
    if(link == "log"){
      y0 <- data[, Y] * as.vector(exp( -(designpsi %*% psi) * data[, X]))
      P0.est <- mean(y0, na.rm = TRUE)
      P.est <- mean(data[, Y], na.rm = TRUE)
      AF.est <- 1 - P0.est / P.est 
      ### Score equations
      ## S(y0)
      score.y0 <- y0 - P0.est
      ## S(y)
      score.y <- data[, Y] - P.est
      ## Score functions of all parameters (y, y0, psi, E(z))
      score <- cbind(score.y, score.y0, object$estfun)
      ### Meat for the sandwich estimator
      meat <- var(score, na.rm = TRUE)
      ### Hessian
      I.y <- c(-1, 0, rep(0, npsi + nZ))
      dy0.dpsi <- colMeans(-designpsi * data[, X] * y0, na.rm = TRUE)
      I.y0 <- c(0, -1, dy0.dpsi, rep(0, nZ))
      ## Bread
      bread <- rbind(I.y, I.y0, cbind(matrix(0, ncol = 2, nrow = npsi + nZ), object$d.estfun))
      ## Variance
      sandwich <- (solve(bread) %*% meat %*% t(solve(bread)) / n)[1:2, 1:2]
      gradient <- as.matrix(c(P0.est / P.est ^ 2, - 1 / P.est), nrow = 2, ncol = 1)
      AF.var <- t(gradient) %*% sandwich %*% gradient
    }
    if(link == "logit"){
      nY <- length(coef(fitY))
      designY <- expand(model.matrix(object=fitY, data = data), rownames(data))
      linear_predY <- predict(object = fitY, newdata = data)
      ## Use only observations with non-missing for Y and X
      linear_predY <- linear_predY
      y0 <- plogis(linear_predY - (designpsi %*% psi) * data[, X])
      P0.est <- mean(y0, na.rm = TRUE)
      P.est <- mean(data[, Y], na.rm = TRUE)
      AF.est <- 1 - P0.est / P.est 
      ### Score equations
      ## S(y0)
      score.y0 <- y0 - mean(y0, na.rm = TRUE)
      ## S(y)
      score.y <- data[, Y] - mean(data[, Y], na.rm = TRUE)
      ## Score functions of all parameters (y, y0, psi, E(z), alpha)
      score <- cbind(score.y, score.y0, object$estfun)
      ### Meat for the sandwich estimator
      meat <- var(score, na.rm = TRUE)
      ### Hessian
      I.y <- c(-1, 0, rep(0, npsi + nZ + nY))
      dy0.dalpha <- colMeans(designY * as.vector((y0 * (1- y0))), na.rm = TRUE)
      dy0.dpsi <- colMeans(-(designpsi * data[, X]) * as.vector((y0 * (1- y0))), na.rm = TRUE)
      I.y0 <- c(0, -1, dy0.dpsi, rep(0, nZ), dy0.dalpha)
      ## Bread
      bread <- rbind(I.y, I.y0, cbind(matrix(0, ncol = 2, nrow = npsi + nZ + nY), object$d.estfun))
      ## Variance
      sandwich <- (solve(bread) %*% meat %*% t(solve(bread)) / n)[1:2, 1:2]
      gradient <- as.matrix(c(P0.est / P.est ^ 2, - 1 / P.est), nrow = 2, ncol = 1)
      AF.var <- t(gradient) %*% sandwich %*% gradient
    }
    #### Output
    out <- list(AF.est = AF.est, AF.var = AF.var, link = link, objectcall = objectcall, call = call, inputcall = inputcall,
                exposure = X, outcome = Y, n = n, n.cases = n.cases, n.cluster = 0, formula = formula,
                psi = psi, fitZ = fitZ, nZ = nZ, fitY = fitY, Z = Z)
  }
  else{
    #### Output
    out <- list(AF.est = NA, AF.var = NA, link = link, objectcall = objectcall, call = call, inputcall = inputcall,
                exposure = X, outcome = Y, n = n, n.cases = n.cases, n.cluster = 0, formula = formula,
                psi = psi, fitZ = fitZ, nZ = nZ, fitY = fitY, Z = Z)
  }
  
  class(out) <- "AF"
  
  return(out)
}
