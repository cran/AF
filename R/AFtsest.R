AFtsest <- function(object, data){
  call <- match.call()
  inputcall <- object$input$estmethod
  objectcall <- object$call
  # Warning if the object is not a tsest object
  if(!(as.character(inputcall) == "ts"))
    stop("The object is not a TS estimate", call. = FALSE)
  psi <- object$est[2]
  fitY <- object$input$fitY.LX
  fitX <- object$input$fitX.LZ
  link <- family(fitY)$link
  Y <- as.character(fitY$formula)[2]
  X <- as.character(fitX$formula)[2]
  ################
  nX <- length(coef(fitX))
  nY <- length(coef(fitY))
  n <- nrow(data)
  n.cases <- sum(data[ , Y], na.rm = TRUE)
  ########### AF estimate #############
  y0 <- data[, Y] * as.vector(exp( -psi * data[, X]))
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
  I.y <- c(-1, 0, rep(0, nY + nX))
  dy0.dpsi <- mean(data[, X] * y0, na.rm = TRUE)
  I.y0 <- c(0, -1, rep(0, nX + nY - 1), dy0.dpsi)
  ## Bread
  bread <- rbind(I.y, I.y0, cbind(matrix(0, ncol = 2, nrow = nX + nY), object$d.estfun))
  ## Variance
  sandwich <- (solve(bread) %*% meat %*% t(solve(bread)) / n)[1:2, 1:2]
  gradient <- as.matrix(c(P0.est / P.est ^ 2, - 1 / P.est), nrow = 2, ncol = 1)
  AF.var <- t(gradient) %*% sandwich %*% gradient
  
  #### Output
  out <- list(AF.est = AF.est, AF.var = AF.var, link = link, objectcall = objectcall, call = call, inputcall = inputcall,
              exposure = X, outcome = Y, n = n, n.cases = n.cases, n.cluster = 0, psi = psi, nX = nX, fitY = fitY, fitX = fitX)
  
  class(out) <- "AF"
  
  return(out)
}