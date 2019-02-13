#' Radial basis function model
#'
#' @param X Matrix witn inputs
#' @param Y Matrix with outputs
#' @param augmented Should the augmented RBF model be used?
#'
#' @return RBF model
#' @export
#'
#' @examples
#' d <- 8
#' n <- 300
#' x <- matrix(runif(d*n), ncol=d)
#' f <- function(x) {cos(x[1]^1.9)}
#' y <- as.matrix(apply(x, 1, f))
#' arb1 <- rbf(X=x, Y = y)
#' np <- 300
#' xp <- matrix(runif(d*np), ncol=d)
#' p1 <- predict(arb1, xp)
rbf <- function(X, Y, augmented=TRUE) {#browser()
  self <- list()
  class(self) <- c("rbf", "list")
  # centers <- X
  n <- nrow(X)
  d <- ncol(X)
  self$BasisFunc <- function(x,y) {
    sqrt(sum((x - y)^2))^3
    # exp(-sum((x - y)^2)/2/.5^2)
  }
  Phi <- outer(1:n, 1:n, Vectorize(
    function(i,j) {
      # sqrt(abs((X[i,] - X[j,])^2))
      self$BasisFunc(X[i,], X[j,])
    }
  ))
  if (augmented) {
  k <- 2
  P <- matrix(0, n, k*d+1)
  P[,1] <- 1
  # browser()
  for (ki in 1:k) {
    P[,1+ki*d - d + 1:d] <- X ^ ki
  }
  AugPhi <- rbind(cbind(Phi, P), cbind(t(P), matrix(0,k*d+1,k*d+1)))
  AugPhiinv_Y <- solve(AugPhi, rbind(Y, matrix(0,k*d+1,1)))
  theta <- AugPhiinv_Y[1:n]
  mu <- AugPhiinv_Y[(n+1):(n+1 + k*d)]
  } else {
    theta <- solve(Phi, Y)
    mu <- NULL
  }

  self$X <- X
  self$Y <- Y
  self$k <- k
  self$AugPhiinv_Y <- AugPhiinv_Y
  self$Phi <- Phi
  self$theta <- theta
  self$mu <- mu
  self$augmented <- augmented
  self
}

#' Predict with RBF
#'
#' @param object RBF model
#' @param xp Matrix of points to predict at
#' @param ... Further arguments passed to predict
#'
#' @return Predictions at xp
#' @export
#' @rdname rbf
predict.rbf <- function(object, xp, ...) {#browser()
  self <- object
  print("using augmented")
  n <- nrow(self$X)
  d <- ncol(self$X)
  np <- nrow(xp)
  phi <- outer(1:n, 1:np, Vectorize(function(i, ip) {#browser()
    # sqrt(abs((self$X[i,] - xp[ip,])^2))
    self$BasisFunc(self$X[i,], xp[ip,])
  }))
  # browser()
  Ypredmean <- t(phi) %*% self$theta
  if (self$augmented) {
    Ypredmean <- Ypredmean + rowSums(matrix(self$mu, nrow(xp), d*self$k+1, byrow=T) * do.call(cbind, lapply(0:self$k, function(ki) {if (ki==0) {matrix(1,nrow(xp),1)} else {xp^ki}})))
  }

  # Get var
  Ypredvar <- self$BasisFunc(self$X[1,], self$X[1,]) - t(phi) %*% solve(self$Phi, phi)

  # Return mean and var
  list(mean=Ypredmean, var=diag(Ypredvar), cov=Ypredvar)
}
