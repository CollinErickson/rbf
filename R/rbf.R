#' Radial basis function model
#'
#' @param X Matrix witn inputs
#' @param Y Matrix with outputs
#' @param augmented Should the augmented RBF model be used?
#' @param kernel Which kernel should be used
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
rbf <- function(X, Y, augmented=TRUE, kernel="cubic") {#browser()
  if (!is.matrix(X)) {stop("X must be matrix")}
  if (!is.matrix(Y)) {
    if (is.numeric(Y) && length(Y) == nrow(X)) {
      Y <- as.matrix(Y)
    } else {stop("Y is wrong length")}
  }

  self <- list()
  class(self) <- c("rbf", "list")
  # centers <- X
  n <- nrow(X)
  d <- ncol(X)
  self$BasisFunc <-
    if (kernel %in% c("cubic")) {
      function(x,y) {
        sqrt(sum((x - y)^2))^3
      }
    } else if (kernel %in% c("gauss")) {
      function(x,y) {
        exp(-sum((x - y)^2)/2/.2^2)
      }
    # } else if (kernel %in% c("tsp")) { # At 0 it is -Inf
    #   function(x,y) {
    #     r2 <- sum((x - y)^2)
    #     r2 * log(r2) / 2
    #   }
    } else if (kernel %in% c("multiquadratic")) {
      function(x,y) {
        # sqrt(1 - min(1,(sum((x - y)^2)) / 3))
        sqrt(1 + (sum((x - y)^2)) / 3)
      }
    } else if (kernel %in% c("inversequadratic")) {
      function(x,y) {
        1/(1 - min(1,(sum((x - y)^2)) / 3))
      }
    } else {stop("kernel no good")}
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
    self$k <- k
    self$AugPhiinv_Y <- AugPhiinv_Y
    self$mu <- mu
  } else {
    theta <- solve(Phi, Y)
    # mu <- NULL
  }

  self$X <- X
  self$Y <- Y
  self$Phi <- Phi
  self$theta <- theta
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
predict.rbf <- function(object, xp, ...) {
  self <- object
  n <- nrow(self$X)
  d <- ncol(self$X)
  np <- nrow(xp)
  phi <- outer(1:n, 1:np, Vectorize(function(i, ip) {
    self$BasisFunc(self$X[i,], xp[ip,])
  }))

  Ypredmean <- t(phi) %*% self$theta
  if (self$augmented) {
    Ypredmean <- Ypredmean + rowSums(matrix(self$mu, nrow(xp), d*self$k+1, byrow=T) * do.call(cbind, lapply(0:self$k, function(ki) {if (ki==0) {matrix(1,nrow(xp),1)} else {xp^ki}})))
  }

  # Get var
  Ypredvar <- self$BasisFunc(self$X[1,], self$X[1,]) - t(phi) %*% solve(self$Phi, phi)
  Ypredvar <- pmax(Ypredvar, 1e-16)

  # Return mean and var
  list(mean=Ypredmean, var=diag(Ypredvar), cov=Ypredvar)
}
