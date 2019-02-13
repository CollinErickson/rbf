test_that("Basics work", {
  d <- 3
  n <- 100
  X <- matrix(runif(n*d), ncol=3)
  f <- function(x) {x[1]*x[2]^1.2 + exp(-x[3]^2)}
  Y <- matrix(apply(X, 1, f))

  np <- 200
  Xp <- matrix(runif(np*d), ncol=3)
  Yp <- matrix(apply(Xp, 1, f))

  for (use_augmented in c(T, F)) {
    for (kerneli in c("cubic", "gauss", "multiquadratic", "inversequadratic")) {
      expect_error(r1 <- rbf(X, Y, augmented=use_augmented, kernel=kerneli), NA, info = paste(use_augmented, kerneli))
      expect_is(r1, "rbf")
      expect_is(r1, "list")
      expect_error(p1 <- predict(r1, Xp), NA)
      if (F) {
        cat(round(unlist(SGGP::valstats(p1$mean, p1$var, Yp)),4), "\t", sep='\t')
        cat(use_augmented, '\t', kerneli, '\n')
      }
    }
  }
})
