test_that("GMC_Y_given_X works correctly", {
  set.seed(123)
  n <- 100
  X <- rnorm(n)
  Y <- 2 * X + rnorm(n, sd = 0.5)
  
  result <- GMC_Y_given_X(X, Y)
  
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0)
  expect_true(result <= 1)
})

test_that("GMC_X_given_Y works correctly", {
  set.seed(123)
  n <- 100
  X <- rnorm(n)
  Y <- X^2 + rnorm(n, sd = 0.5)
  
  result <- GMC_X_given_Y(X, Y)
  
  expect_type(result, "double")
  expect_length(result, 1)
  expect_true(result >= 0)
})

test_that("GMC_feature_ranking works correctly", {
  set.seed(123)
  n <- 100
  X1 <- rnorm(n)
  X2 <- rnorm(n)
  Y <- 2 * X1 + rnorm(n, sd = 0.5)
  X <- cbind(X1, X2)
  
  result <- GMC_feature_ranking(X, Y)
  
  expect_s3_class(result, "data.frame")
  expect_equal(ncol(result), 2)
  expect_equal(nrow(result), 2)
  expect_true(all(c("Variable", "GMC") %in% names(result)))
})
