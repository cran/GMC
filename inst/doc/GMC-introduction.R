## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(GMC)

## -----------------------------------------------------------------------------
# Generate sample data with linear relationship
set.seed(123)
n <- 1000
X <- rnorm(n)
Y <- 2 * X + rnorm(n, sd = 0.5)

# Calculate GMC(Y|X)
gmc_y_given_x <- GMC_Y_given_X(X, Y)
print(paste("GMC(Y|X) =", round(gmc_y_given_x, 3)))

## -----------------------------------------------------------------------------
# Generate sample data with nonlinear relationship
set.seed(123)
X <- rnorm(n)
Y <- X^2 + rnorm(n, sd = 0.5)

# Calculate GMC(X|Y)
gmc_x_given_y <- GMC_X_given_Y(X, Y)
print(paste("GMC(X|Y) =", round(gmc_x_given_y, 3)))

## -----------------------------------------------------------------------------
# Generate sample data with multiple predictors
set.seed(123)
n <- 500
X1 <- rnorm(n)
X2 <- rnorm(n)
X3 <- rnorm(n)
Y <- 2 * X1 + X2^2 + rnorm(n, sd = 0.5)
X <- cbind(X1, X2, X3)

# Rank features by GMC
ranking <- GMC_feature_ranking(X, Y)
print(ranking)

