#' Estimate E[(E[Y|X])^2] using kernel regression
#'
#' This function estimates the squared conditional expectation E[(E[Y|X])^2] using
#' Nadaraya-Watson regression with Gaussian kernel.
#'
#' @param X A numeric vector of predictors.
#' @param Y A numeric vector of responses.
#' @param grid_length Number of grid points for numerical integration (default = 10000).
#' @param kernel Kernel function (default is dnorm).
#'
#' @return A list containing:
#' \describe{
#'   \item{estimate}{Estimated value of E[(E[Y|X])^2]}
#'   \item{bandwidth}{Selected kernel bandwidth}
#'   \item{mean_Y}{Mean of Y}
#'   \item{var_Y}{Variance of Y}
#'   \item{EY_grid}{Grid values of E[Y|X]}
#'   \item{fx_grid}{Estimated marginal density of X}
#'   \item{x_grid}{Grid points used in estimation}
#' }
#' @keywords internal
#' @importFrom ks hscv
#' @importFrom stats dnorm var
#' @references Zheng, S., Shi, N.Z., & Zhang, Z. (2012).
#'   Generalized Measures of Correlation for Asymmetry, Nonlinearity, and Beyond.
#'   Journal of the American Statistical Association, 107(499), 1239-1252.
#'   \doi{10.1080/01621459.2012.710509}

estimate_EY_X_squared <- function(X, Y, grid_length = 10000, kernel = dnorm) {

  # Ensure package is available
  if (!requireNamespace("ks", quietly = TRUE)) {
    stop("Package 'ks' is required. Install it using install.packages('ks')")
  }

  # Convert to matrix (required for ks::hscv)
  X_matrix <- matrix(X, ncol = 1)

  # Step 1: Bandwidth selection via plug-in method
  h_opt <- ks::hscv(x = X_matrix)

  # Step 2: Nadaraya-Watson conditional expectation estimator
  compute_EY_X <- function(X_query, X_data, Y_data, h, K = dnorm) {
    scaled_diff <- outer(X_query, X_data, function(x, xi) (x - xi)/h)
    Kx <- K(scaled_diff) / h
    epsilon <- .Machine$double.eps^0.5
    weights <- Kx / (rowSums(Kx) + epsilon)
    drop(weights %*% Y_data)
  }

  # Step 3: Estimate marginal density f_X(x)
  estimate_fx <- function(x_grid, X_data, h, K = dnorm) {
    scaled_diff <- outer(x_grid, X_data, function(x, xi) (x - xi)/h)
    K_vals <- K(scaled_diff) / h
    rowMeans(K_vals)
  }

  # Step 4: Construct grid for integration
  h <- h_opt
  x_grid <- seq(min(X) - 3 * h, max(X) + 3 * h, length.out = grid_length)
  delta <- diff(x_grid[1:2])

  # Step 5: Compute E[Y|X=x] and f_X(x)
  EY_grid <- compute_EY_X(x_grid, X, Y, h, kernel)
  fx_grid <- estimate_fx(x_grid, X, h, kernel)

  # Step 6: Integral approximation of E[(E[Y|X])^2]
  estimate <- sum(EY_grid^2 * fx_grid) * delta

  return(list(
    estimate = estimate,
    bandwidth = h,
    mean_Y = mean(Y),
    var_Y = var(Y),
    EY_grid = EY_grid,
    fx_grid = fx_grid,
    x_grid = x_grid
  ))
}

#' Generalized Measure of Correlation: GMC(Y | X)
#'
#' @param X Predictor variable
#' @param Y Response variable
#' @param kernel Kernel function (default = dnorm)
#' @importFrom stats dnorm var
#' @return GMC(Y|X) estimate
#' @export
#' @examples
#' # Generate sample data with linear relationship
#' set.seed(123)
#' n <- 1000
#' X <- rnorm(n)
#' Y <- 2 * X + rnorm(n, sd = 0.5)
#' 
#' # Calculate GMC(Y|X)
#' gmc_result <- GMC_Y_given_X(X, Y)
#' print(gmc_result)
GMC_Y_given_X <- function(X, Y, kernel = dnorm) {

  est <- estimate_EY_X_squared(X, Y, kernel = kernel)
  h <- est$bandwidth

  # # Guass kernel E[K] and var(K)
  EK <- 0           # ∫ z K(z) dz = 0
  varK <- 1         # ∫ z^2 K(z) dz = 1

  numerator <- est$estimate - (mean(Y) + h * EK)^2
  denominator <- var(Y) + h^2 * varK

  GMC <- numerator / denominator
  return(GMC)
}

#' Generalized Measure of Correlation: GMC(X | Y)
#'
#' @param X Predictor variable
#' @param Y Response variable
#' @param kernel Kernel function (default = dnorm)
#' @importFrom stats dnorm var
#' @return GMC(X|Y) estimate
#' @export
#' @examples
#' # Generate sample data with nonlinear relationship
#' set.seed(123)
#' n <- 1000
#' X <- rnorm(n)
#' Y <- X^2 + rnorm(n, sd = 0.5)
#' 
#' # Calculate GMC(X|Y)
#' gmc_result <- GMC_X_given_Y(X, Y)
#' print(gmc_result)
GMC_X_given_Y <- function(X, Y, kernel = dnorm) {

  est <- estimate_EY_X_squared(Y, X, kernel = kernel)
  h <- est$bandwidth

  # Guass kernel E[K] and var(K)
  EK <- 0
  varK <- 1

  numerator <- est$estimate - (mean(X) + h * EK)^2
  denominator <- var(X) + h^2 * varK

  GMC <- numerator / denominator
  return(GMC)
}

#' Feature selection using GMC ranking
#'
#' @param X A matrix or data.frame of predictors
#' @param Y A numeric response vector
#' @param kernel Kernel function (default = dnorm)
#' @param sort Logical, whether to sort variables by GMC score
#' @importFrom stats dnorm var
#' @return A data.frame with variable names and GMC scores
#' @export
#' @examples
#' # Generate sample data with multiple predictors
#' set.seed(123)
#' n <- 500
#' X1 <- rnorm(n)
#' X2 <- rnorm(n)
#' X3 <- rnorm(n)
#' Y <- 2 * X1 + X2^2 + rnorm(n, sd = 0.5)
#' X <- cbind(X1, X2, X3)
#' 
#' # Rank features by GMC
#' ranking <- GMC_feature_ranking(X, Y)
#' print(ranking)
GMC_feature_ranking <- function(X, Y, kernel = dnorm, sort = TRUE) {
  #  secret_env$env_lock()
  if (is.vector(X)) X <- matrix(X, ncol = 1)
  p <- ncol(X)
  gmc_scores <- numeric(p)
  for (j in 1:p) {
    gmc_scores[j] <- GMC_Y_given_X(X[, j], Y, kernel = kernel)
  }
  result <- data.frame(
    Variable = paste0("X", 1:p),
    GMC = gmc_scores
  )
  if (sort) {
    result <- result[order(result$GMC, decreasing = TRUE), ]
    rownames(result) <- NULL
  }
  return(result)
}
