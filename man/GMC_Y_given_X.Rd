\name{GMC_Y_given_X}
\alias{GMC_Y_given_X}
\title{Generalized Measure of Correlation: GMC(Y | X)}
\usage{
GMC_Y_given_X(X, Y, kernel = dnorm)
}
\arguments{
\item{X}{Predictor variable}

\item{Y}{Response variable}

\item{kernel}{Kernel function (default = dnorm)}
}
\value{
GMC(Y|X) estimate
}
\description{
Generalized Measure of Correlation: GMC(Y | X)
}
\examples{
# Generate sample data with linear relationship
set.seed(123)
n <- 1000
X <- rnorm(n)
Y <- 2 * X + rnorm(n, sd = 0.5)

# Calculate GMC(Y|X)
gmc_result <- GMC_Y_given_X(X, Y)
print(gmc_result)
}
