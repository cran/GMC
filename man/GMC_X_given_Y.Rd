\name{GMC_X_given_Y}
\alias{GMC_X_given_Y}
\title{Generalized Measure of Correlation: GMC(X | Y)}
\usage{
GMC_X_given_Y(X, Y, kernel = dnorm)
}
\arguments{
\item{X}{Predictor variable}

\item{Y}{Response variable}

\item{kernel}{Kernel function (default = dnorm)}
}
\value{
GMC(X|Y) estimate
}
\description{
Generalized Measure of Correlation: GMC(X | Y)
}
\examples{
# Generate sample data with nonlinear relationship
set.seed(123)
n <- 1000
X <- rnorm(n)
Y <- X^2 + rnorm(n, sd = 0.5)

# Calculate GMC(X|Y)
gmc_result <- GMC_X_given_Y(X, Y)
print(gmc_result)
}
