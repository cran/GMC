\name{GMC_feature_ranking}
\alias{GMC_feature_ranking}
\title{Feature selection using GMC ranking}
\usage{
GMC_feature_ranking(X, Y, kernel = dnorm, sort = TRUE)
}
\arguments{
\item{X}{A matrix or data.frame of predictors}

\item{Y}{A numeric response vector}

\item{kernel}{Kernel function (default = dnorm)}

\item{sort}{Logical, whether to sort variables by GMC score}
}
\value{
A data.frame with variable names and GMC scores
}
\description{
Feature selection using GMC ranking
}
\examples{
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
}
