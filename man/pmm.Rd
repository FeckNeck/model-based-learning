\name{em_poisson}
\alias{em_poisson}
\title{Em algorithm for multivariate Poisson distribution}
\usage{
em_poisson(X, K=2, eps=0.02)
}
\description{
Calculates the maximum likelihood estimates of the parameters of a multivariate Poisson distribution using the EM algorithm.
}
\arguments{
  \item{X}{a matrix of poisson distributed data}
  \item{K}{the number of clusters}
  \item{eps}{the convergence criterion}
}
\value{
  \item{Lambdas}{the estimated means of the clusters}
  \item{pk}{proportions of observations for each clusters}
  \item{clusters}{The affected cluster for each observation}
  \item{tk}{probabilities of each observation to belong to each cluster}
}
\examples{
n <- 200  # Number of rows
lambda1 <- 3  # Poisson parameter for column 1
lambda2 <- 5  # Poisson parameter for column 2
lambda3 <- 7  # Poisson parameter for column 3

# Generate random Poisson-distributed values
x1 <- rpois(n, lambda1)
x2 <- rpois(n, lambda2)
x3 <- rpois(n, lambda2)

# Create the data frame
X <- data.frame(x1, x2,x3)

####EXAMPLE: EM-based Clustering
res <- em_poisson(X)

print(res)
}
