# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Install Package:           'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

em_Poisson <- function(X, k=2, eps=.Machine$double.eps){

  N <- nrow(X) ## number of observations
  error <- eps + 1 #
  iter <- 1

  ## random initialization
  Lambdas <- runif(ncol(X), 1, 20)
  random_pks <- runif(k, 0, 1)
  pks <- random_pks / sum(random_pks)

  while(error > eps){

    tks <- c()
    d.k.s <- 1

    for(i in 1:k){
      d.k = 1
      for(l in 1:ncol(X)){
        d.l <- dpois(X[,l], Lambdas[l])
        d.k <- d.k * d.l
      }
      d.k.s <- d.k.s + d.k
    }

    for(i in 1:k){
      qk <- 1
      for(l in 1:ncol(X)){
        ql <- pks[i] * dpois(X[,l], Lambdas[l])
        qk <- qk * ql
      }
      tk <- qk * pks[i] / d.k.s
      tks <- c(tks, tk)
    }

    print(sum(tks))

    eps = 10000000000000

    ##E-step
    #X$q_x2 <- (1-pk)*dpois(X,Lambda1)
    #X$P1_x <- pk*dpois(X,Lambda0)/(X$q_x1+X$q_x2) ##P=1|X
    #X$P2_x <- (1-pk)*dpois(X,Lambda1)/(X$q_x1+X$q_x2) ##P=2|X
    #Q <- sum(log(X$q_x1)*X$P1_x)+sum(log(X$q_x2)*X$P2_x)

    ##M-step/update parameters
    #Lambda0_k <- sum(X*X$P1_x)/sum(X$P1_x)
    #Lambda1_k <- sum(X*X$P2_x)/sum(X$P2_x)
    #pk_k <- sum(X$P1_x)/N

    ##compare Q
    #X$q_x1_k <- pk_k*dpois(X,Lambda0_k)
    #X$q_x2_k <- (1-pk_k)*dpois(X,Lambda1_k)
    #X$P1_x_k <- X$q_x1/(X$q_x1+X$q_x2)
    #X$P2_x_k <- X$q_x2/(X$q_x1+X$q_x2)
    #Q_k <- sum(log(X$q_x1_k)*X$P1_x_k)+sum(log(X$q_x2_k)*X$P2_x_k)

    ##stop criterion
    #error <- Q_k-Q
    #iter <- iter+1
    #Lambda0 <- Lambda0_k
    #Lambda1 <- Lambda1_k
    #pk <- pk_k
  }

  theta <-c (Lambda0,Lambda1,pk)

  ##cluster assingment based on posterior probability, normalizing term is canceled out
  X$cluster[pk*dpois(X,Lambda0)>=(1-pk)*dpois(X,Lambda1)] <- 1
  X$cluster[pk*dpois(X,Lambda0)<(1-pk)*dpois(X,Lambda1)] <- 2

  return(list(X,theta))
}

##generate simulated data
#set.seed(12345)

##mixture of two poissons
n <- 200  # Number of rows
lambda1 <- 3  # Poisson parameter for column 1
lambda2 <- 5  # Poisson parameter for column 2

# Generate random Poisson-distributed values
x1 <- rpois(n, lambda1)
x2 <- rpois(n, lambda2)

# Create the dataframe
X <- data.frame(x1 = x1, x2 = x2)

options(digits=22)
.Machine$double.eps<-2.220446e-20

####EXAMPLE: EM-based Clustering
a <- em_Poisson(X)

