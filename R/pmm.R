em_poisson <- function(X, clust=c(1,2,3), eps=0.002){

  x <- data.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  max_bic = 0
  df <- list(Lambdas = c(), pk = c(), clusters = c(), tk = c())

  for(K in clust){
    error <- Inf

    # initialization
    Lambdas <- matrix(0, nrow=K, ncol=p)
    for(k in 1:K){
      Lambdas[k,] <- runif(p, 1, 10)
    }
    pk <- rep(1/K, K)
    oldQ <- 0

    # EM algorithm
    while(TRUE){

      # E-STEP
      ## Total density
      ds <- c()
      for(k in 1:K){
        d_row <- c()
        d <- dpois(x, Lambdas[k,])
        for(i in 1:n){
          d_row <- c(d_row,prod(d[i,]))
        }
        ds <- c(ds,sum(d_row))
      }
      ds <- sum(ds)

      ## Calculate Tk
      pkds <- c()
      tk <- c()
      for(k in 1:K){
        d_res <- c()
        d <- dpois(x, Lambdas[k,])
        for(i in 1:n){
          d_res <- c(d_res, prod(d[i,]))
        }
        pkd <- pk[k] * d_res
        pkds <- cbind(pkds, pkd)
        tk <- cbind(tk, pkd/ds)
      }

      # M-step
      for(k in 1:K){
        for(j in 1:p){
          Lambdas[k,j] <- sum(x[,j]*tk[,k])/sum(tk[,k])
        }
        pk[k] <- sum(tk[,k])/n
      }

      # Calculate Likehood
      qs <- c()
      for(k in 1:K){
        q <- sum(log(pkds)*tk)
        qs <- c(qs,q)
      }
      qs <- sum(qs)

      # check convergence
      error <- abs(qs - oldQ)
      cat("error:", error, "\n")
      if(error < eps) break
      oldQ <- qs
    }

    # Determine label for each observation
    labels <- c()
    for(i in 1:n){
      max <- 0
      index <- 0
      for(k in 1:K){
        if(tk[i,k] > max){
          max <- tk[i,k]
          index <- k
        }
      }
      labels <- c(labels,index)
    }

    # BIC
    bic <- -2 * qs + (K*p + K - 1) * log(n)

    # Keep result for best BIC between all k
    if(bic > max_bic){
      max_bic = bic
      df$Lambdas <- Lambdas
      df$pk <- pk
      df$clusters <- labels
      df$tk <- tk
    }

  }

  return(df)
}

set.seed(14052000)

n <- 200  # Number of rows
lambda1 <- 3  # Poisson parameter for column 1
lambda2 <- 5  # Poisson parameter for column 2
lambda3 <- 7  # Poisson parameter for column 3

# Generate random Poisson-distributed values
x1 <- rpois(n, lambda1)
x2 <- rpois(n, lambda2)
x3 <- rpois(n, lambda3)

# Create the data frame
X <- data.frame(x1, x2,x3)

#### EXAMPLE: EM-based Clustering
res <- em_poisson(X)

# print(res)
