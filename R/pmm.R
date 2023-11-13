em_poisson <- function(X, K=2, eps=0.02){

  error <- Inf
  x <- data.matrix(X)

  #initialization
  n <- nrow(X)
  p <- ncol(X)

  Lambdas <- matrix(0, nrow=K, ncol=p)
  for(k in 1:K){
    Lambdas[k,] <- runif(p, 1, 10)
  }
  pk <- rep(1/K, K)
  oldQ <- 0

  #EM algorithm
  while(TRUE){

    #E-step
    ds <- 1
    for(k in 1:K){
      ds <- ds + dpois(x, Lambdas[k,])
    }

    pkds <- c()
    tk <- c()
    for(k in 1:K){
      d <- dpois(x, Lambdas[k,])
      pkd <- pk[k] * d
      pkds <- cbind(pkds, pkd)
      tk <- cbind(tk, pkd/ds)
    }

    #M-step
    for(k in 1:K){
      for(j in 1:p){
        Lambdas[k,j] <- sum(x[,j]*tk[,k])/sum(tk[,k])
      }
      pk[k] <- sum(tk[,k])/n
    }

    Qs <- 1
    for(k in 1:K){
      min <- 1 + (k-1)*p
      max <- k*p
      Q <- sum(log(pkds[,min:max])*tk[,min:max])
      Qs <- Qs + Q
    }

    #check convergence
    error <- abs(Qs - oldQ)
    if(error < eps) break
    oldQ <- Qs
  }

  # assign clusters based on posterior probabilities
  labels <- matrix(0, nrow=n, ncol=p)
  for(i in 1:n){
    for(j in 1:p){
      max = 0
      for(k in 1:K){
        kk <- j + (k-1)*p
        if(tk[i,kk] > max){
          max <- tk[i,kk]
          labels[i,j] <- k
        }
      }
    }
  }

  return(list(Lambdas=Lambdas, pk=pk, clusters=labels, tk=tk))
}
