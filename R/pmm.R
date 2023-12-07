em_poisson <- function(X, clust = c(2, 3, 5), eps = 0.002) {
  x <- data.matrix(X)
  n <- nrow(X)
  p <- ncol(X)

  max_bic <- -Inf
  df <- list(Lambdas = c(), pk = c(), clusters = c(), tk = c(), bic = c())
  bics <- c()
  qss <- c()

  for (K in clust) {
    error <- Inf

    # initialization
    Lambdas <- matrix(0, nrow = K, ncol = p)
    for (k in 1:K) {
      Lambdas[k, ] <- runif(p, 1, 10)
    }
    pk <- rep(1 / K, K)
    oldQ <- 0

    # EM algorithm
    while (TRUE) {
      # E-STEP
      ## Total density
      pkds <- c()
      for (k in 1:K) {
        d_row <- c()
        d_total <- c()
        for (c in 1:p) {
          dp <- dpois(x[, c], Lambdas[k, c])
          d_row <- cbind(d_row, dp)
        }
        for (i in 1:n) {
          d_total <- c(d_total, prod(d_row[i, ]))
        }
        pkd <- pk[k] * d_total
        pkds <- cbind(pkds, pkd)
      }
      tk <- pkds / rowSums(pkds)
      colnames(tk) <- c(1:K)

      # Calculate log Likehood
      qs <- 0
      for (i in 1:n) {
        for (k in 1:K) {
          qs <- qs + sum(tk[i, k] * log(dpois(x[i, ], Lambdas[k, ])))
        }
      }
      qss <- c(qss, qs)

      # Alternative way to calculate log Likehood
      # qs <- 1
      # for (i in 1:n){
      #  q_k <- 0
      #  for(k in 1:K){
      #    q_k <- q_k + (pk[k] * d_total[i,k])
      #  }
      #  qs <- qs * q_k
      # }

      # M-step
      for (k in 1:K) {
        for (j in 1:p) {
          Lambdas[k, j] <- sum(x[, j] * tk[, k]) / sum(tk[, k])
        }
        pk[k] <- sum(tk[, k]) / n
      }

      # check convergence
      error <- abs(qs - oldQ)
      cat("error:", error, "\n")
      if (error < eps) break
      oldQ <- qs
    }

    # Determine label for each observation
    labels <- c()
    for (i in 1:n) {
      max <- 0
      index <- 0
      for (k in 1:K) {
        if (tk[i, k] > max) {
          max <- tk[i, k]
          index <- k
        }
      }
      labels <- c(labels, index)
    }

    # BIC
    bic <- 2 * qs - (K * p + K - 1) * log(n) # (K*p + K - 1) -> v
    bics <- c(bics, bic)
    cat("BIC : ", bic, " for k = ", K, "\n")

    # Keep results for best BIC between all k
    if (bic > max_bic) {
      max_bic <- bic
      df$Lambdas <- Lambdas
      df$pk <- pk
      df$clusters <- labels
      df$tk <- tk
      df$bic <- c(as.integer(K), bic)
    }
  }
  # Plot BIC
  plot(clust, bics, type="b", xlab="Number of clusters", ylab="BIC", main="BIC for Poisson Mixture Model")

  # Plot log Likehood
  plot(qss, type="b", xlab="Number of iterations", ylab="Log Likehood", main="Log Likehood for Poisson Mixture Model")

  return(df)
}
