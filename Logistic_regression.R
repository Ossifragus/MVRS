## Logistic regression

L <- function(par, X, y) {
  p <- 1 - 1 / (1 + exp(X %*% par))
  return(-1 * sum(y * log(p) + (1 - y) * log(1 - p)))
}

Lw <- function(par, X, y, w) {
  p <- 1 - 1 / (1 + exp(X %*% par))
  return(-1 * sum(y * log(p) * w + (1 - y) * log(1 - p) * w))
}

sub_est_sim <- function(data, method, n0 = 200, n = 1000, sim_num = 2, d = 5, par0, m = 10) {
  N <- nrow(data)
  ## uniform subsampling
  if (method == "uniform") {
    par_uni <- matrix(0, nrow = sim_num, ncol = d)
    start_time <- Sys.time()
    for (i in 1:sim_num) {
      index <- sample(1:N, n, replace = TRUE)
      data_uni <- data[index, ]
      result <- optim(par = par00, L, hessian = FALSE, method = "Nelder-Mead", X = as.matrix(cbind(rep(1, length(index)), data_uni[, -1])), y = data_uni[, 1])
      par_uni[i, ] <- result$par
    }
    end_time <- Sys.time()
    MSE <- rep(0, d)
    Bias <- rep(0, d)
    for (i in 1:d) {
      MSE[i] <- mean((par_uni[, i] - par0[i])^2)
      Bias[i] <- mean(par_uni[, i] - par0[i])
    }
  }
  
  ## non-uniform subsampling
  if (method == "osmac") {
    par_sub <- matrix(0, nrow = sim_num, ncol = d)
    start_time <- Sys.time()
    for (i in 1:sim_num) {
      # pilot
      index <- sample(1:N, n0, replace = TRUE)
      data_poi <- data[index, ]
      result <- optim(par = par00, L, hessian = FALSE, method = "Nelder-Mead", X = as.matrix(cbind(rep(1, length(index)), data_poi[, -1])), y = data_poi[, 1])
      par_poi <- result$par
      # prob
      X <- as.matrix(cbind(rep(1, N), data[, -1]))
      p <- 1 - 1 / (1 + exp(X %*% par_poi))
      w <- p * (1 - p)
      w <- solve(t(X[index, ]) %*% (X[index, ] * w[index]))
      sub_pro <- sqrt((data[, 1] - p)^2 * rowSums((X %*% w)^2))
      sub_pro <- sub_pro / sum(sub_pro)
      # subsampling
      index <- sample(1:N, n, replace = TRUE, prob = sub_pro)
      data_sub <- data[index, ]
      result <- optim(par = par00, Lw, hessian = FALSE, method = "Nelder-Mead", X = as.matrix(cbind(rep(1, length(index)), data_sub[, -1])), y = data_sub[, 1], w = sub_pro[1] / sub_pro[index])
      par_sub[i, ] <- result$par
    }
    end_time <- Sys.time()
    MSE <- rep(0, d)
    Bias <- rep(0, d)
    for (i in 1:d) {
      MSE[i] <- mean((par_sub[, i] - par0[i])^2)
      Bias[i] <- mean(par_sub[, i] - par0[i])
    }
  }
  
  
  if (method == "stratified-uni") {
    par_str <- matrix(0, nrow = sim_num, ncol = d)
    start_time <- Sys.time()
    for (i in 1:sim_num) {
      # pilot
      index <- sample(1:N, n0, replace = TRUE)
      data_poi <- data[index, ]
      result <- optim(par = par00, L, hessian = FALSE, method = "Nelder-Mead", X = as.matrix(cbind(rep(1, length(index)), data_poi[, -1])), y = data_poi[, 1])
      par_poi <- result$par
      # stra
      X_poi <- as.matrix(cbind(rep(1, n0), data_poi[, -1]))
      y_poi <- data_poi[, 1]
      p1 <- 1 - 1 / (1 + exp(X_poi %*% par_poi))
      w1 <- p1 * (1 - p1)
      w1 <- solve(t(X_poi) %*% (X_poi * as.numeric(w1)) / n0)
      V <- cov(X_poi %*% w1 * as.numeric(y_poi - p1))
      S <- svd(V, d, d)
      A <- S$u
      sub_pro <- rep(1 / N, N)
      X <- as.matrix(cbind(rep(1, N), data[, -1]))
      p <- 1 - 1 / (1 + exp(X %*% par_poi))
      Y <- as.numeric(data[, 1] - p) * (X %*% w1 %*% A)
      order <- order(Y[, 1])
      
      # subsampling
      order_pro <- sub_pro[order]
      qua <- floor(N * seq(from = 0, to = 1, by = 1 / m))
      index <- c()
      prob1 <- c()
      for (j in 1:m) {
        nj <- round(n * sum(order_pro[(qua[j] + 1):qua[j + 1]]))
        indexj <- sample((qua[j] + 1):qua[j + 1], nj, replace = TRUE, prob = order_pro[(qua[j] + 1):qua[j + 1]])
        index <- c(index, indexj)
        prob1 <- c(prob1, nj * order_pro[indexj] / (sum(order_pro[(qua[j] + 1):qua[j + 1]])))
      }
      data_str <- data[order[index], ]
      result <- optim(par = par00, Lw, hessian = FALSE, method = "Nelder-Mead", X = as.matrix(cbind(rep(1, length(index)), data_str[, -1])), y = data_str[, 1], w = prob1[1] / prob1)
      par_str[i, ] <- result$par
    }
    end_time <- Sys.time()
    MSE <- rep(0, d)
    Bias <- rep(0, d)
    for (i in 1:d) {
      MSE[i] <- mean((par_str[, i] - par0[i])^2)
      Bias[i] <- mean(par_str[, i] - par0[i])
    }
  }
  
  if (method == "stratified-osmac") {
    par_str <- matrix(0, nrow = sim_num, ncol = d)
    start_time <- Sys.time()
    for (i in 1:sim_num) {
      # pilot
      index <- sample(1:N, n0, replace = TRUE)
      data_poi <- data[index, ]
      result <- optim(par = par00, L, hessian = FALSE, method = "Nelder-Mead", X = as.matrix(cbind(rep(1, length(index)), data_poi[, -1])), y = data_poi[, 1])
      par_poi <- result$par
      # stra
      X_poi <- as.matrix(cbind(rep(1, n0), data_poi[, -1]))
      y_poi <- data_poi[, 1]
      p1 <- 1 - 1 / (1 + exp(X_poi %*% par_poi))
      w1 <- p1 * (1 - p1)
      w1 <- solve(t(X_poi) %*% (X_poi * as.numeric(w1)) / n0)
      
      X <- as.matrix(cbind(rep(1, N), data[, -1]))
      p <- 1 - 1 / (1 + exp(X %*% par_poi))
      sub_pro <- sqrt((data[, 1] - p)^2 * rowSums((X %*% w1)^2))
      sub_pro <- sub_pro / sum(sub_pro)
      
      V <- cov(X_poi %*% w1 * as.numeric(y_poi - p1))
      S <- svd(V, d, d)
      A <- S$u
      Y <- as.numeric(data[, 1] - p) * (X %*% w1 %*% A)
      order <- order(Y[, 1])
      
      # subsampling
      order_pro <- sub_pro[order]
      qua <- floor(N * seq(from = 0, to = 1, by = 1 / m))
      index <- c()
      prob1 <- c()
      for (j in 1:m) {
        nj <- round(n * sum(order_pro[(qua[j] + 1):qua[j + 1]]))
        indexj <- sample((qua[j] + 1):qua[j + 1], nj, replace = TRUE, prob = order_pro[(qua[j] + 1):qua[j + 1]])
        index <- c(index, indexj)
        prob1 <- c(prob1, nj * order_pro[indexj] / (sum(order_pro[(qua[j] + 1):qua[j + 1]])))
      }
      data_str <- data[order[index], ]
      result <- optim(par = par00, Lw, hessian = FALSE, method = "Nelder-Mead", X = as.matrix(cbind(rep(1, length(index)), data_str[, -1])), y = data_str[, 1], w = prob1[1] / prob1)
      par_str[i, ] <- result$par
    }
    end_time <- Sys.time()
    MSE <- rep(0, d)
    Bias <- rep(0, d)
    for (i in 1:d) {
      MSE[i] <- mean((par_str[, i] - par0[i])^2)
      Bias[i] <- mean(par_str[, i] - par0[i])
    }
  }
  time <- (end_time - start_time) / sim_num
  result <- list(MSE, sum(MSE), time)
  return(result)
}

for (cas in 1:4) {
  N = 1000000
  library(MASS)
  Sigma <- matrix(c(1, 0.5, 0.5, 0.5,
                    0.5, 1, 0.5, 0.5,
                    0.5, 0.5, 1, 0.5,
                    0.5, 0.5, 0.5, 1), 4, 4)
  set.seed(123)
  
  if (cas == 1) x <- mvrnorm(N, rep(0, 4), Sigma)
  if (cas == 2) x <- mvrnorm(N, rep(1.5, 4), Sigma)
  if (cas == 3) x <- 0.5 * mvrnorm(N, rep(1, 4), Sigma) + 0.5 * mvrnorm(N, rep(-1, 4), Sigma)
  if (cas == 4) x <- matrix(rexp(4 * N, rate = 2), nrow = N, ncol = 4, byrow = TRUE)
  
  par <- 0.5 * c(1, 1, 1, 1, 1)
  p <- exp(par[1] + x %*% par[-1]) / (1 + exp(par[1] + x %*% par[-1]))
  u <- runif(N, 0, 1)
  y <- as.numeric(u < p)
  data <- cbind(y, x)
  data <- as.data.frame(data)
  
  X <- as.matrix(cbind(rep(1, N), data[, -1]))
  start_time <- Sys.time()
  par00 <- c(0, 0, 0, 0, 0)
  result <- optim(par = par00, L, hessian = FALSE, method = "Nelder-Mead", X = X, y = y)
  end_time <- Sys.time()
  end_time - start_time
  par0 <- result$par
  
  result_table <- matrix(rep(0, 16), nrow = 4, ncol = 4)
  col <- c(200, 500, 800, 1000)
  row <- c("uniform", "osmac", "stratified-uni", "stratified-osmac")
  for (i in row) {
    for (j in col) {
      result <- sub_est_sim(data, i, 200, j, 1000, 5, par0)
      i1 <- which(row == i)
      j1 <- which(col == j)
      result_table[i1, j1] <- result[[2]]
    }
  }
  result_table
  write.csv(result_table, paste0("logistic_regression_case", cas, ".csv"))
  
}

