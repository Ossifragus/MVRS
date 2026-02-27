## Poisson regression
L <- function(par, X, y) {
  return(-1 * sum(y * X %*% par - exp(X %*% par)) / length(y))
}

Lw <- function(par, X, y, w) {
  return(-1 * sum((y * X %*% par - exp(X %*% par)) * w) / length(y))
}

library(data.table)
stratify <- function(value, k) {
  dt <- data.table(
    value = value,
    index0 = seq_along(value)  
  )
  breaks <- dt[, quantile(value, probs = seq(0, 1, length.out = k + 1))]
  dt[, layer := cut(value, breaks = breaks, include.lowest = TRUE, labels = FALSE)]
  layers <- split(dt$value, dt$layer)
  index <- split(dt$index0, dt$layer)
  return(list(
    layers = layers,
    index = index
  ))
}

sub_est_sim <- function(data, method, n0 = 200, n = 500, sim_num = 200, d = 5, par0, m = 10) {
  N <- nrow(data)
  ## uniform subsampling
  if (method == "uniform") {
    par_uni <- matrix(0, nrow = sim_num, ncol = d)
    mse_uni <- matrix(0, nrow = sim_num, ncol = d)
    start_time <- Sys.time()
    for (i in 1:sim_num) {
      # subsampling
      index <- sample(1:N, n, replace = TRUE)
      data_uni <- data[index, ]
      # estimation
      result <- optim(par = par00, L, hessian = FALSE, method = "BFGS", X = as.matrix(cbind(rep(1, length(index)), data_uni[, -1])), y = data_uni[, 1])
      par_uni[i, ] <- result$par
      # MSE estimation
      X <- as.matrix(cbind(rep(1, n), data_uni[, -1]))
      w <- solve(t(X) %*% (X * as.numeric(exp(X %*% result$par))) / n)
      inf <- as.numeric(data_uni[, 1] - exp(X %*% result$par)) * (X %*% w) / n
      mse_uni[i, ] <- diag(t(inf) %*% inf) * n / (n - d)
    }
    end_time <- Sys.time()
    MSE <- rep(0, d)
    MSE_est <- rep(0, d)
    Bias <- rep(0, d)
    for (i in 1:d) {
      MSE[i] <- mean((par_uni[, i] - par0[i])^2)
      Bias[i] <- mean(par_uni[, i] - par0[i])
      MSE_est[i] <- mean(mse_uni[, i])
    }
  }
  
  ## non-uniform subsampling
  if (method == "osmac") {
    par_sub <- matrix(0, nrow = sim_num, ncol = d)
    mse_sub <- matrix(0, nrow = sim_num, ncol = d)
    start_time <- Sys.time()
    for (i in 1:sim_num) {
      # pilot
      index <- sample(1:N, n0, replace = TRUE)
      data_poi <- data[index, ]
      result <- optim(par = par00, L, hessian = FALSE, method = "BFGS", X = as.matrix(cbind(rep(1, length(index)), data_poi[, -1])), y = data_poi[, 1])
      par_poi <- result$par
      # prob
      X <- as.matrix(cbind(rep(1, N), data[, -1]))
      w <- solve(t(X[index, ]) %*% (X[index, ] * as.numeric(exp(X[index, ] %*% par_poi))))
      sub_pro <- sqrt((data[, 1] - exp(X %*% par_poi))^2 * rowSums((X %*% w)^2))
      sub_pro <- sub_pro / sum(sub_pro)
      # subsampling
      index <- sample(1:N, n, replace = TRUE, prob = sub_pro)
      data_sub <- data[index, ]
      result <- optim(par = par00, Lw, hessian = FALSE, method = "BFGS", X = as.matrix(cbind(rep(1, length(index)), data_sub[, -1])), y = data_sub[, 1], w = sub_pro[1] / sub_pro[index])
      par_sub[i, ] <- result$par
      # MSE estimation
      X <- as.matrix(cbind(rep(1, n), data_sub[, -1]))
      w <- solve(t(X / sub_pro[index] / N) %*% (X * as.numeric(exp(X %*% par_poi))) / n)
      inf <- as.numeric(data_sub[, 1] - exp(X %*% result$par)) * (X %*% w) / sub_pro[index] / n / N
      mse_sub[i, ] <- diag(t(inf) %*% inf) * n / (n - d)
    }
    end_time <- Sys.time()
    MSE <- rep(0, d)
    MSE_est <- rep(0, d)
    Bias <- rep(0, d)
    for (i in 1:d) {
      MSE[i] <- mean((par_sub[, i] - par0[i])^2)
      Bias[i] <- mean(par_sub[, i] - par0[i])
      MSE_est[i] <- mean(mse_sub[, i])
    }
  }
  
  if (method == "stratified-uni") {
    par_str <- matrix(0, nrow = sim_num, ncol = d)
    mse_str <- matrix(0, nrow = sim_num, ncol = d)
    start_time <- Sys.time()
    for (i in 1:sim_num) {
      # pilot
      index <- sample(1:N, n0, replace = TRUE)
      data_poi <- data[index, ]
      result <- optim(par = par00, L, hessian = FALSE, method = "BFGS", X = as.matrix(cbind(rep(1, length(index)), data_poi[, -1])), y = data_poi[, 1])
      par_poi <- result$par
      # stra
      X_poi <- as.matrix(cbind(rep(1, n0), data_poi[, -1]))
      y_poi <- data_poi[, 1]
      w_poi <- solve(t(X_poi) %*% (X_poi * as.numeric(X_poi %*% par_poi)) / n0)
      V <- cov(X_poi %*% w_poi * as.numeric(y_poi - exp(X_poi %*% par_poi)))
      A <- svd(V, d, d)$u
      sub_pro <- rep(1 / N, N)
      X <- as.matrix(cbind(rep(1, N), data[, -1]))
      Y <- as.numeric(data[, 1] - exp(X %*% result$par)) * (X %*% w_poi %*% A)
      str_result <- stratify(Y[, 1], m)
      
      # subsampling
      index <- c()
      prob1 <- c()
      num <- c(0)
      for (j in 1:m) {
        indexj <- str_result$index[[j]]
        nj <- round(n * sum(sub_pro[indexj]))
        sub_indexj <- sample(indexj, nj, replace = TRUE, prob = sub_pro[indexj])
        index <- c(index, sub_indexj)
        prob1 <- c(prob1, nj * sub_pro[sub_indexj] / (sum(sub_pro[indexj])))
        num <- c(num, nj)
      }
      data_str <- data[index, ]
      result <- optim(par = par00, Lw, hessian = FALSE, method = "BFGS", X = as.matrix(cbind(rep(1, length(index)), data_str[, -1])), y = data_str[, 1], w = prob1[1] / prob1)
      par_str[i, ] <- result$par
      # MSE estimation
      X <- as.matrix(cbind(rep(1, nrow(data_str)), data_str[, -1]))
      w <- solve(t(X / prob1 / N) %*% (X * as.numeric(exp(X %*% result$par))))
      inf <- matrix(0, d, d)
      for (j in 1:m) {
        indexj <- (sum(num[1:j]) + 1):sum(num[1:j + 1])
        nj <- length(indexj)
        Xj <- as.matrix(cbind(rep(1, nj), data_str[indexj, -1]))
        score <- as.numeric(data_str[indexj, 1] - exp(Xj %*% result$par)) * Xj / prob1[indexj] / N * nj
        inf <- inf + cov(score) / sum(sub_pro[str_result$index[[j]]])
      }
      mse_str[i, ] <- diag(w %*% inf %*% w) / n * n / (n - d)
    }
    end_time <- Sys.time()
    MSE <- rep(0, d)
    MSE_est <- rep(0, d)
    Bias <- rep(0, d)
    for (i in 1:d) {
      MSE[i] <- mean((par_str[, i] - par0[i])^2)
      Bias[i] <- mean(par_str[, i] - par0[i])
      MSE_est[i] <- mean(mse_str[, i])
    }
  }
  
  if (method == "stratified-osmac") {
    par_str <- matrix(0, nrow = sim_num, ncol = d)
    mse_str <- matrix(0, nrow = sim_num, ncol = d)
    start_time <- Sys.time()
    for (i in 1:sim_num) {
      # pilot
      index <- sample(1:N, n0, replace = TRUE)
      data_poi <- data[index, ]
      result <- optim(par = par00, L, hessian = FALSE, method = "BFGS", X = as.matrix(cbind(rep(1, length(index)), data_poi[, -1])), y = data_poi[, 1])
      par_poi <- result$par
      # prob
      X <- as.matrix(cbind(rep(1, N), data[, -1]))
      w <- solve(t(X[index, ]) %*% (X[index, ] * as.numeric(exp(X[index, ] %*% par_poi))))
      sub_pro <- sqrt((data[, 1] - exp(X %*% par_poi))^2 * rowSums((X %*% w)^2))
      sub_pro <- sub_pro / sum(sub_pro)
      
      X_poi <- as.matrix(cbind(rep(1, n0), data_poi[, -1]))
      y_poi <- data_poi[, 1]
      w_poi <- solve(t(X_poi) %*% (X_poi * as.numeric(X_poi %*% par_poi)) / n0)
      V <- cov(X_poi %*% w_poi * as.numeric(y_poi - exp(X_poi %*% par_poi)))
      A <- svd(V, d, d)$u
      Y <- as.numeric(data[, 1] - exp(X %*% result$par)) * (X %*% w_poi %*% A)
      str_result <- stratify(Y[, 1], m)
      
      # subsampling
      index <- c()
      prob1 <- c()
      num <- c(0)
      for (j in 1:m) {
        indexj <- str_result$index[[j]]
        nj <- round(n * sum(sub_pro[indexj]))
        sub_indexj <- sample(indexj, nj, replace = TRUE, prob = sub_pro[indexj])
        index <- c(index, sub_indexj)
        prob1 <- c(prob1, nj * sub_pro[sub_indexj] / (sum(sub_pro[indexj])))
        num <- c(num, nj)
      }
      data_str <- data[index, ]
      result <- optim(par = par00, Lw, hessian = FALSE, method = "BFGS", X = as.matrix(cbind(rep(1, length(index)), data_str[, -1])), y = data_str[, 1], w = prob1[1] / prob1)
      par_str[i, ] <- result$par
      # MSE estimation
      X <- as.matrix(cbind(rep(1, nrow(data_str)), data_str[, -1]))
      w <- solve(t(X / prob1 / N) %*% (X * as.numeric(exp(X %*% result$par))))
      inf <- matrix(0, d, d)
      for (j in 1:m) {
        indexj <- (sum(num[1:j]) + 1):sum(num[1:j + 1])
        nj <- length(indexj)
        Xj <- as.matrix(cbind(rep(1, nj), data_str[indexj, -1]))
        score <- as.numeric(data_str[indexj, 1] - exp(Xj %*% result$par)) * Xj / prob1[indexj] / N * nj
        inf <- inf + cov(score) / sum(sub_pro[str_result$index[[j]]])
      }
      mse_str[i, ] <- diag(w %*% inf %*% w) / n * n / (n - d)
    }
    end_time <- Sys.time()
    MSE <- rep(0, d)
    MSE_est <- rep(0, d)
    Bias <- rep(0, d)
    for (i in 1:d) {
      MSE[i] <- mean((par_str[, i] - par0[i])^2)
      Bias[i] <- mean(par_str[, i] - par0[i])
      MSE_est[i] <- mean(mse_str[, i])
    }
  }
  
  time <- (end_time - start_time) / sim_num
  result <- list(MSE, sum(MSE), sum(MSE_est), time)
  return(result)
}

for (cas in 1:4) {
  N <- 1000000
  library(MASS)
  set.seed(123)
  if (cas == 1) { x <- matrix(runif(4 * N, 0, 1), nrow = N, ncol = 4, byrow = TRUE) }
  if (cas == 2) { 
    x1 <- runif(N, 0, 1)
    x <- matrix(c(x1, x1 + runif(N, 0, 0.1), runif(2 * N, 0, 1)), nrow = N, ncol = 4, byrow = TRUE) 
  }
  if (cas == 3) { 
    x1 <- runif(N, 0, 1)
    x <- matrix(c(x1, x1 + runif(N, 0, 1), runif(2 * N, 0, 1)), nrow = N, ncol = 4, byrow = TRUE)
    x1 <- runif(N, 0, 1) 
  }
  if (cas == 4) { 
    x1 <- runif(N, 0, 1)
    x <- matrix(c(x1, x1 + runif(N, 0, 1), runif(2 * N, -1, 1)), nrow = N, ncol = 4, byrow = TRUE) 
  }
  
  par <- 0.5 * c(1, 1, 1, 1, 1)
  y <- rep(0, N)
  for (i in 1:N) {
    y[i] <- rpois(1, exp(par[1] + x[i, ] %*% par[-1]))
  }
  data <- cbind(y, x)
  data <- as.data.frame(data)
  
  X <- as.matrix(cbind(rep(1, N), data[, -1]))
  start_time <- Sys.time()
  par00 <- c(0, 0, 0, 0, 0)
  result <- optim(par = par00, L, hessian = FALSE, method = "BFGS", X = X, y = y)
  end_time <- Sys.time()
  end_time - start_time
  par0 <- result$par
  
  result_table_mse <- matrix(rep(0, 16), nrow = 4, ncol = 4)
  result_table_mse_est <- matrix(rep(0, 16), nrow = 4, ncol = 4)
  col <- c(200, 500, 800, 1000)
  row <- c("uniform", "osmac", "stratified-uni", "stratified-osmac")
  for (i in row) {
    for (j in col) {
      result <- sub_est_sim(data, i, 200, j, 1000, 5, par0)
      i1 <- which(row == i)
      j1 <- which(col == j)
      result_table_mse[i1, j1] <- result[[2]]
      result_table_mse_est[i1, j1] <- result[[3]]
    }
  }
  result_table_mse
  result_table_mse_est
  
  write.csv(result_table_mse, paste0("emp_poisson_regression_case", cas, ".csv"))
  write.csv(result_table_mse_est, paste0("est_poisson_regression_case", cas, ".csv"))
}
