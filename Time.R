## Poisson regression
L <- function(par, X, y) {
  X_par <- X %*% par
  return(-1 * sum(y * X_par - exp(X_par)) / length(y))
}

Lw <- function(par, X, y, w) {
  X_par <- X %*% par  
  return(-1 * sum((y * X_par - exp(X_par)) * w) / length(y))
}


sub_est_sim <- function(data, method, n0 = 200, n = 1000, sim_num = 2, d = 5, par0, m = 100) {
  N <- nrow(data)
  X <- cbind(rep(1, N), data[, -1])
  data <- as.matrix(cbind(data[, 1], X))
  ## uniform subsampling
  if (method == "uniform") {
    par_uni <- matrix(0, nrow = sim_num, ncol = d)
    start_time <- Sys.time()
    for (i in 1:sim_num) {
      index <- sample(1:N, n, replace = TRUE)
      data_uni <- data[index, ]
      result <- optim(par = par00, L, hessian = FALSE, method = "BFGS", X = data_uni[, -1], y = data_uni[, 1])
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
    par_sub <- matrix(0, nrow = sim_num,  ncol = d)
    start_time <- Sys.time()
    for (i in 1:sim_num) {
      # pilot
      index <- sample(1:N, n0, replace = TRUE)
      data_p <- data[index, ]
      result <- optim(par = par00, L, hessian = FALSE, method = "BFGS", X = data_p[, -1], y = data_p[, 1])
      par_p <- result$par
      # prob
      w <- solve(t(data_p[, -1]) %*% (data_p[, -1] * as.numeric(exp(data_p[, -1] %*% par_p))))
      sub_pro <- sqrt((data[, 1] - exp(data[, -1] %*% par_p))^2 * rowSums((data[, -1] %*% w)^2))
      # subsampling
      index <- sample(1:N, n, replace = TRUE, prob = sub_pro)
      data_sub <- data[index, ]
      result <- optim(par = par00, Lw, hessian = FALSE, method = "BFGS", X = data_sub[, -1], y = data_sub[, 1], w = sub_pro[1] / sub_pro[index])
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
      data_p <- data[index, ]
      result <- optim(par = par00, L, hessian = FALSE, method = "BFGS", X = data_p[, -1], y = data_p[, 1])
      par_p <- result$par
      # stra
      b <- as.numeric(data_p[, -1] %*% par_p)
      w_p <- solve(t(data_p[, -1]) %*% (data_p[, -1] * b) / n0) 
      V <- cov(data_p[, -1] %*% w_p * (data_p[, 1] - exp(b)))
      A <- svd(V, d, d)$u
      Y <- as.numeric(data[, 1] - exp(data[, -1] %*% result$par)) * (data[, -1] %*% (w_p %*% A[1, ]))
      order <- order(Y)
      # subsampling
      qua <- floor(N * seq(from = 0, to = 1, by = 1 / m))
      index <- c()
      prob1 <- c()
      for (j in 1:m) {
        nj <- round(n / m)   
        indexj <- sample((qua[j] + 1):qua[j + 1], nj, replace = TRUE)
        index <- c(index, indexj)
      }
      data_str <- data[order[index], ]
      
      result <- optim(par = par00, L, hessian = FALSE, method = "BFGS", X = data_str[, -1], y = data_str[, 1])
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
      data_p <- data[index, ]
      result <- optim(par = par00, L, hessian = FALSE, method = "BFGS", X = data_p[, -1], y = data_p[, 1])
      par_p <- result$par
      # prob
      b <- as.numeric(data_p[, -1] %*% par_p)
      w_p <- solve(t(data_p[, -1]) %*% (data_p[, -1] * b) / n0) 
      inf0 <- as.numeric(data[, 1] - exp(data[, -1] %*% result$par)) * data[, -1]
      sub_pro <- sqrt(rowSums((inf0 %*% w_p)^2))
      V <- cov(data_p[, -1] %*% w_p * (data_p[, 1] - exp(b)))
      A <- svd(V, d, d)$u
      Y <- inf0 %*% (w_p %*% A[1, ])
      order <- order(Y)
      # subsampling
      order_pro <- sub_pro[order]
      
      qua <- floor(N * seq(from = 0, to = 1, by = 1 / m))
      total_sub_pro <- sum(sub_pro)
      cum_pro <- cumsum(order_pro)
      group_sums <- cum_pro[qua[-1]] - c(0, cum_pro[qua[-c(1, m + 1)]])
      njs <- round(n * group_sums / total_sub_pro)    
      index <- c()
      prob1 <- c()
      for (j in 1:m) {
        order_pro_j <- order_pro[(qua[j] + 1):qua[j + 1]]
        nj <- njs[j]
        indexj <- sample((qua[j] + 1):qua[j + 1], nj, replace = TRUE, prob = order_pro_j)
        index <- c(index, indexj)
        prob1 <- c(prob1, nj * order_pro[indexj] / (sum(order_pro_j)))
      }
      data_str <- data[order[index], ]
      result <- optim(par = par00, Lw, hessian = FALSE, method = "BFGS", X = data_str[, -1], y = data_str[, 1], w = prob1[1] / prob1)
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
  
  if (method == "o-stratified-osmac") {
    par_str <- matrix(0, nrow = sim_num, ncol = d)
    start_time <- Sys.time()
    for (i in 1:sim_num) {
      # pilot
      index <- sample(1:N, n0, replace = TRUE)
      data_p <- data[index, ]
      result <- optim(par = par00, L, hessian = FALSE, method = "BFGS", X = data_p[, -1], y = data_p[, 1])
      par_p <- result$par
      # prob
      b <- as.numeric(data_p[, -1] %*% par_p)
      w_p <- solve(t(data_p[, -1]) %*% (data_p[, -1] * b) / n0) 
      inf0 <- as.numeric(data[, 1] - exp(data[, -1] %*% result$par)) * data[, -1]
      sub_pro <- sqrt(rowSums((inf0 %*% w_p)^2))
      sub_pro <- sub_pro / sum(sub_pro)
      V <- cov(data_p[, -1] %*% w_p * (data_p[, 1] - exp(b)))
      A <- svd(V, d, d)$u
      Y <- inf0 %*% (w_p %*% A[1, ])
      order <- order(Y)
      # subsampling
      order_pro <- sub_pro[order]
      pro_sum <- cumsum(order_pro)
      
      thresholds <- (1:(n - 1)) / n
      indices <- findInterval(thresholds, pro_sum) + 1
      qua <- c(0, indices, N)
      
      index <- c()
      prob1 <- c()
      for (j in 1:(n - 1)) {
        nj <- 1
        indexj <- sample((qua[j] + 1):qua[j + 1], nj, replace = TRUE, prob = order_pro[(qua[j] + 1):qua[j + 1]])
        index <- c(index, indexj)
        prob1 <- c(prob1, nj * order_pro[indexj] / (sum(order_pro[(qua[j] + 1):qua[j + 1]])))
      }
      data_str <- data[order[index], ]
      result <- optim(par = par00, Lw, hessian = FALSE, method = "BFGS", X = data_str[, -1], y = data_str[, 1], w = prob1[1] / prob1)
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
  
  if (method == "o-stratified-uni") {
    par_str <- matrix(0, nrow = sim_num, ncol = d)
    start_time <- Sys.time()
    for (i in 1:sim_num) {
      # pilot
      index <- sample(1:N, n0, replace = TRUE)
      data_p <- data[index, ]
      result <- optim(par = par00, L, hessian = FALSE, method = "BFGS", X = data_p[, -1], y = data_p[, 1])
      par_p <- result$par
      # stra
      b <- as.numeric(data_p[, -1] %*% par_p)
      w_p <- solve(t(data_p[, -1]) %*% (data_p[, -1] * b) / n0) 
      V <- cov(data_p[, -1] %*% w_p * (data_p[, 1] - exp(b)))
      A <- svd(V, d, d)$u
      Y <- as.numeric(data[, 1] - exp(data[, -1] %*% result$par)) * (data[, -1] %*% (w_p %*% A[1, ]))
      order <- order(Y)
      sub_pro <- rep(1 / N, N)
      # subsampling
      order_pro <- sub_pro[order]
      pro_sum <- cumsum(order_pro)
      # subsampling
      order_pro <- sub_pro[order]
      pro_sum <- cumsum(order_pro)
      qua <- floor(N * seq(from = 0, to = 1, by = 1 / n))
      index <- c()
      prob1 <- c()
      for (j in 1:(n - 1)) {
        nj <- 1
        indexj <- sample((qua[j] + 1):qua[j + 1], nj, replace = TRUE, prob = order_pro[(qua[j] + 1):qua[j + 1]])
        index <- c(index, indexj)
        prob1 <- c(prob1, nj * order_pro[indexj] / (sum(order_pro[(qua[j] + 1):qua[j + 1]])))
      }
      data_str <- data[order[index], ]
      result <- optim(par = par00, L, hessian = FALSE, method = "BFGS", X = data_str[, -1], y = data_str[, 1])
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


N <- 1000000
library(MASS)
set.seed(123)
x <- matrix(runif(4 * N, 0, 1), nrow = N, ncol = 4, byrow = TRUE)
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

result_table_time <- matrix(rep(0, 28), nrow = 7, ncol = 4)

result <- optim(par = par00, L, hessian = FALSE, method = "BFGS", X = X, y = y)
end_time <- Sys.time()
result_table_time[1, 1] <- end_time - start_time
par0 <- result$par

result <- sub_est_sim(data, "uniform", 200, 1000, 100, 5, par0)
result_table_time[2, 1] <- result[[3]]
result <- sub_est_sim(data, "osmac", 200, 1000, 100, 5, par0)
result_table_time[3, 1] <- result[[3]]

col <- c(5, 10, 50, 100)
row <- c("stratified-uni", "stratified-osmac")
for (i in row) {
  for (m in col) {
    result <- sub_est_sim(data, i, 200, 1000, 100, 5, par0, m)
    i1 <- which(row == i) + 3
    j1 <- which(col == m)
    result_table_time[i1, j1] <- result[[3]]
  }
}

result <- sub_est_sim(data, "o-stratified-osmac", 200, 1000, 100, 5, par0)
result_table_time[6, 1] <- result[[3]]

result <- sub_est_sim(data, "o-stratified-uni", 200, 1000, 100, 5, par0)
result_table_time[7, 1] <- result[[3]]

write.csv(result_table_time, "time.csv")