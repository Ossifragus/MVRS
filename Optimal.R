## Poisson regression
L <- function(par, X, y) {
  return(-1 * sum(y * X %*% par - exp(X %*% par)) / length(y))
}

Lw <- function(par, X, y, w) {
  return(-1 * sum((y * X %*% par - exp(X %*% par)) * w) / length(y))
}

sub_est_sim <- function(data, method, n0 = 200, n = 1000, sim_num = 2, d = 5, par0, m = 100) {
  N <- nrow(data)
  X <- cbind(rep(1, N), data[, -1])
  data <- as.matrix(cbind(data[, 1], X))
  
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
result <- optim(par = par00, L, hessian = FALSE, method = "BFGS", X = X, y = y)
end_time <- Sys.time()
end_time - start_time
par0 <- result$par

result_table <- matrix(rep(0, 8), nrow = 2, ncol = 4)
col <- c(200, 500, 800, 1000)
row <- c("o-stratified-uni", "o-stratified-osmac")
for (i in row) {
  for (j in col) {
    result <- sub_est_sim(data, i, 200, j, 1000, 5, par0)
    i1 <- which(row == i)
    j1 <- which(col == j)
    result_table[i1, j1] <- result[[2]]
  }
}
result_table
write.csv(result_table, "poisson_regression_case1_opt.csv")