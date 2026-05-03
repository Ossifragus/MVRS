getMLE <- function(x, y, w) {
  d <- ncol(x)
  beta <- rep(0, d)
  loop  <- 1
  Loop  <- 1000
  msg <- "NA"
  while (loop <= Loop) {
    pr <- c(1 - 1 / (1 + exp(x %*% beta)))
    H <- t(x) %*% (pr * (1 - pr) * w * x)
    S <- colSums((y - pr) * w * x)
    tryCatch(
      {shs <- NA
      shs <- solve(H, S) },
      error=function(e){
        cat("\n ERROR :", loop, conditionMessage(e), "\n")})
    if (is.na(shs[1])) {
      msg <- "Not converge"
      beta <- loop <- NA
      break
    }
    beta.new <- beta + shs
    tlr  <- sum((beta.new - beta)^2)
    beta  <- beta.new
    if(tlr < 0.000001) {
      msg <- "Successful convergence"
      break
    }
    if (loop == Loop)
      warning("Maximum iteration reached")
    loop  <- loop + 1
  }
  list(par=beta, message=msg, iter=loop)
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

sub_est_sim <- function(data, method, n0 = 200, n = 1000, sim_num = 100, d = 15, par0, m = 10) {
  N <- nrow(data)
  ## uniform subsampling
  if (method == "uniform") {
    par_uni <- matrix(0, nrow = sim_num, ncol = d)
    mse_uni <- matrix(0, nrow = sim_num, ncol = d)
    start_time <- Sys.time()
    for (i in 1:sim_num) {
      index <- sample(1:N, n, replace = TRUE)
      data_uni <- data[index, ]
      X <- as.matrix(cbind(rep(1, n), data_uni[, -1]))
      y <- data_uni[, 1]
      result<-getMLE(X, y, rep(1/n,n))
      par_uni[i, ] <- result$par
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
      data_p <- data[index, ]
      result<-getMLE(as.matrix(cbind(rep(1,length(index)),data_p[,-1])),data_p[,1],rep(1/length(data_p[,1]),length(data_p[,1])))
      par_p <- result$par
      # prob
      X <- as.matrix(cbind(rep(1, N), data[, -1]))
      p <- 1 - 1 / (1 + exp(X %*% par_p))
      w <- p * (1 - p)
      w <- solve(t(X[index, ]) %*% (X[index, ] * w[index]))
      sub_pro <- sqrt((data[, 1] - p)^2 * rowSums((X %*% w)^2))
      sub_pro <- sub_pro / sum(sub_pro)
      # subsampling
      index <- sample(1:N, n, replace = TRUE, prob = sub_pro)
      X_sub <- as.matrix(cbind(rep(1, length(index)), data[index, -1]))
      y_sub <- data[index, 1]
      result<-getMLE(X_sub,y_sub,1/sub_pro[index])
      par_sub[i, ] <- result$par
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
      data_p <- data[index, ]
      result<-getMLE(as.matrix(cbind(rep(1,length(index)),data_p[,-1])),data_p[,1],rep(1/length(data_p[,1]),length(data_p[,1])))
      par_p <- result$par
      # stra
      X_p <- as.matrix(cbind(rep(1, n0), data_p[, -1]))
      y_p <- data_p[, 1]
      p1 <- 1 - 1 / (1 + exp(X_p %*% par_p))
      w1 <- p1 * (1 - p1)
      w1 <- solve(t(X_p) %*% (X_p * as.numeric(w1)) / n0)
      V <- cov(X_p %*% w1* as.numeric(y_p - p1))
      S <- svd(V, d, d)
      A <- S$u

      sub_pro <- rep(1 / N, N)
      X <- as.matrix(cbind(rep(1, N), data[, -1]))
      p <- 1 - 1 / (1 + exp(X %*% par_p))
      Y <- as.numeric(data[, 1] - p) * (X %*% (w1 %*% A[, 1]))
      str_result <- stratify(Y, m)
      
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
      result<-getMLE(as.matrix(cbind(rep(1,length(index)),data_str[,-1])),data_str[,1],prob1[1]/prob1)
      par_str[i, ] <- result$par
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
      data_p <- data[index, ]
      result<-getMLE(as.matrix(cbind(rep(1,length(index)),data_p[,-1])),data_p[,1],rep(1/length(data_p[,1]),length(data_p[,1])))
      par_p <- result$par
      # stra
      X_p <- as.matrix(cbind(rep(1, n0), data_p[, -1]))
      y_p <- data_p[, 1]
      p1 <- 1 - 1 / (1 + exp(X_p %*% par_p))
      w1 <- p1 * (1 - p1)
      w1 <- solve(t(X_p) %*% (X_p * as.numeric(w1)) / n0)
      
      X <- as.matrix(cbind(rep(1, N), data[, -1]))
      p <- 1 - 1 / (1 + exp(X %*% par_p))
      sub_pro <- sqrt((data[, 1] - p)^2 * rowSums((X %*% w1)^2))
      sub_pro <- sub_pro / sum(sub_pro)
      
      V <- cov(X_p %*% w1 * as.numeric(y_p - p1))
      S <- svd(V, d, d)
      A <- S$u
      Y <- as.numeric(data[, 1] - p) * (X %*% (w1 %*% A[, 1]))
      str_result <- stratify(Y, m)
      
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
      result<-getMLE(as.matrix(cbind(rep(1,length(index)),data_str[,-1])),data_str[,1],prob1[1]/prob1)
      par_str[i, ] <- result$par
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

library(MASS)
library(mvtnorm)
for (cas in 1:4) {
  N = 500000
  par <- 0.1 * rep(1, 15)
  d <- length(par)

  set.seed(123)
  if (cas == 1) {
    Sigma <- matrix(0.5, nrow = d-1, ncol = d-1)+ 0.5 * diag(d-1)
    x <- mvrnorm(N, rep(0, d-1), Sigma)
  }
  if (cas == 2) {
    Sigma <- matrix(0.5, nrow = d-1, ncol = d-1)+ 0.5 * diag(d-1)
    x <- mvrnorm(N, rep(1.5, d-1), Sigma)
  }
  if (cas == 3) {
    sds <- 1 / (1:(d-1))
    Corr_mat <- matrix(0.5, nrow = (d-1), ncol = (d-1))
    diag(Corr_mat) <- 1
    Sigma <- diag(sds) %*% Corr_mat %*% diag(sds)
    x <- mvrnorm(N, rep(0, d-1), Sigma)
  }
  if (cas == 4) x <- matrix(rexp((d-1) * N, rate = 2), nrow = N, ncol = d-1, byrow = TRUE)
  p <- exp(par[1] + x %*% par[-1]) / (1 + exp(par[1] + x %*% par[-1]))
  u <- runif(N, 0, 1)
  y <- as.numeric(u < p)
  sum(y)
  data <- cbind(y, x)
  data <- as.data.frame(data)
  
  X <- as.matrix(cbind(rep(1, N), data[, -1]))
  start_time <- Sys.time()
  par00 <- rep(0, d)
  result <- getMLE(X,y,rep(1/N,N))
  end_time <- Sys.time()
  end_time - start_time
  par0 <- result$par
  
  result_table_mse <- matrix(rep(0, 10), nrow = 2, ncol = 5)
  col <- c(1, 5, 10, 50, 100)
  row <- c("stratified-uni", "stratified-osmac")
  for (i in row) {
    for (j in col) {
      result <- sub_est_sim(data, i, 500, 1000, 1000, d, par0, m=j)
      i1 <- which(row == i)
      j1 <- which(col == j)
      result_table_mse[i1, j1] <- result[[2]]
    }
    print(i)
  }
  result_table_mse
  write.csv(result_table_mse, paste0("logistic_regression_case", cas, "(str_num_n=1000).csv"))
}

