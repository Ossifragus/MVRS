## Logistic regression
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



sub_est_sim <- function(data, method, n0 = 1000, n = 1000, sim_num = 20, d = 15, par0, m = 10) {
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
      result<-getMLE(data_uni[, -1],data_uni[,1],rep(1/length(data_uni[,1]),length(data_uni[,1])))
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
      result<-getMLE(data_p[, -1],data_p[, 1],rep(1/length(data_p[,1]),length(data_p[,1])))
      par_p <- result$par
      # prob
      p <- 1 - 1 / (1 + exp(data[, -1] %*% par_p))
      p_p <- p[index]
      w_p <- p_p * (1 - p_p)
      w_p <- solve(t(data_p[, -1]) %*% (data_p[, -1] * w_p))
      sub_pro <- sqrt((data[, 1] - p)^2 * rowSums((data[, -1] %*% w_p)^2))
      #sub_pro <- sub_pro / sum(sub_pro)
      # subsampling
      index <- sample(1:N, n, replace = TRUE, prob = sub_pro)
      data_sub <- data[index, ]
      result<-getMLE(data_sub[, -1],data_sub[,1],1/sub_pro[index])
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
      result<-getMLE(data_p[, -1],data_p[, 1],rep(1/length(data_p[,1]),length(data_p[,1])))
      par_p <- result$par
      # stra
      p <- 1 - 1 / (1 + exp(data[, -1] %*% par_p))
      p_p <- p[index]
      w_p <- p_p * (1 - p_p)
      w_p <- solve(t(data_p[, -1]) %*% (data_p[, -1] * w_p))
      V <- cov(data_p[, -1] %*% w_p * as.numeric(data_p[, 1] - p_p))
      A <- svd(V, d, d)$u
      Y <- as.numeric(data[, 1] - p) * (data[, -1] %*% (w_p %*% A[, 1]))
      order <- order(Y)
      # subsampling
      qua <- floor(N * seq(from = 0, to = 1, by = 1 / m))
      index <- c()
      prob1 <- c()
      for (j in 1:m) {
        nj <- round(n / m)   
        indexj <- sample((qua[j] + 1):qua[j + 1], nj, replace = TRUE)
        index <- c(index, indexj)
        prob1 <- c(prob1, nj / length((qua[j] + 1):qua[j + 1]) * N)
      }
      data_str <- data[order[index], ]
      result<-getMLE(data_str[, -1],data_str[,1],1/prob1)
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
      result<-getMLE(data_p[, -1],data_p[, 1],rep(1/length(data_p[,1]),length(data_p[,1])))
      par_p <- result$par
      # prob
      p <- 1 - 1 / (1 + exp(data[, -1] %*% par_p))
      p_p <- p[index]
      w_p <- p_p * (1 - p_p)
      w_p <- solve(t(data_p[, -1]) %*% (data_p[, -1] * w_p))
      sub_pro <- sqrt((data[, 1] - p)^2 * rowSums((data[, -1] %*% w_p)^2))
      V <- cov(data_p[, -1] %*% w_p * as.numeric(data_p[, 1] - p_p))
      A <- svd(V, d, d)$u
      Y <- as.numeric(data[, 1] - p) * (data[, -1] %*% (w_p %*% A[, 1]))
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
      result<-getMLE(data_str[, -1],data_str[,1],1/prob1)
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
      result<-getMLE(data_p[, -1],data_p[, 1],rep(1/length(data_p[,1]),length(data_p[,1])))
      par_p <- result$par
      # prob
      p <- 1 - 1 / (1 + exp(data[, -1] %*% par_p))
      p_p <- p[index]
      w_p <- p_p * (1 - p_p)
      w_p <- solve(t(data_p[, -1]) %*% (data_p[, -1] * w_p))
      sub_pro <- sqrt((data[, 1] - p)^2 * rowSums((data[, -1] %*% w_p)^2))
      sub_pro <- sub_pro/sum(sub_pro)
      V <- cov(data_p[, -1] %*% w_p * as.numeric(data_p[, 1] - p_p))
      A <- svd(V, d, d)$u
      Y <- as.numeric(data[, 1] - p) * (data[, -1] %*% (w_p %*% A[, 1]))
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
        if(length((qua[j] + 1):qua[j + 1]) == 1) {indexj=qua[j + 1]}
        else{indexj <- sample((qua[j] + 1):qua[j + 1], nj, replace = TRUE, prob = order_pro[(qua[j] + 1):qua[j + 1]])}
        index <- c(index, indexj)
        prob1 <- c(prob1, nj * order_pro[indexj] / (sum(order_pro[(qua[j] + 1):qua[j + 1]])))
      }
      data_str <- data[order[index], ]
      result<-getMLE(data_str[, -1],data_str[,1],1/prob1)
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
      result<-getMLE(data_p[, -1],data_p[, 1],rep(1/length(data_p[,1]),length(data_p[,1])))
      par_p <- result$par
      # prob
      p <- 1 - 1 / (1 + exp(data[, -1] %*% par_p))
      p_p <- p[index]
      w_p <- p_p * (1 - p_p)
      w_p <- solve(t(data_p[, -1]) %*% (data_p[, -1] * w_p))
      V <- cov(data_p[, -1] %*% w_p * as.numeric(data_p[, 1] - p_p))
      A <- svd(V, d, d)$u
      Y <- as.numeric(data[, 1] - p) * (data[, -1] %*% (w_p %*% A[, 1]))
      order <- order(Y)
      sub_pro <- rep(1 / N, N)
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
      result<-getMLE(data_str[, -1],data_str[,1],1/prob1)
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
  duration <- difftime(end_time, start_time, units = "secs")
  time <- as.numeric(duration) / sim_num
  result <- list(MSE, sum(MSE), time)
  return(result)
}


library(MASS)
set.seed(123)

cas = 1
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
    x <- mvrnorm(N, rep(0.5, d-1), Sigma)
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
  data <- cbind(y, x)
  data <- as.data.frame(data)
  
  X <- as.matrix(cbind(rep(1, N), data[, -1]))
  start_time <- Sys.time()
  par00 <- rep(0, d)
  result <- getMLE(X,y,rep(1/N,N))
  end_time <- Sys.time()
  end_time - start_time
result_table_time <- matrix(rep(0, 28), nrow = 7, ncol = 4)
result_table_time[1, 1] <- end_time - start_time
par0 <- result$par

result <- sub_est_sim(data, "uniform", 500, 2500, 50, d, par0)
result_table_time[2, 1] <- result[[3]]
result <- sub_est_sim(data, "osmac", 500, 2500, 50, d, par0)
result_table_time[3, 1] <- result[[3]]

result_table_time

col <- c(5, 10, 50, 100)
row <- c("stratified-uni", "stratified-osmac")
for (i in row) {
  for (m in col) {
    result <- sub_est_sim(data, i, 500, 2500, 50, d, par0, m)
    i1 <- which(row == i) + 3
    j1 <- which(col == m)
    result_table_time[i1, j1] <- result[[3]]
  }
}

result <- sub_est_sim(data, "o-stratified-uni", 500, 2500, 50, d, par0)
result_table_time[6, 1] <- result[[3]]

result <- sub_est_sim(data, "o-stratified-osmac", 500, 2500, 50, d, par0)
result_table_time[7, 1] <- result[[3]]

result_table_time

result_table_time = round(result_table_time, 4)

write.csv(result_table_time, "time.csv")

