library(numDeriv)

shl_loss <- function(y, f_x) 0.5 * pmax(0, 1 - y * f_x)^2

svm_train <- function(X, Y, lambda = 0.1, weights = NULL, init_p = NULL) {
  d <- ncol(X)
  p0 <- if (is.null(init_p)) rep(0, d + 1) else init_p

  obj <- function(p) {
    w <- p[1:d]
    b <- p[d+1]
    f_x <- as.numeric(X %*% w + b)
    losses <- shl_loss(Y, f_x)
    
    if (!is.null(weights)) {
      emp_risk <- mean(losses * weights)
    } else {
      emp_risk <- mean(losses)
    }
    
    return(emp_risk + 0.5 * lambda * sum(w^2))
  }
  
  res <- optim(p0, obj, method = "BFGS", control = list(maxit = 200))
  return(res$par)
}

calc_svm_inf_vecs <- function(X, Y, params, lambda) {
  N <- nrow(X); d <- ncol(X); w <- params[1:d]; b <- params[d+1]
  f_vals <- as.numeric(X %*% w + b)
  diff <- 1 - Y * f_vals
  active <- diff > 0
  X_aug <- cbind(X, 1)
  H <- lambda * diag(c(rep(1, d), 0)) + (t(X_aug[active, ]) %*% X_aug[active, ]) / N
  H_inv <- solve(H + diag(1e-7, d + 1))
  
  inf_vecs <- matrix(0, N, d + 1)
  for(i in which(active)) {
    grad_i <- (-diff[i] * Y[i]) * X_aug[i, ]
    inf_vecs[i, ] <- -H_inv %*% grad_i
  }
  return(inf_vecs)
}


expand_features <- function(X_raw) {
  d <- ncol(X_raw)
  X_poly <- cbind(X_raw, X_raw^2)
  
  interaction_list <- list()
  count <- 1
  for(i in 1:(d-1)) {
    for(j in (i+1):d) {
      interaction_list[[count]] <- X_raw[, i] * X_raw[, j]
      count <- count + 1
    }
  }
  
  X_full <- do.call(cbind, c(list(X_poly), interaction_list))
  return(X_full)
}

stratified_sampling <- function(base_probs, order_idx, m, n_sub) {
  ord_probs <- base_probs[order_idx]
  qua <- floor(N * seq(0, 1, by = 1/m))
  final_idx <- c(); final_ht_w <- c()
  
  for (j in 1:m) {
    range <- (qua[j] + 1):qua[j + 1]
    sum_p_j <- sum(ord_probs[range])
    nj <- round(n_sub * sum_p_j)
    
    if(nj > 0) {
      idx_j <- sample(range, nj, replace = TRUE, prob = ord_probs[range])
      pi_eff <- nj * (ord_probs[idx_j] / sum_p_j)
      
      final_idx <- c(final_idx, order_idx[idx_j])
      final_ht_w <- c(final_ht_w, 1 / pi_eff)
    }
  }
  return(list(idx=final_idx, w=final_ht_w))
}

cumulative_stratified_sampling <- function(base_probs, order_idx, n_sub) {
  N <- length(base_probs)
  ord_probs <- base_probs[order_idx]
  pro_sum <- cumsum(ord_probs) / sum(ord_probs)
  
  thresholds <- (1:(n_sub - 1)) / n_sub

  indices <- findInterval(thresholds, pro_sum)
  qua <- c(0, indices, N)
  
  final_idx <- c()
  final_ht_w <- c()
  
  for (j in 1:n_sub) {
    range_idx <- (qua[j] + 1):qua[j + 1]
    if(length(range_idx) == 0) next 
    nj <- 1
    
    if(length(range_idx) == 1) {
      idx_j_in_range <- 1
    } else {
      idx_j_in_range <- sample(1:length(range_idx), nj, 
                               replace = TRUE, 
                               prob = ord_probs[range_idx])
    }
    
    selected_global_idx <- order_idx[range_idx[idx_j_in_range]]

    sum_p_j <- sum(ord_probs[range_idx])
    pi_eff <- nj * (base_probs[selected_global_idx] / sum_p_j)
    
    final_idx <- c(final_idx, selected_global_idx)
    final_ht_w <- c(final_ht_w, 1 / pi_eff)
  }
  
  return(list(idx = final_idx, w = final_ht_w))
}


sub_est_sim <- function(X_expanded, Y_full, method, n0 = 500, n_sub = 2500, sim_num = 20, d = 66, par0, m = 10) {
  N <- length(Y_full)
  ## uniform subsampling
  if (method == "uniform") {
    par_uni <- matrix(0, nrow = sim_num, ncol = d)
    start_time <- Sys.time()
    for (i in 1:sim_num) {
      idx_u <- sample(1:N, n_sub, replace = TRUE)
      par_uni[i, ] <- svm_train(X_expanded[idx_u, ], Y_full[idx_u], lambda = 0.1)
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
      idx_p <- sample(1:N, n0)
      theta_plt <- svm_train(X_expanded[idx_p, ], Y_full[idx_p], lambda = 0.1)
      inf_vecs <- calc_svm_inf_vecs(X_expanded, Y_full, theta_plt, lambda = 0.1)
      
      inf_norms <- sqrt(rowSums(inf_vecs^2))
      pi_osmac <- 0.9 * (inf_norms/sum(inf_norms)) + 0.1/N
      
      idx_o <- sample(1:N, n_sub, replace = TRUE, prob = pi_osmac)
      par_sub[i, ] <- svm_train(X_expanded[idx_o, ], Y_full[idx_o], 
                                weights = 1/N/pi_osmac[idx_o], lambda = 0.1)
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
      idx_p <- sample(1:N, n0)
      theta_plt <- svm_train(X_expanded[idx_p, ], Y_full[idx_p], lambda = 0.1)
      inf_vecs <- calc_svm_inf_vecs(X_expanded, Y_full, theta_plt, lambda = 0.1)
      
      V_cov <- cov(inf_vecs[idx_p, ])
      A <- svd(V_cov)$u[, 1, drop=FALSE]
      proj_val <- inf_vecs %*% A
      order_idx <- order(proj_val)
      
      s_uni <- stratified_sampling(rep(1/N, N), order_idx, m, n_sub)
      par_str[i, ] <- svm_train(X_expanded[s_uni$idx, ], Y_full[s_uni$idx], weights = s_uni$w/N*n_sub)
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
      idx_p <- sample(1:N, n0)
      theta_plt <- svm_train(X_expanded[idx_p, ], Y_full[idx_p], lambda = 0.1)
      inf_vecs <- calc_svm_inf_vecs(X_expanded, Y_full, theta_plt, lambda = 0.1)
      
      V_cov <- cov(inf_vecs[idx_p, ])
      A <- svd(V_cov)$u[, 1, drop=FALSE]
      proj_val <- inf_vecs %*% A
      order_idx <- order(proj_val)
      
      inf_norms <- sqrt(rowSums(inf_vecs^2))
      pi_osmac <- 0.9 * (inf_norms/sum(inf_norms)) + 0.1/N
      
      s_osmac <- stratified_sampling(pi_osmac, order_idx, m, n_sub)
      par_str[i, ] <- svm_train(X_expanded[s_osmac$idx, ], Y_full[s_osmac$idx], weights = s_osmac$w/N*n_sub)
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
      idx_p <- sample(1:N, n0)
      theta_plt <- svm_train(X_expanded[idx_p, ], Y_full[idx_p], lambda = 0.1)
      inf_vecs <- calc_svm_inf_vecs(X_expanded, Y_full, theta_plt, lambda = 0.1)
      
      V_cov <- cov(inf_vecs[idx_p, ])
      A <- svd(V_cov)$u[, 1, drop=FALSE]
      proj_val <- inf_vecs %*% A
      order_idx <- order(proj_val)
      
      s_uni <- cumulative_stratified_sampling(rep(1/N, N), order_idx, n_sub)
      par_str[i, ] <- svm_train(X_expanded[s_uni$idx, ], Y_full[s_uni$idx], weights = s_uni$w/N*n_sub)
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
      idx_p <- sample(1:N, n0)
      theta_plt <- svm_train(X_expanded[idx_p, ], Y_full[idx_p], lambda = 0.1)
      inf_vecs <- calc_svm_inf_vecs(X_expanded, Y_full, theta_plt, lambda = 0.1)
      
      V_cov <- cov(inf_vecs[idx_p, ])
      A <- svd(V_cov)$u[, 1, drop=FALSE]
      proj_val <- inf_vecs %*% A
      order_idx <- order(proj_val)
      
      inf_norms <- sqrt(rowSums(inf_vecs^2))
      pi_osmac <- 0.9 * (inf_norms/sum(inf_norms)) + 0.1/N
      
      s_osmac <- cumulative_stratified_sampling(pi_osmac, order_idx, n_sub)
      par_str[i, ] <- svm_train(X_expanded[s_osmac$idx, ], Y_full[s_osmac$idx], weights = s_osmac$w/N*n_sub)
    }
    end_time <- Sys.time()
    MSE <- rep(0, d)
    Bias <- rep(0, d)
    for (i in 1:d) {
      MSE[i] <- mean((par_str[, i] - par0[i])^2)
      Bias[i] <- mean(par_str[, i] - par0[i])
    }
  }


  time <- as.numeric(difftime(end_time, start_time, units = "secs"))/ sim_num
  result <- list(MSE, sum(MSE), time)
  return(result)
}

library(MASS)
library(mvtnorm)
for (cas in 1:4) {
N=500000
d=10
set.seed(123)
sds <- 1 / (1:(d))
Corr_mat <- matrix(0.5, nrow = (d), ncol = (d))
diag(Corr_mat) <- 1
Sigma <- diag(sds) %*% Corr_mat %*% diag(sds)
x <- mvrnorm(N, 0.5*(d:1), Sigma)

X_raw <- x
X_expanded <- expand_features(X_raw)
p_dim <- ncol(X_expanded)

if(cas ==1){
  weights <- 0.1 * rep(1,d)
term_lin <- X_raw %*% weights
term_inter <- 0.1 * (X_raw[,1]*X_raw[,2]) 
+ 0.1 * (X_raw[,3]*X_raw[,4]) + 0.1 * (X_raw[,5]*X_raw[,6]) + 0.1 *(X_raw[,7]*X_raw[,8])+ 0.1 *(X_raw[,9]*X_raw[,10])
true_logic <- term_lin + term_inter - mean(term_lin + term_inter)
Y_full <- ifelse(true_logic > 0, 1, -1)
sum(Y_full==1)
}

if(cas ==2){
  weights <- 0.1 *(1:d)
term_lin <- X_raw %*% weights
term_inter <- 0.1 * (X_raw[,1]*X_raw[,2]) 
+ 0.1 * (X_raw[,3]*X_raw[,4]) + 0.1 * (X_raw[,5]*X_raw[,6]) + 0.1 *(X_raw[,7]*X_raw[,8])+ 0.1 *(X_raw[,9]*X_raw[,10])
true_logic <- term_lin + term_inter - mean(term_lin + term_inter)
Y_full <- ifelse(true_logic > 0, 1, -1)
}

if(cas == 3){
 weights <- 0.1 * rep(1,d)
term_lin <- X_raw %*% weights
term_inter <- 0.1 * sin(X_raw[,1]*X_raw[,2]) 
+ 0.1 * sin(X_raw[,3]*X_raw[,4]) + 0.1 * sin(X_raw[,5]*X_raw[,6]) + 0.1 *exp(X_raw[,7]*X_raw[,8])+ 0.1 *exp(X_raw[,9]*X_raw[,10])
true_logic <- term_lin + term_inter - mean(term_lin + term_inter)
Y_full <- ifelse(true_logic > 0, 1, -1)
}

if(cas == 4){
 weights <- 0.1 * rep(1,d)
term_lin <- X_raw %*% weights
term_inter <- 0.1*sign(X_raw[,1]*X_raw[,2]) 
true_logic <- term_lin + term_inter - mean(term_lin + term_inter)
Y_full <- ifelse(true_logic > 0, 1, -1)
}


t_full <- system.time({
  theta_full <- svm_train(X_expanded, Y_full, lambda = 0.1)
})
theta_full

result_table_mse <- matrix(rep(0, 16), nrow = 4, ncol = 4)
result_table_time <- matrix(rep(0, 16), nrow = 4, ncol = 4)
col <- c(1000, 1500, 2000, 2500)
row <- c("uniform", "osmac", "stratified-uni", "stratified-osmac")
for (i in row) {
  for (j in col) {
    result <- sub_est_sim(X_expanded, Y_full, i, 500, j, 1000, p_dim+1, par0=theta_full)
    i1 <- which(row == i)
    j1 <- which(col == j)
    result_table_mse[i1, j1] <- result[[2]]
    result_table_time[i1, j1] <- result[[3]]
  }
  print(i)
}
result_table_mse
#result_table_time
write.csv(result_table_mse, paste0("SVM1_case", cas, ".csv"))
#write.csv(result_table_time, paste0("SVM_time", cas, ".csv"))
}
