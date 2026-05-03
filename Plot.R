library(reshape2)
library(ggplot2)

## logistic model
n <- rep(c(1000, 1500, 2000, 2500), 4)
method <- c(rep("UNIF", 4), rep("OPT", 4), rep("MVRS-U", 4), rep("MVRS-O", 4))

for (j in 1:4) {
  file_name <- paste0("logistic_regression_case", j, ".csv")
  
  data0 <- read.csv(file_name)
  data0 <- data0[, -1]
  mse <- c(as.numeric(data0[1, ]), as.numeric(data0[2, ]),
           as.numeric(data0[3, ]), as.numeric(data0[4, ]))
  data <- data.frame(n = n, method = method, mse = mse)
  
  ggplot(data = data, aes(x = n, y = mse, group = method, color = method)) +
    geom_point(aes(shape = method), size = 3.5) +
    scale_shape_manual(values = c(16, 17, 15, 18),
                       limits = c("UNIF", "MVRS-U", "OPT", "MVRS-O")) +
    scale_color_manual(values = c("#FF6B35", "#D9534F", "#5CB85C", "#2980B9"),
                       limits = c("UNIF", "MVRS-U", "OPT", "MVRS-O")) +
    geom_line() +
    xlab("Subsample Size") +
    ylab("MSE") +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = 15),  
          axis.title.y = element_text(size = 15),  
          axis.text.x = element_text(size = 15),  
          axis.text.y = element_text(size = 15),
          legend.position = if (j == 1) c(0.8, 0.7) else "none",
          legend.box.background = element_rect(color = "black"),
          text = element_text(family = "serif", size = 15)) +
    scale_x_continuous(limits = c(1000, 2500), breaks = seq(1000, 2500, 500))
  
  file_name <- paste0("log_", j, ".pdf")
  ggsave(
    filename = file_name, 
    width = 5.4,              
    height = 3.6,             
    units = "in",           
    dpi = 300               
  )
}


# MSE estimation
n <- rep(c(1000,1500,2000,2500), 4)
method <- c(rep("MVRS-U", 4), rep("MVRS-O", 4), rep("MVRS-U(est)", 4), rep("MVRS-O(est)", 4))

for (j in 1:4) {
  data0 <- read.csv(paste0("emp_logistic_regression_case", j, ".csv"))
  data0 <- data0[, -1]
  data1 <- read.csv(paste0("est_logistic_regression_case", j, ".csv"))
  data1 <- data1[, -1]
  mse <- c(as.numeric(data0[1, ]), as.numeric(data0[2, ]),
           as.numeric(data1[1, ]), as.numeric(data1[2, ]))
  data <- data.frame(n = n, method = method, mse = mse)
  
  ggplot(data = data, aes(x = n, y = mse, group = method, color = method)) +
    geom_point(aes(shape = method), size = 3.5) +
    scale_shape_manual(values = c(17, 18, 15, 16),
                       limits = c("MVRS-U", "MVRS-U(est)", "MVRS-O", "MVRS-O(est)")) +
    scale_color_manual(values = c("#D9534F", "#D9534F", "#2980B9", "#2980B9"),
                       limits = c("MVRS-U", "MVRS-U(est)", "MVRS-O", "MVRS-O(est)")) +
    geom_line(aes(linetype = method)) +
    scale_linetype_manual(
      values = c("solid", "dashed", "solid", "dashed"),
      limits = c("MVRS-U", "MVRS-U(est)", "MVRS-O", "MVRS-O(est)")) +
    xlab("Subsample Size") +
    ylab("MSE") +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = 15),  
          axis.title.y = element_text(size = 15),  
          axis.text.x = element_text(size = 15),  
          axis.text.y = element_text(size = 15),
          legend.position = if (j == 1) c(0.8, 0.7) else "none",
          legend.box.background = element_rect(color = "black"),
          text = element_text(family = "serif", size = 15)) +
    scale_x_continuous(limits = c(1000, 2500), breaks = seq(1000, 2500, 500))
  ggsave(
    filename = paste0("log_", j, "_est.pdf"), 
    width = 5.4,              
    height = 3.6,             
    units = "in",           
    dpi = 300               
  )
}


# number of strata
n <- rep(c(1, 5, 10, 50, 100), 2)
method <- c(rep("MVRS-U", 5), rep("MVRS-O", 5))
for (j in 1:4) {
  file_name <- paste0("logistic_regression_case", j, "(str_num_n=1000).csv")
  
  data0 <- read.csv(file_name)
  data0 <- data0[, -1]
  mse <- c(as.numeric(data0[1, ]), as.numeric(data0[2, ]))
  data <- data.frame(n = n, method = method, mse = mse)


  ggplot(data = data, aes(x = n, y = mse, group = method, color = method)) +
  geom_point(aes(shape = method), size = 3.5) +
  scale_shape_manual(values = c(17, 18),
                     limits = c("MVRS-U", "MVRS-O")) +
  scale_color_manual(values = c("#D9534F", "#2980B9"),
                     limits = c("MVRS-U", "MVRS-O")) +
  geom_line() +
  xlab("Number of strata") +
  ylab("MSE") +
  theme_bw() +
  theme(panel.grid.major = element_line(colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 15),  
        axis.title.y = element_text(size = 15),  
        axis.text.x = element_text(size = 15),  
        axis.text.y = element_text(size = 15),
        legend.position = if (j == 1) c(0.60, 0.85) else "none",
        legend.box.background = element_rect(color = "black"),
        legend.direction = "horizontal",
        text = element_text(family = "serif", size = 17)) +
  scale_x_continuous(breaks = c(1, 10, 50, 100))
  ggsave(
    filename = paste0("num_", j, ".pdf"), 
    width = 5.4,              
    height = 2.2,             
    units = "in",           
    dpi = 300 
  )
}


## SVM
n <- rep(c(1000, 1500, 2000, 2500), 4)
method <- c(rep("UNIF", 4), rep("OPT", 4), rep("MVRS-U", 4), rep("MVRS-O", 4))

for (j in 1:4) {
  file_name <- paste0("SVM1_case", j, ".csv")
  
  data0 <- read.csv(file_name)
  data0 <- data0[, -1]
  mse <- c(as.numeric(data0[1, ]), as.numeric(data0[2, ]),
           as.numeric(data0[3, ]), as.numeric(data0[4, ]))
  data <- data.frame(n = n, method = method, mse = mse)
  
  ggplot(data = data, aes(x = n, y = mse, group = method, color = method)) +
    geom_point(aes(shape = method), size = 3.5) +
    scale_shape_manual(values = c(16, 17, 15, 18),
                       limits = c("UNIF", "MVRS-U", "OPT", "MVRS-O")) +
    scale_color_manual(values = c("#FF6B35", "#D9534F", "#5CB85C", "#2980B9"),
                       limits = c("UNIF", "MVRS-U", "OPT", "MVRS-O")) +
    geom_line() +
    xlab("Subsample Size") +
    ylab("MSE") +
    theme_bw() +
    theme(panel.grid.major = element_line(colour = NA),
          panel.background = element_rect(fill = "transparent", colour = NA),
          plot.background = element_rect(fill = "transparent", colour = NA),
          panel.grid.minor = element_blank(),
          axis.title.x = element_text(size = 15),  
          axis.title.y = element_text(size = 15),  
          axis.text.x = element_text(size = 15),  
          axis.text.y = element_text(size = 15),
          legend.position = if (j == 1) c(0.8, 0.7) else "none",
          legend.box.background = element_rect(color = "black"),
          text = element_text(family = "serif", size = 15)) +
    scale_x_continuous(limits = c(1000, 2500), breaks = seq(1000, 2500, 500))
  
  file_name <- paste0("svm_", j, ".pdf")
  ggsave(
    filename = file_name, 
    width = 5.4,              
    height = 3.6,             
    units = "in",           
    dpi = 300               
  )
}




