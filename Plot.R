library(reshape2)
library(ggplot2)

## logistic model
n <- rep(c(200, 500, 800, 1000), 4)
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
          legend.position = c(.8, .7),
          legend.box.background = element_rect(color = "black"),
          text = element_text(family = "serif", size = 15)) +
    scale_x_continuous(limits = c(200, 1000), breaks = seq(200, 1000, 200))
  
  file_name <- paste0("log_", j, ".pdf")
  ggsave(
    filename = file_name, 
    width = 5.4,              
    height = 3.6,             
    units = "in",           
    dpi = 300               
  )
}



## poisson model
n <- rep(c(200, 500, 800, 1000), 4)
method <- c(rep("UNIF", 4), rep("OPT", 4), rep("MVRS-U", 4), rep("MVRS-O", 4))

for (j in 1:4) {
  file_name <- paste0("poisson_regression_case", j, ".csv")
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
          legend.position = c(.8, .7),
          legend.box.background = element_rect(color = "black"),
          text = element_text(family = "serif", size = 15)) +
    scale_x_continuous(limits = c(200, 1000), breaks = seq(200, 1000, 200))
  
  file_name <- paste0("poi_", j, ".pdf")
  ggsave(
    filename = file_name, 
    width = 5.4,              
    height = 3.6,             
    units = "in",           
    dpi = 300               
  )
}


# number of strata
library(cowplot)
library(grid)
cases <- paste0("case", 1:4)
n <- c(200, 1000)

draw_sub_plot <- function(data_subset, show_x_title = FALSE, label) {
  p <- ggplot(data = data_subset, aes(x = m, y = mse, group = method, color = method)) +
    geom_point(aes(shape = method), size = 3.5) +
    scale_shape_manual(name = label,
                       values = c(17, 18),
                       limits = c("MVRS-U", "MVRS-O")) +
    scale_color_manual(name = label,
                       values = c("#D9534F", "#2980B9"),
                       limits = c("MVRS-U", "MVRS-O")) +
    geom_line() +
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
          legend.position = c(.60, .85),
          legend.box.background = element_rect(color = "black"),
          legend.direction = "horizontal",
          text = element_text(family = "serif", size = 17)) +
    scale_x_continuous(breaks = c(1, 10, 50, 100))
  
  
  if (show_x_title) { p <- p + labs(x = NULL) } 
  else {
    p <- p + labs(x = NULL) + 
      theme(axis.title.x = element_blank(),
            axis.text.x = element_blank())
  }
  return(p)
}

plot_list <- list()
for (i in 1:4) {
  for (j in 1:2) {
    file_name <- paste0("poisson_regression_", cases[i], "(str_num_n=", n[j], ").csv")
    data_plot <- read.csv(file_name)[, -1]
    method <- c(rep("MVRS-U", 5), rep("MVRS-O", 5))
    m <- rep(c(1, 5, 10, 50, 100), 2)
    mse <- c(as.numeric(data_plot[1, ]), as.numeric(data_plot[2, ]))
    data_plot <- data.frame(m = m, mse = mse, method = method)
    label <- if (j == 1) { "n=200" } else { "n=1000" }
    
    p_temp <- draw_sub_plot(data_plot, show_x_title = (j == 2), label = label)
    
    if (j == 1) {
      p_temp <- p_temp + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    }
    
    plot_list[[paste(i, j)]] <- p_temp
  }
  p1 <- plot_list[[paste(i, 1)]]
  p2 <- plot_list[[paste(i, 2)]]
  inner_v <- cowplot::plot_grid(
    p1, p2, 
    ncol = 1, 
    align = 'v',           
    rel_heights = c(0.96, 1.04)  
  )
  strata_grob <- grid::textGrob("Number of Strata", x = 0.55, gp = gpar(fontsize = 17, fontfamily = "serif"))
  cowplot::plot_grid(
    inner_v, 
    strata_grob, 
    ncol = 1, 
    rel_heights = c(1, 0.05)
  )
  
  ggsave(
    filename = paste0("num_", i, ".pdf"), 
    width = 5.4,              
    height = 4.3,             
    units = "in",           
    dpi = 300 
  )
  
}


# MSE estimation
n <- rep(c(200, 500, 800, 1000), 4)
method <- c(rep("MVRS-U", 4), rep("MVRS-O", 4), rep("MVRS-U(est)", 4), rep("MVRS-O(est)", 4))

for (j in 1:4) {
  data0 <- read.csv(paste0("emp_poisson_regression_case", j, ".csv"))
  data0 <- data0[, -1]
  data1 <- read.csv(paste0("est_poisson_regression_case", j, ".csv"))
  data1 <- data1[, -1]
  mse <- c(as.numeric(data0[3, ]), as.numeric(data0[4, ]),
           as.numeric(data1[3, ]), as.numeric(data1[4, ]))
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
          legend.position = c(.8, .7),
          legend.box.background = element_rect(color = "black"),
          text = element_text(family = "serif", size = 15)) +
    scale_x_continuous(limits = c(200, 1000), breaks = seq(200, 1000, 200))
  ggsave(
    filename = paste0("poi_", j, "_est.pdf"), 
    width = 5.4,              
    height = 3.6,             
    units = "in",           
    dpi = 300               
  )
}






