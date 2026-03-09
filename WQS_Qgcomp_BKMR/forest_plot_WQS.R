# Forest plot: WQS clogit results (正向、负向分开画)
# PP_WQS_clogit_results_pos_neg.csv

library(ggplot2)

df <- read.csv("/Users/lixinrui/Desktop/混合暴露分析/output_csv/Mix_results/PP_WQS_clogit_results_pos_neg.csv")

# 分别画正向和负向
for (dir_tag in c("pos", "neg")) {
  sub <- df[df$direction == dir_tag, , drop = FALSE]
  if (nrow(sub) == 0) next
  sub$y_label <- "WQS Index"
  sub$label_text <- sprintf("OR (95%% CI)\n%s\nP value = %.2f", sub$OR_CI, sub$p_value)

  x_lim_hi <- max(4, sub$CI_high * 1.2)
  x_text <- x_lim_hi * 0.85

  p <- ggplot(sub, aes(x = OR, y = y_label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "black", linewidth = 0.6) +
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.2, linewidth = 1, color = "black") +
    geom_point(size = 4, color = "black", shape = 16) +
    scale_x_continuous(
      name = "Adjusted OR (95% CI)",
      breaks = seq(1, ceiling(x_lim_hi), 0.5),
      limits = c(0.8, x_lim_hi * 1.2),
      expand = c(0, 0)
    ) +
    geom_text(aes(x = x_text, label = label_text), hjust = 0, vjust = 0.5, size = 3.5, inherit.aes = TRUE) +
    labs(y = NULL) +
    theme_classic(base_size = 12) +
    theme(
      axis.line = element_line(color = "black"),
      axis.text = element_text(color = "black", size = 12),
      axis.title = element_text(color = "black", size = 12),
      panel.grid = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    ) +
    coord_cartesian(clip = "off")

  fn <- paste0("forest_plot_WQS_", dir_tag)
  ggsave(paste0(fn, ".png"), p, width = 8, height = 3, dpi = 300)
  ggsave(paste0(fn, ".pdf"), p, width = 8, height = 3)
  print(p)
}
