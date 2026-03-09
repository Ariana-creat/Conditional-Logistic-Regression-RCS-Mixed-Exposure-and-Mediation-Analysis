#(1)Mixture Exposure Analysis
#(1.1)WQS
library(gWQS)
library(survival)
library(ggplot2)

  #输出路径与文件前缀
out_dir <- "WQS_results"
prefix  <- "PP_WQS"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

  #数据准备
d <- dat_ana
d$PP <- ifelse(d$ICP == 1 | d$CPP == 1, 1L, 0L)

d$Subclass <- factor(d$Subclass)
if ("Delive_mode" %in% names(d)) d$Delive_mode <- factor(d$Delive_mode)
if ("tobacco" %in% names(d))     d$tobacco     <- factor(d$tobacco)

for (v in exposure_SG) d[[v]][!is.finite(d[[v]])] <- NA

vars_need <- c("PP", "Subclass", exposure_SG, covars)
d_wqs <- d[complete.cases(d[, vars_need]), vars_need]

  #权重图（变量名简写、英文坐标、Times New Roman 12pt、保留颜色渐变）
  exposure_labels <- c(
    "lnSCN_adj"   = "SCN\u207B",
    "lnIodie_adj" = "Iodie",
    "lnCLO3_adj"  = "ClO\u2083\u207B",
    "lnCLO4_adj"  = "ClO\u2084\u207B",
    "lnTCAA_adj"  = "TCAA"
  )
  label_exposure <- function(x) ifelse(x %in% names(exposure_labels), exposure_labels[x], x)

save_weight_plot <- function(w_df, file_png, title = NULL,
                             y_max = NULL, digits = 2) {
  stopifnot(all(c("mix_name", "mean_weight") %in% names(w_df)))
  w_df <- w_df[order(w_df$mean_weight, decreasing = TRUE), ]
  w_df$mix_label <- label_exposure(as.character(w_df$mix_name))
  w_df$mix_label <- factor(w_df$mix_label, levels = unique(w_df$mix_label))
  w_df$rank <- seq_len(nrow(w_df))
  
  if (is.null(y_max)) {
    y_max <- max(w_df$mean_weight, na.rm = TRUE) * 1.25
    if (!is.finite(y_max) || y_max <= 0) y_max <- 1
  }
  w_df$lab <- sprintf(paste0("%.", digits, "f"), w_df$mean_weight)
  
  p <- ggplot(w_df, aes(x = mix_label, y = mean_weight, fill = -rank)) +
    geom_col(width = 0.38) +
    geom_text(aes(label = lab),
              vjust = -0.6, size = 5, fontface = "bold") +
    scale_y_continuous(limits = c(0, y_max),
                       expand = expansion(mult = c(0, 0))) +
    scale_fill_gradient(low = "#DDE4EF", high = "#004C98", guide = "none") +
    labs(x = NULL, y = "Weight", title = title) +
    theme_classic(base_size = 12) +
    theme(
      text = element_text(family = "Times New Roman", size = 12),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 12),
      axis.title.y = element_text(face = "bold"),
      axis.text.x = element_text(face = "bold"),
      axis.line = element_line(linewidth = 1),
      axis.ticks = element_line(linewidth = 1),
      axis.ticks.length = unit(4, "pt")
    )
  
  ggsave(filename = file_png, plot = p, width = 8, height = 6, dpi = 300)
}

  #跑单个方向WQS
run_wqs_one_direction <- function(b1_pos_flag, d_wqs, exposure_SG, covars,
                                  out_dir, prefix,
                                  q = 4, validation = 0.6, b = 200, seed = 20260126) {
  dir_tag <- ifelse(b1_pos_flag, "pos", "neg")

  fml_wqs <- as.formula(paste0("PP ~ wqs + ", paste(covars, collapse = " + ")))
  
  set.seed(seed)
  fit <- gwqs(
    formula    = fml_wqs,
    mix_name   = exposure_SG,
    data       = d_wqs,
    q          = q,
    validation = validation,
    b          = b,
    b1_pos     = b1_pos_flag,
    family     = "binomial",
    seed       = seed
  )
  
    #报错后增
  stopifnot(length(fit$validation_rows) == nrow(d_wqs))
  n_valid <- sum(fit$validation_rows)
  stopifnot(n_valid == length(fit$wqs))
  d_wqs$WQS <- NA_real_
  d_wqs$WQS[fit$validation_rows] <- fit$wqs
  d_valid <- d_wqs[fit$validation_rows, , drop = FALSE]
  
    #条件logistic
  fml_clogit <- as.formula(paste0(
    "PP ~ WQS + ", paste(covars, collapse = " + "),
    " + strata(Subclass)"
  ))
  m <- clogit(fml_clogit, data = d_valid, method = "exact")
  
  b_hat <- coef(m)["WQS"]
  se    <- sqrt(vcov(m)["WQS","WQS"])
  OR    <- exp(b_hat)
  L     <- exp(b_hat - 1.96 * se)
  U     <- exp(b_hat + 1.96 * se)
  pval  <- summary(m)$coefficients["WQS", "Pr(>|z|)"]
  
  res <- data.frame(
    direction = dir_tag,
    OR = OR,
    CI_low = L,
    CI_high = U,
    p_value = pval,
    n_used = m$n,
    events = m$nevent,
    n_validation = nrow(d_valid),
    stringsAsFactors = FALSE
  )
  res$OR_CI <- sprintf("%.2f (%.2f, %.2f)", res$OR, res$CI_low, res$CI_high)
  
    
  #权重表
  w <- fit$final_weights
  w <- w[order(-w$mean_weight), ]
  w$direction <- dir_tag
  
    #保存权重
  write.csv(
    w,
    file = file.path(out_dir, paste0(prefix, "_weights_", dir_tag, ".csv")),
    row.names = FALSE
  )
  
    #保存权重图
  save_weight_plot(
    w_df = w,
    file_png = file.path(out_dir, paste0(prefix, "_weights_", dir_tag, ".png")),
    title = paste0(prefix, " weights (", dir_tag, ")")
  )
  
  list(fit = fit, result = res, weights = w, d_valid = d_valid, model = m)
}

  #跑正向负向
b_boot <- 500
seed0  <- 20260126

pos_out <- run_wqs_one_direction(
  b1_pos_flag = TRUE,
  d_wqs = d_wqs, exposure_SG = exposure_SG, covars = covars,
  out_dir = out_dir, prefix = prefix,
  q = 4, validation = 0.6, b = b_boot, seed = seed0
)

neg_out <- run_wqs_one_direction(
  b1_pos_flag = FALSE,
  d_wqs = d_wqs, exposure_SG = exposure_SG, covars = covars,
  out_dir = out_dir, prefix = prefix,
  q = 4, validation = 0.6, b = b_boot, seed = seed0
)

  #保存主结果
results_wqs <- rbind(pos_out$result, neg_out$result)

write.csv(
  results_wqs,
  file = file.path(out_dir, paste0(prefix, "_clogit_results_pos_neg.csv")),
  row.names = FALSE
)

results_wqs

#(1.2)Q-gcomp
library(qgcomp)
library(survival)

fml_qg <- as.formula(paste0(
  "PP ~ ", paste(exposure_SG, collapse = "+")
))
fit_qg <- qgcomp.noboot(
  f = fml_qg,
  data = d,
  family = binomial(),
  q =4,
  bayes = TRUE #报错（变量共线）后官方推荐的设置贝叶斯处理共线
)
summary(fit_qg)

w_pos <- fit_qg$pos.weights
w_neg <- fit_qg$neg.weights
if (length(w_pos)>0) w_pos <- 
comp_weights <- rbind(
  data.frame(exposure = names(w_pos),
             weight = as.numeric(w_pos),
             direction = "positive"),
  data.frame(exposure = names(w_neg),
             weight = as.numeric(w_neg),
             direction = "negetive")
)
comp_weights
write.csv(
  comp_weights,
  file = "comp_weights_qgcomp.csv",
  row.names = FALSE
)
    #图用python画的

  #提取结果
psi <- fit_qg$psi
se  <- sqrt(fit_qg$var.psi)
OR  <- exp(psi)
L   <- exp(psi - 1.96*se)
U   <- exp(psi + 1.96*se)
p   <- 2*pnorm(-abs(psi/se))

res_qg <- data.frame(
  OR = OR, CI_low = L, CI_high = U, p_value = p,
  n_used = nrow(d), events = sum(d$PP == 1),
  q = 4,
  stringsAsFactors = FALSE
)
res_qg$OR_CI <- sprintf("%.2f (%.2f, %.2f)", res_qg$OR, res_qg$CI_low, res_qg$CI_high)
write.csv(
  res_qg,
  file = "res_qg.csv",
  row.names = FALSE
)

#(1.3)BKMR 
library(aBKMR)
library(bkmr)
library(ggplot2)
library(dplyr)
library(glue)
library(ggtext)

vars_need <- unique(c(outcome, exposure_SG, covars))
d_bkmr <- d[, vars_need, drop = FALSE]
#分类变量处理
d_bkmr$age_c        <- factor(d_bkmr$age_c, levels = c("<28", ">=28"))
d_bkmr$Parity_c     <- factor(d_bkmr$Parity_c, levels = c(1, 2))
d_bkmr$Delive_mode  <- factor(d_bkmr$Delive_mode, levels = c("顺产", "非顺产"))
d_bkmr$GWG          <- factor(d_bkmr$GWG, levels = c(1, 2))
d_bkmr$tobacco      <- factor(d_bkmr$tobacco, levels = c("no", "yes"))
d_bkmr$predclass    <- factor(d_bkmr$predclass, levels = c("low", "high"))
d_bkmr$MENSESAGE_c  <- factor(d_bkmr$MENSESAGE_c, levels = c(1, 2))

#只保留完整病例
d_bkmr <- d_bkmr[complete.cases(d_bkmr), ]
y <- d_bkmr$PP
Z <- as.matrix(d_bkmr[, exposure_SG])
X <- model.matrix(
  ~ age_c + Parity_c + Delive_mode + GWG + tobacco + predclass + MENSESAGE_c,
  data = d_bkmr
)[, -1, drop = FALSE]
k_id <- sam_py_r(R = Z, nd = 100, num_nn = 100, w = FALSE)
km <- a_kmbayes(y = y, Z = Z, X = X, knots = Z[k_id,],
                est.h = TRUE, varsel = TRUE,
                iter = 10000, family = "binomial")
#PIP:后验包含概率
PIP_new <- ExtractPIPs(km) %>% arrange(desc(PIP))
ggplot()+geom_col(data = PIP_new, aes(x = reorder(variable, PIP, decreasing = T), y = PIP), fill = '#48C5B9')+ 
  ylab("Estimated PIPs")+ xlab('Exposures')+ theme_bw()

#联合效应
quants = seq(0.25, 0.75, 0.01)

newz = lapply(1:ncol(Z), function(i){
  qf = data.frame(X = quantile(Z[, i], quants))
  names(qf) = names(data.frame(Z))[i]
  return(qf)
}) %>% bind_cols()

overall_res = a_OverallRiskSummaries_vary(km, newz = newz, data.comps = km$data.comps)
overall_res$risk_overall$quants = seq(0.25, 0.75, 0.01)

je_res = com_je(overall_res$preds_e$postmean, overall_res$preds_e$postvar, X = quants*100)

ggplot() +
  geom_hline(yintercept = 0, linetype = 2, color = "red") +
  geom_pointrange(data = overall_res$risk_overall, 
                  aes(quants, est_p, ymin = est_p - 1.96*sd_p, ymax = est_p + 1.96*sd_p), 
                  color = '#48C5B9', shape = 16) +
  ylab("Relative risk of outcome") + xlab("Quantiles of exposures") +
  theme_bw()

#单变量暴露-反应曲线
singlevar_res = a_PredictorResponseUnivar(km, quants = seq(0, 1, 0.1), method = "approx", data.comps = km$data.comps) 
singlevar_data = lapply(1:length(singlevar_res), function(i){
  df = singlevar_res[[i]]$risk_overall %>% mutate(var_est = glue("{vars}"))
  return(df) 
}) %>% bind_rows() 
ggplot()+ 
  geom_smooth(data = singlevar_data, 
              aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se), 
              color = '#48C5B9', stat = "identity")+
  facet_wrap(~ var_est)+ylab("Relative risk of outcome")+xlab("Exposures")
