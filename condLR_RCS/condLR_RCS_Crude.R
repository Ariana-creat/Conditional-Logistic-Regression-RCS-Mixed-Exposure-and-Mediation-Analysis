#（1）单独建一个分析数据集
dat_ana <- get("ion_cpp_analy")
dim(dat_ana)

#（2）建立变量集合(SG校正后)
exposure_SG <-c("lnSCN_adj","lnIodie_adj","lnCLO3_adj","lnCLO4_adj","lnTCAA_adj")
outcome<-c("PP","ICP","CPP")
covar<-c("age_c","Parity_c","Delive_mode","GWG","tobacco","predclass","MENSESAGE_c")
mediators <- c(
  "lnOHdG_adj","lnOHG_adj","lnHNEMA_adj","lnAlla_adj","lnBY_adj","lndiY_adj",
  "lnE1_adj","lnE2_adj","lnE3_adj","lnP4_adj","lnOH17P4_adj","lnOH21P4_adj",
  "lnTesto_adj","lnA4_adj","lnDHT_adj","lnCOH_adj","lnCTONE_adj",
  "lnCortisoF_adj","lnDeoxyCortisolS_adj","lnALDO_adj"
)
setdiff(c(exposure_SG, outcome, covar, mediators), names(dat_ana))

#（3）条件LR——continuous
dat_ana$PP <- ifelse(dat_ana$ICP == 1 | dat_ana$CPP == 1, 1L, 0L)
dat_ana$Delive_mode <- factor(dat_ana$Delive_mode)
dat_ana$tobacco     <- factor(dat_ana$tobacco)

library(survival)

# 将 infinite 转换成 NA，否则模型报错
for (v in exposure_SG) {
  dat_ana[[v]][!is.finite(dat_ana[[v]])] <- NA
}

# 构造 clogit—_1 函数（不含 covars）
clogit_1_crude <- function(exposure, data,
                           outcome = "PP",
                           strata_var = "Subclass") {
  fomula <- as.formula(paste0(
    outcome, " ~ ", exposure,
    " + strata(", strata_var, ")"
  ))
  
  # 拟合 clogit
  m <- clogit(fomula, data = data)
  
  # 提取 exposure 的 coef(β)，按 OR_×2 = 2^β 转为“每翻倍”OR
  b  <- coef(m)[exposure]
  se <- sqrt(vcov(m)[exposure, exposure])
  OR <- 2^b
  lower <- 2^(b - 1.96 * se)
  upper <- 2^(b + 1.96 * se)
  p  <- summary(m)$coefficients[exposure, "Pr(>|z|)"]
  
  data.frame(
    exposure = exposure,
    OR_doubling = OR,
    CI_low = lower,
    CI_high = upper,
    p_value = p,
    n_used = m$n,
    events = m$nevent,
    stringsAsFactors = FALSE
  )
}

# 循环跑模型
results_PP_clogit_crude <- do.call(
  rbind,
  lapply(exposure_SG, clogit_1_crude, data = dat_ana)
)

# FDR 校正（Benjamini–Hochberg）
results_PP_clogit_crude$p_FDR <- p.adjust(results_PP_clogit_crude$p_value, method = "BH")

results_PP_clogit_crude$OR_CI <- sprintf("%.2f (%.2f, %.2f)",
                                         results_PP_clogit_crude$OR_doubling,
                                         results_PP_clogit_crude$CI_low,
                                         results_PP_clogit_crude$CI_high)

results_PP_clogit_crude <- results_PP_clogit_crude[order(results_PP_clogit_crude$p_value), ]
results_PP_clogit_crude[, c("exposure","OR_CI","p_value","p_FDR","n_used","events")]

# 输出结果（OR 为暴露每翻倍时的 OR_×2 = 2^β）
write.csv(results_PP_clogit_crude, "results_PP_clogit_crude.csv", row.names = FALSE)

#（4）条件logistic——四分位（不做 covars 校正）
library(survival)

clogit_quartile_unadj <- function(exposure, data,
                                  outcome = "PP",
                                  strata_var = "Subclass") {
  
  d <- data
  stopifnot(is.character(exposure), length(exposure) == 1)
  
  x_all <- suppressWarnings(as.numeric(d[[exposure]]))
  idx_ctrl <- d[[outcome]] == 0 & is.finite(x_all)
  x_ctrl <- x_all[idx_ctrl]
  
  # 四分位切点（基于对照组）
  qs <- as.numeric(quantile(x_ctrl, probs = c(0, .25, .5, .75, 1),
                            na.rm = TRUE, type = 7))
  
  # 切点不唯一则返回 NA
  if (length(unique(qs)) < 5) {
    return(data.frame(
      exposure = exposure,
      OR_Q2 = NA, L_Q2 = NA, U_Q2 = NA, p_Q2 = NA,
      OR_Q3 = NA, L_Q3 = NA, U_Q3 = NA, p_Q3 = NA,
      OR_Q4 = NA, L_Q4 = NA, U_Q4 = NA, p_Q4 = NA,
      p_trend = NA, n_used = NA, events = NA,
      note = "Quartile cutpoints not unique",
      stringsAsFactors = FALSE
    ))
  }
  
  # 分四分位（Q1 参照）
  d$Q <- cut(x_all, breaks = qs, include.lowest = TRUE,
             labels = c("Q1","Q2","Q3","Q4"))
  d$Q <- relevel(d$Q, ref = "Q1")
  
  # 趋势变量：对照组各分位中位数映射到全体
  idx_ctrl2 <- d[[outcome]] == 0 & !is.na(d$Q) & is.finite(x_all)
  med_ctrl <- sapply(split(x_all[idx_ctrl2], d$Q[idx_ctrl2]),
                     function(z) median(z, na.rm = TRUE))
  d$Q_med <- med_ctrl[as.character(d$Q)]
  
  # 分类模型（不含 covars）
  fml_cat <- as.formula(paste0(
    outcome, " ~ Q + strata(", strata_var, ")"
  ))
  m_cat <- clogit(fml_cat, data = d)
  
  # 趋势检验（不含 covars）
  fml_trend <- as.formula(paste0(
    outcome, " ~ Q_med + strata(", strata_var, ")"
  ))
  m_trend <- clogit(fml_trend, data = d)
  
  # 提取 summary
  sm <- summary(m_cat)$coefficients
  get_row <- function(rn) {
    b  <- sm[rn, "coef"]
    se <- sm[rn, "se(coef)"]
    OR <- exp(b)
    L  <- exp(b - 1.96 * se)
    U  <- exp(b + 1.96 * se)
    p  <- sm[rn, "Pr(>|z|)"]
    c(OR = OR, L = L, U = U, p = p)
  }
  
  rows <- rownames(sm)
  out_list <- list(Q2 = NA, Q3 = NA, Q4 = NA)
  for (k in c("Q2","Q3","Q4")) {
    rn <- paste0("Q", k)  # "QQ2" "QQ3" "QQ4"
    if (rn %in% rows) out_list[[k]] <- get_row(rn)
  }
  
  p_trend <- summary(m_trend)$coefficients["Q_med", "Pr(>|z|)"]
  
  # 整理结果
  data.frame(
    exposure = exposure,
    OR_Q2 = ifelse(is.na(out_list$Q2[1]), NA, out_list$Q2["OR"]),
    L_Q2  = ifelse(is.na(out_list$Q2[1]), NA, out_list$Q2["L"]),
    U_Q2  = ifelse(is.na(out_list$Q2[1]), NA, out_list$Q2["U"]),
    p_Q2  = ifelse(is.na(out_list$Q2[1]), NA, out_list$Q2["p"]),
    
    OR_Q3 = ifelse(is.na(out_list$Q3[1]), NA, out_list$Q3["OR"]),
    L_Q3  = ifelse(is.na(out_list$Q3[1]), NA, out_list$Q3["L"]),
    U_Q3  = ifelse(is.na(out_list$Q3[1]), NA, out_list$Q3["U"]),
    p_Q3  = ifelse(is.na(out_list$Q3[1]), NA, out_list$Q3["p"]),
    
    OR_Q4 = ifelse(is.na(out_list$Q4[1]), NA, out_list$Q4["OR"]),
    L_Q4  = ifelse(is.na(out_list$Q4[1]), NA, out_list$Q4["L"]),
    U_Q4  = ifelse(is.na(out_list$Q4[1]), NA, out_list$Q4["U"]),
    p_Q4  = ifelse(is.na(out_list$Q4[1]), NA, out_list$Q4["p"]),
    
    p_trend = p_trend,
    n_used = m_cat$n,
    events = m_cat$nevent,
    stringsAsFactors = FALSE
  )
}

# 用 clogit_quartile_unadj 跑 exposure_SG
results_PP_clogit_quartile <- do.call(
  rbind,
  lapply(exposure_SG, clogit_quartile_unadj, data = dat_ana)
)

results_PP_clogit_quartile$FDR_Q2 <- p.adjust(results_PP_clogit_quartile$p_Q2, method = "BH")
results_PP_clogit_quartile$FDR_Q3 <- p.adjust(results_PP_clogit_quartile$p_Q3, method = "BH")
results_PP_clogit_quartile$FDR_Q4 <- p.adjust(results_PP_clogit_quartile$p_Q4, method = "BH")

fmt <- function(or, lo, hi) ifelse(is.na(or), NA, sprintf("%.2f (%.2f, %.2f)", or, lo, hi))
results_PP_clogit_quartile$Q2 <- fmt(results_PP_clogit_quartile$OR_Q2, results_PP_clogit_quartile$L_Q2, results_PP_clogit_quartile$U_Q2)
results_PP_clogit_quartile$Q3 <- fmt(results_PP_clogit_quartile$OR_Q3, results_PP_clogit_quartile$L_Q3, results_PP_clogit_quartile$U_Q3)
results_PP_clogit_quartile$Q4 <- fmt(results_PP_clogit_quartile$OR_Q4, results_PP_clogit_quartile$L_Q4, results_PP_clogit_quartile$U_Q4)

results_PP_clogit_quartile <- results_PP_clogit_quartile[, c("exposure",
                                                             "Q2","p_Q2","FDR_Q2",
                                                             "Q3","p_Q3","FDR_Q3",
                                                             "Q4","p_Q4","FDR_Q4",
                                                             "p_trend",
                                                             "n_used","events")]

write.csv(results_PP_clogit_quartile, "results_PP_clogit_quartile_unadj.csv", row.names = FALSE)


