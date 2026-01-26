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
#将infinite转换成NA，否则模型报错
for (v in exposure_SG) {
     dat_ana[[v]][!is.finite(dat_ana[[v]])] <- NA
}
#构造clogit—_1函数
clogit_1 <- function(exposure, data,
                           outcome = "PP",
                           covars = covar,
                           strata_var = "Subclass") {
    fomula <- as.formula(paste0(
    outcome, " ~ ", exposure, " + ",
    paste(covars, collapse = " + "),
    " + strata(", strata_var, ")"
  ))
  #拟合clogit
  m <- clogit(fomula, data = data)
  #提取exposure的coef(β)，自变量为自然对数形式，按 OR_×2 = exp(β·ln2) = 2^β 转为“每翻倍”OR
  b  <- coef(m)[exposure]
  se <- sqrt(vcov(m)[exposure, exposure])
  OR <- 2^b
  lower <- 2^(b - 1.96 * se)
  upper <- 2^(b + 1.96 * se)
  p  <- summary(m)$coefficients[exposure, "Pr(>|z|)"]
  #返回结果（OR_×2 为原始暴露每翻倍时的OR）
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
#循环跑模型
results_PP_clogit <- do.call(rbind, lapply(exposure_SG, clogit_1, data = dat_ana))
#FDR校正(Benjamini–Hochberg)
results_PP_clogit$p_FDR <- p.adjust(results_PP_clogit$p_value, method = "BH")

results_PP_clogit$OR_CI <- sprintf("%.2f (%.2f, %.2f)",
                                 results_PP_clogit$OR_doubling,
                                 results_PP_clogit$CI_low,
                                 results_PP_clogit$CI_high)
results_PP_clogit <- results_PP_clogit[order(results_PP_clogit$p_value), ]
results_PP_clogit[, c("exposure","OR_CI","p_value","p_FDR","n_used","events")]

#输出结果（OR 为暴露每翻倍时的 OR_×2 = 2^β）
write.csv(results_PP_clogit, "results_PP_clogit.csv", row.names = FALSE)

#（4）条件logistic——四分位
library(survival)

clogit_quartile <- function(exposure, data,
                            outcome = "PP",
                            covars = covar,
                            strata_var = "Subclass") {
  
  d <- data
  #报错后改
  stopifnot(is.character(exposure), length(exposure) == 1)
  x_all <- suppressWarnings(as.numeric(d[[exposure]]))
  idx_ctrl <- d[[outcome]] == 0 & is.finite(x_all)
  x_ctrl <- x_all[idx_ctrl]
  
  #四分位切点
  qs <- as.numeric(quantile(x_ctrl, probs = c(0, .25, .5, .75, 1), na.rm = TRUE, type = 7))
  
  #报错后增
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
  
  #分四分位
  d$Q <- cut(x_all, breaks = qs, include.lowest = TRUE, labels = c("Q1","Q2","Q3","Q4"))
  d$Q <- relevel(d$Q, ref = "Q1")
  
  #趋势变量
  idx_ctrl2 <- d[[outcome]] == 0 & !is.na(d$Q) & is.finite(x_all)
  med_ctrl <- sapply(split(x_all[idx_ctrl2], d$Q[idx_ctrl2]),
                     function(z) median(z, na.rm = TRUE))
  d$Q_med <- med_ctrl[as.character(d$Q)]
  
  #分类
  fml_cat <- as.formula(paste0(
    outcome, " ~ Q + ", paste(covars, collapse = " + "),
    " + strata(", strata_var, ")"
  ))
  m_cat <- clogit(fml_cat, data = d)
  
  #趋势检验
  fml_trend <- as.formula(paste0(
    outcome, " ~ Q_med + ", paste(covars, collapse = " + "),
    " + strata(", strata_var, ")"
  ))
  m_trend <- clogit(fml_trend, data = d)
  
  #提取summary
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
  #整理结果
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

#用clogit_quartile跑5个污染物
results_PP_clogit_quartile <- do.call(rbind, lapply(exposure_SG, clogit_quartile, data = dat_ana))
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
write.csv(results_PP_clogit_quartile, "results_PP_clogit_quartile.csv", row.names = FALSE)

#（4*）RCS
library(rms)
#（5%、27.5%、72.5%、95%）
get_knots <- function(x) {
  quantile(x, probs = c(0.05, 0.275, 0.725, 0.95), na.rm = TRUE)
}
#RCS
run_rcs <- function(exposure, data,
                    outcome = "PP",
                    covars = c("age_c","Parity_c","Delive_mode","GWG","tobacco","predclass","MENSESAGE_c"),
                    strata_var = "Subclass") {
  
  d <- data
  
  #确保暴露是 numeric 且有限
  d[[exposure]] <- suppressWarnings(as.numeric(d[[exposure]]))
  d[[exposure]][!is.finite(d[[exposure]])] <- NA
  
  #knots
  k <- get_knots(d[[exposure]])
  
  f_lin <- as.formula(paste0(
    outcome, " ~ ", exposure, " + ",
    paste(covars, collapse = " + "),
    " + strata(", strata_var, ")"
  ))
  
  f_rcs <- as.formula(paste0(
    outcome, " ~ rcs(", exposure, ", c(",
    paste(sprintf("%.10f", k), collapse = ","),
    ")) + ",
    paste(covars, collapse = " + "),
    " + strata(", strata_var, ")"
  ))
  
  mod_lin <- clogit(f_lin, data = d)
  mod_rcs <- clogit(f_rcs, data = d)
  
  p_linear <- summary(mod_lin)$coefficients[exposure, "Pr(>|z|)"]
  
  #报错后修改
  lrt <- tryCatch(anova(mod_lin, mod_rcs), error = function(e) NULL)
  if (is.null(lrt) || nrow(lrt) < 2) {
    p_nonlinear <- NA_real_
    chisq <- NA_real_
    df <- NA_real_
  } else {
    # anova 输出可能是 matrix / data.frame
    lrt <- as.data.frame(lrt)
    
    # 安全取值：列不存在就返回 NA，不会产生 numeric(0)
    safe_get <- function(tab, row, col) {
      if (!is.null(tab) && nrow(tab) >= row && col %in% colnames(tab)) {
        val <- tab[row, col]
        return(as.numeric(val))
      }
      NA_real_
    }
    
    # 兼容不同 survival 版本的 p 列名
    p_col_candidates <- c("Pr(>Chisq)", "P(>|Chi|)", "Pr(>|Chi|)", "p")
    p_col <- p_col_candidates[p_col_candidates %in% colnames(lrt)][1]
    
    p_nonlinear <- if (!is.na(p_col)) safe_get(lrt, 2, p_col) else NA_real_
    chisq <- safe_get(lrt, 2, "Chisq")
    df    <- safe_get(lrt, 2, "Df")
  }
  
  
  data.frame(
    exposure = exposure,
    p_linear = as.numeric(p_linear),
    p_nonlinear = as.numeric(p_nonlinear),
    chisq = as.numeric(chisq),
    df = as.numeric(df),
    n_used = mod_rcs$n,
    events = mod_rcs$nevent,
    stringsAsFactors = FALSE
  )
}

rcs_p_table <- do.call(rbind, lapply(exposure_SG, run_rcs, data = dat_ana))

# 导出
write.csv(rcs_p_table, "rcs_p_table.csv", row.names = FALSE)

rcs_p_table

#（4**）rcs图
library(Cairo)
library(rcssci)
ls("package:rcssci")

rcssci_logistic(
  data = dat_ana,
  y = "PP",
  x = "lnSCN_adj",
  covs = c("age_c","Parity_c","Delive_mode","GWG","tobacco","predclass","MENSESAGE_c"),
  knots = 4,
  prob = 0.5,
  filepath = "/Users/lixinrui/Desktop/混合暴露分析/Mixture Exposure"
)

