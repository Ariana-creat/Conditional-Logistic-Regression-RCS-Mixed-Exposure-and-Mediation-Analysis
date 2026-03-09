# =============================================================================
# BKMR 双变量交互作用分析（基于 dat_ana）
# 依赖：dat_ana 及 exposure_SG, outcome, covars 已存在；或见下方“变量设定”处定义
# =============================================================================

# -----------------------------------------------------------------------------
# 0. 包与变量设定（若尚未定义则在此定义）
# -----------------------------------------------------------------------------
library(aBKMR)
library(bkmr)
library(ggplot2)
library(dplyr)
library(glue)
library(ggtext)

# 若前面未定义，可在此指定（需与 dat_ana 列名一致）
if (!exists("exposure_SG")) {
  exposure_SG <- c("lnSCN_adj", "lnIodie_adj", "lnCLO3_adj", "lnCLO4_adj", "lnTCAA_adj")
}
if (!exists("outcome")) outcome <- "PP"
if (!exists("covars")) {
  covars <- c("age_c", "Parity_c", "Delive_mode", "GWG", "tobacco", "predclass", "MENSESAGE_c")
}

# 分析用数据：从 dat_ana 取所需变量
d <- dat_ana
d$PP <- ifelse(d$ICP == 1 | d$CPP == 1, 1L, 0L)
for (v in exposure_SG) d[[v]][!is.finite(d[[v]])] <- NA

vars_need <- unique(c(outcome, exposure_SG, covars))
d_bkmr <- d[, vars_need, drop = FALSE]

# 分类变量转为 factor（水平顺序按你数据定）
d_bkmr$age_c       <- factor(d_bkmr$age_c)
d_bkmr$Parity_c    <- factor(d_bkmr$Parity_c)
d_bkmr$Delive_mode <- factor(d_bkmr$Delive_mode)
d_bkmr$GWG         <- factor(d_bkmr$GWG)
d_bkmr$tobacco     <- factor(d_bkmr$tobacco)
d_bkmr$predclass   <- factor(d_bkmr$predclass)
d_bkmr$MENSESAGE_c <- factor(d_bkmr$MENSESAGE_c)

d_bkmr <- d_bkmr[complete.cases(d_bkmr), ]
y <- d_bkmr[[outcome]]
Z <- as.matrix(d_bkmr[, exposure_SG])
X <- model.matrix(
  ~ age_c + Parity_c + Delive_mode + GWG + tobacco + predclass + MENSESAGE_c,
  data = d_bkmr
)[, -1, drop = FALSE]

# -----------------------------------------------------------------------------
# 1. BKMR 拟合（aBKMR：knots + a_kmbayes）
# -----------------------------------------------------------------------------
set.seed(20260126)
k_id <- sam_py_r(R = Z, nd = 100, num_nn = 100, w = FALSE)
km <- a_kmbayes(
  y = y, Z = Z, X = X,
  knots = Z[k_id, ],
  est.h = TRUE,
  varsel = TRUE,
  iter = 10000,
  family = "binomial"
)

# -----------------------------------------------------------------------------
# 2. 双变量交互：a_PredictorResponseBivar + PredictorResponseBivarLevels
# -----------------------------------------------------------------------------
inter <- a_PredictorResponseBivar(
  fit = km,
  data.comps = km$data.comps,
  min.plot.dist = 2 #修改了这个参数），使曲线相连
)
inter2 <- PredictorResponseBivarLevels(inter, Z = Z, qs = c(0.25, 0.75))

# 统一估计/标准误列名（aBKMR 可能输出 est_p/sd_p）
if (!"est" %in% names(inter2) && "est_p" %in% names(inter2)) inter2$est <- inter2$est_p
if (!"se" %in% names(inter2) && "sd_p" %in% names(inter2))  inter2$se  <- inter2$sd_p

# 确保 quantile 为因子且有两水平
inter2$quantile <- factor(as.character(inter2$quantile), levels = c("0.25", "0.75"))

# -----------------------------------------------------------------------------
# 3. 交互作用 P 值：按 (variable1, variable2) 逐对计算
# -----------------------------------------------------------------------------
com_2interaction <- function(x1, x2_group, df, est_col = "est", se_col = "se") {
  stopifnot(all(c(x1, x2_group, est_col, se_col) %in% names(df)))
  d <- df
  d[[x1]]       <- as.numeric(d[[x1]])
  d[[est_col]]  <- as.numeric(d[[est_col]])
  d[[se_col]]   <- as.numeric(d[[se_col]])
  d[[x2_group]] <- as.factor(d[[x2_group]])
  d <- d[is.finite(d[[x1]]) & is.finite(d[[est_col]]) & is.finite(d[[se_col]]) &
         d[[se_col]] > 0 & !is.na(d[[x2_group]]), , drop = FALSE]
  if (nrow(d) < 10L || nlevels(d[[x2_group]]) < 2L) {
    return(data.frame(term = NA_character_, P = NA_real_, stringsAsFactors = FALSE))
  }
  w <- 1 / (d[[se_col]]^2)
  fml <- as.formula(sprintf("%s ~ %s * %s", est_col, x1, x2_group))
  fit_lm <- lm(fml, data = d, weights = w)
  sm <- summary(fit_lm)$coefficients
  inter_rows <- grep(paste0("^", x1, ":"), rownames(sm))
  if (length(inter_rows) == 0L) inter_rows <- grep(paste0(":", x2_group), rownames(sm))
  if (length(inter_rows) == 0L) {
    return(data.frame(term = NA_character_, P = NA_real_, stringsAsFactors = FALSE))
  }
  data.frame(
    term = rownames(sm)[inter_rows],
    P    = sm[inter_rows, "Pr(>|t|)"],
    stringsAsFactors = FALSE
  )
}

# 按 (variable1, variable2) 分组，每组算一个交互 P（不依赖 tidyr）
pairs <- inter2 %>% distinct(variable1, variable2)
inter2_res_list <- vector("list", nrow(pairs))
for (i in seq_len(nrow(pairs))) {
  v1 <- pairs$variable1[i]
  v2 <- pairs$variable2[i]
  sub <- inter2 %>% filter(variable1 == v1, variable2 == v2)
  if (!"est" %in% names(sub) && "est_p" %in% names(sub)) sub$est <- sub$est_p
  if (!"se" %in% names(sub) && "sd_p" %in% names(sub))  sub$se  <- sub$sd_p
  if (!all(c("est", "se") %in% names(sub))) {
    inter2_res_list[[i]] <- data.frame(variable1 = v1, variable2 = v2, P_interaction = NA_real_, stringsAsFactors = FALSE)
  } else {
    res <- com_2interaction("z1", "quantile", sub, est_col = "est", se_col = "se")
    p_val <- if (all(is.na(res$P))) NA_real_ else min(res$P, na.rm = TRUE)
    inter2_res_list[[i]] <- data.frame(variable1 = v1, variable2 = v2, P_interaction = p_val, stringsAsFactors = FALSE)
  }
}
inter2_res_for_plot <- bind_rows(inter2_res_list)
inter2_res_for_plot$P_interaction[!is.finite(inter2_res_for_plot$P_interaction)] <- NA_real_

# -----------------------------------------------------------------------------
# 4. 双变量交互图 + 每格标注 P
# -----------------------------------------------------------------------------
# 用于 geom_richtext 的标注：每个 (variable1, variable2) 一格，需含 variable1/variable2 以便分面
annotations <- inter2_res_for_plot %>%
  mutate(
    x = Inf,
    y = Inf,
    label = ifelse(is.na(P_interaction), "P (NA)",
                   sprintf("<i>P</i><sub>interaction</sub> = %s",
                           ifelse(P_interaction < 0.001, "< 0.001", sprintf("%.3f", P_interaction))))
  )

# 画图：暴露名用 variable1 / variable2（即 Z 的列名）
p_bivar <- ggplot(
  inter2,
  aes(x = z1, y = est, colour = quantile)
) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se, fill = quantile),
              alpha = 0.15, colour = NA) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray50") +
  facet_grid(variable2 ~ variable1, scales = "free_x") +
  labs(
    title = "BKMR 双变量交互：h(暴露1 | 暴露2 分位数)",
    x = "暴露1 (分位数/原始尺度)",
    y = "相对风险估计",
    colour = "暴露2 分位数",
    fill  = "暴露2 分位数"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# 若存在 geom_richtext 的 annotations，可加上（需 ggtext）
if (nrow(annotations) > 0L && "label" %in% names(annotations)) {
  p_bivar <- p_bivar +
    geom_richtext(
      data = annotations,
      aes(x = x, y = y, label = label),
      inherit.aes = FALSE,
      fill = NA, label.color = NA, color = "red", size = 3.5,
      hjust = 1.05, vjust = 1.2
    )
}

print(p_bivar)

# 保存
dir.create("BKMR_results", showWarnings = FALSE)
ggsave("BKMR_results/BKMR_bivar_interaction.png", p_bivar, width = 12, height = 10, dpi = 300)
ggsave("BKMR_results/BKMR_bivar_interaction.pdf", p_bivar, width = 12, height = 10)

# 交互 P 值表
write.csv(inter2_res_for_plot, "BKMR_results/BKMR_bivar_interaction_P.csv", row.names = FALSE)
inter2_res_for_plot
