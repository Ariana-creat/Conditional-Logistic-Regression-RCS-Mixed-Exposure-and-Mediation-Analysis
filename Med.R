# =============================================================================
# 中介效应分析 (Mediation Analysis)
# 变量定义依据: Note.docx
# X: 暴露 | M: 中介 | Y: 结局 (PP/ICP/CPP)
# 若校正 SG：暴露/中介用原始浓度对数(ln*_sub)，协变量含 SG
# 若不校正 SG：暴露/中介用 SG 校正后对数(ln*_adj)，协变量不含 SG
# =============================================================================

#（1）分析数据集
dat_ana <- get("ion_cpp_analy")
dim(dat_ana)

#（2）变量定义（与 Note.docx 一致）
# --- 暴露 Exposure (X) ---
# 原始浓度对数转换
exposure_sub <- c("lnSCN_sub", "lnIodie_sub", "lnCLO3_sub", "lnCLO4_sub", "lnTCAA_sub")
# SG 校正后对数转换
exposure_adj <- c("lnSCN_adj", "lnIodie_adj", "lnCLO3_adj", "lnCLO4_adj", "lnTCAA_adj")

# --- 结局 Outcome (Y) ---
outcomes <- c("PP", "ICP", "CPP")

# --- 协变量 Covariates ---
# 不校正 SG 时使用
covar <- c("age_c", "Parity_c", "Delive_mode", "GWG", "tobacco", "predclass", "MENSESAGE_c")
# 校正 SG 时使用（在 covar 基础上加入 SG）
covar_sg <- c(covar, "SG")

# --- 中介变量 Mediators (M) ---
# 原始浓度对数转换
mediators_sub <- c("lnOHdG_sub", "lnOHG_sub", "lnHNEMA_sub", "lnAlla_sub", "lnBY_sub", "lndiY_sub",
                   "lnE1_sub", "lnE2_sub", "lnE3_sub", "lnP4_sub", "lnOH17P4_sub", "lnOH21P4_sub",
                   "lnTesto_sub", "lnA4_sub", "lnDHT_sub", "lnCOH_sub", "lnCTONE_sub",
                   "lnCortisoF_sub", "lnDeoxyCortisolS_sub", "lnALDO_sub"
)
# SG 校正后对数转换
mediators_adj <- c(
  "lnOHdG_adj", "lnOHG_adj", "lnHNEMA_adj", "lnAlla_adj", "lnBY_adj", "lndiY_adj",
  "lnE1_adj", "lnE2_adj", "lnE3_adj", "lnP4_adj", "lnOH17P4_adj", "lnOH21P4_adj",
  "lnTesto_adj", "lnA4_adj", "lnDHT_adj", "lnCOH_adj", "lnCTONE_adj",
  "lnCortisoF_adj", "lnDeoxyCortisolS_adj", "lnALDO_adj"
)

# 选择分析模式：TRUE = 校正 SG（用 _sub + covar_sg），FALSE = 不校正 SG（用 _adj + covar）
adjust_sg <- FALSE
if (adjust_sg) {
  exposure_set <- exposure_sub
  mediators_set <- mediators_sub
  covar_use <- covar_sg
} else {
  exposure_set <- exposure_adj
  mediators_set <- mediators_adj
  covar_use <- covar
}

# 数据准备
dat_ana$PP <- ifelse(dat_ana$ICP == 1L | dat_ana$CPP == 1L, 1L, 0L)
dat_ana$Delive_mode <- factor(dat_ana$Delive_mode)
dat_ana$tobacco     <- factor(dat_ana$tobacco)

# 暴露与中介中的无穷值置为 NA
all_x_m <- c(exposure_set, mediators_set)
for (v in all_x_m) {
  if (v %in% names(dat_ana))
    dat_ana[[v]][!is.finite(dat_ana[[v]])] <- NA
}

# 检查变量是否都在数据中
setdiff(c(exposure_set, outcomes, covar_use, mediators_set), names(dat_ana))

# =============================================================================
# 性能参数（可按机器性能与精度需求调整）
# =============================================================================
max_cores_use <- 8L
detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
if (!is.finite(detected_cores) || detected_cores < 1L) detected_cores <- 1L
n_cores_use <- min(max_cores_use, max(1L, as.integer(detected_cores) - 1L))
n_boot_main <- 500L            # 原脚本为 2000；可改回 2000 获取更高精度但更慢
boot_screen_mode <- TRUE       # TRUE: 先用 delta 近似筛选，再对候选组合做 bootstrap
boot_p_screen <- 0.20          # delta p_ind <= 该阈值时才做 bootstrap
run_med_analysis <- TRUE       # TRUE 运行传统中介；FALSE 跳过（仅跑 SEM）
run_sem_analysis <- TRUE       # SEM 耗时较高；仅在需要时改为 TRUE
match_strata_var <- "Subclass" # 匹配层变量名
path_sig_alpha <- 0.05         # 路径显著性阈值：用于筛选 X-M、M-Y、X-Y 均显著的组合

# 预筛选：先筛总效应稳定的 (outcome, exposure) 组合，再筛 mediator，减少多重比较
apply_pair_total_screen <- TRUE
screen_total_p_cut <- 0.05
screen_total_min_abs_beta <- 0.08
screen_total_min_abs_z <- 1.96
screen_total_min_n <- 200L
screen_total_min_events <- 40L

apply_mediator_ab_screen <- TRUE
screen_mediator_p_a_cut <- 0.10
screen_mediator_p_b_cut <- 0.10
screen_mediator_top_n <- 5L
screen_mediator_min_keep <- 2L

# SEM 端使用匹配层：在 lavaan 中以 cluster-robust 方式纳入 Subclass
sem_use_subclass_cluster <- TRUE
out_dir <- file.path("..", "output_csv", "Med_results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 中介分析函数（条件 logistic：Y 为二分类，匹配设计 strata = Subclass）
# 路径 a: X -> M (线性回归)
# 路径 c: X -> Y (总效应, clogit)
# 路径 c', b: X + M -> Y (直接效应 c', 中介效应路径 b, clogit)
# 间接效应 ≈ a*b；95% CI 采用 Bootstrap（对匹配层 Subclass 重抽样）
# =============================================================================
library(survival)

join_terms <- function(terms) {
  terms <- terms[!is.na(terms) & nzchar(terms)]
  paste(terms, collapse = " + ")
}

classify_mediation_type <- function(direct_effect, indirect_effect, p_direct, p_indirect, alpha = 0.05) {
  sig_direct <- is.finite(p_direct) && (p_direct < alpha)
  sig_indirect <- is.finite(p_indirect) && (p_indirect < alpha)
  sign_prod <- if (is.finite(direct_effect) && is.finite(indirect_effect)) direct_effect * indirect_effect else NA_real_
  
  if (sig_direct && sig_indirect) {
    if (is.finite(sign_prod) && sign_prod < 0) return("competitive_mediation_suppression")
    if (is.finite(sign_prod) && sign_prod > 0) return("complementary_mediation")
    return("both_significant_direction_unclear")
  }
  if (!sig_direct && sig_indirect) {
    if (is.finite(sign_prod) && sign_prod < 0) return("indirect_only_possible_suppression")
    if (is.finite(sign_prod) && sign_prod > 0) return("indirect_only_mediation")
    return("indirect_only_direction_unclear")
  }
  if (sig_direct && !sig_indirect) return("direct_only_no_mediation")
  "no_evidence_of_mediation"
}

compute_prop_metrics <- function(total_effect, direct_effect, p_total, alpha = 0.05) {
  if (!is.finite(total_effect) || abs(total_effect) <= 1e-8 || !is.finite(direct_effect)) {
    return(list(prop_raw = NA_real_, prop_report = NA_real_, prop_interpretable = FALSE, prop_note = "undefined_total"))
  }
  ind_from_diff <- total_effect - direct_effect
  prop_raw <- ind_from_diff / total_effect
  same_direction <- (total_effect * ind_from_diff) > 0
  total_sig <- is.finite(p_total) && (p_total < alpha)
  interpretable <- isTRUE(total_sig) && isTRUE(same_direction)
  prop_note <- if (!total_sig) {
    "total_not_significant"
  } else if (!same_direction) {
    "suppression_or_inconsistent_direction"
  } else {
    "ok_for_reporting"
  }
  list(
    prop_raw = prop_raw,
    prop_report = if (interpretable) prop_raw else NA_real_,
    prop_interpretable = interpretable,
    prop_note = prop_note
  )
}

safe_extract_term <- function(fit, term, p_col_pref = c("Pr(>|z|)", "Pr(>|t|)")) {
  if (is.null(fit)) return(c(beta = NA_real_, se = NA_real_, p = NA_real_))
  cf <- tryCatch(coef(fit), error = function(e) NULL)
  if (is.null(cf) || !(term %in% names(cf))) return(c(beta = NA_real_, se = NA_real_, p = NA_real_))
  beta <- as.numeric(cf[[term]])
  
  se <- NA_real_
  vcv <- tryCatch(vcov(fit), error = function(e) NULL)
  if (!is.null(vcv) && term %in% rownames(vcv) && term %in% colnames(vcv)) {
    se <- suppressWarnings(as.numeric(sqrt(vcv[term, term])))
  }
  
  p <- NA_real_
  sm <- tryCatch(coef(summary(fit)), error = function(e) NULL)
  if (!is.null(sm) && term %in% rownames(sm)) {
    p_col <- intersect(p_col_pref, colnames(sm))
    if (length(p_col) >= 1L) p <- suppressWarnings(as.numeric(sm[term, p_col[1L]]))
  }
  
  c(beta = beta, se = se, p = p)
}

fit_total_effect_clogit <- function(data, x_var, y_var, covars = covar_use, strata_var = match_strata_var) {
  vars_needed <- unique(c(x_var, y_var, covars, strata_var))
  vars_needed <- vars_needed[vars_needed %in% names(data)]
  if (!all(c(x_var, y_var, strata_var) %in% vars_needed)) {
    return(data.frame(
      outcome = y_var, exposure = x_var,
      coef_total = NA_real_, se_total = NA_real_, p_total = NA_real_, z_total = NA_real_,
      ci_total_low = NA_real_, ci_total_high = NA_real_,
      n_used = 0L, events = 0L, stringsAsFactors = FALSE
    ))
  }
  d <- data[complete.cases(data[, vars_needed, drop = FALSE]), vars_needed, drop = FALSE]
  n_used <- nrow(d)
  events <- if (y_var %in% names(d)) sum(d[[y_var]] == 1L, na.rm = TRUE) else 0L
  if (n_used < 30L || length(unique(d[[y_var]])) < 2L) {
    return(data.frame(
      outcome = y_var, exposure = x_var,
      coef_total = NA_real_, se_total = NA_real_, p_total = NA_real_, z_total = NA_real_,
      ci_total_low = NA_real_, ci_total_high = NA_real_,
      n_used = n_used, events = events, stringsAsFactors = FALSE
    ))
  }
  
  covars_use <- covars[covars %in% names(d)]
  rhs <- join_terms(c(x_var, covars_use))
  fml <- as.formula(paste0(y_var, " ~ ", rhs, " + strata(", strata_var, ")"))
  fit <- tryCatch(clogit(fml, data = d), error = function(e) NULL)
  est <- safe_extract_term(fit, x_var, p_col_pref = "Pr(>|z|)")
  beta <- unname(est["beta"])
  se <- unname(est["se"])
  p <- unname(est["p"])
  z <- if (is.finite(beta) && is.finite(se) && se > 0) beta / se else NA_real_
  ci_low <- if (is.finite(beta) && is.finite(se)) beta - 1.96 * se else NA_real_
  ci_high <- if (is.finite(beta) && is.finite(se)) beta + 1.96 * se else NA_real_
  
  data.frame(
    outcome = y_var, exposure = x_var,
    coef_total = beta, se_total = se, p_total = p, z_total = z,
    ci_total_low = ci_low, ci_total_high = ci_high,
    n_used = n_used, events = events, stringsAsFactors = FALSE
  )
}

screen_total_pairs <- function(data, outcomes, exposures, covars = covar_use, strata_var = match_strata_var,
                               p_cut = screen_total_p_cut,
                               min_abs_beta = screen_total_min_abs_beta,
                               min_abs_z = screen_total_min_abs_z,
                               min_n = screen_total_min_n,
                               min_events = screen_total_min_events) {
  grid <- expand.grid(outcome = outcomes, exposure = exposures, stringsAsFactors = FALSE)
  out <- lapply(seq_len(nrow(grid)), function(i) {
    fit_total_effect_clogit(
      data = data,
      x_var = grid$exposure[i],
      y_var = grid$outcome[i],
      covars = covars,
      strata_var = strata_var
    )
  })
  stats <- do.call(rbind, out)
  stats$pass_total <- with(stats,
                           is.finite(p_total) &
                             is.finite(coef_total) &
                             is.finite(z_total) &
                             (p_total <= p_cut) &
                             (abs(coef_total) >= min_abs_beta) &
                             (abs(z_total) >= min_abs_z) &
                             (n_used >= min_n) &
                             (events >= min_events)
  )
  stats$pass_total <- as.logical(stats$pass_total)
  stats$screen_reason <- ifelse(
    stats$pass_total, "pass",
    ifelse(!is.finite(stats$p_total), "model_or_term_unavailable",
           ifelse(stats$n_used < min_n, "n_below_threshold",
                  ifelse(stats$events < min_events, "events_below_threshold",
                         ifelse(abs(stats$coef_total) < min_abs_beta, "abs_beta_below_threshold",
                                ifelse(abs(stats$z_total) < min_abs_z, "abs_z_below_threshold", "p_total_above_threshold")
                         )
                  )
           )
    )
  )
  stats
}

screen_mediators_for_pair <- function(data, x_var, y_var, mediators, covars = covar_use, strata_var = match_strata_var,
                                      p_a_cut = screen_mediator_p_a_cut,
                                      p_b_cut = screen_mediator_p_b_cut,
                                      top_n = screen_mediator_top_n,
                                      min_keep = screen_mediator_min_keep) {
  meds <- mediators[mediators %in% names(data)]
  if (length(meds) == 0L) {
    return(data.frame(
      outcome = character(0), exposure = character(0), mediator = character(0),
      p_a = numeric(0), p_b = numeric(0), score = numeric(0),
      pass_ab = logical(0), selected = logical(0), selected_by = character(0),
      stringsAsFactors = FALSE
    ))
  }
  
  stats_list <- lapply(meds, function(m_var) {
    vars_needed <- unique(c(x_var, m_var, y_var, covars, strata_var))
    vars_needed <- vars_needed[vars_needed %in% names(data)]
    d <- data[complete.cases(data[, vars_needed, drop = FALSE]), vars_needed, drop = FALSE]
    if (nrow(d) < 30L || length(unique(d[[y_var]])) < 2L) {
      return(data.frame(
        outcome = y_var, exposure = x_var, mediator = m_var,
        p_a = NA_real_, p_b = NA_real_, score = -Inf,
        pass_ab = FALSE, selected = FALSE, selected_by = "insufficient_data",
        stringsAsFactors = FALSE
      ))
    }
    covars_use <- covars[covars %in% names(d)]
    fit_a <- tryCatch(lm(as.formula(paste0(m_var, " ~ ", join_terms(c(x_var, covars_use)))), data = d), error = function(e) NULL)
    fit_b <- tryCatch(
      clogit(as.formula(paste0(y_var, " ~ ", join_terms(c(x_var, m_var, covars_use)), " + strata(", strata_var, ")")), data = d),
      error = function(e) NULL
    )
    est_a <- safe_extract_term(fit_a, x_var, p_col_pref = "Pr(>|t|)")
    est_b <- safe_extract_term(fit_b, m_var, p_col_pref = "Pr(>|z|)")
    p_a <- unname(est_a["p"])
    p_b <- unname(est_b["p"])
    score <- if (is.finite(p_a) && is.finite(p_b)) -log10(max(p_a, 1e-16)) - log10(max(p_b, 1e-16)) else -Inf
    pass_ab <- is.finite(p_a) && is.finite(p_b) && p_a <= p_a_cut && p_b <= p_b_cut
    data.frame(
      outcome = y_var, exposure = x_var, mediator = m_var,
      p_a = p_a, p_b = p_b, score = score,
      pass_ab = pass_ab, selected = FALSE, selected_by = "not_selected",
      stringsAsFactors = FALSE
    )
  })
  stats <- do.call(rbind, stats_list)
  stats <- stats[order(-stats$score), , drop = FALSE]
  pass_idx <- which(stats$pass_ab)
  n_target <- max(as.integer(min_keep), min(as.integer(top_n), nrow(stats)))
  if (length(pass_idx) >= as.integer(min_keep)) {
    keep_idx <- pass_idx[seq_len(min(length(pass_idx), as.integer(top_n)))]
    stats$selected[keep_idx] <- TRUE
    stats$selected_by[keep_idx] <- "pass_ab_screen"
  } else if (nrow(stats) > 0L) {
    keep_idx <- seq_len(n_target)
    stats$selected[keep_idx] <- TRUE
    stats$selected_by[keep_idx] <- "top_score_fallback"
  }
  stats
}

build_screened_combos <- function(data, outcomes, exposures, mediators, covars = covar_use, strata_var = match_strata_var,
                                  apply_pair_screen = apply_pair_total_screen,
                                  apply_mediator_screen = apply_mediator_ab_screen) {
  pair_stats <- screen_total_pairs(data, outcomes, exposures, covars = covars, strata_var = strata_var)
  pair_use <- if (isTRUE(apply_pair_screen)) pair_stats[pair_stats$pass_total, c("outcome", "exposure"), drop = FALSE] else pair_stats[, c("outcome", "exposure"), drop = FALSE]
  if (nrow(pair_use) == 0L) {
    warning("Pair total-effect screening removed all exposure-outcome pairs; fallback to full pair set.")
    pair_use <- unique(pair_stats[, c("outcome", "exposure"), drop = FALSE])
  }
  
  if (!isTRUE(apply_mediator_screen)) {
    combos <- merge(pair_use, data.frame(mediator = mediators, stringsAsFactors = FALSE), by = NULL)
    combos <- combos[combos$mediator %in% names(data), , drop = FALSE]
    mediator_stats <- data.frame(
      outcome = combos$outcome, exposure = combos$exposure, mediator = combos$mediator,
      p_a = NA_real_, p_b = NA_real_, score = NA_real_,
      pass_ab = NA, selected = TRUE, selected_by = "no_mediator_screen",
      stringsAsFactors = FALSE
    )
    return(list(combos = combos, pair_stats = pair_stats, mediator_stats = mediator_stats))
  }
  
  med_stats_list <- lapply(seq_len(nrow(pair_use)), function(i) {
    screen_mediators_for_pair(
      data = data,
      x_var = pair_use$exposure[i],
      y_var = pair_use$outcome[i],
      mediators = mediators,
      covars = covars,
      strata_var = strata_var
    )
  })
  mediator_stats <- do.call(rbind, med_stats_list)
  combos <- mediator_stats[which(!is.na(mediator_stats$selected) & mediator_stats$selected), c("outcome", "exposure", "mediator"), drop = FALSE]
  combos <- unique(combos)
  if (nrow(combos) == 0L) {
    warning("Mediator screening removed all triples; fallback to pair x all mediators.")
    combos <- merge(pair_use, data.frame(mediator = mediators, stringsAsFactors = FALSE), by = NULL)
    combos <- combos[combos$mediator %in% names(data), , drop = FALSE]
  }
  list(combos = combos, pair_stats = pair_stats, mediator_stats = mediator_stats)
}

run_one_mediation <- function(data, x_var, m_var, y_var,
                              covars = covar_use,
                              strata_var = match_strata_var,
                              n_boot = n_boot_main,
                              boot_screen_mode = boot_screen_mode,
                              boot_p_screen = boot_p_screen) {
  # 完整个案
  vars_needed <- c(x_var, m_var, y_var, covars, strata_var)
  vars_needed <- vars_needed[vars_needed %in% names(data)]
  ok <- complete.cases(data[, vars_needed, drop = FALSE])
  d <- data[ok, ]
  if (nrow(d) < 30L || length(unique(d[[y_var]])) < 2L) return(NULL)
  covars <- covars[covars %in% names(d)]
  
  # 路径 a: M ~ X + covars
  rhs_a <- join_terms(c(x_var, covars))
  rhs_c <- join_terms(c(x_var, covars))
  rhs_cb <- join_terms(c(x_var, m_var, covars))
  fml_a <- as.formula(paste0(m_var, " ~ ", rhs_a))
  fit_a <- tryCatch(lm(fml_a, data = d), error = function(e) NULL)
  if (is.null(fit_a)) return(NULL)
  coef_a <- coef(fit_a)[x_var]
  se_a   <- sqrt(vcov(fit_a)[x_var, x_var])
  p_a    <- coef(summary(fit_a))[x_var, "Pr(>|t|)"]
  
  # 路径 c: Y ~ X + covars + strata
  fml_c <- as.formula(paste0(
    y_var, " ~ ", rhs_c,
    " + strata(", strata_var, ")"
  ))
  fit_c <- tryCatch(clogit(fml_c, data = d), error = function(e) NULL)
  if (is.null(fit_c)) return(NULL)
  coef_c <- coef(fit_c)[x_var]
  se_c   <- sqrt(vcov(fit_c)[x_var, x_var])
  p_c    <- coef(summary(fit_c))[x_var, "Pr(>|z|)"]
  
  # 路径 c' 与 b: Y ~ X + M + covars + strata
  fml_cb <- as.formula(paste0(
    y_var, " ~ ", rhs_cb,
    " + strata(", strata_var, ")"
  ))
  fit_cb <- tryCatch(clogit(fml_cb, data = d), error = function(e) NULL)
  if (is.null(fit_cb)) return(NULL)
  coef_cp <- coef(fit_cb)[x_var]
  coef_b  <- coef(fit_cb)[m_var]
  se_cp   <- sqrt(vcov(fit_cb)[x_var, x_var])
  se_b    <- sqrt(vcov(fit_cb)[m_var, m_var])
  p_cp    <- coef(summary(fit_cb))[x_var, "Pr(>|z|)"]
  p_b     <- coef(summary(fit_cb))[m_var, "Pr(>|z|)"]
  
  # 点估计：间接效应 = a * b
  ind_effect <- coef_a * coef_b
  se_ind <- sqrt(coef_a^2 * se_b^2 + coef_b^2 * se_a^2)
  p_ind  <- 2 * pnorm(-abs(ind_effect / se_ind))
  
  # 默认先用 delta 近似；仅对候选组合做 bootstrap（显著提速）
  ci_ind_low  <- ind_effect - 1.96 * se_ind
  ci_ind_high <- ind_effect + 1.96 * se_ind
  ci_ind_method <- "delta"
  ind_boot <- numeric(0)
  
  n_boot_target <- as.integer(max(0L, n_boot))
  if (isTRUE(boot_screen_mode) && (!is.finite(p_ind) || p_ind > boot_p_screen)) {
    n_boot_target <- 0L
  }
  
  if (n_boot_target > 0L) {
    boot_vars <- unique(c(x_var, m_var, y_var, covars, strata_var))
    d_boot <- d[, boot_vars, drop = FALSE]
    strata_idx_list <- split(seq_len(nrow(d_boot)), as.character(d_boot[[strata_var]]))
    n_strata <- length(strata_idx_list)
    if (n_strata >= 2L) {
      fml_cbb_boot <- as.formula(paste0(
        y_var, " ~ ", rhs_cb, " + strata(Subclass_boot)"
      ))
      ind_boot <- rep(NA_real_, n_boot_target)
      for (i in seq_len(n_boot_target)) {
        picked <- sample.int(n_strata, size = n_strata, replace = TRUE)
        picked_rows <- strata_idx_list[picked]
        rows <- unlist(picked_rows, use.names = FALSE)
        d_b <- d_boot[rows, , drop = FALSE]
        d_b$Subclass_boot <- rep.int(seq_along(picked_rows), lengths(picked_rows))
        
        fit_ab <- tryCatch(lm(fml_a, data = d_b), error = function(e) NULL)
        fit_cbb <- tryCatch(clogit(fml_cbb_boot, data = d_b), error = function(e) NULL)
        if (!is.null(fit_ab) && !is.null(fit_cbb)) {
          ind_boot[i] <- coef(fit_ab)[x_var] * coef(fit_cbb)[m_var]
        }
      }
      ind_boot <- ind_boot[is.finite(ind_boot)]
      if (length(ind_boot) >= 20L) {
        ci_ind_low  <- as.numeric(quantile(ind_boot, 0.025))
        ci_ind_high <- as.numeric(quantile(ind_boot, 0.975))
        ci_ind_method <- "bootstrap"
      }
    }
  }
  
  ci_total_low  <- coef_c - 1.96 * se_c
  ci_total_high <- coef_c + 1.96 * se_c
  ci_direct_low  <- coef_cp - 1.96 * se_cp
  ci_direct_high <- coef_cp + 1.96 * se_cp
  ind_from_diff <- coef_c - coef_cp
  prop_info <- compute_prop_metrics(coef_c, coef_cp, p_c)
  mediation_type <- classify_mediation_type(coef_cp, ind_effect, p_cp, p_ind)
  suppressive_role <- is.finite(coef_cp) && is.finite(ind_effect) && (coef_cp * ind_effect < 0)
  
  OR_total   <- 2^coef_c
  OR_direct  <- 2^coef_cp
  OR_ind_approx <- 2^ind_effect
  
  data.frame(
    exposure = x_var,
    mediator = m_var,
    outcome  = y_var,
    coef_a   = coef_a,
    se_a     = se_a,
    p_a      = p_a,
    coef_b   = coef_b,
    se_b     = se_b,
    p_b      = p_b,
    coef_c   = coef_c,
    se_c     = se_c,
    p_c      = p_c,
    ci_total_low  = ci_total_low,
    ci_total_high = ci_total_high,
    total_effect = coef_c,
    se_total = se_c,
    p_total  = p_c,
    coef_cp  = coef_cp,
    se_cp    = se_cp,
    p_cp     = p_cp,
    ci_direct_low  = ci_direct_low,
    ci_direct_high = ci_direct_high,
    ind_effect = ind_effect,
    se_ind   = se_ind,
    ci_ind_low  = ci_ind_low,
    ci_ind_high = ci_ind_high,
    ci_ind_method = ci_ind_method,
    p_ind    = p_ind,
    ind_from_diff = ind_from_diff,
    prop_med = prop_info$prop_raw,
    prop_med_raw = prop_info$prop_raw,
    prop_med_report = prop_info$prop_report,
    prop_med_interpretable = prop_info$prop_interpretable,
    prop_med_note = prop_info$prop_note,
    mediation_type = mediation_type,
    suppressive_role = suppressive_role,
    OR_total = OR_total,
    OR_direct = OR_direct,
    OR_ind_approx = OR_ind_approx,
    n_used   = nrow(d),
    events   = sum(d[[y_var]] == 1L, na.rm = TRUE),
    n_boot   = length(ind_boot),
    stringsAsFactors = FALSE
  )
}

# =============================================================================
# 运行中介分析：暴露 × 中介 × 结局（多核并行，充分利用 CPU）
# =============================================================================
combos_shared <- NULL
if (isTRUE(run_med_analysis) || isTRUE(run_sem_analysis)) {
  combo_pack <- build_screened_combos(
    data = dat_ana,
    outcomes = outcomes,
    exposures = exposure_set,
    mediators = mediators_set,
    covars = covar_use,
    strata_var = match_strata_var,
    apply_pair_screen = apply_pair_total_screen,
    apply_mediator_screen = apply_mediator_ab_screen
  )
  combos_shared <- combo_pack$combos
  write.csv(combo_pack$pair_stats, file.path(out_dir, "screen_total_pairs.csv"), row.names = FALSE)
  write.csv(combo_pack$mediator_stats, file.path(out_dir, "screen_mediator_pairs.csv"), row.names = FALSE)
  cat(
    "筛选完成 | pair总数:", nrow(combo_pack$pair_stats),
    "| pair通过:", sum(combo_pack$pair_stats$pass_total, na.rm = TRUE),
    "| mediator候选记录:", nrow(combo_pack$mediator_stats),
    "| 最终组合数:", nrow(combos_shared), "\n"
  )
}

if (isTRUE(run_med_analysis)) {
  combos_med <- combos_shared
  if (is.null(combos_med) || nrow(combos_med) == 0L) stop("No mediation triples left after screening.")
  cat("传统中介分析组合数:", nrow(combos_med),
      "| cores:", n_cores_use,
      "| n_boot:", n_boot_main,
      "| boot_screen_mode:", boot_screen_mode,
      "| boot_p_screen:", boot_p_screen, "\n")
  if (n_cores_use >= 2L) {
    cl_med <- parallel::makeCluster(n_cores_use)
    parallel::clusterEvalQ(cl_med, library(survival))
    parallel::clusterExport(
      cl_med,
      c("dat_ana", "covar_use", "combos_med", "join_terms", "run_one_mediation",
        "classify_mediation_type", "compute_prop_metrics",
        "n_boot_main", "boot_screen_mode", "boot_p_screen", "match_strata_var"),
      envir = environment()
    )
    res_list_med <- tryCatch(
      parallel::parLapply(cl_med, seq_len(nrow(combos_med)), function(i) {
        run_one_mediation(dat_ana, combos_med$exposure[i], combos_med$mediator[i], combos_med$outcome[i],
                          strata_var = match_strata_var,
                          covars = covar_use, n_boot = n_boot_main,
                          boot_screen_mode = boot_screen_mode, boot_p_screen = boot_p_screen)
      }),
      finally = parallel::stopCluster(cl_med)
    )
  } else {
    res_list_med <- lapply(seq_len(nrow(combos_med)), function(i) {
      run_one_mediation(dat_ana, combos_med$exposure[i], combos_med$mediator[i], combos_med$outcome[i],
                        strata_var = match_strata_var,
                        covars = covar_use, n_boot = n_boot_main,
                        boot_screen_mode = boot_screen_mode, boot_p_screen = boot_p_screen)
    })
  }
  out_list <- Filter(Negate(is.null), res_list_med)
  if (n_cores_use >= 2L) cat("传统中介分析已使用", n_cores_use, "核并行完成。\n")
  
  if (length(out_list) == 0L) stop("No mediation model converged.")
  results_med <- do.call(rbind, out_list)
  
  # FDR 校正（针对间接效应 p 值）
  results_med$p_ind_FDR <- p.adjust(results_med$p_ind, method = "BH")
  # 路径显著性标记：X-M(a), M-Y(b), X-Y(total/c)
  results_med$sig_xm <- is.finite(results_med$p_a) & (results_med$p_a < path_sig_alpha)
  results_med$sig_my <- is.finite(results_med$p_b) & (results_med$p_b < path_sig_alpha)
  results_med$sig_xy <- is.finite(results_med$p_total) & (results_med$p_total < path_sig_alpha)
  results_med$sig_xm_my_xy <- results_med$sig_xm & results_med$sig_my & results_med$sig_xy
  
  # 按 p_ind 排序
  results_med <- results_med[order(results_med$p_ind), ]
  
  # 简洁结果列
  results_med$ind_CI <- sprintf("%.4f (%.4f, %.4f)",
                                results_med$ind_effect,
                                results_med$ci_ind_low,
                                results_med$ci_ind_high)
  
  # 保存完整结果
  write.csv(results_med,
            file.path(out_dir, "mediation_full.csv"),
            row.names = FALSE)
  
  # 保存简洁表（暴露、中介、结局、路径 a/b、总/直接/间接、中介比例、p 及 FDR）
  cols_short <- c("exposure", "mediator", "outcome",
                  "coef_a", "se_a", "p_a", "coef_b", "se_b", "p_b",
                  "coef_c", "se_c", "p_c", "ci_total_low", "ci_total_high", "total_effect", "se_total", "p_total",
                  "coef_cp", "se_cp", "p_cp", "ci_direct_low", "ci_direct_high",
                  "ind_effect", "se_ind", "ind_CI", "ci_ind_method", "p_ind", "p_ind_FDR",
                  "ind_from_diff", "prop_med", "prop_med_raw", "prop_med_report", "prop_med_interpretable", "prop_med_note",
                  "mediation_type", "suppressive_role",
                  "OR_total", "OR_direct", "OR_ind_approx",
                  "sig_xm", "sig_my", "sig_xy", "sig_xm_my_xy",
                  "n_used", "events", "n_boot")
  write.csv(results_med[, cols_short[cols_short %in% names(results_med)]],
            file.path(out_dir, "mediation_summary.csv"),
            row.names = FALSE)

  # 路径关联结果表：用于单独展示 X-M、M-Y、X-Y 的回归关联
  cols_path_trad <- c("exposure", "mediator", "outcome",
                      "coef_a", "se_a", "p_a",
                      "coef_b", "se_b", "p_b",
                      "coef_c", "se_c", "p_c", "ci_total_low", "ci_total_high",
                      "n_used", "events",
                      "sig_xm", "sig_my", "sig_xy", "sig_xm_my_xy")
  write.csv(results_med[, cols_path_trad[cols_path_trad %in% names(results_med)]],
            file.path(out_dir, "association_paths_traditional.csv"),
            row.names = FALSE)

  # 仅保留 X-M、M-Y、X-Y 均显著的中介结果（便于论文主表展示）
  results_med_sig3 <- results_med[results_med$sig_xm_my_xy, , drop = FALSE]
  write.csv(results_med_sig3,
            file.path(out_dir, "mediation_full_traditional_sig_xm_my_xy.csv"),
            row.names = FALSE)
  write.csv(results_med_sig3[, cols_short[cols_short %in% names(results_med_sig3)]],
            file.path(out_dir, "mediation_summary_traditional_sig_xm_my_xy.csv"),
            row.names = FALSE)
  report_trad_sig3 <- data.frame(
    outcome = results_med_sig3$outcome,
    exposure = results_med_sig3$exposure,
    mediator = results_med_sig3$mediator,
    `total effect β (95%CI)` = sprintf("%.3f (%.3f, %.3f)",
                                       results_med_sig3$coef_c,
                                       results_med_sig3$ci_total_low,
                                       results_med_sig3$ci_total_high),
    `total effect P` = results_med_sig3$p_total,
    `direct effect β (95%CI)` = sprintf("%.3f (%.3f, %.3f)",
                                        results_med_sig3$coef_cp,
                                        results_med_sig3$ci_direct_low,
                                        results_med_sig3$ci_direct_high),
    `direct effect P` = results_med_sig3$p_cp,
    `indirect effect β (95%CI)` = sprintf("%.3f (%.3f, %.3f)",
                                          results_med_sig3$ind_effect,
                                          results_med_sig3$ci_ind_low,
                                          results_med_sig3$ci_ind_high),
    `indirect effect P` = results_med_sig3$p_ind,
    `proportion of mediation (%)` = 100 * results_med_sig3$prop_med,
    stringsAsFactors = FALSE
  )
  write.csv(report_trad_sig3,
            file.path(out_dir, "mediation_report_traditional_sig_xm_my_xy.csv"),
            row.names = FALSE)
  cat("传统中介：X-M、M-Y、X-Y 均显著的组合数 =", nrow(results_med_sig3), "\n")
  
  # 在控制台打印有显著间接效应的组合（按 p_ind_FDR < 0.05 或 p_ind < 0.05）
  sig <- results_med[results_med$p_ind < 0.05 | results_med$p_ind_FDR < 0.05, ]
  if (nrow(sig) > 0L) {
    cat("显著间接效应的暴露-中介-结局组合 (p_ind < 0.05 或 p_ind_FDR < 0.05):\n")
    print(sig[, c("exposure", "mediator", "outcome", "ind_effect", "ind_CI", "prop_med", "prop_med_report", "mediation_type", "p_ind", "p_ind_FDR")])
  } else {
    cat("当前无 p_ind 或 p_ind_FDR < 0.05 的显著间接效应。\n")
  }
  
  # 说明：
  # - coef_a: 暴露对中介的效应（路径 a）
  # - coef_b: 中介对结局的效应（路径 b，控制暴露）
  # - coef_c: 暴露对结局的总效应（路径 c）
  # - coef_cp: 暴露对结局的直接效应（路径 c'）
  # - ind_effect: 间接效应 a*b；OR_ind_approx = 2^ind_effect 为近似间接 OR
  # - prop_med: 原始中介比例 (c - c')/c（不再截断到 0~1）
  # - prop_med_report: 仅当总效应显著且方向一致时才建议报告，否则置为 NA
  # - mediation_type: 中介类型（含 competitive/suppression）
} else {
  cat("\n已跳过传统中介分析（run_med_analysis = FALSE）。\n")
}

# =============================================================================
# 结构方程模型 SEM 中介路径分析
# 路径：X -> M (a), X -> Y (c'), M -> Y (b)；间接效应 := a*b，总效应 := c' + a*b
# 说明：SEM 中通过 cluster=Subclass 纳入匹配层，使用 cluster-robust 标准误
#      （与 clogit 的条件似然并不完全等价，但可保持匹配层一致性）
# =============================================================================
if (isTRUE(run_sem_analysis)) {
  if (!requireNamespace("lavaan", quietly = TRUE)) {
    stop("run_sem_analysis=TRUE 但未安装 lavaan，请先安装：install.packages('lavaan')")
  }
  sem_detected_cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
  if (!is.finite(sem_detected_cores) || sem_detected_cores < 1L) {
    ns_parallel <- asNamespace("parallel")
    unlockBinding("detectCores", ns_parallel)
    assign("detectCores", function(all.tests = FALSE, logical = TRUE) 1L, envir = ns_parallel)
    lockBinding("detectCores", ns_parallel)
    cat("检测到 detectCores() = NA，已为 lavaan 临时覆盖为 1 核。\n")
  }
  library(lavaan)
  
  to_numeric_sem <- function(z) {
    if (is.factor(z)) return(as.numeric(z))
    if (is.character(z)) {
      out <- suppressWarnings(as.numeric(z))
      if (all(!is.finite(out))) out <- as.numeric(factor(z))
      return(out)
    }
    suppressWarnings(as.numeric(z))
  }
  
  to_binary01 <- function(y) {
    y_num <- suppressWarnings(as.numeric(as.character(y)))
    if (all(!is.finite(y_num))) y_num <- as.numeric(factor(y)) - 1L
    fin <- is.finite(y_num)
    u <- sort(unique(y_num[fin]))
    if (length(u) != 2L) return(NULL)
    out <- rep(NA_integer_, length(y_num))
    out[fin] <- as.integer(y_num[fin] == u[2L])
    out
  }
  
  has_variation <- function(x) {
    x <- x[is.finite(x)]
    length(unique(x)) >= 2L
  }
  
  sem_try_fit <- function(model_sem, d_sem, estimator_name, cluster_var = NULL) {
    sem_warn <- character(0)
    sem_err <- NA_character_
    sem_args <- list(
      model = model_sem,
      data = d_sem,
      ordered = "Y_use",
      estimator = estimator_name,
      std.lv = FALSE
    )
    if (!is.null(cluster_var) && nzchar(cluster_var) && cluster_var %in% names(d_sem)) {
      sem_args$cluster <- cluster_var
    }
    fit_obj <- withCallingHandlers(
      tryCatch(
        do.call(sem, sem_args),
        error = function(e) {
          sem_err <<- conditionMessage(e)
          e
        }
      ),
      warning = function(w) {
        sem_warn <<- c(sem_warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
    if (inherits(fit_obj, "error")) {
      return(list(fit = NULL, converged = FALSE, estimator = estimator_name, warnings = sem_warn, error = sem_err))
    }
    converged <- isTRUE(tryCatch(lavInspect(fit_obj, "converged"), error = function(e) FALSE))
    list(fit = fit_obj, converged = converged, estimator = estimator_name, warnings = sem_warn, error = sem_err)
  }
  
  run_one_sem_mediation <- function(data, x_var, m_var, y_var, covars = covar_use,
                                    strata_var = match_strata_var,
                                    use_covars = TRUE,
                                    use_subclass_cluster = sem_use_subclass_cluster) {
    if (isTRUE(use_subclass_cluster) && !(strata_var %in% names(data))) return(NULL)
    vars_needed <- unique(c(x_var, m_var, y_var, covars, strata_var))
    vars_needed <- vars_needed[vars_needed %in% names(data)]
    if (!all(c(x_var, m_var, y_var) %in% vars_needed)) return(NULL)
    
    d <- data[, vars_needed, drop = FALSE]
    cc_raw <- complete.cases(d[, c(x_var, m_var, y_var), drop = FALSE])
    d <- d[cc_raw, , drop = FALSE]
    if (nrow(d) < 40L) return(NULL)
    
    x_num <- to_numeric_sem(d[[x_var]])
    m_num <- to_numeric_sem(d[[m_var]])
    y_bin <- to_binary01(d[[y_var]])
    subclass_raw <- if (strata_var %in% names(d)) d[[strata_var]] else NULL
    if (is.null(y_bin)) return(NULL)
    
    keep <- is.finite(x_num) & is.finite(m_num) & !is.na(y_bin)
    if (!is.null(subclass_raw)) keep <- keep & !is.na(subclass_raw)
    if (sum(keep) < 40L) return(NULL)
    
    d_sem <- data.frame(
      X_use = x_num[keep],
      M_use = m_num[keep],
      Y_use = ordered(y_bin[keep], levels = c(0L, 1L)),
      Subclass_use = as.factor(subclass_raw[keep])
    )
    if (!has_variation(d_sem$X_use) || !has_variation(d_sem$M_use)) return(NULL)
    if (length(unique(as.integer(as.character(d_sem$Y_use)))) < 2L) return(NULL)
    if (isTRUE(use_subclass_cluster) && length(unique(d_sem$Subclass_use)) < 2L) return(NULL)
    
    cov_avail <- character(0)
    if (use_covars) {
      for (cv in covars[covars %in% names(d)]) {
        cv_num <- to_numeric_sem(d[[cv]])
        cv_num <- cv_num[keep]
        if (all(is.finite(cv_num)) && has_variation(cv_num)) {
          d_sem[[cv]] <- cv_num
          cov_avail <- c(cov_avail, cv)
        }
      }
    }
    
    rhs_m <- join_terms(c("a*X_use", cov_avail))
    rhs_y <- join_terms(c("b*M_use", "cp*X_use", cov_avail))
    model_sem <- paste0(
      "M_use ~ ", rhs_m, "\n",
      "Y_use ~ ", rhs_y, "\n",
      "indirect := a*b\n",
      "total := cp + a*b\n"
    )
    
    cluster_tag <- if (isTRUE(use_subclass_cluster)) "Subclass_use" else NULL
    used_cluster <- isTRUE(use_subclass_cluster)
    fit_out <- sem_try_fit(model_sem, d_sem, estimator_name = "WLSMV", cluster_var = cluster_tag)
    if (!fit_out$converged) {
      fit_out <- sem_try_fit(model_sem, d_sem, estimator_name = "MLR", cluster_var = cluster_tag)
    }
    # 若 cluster 设定导致不可估，则自动降级为不加 cluster 重试
    if (!fit_out$converged && isTRUE(use_subclass_cluster)) {
      fit_out <- sem_try_fit(model_sem, d_sem, estimator_name = "WLSMV", cluster_var = NULL)
      used_cluster <- FALSE
      if (!fit_out$converged) {
        fit_out <- sem_try_fit(model_sem, d_sem, estimator_name = "MLR", cluster_var = NULL)
      }
    }
    if (!fit_out$converged || is.null(fit_out$fit)) return(NULL)
    fit <- fit_out$fit
    
    pe <- parameterEstimates(fit, ci = TRUE)
    get_par <- function(lab) {
      r <- pe[pe$label == lab | pe$lhs == lab, ]
      if (nrow(r) < 1L) return(c(est = NA, se = NA, p = NA, ci_low = NA, ci_high = NA))
      r <- r[1L, , drop = FALSE]
      c(est = r$est, se = r$se, p = r$pvalue,
        ci_low = if ("ci.lower" %in% names(r)) r$ci.lower else NA,
        ci_high = if ("ci.upper" %in% names(r)) r$ci.upper else NA)
    }
    
    a_   <- get_par("a")
    b_   <- get_par("b")
    cp_  <- get_par("cp")
    ind_ <- get_par("indirect")
    tot_ <- get_par("total")
    
    pe_std <- tryCatch(standardizedSolution(fit), error = function(e) NULL)
    if (!is.null(pe_std) && "std.all" %in% names(pe_std)) {
      get_std <- function(lab) {
        r <- pe_std[pe_std$label == lab | pe_std$lhs == lab, ]
        if (nrow(r) < 1L) return(NA_real_)
        r$std.all[1L]
      }
      std_a  <- get_std("a")
      std_b  <- get_std("b")
      std_cp <- get_std("cp")
      std_ind <- if (!is.na(std_a) && !is.na(std_b)) std_a * std_b else NA_real_
    } else {
      std_a <- std_b <- std_cp <- std_ind <- NA_real_
    }
    
    fm <- tryCatch(fitMeasures(fit, c("chisq", "df", "pvalue", "cfi", "tli", "rmsea", "srmr")), error = function(e) NULL)
    if (is.null(fm)) {
      chisq <- chisq_df <- chisq_pvalue <- cfi <- tli <- rmsea <- srmr <- NA_real_
    } else {
      chisq        <- as.numeric(fm["chisq"])
      chisq_df     <- as.numeric(fm["df"])
      chisq_pvalue <- as.numeric(fm["pvalue"])
      cfi          <- as.numeric(fm["cfi"])
      tli          <- as.numeric(fm["tli"])
      rmsea        <- as.numeric(fm["rmsea"])
      srmr         <- as.numeric(fm["srmr"])
    }
    
    tot_est <- as.numeric(tot_["est"])
    se_tot  <- as.numeric(tot_["se"])
    p_tot   <- as.numeric(tot_["p"])
    ci_tot_low  <- as.numeric(tot_["ci_low"])
    ci_tot_high <- as.numeric(tot_["ci_high"])
    cp_est <- as.numeric(cp_["est"])
    prop_info <- compute_prop_metrics(tot_est, cp_est, p_tot)
    mediation_type <- classify_mediation_type(cp_est, as.numeric(ind_["est"]), as.numeric(cp_["p"]), as.numeric(ind_["p"]))
    suppressive_role <- is.finite(cp_est) && is.finite(as.numeric(ind_["est"])) && (cp_est * as.numeric(ind_["est"]) < 0)
    y_events <- sum(as.integer(as.character(d_sem$Y_use)) == 1L, na.rm = TRUE)
    
    data.frame(
      exposure  = x_var,
      mediator  = m_var,
      outcome   = y_var,
      method    = "SEM",
      with_covars = use_covars,
      strata_var = strata_var,
      use_subclass_cluster = isTRUE(used_cluster),
      n_clusters = length(unique(d_sem$Subclass_use)),
      estimator = fit_out$estimator,
      n_covars_used = length(cov_avail),
      coef_a    = a_["est"],
      se_a      = a_["se"],
      p_a       = a_["p"],
      ci_a_low  = a_["ci_low"],
      ci_a_high = a_["ci_high"],
      std_a     = std_a,
      coef_b    = b_["est"],
      se_b      = b_["se"],
      p_b       = b_["p"],
      ci_b_low  = b_["ci_low"],
      ci_b_high = b_["ci_high"],
      std_b     = std_b,
      coef_cp   = cp_["est"],
      se_cp     = cp_["se"],
      p_cp      = cp_["p"],
      ci_cp_low = cp_["ci_low"],
      ci_cp_high = cp_["ci_high"],
      std_cp    = std_cp,
      ind_effect = ind_["est"],
      se_ind    = ind_["se"],
      ci_ind_low  = ind_["ci_low"],
      ci_ind_high = ind_["ci_high"],
      p_ind     = ind_["p"],
      std_ind   = std_ind,
      total_effect = tot_est,
      se_total  = se_tot,
      p_total   = p_tot,
      ci_total_low  = ci_tot_low,
      ci_total_high = ci_tot_high,
      ind_from_diff = tot_est - cp_est,
      prop_med  = prop_info$prop_raw,
      prop_med_raw = prop_info$prop_raw,
      prop_med_report = prop_info$prop_report,
      prop_med_interpretable = prop_info$prop_interpretable,
      prop_med_note = prop_info$prop_note,
      mediation_type = mediation_type,
      suppressive_role = suppressive_role,
      chisq     = chisq,
      chisq_df  = chisq_df,
      chisq_pvalue = chisq_pvalue,
      cfi       = cfi,
      tli       = tli,
      rmsea     = rmsea,
      srmr      = srmr,
      n_used    = nrow(d_sem),
      events    = y_events,
      stringsAsFactors = FALSE
    )
  }
  
  combos_sem <- combos_shared
  if (is.null(combos_sem) || nrow(combos_sem) == 0L) stop("No SEM triples left after screening.")
  cat("SEM 中介分析组合数:", nrow(combos_sem), "| cores:", n_cores_use, "\n")
  
  run_sem_one_combo <- function(i, use_covars_arg) {
    run_one_sem_mediation(
      dat_ana,
      combos_sem$exposure[i],
      combos_sem$mediator[i],
      combos_sem$outcome[i],
      covars = covar_use,
      strata_var = match_strata_var,
      use_covars = use_covars_arg,
      use_subclass_cluster = sem_use_subclass_cluster
    )
  }
  
  run_sem_batch <- function(use_covars_arg) {
    n_combo <- nrow(combos_sem)
    if (n_cores_use >= 2L) {
      cl_sem <- parallel::makeCluster(n_cores_use)
      parallel::clusterEvalQ(cl_sem, library(lavaan))
      parallel::clusterExport(
        cl_sem,
        c("dat_ana", "covar_use", "combos_sem", "join_terms", "to_numeric_sem", "to_binary01",
          "has_variation", "sem_try_fit", "run_one_sem_mediation", "run_sem_one_combo",
          "classify_mediation_type", "compute_prop_metrics", "match_strata_var", "sem_use_subclass_cluster"),
        envir = environment()
      )
      on.exit(parallel::stopCluster(cl_sem), add = TRUE)
      parallel::parLapply(cl_sem, seq_len(n_combo), function(i) run_sem_one_combo(i, use_covars_arg))
    } else {
      out <- vector("list", n_combo)
      tag <- if (isTRUE(use_covars_arg)) "SEM(带协变量)" else "SEM(无协变量)"
      for (i in seq_len(n_combo)) {
        if (i == 1L || i %% 25L == 0L || i == n_combo) cat(tag, "进度:", i, "/", n_combo, "\n")
        out[[i]] <- run_sem_one_combo(i, use_covars_arg)
      }
      out
    }
  }
  
  res_list_sem <- run_sem_batch(TRUE)
  out_sem <- Filter(Negate(is.null), res_list_sem)
  if (length(out_sem) == 0L) {
    cat("SEM（带协变量）未得到可用模型，尝试无协变量简化模型...\n")
    res_list_sem2 <- run_sem_batch(FALSE)
    out_sem <- Filter(Negate(is.null), res_list_sem2)
  }
  
  if (length(out_sem) > 0L) {
    results_sem <- do.call(rbind, out_sem)
    results_sem$p_ind_FDR <- p.adjust(results_sem$p_ind, method = "BH")
    # 路径显著性标记：X-M(a), M-Y(b), X-Y(total)
    results_sem$sig_xm <- is.finite(results_sem$p_a) & (results_sem$p_a < path_sig_alpha)
    results_sem$sig_my <- is.finite(results_sem$p_b) & (results_sem$p_b < path_sig_alpha)
    results_sem$sig_xy <- is.finite(results_sem$p_total) & (results_sem$p_total < path_sig_alpha)
    results_sem$sig_xm_my_xy <- results_sem$sig_xm & results_sem$sig_my & results_sem$sig_xy
    results_sem <- results_sem[order(results_sem$p_ind), ]
    results_sem$ind_CI <- sprintf("%.4f (%.4f, %.4f)",
                                  results_sem$ind_effect,
                                  results_sem$ci_ind_low,
                                  results_sem$ci_ind_high)
    
    write.csv(results_sem,
              file.path(out_dir, "mediation_SEM_full.csv"),
              row.names = FALSE)
    
    cols_sem <- c("exposure", "mediator", "outcome", "method", "with_covars", "strata_var", "use_subclass_cluster", "n_clusters", "estimator", "n_covars_used",
                  "coef_a", "se_a", "p_a", "ci_a_low", "ci_a_high", "std_a",
                  "coef_b", "se_b", "p_b", "ci_b_low", "ci_b_high", "std_b",
                  "coef_cp", "se_cp", "p_cp", "ci_cp_low", "ci_cp_high", "std_cp",
                  "ind_effect", "ind_CI", "p_ind", "p_ind_FDR", "std_ind",
                  "total_effect", "se_total", "p_total", "ci_total_low", "ci_total_high",
                  "ind_from_diff", "prop_med", "prop_med_raw", "prop_med_report", "prop_med_interpretable", "prop_med_note",
                  "mediation_type", "suppressive_role",
                  "sig_xm", "sig_my", "sig_xy", "sig_xm_my_xy",
                  "chisq", "chisq_df", "chisq_pvalue", "cfi", "tli", "rmsea", "srmr",
                  "n_used", "events")
    write.csv(results_sem[, cols_sem[cols_sem %in% names(results_sem)]],
              file.path(out_dir, "mediation_SEM_summary.csv"),
              row.names = FALSE)

    # 路径关联结果表：用于单独展示 X-M、M-Y、X-Y 的回归关联
    cols_path_sem <- c("exposure", "mediator", "outcome",
                       "coef_a", "se_a", "p_a", "ci_a_low", "ci_a_high",
                       "coef_b", "se_b", "p_b", "ci_b_low", "ci_b_high",
                       "total_effect", "se_total", "p_total", "ci_total_low", "ci_total_high",
                       "n_used", "events",
                       "sig_xm", "sig_my", "sig_xy", "sig_xm_my_xy")
    write.csv(results_sem[, cols_path_sem[cols_path_sem %in% names(results_sem)]],
              file.path(out_dir, "association_paths_SEM.csv"),
              row.names = FALSE)

    # 仅保留 X-M、M-Y、X-Y 均显著的中介结果（便于论文主表展示）
    results_sem_sig3 <- results_sem[results_sem$sig_xm_my_xy, , drop = FALSE]
    write.csv(results_sem_sig3,
              file.path(out_dir, "mediation_SEM_full_sig_xm_my_xy.csv"),
              row.names = FALSE)
    write.csv(results_sem_sig3[, cols_sem[cols_sem %in% names(results_sem_sig3)]],
              file.path(out_dir, "mediation_SEM_summary_sig_xm_my_xy.csv"),
              row.names = FALSE)
    report_sem_sig3 <- data.frame(
      outcome = results_sem_sig3$outcome,
      exposure = results_sem_sig3$exposure,
      mediator = results_sem_sig3$mediator,
      `total effect β (95%CI)` = sprintf("%.3f (%.3f, %.3f)",
                                         results_sem_sig3$total_effect,
                                         results_sem_sig3$ci_total_low,
                                         results_sem_sig3$ci_total_high),
      `total effect P` = results_sem_sig3$p_total,
      `direct effect β (95%CI)` = sprintf("%.3f (%.3f, %.3f)",
                                          results_sem_sig3$coef_cp,
                                          results_sem_sig3$ci_cp_low,
                                          results_sem_sig3$ci_cp_high),
      `direct effect P` = results_sem_sig3$p_cp,
      `indirect effect β (95%CI)` = sprintf("%.3f (%.3f, %.3f)",
                                            results_sem_sig3$ind_effect,
                                            results_sem_sig3$ci_ind_low,
                                            results_sem_sig3$ci_ind_high),
      `indirect effect P` = results_sem_sig3$p_ind,
      `proportion of mediation (%)` = 100 * results_sem_sig3$prop_med,
      stringsAsFactors = FALSE
    )
    write.csv(report_sem_sig3,
              file.path(out_dir, "mediation_report_SEM_sig_xm_my_xy.csv"),
              row.names = FALSE)
    cat("SEM 中介：X-M、M-Y、X-Y 均显著的组合数 =", nrow(results_sem_sig3), "\n")
    
    cat("\nSEM 中介路径分析完成，结果已保存至 output_csv/Med_results/mediation_SEM_*.csv\n")
    cat("说明：SEM 中 Y 为二分类时使用 WLSMV/MLR + cluster(Subclass) 稳健标准误；与 clogit 的条件似然估计不完全等价。\n")
  } else {
    cat("\nSEM 中介模型均未得到可用结果。请重点检查：结局是否为二分类、变量是否有变异、缺失是否过多。\n")
  }
} else {
  cat("\n已跳过 SEM 中介分析（run_sem_analysis = FALSE）。如需运行请改为 TRUE。\n")
}
