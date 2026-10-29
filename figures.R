# figures.R ---------------------------------------------------------------
# Publication-style figures for:
#   "Oral versus transdermal estrogen and triglycerides in women starting
#    menopausal hormone therapy: a target trial emulation study"
#
# SELF-CONTAINED DEMO: every input below is SIMULATED random data (fixed
# seed). The script reproduces the figure *design* used in the article,
# not its results. No real data are read, derived, or shipped.
#
# Figures produced (PDF vector + TIFF, written to ./figures/):
#   figure_participant_flow         Participant flow diagram (140 mm)
#   figure_forest_primary           Primary forest plot (140 mm)
#   figure_forest_full              Extended forest plot (140 mm)
#   figure_tg_baseline_followup     Median (IQR) TG/HDL/LDL, baseline vs follow-up (190 mm)
#   figure_balance_diagnostics      PS overlap, Love plot, weight distribution (190 mm)
#   figure_mice_diagnostics         Imputation trace + density diagnostics (190 mm)
#   figure_supp_diagnostics         Combined diagnostics supplement (190 mm)
#   figure_supp_additional_results  Combined results supplement (190 mm)
#
# Requirements: data.table, ggplot2, patchwork, ragg

library(data.table)
library(ggplot2)
library(patchwork)

set.seed(20250101)  # fixed seed, fully reproducible

# ===================== Output directory =====================
fig_dir <- "figures"
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ===================== Shared style =====================
FONT     <- "Arial"
LINEW    <- 0.35
OKABE_ITO <- c(oral = "#0072B2", tdm = "#D55E00")  # arm colours (colour-blind safe)

theme_pub <- function(base_size = 7) {
  theme_minimal(base_size = base_size, base_family = FONT) %+replace%
    theme(
      panel.grid       = element_blank(),
      panel.background = element_blank(),
      axis.line        = element_line(colour = "black", linewidth = LINEW / .pt),
      axis.ticks       = element_line(colour = "black", linewidth = LINEW / .pt),
      axis.text        = element_text(size = base_size, colour = "black"),
      axis.title       = element_text(size = base_size + 0.5),
      legend.position  = "none",
      strip.text       = element_text(size = base_size, face = "bold"),
      strip.background = element_blank(),
      plot.tag         = element_text(size = 9, face = "bold", family = FONT)
    )
}

W_SINGLE <- 90   # single column, mm
W_1_5    <- 140  # 1.5 column, mm
W_DOUBLE <- 190  # double column, mm

save_pdf <- function(plot, name, w_mm, h_mm) {
  cairo_pdf(file.path(fig_dir, paste0(name, ".pdf")),
            width = w_mm / 25.4, height = h_mm / 25.4, family = FONT)
  print(plot)
  dev.off()
}
save_tiff <- function(plot, name, w_mm, h_mm, dpi = 1000) {
  ragg::agg_tiff(file.path(fig_dir, paste0(name, ".tif")),
                 width = w_mm, height = h_mm, units = "mm", res = dpi,
                 compression = "lzw")
  print(plot)
  dev.off()
}

# ===================== SIMULATED DATA =====================
# All quantities below are randomly generated for demonstration; they carry
# no information about any real cohort.

M_IMP  <- 20   # number of imputed data sets (dimension of weight/PS arrays)
N_ITER <- 20   # imputation iterations shown in trace plots

# ---- Person-level simulated cohort (analysis population) ----
N_POP <- 900
sim <- data.table(
  subject_id      = seq_len(N_POP),
  age_at_index   = round(rnorm(N_POP, 58, 6)),
  index_year     = sample(2012:2022, N_POP, replace = TRUE),
  lab_tg_base_value  = rlnorm(N_POP, log(130), 0.45),   # mg/dL
  lab_hdl_base_value = pmax(20, rnorm(N_POP, 58, 14)),  # mg/dL
  lab_ldl_base_value = pmax(40, rnorm(N_POP, 115, 30))  # mg/dL
)
sim[, log_tg_base := log(lab_tg_base_value)]
bmi_full <- rlnorm(N_POP, log(28), 0.18)
sbp_full <- rnorm(N_POP, 125, 14)

z <- function(x) (x - mean(x)) / sd(x)
# Treatment arm depends weakly on covariates (gives the PS some signal)
sim[, arm := fifelse(runif(.N) < plogis(0.6 * z(bmi_full) + 0.5 * z(log_tg_base) -
                                         0.4 * z(age_at_index) + 0.2),
                     "oral", "transdermal")]
# Follow-up measurement availability (defines the outcome sample)
sim[, in_paired := rbinom(.N, 1, plogis(0.8 - 0.5 * z(age_at_index)))]

# Binary baseline covariates (flags)
mk_flag <- function(p) rbinom(N_POP, 1, p)
sim[, cx_vte_ever         := mk_flag(0.03)]
sim[, cx_chd_ever         := mk_flag(0.04)]
sim[, cx_htn_base         := mk_flag(0.35)]
sim[, cx_t2dm_base        := mk_flag(0.12)]
sim[, rx_statin_base      := mk_flag(0.18)]
sim[, rx_progestogen_base := mk_flag(0.15)]

# Item-level missingness in three covariates (to be imputed below)
sim[, vital_bmi_base_value := bmi_full]
sim[, vital_sbp_base_value := sbp_full]
mis_bmi <- runif(N_POP) < 0.20   # ~20% missing BMI
mis_sbp <- runif(N_POP) < 0.12   # ~12% missing SBP
mis_ltg <- runif(N_POP) < 0.025  # ~2.5% missing log TG (transformed covariate)
sim[mis_bmi, vital_bmi_base_value := NA]
sim[mis_sbp, vital_sbp_base_value := NA]
sim[mis_ltg, log_tg_base          := NA]

# ---- Simulated multiple-imputation frames ----
# Naive stochastic imputation, one draw per imputed data set
bmi_mu <- mean(bmi_full); bmi_sd <- sd(bmi_full)
sbp_mu <- mean(sbp_full); sbp_sd <- sd(sbp_full)
ltg_mu <- mean(sim$log_tg_base, na.rm = TRUE)
ltg_sd <- sd(sim$log_tg_base, na.rm = TRUE)
fr_p <- lapply(seq_len(M_IMP), function(b) {
  x <- copy(sim)
  x[mis_bmi, vital_bmi_base_value := pmax(16, rnorm(.N, bmi_mu, bmi_sd))]
  x[mis_sbp, vital_sbp_base_value := pmax(90, rnorm(.N, sbp_mu, sbp_sd))]
  x[mis_ltg, log_tg_base          := rnorm(.N, ltg_mu, ltg_sd)]
  x
})

# ---- Outcome sample, treatment indicator, arm map ----
ip_vec    <- sim$in_paired == 1
pid_wt    <- sim$subject_id[ip_vec]
arm_map   <- data.table(subject_id = pid_wt, arm = sim$arm[ip_vec])
arm_num_p <- as.integer(sim$arm[ip_vec] == "oral")

# ---- Propensity scores and ATE weights, one fit per imputation ----
# Logistic PS; weights = 1/ps (oral) or 1/(1 - ps) (transdermal)
ps_w_list <- lapply(seq_len(M_IMP), function(b) {
  f  <- fr_p[[b]]
  db <- data.frame(
    armoral = arm_num_p,
    age     = f$age_at_index[ip_vec],
    idxyr   = f$index_year[ip_vec],
    bmi     = f$vital_bmi_base_value[ip_vec],
    sbp     = f$vital_sbp_base_value[ip_vec],
    logtg   = f$log_tg_base[ip_vec],
    cx_htn_base         = f$cx_htn_base[ip_vec],
    cx_t2dm_base        = f$cx_t2dm_base[ip_vec],
    cx_vte_ever         = f$cx_vte_ever[ip_vec],
    cx_chd_ever         = f$cx_chd_ever[ip_vec],
    rx_statin_base      = f$rx_statin_base[ip_vec],
    rx_progestogen_base = f$rx_progestogen_base[ip_vec])
  fit <- glm(armoral ~ ., data = db, family = binomial())
  ps  <- as.numeric(fitted(fit))
  w   <- ifelse(arm_num_p == 1, 1 / ps, 1 / (1 - ps))
  list(ps = ps, w = w)
})
ps_main <- ps_w_list[[1]]$ps
wC_main <- as.matrix(sapply(ps_w_list, function(l) l$w))  # outcome sample x M_IMP

# Absolute SMD of one covariate under a given weight vector
calc_smd_fig <- function(x, arm_num, w) {
  use <- !is.na(x) & !is.na(w) & !is.na(arm_num)
  x <- as.numeric(x[use]); arm_num <- arm_num[use]; w <- w[use]
  if (length(x) < 10 || !any(arm_num == 1) || !any(arm_num == 0)) return(NA_real_)
  m1 <- sum(w[arm_num == 1] * x[arm_num == 1]) / sum(w[arm_num == 1])
  m0 <- sum(w[arm_num == 0] * x[arm_num == 0]) / sum(w[arm_num == 0])
  sd1 <- sqrt(sum(w[arm_num == 1] * (x[arm_num == 1] - m1)^2) / sum(w[arm_num == 1]))
  sd0 <- sqrt(sum(w[arm_num == 0] * (x[arm_num == 0] - m0)^2) / sum(w[arm_num == 0]))
  sp <- sqrt((sd1^2 + sd0^2) / 2); if (sp == 0 || !is.finite(sp)) return(NA_real_)
  abs(m1 - m0) / sp
}
# Mean SMD across imputations, weighted or unweighted
smd_imp_mean <- function(v, weighted) {
  mean(sapply(seq_len(M_IMP), function(b) {
    w <- if (weighted) wC_main[, b] else rep(1, sum(ip_vec))
    calc_smd_fig(fr_p[[b]][[v]][ip_vec], arm_num_p, w)
  }))
}

# ---- Derived sample counts (flow diagram and forest plots) ----
n_base_oral <- sim[arm == "oral", .N]
n_base_tdm  <- sim[arm == "transdermal", .N]
n_oral <- sum(arm_num_p == 1)   # analysed, oral arm
n_tdm  <- sum(arm_num_p == 0)   # analysed, transdermal arm
n_fu_oral <- n_base_oral - n_oral
n_fu_tdm  <- n_base_tdm  - n_tdm
n_base_tg <- N_POP
n_all     <- 3000               # simulated eligible starters
n_no_base <- n_all - n_base_tg

# ---- Simulated follow-up lipid values (outcome sample only) ----
fu_mult <- ifelse(arm_map$arm == "oral", 1.06, 1.09)  # arbitrary simulated change
sim_p <- copy(sim)[ip_vec]
sim_p[, lab_tg_fu_value  := lab_tg_base_value * fu_mult * rlnorm(.N, 0, 0.30)]
sim_p[, lab_hdl_fu_value := pmax(20, lab_hdl_base_value + rnorm(.N, 0, 8))]
sim_p[, lab_ldl_fu_value := pmax(40, lab_ldl_base_value * (1 + rnorm(.N, 0, 0.10)))]

qtl <- function(x) quantile(x, c(0.25, 0.5, 0.75), na.rm = TRUE)
mk_med <- function(an, base_col, fu_col) {
  rbindlist(lapply(c("oral", "transdermal"), function(a) {
    sb <- sim_p[arm == a]
    qb <- qtl(sb[[base_col]]); qf <- qtl(sb[[fu_col]])
    data.table(
      analyte = an, arm = a, time = c("baseline", "follow-up"),
      median = c(qb[2], qf[2]), q1 = c(qb[1], qf[1]), q3 = c(qb[3], qf[3]))
  }))
}
med_fu <- rbind(
  mk_med("TG",  "lab_tg_base_value", "lab_tg_fu_value"),
  mk_med("HDL", "lab_hdl_base_value", "lab_hdl_fu_value"),
  mk_med("LDL", "lab_ldl_base_value", "lab_ldl_fu_value"))

# ---- Simulated imputation chain means (trace plots) ----
mice_pub <- c(vital_bmi_base_value = "BMI",
              vital_sbp_base_value = "Systolic BP",
              log_tg_base          = "log TG baseline")
mice_imp_vars <- names(mice_pub)
mice_pct_miss <- sapply(mice_imp_vars, function(v) 100 * mean(is.na(sim[[v]])))
mice_lab <- sprintf("%s (%.1f%% missing)", unname(mice_pub), mice_pct_miss)

# Geometric-decay chains converging to the imputation target, per chain
mk_chain <- function(target, start) {
  sapply(seq_len(M_IMP), function(j)
    target + (start - target) * exp(-(seq_len(N_ITER) - 1) / 5) +
      rnorm(N_ITER, 0, abs(start - target) / 12 + 0.002))
}
chain_targets <- c(vital_bmi_base_value = bmi_mu,
                   vital_sbp_base_value = sbp_mu,
                   log_tg_base          = ltg_mu)
chain_starts  <- chain_targets + c(0.8, 2.0, 0.05)
trace_list <- list()
for (v in mice_imp_vars) {
  mat <- mk_chain(chain_targets[v], chain_starts[v])  # N_ITER x M_IMP
  for (j in seq_len(ncol(mat))) {
    trace_list[[length(trace_list) + 1]] <- data.table(
      variable = v, iteration = seq_len(nrow(mat)),
      chain = j, mean_val = mat[, j])
  }
}
trace_dt <- rbindlist(trace_list)
trace_dt[, variable_f := factor(variable, levels = mice_imp_vars, labels = mice_lab)]

# Observed vs imputed values (density overlays), imputation 1
dens_list <- list()
for (v in mice_imp_vars) {
  obs <- sim[[v]][!is.na(sim[[v]])]
  imp <- fr_p[[1]][[v]][is.na(sim[[v]])]
  dens_list[[length(dens_list) + 1]] <- data.table(variable = v, value = obs, type = "Observed")
  dens_list[[length(dens_list) + 1]] <- data.table(variable = v, value = imp, type = "Imputed")
}
dens_dt <- rbindlist(dens_list)
dens_dt[, variable_f := factor(variable, levels = mice_imp_vars, labels = mice_lab)]
dens_dt[, type_f := factor(type, levels = c("Observed", "Imputed"))]

# ---- Simulated analysis results (forest plots) ----
# Random effects on the percentage-difference scale, purely illustrative
mk_est <- function(n_rows = 1) {
  est <- round(runif(n_rows, -4, 16), 1)
  hw  <- round(runif(n_rows, 3, 6), 1)
  data.table(est = est, ci_lo = est - hw, ci_hi = est + hw,
             p = round(2 * pnorm(-abs(est) / (hw / 1.96)), 3))
}
# Simulated per-analysis sample counts (fractions of the outcome sample)
cnt <- function(frac_oral, frac_tdm) {
  c(oral = round(frac_oral * n_oral), tdm = round(frac_tdm * n_tdm))
}
m2_oral   <- cnt(0.55, 1)          # oral formulation subgroups vs full transdermal
m2_e2d    <- cnt(1 - 0.55, 1)
dose_low  <- cnt(0.52, 0.47)
dose_std  <- cnt(0.48, 0.53)
pp2       <- cnt(0.55, 0.58)       # on-treatment subset
pp3       <- cnt(0.48, 0.56)       # on-treatment, no new lipid-lowering
m1_cnt    <- cnt(0.70, 0.87)       # estrogen-only regimen
m1b_cnt   <- cnt(0.40, 0.53)       # on-treatment, estrogen-only

res_rows <- list(
  P1_TG                     = list(n = cnt(1, 1)),
  B4_TG                     = list(n = pp2),
  P2_TG                     = list(n = pp2),
  P3_TG                     = list(n = pp3),
  M2_cnj_vs_tdm_TG          = list(n = m2_oral),
  M2_e2d_vs_tdm_TG          = list(n = m2_e2d),
  Dose_low_TG               = list(n = dose_low),
  Dose_standard_TG          = list(n = dose_std),
  S1_TG                     = list(n = cnt(1, 1)),
  SENS_rawbase_P1_TG        = list(n = cnt(1, 1)),
  M1_TG                     = list(n = m1_cnt),
  M1b_TG                    = list(n = m1b_cnt),
  B4P3_TG                   = list(n = pp3),
  SENS_persist_noarm_B4_TG  = list(n = pp2),
  SENS_persist_noarm_B4P3_TG = list(n = pp3),
  M2_cnj_vs_e2d_TG          = list(n = m2_oral)
)
res <- rbindlist(lapply(names(res_rows), function(a) {
  e <- mk_est()
  data.table(analysis = a, est = e$est, ci_lo = e$ci_lo, ci_hi = e$ci_hi, p = e$p,
             n = sum(res_rows[[a]]$n))
}))
# Alternative-weighting variants of the primary analysis (same simulated sample)
p1_v2 <- mk_est(); p1_v3 <- mk_est()
res[analysis == "P1_TG", `:=`(est_v2 = p1_v2$est, ci_lo_v2 = p1_v2$ci_lo,
                              ci_hi_v2 = p1_v2$ci_hi, p_v2 = p1_v2$p,
                              est_v3 = p1_v3$est, ci_lo_v3 = p1_v3$ci_lo,
                              ci_hi_v3 = p1_v3$ci_hi, p_v3 = p1_v3$p)]

# ===================== figure_participant_flow (140 mm) =====================
# Vertical main column with side exclusion boxes; grey fill for exclusions
stopifnot(n_base_tg == n_all - n_no_base)
stopifnot(n_base_tg == n_base_oral + n_base_tdm)
stopifnot(n_oral + n_fu_oral == n_base_oral)
stopifnot(n_tdm  + n_fu_tdm  == n_base_tdm)
# Disclosure-avoidance check: no small counts anywhere in the diagram
stopifnot(all(c(n_all, n_no_base, n_base_tg, n_base_oral, n_base_tdm,
                n_fu_oral, n_fu_tdm, n_oral, n_tdm) > 20))

fmt_n <- function(n) format(n, big.mark = ",")
flow_boxes <- data.table(
  id = c("all", "base",
         "arm_oral", "arm_tdm",
         "an_oral", "an_tdm",
         "no_base", "fu_oral", "fu_tdm"),
  label = c(
    sprintf("Women initiating oral or transdermal\nestrogen who met eligibility criteria\nN = %s", fmt_n(n_all)),
    sprintf("Baseline TG available\nN = %s", fmt_n(n_base_tg)),
    sprintf("Oral\nn = %s", fmt_n(n_base_oral)),
    sprintf("Transdermal\nn = %s", fmt_n(n_base_tdm)),
    sprintf("Analysed\n(baseline and\nfollow-up TG)\nn = %s", fmt_n(n_oral)),
    sprintf("Analysed\n(baseline and\nfollow-up TG)\nn = %s", fmt_n(n_tdm)),
    sprintf("No valid TG measurement\nwithin 365 days before\nthe index date\n(n = %s)", fmt_n(n_no_base)),
    sprintf("No follow-up TG\n(days 30\u2013365)\nn = %s", fmt_n(n_fu_oral)),
    sprintf("No follow-up TG\n(days 30\u2013365)\nn = %s", fmt_n(n_fu_tdm))
  ),
  x = c(0.50, 0.50,
        0.30, 0.70,
        0.30, 0.70,
        0.89, 0.08, 0.92),
  y = c(0.91, 0.72,
        0.52, 0.52,
        0.17, 0.17,
        0.805, 0.35, 0.35),
  w = c(0.23, 0.17,
        0.13, 0.13,
        0.13, 0.13,
        0.14, 0.08, 0.08),
  h = c(0.075, 0.055,
        0.055, 0.055,
        0.065, 0.065,
        0.09, 0.055, 0.055),
  is_main = c(TRUE, TRUE,
              TRUE, TRUE,
              TRUE, TRUE,
              FALSE, FALSE, FALSE)
)

draw_flow <- function() {
  ggplot() +
    # main boxes (white) and exclusion boxes (grey95)
    geom_rect(data = flow_boxes[is_main == TRUE],
              aes(xmin = x - w, xmax = x + w, ymin = y - h, ymax = y + h),
              fill = "white", colour = "black", linewidth = LINEW / .pt) +
    geom_rect(data = flow_boxes[is_main == FALSE],
              aes(xmin = x - w, xmax = x + w, ymin = y - h, ymax = y + h),
              fill = "grey95", colour = "black", linewidth = LINEW / .pt) +
    geom_text(data = flow_boxes, aes(x = x, y = y, label = label),
              family = FONT, size = 5.5 / .pt, lineheight = 1.0) +
    # all -> base (down)
    geom_segment(aes(x = 0.5, y = 0.91 - 0.075, xend = 0.5, yend = 0.72 + 0.055),
                 arrow = arrow(length = unit(1.5, "mm")), linewidth = LINEW / .pt) +
    # all->base vertical -> no_base exclusion (right, head at exclusion box)
    geom_segment(aes(x = 0.5, y = 0.805, xend = 0.89 - 0.14, yend = 0.805),
                 arrow = arrow(length = unit(1.5, "mm")), linewidth = LINEW / .pt) +
    # base -> split elbow
    geom_segment(aes(x = 0.5, y = 0.72 - 0.055, xend = 0.5, yend = 0.615),
                 linewidth = LINEW / .pt) +
    geom_segment(aes(x = 0.5, y = 0.615, xend = 0.30, yend = 0.615),
                 linewidth = LINEW / .pt) +
    geom_segment(aes(x = 0.5, y = 0.615, xend = 0.70, yend = 0.615),
                 linewidth = LINEW / .pt) +
    # split -> arms (down, head at arm box)
    geom_segment(aes(x = 0.30, y = 0.615, xend = 0.30, yend = 0.52 + 0.055),
                 arrow = arrow(length = unit(1.5, "mm")), linewidth = LINEW / .pt) +
    geom_segment(aes(x = 0.70, y = 0.615, xend = 0.70, yend = 0.52 + 0.055),
                 arrow = arrow(length = unit(1.5, "mm")), linewidth = LINEW / .pt) +
    # arms -> analysed (down, head at analysed box)
    geom_segment(aes(x = 0.30, y = 0.52 - 0.055, xend = 0.30, yend = 0.17 + 0.065),
                 arrow = arrow(length = unit(1.5, "mm")), linewidth = LINEW / .pt) +
    geom_segment(aes(x = 0.70, y = 0.52 - 0.055, xend = 0.70, yend = 0.17 + 0.065),
                 arrow = arrow(length = unit(1.5, "mm")), linewidth = LINEW / .pt) +
    # oral vertical -> fu_oral exclusion (left, head at exclusion box)
    geom_segment(aes(x = 0.30, y = 0.35, xend = 0.08 + 0.08, yend = 0.35),
                 arrow = arrow(length = unit(1.5, "mm")), linewidth = LINEW / .pt) +
    # tdm vertical -> fu_tdm exclusion (right, head at exclusion box)
    geom_segment(aes(x = 0.70, y = 0.35, xend = 0.92 - 0.08, yend = 0.35),
                 arrow = arrow(length = unit(1.5, "mm")), linewidth = LINEW / .pt) +
    xlim(-0.03, 1.05) + ylim(0.06, 1.02) +
    theme_void(base_family = FONT) +
    theme(plot.margin = margin(2, 2, 2, 2))
}

p_flow <- draw_flow()
save_pdf(p_flow, "figure_participant_flow", W_1_5, 130)
save_tiff(p_flow, "figure_participant_flow", W_1_5, 130, dpi = 1000)

# ===================== figure_forest_primary (140 mm) =====================
# Groups: ITT, per-protocol, on-treatment, by formulation, by dose, sensitivity
# X axis: log-ratio scale labelled as percentages, reference line at 0%

fmt_pct <- function(x) {  # minus sign uses U+2212
  neg <- x < 0
  s <- sprintf("%.1f", abs(x))
  ifelse(neg, paste0("\u2212", s), s)
}

fp_def <- data.table(
  analysis = c("P1_TG", "B4_TG", "P2_TG", "P3_TG",
               "M2_cnj_vs_tdm_TG", "M2_e2d_vs_tdm_TG",
               "Dose_low_TG", "Dose_standard_TG",
               "S1_TG", "SENS_rawbase_P1_TG", "P1_TG_v2", "P1_TG_v3"),
  label = c("Primary analysis (selection-weighted)",
            "Persistence-weighted",
            "On-treatment at follow-up",
            "On-treatment, no new lipid-lowering",
            "Conjugated estrogens",
            "Oral estradiol",
            "Low dose (<1 mg/day estradiol-eq)",
            "Standard dose (\u22651 mg/day estradiol-eq)",
            "Parsimonious propensity score",
            "Original-scale baseline TG",
            "Overlap weights",
            "CBPS weights only (no selection weighting)"),
  group = c("Intention-to-treat", "Per-protocol",
            "On-treatment subsets", "On-treatment subsets",
            "By oral formulation (vs transdermal)", "By oral formulation (vs transdermal)",
            "By dose stratum", "By dose stratum",
            "Sensitivity analysis", "Sensitivity analysis",
            "Sensitivity analysis", "Sensitivity analysis"),
  n_oral = c(n_oral, pp2[["oral"]], pp2[["oral"]], pp3[["oral"]],
             m2_oral[["oral"]], m2_e2d[["oral"]],
             dose_low[["oral"]], dose_std[["oral"]],
             n_oral, n_oral, n_oral, n_oral),
  n_tdm  = c(n_tdm, pp2[["tdm"]], pp2[["tdm"]], pp3[["tdm"]],
             m2_oral[["tdm"]], m2_e2d[["tdm"]],
             dose_low[["tdm"]], dose_std[["tdm"]],
             n_tdm, n_tdm, n_tdm, n_tdm)
)

# Point estimates and intervals from the simulated results table
fp_data <- merge(fp_def[!grepl("^P1_TG_v[23]$", analysis)],
                 res[, .(analysis, est, ci_lo, ci_hi, p, n)], by = "analysis")
# Alternative-weighting variants (columns *_v2 / *_v3 of the primary row)
p1_var <- res[analysis == "P1_TG"]
fp_v2v3 <- data.table(
  analysis = c("P1_TG_v2", "P1_TG_v3"),
  est    = unlist(p1_var[, .(est_v2, est_v3)]),
  ci_lo  = unlist(p1_var[, .(ci_lo_v2, ci_lo_v3)]),
  ci_hi  = unlist(p1_var[, .(ci_hi_v2, ci_hi_v3)]),
  p      = unlist(p1_var[, .(p_v2, p_v3)]),
  n      = c(n_oral + n_tdm, n_oral + n_tdm)
)
fp_data <- rbind(fp_data, merge(fp_def[grepl("^P1_TG_v[23]$", analysis)],
                                fp_v2v3, by = "analysis"))
stopifnot(nrow(fp_data) == 12)

# Rows with group headers stacked on top
group_order <- c("Intention-to-treat", "Per-protocol", "On-treatment subsets",
                 "By oral formulation (vs transdermal)", "By dose stratum",
                 "Sensitivity analysis")
fp_rows <- list()
y_pos <- nrow(fp_data) + length(group_order)  # start from the top
for (g in group_order) {
  fp_rows[[length(fp_rows) + 1]] <- data.table(
    is_header = TRUE, label = g, y = y_pos,
    est = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
    n_oral = NA_integer_, n_tdm = NA_integer_
  )
  y_pos <- y_pos - 1
  sub <- fp_data[group == g]
  for (i in seq_len(nrow(sub))) {
    r <- sub[i]
    fp_rows[[length(fp_rows) + 1]] <- data.table(
      is_header = FALSE, label = r$label, y = y_pos,
      est = r$est, ci_lo = r$ci_lo, ci_hi = r$ci_hi,
      n_oral = r$n_oral, n_tdm = r$n_tdm
    )
    y_pos <- y_pos - 1
  }
}
fp_plot <- rbindlist(fp_rows)

# Ratio conversion: ratio = 1 + est / 100
fp_plot[, ratio := 1 + est / 100]
fp_plot[, ratio_lo := 1 + ci_lo / 100]
fp_plot[, ratio_hi := 1 + ci_hi / 100]
fp_plot[is_header == FALSE, n_txt := sprintf("%d / %d", n_oral, n_tdm)]
fp_plot[is_header == FALSE, est_txt := paste0(fmt_pct(est), "%")]
fp_plot[is_header == FALSE, ci_txt := sprintf("(%s%%, %s%%)", fmt_pct(ci_lo), fmt_pct(ci_hi))]

draw_forest_primary <- function() {
  # Column proportions: label 0.38 | N 0.10 | plot 0.32 | estimate 0.20
  # X range -10% to 35%; intervals must not enter the estimate column
  r_label <- 0.38; r_n <- 0.10; r_plot <- 0.32; r_est <- 0.20
  plot_x_lo <- 0.9; plot_x_hi <- 1.35  # -10% to 35%
  log_plot_lo <- log10(plot_x_lo); log_plot_hi <- log10(plot_x_hi)
  log_plot_range <- log_plot_hi - log_plot_lo
  log_total_range <- log_plot_range / r_plot
  log_min <- log_plot_lo - (r_label + r_n) * log_total_range
  log_max <- log_plot_hi + r_est * log_total_range
  x_min <- 10^log_min; x_max <- 10^log_max
  x_lab <- x_min
  x_n <- 10^(log_min + (r_label + r_n / 2) * log_total_range)
  x_est <- x_max * 0.98
  # Direction annotations sit 2 mm either side of the reference line
  # (panel width about 136 mm = 140 - 2 x 2 mm margins)
  log_gap <- (2 / 136) * (log_max - log_min)
  x_gap_l <- 10^(log10(1) - log_gap)  # left annotation right edge
  x_gap_r <- 10^(log10(1) + log_gap)  # right annotation left edge

  ggplot(fp_plot, aes(y = y)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = LINEW / .pt) +
    # data-row intervals
    geom_segment(data = fp_plot[is_header == FALSE],
                 aes(x = ratio_lo, xend = ratio_hi, yend = y),
                 linewidth = LINEW / .pt, colour = "black") +
    geom_point(data = fp_plot[is_header == FALSE],
               aes(x = ratio), size = 2, shape = 16, colour = "black") +
    # group headers (left-aligned, bold)
    geom_text(data = fp_plot[is_header == TRUE],
              aes(x = x_lab, label = label), hjust = 0,
              family = FONT, size = 6.5 / .pt, fontface = "bold") +
    # row labels
    geom_text(data = fp_plot[is_header == FALSE],
              aes(x = x_lab, label = label), hjust = 0,
              family = FONT, size = 6 / .pt) +
    # N column
    geom_text(data = fp_plot[is_header == FALSE],
              aes(x = x_n, label = n_txt), hjust = 0.5,
              family = FONT, size = 6 / .pt) +
    # estimates (right-aligned, right of the plot area)
    geom_text(data = fp_plot[is_header == FALSE],
              aes(x = x_est, label = paste0(est_txt, " ", ci_txt)), hjust = 1,
              family = FONT, size = 6 / .pt) +
    scale_x_log10(limits = c(x_min, x_max),
                  breaks = c(0.9, 1, 1.1, 1.2, 1.3),
                  labels = c("\u221210%", "0%", "10%", "20%", "30%")) +
    scale_y_continuous(limits = c(0, max(fp_plot$y) + 1), expand = c(0, 0)) +
    labs(x = "Difference in follow-up TG, oral vs transdermal (%)", y = NULL) +
    # direction annotations: 2 mm gap to the reference line, 7 pt, grey
    annotate("text", x = x_gap_l, y = max(fp_plot$y) + 0.7,
             label = "\u2190 Lower with oral", hjust = 1,
             family = FONT, size = 7 / .pt, colour = "grey40") +
    annotate("text", x = x_gap_r, y = max(fp_plot$y) + 0.7,
             label = "Higher with oral \u2192", hjust = 0,
             family = FONT, size = 7 / .pt, colour = "grey40") +
    theme_pub() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25 / .pt),
      plot.margin = margin(2, 2, 2, 2)
    )
}

p_forest <- draw_forest_primary()
save_pdf(p_forest, "figure_forest_primary", W_1_5, 135)
save_tiff(p_forest, "figure_forest_primary", W_1_5, 135, dpi = 1000)

# ===================== figure_forest_full (140 mm) =====================
# Extended forest plot: additional sensitivity analyses

ff_def <- data.table(
  analysis = c("M1_TG", "M1b_TG",
               "B4P3_TG", "SENS_persist_noarm_B4_TG", "SENS_persist_noarm_B4P3_TG",
               "M2_cnj_vs_e2d_TG"),
  label = c("Estrogen only (no progestogen)",
            "On-treatment, estrogen only",
            "No new lipid-lowering, persistence-weighted",
            "Persistence model not fitted by arm",
            "No new lipid-lowering, persistence model not fitted by arm",
            "Conjugated estrogens vs oral estradiol"),
  group = c("Oral vs transdermal: modified regimen", "Oral vs transdermal: modified regimen",
            "Oral vs transdermal: per-protocol variants", "Oral vs transdermal: per-protocol variants", "Oral vs transdermal: per-protocol variants",
            "Between oral formulations"),
  n_oral = c(m1_cnt[["oral"]], m1b_cnt[["oral"]], pp3[["oral"]], pp2[["oral"]],
             pp3[["oral"]], m2_oral[["oral"]]),
  n_tdm  = c(m1_cnt[["tdm"]], m1b_cnt[["tdm"]], pp3[["tdm"]], pp2[["tdm"]],
             pp3[["tdm"]], m2_e2d[["oral"]]),
  n_total = c(sum(m1_cnt), sum(m1b_cnt), sum(pp3), sum(pp2), sum(pp3),
              m2_oral[["oral"]] + m2_e2d[["oral"]])
)
stopifnot(all(ff_def$analysis %in% res$analysis))
ff_data <- merge(ff_def, res[, .(analysis, est, ci_lo, ci_hi, p, n)], by = "analysis")
stopifnot(nrow(ff_data) == 6)

ff_group_order <- c("Oral vs transdermal: modified regimen",
                    "Oral vs transdermal: per-protocol variants",
                    "Between oral formulations")
ff_rows <- list()
y_pos <- nrow(ff_data) + length(ff_group_order)
for (g in ff_group_order) {
  ff_rows[[length(ff_rows) + 1]] <- data.table(
    is_header = TRUE, label = g, y = y_pos,
    est = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_,
    n_oral = NA_integer_, n_tdm = NA_integer_, n_total = NA_integer_
  )
  y_pos <- y_pos - 1
  sub <- ff_data[group == g]
  for (i in seq_len(nrow(sub))) {
    r <- sub[i]
    ff_rows[[length(ff_rows) + 1]] <- data.table(
      is_header = FALSE, label = r$label, y = y_pos,
      est = r$est, ci_lo = r$ci_lo, ci_hi = r$ci_hi,
      n_oral = r$n_oral, n_tdm = r$n_tdm, n_total = r$n_total
    )
    y_pos <- y_pos - 1
  }
}
ff_plot <- rbindlist(ff_rows)

ff_plot[, ratio := 1 + est / 100]
ff_plot[, ratio_lo := 1 + ci_lo / 100]
ff_plot[, ratio_hi := 1 + ci_hi / 100]
ff_plot[is_header == FALSE, est_txt := paste0(fmt_pct(est), "%")]
ff_plot[is_header == FALSE, ci_txt := sprintf("(%s%%, %s%%)", fmt_pct(ci_lo), fmt_pct(ci_hi))]
ff_plot[is_header == FALSE & !is.na(n_oral), n_txt := sprintf("%d / %d", n_oral, n_tdm)]
ff_plot[is_header == FALSE & is.na(n_oral), n_txt := sprintf("N = %d", n_total)]

draw_forest_full <- function() {
  x_lab <- 0.16; x_n <- 0.72; x_est <- 2.00; x_max <- 2.15
  ggplot(ff_plot, aes(y = y)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = LINEW / .pt) +
    geom_segment(data = ff_plot[is_header == FALSE],
                 aes(x = ratio_lo, xend = ratio_hi, yend = y),
                 linewidth = LINEW / .pt, colour = "black") +
    geom_point(data = ff_plot[is_header == FALSE],
               aes(x = ratio), size = 2, shape = 16, colour = "black") +
    geom_text(data = ff_plot[is_header == TRUE],
              aes(x = x_lab, label = label), hjust = 0,
              family = FONT, size = 6.5 / .pt, fontface = "bold") +
    geom_text(data = ff_plot[is_header == FALSE],
              aes(x = x_lab, label = label), hjust = 0,
              family = FONT, size = 6 / .pt) +
    geom_text(data = ff_plot[is_header == FALSE],
              aes(x = x_n, label = n_txt), hjust = 0.5,
              family = FONT, size = 6 / .pt) +
    geom_text(data = ff_plot[is_header == FALSE],
              aes(x = x_est, label = paste0(est_txt, " ", ci_txt)), hjust = 1,
              family = FONT, size = 6 / .pt) +
    scale_x_log10(limits = c(0.12, x_max),
                  breaks = c(0.9, 1, 1.1, 1.2, 1.3),
                  labels = c("\u221210%", "0%", "10%", "20%", "30%")) +
    scale_y_continuous(limits = c(0, max(ff_plot$y) + 1), expand = c(0, 0)) +
    labs(x = "Difference in follow-up TG (%)", y = NULL) +
    theme_pub() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25 / .pt),
          plot.margin = margin(2, 2, 2, 2))
}

p_forest_full <- draw_forest_full()
save_pdf(p_forest_full, "figure_forest_full", W_1_5, 120)
save_tiff(p_forest_full, "figure_forest_full", W_1_5, 120, dpi = 1000)

# ===================== figure_tg_baseline_followup (190 x 70 mm) =====================
# Three panels: TG | HDL | LDL, median (IQR); oral solid circle,
# transdermal hollow triangle
tg_data <- copy(med_fu)
tg_data[, time_f := factor(time, levels = c("baseline", "follow-up"),
                           labels = c("Baseline", "Follow-up"))]
tg_data[, arm_f := factor(arm, levels = c("oral", "transdermal"),
                          labels = c("Oral", "Transdermal"))]
tg_data[, x_n := ifelse(arm == "oral", -0.15, 0.15)]  # dodge left/right
tg_data[, x_pos := as.numeric(time_f) + x_n]

panel_list <- list()
for (an in c("TG", "HDL", "LDL")) {
  dd <- tg_data[analyte == an]
  y_lab <- if (an == "TG") "Median (IQR), mg/dL" else NULL
  p <- ggplot(dd, aes(x = x_pos, y = median, colour = arm_f, shape = arm_f)) +
    geom_errorbar(aes(ymin = q1, ymax = q3), width = 0, linewidth = LINEW / .pt) +
    geom_point(size = 2.2) +
    geom_line(aes(group = arm_f), linewidth = LINEW / .pt) +
    scale_colour_manual(values = c(Oral = unname(OKABE_ITO["oral"]),
                                   Transdermal = unname(OKABE_ITO["tdm"]))) +
    scale_shape_manual(values = c(Oral = 16, Transdermal = 2)) +
    scale_x_continuous(breaks = c(1, 2), labels = c("Baseline", "Follow-up"),
                       limits = c(0.5, 2.5)) +
    labs(x = NULL, y = y_lab, title = an) +
    theme_pub() +
    theme(plot.title = element_text(hjust = 0.5, size = 7, face = "bold"))
  if (an == "TG") {
    tg_pan <- tg_data[analyte == "TG"]
    p <- p + scale_y_log10(breaks = c(75, 100, 150, 200),
                           limits = c(min(tg_pan$q1) * 0.85, max(tg_pan$q3) * 1.08))
  }
  panel_list[[an]] <- p
}
p_tg_fu <- wrap_plots(panel_list, nrow = 1, guides = "collect") &
  theme(legend.position = "bottom",
        legend.key.height = unit(4, "mm"),
        legend.text = element_text(size = 7),
        legend.title = element_blank())
save_pdf(p_tg_fu, "figure_tg_baseline_followup", W_DOUBLE, 70)
save_tiff(p_tg_fu, "figure_tg_baseline_followup", W_DOUBLE, 70, dpi = 1000)

# ===================== figure_balance_diagnostics (190 x 140 mm) =====================
# A: PS overlap density (with legend); C: combined-weight boxplots (top row)
# B: Love plot, unweighted vs weighted SMD (bottom row, full width)

# --- A. PS overlap ---
ps_df <- data.table(subject_id = pid_wt, ps = ps_main)
ps_df <- merge(ps_df, arm_map, by = "subject_id")
ps_df[, arm_f := factor(arm, levels = c("oral", "transdermal"),
                        labels = c("Oral", "Transdermal"))]

p_ps_overlap <- ggplot(ps_df, aes(x = ps, colour = arm_f, fill = arm_f)) +
  geom_density(alpha = 0.3, linewidth = LINEW / .pt) +
  scale_colour_manual(values = c(Oral = unname(OKABE_ITO["oral"]),
                                 Transdermal = unname(OKABE_ITO["tdm"]))) +
  scale_fill_manual(values = c(Oral = unname(OKABE_ITO["oral"]),
                               Transdermal = unname(OKABE_ITO["tdm"]))) +
  labs(x = "Propensity score", y = "Density", tag = "A") +
  theme_pub() +
  theme(legend.position = c(0.6, 0.8),
        legend.title = element_blank(),
        legend.key.height = unit(3, "mm"),
        legend.text = element_text(size = 6),
        legend.background = element_blank())

# --- B. Love plot: 11 PS-model covariates plus HDL/LDL baseline ---
love_covars <- c("age_at_index", "index_year",
                 "vital_bmi_base_value", "vital_sbp_base_value",
                 "cx_vte_ever", "cx_chd_ever", "cx_htn_base", "cx_t2dm_base",
                 "rx_statin_base", "rx_progestogen_base", "log_tg_base")
love_dt <- rbindlist(lapply(love_covars, function(v) {
  data.table(cov = v,
             smd_unw = smd_imp_mean(v, weighted = FALSE),
             smd_wt  = smd_imp_mean(v, weighted = TRUE))
}))
# HDL/LDL baseline: not in the PS model, computed directly from the cohort
love_dt <- rbind(love_dt, rbindlist(lapply(c("lab_hdl_base_value", "lab_ldl_base_value"),
  function(v) {
    data.table(cov = v,
               smd_unw = calc_smd_fig(sim[[v]][ip_vec], arm_num_p, rep(1, sum(ip_vec))),
               smd_wt  = mean(sapply(seq_len(M_IMP), function(b)
                 calc_smd_fig(sim[[v]][ip_vec], arm_num_p, wC_main[, b]))))
  })))
stopifnot(nrow(love_dt) == 13)

pub_name <- c(
  "age_at_index" = "Age at index",
  "index_year" = "Index year",
  "vital_bmi_base_value" = "BMI",
  "vital_sbp_base_value" = "Systolic BP",
  "cx_vte_ever" = "Venous thromboembolism",
  "cx_chd_ever" = "Coronary heart disease",
  "cx_htn_base" = "Hypertension",
  "cx_t2dm_base" = "Type 2 diabetes",
  "rx_statin_base" = "Statin",
  "rx_progestogen_base" = "Progestogen",
  "log_tg_base" = "Triglycerides (log)",
  "lab_hdl_base_value" = "HDL cholesterol",
  "lab_ldl_base_value" = "LDL cholesterol"
)
love_dt[, pub := pub_name[cov]]
love_dt <- love_dt[order(smd_unw, decreasing = TRUE)]
love_dt[, pub_f := factor(pub, levels = pub)]

love_long <- melt(love_dt, id.vars = c("cov", "pub", "pub_f"),
                  measure.vars = c("smd_unw", "smd_wt"),
                  variable.name = "weighting", value.name = "smd")
love_long[, weighting := factor(weighting,
                                levels = c("smd_unw", "smd_wt"),
                                labels = c("Unweighted", "Weighted"))]

p_love <- ggplot(love_long, aes(x = smd, y = pub_f, colour = weighting, shape = weighting)) +
  geom_vline(xintercept = 0.1, linetype = "dashed", colour = "grey50", linewidth = LINEW / .pt) +
  geom_vline(xintercept = 0, colour = "black", linewidth = LINEW / .pt) +
  geom_point(size = 1.8) +
  geom_line(aes(group = pub_f), colour = "grey70", linewidth = 0.2 / .pt) +
  scale_colour_manual(values = c("Unweighted" = "grey40", "Weighted" = "black")) +
  scale_shape_manual(values = c("Unweighted" = 1, "Weighted" = 16)) +
  labs(x = "Absolute SMD", y = NULL, tag = "B") +
  scale_x_continuous(limits = c(0, max(love_long$smd, na.rm = TRUE) * 1.05)) +
  theme_pub() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        axis.text.y = element_text(size = 6),
        legend.key.height = unit(3, "mm"))

# --- C. Weight distribution (first imputation) ---
w_df <- data.table(subject_id = pid_wt, weight = as.numeric(wC_main[, 1]))
w_df <- merge(w_df, arm_map, by = "subject_id")
w_df[, arm_f := factor(arm, levels = c("oral", "transdermal"))]

p_wt <- ggplot(w_df, aes(x = arm_f, y = weight, fill = arm_f)) +
  geom_boxplot(outlier.shape = NA, linewidth = LINEW / .pt, width = 0.5, alpha = 0.3) +
  scale_fill_manual(values = c(oral = unname(OKABE_ITO["oral"]),
                               transdermal = unname(OKABE_ITO["tdm"]))) +
  scale_x_discrete(labels = c("Oral", "Transdermal")) +
  labs(x = NULL, y = "Combined weight", tag = "C") +
  theme_pub() +
  theme(legend.position = "none")

# Layout: A + C top row, B bottom row (full width)
p_balance <- wrap_plots(
  wrap_plots(p_ps_overlap, p_wt, nrow = 1),
  p_love,
  ncol = 1, heights = c(1, 1.3)
) & theme(plot.tag = element_text(size = 9, face = "bold", family = FONT))
save_pdf(p_balance, "figure_balance_diagnostics", W_DOUBLE, 140)
save_tiff(p_balance, "figure_balance_diagnostics", W_DOUBLE, 140, dpi = 1000)

stopifnot(nrow(ps_df) == sum(ip_vec))
stopifnot(nrow(w_df) == sum(ip_vec))

# ===================== figure_mice_diagnostics (190 x 90 mm) =====================
# A: chain means by iteration (trace); B: observed vs imputed density
p_trace <- ggplot(trace_dt, aes(x = iteration, y = mean_val, colour = factor(chain))) +
  geom_line(linewidth = 0.2 / .pt) +
  facet_wrap(~ variable_f, scales = "free_y", ncol = 3) +
  scale_colour_grey(start = 0.3, end = 0.7) +
  labs(x = "Iteration", y = "Mean", tag = "A") +
  theme_pub() +
  theme(legend.position = "none",
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 5),
        panel.spacing = unit(1.5, "mm"))

p_dens <- ggplot(dens_dt, aes(x = value, colour = type_f, fill = type_f)) +
  geom_density(alpha = 0.2, linewidth = LINEW / .pt) +
  facet_wrap(~ variable_f, scales = "free", ncol = 3) +
  scale_colour_manual(values = c("Observed" = "black", "Imputed" = unname(OKABE_ITO["oral"]))) +
  scale_fill_manual(values = c("Observed" = "black", "Imputed" = unname(OKABE_ITO["oral"]))) +
  labs(x = "Value", y = "Density", tag = "B") +
  theme_pub() +
  theme(legend.position = "bottom",
        legend.title = element_blank(),
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 5),
        panel.spacing = unit(1.5, "mm"))

p_mice <- wrap_plots(p_trace, p_dens, ncol = 1, heights = c(1, 1)) &
  theme(plot.tag = element_text(size = 9, face = "bold", family = FONT))
save_pdf(p_mice, "figure_mice_diagnostics", W_DOUBLE, 90)
save_tiff(p_mice, "figure_mice_diagnostics", W_DOUBLE, 90, dpi = 1000)

# ===================== figure_supp_diagnostics (190 x 230 mm) =====================
# Row 1: A PS overlap | B weight distribution (equal width, shared arm legend)
# Row 2: C Love plot (full width); Row 3: D trace (full width);
# Row 4: E observed vs imputed density (full width)
# Legends are collected at the bottom via patchwork guides = "collect"
# --- A: rebuilt with fill layer (no guide) + colour layer (Arm guide) ---
p_A <- ggplot(ps_df, aes(x = ps)) +
  geom_density(aes(fill = arm_f), alpha = 0.3, colour = NA, linewidth = LINEW / .pt) +
  geom_density(aes(colour = arm_f), fill = NA, linewidth = LINEW / .pt) +
  scale_fill_manual(values = c(Oral = unname(OKABE_ITO["oral"]),
                               Transdermal = unname(OKABE_ITO["tdm"])), guide = "none") +
  scale_colour_manual(name = "Arm",
                      values = c(Oral = unname(OKABE_ITO["oral"]),
                                 Transdermal = unname(OKABE_ITO["tdm"])),
                      guide = guide_legend(override.aes = list(
                        fill = unname(OKABE_ITO), alpha = 0.3))) +
  labs(x = "Propensity score", y = "Density", tag = "A") +
  theme_pub() + theme(legend.position = "right")
p_B <- p_wt + labs(tag = "B") + theme(legend.position = "right") +
  scale_fill_manual(name = "Arm",
                    values = c(oral = unname(OKABE_ITO["oral"]),
                               transdermal = unname(OKABE_ITO["tdm"])), guide = "none")
p_C <- p_love + labs(tag = "C") +
  scale_colour_manual(values = c("Unweighted" = "grey40", "Weighted" = "black"), guide = "none") +
  scale_shape_manual(name = "Weighting", values = c("Unweighted" = 1, "Weighted" = 16))
# --- D/E panels: drop the log TG baseline column, keep BMI and Systolic BP ---
supp_vars <- setdiff(mice_imp_vars, "log_tg_base")
stopifnot(setequal(supp_vars, c("vital_bmi_base_value", "vital_sbp_base_value")))
supp_lab <- sapply(supp_vars, function(v)
  sprintf("%s (%.1f%% missing)", mice_pub[v], 100 * sum(is.na(sim[[v]])) / nrow(sim)))
trace_dt_s <- trace_dt[variable %in% supp_vars]
trace_dt_s[, variable_f := factor(variable, levels = supp_vars, labels = unname(supp_lab))]
p_D <- ggplot(trace_dt_s, aes(x = iteration, y = mean_val, colour = factor(chain))) +
  geom_line(linewidth = 0.2 / .pt) +
  facet_wrap(~ variable_f, scales = "free_y", ncol = 2) +
  scale_colour_grey(start = 0.3, end = 0.7, guide = "none") +  # keep chain guide out of collected legends
  labs(x = "Iteration", y = "Mean", tag = "D") +
  theme_pub() +
  theme(legend.position = "none",
        strip.text = element_text(size = 6),
        axis.text = element_text(size = 5),
        panel.spacing = unit(1.5, "mm"))
# --- E: rebuilt with fill layer (Type guide) + colour layer (no guide) ---
dens_dt_s <- dens_dt[variable %in% supp_vars]
dens_dt_s[, variable_f := factor(variable, levels = supp_vars, labels = unname(supp_lab))]
p_E <- ggplot(dens_dt_s, aes(x = value)) +
  geom_density(aes(fill = type_f), alpha = 0.2, colour = NA, linewidth = LINEW / .pt) +
  geom_density(aes(colour = type_f), fill = NA, linewidth = LINEW / .pt) +
  scale_fill_manual(name = "Type",
                    values = c("Observed" = "black", "Imputed" = unname(OKABE_ITO["oral"])),
                    guide = guide_legend(override.aes = list(alpha = 1))) +
  scale_colour_manual(values = c("Observed" = "black", "Imputed" = unname(OKABE_ITO["oral"])),
                      guide = "none") +
  facet_wrap(~ variable_f, scales = "free", ncol = 2) +
  labs(x = "Value", y = "Density", tag = "E") +
  theme_pub() +
  theme(legend.position = "bottom", legend.title = element_blank(),
        strip.text = element_text(size = 6), axis.text = element_text(size = 5),
        panel.spacing = unit(1.5, "mm"))

design_diag <- "
AB
CC
DD
EE
"
p_supp_diag <- wrap_plots(
  A = p_A, B = p_B, C = p_C, D = p_D, E = p_E,
  design = design_diag
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.box = "horizontal", legend.title = element_blank(),
        plot.tag = element_text(size = 9, face = "bold", family = FONT))
save_pdf(p_supp_diag, "figure_supp_diagnostics", W_DOUBLE, 230)
save_tiff(p_supp_diag, "figure_supp_diagnostics", W_DOUBLE, 230, dpi = 600)

# Legend audit: collected guides must be exactly Arm / Weighting / Type
# (2 keys each); no chain-level keys may leak from panel D
gt_supp <- patchwork::patchworkGrob(p_supp_diag)
guide_txt <- character(0)
collect_guide_text <- function(g) {
  if (inherits(g, "text")) {
    guide_txt <<- c(guide_txt, g$label)
  } else if (inherits(g, "gtable")) {
    hit_idx <- which(grepl("guide", g$layout$name))
    if (length(hit_idx)) {  # guide-box subtree: collect all text within
      for (i in hit_idx) for (k in g$grobs[i]) collect_guide_text(k)
    } else {
      for (k in g$grobs) collect_guide_text(k)
    }
  } else if (inherits(g, "gTree")) {  # ggplot2 >= 3.5 wraps labels in gTree
    for (k in g$children) collect_guide_text(k)
  }
  invisible()
}
collect_guide_text(gt_supp)
stopifnot(any(guide_txt %in% c("Oral", "Transdermal")))     # Arm legend
stopifnot(any(guide_txt %in% c("Unweighted", "Weighted")))  # Weighting legend
stopifnot(any(guide_txt %in% c("Observed", "Imputed")))     # Type legend
stopifnot(!any(guide_txt %in% as.character(1:N_ITER)))      # no chain keys

# ===================== figure_supp_additional_results (190 x 170 mm) =====================
# Row 1: A TG | B HDL | C LDL (shared arm legend); Row 2: D extended forest (full width)
p_TG  <- panel_list[["TG"]] + labs(tag = "A") + theme(legend.position = "right") +
  scale_colour_manual(name = "Arm",
                      values = c(Oral = unname(OKABE_ITO["oral"]),
                                 Transdermal = unname(OKABE_ITO["tdm"]))) +
  scale_shape_manual(name = "Arm", values = c(Oral = 16, Transdermal = 2))
p_HDL <- panel_list[["HDL"]] + labs(tag = "B") + theme(legend.position = "right") +
  scale_colour_manual(name = "Arm",
                      values = c(Oral = unname(OKABE_ITO["oral"]),
                                 Transdermal = unname(OKABE_ITO["tdm"]))) +
  scale_shape_manual(name = "Arm", values = c(Oral = 16, Transdermal = 2))
p_LDL <- panel_list[["LDL"]] + labs(tag = "C") + theme(legend.position = "right") +
  scale_colour_manual(name = "Arm",
                      values = c(Oral = unname(OKABE_ITO["oral"]),
                                 Transdermal = unname(OKABE_ITO["tdm"]))) +
  scale_shape_manual(name = "Arm", values = c(Oral = 16, Transdermal = 2))

# D forest: same column proportions as the primary forest plot
draw_forest_combined <- function() {
  r_label <- 0.38; r_n <- 0.10; r_plot <- 0.32; r_est <- 0.20
  plot_x_lo <- 0.9; plot_x_hi <- 1.35  # -10% to 35%
  log_plot_lo <- log10(plot_x_lo); log_plot_hi <- log10(plot_x_hi)
  log_plot_range <- log_plot_hi - log_plot_lo
  log_total_range <- log_plot_range / r_plot
  log_min <- log_plot_lo - (r_label + r_n) * log_total_range
  log_max <- log_plot_hi + r_est * log_total_range
  x_min <- 10^log_min; x_max <- 10^log_max
  x_lab <- x_min
  x_n <- 10^(log_min + (r_label + r_n / 2) * log_total_range)
  x_est <- x_max * 0.98
  ggplot(ff_plot, aes(y = y)) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50", linewidth = LINEW / .pt) +
    geom_segment(data = ff_plot[is_header == FALSE],
                 aes(x = ratio_lo, xend = ratio_hi, yend = y),
                 linewidth = LINEW / .pt, colour = "black") +
    geom_point(data = ff_plot[is_header == FALSE],
               aes(x = ratio), size = 2, shape = 16, colour = "black") +
    geom_text(data = ff_plot[is_header == TRUE],
              aes(x = x_lab, label = label), hjust = 0,
              family = FONT, size = 6 / .pt, fontface = "bold") +
    geom_text(data = ff_plot[is_header == FALSE],
              aes(x = x_lab, label = label), hjust = 0,
              family = FONT, size = 5.5 / .pt) +
    geom_text(data = ff_plot[is_header == FALSE],
              aes(x = x_n, label = n_txt), hjust = 0.5,
              family = FONT, size = 5.5 / .pt) +
    geom_text(data = ff_plot[is_header == FALSE],
              aes(x = x_est, label = paste0(est_txt, " ", ci_txt)), hjust = 1,
              family = FONT, size = 5.5 / .pt) +
    scale_x_log10(limits = c(x_min, x_max),
                  breaks = c(0.9, 1, 1.1, 1.2, 1.3),
                  labels = c("\u221210%", "0%", "10%", "20%", "30%")) +
    scale_y_continuous(limits = c(0, max(ff_plot$y) + 1), expand = c(0, 0)) +
    labs(x = "Difference in follow-up TG (%)", y = NULL, tag = "D") +
    theme_pub() +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.25 / .pt),
          plot.margin = margin(2, 2, 2, 2))
}
p_D_forest <- draw_forest_combined()

design_res <- "
ABC
DDD
"
p_supp_res <- wrap_plots(
  A = p_TG, B = p_HDL, C = p_LDL, D = p_D_forest,
  design = design_res
) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.box = "horizontal", legend.title = element_blank(),
        plot.tag = element_text(size = 9, face = "bold", family = FONT))
save_pdf(p_supp_res, "figure_supp_additional_results", W_DOUBLE, 170)
save_tiff(p_supp_res, "figure_supp_additional_results", W_DOUBLE, 170, dpi = 600)

# ===================== Final structural checks =====================
stopifnot(nrow(fp_data) == 12)
stopifnot(nrow(ff_data) == 6)
stopifnot(nrow(love_dt) == 13)
stopifnot(nrow(med_fu) == 12)
stopifnot(nrow(trace_dt) == N_ITER * M_IMP * length(mice_imp_vars))
