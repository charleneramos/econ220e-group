###############################################################################
## Monte Carlo Simulations for Roth & Sant'Anna (2023)
## "When Is Parallel Trends Sensitive to Functional Form?"
##
## Authors: Roberto Diaz Hernandez, Riana Ramonjamanana, Charlene Ramos
## Date:    March 2026
##
###############################################################################

# ---- 0. Setup ----------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
library(xtable)

set.seed(20260318)

B     <- 2000    # number of Monte Carlo replications
n     <- 1000    # observations per group-period cell
alpha <- 0.05    # significance level

# ---- Shared ggplot theme for all figures -------------------------------------
theme_paper <- theme_minimal(base_size = 11) +
  theme(
    legend.position   = "bottom",
    plot.title         = element_text(size = 11, face = "bold", hjust = 0.5),
    plot.subtitle      = element_text(size = 9, hjust = 0.5, color = "grey40"),
    panel.grid.minor   = element_blank(),
    axis.title         = element_text(size = 10),
    legend.text        = element_text(size = 9),
    legend.title       = element_blank(),
    strip.text         = element_text(size = 10, face = "bold")
  )

# Color scheme
col_ctrl <- "#D55E00"  
col_trt  <- "#0072B2"   
col_imp  <- "#E69F00"   

# ---- 1. Helper: draw from a 50/50 lognormal mixture -------------------------
rmix_lognorm <- function(n, mu1, sigma1, mu2, sigma2, theta = 0.5) {
  u <- rbinom(n, 1, theta)
  u * rlnorm(n, mu1, sigma1) + (1 - u) * rlnorm(n, mu2, sigma2)
}

# ---- 2. DGP 1 — Favorable Case -----------
run_dgp1 <- function(n) {
  Y_00 <- rmix_lognorm(n, mu1 = 2, sigma1 = 1, mu2 = 3, sigma2 = 1)
  Y_01 <- rmix_lognorm(n, mu1 = 3, sigma1 = 1, mu2 = 3, sigma2 = 1)
  Y_10 <- rmix_lognorm(n, mu1 = 2, sigma1 = 1, mu2 = 4, sigma2 = 1)
  Y_11 <- rmix_lognorm(n, mu1 = 3, sigma1 = 1, mu2 = 4, sigma2 = 1)

  did_lev <- (mean(Y_11) - mean(Y_10)) - (mean(Y_01) - mean(Y_00))
  se_lev  <- sqrt(var(Y_11)/n + var(Y_10)/n + var(Y_01)/n + var(Y_00)/n)

  did_log <- (mean(log(Y_11)) - mean(log(Y_10))) -
             (mean(log(Y_01)) - mean(log(Y_00)))
  se_log  <- sqrt(var(log(Y_11))/n + var(log(Y_10))/n +
                  var(log(Y_01))/n + var(log(Y_00))/n)

  c(did_lev = did_lev, se_lev = se_lev,
    did_log = did_log, se_log = se_log)
}

# ---- 3. DGP 2 — Unfavorable Case -------------------------
run_dgp2 <- function(n) {
  Y_00 <- rlnorm(n, meanlog = 2.5, sdlog = 1.0)
  Y_01 <- rlnorm(n, meanlog = 3.0, sdlog = 1.0)
  Y_10 <- rlnorm(n, meanlog = 3.0, sdlog = 0.5)
  Y_11 <- rlnorm(n, meanlog = 3.5, sdlog = 0.5)

  did_lev <- (mean(Y_11) - mean(Y_10)) - (mean(Y_01) - mean(Y_00))
  se_lev  <- sqrt(var(Y_11)/n + var(Y_10)/n + var(Y_01)/n + var(Y_00)/n)

  did_log <- (mean(log(Y_11)) - mean(log(Y_10))) -
             (mean(log(Y_01)) - mean(log(Y_00)))
  se_log  <- sqrt(var(log(Y_11))/n + var(log(Y_10))/n +
                  var(log(Y_01))/n + var(log(Y_00))/n)

  c(did_lev = did_lev, se_lev = se_lev,
    did_log = did_log, se_log = se_log)
}

# ---- 4. Run Monte Carlo Simulations -----------------------------------------
cat("Running Monte Carlo simulations (B =", B, ", n =", n, ")...\n")

mc1 <- t(replicate(B, run_dgp1(n)))
mc2 <- t(replicate(B, run_dgp2(n)))

cat("Done.\n\n")

# ---- 5. Summarize Results ----------------------------------------------------
summarize_mc <- function(mc, true_lev = 0, true_log = 0) {
  data.frame(
    Transformation = c("Levels", "Logs"),
    True_DiD       = c(true_lev, true_log),
    Mean_DiD       = c(mean(mc[, "did_lev"]), mean(mc[, "did_log"])),
    Bias           = c(mean(mc[, "did_lev"]) - true_lev,
                       mean(mc[, "did_log"]) - true_log),
    SD             = c(sd(mc[, "did_lev"]), sd(mc[, "did_log"])),
    Mean_SE        = c(mean(mc[, "se_lev"]), mean(mc[, "se_log"])),
    Reject_Rate    = c(
      mean(abs(mc[, "did_lev"] / mc[, "se_lev"]) > qnorm(1 - alpha/2)),
      mean(abs(mc[, "did_log"] / mc[, "se_log"]) > qnorm(1 - alpha/2))
    )
  )
}

tab1 <- summarize_mc(mc1, true_lev = 0, true_log = 0)
tab2 <- summarize_mc(mc2, true_lev = 0, true_log = 0)

cat("=== DGP 1: Favorable Case (Mixture Satisfied) ===\n")
print(tab1, digits = 4, row.names = FALSE)
cat("\n=== DGP 2: Unfavorable Case (Mixture Violated) ===\n")
print(tab2, digits = 4, row.names = FALSE)

# ---- 6. Export LaTeX Tables --------------------------------------------------
make_latex_table <- function(tab, caption, label) {
  tab_print <- tab
  tab_print$True_DiD    <- sprintf("%.3f", tab$True_DiD)
  tab_print$Mean_DiD    <- sprintf("%.3f", tab$Mean_DiD)
  tab_print$Bias        <- sprintf("%.3f", tab$Bias)
  tab_print$SD          <- sprintf("%.3f", tab$SD)
  tab_print$Mean_SE     <- sprintf("%.3f", tab$Mean_SE)
  tab_print$Reject_Rate <- sprintf("%.3f", tab$Reject_Rate)
  names(tab_print) <- c("$g(\\cdot)$", "True DiD", "Mean $\\hat{\\tau}$",
                         "Bias", "SD", "Mean SE", "Rej.\\ Rate")
  xt <- xtable(tab_print, align = "llcccccc")
  print(xt, file = file.path("output/tables", paste0(label, ".tex")),
        sanitize.text.function = identity,
        include.rownames = FALSE,
        booktabs = TRUE,
        floating = FALSE)
}

make_latex_table(tab1,
  caption = "Monte Carlo results for DGP~1 (favorable case). $B = 2{,}000$ replications, $n = 1{,}000$ per cell. True ATT $= 0$ under all transformations.",
  label = "tab_dgp1")

make_latex_table(tab2,
  caption = "Monte Carlo results for DGP~2 (unfavorable case). $B = 2{,}000$ replications, $n = 1{,}000$ per cell. True ATT $= 0$; PT holds in logs but not levels.",
  label = "tab_dgp2")

cat("Tables saved to tables/\n")

###############################################################################
# ---- 7. Figures                                                           ----
###############################################################################

## Mixture lognormal density helper
dmix_lognorm <- function(x, mu1, sigma1, mu2, sigma2, theta = 0.5) {
  theta * dlnorm(x, mu1, sigma1) + (1 - theta) * dlnorm(x, mu2, sigma2)
}

x_grid <- seq(0.01, 150, length.out = 2000)

# =============================================================================
# FIGURE 1 — Density plots 
# =============================================================================

dens1 <- data.frame(
  y = rep(x_grid, 4),
  density = c(
    dmix_lognorm(x_grid, 2, 1, 3, 1),
    dmix_lognorm(x_grid, 3, 1, 3, 1),
    dmix_lognorm(x_grid, 2, 1, 4, 1),
    dmix_lognorm(x_grid, 3, 1, 4, 1)
  ),
  Group  = rep(c("Control", "Control", "Treated", "Treated"), each = length(x_grid)),
  Period = rep(c("Pre", "Post", "Pre", "Post"), each = length(x_grid))
)
dens1$DGP <- "DGP 1"

dens2 <- data.frame(
  y = rep(x_grid, 4),
  density = c(
    dlnorm(x_grid, 2.5, 1.0),
    dlnorm(x_grid, 3.0, 1.0),
    dlnorm(x_grid, 3.0, 0.5),
    dlnorm(x_grid, 3.5, 0.5)
  ),
  Group  = rep(c("Control", "Control", "Treated", "Treated"), each = length(x_grid)),
  Period = rep(c("Pre", "Post", "Pre", "Post"), each = length(x_grid))
)
dens2$DGP <- "DGP 2"

dens_all <- rbind(dens1, dens2)

p_dens <- ggplot(dens_all, aes(x = y, y = density, color = Group, linetype = Period)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ DGP, ncol = 2) +
  scale_color_manual(values = c("Control" = col_ctrl, "Treated" = col_trt)) +
  scale_linetype_manual(values = c("Pre" = "dashed", "Post" = "solid")) +
  labs(x = expression(italic(y)), y = "Density") +
  coord_cartesian(xlim = c(0, 120)) +
  theme_paper +
  guides(color = guide_legend(order = 1), linetype = guide_legend(order = 2))

ggsave("output/figures/fig_densities.png", p_dens, width = 7, height = 3.5)
cat("Saved: figures/fig_densities.png\n")

# =============================================================================
# FIGURE 2 — Change-in-density diagnostic 
# =============================================================================

delta_ctrl_1 <- function(y) dmix_lognorm(y, 3, 1, 3, 1) - dmix_lognorm(y, 2, 1, 3, 1)
delta_trt_1  <- function(y) dmix_lognorm(y, 3, 1, 4, 1) - dmix_lognorm(y, 2, 1, 4, 1)

delta_ctrl_2 <- function(y) dlnorm(y, 3.0, 1.0) - dlnorm(y, 2.5, 1.0)
delta_trt_2  <- function(y) dlnorm(y, 3.5, 0.5) - dlnorm(y, 3.0, 0.5)

delta_df <- rbind(
  data.frame(y = rep(x_grid, 2),
             delta = c(delta_ctrl_1(x_grid), delta_trt_1(x_grid)),
             Group = rep(c("Control", "Treated"), each = length(x_grid)),
             DGP = "DGP 1"),
  data.frame(y = rep(x_grid, 2),
             delta = c(delta_ctrl_2(x_grid), delta_trt_2(x_grid)),
             Group = rep(c("Control", "Treated"), each = length(x_grid)),
             DGP = "DGP 2")
)

p_delta <- ggplot(delta_df, aes(x = y, y = delta, color = Group, linetype = Group)) +
  geom_line(linewidth = 0.7) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50", linetype = "dotted") +
  facet_wrap(~ DGP, ncol = 2) +
  scale_color_manual(values = c("Control" = col_ctrl, "Treated" = col_trt)) +
  labs(x = expression(italic(y)), y = expression(Delta ~ "density")) +
  coord_cartesian(xlim = c(0, 120)) +
  theme_paper

ggsave("output/figures/fig_delta_density.png", p_delta, width = 7, height = 3.5)
cat("Saved: output/figures/fig_delta_density.png\n")

# =============================================================================
# FIGURE 3 — Implied vs. True Counterfactual CDF 
# =============================================================================

y_cdf <- seq(0.5, 150, length.out = 2000)

# --- DGP 1 CDF components (mixture) ---
pmix <- function(y, mu1, sigma1, mu2, sigma2, theta = 0.5) {
  theta * plnorm(y, mu1, sigma1) + (1 - theta) * plnorm(y, mu2, sigma2)
}

F1_ctrl_pre  <- function(y) pmix(y, 2, 1, 3, 1)
F1_ctrl_post <- function(y) pmix(y, 3, 1, 3, 1)
F1_trt_pre   <- function(y) pmix(y, 2, 1, 4, 1)
F1_trt_post_true <- function(y) pmix(y, 3, 1, 4, 1)  # true counterfactual
F1_implied   <- function(y) F1_trt_pre(y) + F1_ctrl_post(y) - F1_ctrl_pre(y)

# --- DGP 2 CDF components ---
F2_ctrl_pre  <- function(y) plnorm(y, 2.5, 1.0)
F2_ctrl_post <- function(y) plnorm(y, 3.0, 1.0)
F2_trt_pre   <- function(y) plnorm(y, 3.0, 0.5)
F2_trt_post_true <- function(y) plnorm(y, 3.5, 0.5)
F2_implied   <- function(y) F2_trt_pre(y) + F2_ctrl_post(y) - F2_ctrl_pre(y)

cdf_df <- rbind(
  data.frame(y = rep(y_cdf, 2),
             CDF = c(F1_implied(y_cdf), F1_trt_post_true(y_cdf)),
             Type = rep(c("Implied Counterfactual", "True Counterfactual"), each = length(y_cdf)),
             DGP = "DGP 1"),
  data.frame(y = rep(y_cdf, 2),
             CDF = c(F2_implied(y_cdf), F2_trt_post_true(y_cdf)),
             Type = rep(c("Implied Counterfactual", "True Counterfactual"), each = length(y_cdf)),
             DGP = "DGP 2")
)

p_cdf <- ggplot(cdf_df, aes(x = y, y = CDF, color = Type, linetype = Type)) +
  geom_line(linewidth = 0.7) +
  facet_wrap(~ DGP, ncol = 2) +
  scale_color_manual(values = c("Implied Counterfactual" = col_imp,
                                "True Counterfactual" = col_trt)) +
  scale_linetype_manual(values = c("Implied Counterfactual" = "dashed",
                                   "True Counterfactual" = "solid")) +
  geom_hline(yintercept = c(0, 1), linewidth = 0.3, color = "grey50", linetype = "dotted") +
  labs(x = expression(italic(y)), y = "CDF") +
  coord_cartesian(xlim = c(0, 120), ylim = c(-0.05, 1.05)) +
  theme_paper

ggsave("output/figures/fig_counterfactual_cdf.png", p_cdf, width = 7, height = 3.5)
cat("Saved: output/figures/fig_counterfactual_cdf.png\n")

# Check DGP 2 non-monotonicity
implied2_vals <- F2_implied(y_cdf)
non_mono <- any(diff(implied2_vals) < -1e-10)
cat("\nDGP 2 implied CDF is",
    ifelse(non_mono, "NON-MONOTONE (violation detected)\n",
                     "monotone (no violation detected)\n"))
cat("  Min implied CDF:", min(implied2_vals), "\n")
cat("  Max implied CDF:", max(implied2_vals), "\n")

# Also check DGP 1
implied1_vals <- F1_implied(y_cdf)
non_mono1 <- any(diff(implied1_vals) < -1e-10)
cat("DGP 1 implied CDF is",
    ifelse(non_mono1, "NON-MONOTONE (violation detected)\n",
                      "monotone (no violation)\n"))

# =============================================================================
# FIGURE 4 — Sampling distribution of DiD estimator (2x2 faceted)
# =============================================================================

hist_df <- data.frame(
  DiD = c(mc1[, "did_lev"], mc2[, "did_lev"]),
  DGP = rep(c("DGP 1 (Favorable)", "DGP 2 (Unfavorable)"), each = B)
)

# Add vertical lines for the true DiD values
vline_df <- data.frame(
  DGP = c("DGP 1 (Favorable)", "DGP 2 (Unfavorable)"),
  true_val = c(0, 0)  # H0 value under test
)

p_hist <- ggplot(hist_df, aes(x = DiD)) +
  geom_histogram(aes(y = after_stat(density)), bins = 50,
                 fill = col_trt, alpha = 0.55, color = "white", linewidth = 0.2) +
  geom_vline(data = vline_df, aes(xintercept = true_val),
             linewidth = 0.6, linetype = "dashed", color = col_ctrl) +
  facet_wrap(~ DGP, ncol = 2, scales = "free") +
  labs(x = expression(hat(tau)[DiD]), y = "Density") +
  theme_paper +
  theme(strip.text = element_text(size = 9))

ggsave("output/figures/fig_sampling_dist.png", p_hist, width = 6.5, height = 3.5)
cat("Saved: output/figures/fig_sampling_dist.png\n")

# =============================================================================
# FIGURE 5 — Sample-size sensitivity: Bias and Rejection Rate vs. n
# =============================================================================

cat("\nRunning sample-size sensitivity analysis...\n")

n_vec   <- c(100, 250, 500, 1000, 2500, 5000)
B_sens  <- 1000  # fewer reps for speed

sens_results <- list()

for (nn in n_vec) {
  cat("  n =", nn, "...")

  mc1_s <- t(replicate(B_sens, run_dgp1(nn)))
  mc2_s <- t(replicate(B_sens, run_dgp2(nn)))

  sens_results[[as.character(nn)]] <- data.frame(
    n = nn,
    DGP = rep(c("DGP 1 (Favorable)", "DGP 2 (Unfavorable)"), each = 2),
    Transformation = rep(c("Levels", "Logs"), 2),
    Bias = c(
      mean(mc1_s[, "did_lev"]),
      mean(mc1_s[, "did_log"]),
      mean(mc2_s[, "did_lev"]),
      mean(mc2_s[, "did_log"])
    ),
    Reject_Rate = c(
      mean(abs(mc1_s[, "did_lev"] / mc1_s[, "se_lev"]) > qnorm(1 - alpha/2)),
      mean(abs(mc1_s[, "did_log"] / mc1_s[, "se_log"]) > qnorm(1 - alpha/2)),
      mean(abs(mc2_s[, "did_lev"] / mc2_s[, "se_lev"]) > qnorm(1 - alpha/2)),
      mean(abs(mc2_s[, "did_log"] / mc2_s[, "se_log"]) > qnorm(1 - alpha/2))
    ),
    RMSE = c(
      sqrt(mean(mc1_s[, "did_lev"]^2)),
      sqrt(mean(mc1_s[, "did_log"]^2)),
      sqrt(mean(mc2_s[, "did_lev"]^2)),
      sqrt(mean(mc2_s[, "did_log"]^2))
    )
  )
  cat(" done\n")
}

sens_df <- do.call(rbind, sens_results)

# Panel A: Rejection rate vs n
p_sens_rej <- ggplot(sens_df, aes(x = n, y = Reject_Rate,
                                   color = Transformation, shape = Transformation)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 0.05, linewidth = 0.4, linetype = "dashed", color = "grey40") +
  facet_wrap(~ DGP, ncol = 2) +
  scale_color_manual(values = c("Levels" = col_ctrl, "Logs" = col_trt)) +
  scale_x_log10(breaks = n_vec, labels = scales::comma) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = expression(italic(n) ~ "(per cell)"), y = "Rejection Rate (5% level)") +
  theme_paper 

ggsave("output/figures/fig_sensitivity_rejection.png", p_sens_rej, width = 7, height = 3.5)
cat("Saved: output/figures/fig_sensitivity_rejection.png\n")

# Panel B: RMSE vs n
p_sens_rmse <- ggplot(sens_df, aes(x = n, y = RMSE,
                                    color = Transformation, shape = Transformation)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 2.5) +
  facet_wrap(~ DGP, ncol = 2, scales = "free_y") +
  scale_color_manual(values = c("Levels" = col_ctrl, "Logs" = col_trt)) +
  scale_x_log10(breaks = n_vec, labels = scales::comma) +
  labs(x = expression(italic(n) ~ "(per cell)"), y = "RMSE") +
  theme_paper

ggsave("output/figures/fig_sensitivity_rmse.png", p_sens_rmse, width = 7, height = 3.5)
cat("Saved: output/figures/fig_sensitivity_rmse.png\n")



