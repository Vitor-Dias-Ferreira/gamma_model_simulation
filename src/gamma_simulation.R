# ==========================================
# Monte Carlo study: Gamma model, MLE vs. method of moments
# Run from the repository root: Rscript src/gamma_simulation.R
# ==========================================

# 0. Setup
set.seed(42) # Fixed seed for reproducible results
alpha_0 <- 2 # True shape
beta_0 <- 3  # True scale
sample_sizes <- c(20, 50, 100, 200, 500)
R <- 10000   # Replications per sample size

# 1. Method of moments (MM): closed form
estimate_moments <- function(x) {
  n <- length(x)
  x_bar <- mean(x)

  # var() uses the n - 1 divisor; the assignment defines S^2 with divisor n
  s2 <- var(x) * ((n - 1) / n)

  alpha_mm <- (x_bar^2) / s2
  beta_mm <- s2 / x_bar

  return(c(alpha = alpha_mm, beta = beta_mm))
}

# Quick check on a single sample
test_data <- rgamma(n = 100, shape = alpha_0, scale = beta_0)
print("Method of moments estimates:")
print(estimate_moments(test_data))


# ==========================================
# 2. Maximum likelihood (Newton-Raphson)
# ==========================================

estimate_mle <- function(x, tol = 1e-8, max_iter = 100) {
  n <- length(x)
  x_bar <- mean(x)

  # Constant term of the profile score equation
  constant <- log(x_bar) - mean(log(x))

  # Method of moments as the starting value
  s2 <- var(x) * ((n - 1) / n)
  alpha_init <- (x_bar^2) / s2

  # Work on the log scale to keep alpha positive (eta = ln(alpha))
  eta <- log(alpha_init)

  for (i in 1:max_iter) {
    alpha_current <- exp(eta)

    # Equation to solve
    g_val <- eta - digamma(alpha_current) - constant

    # Derivative with respect to eta
    g_deriv <- 1 - trigamma(alpha_current) * alpha_current

    eta_new <- eta - (g_val / g_deriv)

    if (abs(eta_new - eta) < tol) {
      eta <- eta_new
      break
    }

    eta <- eta_new
  }

  # Back to the original scale
  alpha_mle <- exp(eta)
  beta_mle <- x_bar / alpha_mle

  return(c(alpha = alpha_mle, beta = beta_mle))
}

print("Maximum likelihood estimates:")
print(estimate_mle(test_data))


# ==========================================
# 3. Standard errors (delta method, expected and observed information)
# ==========================================

estimate_standard_errors <- function(x, est_mm, est_mle) {
  n <- length(x)

  # --- 1. Method of moments (MM) ---
  a_mm <- est_mm["alpha"]
  b_mm <- est_mm["beta"]

  # Delta-method asymptotic covariance for the shape-scale parameterization.
  # This replaces the matrix printed in the assignment; see README.md.
  Sigma_MM <- matrix(c(
    2 * a_mm * (a_mm + 1),     -2 * b_mm * (a_mm + 1),
    -2 * b_mm * (a_mm + 1),    b_mm^2 * (2 * a_mm + 3) / a_mm
  ), nrow = 2, byrow = TRUE)

  se_mm <- sqrt(diag(Sigma_MM) / n)

  # --- 2. MLE, expected Fisher information (MLE-IF) ---
  a_mle <- est_mle["alpha"]
  b_mle <- est_mle["beta"]

  I_F <- matrix(c(
    trigamma(a_mle),   1 / b_mle,
    1 / b_mle,         a_mle / (b_mle^2)
  ), nrow = 2, byrow = TRUE)

  se_mle_if <- sqrt(diag(solve(n * I_F)))

  # --- 3. MLE, observed information (MLE-IO) ---
  H11 <- -n * trigamma(a_mle)
  H12 <- -n / b_mle
  H22 <- (n * a_mle / (b_mle^2)) - (2 * sum(x) / (b_mle^3))
  Hessian <- matrix(c(H11, H12, H12, H22), nrow = 2, byrow = TRUE)

  se_mle_io <- sqrt(diag(solve(-Hessian)))

  return(list(
    SE_MM = setNames(se_mm, c("alpha", "beta")),
    SE_MLE_IF = setNames(se_mle_if, c("alpha", "beta")),
    SE_MLE_IO = setNames(se_mle_io, c("alpha", "beta"))
  ))
}

print("Standard errors:")
est_mm_test <- estimate_moments(test_data)
est_mle_test <- estimate_mle(test_data)
print(estimate_standard_errors(test_data, est_mm_test, est_mle_test))


# ==========================================
# 4. Monte Carlo loop
# ==========================================

z_val <- qnorm(0.975) # 95% confidence
result_list <- list()

compute_metrics <- function(n_current, est, se, theta_0, method_name, param) {
  lower <- est - z_val * se
  upper <- est + z_val * se

  data.frame(
    n = n_current,
    Parameter = param,
    Method = method_name,
    Mean = mean(est),
    Mean_SE = mean(se),
    SD = sd(est),
    Bias = mean(est - theta_0),
    MSE = mean((est - theta_0)^2),
    Width = mean(upper - lower),
    Coverage = mean(theta_0 >= lower & theta_0 <= upper),
    Alpha1 = mean(theta_0 < lower),  # Interval entirely above the true value
    Alpha2 = mean(theta_0 > upper)   # Interval entirely below the true value
  )
}

for (n in sample_sizes) {
  cat(sprintf("Running %d replications for n = %d...\n", R, n))

  # R x 2 matrices for this sample size
  est_mm_mat <- matrix(NA, nrow = R, ncol = 2)
  est_mle_mat <- matrix(NA, nrow = R, ncol = 2)
  se_mm_mat <- matrix(NA, nrow = R, ncol = 2)
  se_mle_if_mat <- matrix(NA, nrow = R, ncol = 2)
  se_mle_io_mat <- matrix(NA, nrow = R, ncol = 2)

  for (j in 1:R) {
    x <- rgamma(n, shape = alpha_0, scale = beta_0)

    est_mm <- estimate_moments(x)
    est_mle <- estimate_mle(x)

    errors <- estimate_standard_errors(x, est_mm, est_mle)

    est_mm_mat[j, ] <- est_mm
    est_mle_mat[j, ] <- est_mle
    se_mm_mat[j, ] <- errors$SE_MM
    se_mle_if_mat[j, ] <- errors$SE_MLE_IF
    se_mle_io_mat[j, ] <- errors$SE_MLE_IO
  }

  table_n <- rbind(
    # --- Alpha ---
    compute_metrics(n, est_mm_mat[,1], se_mm_mat[,1], alpha_0, "MM", "Alpha"),
    compute_metrics(n, est_mle_mat[,1], se_mle_if_mat[,1], alpha_0, "MLE-IF", "Alpha"),
    compute_metrics(n, est_mle_mat[,1], se_mle_io_mat[,1], alpha_0, "MLE-IO", "Alpha"),

    # --- Beta ---
    compute_metrics(n, est_mm_mat[,2], se_mm_mat[,2], beta_0, "MM", "Beta"),
    compute_metrics(n, est_mle_mat[,2], se_mle_if_mat[,2], beta_0, "MLE-IF", "Beta"),
    compute_metrics(n, est_mle_mat[,2], se_mle_io_mat[,2], beta_0, "MLE-IO", "Beta")
  )

  result_list[[as.character(n)]] <- table_n
}

final_table <- do.call(rbind, result_list)
rownames(final_table) <- NULL

print("Simulation complete. First rows of the consolidated table:")
print(head(final_table))


# ==========================================
# 5. Export table and figures
# ==========================================
library(ggplot2)

results_dir <- "results"
if(!dir.exists(results_dir)) dir.create(results_dir)

table_path <- file.path(results_dir, "full_results.csv")
write.csv(final_table, table_path, row.names = FALSE)
cat("Table saved to:", table_path, "\n")

# MLE-IF and MLE-IO coincide, so the figures show MM and MLE-IF
plot_data <- subset(final_table, Method %in% c("MM", "MLE-IF"))

g_mse <- ggplot(plot_data, aes(x = n, y = MSE, color = Method, group = Method)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  facet_wrap(~Parameter, scales = "free_y") +
  theme_minimal() +
  labs(title = "Mean squared error vs. sample size",
       y = "MSE", x = "Sample size (n)") +
  scale_color_brewer(palette = "Set1")

g_coverage <- ggplot(plot_data, aes(x = n, y = Coverage, color = Method, group = Method)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") +
  facet_wrap(~Parameter) +
  theme_minimal() +
  labs(title = "Coverage probability vs. sample size",
       subtitle = "Dashed line: nominal 95%",
       y = "Coverage", x = "Sample size (n)") +
  scale_color_brewer(palette = "Set1")

g_width <- ggplot(plot_data, aes(x = n, y = Width, color = Method, group = Method)) +
  geom_line(linewidth = 1) + geom_point(size = 3) +
  facet_wrap(~Parameter, scales = "free_y") +
  theme_minimal() +
  labs(title = "Mean interval width vs. sample size",
       y = "Interval width", x = "Sample size (n)") +
  scale_color_brewer(palette = "Set1")

ggsave(file.path(results_dir, "mse.png"), plot = g_mse, width = 8, height = 5, bg = "white")
ggsave(file.path(results_dir, "coverage.png"), plot = g_coverage, width = 8, height = 5, bg = "white")
ggsave(file.path(results_dir, "width.png"), plot = g_width, width = 8, height = 5, bg = "white")

cat("Figures saved to:", results_dir, "\n")
