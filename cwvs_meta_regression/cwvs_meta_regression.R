#!/usr/bin/env Rscript
# -*- coding:utf-8 -*-

##############################################################
# Critical Window Variable Selection (CWVS) Meta-Regression
# MCMC analysis pipeline
##############################################################

# 1. Utility functions -----------------------------------------------------------

inv_logit <- plogis   # built‑in logistic inverse
logit     <- qlogis   # built‑in logit

ci95 <- function(x) quantile(x, c(.025, .5, .975), names = FALSE)

drop_burn <- function(x, burn) {
  n <- if (is.matrix(x)) nrow(x) else length(x)
  k <- if (burn < 1) ceiling(n * burn) else min(burn, n - 1)
  if (is.matrix(x)) x[-seq_len(k), , drop = FALSE] else x[-seq_len(k)]
}

trace_hist <- function(draws, name, breaks = "FD") {
  op <- par(mfrow = c(1, 2)); on.exit(par(op), add = TRUE)
  plot(draws, type = "l", xlab = "iter", ylab = name,
       main = paste("Trace:", name))
  hist(draws, breaks = breaks, freq = FALSE,
       main = paste("Histogram:", name),
       xlab = name, ylab = "density", border = "grey")
}

rw_update <- function(current, sd, logpost) {
  prop <- current + rnorm(1, 0, sd)
  if (log(runif(1)) < logpost(prop) - logpost(current)) prop else current
}

log_ar1_prior <- function(delta, rho) {
  n <- length(delta)
  lp <- dnorm(delta[1], 0, sqrt(1 - rho^2), log = TRUE)
  if (n > 1)
    lp <- lp + sum(dnorm(delta[-1], rho * delta[-n], log = TRUE))
  lp
}

neigh_lp <- function(k, val, vec, rho) {
  n <- length(vec)
  if (k == 1) {
    lp <- dnorm(val, 0, sqrt(1 - rho^2), log = TRUE)
    if (n > 1) lp <- lp + dnorm(vec[2], rho * val, log = TRUE)
  } else if (k == n) {
    lp <- dnorm(val, rho * vec[n - 1], log = TRUE)
  } else {
    lp <- dnorm(val,           rho * vec[k - 1], log = TRUE) +
          dnorm(vec[k + 1], rho * val, log = TRUE)
  }
  lp
}

logbern <- function(g, d1, d2, A21, A22) {
  eta <- A21 * d1 + A22 * d2
  dbinom(g, 1, pnorm(eta), log = TRUE)
}

# 2. Sampler (Gibbs in Metropolis-Hastings) -----------------------------------------------------------

cwvs_mcmc <- function(alpha_hat, sigma_hat,
                      n_iter = 2e4, thin = 20,
                      prop = list(),
                      init = NULL,
                      seed = NULL,
                      verbose = TRUE,
                      save_all = TRUE) {
  if (!is.null(seed)) set.seed(seed)

  # ---- proposal SD ----
  prop <- modifyList(list(sd_delta = 0.05,
                          sd_A     = 0.02,
                          sd_A21   = 0.05,
                          sd_rho   = 0.10), prop)

  # ---- data ----
  alpha_hat <- as.numeric(alpha_hat)
  sigma_hat <- pmax(as.numeric(sigma_hat), 1e-6)  # avoid 0
  N <- length(alpha_hat)

  # ---- initial values ----
  if (is.null(init)) {
    delta1 <- delta2 <- gamma <- numeric(N)
    A11 <- 0.1;  A22 <- 0.1; A21 <- 0
    rho1 <- 0.2; rho2 <- 0.2
  } else list2env(init, environment())

  # ---- storage ----
  n_save <- floor(n_iter / thin)
  gamma_save <- theta_save <- matrix(NA_real_, n_save, N)
  if (save_all) {
    A11_save <- A22_save <- A21_save <- rho1_save <- rho2_save <- numeric(n_save)
  }
  save_idx <- 1L

  # ---- MCMC ----
  for (iter in seq_len(n_iter)) {
    ## 1.  Gibbs for gamma_t ------------------------------
    theta_vec <- A11 * delta1
    eta_vec   <- A21 * delta1 + A22 * delta2
    pi_vec    <- pnorm(eta_vec)

    loglik1 <- dnorm(alpha_hat, theta_vec, sigma_hat, log = TRUE)
    loglik0 <- dnorm(alpha_hat, 0,          sigma_hat, log = TRUE)
    logp1   <- log(pi_vec)    + loglik1
    logp0   <- log1p(-pi_vec) + loglik0
    p1      <- plogis(logp1 - logp0)  # numerically stable
    gamma   <- rbinom(N, 1, p1)

    ## 2.  Random‑walk MH for delta1_t, delta2_t ----------
    for (t in seq_len(N)) {
      # delta1_t
      prop_d1 <- delta1[t] + rnorm(1, 0, prop$sd_delta)
      logratio <- dnorm(alpha_hat[t], A11 * prop_d1 * gamma[t], sigma_hat[t], log = TRUE) -
                  dnorm(alpha_hat[t], A11 * delta1[t] * gamma[t], sigma_hat[t], log = TRUE) +
                  logbern(gamma[t], prop_d1, delta2[t], A21, A22) -
                  logbern(gamma[t], delta1[t], delta2[t], A21, A22)
      lp_new <- neigh_lp(t, prop_d1, delta1, rho1)
      lp_old <- neigh_lp(t, delta1[t], delta1, rho1)
      if (log(runif(1)) < logratio + lp_new - lp_old)
        delta1[t] <- prop_d1

      # delta2_t
      prop_d2 <- delta2[t] + rnorm(1, 0, prop$sd_delta)
      logratio <- logbern(gamma[t], delta1[t], prop_d2, A21, A22) -
                  logbern(gamma[t], delta1[t], delta2[t], A21, A22)
      lp_new <- neigh_lp(t, prop_d2, delta2, rho2)
      lp_old <- neigh_lp(t, delta2[t], delta2, rho2)
      if (log(runif(1)) < logratio + lp_new - lp_old)
        delta2[t] <- prop_d2
    }

    ## 3.  Update A11, A22, A21 ---------------------------
    A11 <- rw_update(A11, prop$sd_A, function(a) {
      ifelse(a <= 0 || a >= 100, -Inf,
             sum(dnorm(alpha_hat, a * delta1 * gamma, sigma_hat, log = TRUE)))
    })

    A22 <- rw_update(A22, prop$sd_A, function(a) {
      ifelse(a <= 0 || a >= 100, -Inf,
             sum(dbinom(gamma, 1, pnorm(A21 * delta1 + a * delta2), log = TRUE)))
    })

    A21 <- rw_update(A21, prop$sd_A21, function(a) {
      sum(dbinom(gamma, 1, pnorm(a * delta1 + A22 * delta2), log = TRUE)) +
        dnorm(a, 0, 100, log = TRUE)
    })

    ## 4.  Update rho1, rho2 (logit space) ----------------
    fisher_rw <- function(rho, sd, delta) {
      z_curr <- logit(rho)
      z_new  <- rw_update(z_curr, sd, function(z) {
        r <- inv_logit(z)
        if (r <= 0 || r >= 1) return(-Inf)
        log_ar1_prior(delta, r) + log(r) + log(1 - r)  # Jacobian
      })
      inv_logit(z_new)
    }
    rho1 <- fisher_rw(rho1, prop$sd_rho, delta1)
    rho2 <- fisher_rw(rho2, prop$sd_rho, delta2)

    ## 5.  Save draws ------------------------------------
    if (iter %% thin == 0L) {
      gamma_save[save_idx, ] <- gamma
      theta_save[save_idx, ] <- A11 * delta1
      if (save_all) {
        A11_save[save_idx] <- A11; A22_save[save_idx] <- A22; A21_save[save_idx] <- A21
        rho1_save[save_idx] <- rho1; rho2_save[save_idx] <- rho2
      }
      if (verbose && (save_idx %% 50L == 0L))
        cat(sprintf("saved %d / %d\n", save_idx, n_save))
      save_idx <- save_idx + 1L
    }
  }

  out <- list(gamma = gamma_save, theta = theta_save)
  if (save_all)
    out <- c(out, list(A11 = A11_save, A22 = A22_save, A21 = A21_save,
                       rho1 = rho1_save, rho2 = rho2_save))
  class(out) <- c("cwvs_mcmc", class(out))
  out
}

# 3. Perform analysis --------------------------------------------------------------

dat <- read.csv("synthetic_cwvs_input.csv")
alpha_hat <- log(dat$OR)
sigma_hat <- (log(dat$UCI) - log(dat$LCI)) / (2 * 1.96)

fit <- cwvs_mcmc(alpha_hat, sigma_hat,
                 n_iter = 1e5, thin = 20,
                 prop = list(sd_delta = 0.05,
                             sd_A     = 0.02,
                             sd_A21   = 0.05,
                             sd_rho   = 0.10),
                 verbose = TRUE)

# 4. Posterior summaries & plots ------------------------------------------

burn <- 200
week_start <- -15
weeks <- ncol(fit$theta)
week_lab <- seq(week_start, by = 1, length.out = weeks)

## ---- scalar parameters ----

scalar_params <- c("A11", "A22", "A21", "rho1", "rho2")
scalar_summ <- sapply(scalar_params, function(p) {
  vec <- drop_burn(fit[[p]], burn)
  c(mean(vec), sd(vec), ci95(vec))
})

cat("\n===== Scalar parameter summary =====\n")
print(round(t(scalar_summ), 4), row.names = FALSE)

## ---- theta_t ----

theta_mat <- drop_burn(as.matrix(fit$theta), burn)
theta_q <- t(apply(theta_mat, 2, ci95))
theta_tab <- data.frame(
  week   = week_lab,
  mean   = colMeans(theta_mat),
  sd     = apply(theta_mat, 2, sd),
  `2.5%` = theta_q[, 1],
  median = theta_q[, 2],
  `97.5%`= theta_q[, 3]
)
cat("\n===== theta_t summary =====\n")
print(round(theta_tab, 4), row.names = FALSE)

## ---- gamma_t ----

gamma_mat <- drop_burn(as.matrix(fit$gamma), burn)
gamma_bar <- colMeans(gamma_mat)
gamma_tab <- data.frame(week = week_lab, P_incl = round(gamma_bar, 4))
cat("\n===== gamma_t inclusion probability =====\n")
print(gamma_tab, row.names = FALSE)

## ---- diagnostics ----

for (w in seq_along(week_lab)) {
  trace_hist(theta_mat[, w], sprintf("theta_%d", week_lab[w]))
  trace_hist(gamma_mat[, w], sprintf("gamma_%d", week_lab[w]), breaks = 2)
}

## ---- Point‑wise OR & P(gamma_t = 1) ----

alpha_mat <- theta_mat * gamma_mat
OR_mat    <- exp(alpha_mat)
OR_mean <- colMeans(OR_mat)
OR_ci   <- t(apply(OR_mat, 2, function(x) quantile(x, c(.025, .975))))

pdf("pointwise_or_cwvs.pdf", width = 7, height = 4)
plot(week_lab, OR_mean, type = "n",
     ylim = range(OR_ci),
     xlab = "week", ylab = "OR")
lines(week_lab, OR_mean, lty = 2)
points(week_lab, OR_mean, pch = 19, cex = 0.7)
segments(x0 = week_lab, y0 = OR_ci[, 1],
         x1 = week_lab, y1 = OR_ci[, 2])
abline(h = 1, lty = 3, col = "grey60")
dev.off()

pdf("posterior_inclusion_probability.pdf", width = 7, height = 4)
plot(week_lab, gamma_bar, type = "h", lwd = 3, ylim = c(0, 1),
     xlab = "week", ylab = "posterior P(gamma_t = 1)")
points(week_lab, gamma_bar, pch = 19)
abline(h = 0.5, lty = 2, col = "grey60")
dev.off()
