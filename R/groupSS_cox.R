## Cox model with grouped spike-and-slab prior
##
## Key fixes (relative to an earlier draft, kept here for context):
##   1. Zeta update now includes Cox likelihood
##   2. Group-adaptive MH proposal scaling
##   3. Tempering during burnin for better mixing
##   4. Sensible defaults for simulation (faster + correct)

#' MCMC for Cox model with grouped spike-and-slab prior
#'
#' @param time survival times
#' @param status event indicator (1=event, 0=censored)
#' @param X predictor matrix (n x p), should be standardized
#' @param g_id group membership vector, values in 1,...,G
#' @param beta_init initial beta. Defaults to small random values.
#' @param zeta_init initial zeta. Defaults to Bernoulli(0.5).
#' @param w_init initial w. Defaults to mean(zeta_init).
#' @param nusq_init initial nu^2. Defaults to 1.
#' @param v0 spike variance. Default 0.01.
#' @param v1 slab variance. Default 10.
#' @param a1 shape hyperparameter for nu^2. Default 0.5.
#' @param a2 rate hyperparameter for nu^2. Default 0.5.
#' @param alpha_w Beta hyperparameter for w. Default 1.
#' @param beta_w Beta hyperparameter for w. Default 1.
#' @param Nmcmc total MCMC iterations. Default 50000.
#' @param burnin burnin length. Default 20000.
#' @param thin thinning interval. Default 5.
#' @param mh_sd initial MH proposal SD per group. Default 0.15.
#' @param adapt_mh adaptive MH during burnin. Default TRUE.
#' @param temper use likelihood tempering during early burnin. Default TRUE.
#' @param verbose print progress. Default TRUE.
#' @param track_logpost track log posterior. Default TRUE.
#' @param seed random seed. Default NULL.
#' @return List with MCMC samples and diagnostics.
#' @export
groupSS_cox_final = function(time, status, X, g_id,
                             beta_init  = NULL, zeta_init = NULL,
                             w_init     = NULL, nusq_init = NULL,
                             v0         = 0.01, v1 = 10,
                             a1         = 0.5,  a2 = 0.5,
                             alpha_w    = 1,    beta_w = 1,
                             Nmcmc      = 50000, burnin = 20000, thin = 5,
                             mh_sd      = 0.15,
                             adapt_mh   = TRUE,
                             temper     = TRUE,
                             verbose    = TRUE,
                             track_logpost = TRUE,
                             seed       = NULL) {

  if (!is.null(seed)) set.seed(seed)

  nc = ncol(X)
  nr = nrow(X)
  G  = max(g_id)

  # --- Initialisation --------------------------------------------------------
  # Start beta near zero — flat start avoids biasing toward any particular
  # local mode under a low-information likelihood
  if (is.null(beta_init))  beta_init  = rnorm(nc, 0, 0.1)
  if (is.null(zeta_init))  zeta_init  = rbinom(G, 1, 0.5)
  if (is.null(w_init))     w_init     = min(max(mean(zeta_init), 1e-3), 1 - 1e-3)
  if (is.null(nusq_init))  nusq_init  = rep(1, nc)

  b    = beta_init
  zeta = zeta_init
  w    = w_init
  nusq = nusq_init

  # --- Storage ---------------------------------------------------------------
  n_save     = floor((Nmcmc - burnin) / thin)
  zeta_store = matrix(0L,  nrow = n_save, ncol = G)
  b_store    = matrix(0.0, nrow = n_save, ncol = nc)
  w_store    = numeric(n_save)
  log_post_store = if (track_logpost) numeric(n_save) else NULL

  # --- MH bookkeeping --------------------------------------------------------
  mh_sd_vec              = rep(mh_sd, G)
  accept_count_beta      = rep(0, G)
  accept_count_prev_beta = rep(0, G)
  accept_count_zeta      = rep(0, G)

  # --- Precompute sorted order -----------------------------------------------
  ord        = order(time)
  time_ord   = time[ord]
  status_ord = status[ord]
  X_ord      = X[ord, , drop = FALSE]

  save_idx = 1

  for (iter in 1:Nmcmc) {

    # ----- Likelihood tempering during early burnin --------------------------
    # Flattens the beta MH acceptance in the first 25% of burnin so the chain
    # explores group activation patterns before the full likelihood locks in.
    # For the zeta step, "tempering" scales the prior-only acceptance ratio
    # (there is no likelihood term in the zeta conditional given beta).
    if (temper && iter <= burnin) {
      temp = 0.1 + 0.9 * (iter / burnin)
    } else {
      temp = 1.0
    }

    # ----- 1. Update beta (group-wise Random-Walk MH) ------------------------
    for (g_idx in 1:G) {
      j_idx = which(g_id == g_idx)
      p_g   = length(j_idx)
      var_g = ifelse(zeta[g_idx] == 1, v1, v0)

      b_prop          = b
      b_prop[j_idx]   = b[j_idx] + rnorm(p_g, 0, mh_sd_vec[g_idx])

      log_prior_curr  = -0.5 * sum(b[j_idx]^2     / (nusq[j_idx] * var_g))
      log_prior_prop  = -0.5 * sum(b_prop[j_idx]^2 / (nusq[j_idx] * var_g))

      log_lik_curr    = cox_loglik_fast(b,      X_ord, time_ord, status_ord)
      log_lik_prop    = cox_loglik_fast(b_prop, X_ord, time_ord, status_ord)

      # Apply temperature to likelihood component only
      log_accept = temp * (log_lik_prop - log_lik_curr) +
                   (log_prior_prop - log_prior_curr)

      if (is.finite(log_accept) && log(runif(1)) < log_accept) {
        b = b_prop
        accept_count_beta[g_idx] = accept_count_beta[g_idx] + 1
      }
    }

    # ----- Adaptive MH tuning (every 200 iters during burnin) ----------------
    if (adapt_mh && iter <= burnin && iter %% 200 == 0) {
      for (g_idx in 1:G) {
        recent = (accept_count_beta[g_idx] - accept_count_prev_beta[g_idx]) / 200

        # Target 20-35% acceptance (standard for block MH)
        if      (recent < 0.10) mh_sd_vec[g_idx] = mh_sd_vec[g_idx] * 0.60
        else if (recent < 0.20) mh_sd_vec[g_idx] = mh_sd_vec[g_idx] * 0.85
        else if (recent > 0.45) mh_sd_vec[g_idx] = mh_sd_vec[g_idx] * 1.50
        else if (recent > 0.35) mh_sd_vec[g_idx] = mh_sd_vec[g_idx] * 1.20

        # Hard bounds
        mh_sd_vec[g_idx] = pmin(pmax(mh_sd_vec[g_idx], 0.005), 3.0)
      }
      accept_count_prev_beta = accept_count_beta
    }

    # ----- 2. Update zeta (deterministic flip MH) ----------------------------
    # Acceptance ratio is purely prior-driven; the Cox likelihood p(data|b)
    # cancels because beta does not change during this step (non-collapsed
    # sampler). The normalising constant of the Gaussian prior MUST be
    # included: omitting it biases the sampler toward the spike.
    #
    #   log_accept = [log p(b|zeta*) + log p(zeta*|w)]
    #              - [log p(b|zeta)  + log p(zeta|w)]
    for (g_idx in 1:G) {
      j_idx   = which(g_id == g_idx)
      p_g     = length(j_idx)

      var_curr = if (zeta[g_idx] == 1L) v1 else v0
      var_prop = if (zeta[g_idx] == 1L) v0 else v1   # flipped

      # Log p(b_g | zeta) including Gaussian normalising constant
      lp_beta_curr = -0.5 * p_g * log(2 * pi * var_curr) -
                     0.5 * sum(log(nusq[j_idx])) -
                     0.5 * sum(b[j_idx]^2 / (nusq[j_idx] * var_curr))

      lp_beta_prop = -0.5 * p_g * log(2 * pi * var_prop) -
                     0.5 * sum(log(nusq[j_idx])) -
                     0.5 * sum(b[j_idx]^2 / (nusq[j_idx] * var_prop))

      # Log p(zeta_g | w) — evaluate at current and proposed states directly;
      # do NOT use 1 - lp_zeta_curr (that is arithmetic on log-probs, not valid)
      lp_zeta_curr = if (zeta[g_idx] == 1L) log(w + 1e-10) else log(1 - w + 1e-10)
      lp_zeta_prop = if (zeta[g_idx] == 0L) log(w + 1e-10) else log(1 - w + 1e-10)

      log_accept = (lp_beta_prop + lp_zeta_prop) -
                   (lp_beta_curr + lp_zeta_curr)

      # During tempered burnin, also scale zeta moves to aid exploration
      if (temper && iter <= burnin) log_accept = log_accept * temp

      if (is.finite(log_accept) && log(runif(1)) < log_accept) {
        zeta[g_idx] = 1L - zeta[g_idx]
        accept_count_zeta[g_idx] = accept_count_zeta[g_idx] + 1
      }
    }

    # ----- 3. Update nu^2 (Gibbs, conjugate inverse-Gamma) ------------------
    var_vec = ifelse(zeta[g_id] == 1, v1, v0)
    nusq    = 1 / rgamma(nc,
                         shape = a1 + 0.5,
                         rate  = a2 + b^2 / (2 * var_vec))

    # ----- 4. Update w (Gibbs, conjugate Beta) -------------------------------
    n1 = sum(zeta)
    w  = rbeta(1, alpha_w + n1, beta_w + G - n1)

    # ----- 5. Store post-burnin samples --------------------------------------
    if (iter > burnin && (iter - burnin) %% thin == 0) {
      zeta_store[save_idx, ] = zeta
      b_store[save_idx, ]    = b
      w_store[save_idx]      = w

      if (track_logpost) {
        # lp_beta includes normalising constant for consistency with zeta MH
        lp_beta   = sum(-0.5 * log(2 * pi * var_vec) - 0.5 * b^2 / (nusq * var_vec))
        # Jacobian: if u=1/nu^2 ~ Gamma(a1,a2), then log p(nu^2) = log p(u) - 2*log(nu^2)
        lp_nusq   = sum(dgamma(1 / nusq, shape = a1, rate = a2, log = TRUE) - 2 * log(nusq))
        lp_w      = dbeta(w, alpha_w, beta_w, log = TRUE)
        lp_zeta   = sum(zeta * log(w + 1e-10) + (1 - zeta) * log(1 - w + 1e-10))
        log_lik   = cox_loglik_fast(b, X_ord, time_ord, status_ord)
        log_post_store[save_idx] = log_lik + lp_beta + lp_w + lp_zeta + lp_nusq
      }

      save_idx = save_idx + 1
    }

    if (verbose && iter %% 500 == 0) {
      if (iter > burnin) {
        pip_str = paste(sprintf("%.2f", colMeans(zeta_store[1:(save_idx - 1), , drop = FALSE])),
                        collapse = ", ")
      } else {
        pip_str = "burnin"
      }
      cat(sprintf("Iter %d/%d | temp=%.2f | w=%.3f | Active: %d | PIPs: %s\n",
                  iter, Nmcmc, temp, w, sum(zeta), pip_str))
    }
  }

  if (verbose) {
    cat("\n=== MCMC Complete ===\n")
    cat(sprintf("Final PIPs: %s\n",
                paste(sprintf("%.3f", colMeans(zeta_store)), collapse = ", ")))
    cat(sprintf("Mean w: %.3f\n", mean(w_store)))
    cat(sprintf("Beta acceptance rates: %s\n",
                paste(sprintf("%.2f", accept_count_beta / Nmcmc), collapse = ", ")))
    cat(sprintf("Zeta acceptance rates: %s\n",
                paste(sprintf("%.2f", accept_count_zeta / Nmcmc), collapse = ", ")))
    cat(sprintf("Final MH SDs: %s\n",
                paste(sprintf("%.4f", mh_sd_vec), collapse = ", ")))
  }

  list(
    beta              = b_store,
    zeta              = zeta_store,
    w                 = w_store,
    log_posterior     = log_post_store,
    x_cnames          = if (is.null(colnames(X))) paste0("X", seq_len(nc)) else colnames(X),
    g_id              = g_id,
    acceptance_rates_beta  = accept_count_beta  / Nmcmc,
    acceptance_rates_zeta  = accept_count_zeta  / Nmcmc,
    mh_sd_final       = mh_sd_vec,
    beta_init         = beta_init,
    zeta_init         = zeta_init,
    w_init            = w_init,
    seed              = seed,
    hyperparameters   = list(v0 = v0, v1 = v1, a1 = a1, a2 = a2,
                             alpha_w = alpha_w, beta_w = beta_w)
  )
}

#' Cox partial log-likelihood (Breslow approximation, sorted input)
#'
#' @param beta coefficient vector
#' @param X_ord predictor matrix, pre-sorted by \code{time_ord}
#' @param time_ord survival times, sorted ascending
#' @param status_ord event indicators, in the same order as \code{time_ord}
#' @return The Cox partial log-likelihood (Breslow tie handling).
#' @keywords internal
#' @noRd
cox_loglik_fast = function(beta, X_ord, time_ord, status_ord) {
  eta      = as.vector(X_ord %*% beta)
  # The partial likelihood is unchanged when a constant is added to the linear
  # predictor, so the maximum is subtracted for numerical stability. (Version
  # 0.1.0 clamped eta to [-100, 100] instead, which returned a wrong, flat
  # likelihood whenever predictors were not centred and scaled.)
  eta      = eta - max(eta)
  risk_set = rev(cumsum(rev(exp(eta))))
  # Breslow risk set for tied times: all subjects with time >= t
  if (anyDuplicated(time_ord)) risk_set = risk_set[match(time_ord, time_ord)]
  sum(status_ord * (eta - log(risk_set)))
}

#' Posterior inclusion probabilities
#'
#' @param fit Output of \code{\link{groupSS_cox_final}}.
#' @return A named numeric vector of posterior inclusion probabilities, one
#'   per group.
#' @export
compute_PIP = function(fit) {
  pip = colMeans(fit$zeta)
  names(pip) = paste0("Group", seq_len(ncol(fit$zeta)))
  pip
}

#' Posterior summary for beta
#'
#' @param fit Output of \code{\link{groupSS_cox_final}}.
#' @param prob Credible interval coverage. Default 0.95.
#' @return A data frame with one row per predictor: posterior mean, median,
#'   SD, and credible interval.
#' @export
summarize_beta = function(fit, prob = 0.95) {
  alpha = 1 - prob
  data.frame(
    Variable    = fit$x_cnames,
    Group       = fit$g_id,
    Post_Mean   = colMeans(fit$beta),
    Post_Median = apply(fit$beta, 2, median),
    Post_SD     = apply(fit$beta, 2, sd),
    CI_Lower    = apply(fit$beta, 2, quantile, probs = alpha / 2),
    CI_Upper    = apply(fit$beta, 2, quantile, probs = 1 - alpha / 2),
    stringsAsFactors = FALSE
  )
}

#' Trace plot for zeta (group inclusion indicators)
#'
#' @param fit Output of \code{\link{groupSS_cox_final}}.
#' @param group_names Optional character vector of group labels.
#' @return Invisibly \code{NULL}; called for its plotting side effect.
#' @export
plot_zeta_trace = function(fit, group_names = NULL) {
  G = ncol(fit$zeta)
  if (is.null(group_names)) group_names = paste0("Group ", seq_len(G))
  op = par(mfrow = c(ceiling(G / 2), 2), mar = c(3, 4, 2, 1))
  on.exit(par(op), add = TRUE)
  for (g in seq_len(G)) {
    plot(fit$zeta[, g], type = "l",
         main = group_names[g], ylab = "Indicator",
         xlab = "Post-burnin iteration", ylim = c(-0.1, 1.1))
    abline(h = 0.5, col = "red", lty = 2)
  }
}

#' Convergence diagnostics (effective sample size, log-posterior trace)
#'
#' Requires the \pkg{coda} package.
#'
#' @param fit Output of \code{\link{groupSS_cox_final}}.
#' @return Invisibly, a list with \code{beta_ess}, \code{pip_ess}, and
#'   \code{w_ess}. Also plots the log-posterior trace if
#'   \code{fit$log_posterior} is present, and warns if any effective sample
#'   size is below 100.
#' @export
check_convergence = function(fit) {
  if (!requireNamespace("coda", quietly = TRUE)) {
    stop("Package 'coda' is required for check_convergence(). ",
         "Install it with install.packages(\"coda\").", call. = FALSE)
  }
  ess = list(
    beta_ess = apply(fit$beta, 2, function(x) coda::effectiveSize(coda::mcmc(x))),
    pip_ess  = apply(fit$zeta, 2, function(x) coda::effectiveSize(coda::mcmc(x))),
    w_ess    = coda::effectiveSize(coda::mcmc(fit$w))
  )
  if (!is.null(fit$log_posterior)) {
    plot(fit$log_posterior, type = "l",
         main = "Log Posterior Trace",
         xlab = "Post-burnin iteration", ylab = "Log Posterior")
  }
  cat("\n=== Effective Sample Sizes ===\n")
  cat(sprintf("Min beta ESS : %.1f\n", min(ess$beta_ess)))
  cat(sprintf("Min PIP ESS  : %.1f\n", min(ess$pip_ess)))
  cat(sprintf("w ESS        : %.1f\n", ess$w_ess))
  if (min(c(ess$beta_ess, ess$pip_ess, ess$w_ess)) < 100)
    warning("ESS < 100. Consider longer chains or less thinning.")
  invisible(ess)
}
