## rismbma: group selection for the Cox model by exact enumeration of the group-inclusion models.
##
## Design (each element follows a published source; see README):
##   * point-mass spike and Gaussian slab on whole groups (Xu and Ghosh 2015);
##   * the 2^G inclusion models are enumerated and their marginal likelihoods are obtained by the
##     Laplace approximation at the posterior mode (Raftery 1996; Bayesian model averaging for the
##     Cox model with the partial likelihood: Volinsky et al. 1997, Volinsky and Raftery 2000);
##   * beta-binomial prior on the number of included groups (Scott and Berger 2010);
##   * slab variance fixed by the user or chosen by empirical Bayes, i.e. by maximising the
##     marginal likelihood summed over models (George and Foster 2000);
##   * functional groups (principal component scores of one functional object) are scaled by one
##     common factor per group, so that high-order components with little variance are shrunk
##     more than the leading ones.

#' Centre the predictors and scale each group by one common factor
#'
#' Every column is centred. All columns of a group are divided by the same number, the root mean
#' variance of the columns of the group, so that the average variance within a group is one and
#' the relative variances of the columns are kept. For a group with one column this is ordinary
#' standardisation.
#' @param X numeric matrix (n x p).
#' @param g_id integer vector of length p with values 1..G.
#' @param center,scale optional values from a previous call, to apply a training scaling to new data.
#' @return The scaled matrix with attributes \code{center} and \code{scale}.
#' @export
group_scale = function(X, g_id, center = NULL, scale = NULL) {
  X = as.matrix(X)
  if (is.null(center)) center = colMeans(X)
  if (is.null(scale)) {
    v = apply(X, 2, stats::var)
    scale = sqrt(tapply(v, g_id, mean))[as.character(g_id)]
    scale = as.numeric(scale)
  }
  Z = sweep(sweep(X, 2, center, "-"), 2, scale, "/")
  attr(Z, "center") = center; attr(Z, "scale") = scale
  Z
}

# Cox partial log-likelihood (Breslow), gradient and Hessian, for data sorted by time.
# 'first' maps each subject to the first subject with the same time (Breslow risk sets).
.cox_parts = function(beta, X, status, first, hess = TRUE) {
  n = nrow(X); d = ncol(X)
  eta = if (d) as.vector(X %*% beta) else numeric(n)
  eta = eta - max(eta)
  w   = exp(eta)
  S0  = rev(cumsum(rev(w)))[first]
  ll  = sum(status * (eta - log(S0)))
  if (!d) return(list(ll = ll, grad = numeric(0), hess = matrix(0, 0, 0)))
  S1  = apply(w * X, 2, function(z) rev(cumsum(rev(z))))
  if (!is.matrix(S1)) S1 = matrix(S1, n, d)
  S1  = S1[first, , drop = FALSE]
  ev  = which(status == 1)
  E   = S1[ev, , drop = FALSE] / S0[ev]
  grad = colSums(X[ev, , drop = FALSE]) - colSums(E)
  H = NULL
  if (hess) {
    H = matrix(0, d, d)
    for (i in ev) {                      # few events, so the loop is cheap
      r  = first[i]:n
      Xr = X[r, , drop = FALSE]
      H  = H + crossprod(Xr, w[r] * Xr) / S0[i]
    }
    H = H - crossprod(E)                 # information matrix (positive semi-definite)
  }
  list(ll = ll, grad = grad, hess = H)
}

# Posterior mode and Laplace log marginal likelihood of one model with slab N(0, v I).
.laplace_model = function(X, status, first, v, maxit = 100, tol = 1e-8) {
  d = ncol(X)
  if (d == 0) return(list(logml = .cox_parts(numeric(0), X, status, first)$ll, beta = numeric(0), cov = matrix(0, 0, 0)))
  b = numeric(d); obj_old = -Inf
  for (it in seq_len(maxit)) {
    cp  = .cox_parts(b, X, status, first)
    obj = cp$ll - sum(b^2) / (2 * v)
    g   = cp$grad - b / v
    A   = cp$hess + diag(1 / v, d)
    step = solve(A, g)
    t = 1
    repeat {                              # step halving keeps the objective increasing
      bn   = b + t * step
      objn = .cox_parts(bn, X, status, first, hess = FALSE)$ll - sum(bn^2) / (2 * v)
      if (objn >= obj - 1e-12 || t < 1e-8) break
      t = t / 2
    }
    b = bn
    if (abs(objn - obj) < tol * (abs(obj) + 1)) break
  }
  cp = .cox_parts(b, X, status, first)
  A  = cp$hess + diag(1 / v, d)
  logml = cp$ll - sum(b^2) / (2 * v) - (d / 2) * log(v) - 0.5 * as.numeric(determinant(A, logarithm = TRUE)$modulus)
  list(logml = logml, beta = b, cov = solve(A))
}

#' Bayesian model averaging over groups for the Cox model
#'
#' Enumerates all 2^G models defined by including or excluding each group of predictors,
#' approximates the marginal likelihood of each by the Laplace method, and returns posterior
#' model probabilities, posterior inclusion probabilities of the groups and model-averaged
#' coefficients.
#'
#' @param time survival or censoring times.
#' @param status event indicator (1 = event, 0 = censored).
#' @param X predictor matrix (n x p). Use \code{\link{group_scale}} first.
#' @param g_id group of each column of \code{X}, integers 1..G (G at most 15).
#' @param v slab variance of a coefficient. Either one number, or a vector of candidate values
#'   among which the value with the largest marginal likelihood (summed over models) is chosen
#'   (empirical Bayes).
#' @param group_names optional names of the groups.
#' @param always indices of groups that are included in every model (e.g. clinical covariates).
#' @return A list with \code{pip} (posterior inclusion probability per group), \code{models}
#'   (one row per model: inclusion pattern, log marginal likelihood, posterior probability),
#'   \code{beta} (model-averaged posterior mean), \code{beta_sd}, \code{v}, \code{v_path}
#'   (empirical-Bayes criterion for each candidate), and \code{fits} (mode and covariance of
#'   each model).
#' @export
rism_bma = function(time, status, X, g_id, v = c(0.05, 0.1, 0.25, 0.5, 1, 2.5), group_names = NULL, always = integer(0)) {
  X = as.matrix(X); n = nrow(X); p = ncol(X); G = max(g_id)
  stopifnot(length(g_id) == p, length(time) == n, length(status) == n, G <= 15, all(sort(unique(g_id)) == seq_len(G)))
  if (is.null(group_names)) group_names = paste0("Group", seq_len(G))
  o = order(time); X = X[o, , drop = FALSE]; time_o = time[o]; status_o = status[o]
  first = match(time_o, time_o)
  free = setdiff(seq_len(G), always)
  pat  = as.matrix(expand.grid(rep(list(0:1), length(free))))
  gam  = matrix(0L, nrow(pat), G); gam[, free] = pat; gam[, always] = 1L
  colnames(gam) = group_names
  k = rowSums(gam[, free, drop = FALSE]); Gf = length(free)
  logprior = -log(Gf + 1) - lchoose(Gf, k)          # beta-binomial(1, 1) on the number of free groups
  run = function(vv) lapply(seq_len(nrow(gam)), function(m) {
    cols = which(gam[m, g_id] == 1L)
    .laplace_model(X[, cols, drop = FALSE], status_o, first, vv)
  })
  lse = function(a) { mx = max(a); mx + log(sum(exp(a - mx))) }
  crit = numeric(length(v)); allfits = vector("list", length(v))
  for (i in seq_along(v)) { allfits[[i]] = run(v[i]); crit[i] = lse(sapply(allfits[[i]], `[[`, "logml") + logprior) }
  best = which.max(crit); fits = allfits[[best]]
  logml = sapply(fits, `[[`, "logml")
  post  = exp(logml + logprior - lse(logml + logprior))
  pip   = colSums(gam * post); names(pip) = group_names
  bm = matrix(0, nrow(gam), p); bv = matrix(0, nrow(gam), p)
  for (m in seq_len(nrow(gam))) { cols = which(gam[m, g_id] == 1L)
    if (length(cols)) { bm[m, cols] = fits[[m]]$beta; bv[m, cols] = diag(fits[[m]]$cov) } }
  beta = colSums(bm * post)
  beta_sd = sqrt(pmax(colSums((bv + bm^2) * post) - beta^2, 0))
  names(beta) = names(beta_sd) = colnames(X)
  models = data.frame(gam, logml = logml, posterior = post, check.names = FALSE)
  models = models[order(-models$posterior), ]
  list(pip = pip, models = models, beta = beta, beta_sd = beta_sd, v = v[best],
       v_path = data.frame(v = v, log_marginal = crit), fits = fits, g_id = g_id, group_names = group_names)
}

#' Pre-validated risk scores
#'
#' Each subject's risk score is computed from a model fitted without that subject's fold; the
#' scaling of the predictors and the slab variance are chosen within the training folds.
#' @inheritParams rism_bma
#' @param X unscaled predictor matrix; \code{\link{group_scale}} is applied within each training set.
#' @param folds integer vector of fold labels.
#' @return Numeric vector of pre-validated risk scores.
#' @export
rism_prevalidate = function(time, status, X, g_id, folds, v = c(0.05, 0.1, 0.25, 0.5, 1, 2.5), always = integer(0)) {
  eta = numeric(length(time))
  for (k in sort(unique(folds))) {
    tr = folds != k
    Z  = group_scale(X[tr, , drop = FALSE], g_id)
    Zt = group_scale(X[!tr, , drop = FALSE], g_id, attr(Z, "center"), attr(Z, "scale"))
    f  = rism_bma(time[tr], status[tr], Z, g_id, v = v, always = always)
    eta[!tr] = as.vector(Zt %*% f$beta)
  }
  eta
}

#' Permutation reference for the inclusion probabilities
#'
#' The rows of the columns in \code{perm_cols} are permuted across subjects, which removes any
#' association between those predictors and survival while keeping everything else, and the
#' inclusion probabilities are recomputed. The result shows what PIPs the design produces when
#' the permuted predictors are unrelated to survival.
#' @inheritParams rism_bma
#' @param perm_cols columns of \code{X} to permute jointly.
#' @param B number of permutations.
#' @param seed random seed.
#' @return A B x G matrix of inclusion probabilities under permutation.
#' @export
rism_permute = function(time, status, X, g_id, perm_cols, B = 200, v = c(0.05, 0.1, 0.25, 0.5, 1, 2.5), always = integer(0), seed = 1) {
  set.seed(seed); n = nrow(X)
  t(sapply(seq_len(B), function(b) { Xp = X; Xp[, perm_cols] = X[sample(n), perm_cols]
    rism_bma(time, status, Xp, g_id, v = v, always = always)$pip }))
}

#' Stochastic search over group-inclusion models (MC3)
#'
#' For more than about 15 groups the 2^G models cannot be enumerated. This function runs the
#' Markov chain Monte Carlo model composition sampler of Madigan and York (1995): at each step
#' one group is proposed to enter or leave the model, and the move is accepted with the ratio of
#' the posterior model probabilities, which use the same Laplace marginal likelihoods and the
#' same beta-binomial model prior as \code{\link{rism_bma}}. The coefficients are integrated
#' out, so the chain moves on the indicators only and cannot get stuck through a mismatch
#' between an indicator and its coefficients. Marginal likelihoods are cached.
#'
#' @inheritParams rism_bma
#' @param v slab variance (one number).
#' @param n_iter number of iterations after burn-in; \code{burnin} iterations are discarded.
#' @param burnin number of initial iterations to discard.
#' @param init starting inclusion pattern (0/1 vector of length G); default all excluded.
#' @param seed random seed.
#' @return A list with \code{pip}, \code{visited} (distinct models visited, with log marginal
#'   likelihood and visit frequency), \code{accept} (acceptance rate) and \code{gamma} (the chain).
#' @export
rism_mc3 = function(time, status, X, g_id, v = 1, n_iter = 20000, burnin = 2000, init = NULL, group_names = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  X = as.matrix(X); G = max(g_id)
  if (is.null(group_names)) group_names = paste0("Group", seq_len(G))
  o = order(time); X = X[o, , drop = FALSE]; status_o = status[o]; first = match(time[o], time[o])
  cache = new.env()
  lpost = function(g) {
    key = paste(g, collapse = "")
    if (is.null(cache[[key]])) {
      cols = which(g[g_id] == 1L)
      cache[[key]] = .laplace_model(X[, cols, drop = FALSE], status_o, first, v)$logml - log(G + 1) - lchoose(G, sum(g))
    }
    cache[[key]]
  }
  g = if (is.null(init)) integer(G) else as.integer(init)
  cur = lpost(g); acc = 0; total = n_iter + burnin
  chain = matrix(0L, n_iter, G)
  for (it in seq_len(total)) {
    j = sample.int(G, 1); gp = g; gp[j] = 1L - gp[j]
    new = lpost(gp)
    if (log(stats::runif(1)) < new - cur) { g = gp; cur = new; if (it > burnin) acc = acc + 1 }
    if (it > burnin) chain[it - burnin, ] = g
  }
  pip = colMeans(chain); names(pip) = group_names
  keys = apply(chain, 1, paste, collapse = ""); tab = sort(table(keys), decreasing = TRUE)
  visited = data.frame(model = names(tab), frequency = as.numeric(tab) / n_iter,
                       log_posterior_unnormalised = sapply(names(tab), function(k) cache[[k]]), row.names = NULL)
  list(pip = pip, visited = visited, accept = acc / n_iter, gamma = chain, v = v)
}
