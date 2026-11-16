# rism

**R**adiomic feature **I**ntegration with **S**pike-and-slab for **M**ultimodal
survival analysis.

`rism` selects *groups* of predictors in a Cox proportional hazards model. A
group is a set of predictors that belong together, for example the principal
component (PC) scores of a tumor boundary, or the PC scores of the
voxel-intensity density of one MRI sequence. For every group the package
returns a posterior inclusion probability (PIP), and for every predictor a
posterior summary of its coefficient.

Two methods are provided.

| | **Group BMA** (`rism_bma`) | **Grouped SS** (`groupSS_cox_final`) |
|---|---|---|
| Idea | Bayesian model averaging over the 2^G models that include or exclude each group | MCMC with a continuous grouped spike-and-slab prior |
| Spike | point mass: an excluded group contributes exactly zero | narrow Gaussian: an excluded group keeps small coefficients |
| Computation | exact enumeration with Laplace marginal likelihoods; stochastic search (`rism_mc3`) when there are many groups | random-walk Metropolis-Hastings within Gibbs |
| Tuning | slab variance chosen by empirical Bayes | spike and slab variances set by the user |
| Speed (61 subjects, 42 predictors, 7 groups) | about 3 seconds | about 45 seconds per chain |
| Recommended | **yes, as the default** | for comparison, or when a fully sampled posterior of the coefficients is wanted |

The two methods belong to the same class of models: a Cox proportional hazards
model with a spike-and-slab prior on groups of coefficients. They differ in the
spike (a point mass or a narrow Gaussian) and in the computation (enumeration
with the coefficients integrated out, or MCMC).

In the simulation study of the accompanying manuscript (450 data sets, six
methods) Group BMA had the highest F1 for group selection in 10 of 15 scenarios
and the smallest coefficient error in 14.

## Contents

- [Installation](#installation)
- [Preparing your data](#preparing-your-data)
- [Group BMA](#group-bma)
- [Grouped SS](#grouped-ss)
- [Example data: Wisconsin Prognostic Breast Cancer](#example-data-wisconsin-prognostic-breast-cancer)
- [Results in the accompanying manuscript](#results-in-the-accompanying-manuscript)
- [How much data the methods need](#how-much-data-the-methods-need)
- [Function reference](#function-reference)
- [Version notes](#version-notes)
- [References](#references)

## Installation

```r
# install.packages("remotes")
remotes::install_github("sagnikbhadury/rism")
```

The package depends only on base R (`stats`, `graphics`). `coda` is needed for
`check_convergence()`.

## Preparing your data

1. **One row per subject.** `time` is the follow-up time and `status` is 1 for
   an event and 0 for a censored subject. Missing values are not handled;
   impute or remove them first.
2. **Group vector.** `g_id` has one entry per column of `X`, with integer
   values `1, ..., G` and no gaps. A single predictor can be its own group.
3. **Scaling.**
   - For Group BMA use `group_scale(X, g_id)`. Every column is centred and all
     columns of a group are divided by one common number, so the relative
     variances of the PC scores of a functional predictor are kept and
     high-order components, which carry little variance, are shrunk more.
   - For Grouped SS use `scale(X)` (every column to mean 0 and variance 1).
   - Always centre and scale clinical covariates such as age. A coefficient is
     then a log hazard ratio per standard deviation.
4. **Keep groups small relative to the number of events.** A group of 23 scores
   cannot be assessed with 14 events (see below).

## Group BMA

### The model

For subject `i` the hazard is `lambda_0(t) exp(x_i' beta)` and inference uses
the Cox partial likelihood (Breslow form for tied times). Let `gamma_g` be 1 if
group `g` is in the model and 0 otherwise.

```
beta_j | gamma_g = 1  ~ N(0, v)        for j in group g   (slab)
beta_j | gamma_g = 0  = 0                                  (point-mass spike)
number of included groups ~ Beta-binomial(1, 1)
```

For each of the 2^G inclusion patterns the marginal likelihood is obtained by
the Laplace approximation at the posterior mode. The posterior probability of a
model is proportional to its marginal likelihood times its prior, and the PIP of
a group is the total probability of the models that contain it. A group is
selected when its PIP exceeds 0.5 (the median probability model). The slab
variance `v` is either given or chosen among candidate values by empirical
Bayes, i.e. by maximising the marginal likelihood summed over the models.

### The algorithm

`rism_bma` carries out the following steps (Algorithm 1 of the manuscript).

1. Sort the subjects by time. (Scale the predictors beforehand, see
   [Preparing your data](#preparing-your-data).)
2. For each candidate slab variance `v` and each of the 2^G models: find the
   posterior mode of the included coefficients by Newton's method (maximising
   the log partial likelihood minus `sum(beta^2) / (2 v)`), and compute the
   Laplace approximation to the log marginal likelihood of the model.
3. For each `v`, sum the marginal likelihoods over the models, weighted by the
   model prior. Keep the `v` with the largest sum (empirical Bayes).
4. At that `v`, normalise marginal likelihood times prior over the models to
   obtain the posterior model probabilities.
5. The PIP of a group is the sum of the probabilities of the models that
   contain it; the coefficients are averaged over the models, with zeros for
   excluded groups.

No step uses random numbers, so the result is the same on every run and there
is no chain to monitor. The cost is (number of candidate `v`) x 2^G Newton fits.

### Quick start

```r
library(rism)

set.seed(1)
n    <- 80
X    <- matrix(rnorm(n * 6), n)
g_id <- c(1, 1, 2, 3, 3, 3)                       # three groups
beta <- c(0.8, -0.6, 0.4, 0, 0, 0)                # groups 1 and 2 active
t_event <- rexp(n, exp(scale(X) %*% beta))
t_cens  <- rexp(n, 0.4)
time    <- pmin(t_event, t_cens)
status  <- as.numeric(t_event <= t_cens)          # 59 events

Z   <- group_scale(X, g_id)
fit <- rism_bma(time, status, Z, g_id)

fit$pip                 # posterior inclusion probability of each group
#> Group1 Group2 Group3
#>  1.000  0.915  0.072
head(fit$models)        # models ordered by posterior probability
fit$beta                # model-averaged coefficients (posterior means)
fit$v                   # slab variance chosen by empirical Bayes
```

### Reading the result

- `fit$pip`: PIP per group.
- `fit$models`: one row per model with its inclusion pattern, log marginal
  likelihood and posterior probability, most probable first.
- `fit$beta`, `fit$beta_sd`: model-averaged posterior mean and standard
  deviation of each coefficient (zero where a group is always excluded).
- `fit$v_path`: the empirical-Bayes criterion for each candidate slab variance.
- A risk score for each subject is `Z %*% fit$beta`.

### Options

```r
rism_bma(time, status, Z, g_id, v = 1)              # fixed slab variance
rism_bma(time, status, Z, g_id, always = c(1, 2))   # groups 1 and 2 in every model
rism_mc3(time, status, Z, g_id, v = 1, n_iter = 20000, burnin = 2000, seed = 1)
```

`rism_bma` enumerates up to 15 groups (32,768 models). For more groups use
`rism_mc3`, which runs a stochastic search over the same models with the same
marginal likelihoods; run it from several starting patterns (`init`) and check
that the PIPs agree.

### Checking the result

- **Out-of-sample prediction.** `rism_prevalidate()` returns a risk score for
  each subject from a model fitted without that subject's fold, with the scaling
  and the slab variance chosen inside the training folds. Compare its
  concordance with that of a simple reference model.

  ```r
  folds <- sample(rep_len(1:10, n))
  eta   <- rism_prevalidate(time, status, X, g_id, folds)
  survival::concordance(survival::Surv(time, status) ~ eta, reverse = TRUE)
  ```

- **Permutation reference.** `rism_permute()` permutes the rows of chosen
  columns, which removes any association between them and survival, and returns
  the PIPs for each permutation. A PIP is evidence only if it is large relative
  to this reference.

  ```r
  ref <- rism_permute(time, status, Z, g_id, perm_cols = 3:6, B = 200)
  apply(ref, 2, quantile, 0.95)      # 95th percentile of the PIP under no association
  ```

- **Stochastic-search check and trace plots.** `rism_bma` has no chain, but the
  same posterior over models can be sampled with `rism_mc3`. Running a few
  chains from different starting models gives trace plots and an independent
  check of the enumerated PIPs.

  ```r
  G      <- max(g_id)
  starts <- list(integer(G), rep(1L, G))            # no group, all groups
  chains <- lapply(seq_along(starts), function(k)
    rism_mc3(time, status, Z, g_id, v = fit$v, n_iter = 20000, burnin = 2000,
             init = starts[[k]], seed = k))

  sapply(chains, `[[`, "pip")                       # compare with fit$pip
  matplot(sapply(chains, `[[`, "logpost"), type = "l", lty = 1,
          xlab = "Iteration", ylab = "Log posterior of the model")
  running <- sapply(chains, function(f) cumsum(f$gamma[, 1]) / seq_len(nrow(f$gamma)))
  matplot(running, type = "l", lty = 1, ylim = c(0, 1),
          xlab = "Iteration", ylab = "Running PIP of group 1")
  abline(h = fit$pip[1], lty = 2)                   # enumerated value
  ```

- **Accuracy of the Laplace approximation.** It is accurate to a few tenths of a
  log unit for models with up to about 20 coefficients at 14 events, and
  underestimates the marginal likelihood of larger models by 1 to 2 log units.

## Grouped SS

### The model

```
beta_j | nu_j^2, zeta_g ~ N(0, nu_j^2 * (zeta_g * v1 + (1 - zeta_g) * v0))
nu_j^2                  ~ Inverse-Gamma(a1, a2)
zeta_g | w              ~ Bernoulli(w)
w                       ~ Beta(alpha_w, beta_w)
```

`zeta_g = 1` puts the whole group in the slab (variance scale `v1`) and
`zeta_g = 0` in the spike (`v0 << v1`). The sampler updates the coefficients of
one group at a time by random-walk Metropolis-Hastings, flips each `zeta_g` by a
Metropolis-Hastings step and draws `nu_j^2` and `w` from their full
conditionals. The proposal scale is adapted and the likelihood is tempered
during burn-in only.

### Quick start

```r
Xs  <- scale(X)
fit <- groupSS_cox_final(time, status, Xs, g_id,
                         v0 = 0.01, v1 = 2.5, a1 = 5, a2 = 5,
                         Nmcmc = 30000, burnin = 10000, thin = 5,
                         seed = 1, verbose = FALSE)

compute_PIP(fit)          # posterior inclusion probability per group
summarize_beta(fit)       # posterior mean, median, SD and credible interval
plot_zeta_trace(fit)      # trace of the group indicators
check_convergence(fit)    # effective sample sizes (needs coda)
```

### Using it responsibly

- **The PIPs depend on `v0` and `v1`**, strongly so when events are few. Fix the
  slab, fit a decreasing sequence of spike variances and keep the value with
  the best pre-validated partial log-likelihood (Tang et al., 2017; Verweij and
  van Houwelingen, 1993). Do not choose them by which groups come out selected.
- **Run at least four chains** from dispersed starting values (`beta_init`,
  `zeta_init`, different seeds), including all groups in and all groups out, and
  require rank-normalised R-hat below 1.01 and effective sample size above 400
  (Vehtari et al., 2021).
- **Always centre and scale `X`.**

## Example data: Wisconsin Prognostic Breast Cancer

The package includes the Wisconsin Prognostic Breast Cancer data (`wpbc`): 198
patients with invasive breast cancer, time to recurrence (47 recurrences), 30
features of the cell nuclei on a digitised image of a fine needle aspirate, tumor
size and the number of positive lymph nodes. The image features come in ten
natural groups (radius, texture, perimeter, area, smoothness, compactness,
concavity, concave points, symmetry, fractal dimension), each with its mean,
standard error and worst value. `wpbc_groups()` returns the complete cases with
the predictors arranged by group.

```r
library(rism)
data(wpbc)                       # the data frame, see ?wpbc
w <- wpbc_groups()               # time, status, X (194 x 32), g_id, group_names

fit <- rism_bma(w$time, w$status, scale(w$X), w$g_id, group_names = w$group_names)
round(fit$pip, 2)                # 12 groups, 4,096 models; about a minute per candidate
                                 # slab variance (six by default); use v = 0.25 for one

ss <- groupSS_cox_final(w$time, w$status, scale(w$X), w$g_id,
                        v0 = 0.05, v1 = 2.5, a1 = 5, a2 = 5, seed = 1, verbose = FALSE)
round(compute_PIP(ss), 2)        # about 1.5 minutes; choose v0 as described under Grouped SS
```

Here the three summaries of a characteristic are on different scales, so the
columns are standardised with `scale()`; `group_scale()` is meant for groups of
principal component scores.

**Source and licence of the data.** Wolberg, W., Street, W. and Mangasarian, O.
(1995). *Breast Cancer Wisconsin (Prognostic)* [Dataset]. UCI Machine Learning
Repository. https://doi.org/10.24432/C5GK50. The data are distributed under the
Creative Commons Attribution 4.0 International (CC BY 4.0) licence; the copy in
this package was taken from the R package `TH.data`. If you use the data, cite
the data set above and Street, Mangasarian and Wolberg (1995), "An inductive
learning approach to prognostic prediction", *Proceedings of the Twelfth
International Conference on Machine Learning*, 522-530. The MIT licence of this
package applies to the code, not to the data.

## Results in the accompanying manuscript

The manuscript applies both methods to two cohorts.

**Lower-grade glioma (61 patients, 14 deaths; 42 predictors in 7 groups).** The
groups are age, the Karnofsky performance score (KPS), the shape of the tumor
boundary (23 PC scores) and the intensity densities of four MRI sequences. Group
BMA retains KPS (PIP 0.99); age is borderline (0.46). No imaging group reaches
0.5 (T1Gd 0.24, T1 0.10, the others below 0.03), and under permutation of the
imaging scores the observed PIP of T1Gd is exceeded 18% of the time. With 14
deaths, a prognostic contribution of the imaging groups is not established. The
image-derived scores of this cohort are not distributed with the package.

**Breast cancer recurrence (`wpbc`; 194 complete cases, 46 recurrences; 32
predictors in 12 groups).** Group BMA gives the largest PIPs among the image
features to the three groups that measure the size of the nuclei (perimeter
0.59, area 0.50, radius 0.49), the only groups above the 95th percentile of the
permutation reference. These three are nearly collinear, and the probability
that at least one of them is in the model is 0.97 (prior 0.75). Grouped SS, with
the spike variance chosen by pre-validation, gives every group a PIP below 0.04.

**Convergence.** Group BMA is computed by enumeration. As a check, four
`rism_mc3` chains per cohort (20,000 iterations each) gave PIPs within 0.02
(glioma) and 0.04 (breast cancer) of the enumerated values. The four Grouped SS
chains per cohort had rank-normalised R-hat of at most 1.002. Trace plots of
both are in the supplement of the manuscript.

## How much data the methods need

Group selection needs events. In the simulation study (22 predictors in six
groups, 5 to 27 events):

- Group BMA selected an inactive group in 8% to 42% of cases across the design.
- Grouped SS with a diffuse slab selected an inactive group in 66% to 97% of
  cases at 5 to 9 events, and in 1% to 14% once about one event per predictor
  was available.

Before interpreting a PIP from a small cohort, look at the permutation
reference, and simulate survival times on your own design matrix with a planted
effect to see whether it is recovered at your number of events.

## Function reference

| Function | Purpose |
|---|---|
| `group_scale(X, g_id)` | centre the columns; one common scale per group |
| `rism_bma(time, status, X, g_id, v, group_names, always)` | Group BMA by enumeration |
| `rism_mc3(time, status, X, g_id, v, n_iter, burnin, init, seed)` | Group BMA by stochastic search |
| `rism_prevalidate(time, status, X, g_id, folds, v, always)` | pre-validated risk scores |
| `rism_permute(time, status, X, g_id, perm_cols, B, v, always, seed)` | permutation reference for the PIPs |
| `groupSS_cox_final(time, status, X, g_id, ...)` | Grouped SS sampler |
| `wpbc`, `wpbc_groups()` | example data and its group structure |
| `compute_PIP(fit)`, `summarize_beta(fit)`, `plot_zeta_trace(fit)`, `check_convergence(fit)` | summaries and diagnostics of a Grouped SS fit |

Each function has a help page (`?rism_bma`, `?groupSS_cox_final`, ...).

## Version notes

- **0.2.2.** `rism_mc3` also returns `logpost`, the log posterior (up to a
  constant) of the model at each iteration, for trace plots. README: the
  algorithm of Group BMA, the stochastic-search check, and the results of the
  two applications of the manuscript.
- **0.2.1.** Added the Wisconsin Prognostic Breast Cancer data (`wpbc`, `wpbc_groups()`).
- **0.2.0.** Added Group BMA (`rism_bma`, `rism_mc3`, `group_scale`,
  `rism_prevalidate`, `rism_permute`).
- **0.1.1.** The partial likelihood of Grouped SS no longer clamps the linear
  predictor to [-100, 100]; it subtracts the maximum instead, which is exact.
  Version 0.1.0 returned a wrong (flat) likelihood when predictors were not
  centred and scaled, for example age in years. Tied times now use the Breslow
  risk set. Results obtained with 0.1.0 on unscaled predictors should be rerun.
- **0.1.0.** First version (Grouped SS only).

## References

- George, E. I. and Foster, D. P. (2000). Calibration and empirical Bayes
  variable selection. *Biometrika* 87, 731-747.
- Madigan, D. and York, J. (1995). Bayesian graphical models for discrete data.
  *International Statistical Review* 63, 215-232.
- Raftery, A. E. (1996). Approximate Bayes factors and accounting for model
  uncertainty in generalised linear models. *Biometrika* 83, 251-266.
- Scott, J. G. and Berger, J. O. (2010). Bayes and empirical-Bayes multiplicity
  adjustment in the variable-selection problem. *Annals of Statistics* 38,
  2587-2619.
- Tang, Z., Shen, Y., Zhang, X. and Yi, N. (2017). The spike-and-slab lasso Cox
  model for survival prediction and associated genes detection.
  *Bioinformatics* 33, 2799-2807.
- Vehtari, A., Gelman, A., Simpson, D., Carpenter, B. and Buerkner, P.-C.
  (2021). Rank-normalization, folding, and localization: an improved R-hat for
  assessing convergence of MCMC. *Bayesian Analysis* 16, 667-718.
- Verweij, P. J. M. and van Houwelingen, H. C. (1993). Cross-validation in
  survival analysis. *Statistics in Medicine* 12, 2305-2314.
- Volinsky, C. T. and Raftery, A. E. (2000). Bayesian information criterion for
  censored survival models. *Biometrics* 56, 256-262.
- Xu, X. and Ghosh, M. (2015). Bayesian variable selection and estimation for
  group lasso. *Bayesian Analysis* 10, 909-936.

## License

MIT.
