# rism

**R**adiomic feature **I**ntegration with **S**pike-and-slab for **M**ultimodal
survival analysis.

An R package fitting a Cox proportional hazards model with a grouped
spike-and-slab prior over predictor groups (e.g. radiomic feature families,
or feature groups from different data modalities), via a
Metropolis-Hastings / Gibbs MCMC sampler with adaptive proposal scaling and
likelihood tempering during burn-in.

## Installation

```r
# install.packages("remotes")
remotes::install_github("sagnikbhadury/rism")
```

## Usage

```r
library(rism)

fit <- groupSS_cox_final(time = surv_time, status = event, X = X_std, g_id = group_ids,
                          Nmcmc = 50000, burnin = 20000, thin = 5, seed = 1)

compute_PIP(fit)          # posterior inclusion probability per group
summarize_beta(fit)       # posterior mean/median/SD/CI per predictor
plot_zeta_trace(fit)      # trace plots of group inclusion indicators
check_convergence(fit)    # effective sample sizes (requires the coda package)
```

See `?groupSS_cox_final` for full documentation.
