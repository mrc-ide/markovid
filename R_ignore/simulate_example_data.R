
# simulate_example_data.R
#
# Author: Bob Verity
# Date: 2020-11-05
#
# Purpose:
# Produce simulated individual-level data. Save to file inside package.

# ------------------------------------------------------------------

set.seed(1)


# define scalar parameters
params_scalar <- c(s_AI = 0.9,
                   s_AD = 0.9,
                   s_AC = 0.9,
                   s_ID = 0.9,
                   s_I1S = 0.9,
                   s_I2S = 0.9,
                   s_SD = 0.9,
                   s_SC = 0.9)

# define age-varying params
params_age <- list(p_AI = 0.2 * dbeta(seq(0, 1, l = 101), 4, 4),
                   p_AD = seq(0, 0.5, l = 101),
                   p_ID = seq(0, 0.8, l = 101),
                   p_SD = seq(0, 0.5, l = 101),
                   m_AI = rep(3, 101),
                   m_AD = rep(10, 101),
                   m_AC = seq(5, 15, l = 101),
                   m_ID = seq(15, 5, l = 101),
                   m_I1S = 10*dbeta(seq(0.1, 0.9, l = 101), 2, 2),
                   m_I2S = 5*dbeta(seq(0.1, 0.9, l = 101), 2, 2),
                   m_SD = seq(10, 5, l = 101),
                   m_SC = seq(5, 10, l = 101))

# draw dates of admission and ages, and define censoring date
n <- 1e4
df_sim <- data.frame(date_admission = floor(rlnorm(n, meanlog = 4, sdlog = 0.2)),
                     date_censor = 100,
                     age = sample(0:100, n, replace = TRUE))

# simulate from basic hospital progression model
df_sim <- sim_indlevel(params_scalar = params_scalar,
                       params_age = params_age,
                       age_vec = 0:100,
                       df_sim = df_sim)

# save example data to package
saveRDS(df_sim, "inst/extdata/dummy_indlevel.rds")
