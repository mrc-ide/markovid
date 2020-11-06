
# deploy.R
#
# Author: Bob Verity
# Date: 2020-11-05
#
# Purpose:
# Test package functions.

# ------------------------------------------------------------------

# load packages
library(ggplot2)
library(epitools)

set.seed(1)

# ------------------------------------------------------------------
# prepare data

# load dummy individual-level line list data
data_linelist <- readRDS(system.file("extdata", "dummy_indlevel.rds",
                                     package = "markovid",
                                     mustWork = TRUE))

# aggregate line list data
age_vec <- 0:100
data_agg <- aggregate_indlevel(df_data = data_linelist,
                               age_vec = age_vec)

# define age spline nodes
node_x <- seq(0, 100, 20)
n_node <- length(node_x)

# create final data list
data_list <- list(indlevel = data_agg,
                  max_indlevel_age = max(age_vec),
                  node_x = node_x)


# ------------------------------------------------------------------
# MCMC parameters

# create parameters dataframe
df_params <- rbind(data.frame(name = sprintf("p_AI_node%s", 1:n_node), min = -5, max = 5, init = 0),
                   data.frame(name = sprintf("p_AD_node%s", 1:n_node), min = -5, max = 5, init = 0),
                   data.frame(name = sprintf("p_ID_node%s", 1:n_node), min = -5, max = 5, init = 0),
                   data.frame(name = "m_AI", min = 0, max = 20, init = 5),
                   data.frame(name = "m_AD", min = 0, max = 20, init = 5),
                   data.frame(name = "m_AC", min = 0, max = 20, init = 5),
                   data.frame(name = "m_ID", min = 0, max = 20, init = 5),
                   data.frame(name = "m_IS", min = 0, max = 20, init = 5),
                   data.frame(name = "m_SC", min = 0, max = 20, init = 5),
                   data.frame(name = "s_AI", min = 0, max = 1, init = 0.5),
                   data.frame(name = "s_AD", min = 0, max = 1, init = 0.5),
                   data.frame(name = "s_AC", min = 0, max = 1, init = 0.5),
                   data.frame(name = "s_ID", min = 0, max = 1, init = 0.5),
                   data.frame(name = "s_IS", min = 0, max = 1, init = 0.5),
                   data.frame(name = "s_SC", min = 0, max = 1, init = 0.5)
)


# ------------------------------------------------------------------
# Run MCMC

# define MCMC parameters
burnin <- 1e3
samples <- 1e4
chains <- 1

# run MCMC
set.seed(1)
t0 <- Sys.time()
mcmc <- markovid::run_mcmc(data_list = data_list,
                           df_params = df_params,
                           burnin = burnin,
                           samples = samples,
                           chains = chains)

print(Sys.time() - t0)


# basic exploratory plots
plot_credible(mcmc)
#plot_par(mcmc, show = "p_AI", phase = "both")


# ------------------------------------------------------------------
# Get posterior splines and quantiles

# subset to sampling iterations of desired parameters only
mcmc_samples <- subset(mcmc$output,
                       stage == "sampling" & rung == "rung1",
                       select = -c(chain, rung, iteration, stage, logprior, loglikelihood))


# get spline quantiles for transition probabilities
mcmc_samples_p_AI <- as.matrix(mcmc_samples[, sprintf("p_AI_node%s", seq_along(node_x)), drop = FALSE])
p_AI_spline <- get_spline(mcmc_samples_p_AI, node_x, age_vec)
p_AI_quantile <- get_spline_quantiles(p_AI_spline)

mcmc_samples_p_AD <- as.matrix(mcmc_samples[, sprintf("p_AD_node%s", seq_along(node_x)), drop = FALSE])
p_AD_spline <- get_spline(mcmc_samples_p_AD, node_x, age_vec)
p_AD_quantile <- get_spline_quantiles(p_AD_spline)

mcmc_samples_p_ID <- as.matrix(mcmc_samples[, sprintf("p_ID_node%s", seq_along(node_x)), drop = FALSE])
p_ID_spline <- get_spline(mcmc_samples_p_ID, node_x, age_vec)
p_ID_quantile <- get_spline_quantiles(p_ID_spline)

p_OD_spline <- (1 - p_AI_spline)*p_AD_spline + p_AI_spline*p_ID_spline
p_OD_quantile <- get_spline_quantiles(p_OD_spline)


# ------------------------------------------------------------------
# Get 95% exact binomial intervals from raw data

# prob. admission to ICU
df_quantiles_p_AI <- get_data_quantiles_p(age = data_linelist$age,
                                          max_age = max(age_vec),
                                          outcome1 = (data_linelist$icu == TRUE),
                                          outcome2 = (data_linelist$icu == FALSE))

# prob. general ward to death
df_quantiles_p_AD <- get_data_quantiles_p(age = data_linelist$age,
                                          max_age = max(age_vec),
                                          outcome1 = (data_linelist$icu == FALSE) & (data_linelist$final_outcome == "death"),
                                          outcome2 = (data_linelist$icu == FALSE) & (data_linelist$final_outcome == "discharge"))

# prob. ICU ward to death
df_quantiles_p_ID <- get_data_quantiles_p(age = data_linelist$age,
                                          max_age = max(age_vec),
                                          outcome1 = (data_linelist$icu == TRUE) & (data_linelist$final_outcome == "death"),
                                          outcome2 = (data_linelist$icu == TRUE) & (data_linelist$final_outcome == "discharge"))

# overall prob. death
df_quantiles_p_OD <- get_data_quantiles_p(age = data_linelist$age,
                                          max_age = max(age_vec),
                                          outcome1 = (data_linelist$final_outcome == "death"),
                                          outcome2 = (data_linelist$final_outcome == "discharge"))

# ------------------------------------------------------------------
# Plot posteriors


# produce spline plots
plot_p_AI <- plot_spline_quantiles(df_spline_quantile = p_AI_quantile,
                                   df_data_quantile = df_quantiles_p_AI,
                                   title = "p_AI", ylim = c(0,1))
plot_p_AD <- plot_spline_quantiles(df_spline_quantile = p_AD_quantile,
                                   df_data_quantile = df_quantiles_p_AD,
                                   title = "p_AD", ylim = c(0,1))
plot_p_ID <- plot_spline_quantiles(df_spline_quantile = p_ID_quantile,
                                   df_data_quantile = df_quantiles_p_ID,
                                   title = "p_ID", ylim = c(0,1))
plot_p_OD <- plot_spline_quantiles(df_spline_quantile = p_OD_quantile,
                                   df_data_quantile = df_quantiles_p_OD,
                                   title = "p_OD", ylim = c(0,1))

# produce combined plot of transition splines
cowplot::plot_grid(plot_p_AI, plot_p_AD, plot_p_ID, plot_p_OD)

# plot durations
plot_credible(mcmc, show = c("m_AI", "m_AD", "m_AC", "m_ID", "m_IS", "m_SC")) +
  ggplot2::ylim(0, 20)




