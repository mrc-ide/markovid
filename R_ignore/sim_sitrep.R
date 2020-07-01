
set.seed(1)

# logit transform
logit <- function(x, z = 1) {
  -log((z - x) / x)
}

# define true parameters
true_params <- c(m_AL = 0.1,
                 p_AI = 0.5,
                 p_AD = 0.5,
                 m_AD = 10,
                 s_AD = 0.8,
                 m_AI = 3,
                 s_AI = 1,
                 m_AC = 10,
                 s_AC = 0.8,
                 p_ID = 0.8,
                 m_ID = 10,
                 s_ID = 0.8,
                 m_IS = 10,
                 s_IS = 0.8,
                 m_SC = 5,
                 s_SC = 1,
                 c_AD = 1,
                 c_AC = 1,
                 c_ID = 1,
                 c_IS = 1,
                 scale_p_AI1 = 1,
                 scale_p_AD1 = 1,
                 scale_p_ID1 = 1)

# convert to node values
n_node <- 7

p_AI_node <- rep(logit(true_params["p_AI"]), n_node)
names(p_AI_node) <- sprintf("p_AI_node%s", seq_len(n_node))

p_AD_node <- rep(logit(true_params["p_AD"]), n_node)
names(p_AD_node) <- sprintf("p_AD_node%s", seq_len(n_node))

p_ID_node <- rep(logit(true_params["p_ID"]), n_node)
names(p_ID_node) <- sprintf("p_ID_node%s", seq_len(n_node))

m_AC_node <- rep(logit(true_params["m_AC"], 20), n_node)
names(m_AC_node) <- sprintf("m_AC_node%s", seq_len(n_node))

true_params <- c(true_params, p_AI_node, p_AD_node, p_ID_node, m_AC_node)


# define true admissions curve
t_vec <- -15:40
true_admissions_curve <- ceiling(1e3 * dnorm(t_vec, mean = 12, sd = 5))

# simulate progression
sim <- data.frame(date_admission = rep(t_vec, times = true_admissions_curve))
n_sim <- nrow(sim)
sim$date_test <- sim$date_admission + floor(rgamma(n_sim, shape = 1, scale = true_params["m_AL"]))
sim$icu <- runif(n_sim) < true_params["p_AI"]
w <- which(!sim$icu)
sim$final_outcome <- NA
sim$final_outcome[w] <- sample(c("death", "discharge"), length(w), replace = TRUE, prob = c(true_params["p_AD"], 1 - true_params["p_AD"]))
w2 <- which(!sim$icu & sim$final_outcome == "death")
sim$date_final_outcome <- NA
sim$date_final_outcome[w2] <- sim$date_admission[w2] + 
                              floor(rgamma(length(w2), shape = 1/true_params["s_AD"]^2, scale = true_params["m_AD"]*true_params["s_AD"]^2))
w2 <- which(!sim$icu & sim$final_outcome == "discharge")
sim$date_final_outcome[w2] <- sim$date_admission[w2] + 
                              floor(rgamma(length(w2), shape = 1/true_params["s_AC"]^2, scale = true_params["m_AC"]*true_params["s_AC"]^2))
w <- which(sim$icu)
sim$date_icu <- NA
sim$date_icu[w] <- sim$date_admission[w] +
                   floor(rgamma(length(w), shape = 1/true_params["s_AI"]^2, scale = true_params["m_AI"]*true_params["s_AI"]^2))
sim$final_outcome[w] <- sample(c("death", "discharge"), length(w), replace = TRUE, prob = c(true_params["p_ID"], 1 - true_params["p_ID"]))
w2 <- which(sim$icu & sim$final_outcome == "death")
sim$date_final_outcome[w2] <- sim$date_icu[w2] + 
                              floor(rgamma(length(w2), shape = 1/true_params["s_ID"]^2, scale = true_params["m_ID"]*true_params["s_ID"]^2))
w2 <- which(sim$icu & sim$final_outcome == "discharge")
sim$date_stepdown <- NA
sim$date_stepdown[w2] <- sim$date_icu[w2] +
                         floor(rgamma(length(w2), shape = 1/true_params["s_IS"]^2, scale = true_params["m_IS"]*true_params["s_IS"]^2))
sim$date_final_outcome[w2] <- sim$date_stepdown[w2] + 
                              floor(rgamma(length(w2), shape = 1/true_params["s_SC"]^2, scale = true_params["m_SC"]*true_params["s_SC"]^2))

sim$age <- sample(0:105, n_sim, replace = TRUE)

# initialise sitrep object over all age groups
age_lower <- c(0, 6, 18, 65, 85)
age_upper <- c(age_lower[-1] - 1, 110)
nt <- length(t_vec)
sitrep <- replicate(5, data.frame(t = t_vec,
                                  obs_admissions = NA,
                                  deaths_general = NA,
                                  discharges_general = NA,
                                  stepup = NA,
                                  deaths_critical = NA,
                                  stepdown = NA,
                                  discharges_stepdown = NA,
                                  prevalence_general = NA,
                                  prevalence_critical = NA),
                    simplify = FALSE)

# convert individual-level to aggregate level
for (age_i in seq_along(age_lower)) {
  sim_age <- subset(sim, age >= age_lower[age_i] & age < age_upper[age_i])
  for (i in seq_along(t_vec)) {
    j <- t_vec[i]
    sitrep[[age_i]]$obs_admissions[i] <- sum(sim_age$date_test == j)
    sitrep[[age_i]]$deaths_general[i] <- sum(sim_age$date_test <= j & !sim_age$icu & sim_age$final_outcome == "death" & sim_age$date_final_outcome == j)
    sitrep[[age_i]]$discharges_general[i] <- sum(sim_age$date_test <= j & !sim_age$icu & sim_age$final_outcome == "discharge" & sim_age$date_final_outcome == j)
    sitrep[[age_i]]$stepup[i] <- sum(sim_age$date_test <= j & sim_age$date_icu == j, na.rm = TRUE)
    sitrep[[age_i]]$deaths_critical[i] <- sum(sim_age$date_test <= j & sim_age$icu & sim_age$final_outcome == "death" & sim_age$date_final_outcome == j)
    sitrep[[age_i]]$stepdown[i] <- sum(sim_age$date_test <= j & sim_age$icu & sim_age$final_outcome == "discharge" & sim_age$date_stepdown == j)
    sitrep[[age_i]]$discharges_stepdown[i] <- sum(sim_age$date_test <= j & sim_age$icu & sim_age$final_outcome == "discharge" & sim_age$date_final_outcome == j)
    sitrep[[age_i]]$prevalence_general[i] <- sum(sim_age$date_test <= j & !sim_age$icu & sim_age$date_final_outcome >= j) +
                                             sum(sim_age$date_test <= j & sim_age$icu & sim_age$final_outcome == "discharge" & sim_age$date_stepdown <= j & sim_age$date_final_outcome >= j)
    sitrep[[age_i]]$prevalence_critical[i] <- sum(sim_age$date_test <= j & sim_age$icu & sim_age$final_outcome == "death" & sim_age$date_final_outcome > j) +
                                              sum(sim_age$date_test <= j & sim_age$icu & sim_age$final_outcome == "discharge" & sim_age$date_stepdown > j)
  }
  
  # clean up
  sitrep[[age_i]] <- subset(sitrep[[age_i]], t > 0)
  sitrep[[age_i]]$date_numeric <- sitrep[[age_i]]$t
  sitrep[[age_i]]$daily_influx <- sitrep[[age_i]]$obs_admissions
  sitrep[[age_i]]$total_hdu_icu <- sitrep[[age_i]]$prevalence_critical
  sitrep[[age_i]]$total_general <- sitrep[[age_i]]$prevalence_general
  sitrep[[age_i]]$deaths <- sitrep[[age_i]]$deaths_general + sitrep[[age_i]]$deaths_critical
  sitrep[[age_i]]$new_discharges <- sitrep[[age_i]]$discharges_general + sitrep[[age_i]]$discharges_stepdown
  
  sitrep[[age_i]] <- subset(sitrep[[age_i]], select = c(date_numeric, daily_influx, total_hdu_icu,
                                                        total_general, deaths, new_discharges))
  
}

# make list over regions
sitrep_list <- list(sitrep)

# ------------------------------------------------------------------------------------------------

# get age weights
age_group <- cut(0:110, breaks = c(age_lower, Inf), right = FALSE)
age_tab <- rep(1 / 111, 111)
age_weights <- split(age_tab, f = age_group)
age_values <- split(seq_along(age_tab) - 1, f = age_group)

# get number of rows in each sitrep
n_date_sitrep <- length(sitrep_list[[1]][[1]][[1]])

# specify x-positions (numeric dates) in admissions spline
n_node <- 7
node_x <- round(seq(-15, n_date_sitrep, length.out = n_node))








plot(t_vec, true_admissions_curve)
abline(v = node_x)
points(node_x, true_admissions_curve[m], col = 2)





# get longest interval that could possibly be required for lookup table
lookup_max <-  n_date_sitrep - min(node_x) + 1

# define age spline nodes
max_indlevel_age <- 110
age_seq <- c(seq(0, 100, 20), max_indlevel_age)

p_AI_nodex <- age_seq
p_AI_noden <- length(p_AI_nodex)

p_AD_nodex <- age_seq
p_AD_noden <- length(p_AD_nodex)

p_ID_nodex <- age_seq
p_ID_noden <- length(p_ID_nodex)

m_AC_nodex <- age_seq
m_AC_noden <- length(m_AC_nodex)

names(sitrep_list[[1]][[1]])

# create final data list
data_list <- list(sitrep = sitrep_list,
                  indlevel = indlevel_list,
                  node_x = node_x,
                  lookup_max = lookup_max,
                  n_region = length(sitrep_list),
                  n_age_sitrep = length(sitrep_list[[1]]),
                  n_date_sitrep = n_date_sitrep,
                  max_indlevel_age = max_indlevel_age,
                  age_weights = age_weights,
                  age_values = age_values,
                  p_AI_nodex = p_AI_nodex,
                  p_AD_nodex = p_AD_nodex,
                  p_ID_nodex = p_ID_nodex,
                  m_AC_nodex = m_AC_nodex)

# ------------------------------------------------------------------
# MCMC parameters

# create parameters dataframe
n_region <- length(sitrep_list)
eg <- expand.grid(1:n_node, 1:n_region)
df_params <- rbind(data.frame(name = sprintf("region%s_node%s", eg[,2], eg[,1]), min = 0, max = 10, init = 1, region = eg[,2]),
                   data.frame(name = sprintf("p_AI_node%s", 1:p_AI_noden), min = -5, max = 5, init = 0, region = 0),
                   data.frame(name = sprintf("p_AD_node%s", 1:p_AD_noden), min = -5, max = 5, init = 0, region = 0),
                   data.frame(name = sprintf("p_ID_node%s", 1:p_ID_noden), min = -5, max = 5, init = 0, region = 0),
                   data.frame(name = "m_AI", min = 0, max = 20, init = 1, region = 0),
                   data.frame(name = "m_AD", min = 0, max = 20, init = 1, region = 0),
                   data.frame(name = sprintf("m_AC_node%s", 1:m_AC_noden), min = -5, max = 5, init = 0, region = 0),
                   data.frame(name = "m_ID", min = 0, max = 20, init = 1, region = 0),
                   data.frame(name = "m_IS", min = 0, max = 20, init = 1, region = 0),
                   data.frame(name = "m_SC", min = 0, max = 20, init = 1, region = 0),
                   data.frame(name = "s_AI", min = 0, max = 1, init = 0.9, region = 0),
                   data.frame(name = "s_AD", min = 0, max = 1, init = 0.7, region = 0),
                   data.frame(name = "s_AC", min = 0, max = 1, init = 0.7, region = 0),
                   data.frame(name = "s_ID", min = 0, max = 1, init = 0.9, region = 0),
                   data.frame(name = "s_IS", min = 0, max = 1, init = 0.9, region = 0),
                   data.frame(name = "s_SC", min = 0, max = 1, init = 0.9, region = 0),
                   data.frame(name = "c_AD", min = 0, max = 1, init = 0.5, region = 0),
                   data.frame(name = "c_AC", min = 0, max = 1, init = 0.5, region = 0),
                   data.frame(name = "c_ID", min = 0, max = 1, init = 0.5, region = 0),
                   data.frame(name = "c_IS", min = 0, max = 1, init = 0.5, region = 0),
                   data.frame(name = "m_AL", min = 0, max = 10, init = 1, region = 0),
                   data.frame(name = sprintf("scale_p_AI%s", 1:n_region), min = 0.1, max = 10, init = 1, region = 0),
                   data.frame(name = sprintf("scale_p_AD%s", 1:n_region), min = 0.1, max = 10, init = 1, region = 0),
                   data.frame(name = sprintf("scale_p_ID%s", 1:n_region), min = 1, max = 1, init = 1, region = 0)
)

# set initial values
m <- match(names(true_params), df_params$name)
df_params$init[m[!is.na(m)]] <- true_params[which(!is.na(m))]

# set fixed parameters
fixed_params <- setdiff(names(true_params), c("m_AD")[-1])
m <- match(fixed_params, df_params$name)
df_params$min[m[!is.na(m)]] <- df_params$max[m[!is.na(m)]] <- df_params$init[m[!is.na(m)]]

# fix regional node y values
m <- match(node_x, t_vec)
region_fitted <- log(true_admissions_curve[m])
for (i in seq_along(region_fitted)) {
  df_params$min[i] <- df_params$max[i] <- df_params$init[i] <- region_fitted[i]
}

# append update rules to data list
data_list$update_region <- df_params$region

# ------------------------------------------------------------------
# Run MCMC

# define MCMC parameters
burnin <- 1e1
samples <- 1e1
chains <- 1
beta_vec <- 1
pb_markdown <- FALSE

beta_vec <- rev(c(seq(1e-5, 9.5e-5, 5e-6),
                  seq(1e-4, 9e-4, 1e-4),
                  1e-3, 3e-3, 6e-3, 1e-2, 3e-2, 6e-2, 0.1, 0.5, 1))
beta_vec <- c(1)
length(beta_vec)

# run MCMC
set.seed(1)
t0 <- Sys.time()
mcmc <- run_mcmc(data_list = data_list,
                 df_params = df_params,
                 burnin = burnin,
                 samples = samples,
                 chains = chains,
                 beta_vec = beta_vec,
                 pb_markdown = pb_markdown)

print(Sys.time() - t0)

#as.data.frame(mcmc$diagnostics$mc_accept)
plot_credible(mcmc)
#plot_par(mcmc, show = "m_AD", phase = "both")
#plot_par(mcmc, show = "scale_p_AD1", phase = "both")
#subset(mcmc$diagnostics$mc_accept, stage == "sampling")



# subset to sampling iterations and correct rung
mcmc_samples_raw <- subset(mcmc$output, stage == "sampling" & rung == "rung1")

# subset to desired parameter columns only
mcmc_samples <- subset(mcmc_samples_raw, select = -c(chain, rung, iteration, stage, logprior, loglikelihood))

# create parameter descriptions dataframe
param_desc <- setNames(rep(NA, ncol(mcmc_samples)), names(mcmc_samples))
param_desc[sprintf("region%s_node%s", eg[,2], eg[,1])] <- sprintf("spline t%s", 1:n_node)
param_desc[sprintf("p_AI_node%s", 1:p_AI_noden)] <- "probability\nadmission\nto ICU\nspline node"
param_desc[sprintf("p_AD_node%s", 1:p_AD_noden)] <- "probability\ngeneral ward\nto death\nspline node"
param_desc[sprintf("p_ID_node%s", 1:p_ID_noden)] <- "probability\nICU\nto death\nspline node"
param_desc["m_AI"] <- "mean\nadmission\nto ICU"
param_desc["m_AD"] <- "mean\nadmission\nto death"
param_desc[sprintf("m_AC_node%s", 1:m_AC_noden)] <- "mean\nadmission\nto discharge\nspline node"
param_desc["m_ID"] <- "mean\nICU\nto death"
param_desc["m_IS"] <- "mean\nICU\nto stepdown"
param_desc["m_SC"] <- "mean\nstepdown\nto discharge"
param_desc["s_AI"] <- "CV\nadmission\nto ICU"
param_desc["s_AD"] <- "CV\nadmission\nto death"
param_desc["s_AC"] <- "CV\nadmission\nto discharge"
param_desc["s_ID"] <- "CV\nICU\nto death"
param_desc["s_IS"] <- "CV\nICU\nto stepdown"
param_desc["s_SC"] <- "CV\nstepdown\nto discharge"
param_desc["c_AD"] <- "censoring admission to death"
param_desc["c_AC"] <- "censoring admission to discharge"
param_desc["c_ID"] <- "censoring ICU to death"
param_desc["c_IS"] <- "censoring ICU to stepdown"
param_desc["m_AL"] <- "mean\nadmission\nto lab result"
param_desc[sprintf("scale_p_AI%s", 1:n_region)] <- "scale\nprobability\nadmission\nto ICU"
param_desc[sprintf("scale_p_AD%s", 1:n_region)] <- "scale\nprobability\nadmission\nto death"
param_desc[sprintf("scale_p_ID%s", 1:n_region)] <- "scale\nprobability\nICU\nto death"


df_param_desc <- data.frame(param = names(param_desc),
                            description = param_desc)



# ------------------------------------------------------------------
# Save diagnostics to file

# save trace plots
save_diag <- FALSE
if (save_diag) {
  trace_both <- plot_par(mcmc, phase = "both", display = FALSE)
  trace_sampling <- NULL#plot_par(mcmc, phase = "sampling", display = FALSE)
} else {
  trace_both <- NULL
  trace_sampling <- NULL
}

# save credible interval plot of all parameters
credible_all <- plot_credible(mcmc)

# save diagnostics to file
diags <- list(df_params = df_params,
              df_param_desc = df_param_desc,
              region_list = region_list,
              trace_both = trace_both,
              trace_sampling = trace_sampling,
              credible_all = credible_all,
              mcmc_diagnostics = mcmc$diagnostics)
saveRDS(diags, "ignore/output/diagnostics_mcmc4.rds")


# ------------------------------------------------------------------
# Save posterior summaries to file

# get 95% credible intervals
df_summary <- as.data.frame(t(apply(mcmc_samples, 2, quantile_95)))

# get posterior means and variances
df_summary$mean <- colMeans(mcmc_samples)
df_summary$var <- apply(mcmc_samples, 2, var)

# get maximum a-posteriori
w <- which.max(mcmc_samples_raw$logprior + mcmc_samples_raw$loglikelihood)
df_summary$MAP <- as.vector(unlist(mcmc_samples[w,]))

# append parameter names and descriptions
df_summary <- cbind(param = df_param_desc$param,
                    description = stringr::str_replace_all(df_param_desc$description, "\n", " "),
                    df_summary)

# save results to file
write.csv(df_summary, "ignore/output/summary_mcmc4.csv", row.names = FALSE)
#df_summary <- read.csv("ignore/output/summary_mcmc4.csv")

# get best posterior parameter estimates
params <- df_summary$mean
#params <- df_summary$MAP
names(params) <- df_summary$param


# ------------------------------------------------------------------
# Save spline posterior summaries to file

# function for getting transformed spline from raw nodes
get_spline <- function(mcmc_samples, nodex, age_vec, scale = 1) {
  mapply(function(i) {
    ret <- cubic_spline(nodex, mcmc_samples[i,], age_vec)
    scale / (1 + exp(-ret))
  }, seq_len(nrow(mcmc_samples)))
}

# function for getting 95% CrI from mcmc draws
get_spline_quantiles <- function(spline_draws) {
  ret <- as.data.frame(t(apply(spline_draws, 1, quantile_95)))
  ret$age = seq_len(nrow(ret)) - 1
  ret
}

# get quantiles for transition probabilities
mcmc_samples_p_AI <- as.matrix(mcmc_samples[, sprintf("p_AI_node%s", 1:p_AI_noden), drop = FALSE])
p_AI_spline <- get_spline(mcmc_samples_p_AI, p_AI_nodex, 0:max_indlevel_age)
p_AI_quantile <- get_spline_quantiles(p_AI_spline)

mcmc_samples_p_AD <- as.matrix(mcmc_samples[, sprintf("p_AD_node%s", 1:p_AD_noden), drop = FALSE])
p_AD_spline <- get_spline(mcmc_samples_p_AD, p_AD_nodex, 0:max_indlevel_age)
p_AD_quantile <- get_spline_quantiles(p_AD_spline)

mcmc_samples_p_ID <- as.matrix(mcmc_samples[, sprintf("p_ID_node%s", 1:p_ID_noden), drop = FALSE])
p_ID_spline <- get_spline(mcmc_samples_p_ID, p_ID_nodex, 0:max_indlevel_age)
p_ID_quantile <- get_spline_quantiles(p_ID_spline)

p_OD_spline <- (1 - p_AI_spline)*p_AD_spline + p_AI_spline*p_ID_spline
p_OD_quantile <- get_spline_quantiles(p_OD_spline)

mcmc_samples_m_AC <- as.matrix(mcmc_samples[, sprintf("m_AC_node%s", 1:m_AC_noden), drop = FALSE])
m_AC_spline <- get_spline(mcmc_samples_m_AC, m_AC_nodex, 0:max_indlevel_age, scale = 20)
m_AC_quantile <- get_spline_quantiles(m_AC_spline)


# save to file
age_spline_list <- list(p_AI = p_AI_quantile,
                        p_AD = p_AD_quantile,
                        p_ID = p_ID_quantile,
                        p_OD = p_OD_quantile,
                        m_AC = m_AC_quantile)
saveRDS(age_spline_list, file = "ignore/output/age_spline.rds")



# ------------------------------------------------------------------
# Save data quantiles to file

# admission to ICU
df_quantiles_p_AI <- as.data.frame(t(mapply(function(i) {
  w <- which(indlevel$age == i-1)
  n1 <- sum(indlevel$icu[w] == TRUE, na.rm = TRUE)
  n2 <- sum(indlevel$icu[w] == FALSE, na.rm = TRUE)
  if (n1 == 0 & n2 == 0) {
    return(c(i-1,rep(NA, 3)))
  }
  tmp <- epitools::binom.exact(n1, n1 + n2)
  c(i-1, tmp$lower, tmp$upper, n1 / (n1 + n2))
}, 1:(max_indlevel_age + 1))))
names(df_quantiles_p_AI) <- c("age", "lower", "upper", "mean")

# general ward to death
df_quantiles_p_AD <- as.data.frame(t(mapply(function(i) {
  w <- which(indlevel$age == i-1 & indlevel$icu == FALSE)
  n1 <- sum(indlevel$final_outcome[w] == "death", na.rm = TRUE)
  n2 <- sum(indlevel$final_outcome[w] == "discharge", na.rm = TRUE)
  if (n1 == 0 & n2 == 0) {
    return(c(i-1,rep(NA, 3)))
  }
  tmp <- epitools::binom.exact(n1, n1 + n2)
  c(i-1, tmp$lower, tmp$upper, n1 / (n1 + n2))
}, 1:(max_indlevel_age + 1))))
names(df_quantiles_p_AD) <- c("age", "lower", "upper", "mean")

# ICU to death
df_quantiles_p_ID <- as.data.frame(t(mapply(function(i) {
  w <- which(indlevel$age == i-1 & indlevel$icu == TRUE)
  n1 <- sum(indlevel$final_outcome[w] == "death", na.rm = TRUE)
  n2 <- sum(indlevel$final_outcome[w] == "discharge", na.rm = TRUE)
  if (n1 == 0 & n2 == 0) {
    return(c(i-1,rep(NA, 3)))
  }
  tmp <- epitools::binom.exact(n1, n1 + n2)
  c(i-1, tmp$lower, tmp$upper, n1 / (n1 + n2))
}, 1:(max_indlevel_age + 1))))
names(df_quantiles_p_ID) <- c("age", "lower", "upper", "mean")

# overall death
df_quantiles_p_OD <- as.data.frame(t(mapply(function(i) {
  w <- which(indlevel$age == i-1)
  n1 <- sum(indlevel$final_outcome[w] == "death", na.rm = TRUE)
  n2 <- sum(indlevel$final_outcome[w] == "discharge", na.rm = TRUE)
  if (n1 == 0 & n2 == 0) {
    return(c(i-1,rep(NA, 3)))
  }
  tmp <- epitools::binom.exact(n1, n1 + n2)
  c(i-1, tmp$lower, tmp$upper, n1 / (n1 + n2))
}, 1:(max_indlevel_age + 1))))
names(df_quantiles_p_OD) <- c("age", "lower", "upper", "mean")

# time general admission to death
df_quantiles_m_AC <- as.data.frame(t(mapply(function(i) {
  w <- which(indlevel$age == i-1 & indlevel$icu == FALSE & indlevel$final_outcome == "discharge")
  delta <- indlevel$date_final_outcome[w] - indlevel$date_admission[w]
  delta <- delta[!is.na(delta)]
  if (length(delta) == 0) {
    return(c(i-1,rep(NA, 3)))
  }
  m <- mean(delta)
  s <- sd(delta) / sqrt(length(delta))
  c(i-1, m - 1.96*s, m + 1.96*s, m)
}, 1:(max_indlevel_age + 1))))
names(df_quantiles_m_AC) <- c("age", "lower", "upper", "mean")


# save list to file
data_quantiles_list <- list(p_AI = df_quantiles_p_AI,
                            p_AD = df_quantiles_p_AD,
                            p_ID = df_quantiles_p_ID,
                            p_OD = df_quantiles_p_OD,
                            m_AC = df_quantiles_m_AC)
saveRDS(data_quantiles_list, file = "ignore/output/data_quantiles.rds")


# ------------------------------------------------------------------
# Save ccdf to file

# get ccdf
s_name <- c("s_AI", "s_AD", "s_AC", "s_ID", "s_IS", "s_SC")
s_desc <- c("admission to ICU",
            "admission to death",
            "admission to discharge",
            "ICU to death",
            "ICU to stepdown",
            "stepdown to discharge")
m <- rep(1,6)
s <- df_summary$mean[match(s_name, df_summary$param)]
df_ccdf <- cbind(description = s_desc,
                 get_ccdf(m, s))

# write ccdf to file
write.csv(df_ccdf, "ignore/output/ccdf_mcmc4.csv", row.names = FALSE)



# ------------------------------------------------------------------
# Save regional proportions to file

# get proportion in each 5-year age band
weight_mat <- matrix(unlist(age_weights)[1:110], nrow = 5)
df_prop <- data.frame(age = sprintf("%s-%s", seq(0,105,5), seq(4,109,5)),
                      prop = colSums(weight_mat))

# sum splines by 5-year age bands
p_AI_mat <- matrix(p_AI_quantile$Q50[1:110], nrow = 5)
df_prop$p_AI <- colSums(weight_mat * p_AI_mat) / colSums(weight_mat)

p_AD_mat <- matrix(p_AD_quantile$Q50[1:110], nrow = 5)
df_prop$p_AD <- colSums(weight_mat * p_AD_mat) / colSums(weight_mat)

p_ID_mat <- matrix(p_ID_quantile$Q50[1:110], nrow = 5)
df_prop$p_ID <- colSums(weight_mat * p_ID_mat) / colSums(weight_mat)

m_AC_mat <- matrix(m_AC_quantile$Q50[1:110], nrow = 5)
df_prop$m_AC <- colSums(weight_mat * m_AC_mat) / colSums(weight_mat)

# write to file
write.csv(df_prop, "ignore/output/age_prop_mcmc4.csv", row.names = FALSE)


# ------------------------------------------------------------------
# Save individual-level fits to file

# simulate individual-level data from model
sim_data <- sim_indlevel(params = params,
                         p_AI = p_AI_quantile$Q50,
                         p_AD = p_AD_quantile$Q50,
                         p_ID = p_ID_quantile$Q50,
                         m_AC = m_AC_quantile$Q50,
                         date_admission = indlevel_list$date_admission,
                         date_censor = indlevel_list$date_censor,
                         age = indlevel_list$age,
                         n_samp = 1e6)


# bin delay distributions
bin_real <- bin_indlevel(indlevel, max_days = 30)
bin_sim <- bin_indlevel(sim_data, max_days = 30)

# combine data and model fit
bin_real$type <- "data"
bin_sim$type <- "model"
bin_combined <- rbind(bin_real, bin_sim)

# save results to file
write.csv(bin_combined, "ignore/output/indlevel_fit_mcmc4.csv", row.names = FALSE)



# ------------------------------------------------------------------
# Save SitRep fits to file

# get model fit
df_fit <- df_params
df_fit$init <- params[match(df_fit$name, names(params))]
df_model_fit <- run_mcmc(data_list = data_list,
                         df_params = df_fit,
                         return_fit = TRUE)
df_model_fit$by_age$x <- df_model_fit$by_age$x + node_x[1]

# get sitrep into same format
df_sitrep <- nested_to_long(sitrep_list)
names(df_sitrep) <- c("x", "value", "metric", "age", "region")
df_sitrep$metric <- names(sitrep_list[[1]][[1]])[df_sitrep$metric]
df_sitrep <- subset(df_sitrep, metric != "date_numeric")
metric_rename <- c(daily_influx = "admission_incidence",
                   deaths = "deaths_incidence",
                   new_discharges = "discharges_incidence",
                   total_general = "general_prevalence",
                   total_hdu_icu = "critical_prevalence")
df_sitrep$metric <- metric_rename[df_sitrep$metric]
df_sitrep$value[df_sitrep$value == -1] <- NA

# combine model fit with sitrep
df_sitrep_fit <- rbind(cbind(df_model_fit$by_age[,names(df_sitrep)], type = "model"),
                       cbind(df_sitrep, type = "data"))

# get combined fit
df_sitrep_fit_combined <- subset(df_sitrep_fit, region == 1)
for (i in seq_len(n_region - 1)) {
  tmp <- subset(df_sitrep_fit, region == i + 1)
  df_sitrep_fit_combined$value <- df_sitrep_fit_combined$value + tmp$value
}
df_sitrep_fit_combined$region <- 0

# append to fit dataframe
df_sitrep_fit <- rbind(df_sitrep_fit, df_sitrep_fit_combined)

# save combined sitrep + fit to file
write.csv(df_sitrep_fit, "ignore/output/sitrep_fit_mcmc4.csv", row.names = FALSE)



#rmarkdown::render("R_ignore/diagnostics_mcmc4.Rmd")
#rmarkdown::render("R_ignore/summary_mcmc4.Rmd")
