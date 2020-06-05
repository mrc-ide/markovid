
# deploy.R
#
# Author: Bob Verity
# Date: 2019-05-02
#
# Purpose:
# Test package functions.

# RStudio shortcuts:
#    cmd+shift+L     : load package from local version
#    cmd+shift+D     : document (NB, option must be activated in Build Tools)
#    cmd+shift+E     : check
#    cmd+shift+T     : test

# Useful commands:
# devtools::install()  # install package
# pkgdown::build_site() # build all pages of pkgdown website
# pkgdown::build_article('index')  # build single vignette
# check('.', args = '--no-examples')  # run checks without examples
# covr::report()    # interactive coverage report
# devtools::build_vignettes()
# ------------------------------------------------------------------

#devtools::install_github("mrc-ide/drjacoby", ref = "develop")
#library(drjacoby)

library(ggplot2)

set.seed(1)


# ------------------------------------------------------------------
# Prepare individual-level data

# load raw data
#cocin_raw <- read.csv("C:/Users/rverity/Desktop/ncov/ncov-outputs/src/mcmc_uk_progression/raw_data/cocin_cleaned.csv",
#                      stringsAsFactors = FALSE)
cocin_raw <- read.csv("C:/Users/rverity/Desktop/cocin_cleaned_20200421_region.csv",
                      stringsAsFactors = FALSE)

# prepare data
cocin <- prepare_cocin_v2(cocin_raw)
indlevel <- prepare_indlevel(cocin)

# create list of regions
region_names <- c("east of england", "london", "midlands",
                  "north east and yorkshire", "north west",
                  "south east", "south west")
region_list <- as.list(region_names)
names(region_list) <- region_list

region_list <- list("all England" = region_names)

# subset regions
indlevel$region <- as.character(indlevel$region)
n_region <- length(region_list)
indlevel <- subset(indlevel, region %in% unlist(region_list))

# must have an admission date within specified range. Should be chosen
# pragmatically on the basis of admissions curve with the aim of removing backfill issues
#tab1 <- table(indlevel$date_admission)
#plot(tab1); abline(v = which(names(tab1) == "2020-04-18"), col = 2, lwd = 2)
indlevel <- subset(indlevel, date_admission >= as.Date("2020-03-09") &
                     date_admission <= as.Date("2020-04-18"))

# must have an age
indlevel <- subset(indlevel, !is.na(age))

# make age integer
indlevel$age <- floor(indlevel$age)

# ------------------------------------------------------------------
# Define stratification

# define sitrep age grouping
sitrep_lower = c(0, 6, 18, 65, 85)
sitrep_upper = c(sitrep_lower[-1] - 1, Inf)
sitrep_name <- paste(sitrep_lower, sitrep_upper, sep = "-")
sitrep_name[length(sitrep_name)] <- sprintf("%s+", sitrep_lower[length(sitrep_name)])
n_age_sitrep <- length(sitrep_lower)

df_age_sitrep <- data.frame(sitrep_group = seq_along(sitrep_lower),
                            sitrep_lower = sitrep_lower,
                            sitrep_upper = sitrep_upper,
                            sitrep_name = sitrep_name,
                            stringsAsFactors = FALSE)

# define individual-level age grouping
indlevel_lower <- c(0, 6, 18, seq(45, 85, 5))
#indlevel_lower <- seq(0, 85, 5)
#indlevel_lower <- c(0, 6, 18, 65, 70, 85)
#indlevel_lower <- c(0, 6, 18, 65, 85)
indlevel_upper <- c(indlevel_lower[-1] - 1, Inf)
indlevel_name <- paste(indlevel_lower, indlevel_upper, sep = "-")
indlevel_name[length(indlevel_name)] <- sprintf("%s+", indlevel_lower[length(indlevel_name)])
n_age_indlevel <- length(indlevel_lower)

df_age_indlevel <- data.frame(indlevel_group = seq_along(indlevel_lower),
                              indlevel_lower = indlevel_lower,
                              indlevel_upper = indlevel_upper,
                              indlevel_name = indlevel_name,
                              stringsAsFactors = FALSE)

# add age groups to individual-level line list
indlevel$age_group <- as.numeric(cut(indlevel$age,
                                     breaks = c(df_age_indlevel$indlevel_lower, Inf),
                                     right = FALSE))

# get expected proportions from individual-level line list
tab1 <- table(indlevel$age_group)
df_age_indlevel$prop <- as.vector(tab1 / sum(tab1))

# convert to proportions relative to oldest age group
df_age_indlevel$rel_prop <- df_age_indlevel$prop / df_age_indlevel$prop[n_age_indlevel]

# map sitrep groups to indlevel groups
age_cut <- as.numeric(cut(df_age_indlevel$indlevel_lower,
                          breaks = c(df_age_sitrep$sitrep_lower, Inf),
                          right = FALSE))
df_age_indlevel$sitrep_group = age_cut

# map indlevel groups to sitrep groups
df_age_sitrep$indlevel_group <- split(1:n_age_indlevel, f = age_cut)


# ------------------------------------------------------------------
# Prepare SitRep data

# load raw SitRep and deaths data
sitrep_raw <- readRDS("C:/Users/rverity/Desktop/ncov/ncov-outputs/src/mcmc_uk_progression/raw_data/combined_covid_sitrep_by_region.rds")
deaths_raw <- readRDS("C:/Users/rverity/Desktop/deaths_full_linelist.rds")

# fix region format
deaths_raw$region <- tolower(deaths_raw$nhser_name)

# filter sitrep on dates
sitrep_raw <- subset(sitrep_raw, sitrep_raw$date <= as.Date("2020-04-26"))

# prepare deaths
deaths <- prepare_deaths(deaths_raw)

# specify max date at which death data should be included. Should be chosen
# pragmatically on the basis of deaths curve with the aim of removing backfill issues
#tab1 <- table(deaths$date_admission)
#plot(tab1); abline(v = which(names(tab1) == "2020-05-01"), col = 2, lwd = 2)
deaths_max_date <- as.Date("2020-05-01")

# prepare each region separately
sitrep_list <- list()
for (i in seq_along(region_list)) {
  
  # sum sitrep over regions in this element of region_list
  sitrep_i <- data.frame(date = sitrep_raw$date,
                         metric_name = sitrep_raw$metric_name,
                         value = rowSums(sitrep_raw[, region_list[[i]], drop = FALSE]),
                         stringsAsFactors = FALSE)
  
  # subset deaths to regions in this element of region_list
  deaths_i <- subset(deaths, region %in% region_list[[i]])
  
  # prepare sitrep and merge deaths
  sitrep_i <- prepare_sitrep_age(sitrep_i, deaths_i)
  
  # apply death date cutoff
  sitrep_i$deaths[sitrep_i$date > deaths_max_date] <- NA
  
  # create new field for total daily influx
  sitrep_i$daily_influx <- sitrep_i$new_admissions + sitrep_i$new_inpatients_diagnosed
  
  # split by age group
  sitrep_list[[i]] <- split(sitrep_i, f = sitrep_i$age_group)
  
  # subset fields
  for (j in seq_along(sitrep_list[[i]])) {
    sitrep_list[[i]][[j]] <- as.list(subset(sitrep_list[[i]][[j]],
                                            select = c(date_numeric,
                                                       daily_influx,
                                                       total_hdu_icu,
                                                       total_general,
                                                       deaths,
                                                       new_discharges)))
  }
}


# ------------------------------------------------------------------
# Finalise individual-level data

# make dates numeric
for (i in c("date_admission", "date_labtest", "date_icu", "date_leave_icu",
            "date_stepdown", "date_final_outcome", "date_censor")) {
  indlevel[[i]] <- as.numeric(indlevel[[i]] - as.Date("2020-03-09"))
}

# get numeric outcome
indlevel$final_outcome_numeric <- as.numeric(indlevel$final_outcome)

# subset to list of values to pass into MCMC
indlevel_list <- as.list(subset(indlevel, select = c(age_group,
                                                     age,
                                                     icu,
                                                     stepdown,
                                                     date_admission,
                                                     date_labtest,
                                                     date_icu,
                                                     date_stepdown,
                                                     date_final_outcome,
                                                     final_outcome_numeric,
                                                     date_censor)))


# ------------------------------------------------------------------
# Final data prep

# get number of rows in each sitrep
n_date_sitrep <- length(sitrep_list[[1]][[1]][[1]])

# specify x-positions (numeric dates) in admissions spline
n_node <- 6
node_x <- round(seq(-15, n_date_sitrep, length.out = n_node))

# ensure no duplicated spline points
if (any(duplicated(node_x))) {
  stop("cannot have duplicated spline points")
}

# get longest interval that could possibly be required for lookup table
max_indlevel <- max(indlevel_list$date_final_outcome - indlevel_list$date_admission, na.rm = TRUE)
max_sitrep <- n_date_sitrep - min(node_x)
lookup_max <- max(max_indlevel, max_sitrep) + 1

# replace NAs in lists with -1 (understood by C++)
for (i in seq_along(sitrep_list)) {
  for (j in seq_along(sitrep_list[[i]])) {
    for (k in seq_along(sitrep_list[[i]][[j]])) {
      sitrep_list[[i]][[j]][[k]][is.na(sitrep_list[[i]][[j]][[k]])] <- -1
    }
  }
}
for (i in seq_along(indlevel_list)) {
  indlevel_list[[i]][is.na(indlevel_list[[i]])] <- -1
}

# define age spline nodes
#print(max(indlevel$age))
max_indlevel_age <- 110

p_AI_nodex <- c(0, 25, 50, 75, max_indlevel_age)
p_AI_noden <- length(p_AI_nodex)

p_AD_nodex <- c(0, 25, 50, 75, max_indlevel_age)
p_AD_noden <- length(p_AD_nodex)

# create final data list
data_list <- list(sitrep = sitrep_list,
                  indlevel = indlevel_list,
                  n_region = n_region,
                  node_x = node_x,
                  lookup_max = lookup_max,
                  n_date_sitrep = n_date_sitrep,
                  n_age_sitrep = n_age_sitrep,
                  n_age_indlevel = n_age_indlevel,
                  rel_prop = df_age_indlevel$rel_prop,
                  map_age_indlevel = df_age_indlevel$sitrep_group,
                  map_age_sitrep = df_age_sitrep$indlevel_group,
                  max_indlevel_age = max_indlevel_age,
                  p_AI_nodex = p_AI_nodex,
                  p_AD_nodex = p_AD_nodex)


# ------------------------------------------------------------------
# MCMC parameters

# define MCMC parameters
burnin <- 3e2
samples <- 3e2
chains <- 1
beta_vec <- 1
pb_markdown <- FALSE

# create parameters dataframe
create_param_df <- function(name, min, max, init, density) {
  data.frame(name = paste0(name, 1:n_age_indlevel), min = min, max = max, init = init,
             density = density,
             region = 0,
             indlevel_age = 1:n_age_indlevel,
             sitrep_age = df_age_indlevel$sitrep_group)
}
eg <- expand.grid(1:n_node, 1:n_region)
df_params <- rbind(data.frame(name = sprintf("region%s_node%s", eg[,2], eg[,1]), min = 0, max = 10, init = 1, density = -1, region = eg[,2], indlevel_age = -1, sitrep_age = 0),
                   data.frame(name = paste0("scale_rel_prop", 1:n_age_sitrep), min = 0, max = 5, init = 1, density = -1, region = 0, indlevel_age = -1, sitrep_age = 1:n_age_sitrep),
                   data.frame(name = sprintf("p_AI_node%s", 1:p_AI_noden), min = -5, max = 5, init = 0, density = -1, region = -1, indlevel_age = -1, sitrep_age = -1),
                   data.frame(name = sprintf("p_AD_node%s", 1:p_AD_noden), min = -5, max = 5, init = 0, density = -1, region = -1, indlevel_age = -1, sitrep_age = -1),
                   #create_param_df("p_AI", min = 0, max = 1, init = 0.5, density = -1),
                   #create_param_df("p_AD", min = 0, max = 1, init = 0.5, density = -1),
                   create_param_df("p_ID", min = 0, max = 1, init = 0.5, density = -1),
                   data.frame(name = "m_AI", min = 0, max = 20, init = 5, density = 1, region = 0, indlevel_age = 0, sitrep_age = 0),
                   data.frame(name = "m_AD", min = 0, max = 20, init = 5, density = 2, region = 0, indlevel_age = 0, sitrep_age = 0),
                   create_param_df("m_AC", min = 0, max = 20, init = 5, density = 3),
                   data.frame(name = "m_ID", min = 0, max = 20, init = 5, density = 4, region = 0, indlevel_age = 0, sitrep_age = 0),
                   data.frame(name = "m_IS", min = 0, max = 20, init = 5, density = 5, region = 0, indlevel_age = 0, sitrep_age = 0),
                   data.frame(name = "m_SC", min = 0, max = 20, init = 5, density = 6, region = 0, indlevel_age = 0, sitrep_age = 0),
                   data.frame(name = "s_AI", min = 0, max = 1, init = 0.9, density = 1, region = 0, indlevel_age = 0, sitrep_age = 0),
                   data.frame(name = "s_AD", min = 0, max = 1, init = 0.9, density = 2, region = 0, indlevel_age = 0, sitrep_age = 0),
                   data.frame(name = "s_AC", min = 0, max = 1, init = 0.9, density = 3, region = 0, indlevel_age = 0, sitrep_age = 0),
                   data.frame(name = "s_ID", min = 0, max = 1, init = 0.9, density = 4, region = 0, indlevel_age = 0, sitrep_age = 0),
                   data.frame(name = "s_IS", min = 0, max = 1, init = 0.9, density = 5, region = 0, indlevel_age = 0, sitrep_age = 0),
                   data.frame(name = "s_SC", min = 0, max = 1, init = 0.9, density = 6, region = 0, indlevel_age = 0, sitrep_age = 0),
                   data.frame(name = "m_AL", min = 3.5, max = 3.5, init = 3.5, density = -1, region = 0, indlevel_age = -1, sitrep_age = 0),
                   data.frame(name = sprintf("scale_p_AI%s", 1:n_region), min = 0.1, max = 10, init = 1, density = -1, region = 1:n_region, indlevel_age = -1, sitrep_age = 0),
                   data.frame(name = sprintf("scale_p_AD%s", 1:n_region), min = 0.1, max = 10, init = 1, density = -1, region = 1:n_region, indlevel_age = -1, sitrep_age = 0),
                   data.frame(name = sprintf("scale_p_ID%s", 1:n_region), min = 1, max = 1, init = 1, density = -1, region = 1:n_region, indlevel_age = -1, sitrep_age = 0)
)

# fix scaling in highest age group
w <- which(df_params$name == sprintf("scale_rel_prop%s", n_age_sitrep))
df_params$min[w] <- df_params$max[w] <- df_params$init[w] <- 1.0

# append update rules to data list
data_list$update_density <- df_params$density
data_list$update_region <- df_params$region
data_list$update_indlevel_age <- df_params$indlevel_age
data_list$update_sitrep_age <- df_params$sitrep_age


# ------------------------------------------------------------------
# Run MCMC

beta_vec <- rev(c(1e-5, 2e-5, 3e-5, 4e-5, 5e-5, 6e-5, 7e-5, 8e-5, 9e-5,
                  1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 6e-4, 7e-4, 8e-4, 9e-4,
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

as.data.frame(mcmc$diagnostics$mc_accept)

#plot_credible(mcmc)
#plot_par(mcmc, show = "p_AD_node2", phase = "both")
#subset(mcmc$diagnostics$mc_accept, stage == "sampling")

# subset to sampling iterations and correct rung
mcmc_samples_raw <- subset(mcmc$output, stage == "sampling" & rung == "rung1")

# subset to desired parameter columns only
mcmc_samples <- subset(mcmc_samples_raw, select = -c(chain, rung, iteration, stage, logprior, loglikelihood))

# create parameter descriptions dataframe
param_desc <- setNames(rep(NA, ncol(mcmc_samples)), names(mcmc_samples))
param_desc[sprintf("region%s_node%s", eg[,2], eg[,1])] <- sprintf("spline t%s", 1:n_node)
param_desc[sprintf("scale_rel_prop%s", 1:n_age_sitrep)] <- "relative\nproportion\nadmissions"
param_desc[sprintf("p_AI_node%s", 1:p_AI_noden)] <- "probability\nadmission\nto ICU\nspline node"
param_desc[sprintf("p_AD_node%s", 1:p_AD_noden)] <- "probability\ngeneral ward\nto death\nspline node"
param_desc[sprintf("p_ID%s", 1:n_age_indlevel)] <- "probability\nICU\nto death"
param_desc["m_AI"] <- "mean\nadmission\nto ICU"
param_desc["m_AD"] <- "mean\nadmission\nto death"
param_desc[sprintf("m_AC%s", 1:n_age_indlevel)] <- "mean\nadmission\nto discharge"
param_desc["m_ID"] <- "mean\nICU\nto death"
param_desc["m_IS"] <- "mean\nICU\nto stepdown"
param_desc["m_SC"] <- "mean\nstepdown\nto discharge"
param_desc["s_AI"] <- "CV\nadmission\nto ICU"
param_desc["s_AD"] <- "CV\nadmission\nto death"
param_desc["s_AC"] <- "CV\nadmission\nto discharge"
param_desc["s_ID"] <- "CV\nICU\nto death"
param_desc["s_IS"] <- "CV\nICU\nto stepdown"
param_desc["s_SC"] <- "CV\nstepdown\nto discharge"
param_desc["m_AL"] <- "mean\nadmission\nto lab result"
param_desc[sprintf("scale_p_AI%s", 1:n_region)] <- "scale\nprobability\nadmission\nto ICU"
param_desc[sprintf("scale_p_AD%s", 1:n_region)] <- "scale\nprobability\nadmission\nto death"
param_desc[sprintf("scale_p_ID%s", 1:n_region)] <- "scale\nprobability\nICU\nto death"


df_param_desc <- data.frame(param = names(param_desc),
                            description = param_desc)


# ------------------------------------------------------------------
# Save diagnostics to file

# save trace plots
save_diag <- TRUE
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
diags <- list(df_age_indlevel = df_age_indlevel,
              df_age_sitrep = df_age_sitrep,
              df_params = df_params,
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
#params <- df_summary$mean
params <- df_summary$MAP
names(params) <- df_summary$param

# ------------------------------------------------------------------
# Save spline posterior summaries to file

# function for getting 95% CrI from mcmc draws
get_spline_quantiles <- function(mcmc_samples, name, nodex, age_vec) {
  mcmc_samples <- as.matrix(mcmc_samples[,name])
  p_mat <- mapply(function(i) {
    ret <- cubic_spline(nodex, mcmc_samples[i,], age_vec)
    1 / (1 + exp(-ret))
  }, seq_len(nrow(mcmc_samples)))
  ret <- as.data.frame(t(apply(p_mat, 1, quantile_95)))
  ret$age = age_vec
  ret
}

# get quantiles for transition probabilities
p_AI_quantile <- get_spline_quantiles(mcmc_samples, sprintf("p_AI_node%s", 1:p_AI_noden), p_AD_nodex, 0:max_indlevel_age)
p_AD_quantile <- get_spline_quantiles(mcmc_samples, sprintf("p_AD_node%s", 1:p_AD_noden), p_AD_nodex, 0:max_indlevel_age)

# save to file
age_spline_list <- list(p_AI = p_AI_quantile,
                        p_AD = p_AD_quantile)
saveRDS(age_spline_list, file = "ignore/output/age_spline.rds")

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

# get final proportion in each indlevel age group from fitted scaling factors
scale_rel_prop <- params[sprintf("scale_rel_prop%s", 1:n_age_sitrep)]
final_prop <- df_age_indlevel$rel_prop * scale_rel_prop[df_age_indlevel$sitrep_group]
final_prop <- final_prop / sum(final_prop)

# write to file
write.csv(data.frame(age_group = sprintf(" %s", df_age_indlevel$indlevel_name),
                     proportion = final_prop),
          "ignore/output/age_prop_mcmc4.csv", row.names = FALSE)

# ------------------------------------------------------------------
# Save individual-level fits to file

# simulate individual-level data from model
sim_data <- sim_indlevel(params = params,
                         date_admission = indlevel_list$date_admission,
                         date_censor = indlevel_list$date_censor,
                         age_group = indlevel_list$age_group,
                         n_age = n_age_indlevel,
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
df_model_fit_indlevel <- df_model_fit$indlevel
df_model_fit <- df_model_fit$sitrep
df_model_fit$x <- df_model_fit$x + node_x[1]

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
df_sitrep_fit <- rbind(cbind(df_model_fit[,names(df_sitrep)], type = "model"),
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


#rmarkdown::render("ignore/Rmarkdown/diagnostics_mcmc4.Rmd")
#rmarkdown::render("ignore/Rmarkdown/summary_mcmc4.Rmd")
