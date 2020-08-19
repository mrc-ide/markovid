
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
library(epitools)

set.seed(4)


# ------------------------------------------------------------------
# Prepare individual-level data

# load dummy individual-level data
indlevel_list <- readRDS(system.file("extdata", "dummy_indlevel.rds", package = "markovid", mustWork = TRUE))


# ------------------------------------------------------------------
# Prepare SitRep data

# load dummy aggregate data
sitrep_list <- readRDS(system.file("extdata", "dummy_aggregate.rds", package = "markovid", mustWork = TRUE))
n_region <- length(sitrep_list)

# ------------------------------------------------------------------
# Define age weights within sitrep groups based on individual-level data

# tabulate indlevel age and split based on sitrep age bands
sitrep_lower = c(0, 6, 18, 65, 85, Inf)
sitrep_upper = c(sitrep_lower[-1] - 1, 110)
age_group <- cut(0:110, breaks = sitrep_lower, right = FALSE)
age_tab <- tabulate(indlevel$age + 1, nbins = 111)
age_weights <- split(age_tab, f = age_group)
age_values <- split(seq_along(age_tab) - 1, f = age_group)

# normalise weights to sum to 1 within sitrep groups
for (i in seq_along(age_weights)) {
  age_weights[[i]] <- age_weights[[i]] / sum(age_weights[[i]])
}


# ------------------------------------------------------------------
# Final data prep

# get number of rows in each sitrep
n_date_sitrep <- length(sitrep_list[[1]][[1]][[1]])

# get longest interval that could possibly be required for lookup table
max_indlevel1 <- max(indlevel_list$date_final_outcome - indlevel_list$date_admission, na.rm = TRUE)
max_indlevel2 <- max(indlevel_list$date_censor - indlevel_list$date_admission, na.rm = TRUE)
max_sitrep <- n_date_sitrep
lookup_max <- max(max_indlevel1, max_indlevel2, max_sitrep) + 1

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
#print(max(indlevel$age)) # max_indlevel_age should exceed this value
max_indlevel_age <- 120
age_seq <- seq(0, 120, 20)

p_AI_nodex <- age_seq
p_AI_noden <- length(p_AI_nodex)

p_AD_nodex <- age_seq
p_AD_noden <- length(p_AD_nodex)

p_ID_nodex <- age_seq
p_ID_noden <- length(p_ID_nodex)

m_AC_nodex <- age_seq
m_AC_noden <- length(m_AC_nodex)

# create final data list
data_list <- list(sitrep = sitrep_list,
                  indlevel = indlevel_list,
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
df_params <- rbind(data.frame(name = sprintf("p_AI_node%s", 1:p_AI_noden), min = -5, max = 5, init = 0),
                   data.frame(name = sprintf("p_AD_node%s", 1:p_AD_noden), min = -5, max = 5, init = 0),
                   data.frame(name = sprintf("p_ID_node%s", 1:p_ID_noden), min = -5, max = 5, init = 0),
                   data.frame(name = "m_AI", min = 0, max = 20, init = 1),
                   data.frame(name = "m_AD", min = 0, max = 20, init = 1),
                   data.frame(name = sprintf("m_AC_node%s", 1:m_AC_noden), min = -5, max = 5, init = 0),
                   data.frame(name = "m_ID", min = 0, max = 20, init = 1),
                   data.frame(name = "m_IS", min = 0, max = 20, init = 1),
                   data.frame(name = "m_SC", min = 0, max = 20, init = 1),
                   data.frame(name = "s_AI", min = 0, max = 1, init = 0.9),
                   data.frame(name = "s_AD", min = 0, max = 1, init = 0.7),
                   data.frame(name = "s_AC", min = 0, max = 1, init = 0.7),
                   data.frame(name = "s_ID", min = 0, max = 1, init = 0.9),
                   data.frame(name = "s_IS", min = 0, max = 1, init = 0.9),
                   data.frame(name = "s_SC", min = 0, max = 1, init = 0.9),
                   data.frame(name = "c_AD", min = 0, max = 1, init = 0.5),
                   data.frame(name = "c_AC", min = 0, max = 1, init = 0.5),
                   data.frame(name = "c_ID", min = 0, max = 1, init = 0.5),
                   data.frame(name = "c_IS", min = 0, max = 1, init = 0.5),
                   data.frame(name = sprintf("scale_p_AI_region%s", 1:n_region), min = 0.1, max = 10, init = 1),
                   data.frame(name = sprintf("scale_p_AD_region%s", 1:n_region), min = 0.1, max = 10, init = 1),
                   data.frame(name = sprintf("scale_p_ID_region%s", 1:n_region), min = 0.1, max = 10, init = 1)
)


# ------------------------------------------------------------------
# Run MCMC

# define MCMC parameters
burnin <- 20
samples <- 20
chains <- 1
beta_vec <- 1
pb_markdown <- FALSE

# run MCMC
set.seed(1)
t0 <- Sys.time()
mcmc <- markovid::run_mcmc(data_list = data_list,
                           df_params = df_params,
                           burnin = burnin,
                           samples = samples,
                           chains = chains,
                           beta_vec = 1,
                           pb_markdown = pb_markdown,
                           n_threads = 8)

print(Sys.time() - t0)


