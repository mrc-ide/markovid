
# standalone script, related to sim_sitrep.R. Simulates data and compares model
# fit under true parameter values (calculated directly within R) with simulated
# distributions.

#set.seed(1)

# logit transform
logit <- function(x, z = 1) {
  -log((z - x) / x)
}

# define true parameters
true_params <- c(m_AL = 0.001,
                 p_AI = 0.2,
                 p_AD = 0.15,
                 m_AD = 5,
                 s_AD = 0.8,
                 m_AI = 10,
                 s_AI = 1,
                 m_AC = 10,
                 s_AC = 0.8,
                 p_ID = 0.99,
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
t_vec <- -15:100
true_admissions_curve <- ceiling(1e4 * dnorm(t_vec, mean = 15, sd = 5))

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
    sitrep[[age_i]]$prevalence_general[i] <- sum(sim_age$date_test <= j & !sim_age$icu & sim_age$date_final_outcome > j) +
                                             sum(sim_age$date_test <= j & sim_age$icu & sim_age$date_icu >= j) +
                                             sum(sim_age$date_test <= j & sim_age$icu & sim_age$final_outcome == "discharge" & sim_age$date_stepdown <= j & sim_age$date_final_outcome >= j)
    sitrep[[age_i]]$prevalence_critical[i] <- sum(sim_age$date_test <= j & sim_age$icu & sim_age$date_icu <= j & sim_age$final_outcome == "death" & sim_age$date_final_outcome > j) +
                                              sum(sim_age$date_test <= j & sim_age$icu & sim_age$date_icu <= j & sim_age$final_outcome == "discharge" & sim_age$date_stepdown > j)
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

# ----------------------------------------------------------------

params <- true_params

age_this <- 3

get_delay_density <- function(x, m, s) {
  pgamma(x + 1, shape = 1/s^2, scale = m*s^2) - pgamma(x, shape = 1/s^2, scale = m*s^2)
}
get_delay_tail <- function(x, m, s) {
  pgamma(x + 1, shape = 1/s^2, scale = m*s^2, lower.tail = FALSE)
}

# get true admissions spline
s <- 0
for (i in seq_along(sitrep)) {
  s <- s + sitrep[[i]]$daily_influx
}
s <- sitrep[[age_this]]$daily_influx


# get age aplines
p_AI <- get_spline(matrix(params[sprintf("p_AI_node%s", 1:p_AI_noden)], nrow = 1),
                   p_AI_nodex, 0:max_indlevel_age)[,1]

p_AD <- get_spline(matrix(params[sprintf("p_AD_node%s", 1:p_AD_noden)], nrow = 1),
                   p_AD_nodex, 0:max_indlevel_age)[,1]

p_ID <- get_spline(matrix(params[sprintf("p_ID_node%s", 1:p_ID_noden)], nrow = 1),
                   p_ID_nodex, 0:max_indlevel_age)[,1]

m_AC <- get_spline(matrix(params[sprintf("m_AC_node%s", 1:m_AC_noden)], nrow = 1),
                   m_AC_nodex, 0:max_indlevel_age, scale = 20)[,1]


# misc
w <- mapply(sum, data_list$age_weights)
#w <- c(6, 12, 47, 20, 26) / 111
t_max <- length(s)
t <- seq_len(t_max) - 1
df_check <- data.frame(true_admissions = s,
                       open_general = 0,
                       open_critical = 0,
                       open_stepdown = 0,
                       admission_incidence = 0,
                       deaths_incidence = 0,
                       discharges_incidence = 0)

# define testing delays
pos_on_day <- get_delay_density(t, params["m_AL"], 1.0)
pos_by_day <- 1 - get_delay_tail(t, params["m_AL"], 1.0)

pos_on_day <- pos_on_day * 0
pos_on_day[1] <- 1

# first pass
new_stepup <- rep(0, length(s))
for (i in seq_along(s)) {
  j <- 1:(length(s) + 1 - i)
  
  # get open cases in general
  open_general <- s[i] * ( (1 - p_AI[1]) * p_AD[1] * get_delay_tail(j - 1, params["m_AD"], params["s_AD"]) +
                                           (1 - p_AI[1]) * (1 - p_AD[1]) * get_delay_tail(j - 1, m_AC[1], params["s_AC"]) +
                                           p_AI[1] * get_delay_tail(j - 1, params["m_AI"], params["s_AI"]) )
  
  # new stepup
  new_stepup[i:length(s)] <- new_stepup[i:length(s)] +
    s[i] * p_AI[1] * get_delay_density(j - 1, params["m_AI"], params["s_AI"])
  
  # update observed admissions
  df_check$admission_incidence[i:length(s)] <- df_check$admission_incidence[i:length(s)] +
    pos_on_day[j] * s[i]
  
  # deaths in general ward
  df_check$deaths_incidence[i:length(s)] <- df_check$deaths_incidence[i:length(s)] +
    w[age_this] * s[i] * (1 - p_AI[1]) * p_AD[1] * get_delay_density(j - 1, params["m_AD"], params["s_AD"])
  
  # discharges in general ward
  df_check$discharges_incidence[i:length(s)] <- df_check$discharges_incidence[i:length(s)] +
    w[age_this] * s[i] * (1 - p_AI[1]) * (1 - p_AD[1]) * get_delay_density(j - 1, m_AC[1], params["s_AC"])
  
  # store open cases in general ward
  df_check$open_general[i:length(s)] <- df_check$open_general[i:length(s)] + open_general
  
}

# second pass
new_stepdown <- rep(0, length(s))
for (i in seq_along(s)) {
  j <- 1:(length(s) + 1 - i)
  
  # open cases in critical
  open_critical <- new_stepup[i] * ( p_ID[1] * get_delay_tail(j - 1, params["m_ID"], params["s_ID"]) +
                                    (1 - p_ID[1]) * get_delay_tail(j - 1, params["m_IS"], params["s_IS"]) )
  
  # new stepdown
  new_stepdown[i:length(s)] <- new_stepdown[i:length(s)] +
    new_stepup[i] * (1 - p_ID[1]) * get_delay_density(j - 1, params["m_IS"], params["s_IS"]);
  
  # update observed admissions
  #df_check$admission_incidence[i:length(s)] <- df_check$admission_incidence[i:length(s)] +
  #  pos_on_day[j] * open_critical
  
  # deaths in critical
  df_check$deaths_incidence[i:length(s)] <- df_check$deaths_incidence[i:length(s)] +
    new_stepup[i] * p_ID[1] * get_delay_density(j - 1, params["m_ID"], params["s_ID"])
  
  # store open cases in critical ward
  df_check$open_critical[i:length(s)] <- df_check$open_critical[i:length(s)] + open_critical
}

# third pass
for (i in seq_along(s)) {
  j <- 1:(length(s) + 1 - i)
  
  # open cases in stepdown
  open_stepdown <- new_stepdown[i] * get_delay_tail(j - 1, params["m_SC"], params["s_SC"])
  
  # update observed admissions
  #df_check$admission_incidence[i:length(s)] <- df_check$admission_incidence[i:length(s)] +
  #  pos_on_day[j] * open_stepdown
  
  # discharges in stepdown
  df_check$discharges_incidence[i:length(s)] <- df_check$discharges_incidence[i:length(s)] +
    new_stepdown[i] * get_delay_density(j - 1, params["m_SC"], params["s_SC"])
  
  # store open cases in stepdown
  df_check$open_stepdown[i:length(s)] <- df_check$open_stepdown[i:length(s)] + open_stepdown
}

# final fields
df_check$general_prevalence <- df_check$open_general + df_check$open_stepdown
df_check$critical_prevalence <- df_check$open_critical

# plot admission_incidence
plot(df_check$admission_incidence)
lines(sitrep[[1]]$date_numeric, sitrep[[age_this]]$daily_influx, pch = 20)

# plot deaths
plot(df_check$deaths_incidence)
points(sitrep[[1]]$date_numeric, sitrep[[age_this]]$deaths, pch = 20)

# plot discharges
plot(df_check$discharges_incidence)
points(sitrep[[1]]$date_numeric, sitrep[[age_this]]$new_discharges, pch = 20)

# plot open cases in general ward
plot(df_check$general_prevalence)
points(sitrep[[1]]$date_numeric, sitrep[[age_this]]$total_general, pch = 20)

# plot open cases in critical ward
plot(df_check$critical_prevalence)
points(sitrep[[1]]$date_numeric, sitrep[[age_this]]$total_hdu_icu, pch = 20)

#abline(h = 47/111 * 0.2 * 1e4)


