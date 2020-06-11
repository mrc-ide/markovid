
# model_functions.R
#
# Author: Bob Verity
# Date: 2020-04-03
#
# Purpose:
# Set of functions related to the progression model.
#
# ------------------------------------------------------------------

# ------------------------------------------------------------------
# get complementary cumulative distribution quantiles from parameters
get_ccdf <- function(m, s) {
  
  # define quantiles
  q <- seq(0, 1, l = 21)
  q[21] <- 0.99
  
  # get quantiles of all distributions
  ret <- as.data.frame(t(mapply(function(i) {
    qgamma(q, shape = 1/s[i]^2, rate = 1/(m[i]*s[i]^2))
  }, seq_along(m))))
  names(ret) <- q
  
  return(ret)
}

# ------------------------------------------------------------------
# simulate individual-level data from model, using real dates of admission and censoring
sim_indlevel <- function(params, date_admission, date_censor, age_group,
                             n_age, n_samp = 1e3) {
  
  # extract parameters
  p_AI <- params[sprintf("p_AI%s", 1:n_age)]
  p_AD <- params[sprintf("p_AD%s", 1:n_age)]
  p_ID <- params[sprintf("p_ID%s", 1:n_age)]
  m_AI <- params["m_AI"]
  m_AD <- params["m_AD"]
  m_AC <- params[sprintf("m_AC%s", 1:n_age)]
  #m_AC <- params["m_AC"]
  m_ID <- params["m_ID"]
  m_IS <- params["m_IS"]
  m_SC <- params["m_SC"]
  s_AI <- params["s_AI"]
  s_AD <- params["s_AD"]
  s_AC <- params["s_AC"]
  s_ID <- params["s_ID"]
  s_IS <- params["s_IS"]
  s_SC <- params["s_SC"]
  
  # sample from data with replacement
  n <- length(date_admission)
  s <- sample(n, n_samp, replace = TRUE)
  df_sim <- data.frame(date_admission = date_admission[s],
                       date_censor = date_censor[s],
                       age_group = age_group[s])
  
  # draw interval times
  t_AI = rgamma(n_samp, shape = 1/s_AI^2, scale = m_AI*s_AI^2)
  t_AD = rgamma(n_samp, shape = 1/s_AD^2, scale = m_AD*s_AD^2)
  t_AC = rgamma(n_samp, shape = 1/s_AC^2, scale = m_AC[df_sim$age_group]*s_AC^2)
  #t_AC = rgamma(n_samp, shape = 1/s_AC^2, scale = m_AC*s_AC^2)
  t_ID = rgamma(n_samp, shape = 1/s_ID^2, scale = m_ID*s_ID^2)
  t_IS = rgamma(n_samp, shape = 1/s_IS^2, scale = m_IS*s_IS^2)
  t_SC = rgamma(n_samp, shape = 1/s_SC^2, scale = m_SC*s_SC^2)
  
  # draw progression routes
  route_AI = as.logical(rbinom(n_samp, 1, p_AI[df_sim$age_group]))
  route_AD = as.logical(rbinom(n_samp, 1, p_AD[df_sim$age_group]))
  route_ID = as.logical(rbinom(n_samp, 1, p_ID[df_sim$age_group]))
  
  w_AI <- which(route_AI)
  w_AD <- which(!route_AI & route_AD)
  w_AC <- which(!route_AI & !route_AD)
  w_ID <- which(route_AI & route_ID)
  w_IC <- which(route_AI & !route_ID)
  
  # fill in progression routes
  df_sim$icu <- route_AI
  df_sim$stepdown <- NA
  df_sim$stepdown[w_IC] <- TRUE
  df_sim$stepdown[w_ID] <- FALSE
  
  # fill in final outcomes
  df_sim$final_outcome <- NA
  df_sim$final_outcome[c(w_AD, w_ID)] <- "death"
  df_sim$final_outcome[c(w_AC, w_IC)] <- "discharge"
  
  # fill in times
  df_sim$date_icu <- NA
  df_sim$date_icu[w_AI] <- df_sim$date_admission[w_AI] + t_AI[w_AI]
  
  df_sim$date_leave_icu <- NA
  df_sim$date_leave_icu[w_IC] <- df_sim$date_admission[w_IC] + t_AI[w_IC] + t_IS[w_IC]
  
  df_sim$date_final_outcome <- NA
  df_sim$date_final_outcome[w_AD] <- df_sim$date_admission[w_AD] + t_AD[w_AD]
  df_sim$date_final_outcome[w_AC] <- df_sim$date_admission[w_AC] + t_AC[w_AC]
  df_sim$date_final_outcome[w_ID] <- df_sim$date_icu[w_ID] + t_ID[w_ID]
  df_sim$date_final_outcome[w_IC] <- df_sim$date_leave_icu[w_IC] + t_SC[w_IC]
  
  # apply censoring
  w <- which(df_sim$date_icu > df_sim$date_censor)
  df_sim$icu[w] <- NA
  df_sim$date_icu[w] <- NA
  
  w <- which(df_sim$date_leave_icu > df_sim$date_censor)
  df_sim$stepdown[w] <- NA
  df_sim$date_leave_icu[w] <- NA
  
  w <- which(df_sim$date_final_outcome > df_sim$date_censor)
  df_sim$final_outcome[w] <- NA
  df_sim$date_final_outcome[w] <- NA
  
  return(df_sim)
}

# ------------------------------------------------------------------
# get all delay distributions binned into daily proportions
bin_indlevel <- function(indlevel, max_days = 30) {
  
  ret <- NULL
  
  # admission to ICU
  tab1 <- tabulate(indlevel$date_icu - indlevel$date_admission + 1,
                   nbins = max_days + 1)
  ret <- rbind(ret, data.frame(day = 0:max_days,
                               metric = "admission to ICU",
                               proportion = tab1 / sum(tab1)))
  
  # general ward to death
  w <- which(!indlevel$icu & indlevel$final_outcome == "death")
  tab1 <- tabulate(indlevel$date_final_outcome[w] - indlevel$date_admission[w] + 1,
                   nbins = max_days + 1)
  ret <- rbind(ret, data.frame(day = 0:max_days,
                               metric = "general ward to death",
                               proportion = tab1 / sum(tab1)))
  
  # general ward to discharge
  w <- which(!indlevel$icu & indlevel$final_outcome == "discharge")
  tab1 <- tabulate(indlevel$date_final_outcome[w] - indlevel$date_admission[w] + 1,
                   nbins = max_days + 1)
  ret <- rbind(ret, data.frame(day = 0:max_days,
                               metric = "general ward to discharge",
                               proportion = tab1 / sum(tab1)))
  
  # ICU to death
  w <- which(indlevel$icu & indlevel$final_outcome == "death")
  tab1 <- tabulate(indlevel$date_final_outcome[w] - indlevel$date_icu[w] + 1,
                   nbins = max_days + 1)
  ret <- rbind(ret, data.frame(day = 0:max_days,
                               metric = "ICU to death",
                               proportion = tab1 / sum(tab1)))
  
  # ICU to stepdown
  w <- which(indlevel$icu & indlevel$stepdown)
  tab1 <- tabulate(indlevel$date_leave_icu[w] - indlevel$date_icu[w] + 1,
                   nbins = max_days + 1)
  ret <- rbind(ret, data.frame(day = 0:max_days,
                               metric = "ICU to stepdown",
                               proportion = tab1 / sum(tab1)))
  
  # stepdown to discharge
  w <- which(indlevel$icu & indlevel$final_outcome == "discharge")
  tab1 <- tabulate(indlevel$date_final_outcome[w] - indlevel$date_leave_icu[w] + 1,
                   nbins = max_days + 1)
  ret <- rbind(ret, data.frame(day = 0:max_days,
                               metric = "stepdown to discharge",
                               proportion = tab1 / sum(tab1)))
  
  return(ret)
}

# ------------------------------------------------------------------
# convert nested list to long dataframe
nested_to_long <- function(nl) {
  
  do.call(rbind, mapply(function(k) {
    cbind(do.call(rbind, mapply(function(j) {
      cbind(do.call(rbind, mapply(function(i) {
        data.frame(x = seq_along(nl[[k]][[j]][[i]]),
                   value = nl[[k]][[j]][[i]],
                   i = i)
      }, seq_along(nl[[k]][[j]]), SIMPLIFY = FALSE)), j = j)
    }, seq_along(nl[[k]]), SIMPLIFY = FALSE)), k = k)
  }, seq_along(nl), SIMPLIFY = FALSE))
  
}
