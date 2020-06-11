
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
sim_indlevel <- function(params, p_AI, p_AD, p_ID, m_AC,
                         date_admission, date_censor, age, n_samp = 1e3) {
  
  # extract parameters
  m_AI <- params["m_AI"]
  m_AD <- params["m_AD"]
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
                       age = age[s])
  
  # draw interval times
  t_AI = rgamma(n_samp, shape = 1/s_AI^2, scale = m_AI*s_AI^2)
  t_AD = rgamma(n_samp, shape = 1/s_AD^2, scale = m_AD*s_AD^2)
  t_AC = rgamma(n_samp, shape = 1/s_AC^2, scale = m_AC[df_sim$age]*s_AC^2)
  t_ID = rgamma(n_samp, shape = 1/s_ID^2, scale = m_ID*s_ID^2)
  t_IS = rgamma(n_samp, shape = 1/s_IS^2, scale = m_IS*s_IS^2)
  t_SC = rgamma(n_samp, shape = 1/s_SC^2, scale = m_SC*s_SC^2)
  
  # draw progression routes
  route_AI = as.logical(rbinom(n_samp, 1, p_AI[df_sim$age]))
  route_AD = as.logical(rbinom(n_samp, 1, p_AD[df_sim$age]))
  route_ID = as.logical(rbinom(n_samp, 1, p_ID[df_sim$age]))
  
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

# ------------------------------------------------------------------
# get cubic spline
cubic_spline <- function (x, y, x_pred) {
  assert_vector_numeric(x)
  assert_increasing(x)
  assert_vector_numeric(y)
  assert_same_length(x, y)
  assert_vector_numeric(x_pred)
  assert_increasing(x_pred)
  assert_greq(min(x_pred), min(x))
  assert_leq(max(x_pred), max(x))
  n <- length(x) - 1
  c <- l <- mu <- z <- rep(0, n + 1)
  h <- b <- d <- alpha <- rep(NA, n)
  for (i in 1:n) {
    h[i] <- x[i + 1] - x[i]
  }
  for (i in 2:n) {
    alpha[i] <- 3/h[i] * (y[i + 1] - y[i]) - 3/h[i - 1] * 
      (y[i] - y[i - 1])
  }
  l[1] <- 1
  for (i in 2:n) {
    l[i] <- 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 
                                                        1]
    mu[i] <- h[i]/l[i]
    z[i] <- (alpha[i] - h[i - 1] * z[i - 1])/l[i]
  }
  l[n + 1] <- 1
  for (i in n:1) {
    c[i] <- z[i] - mu[i] * c[i + 1]
    b[i] <- (y[i + 1] - y[i])/h[i] - h[i] * (c[i + 1] + 2 * 
                                               c[i])/3
    d[i] <- (c[i + 1] - c[i])/(3 * h[i])
  }
  s <- rep(NA, length(x_pred))
  j <- 1
  for (i in seq_along(x_pred)) {
    if (x_pred[i] > x[j + 1]) {
      j <- j + 1
    }
    s[i] <- y[j] + b[j] * (x_pred[i] - x[j]) + c[j] * (x_pred[i] - 
                                                         x[j])^2 + d[j] * (x_pred[i] - x[j])^3
  }
  return(s)
}

