
#------------------------------------------------
#' @title Simulate individual-level data
#'
#' @description Draw individual-level data at random from the basic hospital
#'   progression model.
#'
#' @param params_scalar vector of scalar parameters.
#' @param params_age list of age-varying parameters, each of which is a vector
#'   over integer ages.
#' @param age_vec the vector of ages corresponding to the age-varying
#'   parameters.
#' @param df_sim dataframe of basic data properties. Must include the columns
#'   "date_admission", "date_censor" and "age".
#'
#' @importFrom stats rgamma rbinom
#' @export

sim_indlevel <- function(params_scalar,
                         params_age,
                         age_vec,
                         df_sim) {
  
  # check inputs
  assert_vector_numeric(params_scalar)
  assert_in(c("s_AI", "s_AD", "s_AC", "s_ID", "s_I1S", "s_I2S", "s_SC", "s_SD"),
            names(params_scalar))
  assert_list(params_age)
  assert_in(c("p_AI", "p_AD", "p_ID", "p_SD", "m_AI", "m_AC", "m_AD", "m_ID", "m_I1S", "m_I2S", "m_SC", "m_SD"),
            names(params_age))
  assert_vector_pos_int(age_vec, zero_allowed = TRUE)
  assert_dataframe(df_sim)
  assert_in(c("date_admission", "date_censor", "age"),
            names(df_sim))
  
  # extract scalar parameters
  s_AI <- params_scalar["s_AI"]
  s_AD <- params_scalar["s_AD"]
  s_AC <- params_scalar["s_AC"]
  s_ID <- params_scalar["s_ID"]
  s_I1S <- params_scalar["s_I1S"]
  s_I2S <- params_scalar["s_I2S"]
  s_SC <- params_scalar["s_SC"]
  s_SD <- params_scalar["s_SD"]
  
  # extract age-varying parameters
  p_AI <- params_age$p_AI
  p_AD <- params_age$p_AD
  p_ID <- params_age$p_ID
  p_SD <- params_age$p_SD
  m_AI <- params_age$m_AI
  m_AC <- params_age$m_AC
  m_AD <- params_age$m_AD
  m_ID <- params_age$m_ID
  m_I1S <- params_age$m_I1S
  m_I2S <- params_age$m_I2S
  m_SC <- params_age$m_SC
  m_SD <- params_age$m_SD
  
  # further checks on inputs
  assert_single_pos(s_AI)
  assert_single_pos(s_AD)
  assert_single_pos(s_AC)
  assert_single_pos(s_ID)
  assert_single_pos(s_I1S)
  assert_single_pos(s_I2S)
  assert_single_pos(s_SC)
  assert_single_pos(s_SD)
  assert_vector_bounded(p_AI)
  assert_vector_bounded(p_AD)
  assert_vector_bounded(p_ID)
  assert_vector_bounded(p_SD)
  assert_vector_pos(m_AI)
  assert_vector_pos(m_AC)
  assert_vector_pos(m_AD)
  assert_vector_pos(m_ID)
  assert_vector_pos(m_I1S)
  assert_vector_pos(m_I2S)
  assert_vector_pos(m_SC)
  assert_vector_pos(m_SD)
  
  assert_int(df_sim$date_admission)
  assert_int(df_sim$date_censor)
  assert_int(df_sim$age)
  
  # --------------------------------
  
  # get data size
  n <- nrow(df_sim)
  
  # match ages against the age_vec
  w_age <- match(df_sim$age, age_vec)
  
  # draw interval times
  t_AI = round(rgamma(n, shape = 1/s_AI^2, scale = m_AI[w_age]*s_AI^2))
  t_AD = round(rgamma(n, shape = 1/s_AD^2, scale = m_AD[w_age]*s_AD^2))
  t_AC = round(rgamma(n, shape = 1/s_AC^2, scale = m_AC[w_age]*s_AC^2))
  t_ID = round(rgamma(n, shape = 1/s_ID^2, scale = m_ID[w_age]*s_ID^2))
  t_I1S = round(rgamma(n, shape = 1/s_I1S^2, scale = m_I1S[w_age]*s_I1S^2))
  t_I2S = round(rgamma(n, shape = 1/s_I2S^2, scale = m_I2S[w_age]*s_I2S^2))
  t_SC = round(rgamma(n, shape = 1/s_SC^2, scale = m_SC[w_age]*s_SC^2))
  t_SD = round(rgamma(n, shape = 1/s_SD^2, scale = m_SD[w_age]*s_SD^2))
  
  # draw progression routes
  route_AI = as.logical(rbinom(n, 1, p_AI[w_age]))
  route_AD = as.logical(rbinom(n, 1, p_AD[w_age]))
  route_ID = as.logical(rbinom(n, 1, p_ID[w_age]))
  route_SD = as.logical(rbinom(n, 1, p_SD[w_age]))
  
  w_AI <- which(route_AI)
  w_AD <- which(!route_AI & route_AD)
  w_AC <- which(!route_AI & !route_AD)
  w_ID <- which(route_AI & route_ID)
  w_IS <- which(route_AI & !route_ID)
  w_I1S <- which(route_AI & !route_ID & route_SD)
  w_I2S <- which(route_AI & !route_ID & !route_SD)
  
  # fill in progression routes
  df_sim$icu <- route_AI
  df_sim$stepdown <- NA
  df_sim$stepdown[w_IS] <- TRUE
  df_sim$stepdown[w_ID] <- FALSE
  
  # fill in final outcomes
  df_sim$final_outcome <- NA
  df_sim$final_outcome[c(w_AD, w_ID, w_I1S)] <- "death"
  df_sim$final_outcome[c(w_AC, w_I2S)] <- "discharge"
  
  # fill in times
  df_sim$date_icu <- NA
  df_sim$date_icu[w_AI] <- df_sim$date_admission[w_AI] + t_AI[w_AI]
  
  df_sim$date_stepdown <- NA
  df_sim$date_stepdown[w_I1S] <- df_sim$date_admission[w_I1S] + t_AI[w_I1S] + t_I1S[w_I1S]
  df_sim$date_stepdown[w_I2S] <- df_sim$date_admission[w_I2S] + t_AI[w_I2S] + t_I2S[w_I2S]
  
  df_sim$date_final_outcome <- NA
  df_sim$date_final_outcome[w_AD] <- df_sim$date_admission[w_AD] + t_AD[w_AD]
  df_sim$date_final_outcome[w_AC] <- df_sim$date_admission[w_AC] + t_AC[w_AC]
  df_sim$date_final_outcome[w_ID] <- df_sim$date_icu[w_ID] + t_ID[w_ID]
  df_sim$date_final_outcome[w_I1S] <- df_sim$date_stepdown[w_I1S] + t_SD[w_I1S]
  df_sim$date_final_outcome[w_I2S] <- df_sim$date_stepdown[w_I2S] + t_SC[w_I2S]
  
  # apply censoring
  w <- which(df_sim$date_icu > df_sim$date_censor)
  df_sim$icu[w] <- NA
  df_sim$date_icu[w] <- NA
  
  w <- which(df_sim$date_stepdown > df_sim$date_censor)
  df_sim$stepdown[w] <- NA
  df_sim$date_stepdown[w] <- NA
  
  w <- which(df_sim$date_final_outcome > df_sim$date_censor)
  df_sim$final_outcome[w] <- NA
  df_sim$date_final_outcome[w] <- NA
  
  return(df_sim)
}

#------------------------------------------------
#' @title Aggregate individual-level data
#'
#' @description Given a line list of individual-level data (for example
#'   generated by \code{sim_indlevel()}), aggregate data - in some cases broken
#'   down by age.
#'
#' @details Transmision probabilities can be derived from the number of patients
#'   going down each pathway (numerator) and the total number of patients
#'   (denominator) at each stage. These two outputs are calculated and returned
#'   for each transition probability. For lengths of stay in each state, we
#'   require a different tabulation in terms of the number of patients that were
#'   in a given state for x days. Transition probabilities are returned broken
#'   down by age, but in this version of the package durations are not broken
#'   down by age.
#'   
#' @param df_data dataframe of individual-level data, in the format output by
#'   \code{sim_indlevel()}.
#' @param age_vec an integer sequence of ages over which to aggregate.
#' @param t_max the maximum time considered when tabulating durations in each
#'   state.
#'
#' @export

aggregate_indlevel <- function(df_data,
                               age_vec = 0:100,
                               t_max = 100) {
  
  # check inputs
  assert_dataframe(df_data)
  assert_in(c("age", "icu", "stepdown", "final_outcome", "date_admission", "date_icu", "date_stepdown", "date_final_outcome"), names(df_data))
  assert_vector_pos_int(age_vec, zero_allowed = TRUE)
  assert_greq(length(age_vec), 5)
  assert_single_pos_int(t_max, zero_allowed = FALSE)
  
  # initialise objects for storing aggregate values
  n_age <- length(age_vec)
  p_AI_numer <- p_AI_denom <- p_AD_numer <- p_AD_denom <- p_ID_numer <- p_ID_denom <- p_SD_numer <- p_SD_denom <- rep(NA, n_age)
  
  # get aggregate values for each 1-year age band
  for (i in seq_len(n_age)) {
    
    # ICU counts
    w <- which(df_data$age == age_vec[i])
    p_AI_numer[i] <- sum(df_data$icu[w] == TRUE, na.rm = TRUE)
    p_AI_denom[i] <- sum(df_data$icu[w] %in% c(TRUE, FALSE), na.rm = TRUE)
    
    # death in general ward
    w <- which((df_data$age == age_vec[i]) & (df_data$icu == FALSE))
    p_AD_numer[i] <- sum(df_data$final_outcome[w] == "death", na.rm = TRUE)
    p_AD_denom[i] <- sum(df_data$final_outcome[w] %in% c("death", "discharge"), na.rm = TRUE)
    
    # death in ICU
    w <- which((df_data$age == age_vec[i]) & (df_data$icu == TRUE))
    p_ID_numer[i] <- sum(df_data$stepdown[w] == FALSE, na.rm = TRUE)
    p_ID_denom[i] <- sum(df_data$stepdown[w] %in% c(TRUE, FALSE), na.rm = TRUE)
    
    # death in stepdown
    w <- which((df_data$age == age_vec[i]) & (df_data$stepdown == TRUE))
    p_SD_numer[i] <- sum(df_data$final_outcome[w] == "death", na.rm = TRUE)
    p_SD_denom[i] <- sum(df_data$final_outcome[w] %in% c("death", "discharge"), na.rm = TRUE)
    
  }
  
  # time admission to ICU
  w <- which(df_data$icu == TRUE)
  m_AI_count <- tabulate(df_data$date_icu[w] - df_data$date_admission[w] + 1, nbins = 100)
  
  # time admission to death in general ward
  w <- which((df_data$icu == FALSE) & (df_data$final_outcome == "death"))
  m_AD_count <- tabulate(df_data$date_final_outcome[w] - df_data$date_admission[w] + 1, nbins = 100)
  
  # time admission to discharge in general ward
  w <- which((df_data$icu == FALSE) & (df_data$final_outcome == "discharge"))
  m_AC_count <- tabulate(df_data$date_final_outcome[w] - df_data$date_admission[w] + 1, nbins = 100)
  
  # time admission to death in ICU
  w <- which((df_data$icu == TRUE) & (df_data$final_outcome == "death"))
  m_ID_count <- tabulate(df_data$date_final_outcome[w] - df_data$date_icu[w] + 1, nbins = 100)
  
  # time admission to stepdown (to death) from ICU
  w <- which((df_data$icu == TRUE) & (df_data$stepdown == TRUE) & (df_data$final_outcome == "death"))
  m_I1S_count <- tabulate(df_data$date_stepdown[w] - df_data$date_icu[w] + 1, nbins = 100)
  
  # time admission to stepdown (to discharge) from ICU
  w <- which((df_data$icu == TRUE) & (df_data$stepdown == TRUE) & (df_data$final_outcome == "discharge"))
  m_I2S_count <- tabulate(df_data$date_stepdown[w] - df_data$date_icu[w] + 1, nbins = 100)
  
  # time stepdown to death
  w <- which((df_data$icu == TRUE) & (df_data$stepdown == TRUE) & (df_data$final_outcome == "death"))
  m_SD_count <- tabulate(df_data$date_final_outcome[w] - df_data$date_stepdown[w] + 1, nbins = 100)
  
  # time stepdown to discharge
  w <- which((df_data$icu == TRUE) & (df_data$stepdown == TRUE) & (df_data$final_outcome == "discharge"))
  m_SC_count <- tabulate(df_data$date_final_outcome[w] - df_data$date_stepdown[w] + 1, nbins = 100)
  
  
  # return as list
  ret <- list(p_AI_numer = p_AI_numer,
              p_AI_denom = p_AI_denom,
              p_AD_numer = p_AD_numer,
              p_AD_denom = p_AD_denom,
              p_ID_numer = p_ID_numer,
              p_ID_denom = p_ID_denom,
              p_SD_numer = p_SD_numer,
              p_SD_denom = p_SD_denom,
              m_AI_count = m_AI_count,
              m_AD_count = m_AD_count,
              m_AC_count = m_AC_count,
              m_ID_count = m_ID_count,
              m_I1S_count = m_I1S_count,
              m_I2S_count = m_I2S_count,
              m_SD_count = m_SD_count,
              m_SC_count = m_SC_count)
  
  return(ret)
}

#------------------------------------------------
#' @title Get cubic spline from node positions
#'
#' @description Given a series of x,y node positions, return a cubic spline
#'   through these points.
#'
#' @param x,y coordinates of spline nodes.
#' @param x_pred x-coordinates at which to evaluate and return cubic spline.
#'
#' @export

cubic_spline <- function (x, y, x_pred) {
  
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


#------------------------------------------------
#' @title Get cubic splines from posterior draws
#'
#' @description Given a matrix of posterior draws of node locations, calculate
#'   cubic splines over all draws, then apply logistic transform to get
#'   transformed spline in a finite range.
#'
#' @param mcmc_samples matrix of posterior draws.
#' @param nodex x-coordinates of nodes.
#' @param age_vec x-coordinates at which to evaluate cubic splines.
#' @param scale posterior splines are logistic transformed and scaled to the
#'   interval \code{[0,scale]}.
#'
#' @export

# function for getting transformed spline from raw nodes
get_spline <- function(mcmc_samples, nodex, age_vec, scale = 1) {
  mapply(function(i) {
    ret <- cubic_spline(nodex, mcmc_samples[i,], age_vec)
    scale / (1 + exp(-ret))
  }, seq_len(nrow(mcmc_samples)))
}

#------------------------------------------------
#' @title Get 95\% quantiles over posterior splines
#'
#' @description Given output from \code{get_spline()}, produce 95\% quantiles
#'   over all draws.
#'
#' @param spline_draws series of draws of posterior splines, as returned by the
#'   \code{get_spline()} function.
#'
#' @export

get_spline_quantiles <- function(spline_draws) {
  ret <- as.data.frame(t(apply(spline_draws, 1, quantile_95)))
  ret$age = seq_len(nrow(ret)) - 1
  ret
}

#------------------------------------------------
#' @title Get 95\% exact binomial intervals from data
#'
#' @description Given line list data, return exact binomial 95\% confidence
#'   intervals of a given outcome vs the alternative. The numerator in the
#'   interval calculation is the number of patients outcome1, and the
#'   denominator is the number of patients in either outcome1 or outcome2.
#'
#' @param age vector of ages for each patient.
#' @param max_age maximum age used in tabulation.
#' @param outcome1,outcome2 for each patient, a logical TRUE/FALSE as to whether
#'   they satisfy outcome1 or outcome2.
#'
#' @export

get_data_quantiles_p <- function(age, max_age, outcome1, outcome2) {
  
  # check inputs
  assert_vector_pos_int(age, zero_allowed = TRUE)
  assert_single_pos_int(max_age, zero_allowed = FALSE)
  assert_vector(outcome1)
  assert_logical(outcome1)
  assert_vector(outcome2)
  assert_logical(outcome2)
  
  # get dataframe of quantiles by age
  ret <- as.data.frame(t(mapply(function(i) {
    w <- which(age == i-1)
    n1 <- sum(outcome1[w] == TRUE, na.rm = TRUE)
    n2 <- sum(outcome2[w] == TRUE, na.rm = TRUE)
    if (n1 == 0 & n2 == 0) {
      return(c(i-1,rep(NA, 3)))
    }
    tmp <- epitools::binom.exact(n1, n1 + n2, conf.level = 0.95)
    c(i-1, tmp$lower, tmp$upper, n1 / (n1 + n2))
  }, 1:(max_age + 1))))
  names(ret) <- c("age", "lower", "upper", "mean")
  
  ret
}

#------------------------------------------------
#' @title Get 95\% confidence intervals from data
#'
#' @description Given line list data, return 95\% confidence intervals of time
#'   to a given outcome.
#'
#' @param age vector of ages for each patient.
#' @param max_age maximum age used in tabulation.
#' @param delta for each patient, for a given outcome, the number of days
#'   between start and end dates.
#' 
#' @importFrom stats sd
#' @export

get_data_quantiles_m <- function(age, max_age, delta) {
  
  # check inputs
  assert_vector_pos_int(age, zero_allowed = TRUE)
  assert_single_pos_int(max_age, zero_allowed = FALSE)
  assert_vector_pos_int(delta, zero_allowed = TRUE)
  
  # get dataframe of quantiles by age
  ret <- as.data.frame(t(mapply(function(i) {
    w <- which(age == i-1)
    delta <- delta[w][!is.na(delta[w])]
    if (length(delta) == 0) {
      return(c(i-1,rep(NA, 3)))
    }
    m <- mean(delta)
    s <- sd(delta) / sqrt(length(delta))
    c(i-1, m - 1.96*s, m + 1.96*s, m)
  }, 1:(max_age + 1))))
  names(ret) <- c("age", "lower", "upper", "mean")
  
  ret
}
