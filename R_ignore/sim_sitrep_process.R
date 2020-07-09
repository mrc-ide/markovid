
# run after running sim_sitrep.R. Takes the output of the model and compares C++
# model prediction with manual re-implementation in R.

age_this <- 3

get_delay_density <- function(x, m, s) {
  pgamma(x + 1, shape = 1/s^2, scale = m*s^2) - pgamma(x, shape = 1/s^2, scale = m*s^2)
}
get_delay_tail <- function(x, m, s) {
  pgamma(x + 1, shape = 1/s^2, scale = m*s^2, lower.tail = FALSE)
}

# get true admissions spline
#s <- df_model_fit$admissions_spline[[1]]
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
  open_general <- w[age_this] * s[i] * ( (1 - p_AI[1]) * p_AD[1] * get_delay_tail(j - 1, params["m_AD"], params["s_AD"]) +
                    (1 - p_AI[1]) * (1 - p_AD[1]) * get_delay_tail(j - 1, m_AC[1], params["s_AC"]) +
                    p_AI[1] * get_delay_tail(j - 1, params["m_AI"], params["s_AI"]) )
  
  # new stepup
  new_stepup[i:length(s)] <- new_stepup[i:length(s)] +
    w[age_this] * s[i] * p_AI[1] * get_delay_density(j - 1, params["m_AI"], params["s_AI"])
  
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
tmp <- subset(df_model_fit$by_age, age == age_this & metric == "admission_incidence")
plot(tmp$x, tmp$value)
lines(tmp$x, df_check$admission_incidence)
points(sitrep[[1]]$date_numeric, sitrep[[age_this]]$daily_influx, pch = 20)

# plot deaths
tmp <- subset(df_model_fit$by_age, age == age_this & metric == "deaths_incidence")
plot(tmp$x, tmp$value)
lines(tmp$x, df_check$deaths_incidence)
points(sitrep[[1]]$date_numeric, sitrep[[age_this]]$deaths, pch = 20)

# plot discharges
tmp <- subset(df_model_fit$by_age, age == age_this & metric == "discharges_incidence")
plot(tmp$x, tmp$value)
lines(tmp$x, df_check$discharges_incidence)
points(sitrep[[1]]$date_numeric, sitrep[[age_this]]$new_discharges, pch = 20)

# plot open cases in general ward
tmp <- subset(df_model_fit$by_age, age == age_this & metric == "general_prevalence")
plot(tmp$x, tmp$value)
lines(tmp$x, df_check$general_prevalence)
points(sitrep[[1]]$date_numeric, sitrep[[age_this]]$total_general, pch = 20)

# plot open cases in critical ward
tmp <- subset(df_model_fit$by_age, age == age_this & metric == "critical_prevalence")
plot(tmp$x, tmp$value, ylim = c(0,300))
lines(tmp$x, df_check$critical_prevalence)
points(sitrep[[1]]$date_numeric, sitrep[[age_this]]$total_hdu_icu, pch = 20)


