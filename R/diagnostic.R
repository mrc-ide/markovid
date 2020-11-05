#------------------------------------------------
# get Gelman-Rubin diagnostic value from parameter matrix
gelman_rubin <- function(par_draws, par_chain, chains) {
  
  # get number of samples
  samples <- length(par_draws)
  
  # get mean over all samples
  all_mean <- mean(par_draws, na.rm = TRUE)
  
  # get mean of each chain
  chain_mean <- tapply(par_draws, par_chain, mean)
  
  # get variance of each chain
  chain_var <- tapply(par_draws, par_chain, stats::var)
  W <- (1 / chains) * sum(chain_var)
  B <- samples / (chains - 1) * sum((chain_mean - all_mean)^2)
  V <- (1 - 1 / samples) * W + (1 / samples) * B
  
  # compute Gelman-Rubin diagnostic
  gr <- round(sqrt(V / W), 4)
  
  return(gr)
}

#------------------------------------------------
# return autocorrelation for range of lags without plotting
acf_data <- function(x, lag) {
  stats::acf(x, plot = FALSE, lag.max = lag)$acf
}
