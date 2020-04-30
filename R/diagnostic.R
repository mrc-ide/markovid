#------------------------------------------------
gelman_rubin <- function(par_matrix, chains, samples){
  par_matrix <- as.data.frame(par_matrix)
  # Mean over all samples
  all_mean <- mean(par_matrix[,2])
  
  # Mean of each chain
  chain_mean <- tapply(par_matrix[,2], par_matrix[,1], mean)
  
  # Variance of each chain
  chain_var <- tapply(par_matrix[,2], par_matrix[,1], stats::var)
  W <- (1 / chains) * sum(chain_var)
  B <- samples / (chains - 1) * sum((chain_mean - all_mean)^2)
  V <- (1 - 1 / samples) * W + (1 / samples) * B
  round(sqrt(V / W), 4)
}

#------------------------------------------------
acf_data <- function(x, lag){
  stats::acf(x, plot = FALSE, lag.max = lag)$acf
}
