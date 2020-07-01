
#------------------------------------------------
run_mcmc <- function(data_list,
                     df_params,
                     burnin = 1e3,
                     samples = 1e4,
                     beta_vec = 1,
                     chains = 1,
                     pb_markdown = FALSE,
                     silent = FALSE,
                     return_fit = FALSE) {
  
  # ---------- check inputs ----------
  
  # check df_params
  assert_dataframe(df_params)
  assert_in(c("name", "min", "max", "init"), names(df_params),
            message = "df_params must contain the columns 'name', 'min', 'max', 'init")
  assert_numeric(df_params$min)
  assert_numeric(df_params$max)
  assert_leq(df_params$min, df_params$max)
  assert_numeric(df_params$init)
  assert_greq(df_params$init, df_params$min)
  assert_leq(df_params$init, df_params$max)
  
  # check MCMC parameters
  assert_single_pos_int(burnin, zero_allowed = FALSE)
  assert_single_pos_int(samples, zero_allowed = FALSE)
  assert_single_pos_int(chains, zero_allowed = FALSE)
  
  # check misc parameters
  assert_single_logical(pb_markdown)
  assert_single_logical(silent)
  
  
  # ---------- pre-processing ----------
  
  # get number of rungs
  rungs <- length(beta_vec)
  
  # calculate transformation type for each parameter
  # 0 = [-Inf,Inf] -> phi = theta
  # 1 = [-Inf,b]   -> phi = log(b - theta)
  # 2 = [a,Inf]    -> phi = log(theta - a)
  # 3 = [a,b]      -> phi = log((theta - a)/(b - theta))
  df_params$trans_type <- 2*is.finite(df_params$min) + is.finite(df_params$max)
  
  # flag to skip over fixed parameters
  skip_param <- (df_params$min == df_params$max)
  
  
  # ---------- define argument lists ----------
  
  # parameters to pass to C++
  args_params <- list(data_list = data_list,
                      theta_min = df_params$min,
                      theta_max = df_params$max,
                      theta_init = df_params$init,
                      trans_type = df_params$trans_type,
                      skip_param = skip_param,
                      burnin = burnin,
                      samples = samples,
                      beta_vec = beta_vec,
                      pb_markdown = pb_markdown,
                      silent = silent,
                      return_fit = return_fit)
  
  # functions to pass to C++
  args_functions <- list(test_convergence = test_convergence,
                         update_progress = update_progress)
  
  # complete list of arguments
  args <- list(args_params = args_params,
               args_functions = args_functions)
  
  # replicate arguments over chains
  parallel_args <- replicate(chains, args, simplify = FALSE)
  for (i in 1:chains) {
    parallel_args[[i]]$args_params$chain <- i
  }
  
  
  # ---------- run MCMC ----------
  
  # run in serial
  output_raw <- lapply(parallel_args, deploy_chain)
  
  
  # ---------- option to return model fit ----------
  
  if (return_fit) {
    
    # extract elements that are not defined over ages
    admissions_spline <- output_raw[[1]]$admissions_spline
    
    # drop elements that are not defined over ages
    raw_names <- names(output_raw[[1]])
    output_raw[[1]] <- output_raw[[1]][setdiff(raw_names, "admissions_spline")]
    
    # get remaining into long form
    by_age <- nested_to_long(output_raw[[1]])
    names(by_age) <- c("x", "value", "age", "region", "metric")
    by_age$metric <- c("admission_incidence",
                       "deaths_incidence",
                       "discharges_incidence",
                       "general_prevalence",
                       "critical_prevalence")[by_age$metric]
    
    ret <- list(admissions_spline = admissions_spline,
                by_age = by_age)
    return(ret)
    
  }
  
  # ---------- process output ----------
  
  # define names
  chain_names <- sprintf("chain%s", 1:chains)
  rung_names <- sprintf("rung%s", 1:rungs)
  param_names <- df_params$name
  
  # get raw output into dataframe
  df_output <- do.call(rbind, mapply(function(j) {
    do.call(rbind, mapply(function(i) {
      
      # concatenate burn-in and sampling loglike and logprior
      logprior <- c(output_raw[[j]]$logprior_burnin[[i]], output_raw[[j]]$logprior_sampling[[i]])
      loglike <- c(output_raw[[j]]$loglike_burnin[[i]], output_raw[[j]]$loglike_sampling[[i]])
      
      # create dataframe of loglike and logprior
      ret <- data.frame(chain = chain_names[j],
                        rung = rung_names[i],
                        iteration = seq_along(loglike),
                        stage = rep(c("burnin", "sampling"), times = c(burnin, samples)),
                        logprior = logprior,
                        loglikelihood = loglike)
      
      # concatenate theta into dataframe and append
      theta <- as.data.frame(do.call(rbind, c(output_raw[[j]]$theta_burnin[[i]], output_raw[[j]]$theta_sampling[[i]])))
      names(theta) <- param_names
      ret <- cbind(ret, theta)
      
      return(ret)
    }, seq_along(output_raw[[j]]$loglike_burnin), SIMPLIFY = FALSE))
  }, seq_along(output_raw), SIMPLIFY = FALSE))
  
  # append to output list
  output_processed <- list(output = df_output)
  output_processed$diagnostics <- list()
  
  ## Diagnostics
  # Rhat
  if (chains > 1) {
    rhat_est <- c()
    for(p in seq_along(param_names)){
      pm <- output_processed$output[output_processed$output$stage == "sampling",c("chain", as.character(param_names[p]))]
      rhat_est[p] <- gelman_rubin(pm, chains, samples)
    }
    rhat_est[skip_param] <- NA
    output_processed$diagnostics$rhat <- rhat_est
  }
  
  # MC
  if (rungs > 1) {
    
    # Beta raised
    output_processed$diagnostics$beta <- tidyr::expand_grid(chain = chain_names, rung = rung_names)
    output_processed$diagnostics$beta$value <- unlist(lapply(output_raw, function(x) x$beta))
    
    # MC accept
    mc_accept <- tidyr::expand_grid(chain = chain_names, link = 1:(length(rung_names) - 1))
    mc_accept$burnin <- unlist(lapply(output_raw, function(x){x$mc_accept_burnin})) / burnin
    mc_accept$sampling <- unlist(lapply(output_raw, function(x){x$mc_accept_sampling})) / samples
    mc_accept <- tidyr::gather(mc_accept, stage, value, -chain, -link)
    
    output_processed$diagnostics$mc_accept <- mc_accept
  }
  
  ## Parameters
  output_processed$parameters <- list(data = data,
                                      df_params = df_params,
                                      burnin = burnin,
                                      samples = samples,
                                      rungs = rungs,
                                      chains = chains)

  # save output as custom class
  class(output_processed) <- "drjacoby_output"
  
  # return
  return(output_processed)
}

#------------------------------------------------
# deploy main_mcmc for this chain
#' @noRd
deploy_chain <- function(args) {
  
  # get parameters
  burnin <- args$args_params$burnin
  samples <- args$args_params$samples
  
  # make progress bars
  pb_burnin <- txtProgressBar(min = 0, max = burnin, initial = NA, style = 3)
  pb_samples <- txtProgressBar(min = 0, max = samples, initial = NA, style = 3)
  args$args_progress <- list(pb_burnin = pb_burnin,
                             pb_samples = pb_samples)
  
  # run C++ function
  ret <- run_mcmc_cpp(args)
  
  return(ret)
}

