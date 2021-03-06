
#------------------------------------------------
#' @title Plot autocorrelation
#'   
#' @description Plot loglikelihood 95\% credible intervals.
#'   
#' @param x an object of class \code{drjacoby_output}
#' @param lag calculate autocorrelation up to this many lags.
#' @param par which parameter to plot.
#' @param chain which chain to plot.
#' @param phase whether to plot from burnin or sampling phase.
#' @param rung which thermodynamic rung to plot.
#'
#' @export

plot_autocorrelation <- function(x, lag = 20, par = NULL, chain = 1, phase = "sampling", rung = 1) {
  
  # check inputs
  assert_custom_class(x, "drjacoby_output")
  assert_single_pos_int(chain)
  assert_leq(chain, length(x))
  assert_in(phase, c("burnin", "sampling"))
  assert_single_bounded(lag, 1, 500)
  
  # declare variables to avoid "no visible binding" issues
  stage <- iteration <- loglikelihood <- NULL
  
  # get values
  chain_get <- paste0("chain", chain)
  rung_get <- paste0("rung", rung)
  data <- dplyr::filter(x$output, chain == chain_get, rung == rung_get, stage == phase) %>%
    dplyr::select(-chain, -rung, -iteration, -stage, -loglikelihood) %>%
    as.data.frame()
  
  # Select parameters
  if(!is.null(par)){
    data <- data[, colnames(data) %in% par, drop = FALSE]
  }
  
  # Estimate autocorrelation
  out <- as.data.frame(apply(data, 2, acf_data, lag = lag))
  
  # Format data for plotting
  out$lag <- 0:lag
  out <- tidyr::gather(out, "parameter", "Autocorrelation", -lag)
  
  ggplot2::ggplot(data = out,
                  ggplot2::aes(x = .data$lag, y = 0, xend = .data$lag, yend =.data$Autocorrelation)) + 
    ggplot2::geom_hline(yintercept = 0, lty = 2, col = "red") + 
    ggplot2::geom_segment(size = 1.5) +
    ggplot2::theme_bw() +
    ggplot2::ylab("Autocorrelation") +
    ggplot2::xlab("Lag") +
    ggplot2::ylim(min(0, min(out$Autocorrelation)), 1) +
    ggplot2::facet_wrap(~ parameter)
}

#------------------------------------------------
#' @title Plot parameter summary
#'   
#' @description Combined plot of posterior summaries of a given parameter or set of parameters.
#'   
#' @param x an object of class \code{drjacoby_output}.
#' @param show which parameters to plot.
#' @param hide which parameters to omit from plot.
#' @param lag calculate autocorrelation up to this many lags.
#' @param downsample whether to downsample iterations.
#' @param phase whether to plot from burnin or sampling phase.
#' @param rung which thermodynamic rung to plot.
#' @param display if FALSE the plot is produced but not saved.
#'
#' @export

plot_par <- function(x, show = NULL, hide = NULL, lag = 20,
                     downsample = TRUE, phase = "sampling", rung = 1,
                     display = TRUE) {
  
  # check inputs
  assert_custom_class(x, "drjacoby_output")
  assert_single_bounded(lag, 1, 500)
  assert_single_logical(downsample)
  assert_in(phase, c("burnin", "sampling", "both"))
  assert_single_pos_int(rung)
  assert_single_logical(display)
  
  # declare variables to avoid "no visible binding" issues
  stage <- chain <- NULL
  
  # deal with phase = "both" situation
  if (phase == "both") {
    phase <- c("burnin", "sampling")
  }
  
  # get basic properties
  rung_get <- paste0("rung", rung)
  data <- dplyr::filter(x$output, rung == rung_get, stage %in% phase) 
  
  # choose which parameters to plot
  parameter <- setdiff(names(data), c("chain", "rung", "iteration", "stage", "logprior", "loglikelihood"))
  if(!is.null(show)){
    stopifnot(is.character(show))
    parameter <- parameter[grepl(paste(show, collapse = "|"), parameter)]
  }
  if(!is.null(hide)){
    stopifnot(is.character(hide))
    parameter <- parameter[!grepl(paste(hide, collapse = "|"), parameter)]
  }
  if(length(parameter) > 10){
    message("More than 10 parameters to summarise, consider using the show or hide arguments 
            to select parameters and reduce computation time.")
  }
  
  # Downsample
  if(downsample & nrow(data) > 2000){
    data <- data[seq.int(1, nrow(data), length.out = 2000),]
  }
  
  data <- dplyr::group_by(data, chain)
  data <- dplyr::mutate(data, plot_par_x = 1:dplyr::n())
  data <- dplyr::ungroup(data)
  
  # Autocorrealtion (on downsample)
  ac_data <- as.data.frame(apply(data[,parameter], 2, acf_data, lag = lag))
  ac_data$lag <- 0:lag
  
  # Set minimum bin number
  b <- min(nrow(data) / 4, 40)
  
  # produce plots over all parameters
  plot_list <- c()
  for(j in 1:length(parameter)){
    pd <- data[, c("chain", "plot_par_x", parameter[j])]
    names(pd) <- c("chain", "plot_par_x", "y")
    pd2 <- ac_data[, c("lag", parameter[j])]
    names(pd2) <- c("lag", "Autocorrelation")
    
    # Histogram
    p1 <- ggplot2::ggplot(pd, ggplot2::aes(x = .data$y)) + 
      ggplot2::geom_histogram(bins = b, fill = "deepskyblue3", col = "darkblue") + 
      ggplot2::ylab("Count") + 
      ggplot2::xlab(parameter[j]) +
      ggplot2::theme_bw()
    
    # Chains
    p2 <- ggplot2::ggplot(pd, ggplot2::aes(x = .data$plot_par_x, y = .data$y, col = .data$chain)) + 
      ggplot2::geom_line() +
      scale_color_discrete(name = "Chain") +
      ggplot2::xlab("Iteration") +
      ggplot2::ylab(parameter[j]) +
      ggplot2::theme_bw()
    
    # Autocorrealtion
    p3 <- ggplot2::ggplot(data = pd2,
                          ggplot2::aes(x = .data$lag, y = 0, xend = .data$lag, yend =.data$Autocorrelation)) + 
      ggplot2::geom_hline(yintercept = 0, lty = 2, col = "red") + 
      ggplot2::geom_segment(size = 1.5) +
      ggplot2::theme_bw() +
      ggplot2::ylab("Autocorrelation") +
      ggplot2::xlab("Lag") +
      ggplot2::ylim(min(0, min(pd2$Autocorrelation)), 1)
    
    # Arrange
    pc1 <- cowplot::plot_grid(p1, p3, ncol = 2)
    pc2 <- cowplot::plot_grid(p2, pc1, nrow = 2)
    plot_list[[j]] <- list(trace = p2,
                           hist = p1,
                           acf = p3,
                           combined = pc2)
  }
  names(plot_list) <- paste0("Plot_", parameter)
  
  if(!display){
    return(invisible(plot_list))
  } else {
    # Display plots, asking user for next page if multiple parameters
    for(j in 1:length(parameter)){
      graphics::plot(plot_list[[j]]$combined)
      if(j == 1){
        default_ask <- grDevices::devAskNewPage()
        on.exit(grDevices::devAskNewPage(default_ask))
        grDevices::devAskNewPage(ask = TRUE)
      }
    }
  }
  
  return(invisible(plot_list))
}

#------------------------------------------------
#' @title Plot parameter autocorrelation
#'   
#' @description Plot parameter autocorrelation.
#'   
#' @param x an object of class \code{drjacoby_output}.
#' @param parameter1,parameter2 which parameters to plot.
#' @param downsample whether to downsample iterations.
#' @param phase whether to plot from burnin or sampling phase.
#' @param rung which thermodynamic rung to plot.
#'
#' @export

plot_cor <- function(x, parameter1, parameter2,
                     downsample = TRUE, phase = "sampling",
                     rung = 1) {
  
  # check inputs
  assert_custom_class(x, "drjacoby_output")
  assert_string(parameter1)
  assert_string(parameter2)
  assert_in(parameter1, names(x$output))
  assert_in(parameter2, names(x$output))
  assert_single_logical(downsample)
  assert_in(phase, c("burnin", "sampling"))
  assert_single_pos_int(rung)
  
  # declare variables to avoid "no visible binding" issues
  stage <- NULL
  
  # get basic quantities
  rung_get <- paste0("rung", rung)
  data <- dplyr::filter(x$output, rung == rung_get, stage == phase) 
  data <- data[,c("chain", parameter1, parameter2)]  
  colnames(data) <- c("chain", "x", "y")
  
  # Downsample
  if(downsample & nrow(data) > 2000){
    data <- data[seq.int(1, nrow(data), length.out = 2000),]
  }
  
  # produce plot
  ggplot2::ggplot(data = data,
                  ggplot2::aes(x = .data$x, y = .data$y, col = as.factor(.data$chain))) + 
    ggplot2::geom_point(alpha = 0.5) + 
    ggplot2::xlab(parameter1) +
    ggplot2::ylab(parameter2) +
    scale_color_discrete(name = "Chain") +
    ggplot2::theme_bw()
  
}

#------------------------------------------------
#' @title Credible interval plot
#'   
#' @description Produce credible interval plot.
#'   
#' @param x an object of class \code{drjacoby_output}.
#' @param show which parameters to plot.
#' @param phase whether to plot from burnin or sampling phase.
#' @param rung which thermodynamic rung to plot.
#' @param param_names manually define parameter names on plot.
#'
#' @export

plot_credible <- function(x, show = NULL, phase = "sampling", rung = 1, param_names = NULL) {
  
  # check inputs
  assert_custom_class(x, "drjacoby_output")
  if (!is.null(show)) {
    assert_string(show)
    assert_in(show, names(x$output))
  }
  assert_in(phase, c("burnin", "sampling", "both"))
  assert_single_pos_int(rung)
  
  # declare variables to avoid "no visible binding" issues
  stage <- NULL
  
  # define defaults
  if (is.null(show)) {
    show <- setdiff(names(x$output), c("chain", "rung", "iteration", "stage", "logprior", "loglikelihood"))
  }
  if (is.null(param_names)) {
    param_names <- show
  }
  
  # deal with phase = "both" situation
  if (phase == "both") {
    phase <- c("burnin", "sampling")
  }
  
  # subset based on phase and rung
  rung_get <- paste0("rung", rung)
  data <- dplyr::filter(x$output, rung == rung_get, stage %in% phase) 
  data <- data[, show, drop = FALSE]
  
  # get quantiles
  df_plot <- as.data.frame(t(apply(data, 2, quantile_95)))
  df_plot$param <- factor(param_names, levels = param_names)
  
  # produce plot
  ggplot2::ggplot(data = df_plot) + ggplot2::theme_bw() +
    ggplot2::geom_point(ggplot2::aes(x = .data$param, y = .data$Q50)) +
    ggplot2::geom_segment(ggplot2::aes(x = .data$param, y = .data$Q2.5, xend = .data$param, yend = .data$Q97.5)) +
    ggplot2::xlab("") +
    ggplot2::ylab("95% CrI")
  
}


#------------------------------------------------
#' @title Produce credible interval plot from dataframe
#'   
#' @description Combined plot of posterior summaries of a given parameter or set of parameters.
#'   
#' @param x a dataframe.
#' @param show which parameters to plot.
#' @param param_names manually define parameter names.
#' @param rotate If TRUE, parameter names are rotated.
#'
#' @export

plot_CrI <- function (x, show = NULL, param_names = NULL, rotate = FALSE) {
  
  # avoid "no visible binding" note
  param <- NULL
  
  # check inputs
  assert_dataframe(x)
  assert_in(c("param", "Q2.5", "Q50", "Q97.5"), names(x))
  
  # subset parameters
  if (is.null(show)) {
    show <- x$param
  }
  df_plot <- subset(x, param %in% show)
  
  # default param names
  if (is.null(param_names)) {
    param_names <- unique(show)
  }
  
  # ensure factor levels ordered correctly
  df_plot$param <- factor(df_plot$param, levels = show)
  
  # create plot
  ret <- ggplot2::ggplot(df_plot) + ggplot2::theme_bw() +
    ggplot2::geom_pointrange(ggplot2::aes_(x = ~param, y = ~Q50, ymin = ~Q2.5, ymax = ~Q97.5)) +
    ggplot2::scale_x_discrete(labels = param_names) +
    ggplot2::xlab("") + ggplot2::ylab("95% CrI")
  
  # optionally rotate axis labels
  if (rotate) {
    ret <- ret + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  }
  
  return(ret)
}

#------------------------------------------------
#' @title Plot loglikelihood 95\% credible intervals
#'   
#' @description Plot loglikelihood 95\% credible intervals.
#'   
#' @param x an object of class \code{drjacoby_output}
#' @param chain which chain to plot.
#' @param phase which phase to plot. Must be either "burnin" or "sampling".
#' @param x_axis_type how to format the x-axis. 1 = integer rungs, 2 = values of
#'   the thermodynamic power.
#' @param y_axis_type how to format the y-axis. 1 = raw values, 2 = truncated at
#'   auto-chosen lower limit. 3 = double-log scale.
#'
#' @export

plot_rung_loglike <- function(x, chain = 1, phase = "sampling", x_axis_type = 1, y_axis_type = 1) {
  
  # check inputs
  assert_custom_class(x, "drjacoby_output")
  assert_single_pos_int(chain)
  assert_leq(chain, length(x))
  assert_in(phase, c("burnin", "sampling"))
  assert_single_pos_int(x_axis_type)
  assert_in(x_axis_type, 1:2)
  assert_single_pos_int(y_axis_type)
  assert_in(y_axis_type, 1:3)
  
  # declare variables to avoid "no visible binding" issues
  stage <- rung <- value <- loglikelihood <- NULL
  
  # get useful quantities
  chain_get <- paste0("chain", chain)
  thermo_power <- x$diagnostics$rung_details$thermodynamic_power
  rungs <- length(thermo_power)
  
  # define x-axis type
  if (x_axis_type == 1) {
    x_vec <- rungs:1
    x_lab <- "rung"
  } else {
    x_vec <- thermo_power
    x_lab <- "thermodynamic power"
  }
  
  # get plotting values (loglikelihoods)
  data <- dplyr::filter(x$output, chain == chain_get, stage == phase)
  y_lab <- "log-likelihood"
  
  # move to plotting deviance if specified
  if (y_axis_type == 3) {
    data$loglikelihood <- -2 * data$loglikelihood
    y_lab <- "deviance"
    
    # if needed, scale by adding/subtracting a power of ten until all values are
    # positive
    if (min(data$loglikelihood) < 0) {
      dev_scale_power <- ceiling(log(abs(min(data$loglikelihood)))/log(10))
      dev_scale_sign <- -sign(min(data$loglikelihood))
      data$loglikelihood <- data$loglikelihood + dev_scale_sign*10^dev_scale_power
      
      dev_scale_base <- ifelse(dev_scale_power == 0, 1, 10)
      dev_scale_power_char <- ifelse(dev_scale_power <= 1, "", paste("^", dev_scale_power))
      dev_scale_sign_char <- ifelse(dev_scale_sign < 0, "-", "+")
      y_lab <- parse(text = paste("deviance", dev_scale_sign_char, dev_scale_base, dev_scale_power_char))
    }
  }
  
  # get 95% credible intervals over plotting values
  y_intervals <- data %>%
    dplyr::group_by(rung) %>%
    dplyr::summarise(Q2.5 = quantile(loglikelihood, 0.025),
                     Q50 =  quantile(loglikelihood, 0.5),
                     Q97.5 = quantile(loglikelihood, 0.975))
  
  # get data into ggplot format and define temperature colours
  df <- y_intervals
  df$col <- thermo_power
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw() + theme(panel.grid.minor.x = element_blank(),
                                           panel.grid.major.x = element_blank())
  plot1 <- plot1 + geom_vline(aes(xintercept = x_vec), col = grey(0.9))
  plot1 <- plot1 + geom_segment(aes_(x = ~x_vec, y = ~Q2.5, xend = ~x_vec, yend = ~Q97.5))
  plot1 <- plot1 + geom_point(aes_(x = ~x_vec, y = ~Q50, color = ~col))
  plot1 <- plot1 + xlab(x_lab) + ylab(y_lab)
  plot1 <- plot1 + scale_colour_gradientn(colours = c("red", "blue"), name = "thermodynamic\npower", limits = c(0,1))
  
  # define y-axis
  if (y_axis_type == 2) {
    y_min <- quantile(df$Q2.5, probs = 0.5)
    y_max <- max(df$Q97.5)
    plot1 <- plot1 + coord_cartesian(ylim = c(y_min, y_max))
  } else if (y_axis_type == 3) {
    plot1 <- plot1 + scale_y_continuous(trans = "log10")
  }
  
  # return plot object
  return(plot1)
}

#------------------------------------------------
#' @title Plot Metropolis coupling acceptance rates
#'
#' @description Plot Metropolis coupling acceptance rates between all rungs.
#'
#' @inheritParams plot_rung_loglike
#'
#' @import ggplot2
#' @importFrom grDevices grey
#' @export

plot_mc_acceptance <- function(x, chain = 1, phase = "sampling", x_axis_type = 1) {
  
  # check inputs
  assert_custom_class(x, "drjacoby_output")
  assert_single_pos_int(chain)
  assert_in(phase, c("burnin", "sampling"))
  assert_single_pos_int(x_axis_type)
  assert_in(x_axis_type, 1:2)
  
  # declare variables to avoid "no visible binding" issues
  stage <- value <- NULL
  
  # get useful quantities
  chain_get <- paste0("chain", chain)
  thermo_power <- x$diagnostics$rung_details$thermodynamic_power
  thermo_power_mid <- thermo_power[-1] - diff(thermo_power)/2
  rungs <- length(thermo_power)
  
  # exit if rungs = 1
  if (rungs == 1) {
    stop("no metropolis coupling when rungs = 1")
  }
  
  # define x-axis type
  if (x_axis_type == 1) {
    breaks_vec <- rungs:2
    x_vec <- (rungs:2) - 0.5
    x_lab <- "rung"
  } else {
    breaks_vec <- thermo_power
    x_vec <- thermo_power_mid
    x_lab <- "thermodynamic power"
  }
  
  # get acceptance rates
  mc_accept <- dplyr::filter(x$diagnostics$mc_accept, stage == phase, chain == chain_get) %>%
    dplyr::pull(value)
  
  # get data into ggplot format and define temperature colours
  df <- as.data.frame(mc_accept)
  df$col <- thermo_power_mid
  
  # produce plot
  plot1 <- ggplot(df) + theme_bw() + theme(panel.grid.minor.x = element_blank(),
                                           panel.grid.major.x = element_blank())
  plot1 <- plot1 + geom_vline(aes(xintercept = breaks_vec), col = grey(0.9))
  plot1 <- plot1 + scale_y_continuous(limits = c(0,1), expand = c(0,0))
  plot1 <- plot1 + geom_point(aes(x = x_vec, y = mc_accept, color = col))
  plot1 <- plot1 + xlab(x_lab) + ylab("coupling acceptance rate")
  plot1 <- plot1 + scale_colour_gradientn(colours = c("red", "blue"), name = "thermodynamic\npower", limits = c(0,1))
  
  return(plot1)
}

#------------------------------------------------
#' @title Plot spline quantiles
#'
#' @description Given age-dependent quantile information, as produced by
#'   \code{get_spline_quantiles()}, produce a ribbon plot over ages. Optionally
#'   overlay
#'
#' @param df_spline_quantile dataframe of quantiles broken down by age, as
#'   produced by \code{get_spline_quantiles()}.
#' @param df_data_quantile dataframe of exact data quantiles, as produced by
#'   \code{get_data_quantiles_p()} or \code{get_data_quantiles_m()}. Ignored if
#'   \code{NULL}.
#' @param title plot title.
#' @param ylim y-limits of the plot.
#'
#' @export

# function for plotting spline quantiles
plot_spline_quantiles <- function(df_spline_quantile,
                                  df_data_quantile = NULL,
                                  title = "",
                                  ylim = c(0,1)) {
  
  # avoid no visible binding error
  age <- Q2.5 <- Q97.5 <- Q50 <- lower <- upper <- NULL
    
  # check inputs
  assert_dataframe(df_spline_quantile)
  assert_in(c("age", "Q2.5", "Q50", "Q97.5"), names(df_spline_quantile))
  if (!is.null(df_data_quantile)) {
    assert_dataframe(df_data_quantile)
    assert_in(c("age", "lower", "upper", "mean"), names(df_data_quantile))
  }
  assert_single_string(title)
  assert_limit(ylim)
  
  # produce plotting object
  ret <- ggplot(df_spline_quantile) + theme_bw() + 
    geom_ribbon(aes(x = age, ymin = Q2.5, ymax = Q97.5), fill = "blue", alpha = 0.5) +
    geom_line(aes(x = age, y = Q50)) +
    ggplot2::ggtitle(title) + ggplot2::ylim(ylim)
  
  # add raw data confidence intervals
  if (!is.null(df_data_quantile)) {
    ret <- ret + geom_pointrange(aes(x = age, y = mean, ymin = lower, ymax = upper),
                                 size = 0.2, data = df_data_quantile)
  }
  
  return(ret)
}
