
#------------------------------------------------
# Plot Metropolis coupling acceptance rates
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
# Plot autocorrelation
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
# Plot parameter estimates
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
# Plot parameter correlation
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
# Plot 95\% credible intervals
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
# plot credible intervals
plot_CrI <- function (x, show = NULL, param_names = NULL, rotate = FALSE) {
  
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
    ggplot2::geom_pointrange(ggplot2::aes(x = param, y = Q50, ymin = Q2.5, ymax = Q97.5)) +
    ggplot2::geom_point(ggplot2::aes(x = param, y = MAP, col = "red")) +
    ggplot2::scale_x_discrete(labels = param_names) +
    ggplot2::xlab("") + ggplot2::ylab("95% CrI")
  
  # optionally rotate axis labels
  if (rotate) {
    ret <- ret + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  }
  
  return(ret)
}

#------------------------------------------------
# plot credible intervals coloured by a given variable
plot_CrI_colour <- function (x, col = 1, col_names = NULL,
                             show = NULL, param_names = NULL, rotate = FALSE) {
  #x <- df_summary
  
  # check inputs
  assert_dataframe(x)
  assert_in(c("param", "Q2.5", "Q50", "Q97.5"), names(x))
  
  # add colour column
  x$col <- as.factor(col)
  
  # subset parameters
  if (is.null(show)) {
    show <- x$param
  }
  df_plot <- subset(x, param %in% show)
  
  # default param names
  if (is.null(param_names)) {
    param_names <- unique(show)
  }
  
  # default colour names
  if (is.null(col_names)) {
    col_names <- as.character(col)
  }
  
  # ensure factor levels ordered correctly
  df_plot$param <- factor(df_plot$param, levels = show)
  
  # create plot
  ret <- ggplot2::ggplot(df_plot) + ggplot2::theme_bw() +
    ggplot2::geom_pointrange(ggplot2::aes(x = param, y = Q50, ymin = Q2.5, ymax = Q97.5, col = col),
                          position = ggplot2::position_dodge(width = 0.5)) + 
    ggplot2::scale_x_discrete(labels = param_names) +
    ggplot2::xlab("") + ggplot2::ylab("95% CrI") +
    ggplot2::scale_color_discrete(labels = col_names) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank()) +
    ggplot2::geom_vline(xintercept = 1:6 + 0.5, col = grey(0.75))
  
  # optionally rotate axis labels
  if (rotate) {
    ret <- ret + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  }
  
  return(ret)
}
