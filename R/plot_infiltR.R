### Plot infiltR figures -------------------------------------------------------
#'
#' @author Andrea Del Cortona <andrea.delcortona@gmail.com>
#' @note wrapper for plotting infiltR summary figures
#' @note 2024-11-13
#'
#' @param infiltr_out infiltR output object
#' @param save_plots save plots in pdf and png format
#'   Default: TRUE
#' @param outdir Output directory of the plots
#'   Default: "default", it prints in the current working directory
plot_infiltR <- function(
    infiltr_out,
    save_plots = TRUE,
    outdir = "default",
){


  mcp_absolute_infil = infiltr_out[["mcp_counter"]]

  # t-test
  tests = mcp_absolute_infil %>%
    rstatix::t_test(all_infiltrating ~ group)

  p_abs_all = ggplot(mcp_absolute_infil, aes(group, all_infiltrating, fill = group)) +
    geom_boxplot() +
    ggpubr::stat_compare_means(method = "anova", label.x = 0.75, label.y = 1) +
    scale_y_log10() +
    xlab("") +
    ylab("cell fractions") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  for(k in 1:(nrow(tests)-1)){
    p_abs_all = p_abs_all +
      ggplot2::annotate(
        "text",
        x = k+1,
        y = max(mcp_absolute_infil$all_infiltrating)*1.2,
        label = c(tests$p.adj.signif)[[k]]
      )
  }
  p_abs_all

}
