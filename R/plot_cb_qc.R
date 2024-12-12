### Plot CIBERSORT QC ----------------------------------------------------------
#'
#' @author Andrea Del Cortona <andrea.delcortona@gmail.com>
#' @note print CIBERSORT QC
#' @note 2024-12-09
#' @title plot_cb_qc
#' @details
#'
#' This function plots QC statistics for CIBERSORT.
#'
#' @param infiltr_out infiltr_out output object
#' @param sample_groups The columns of the metadata matrix to be used for the
#'   comparisons between sample groups. E.g.: a column indicating which sample
#'   was treated and which one not.
#' @param metadata A metadata matrix with samples as rows
#' @param save_plots save plots in pdf and png format
#'   Default: TRUE
#' @param outdir Output directory of the plots
#'   Default: "default", it prints in the current working directory
#' @export
plot_cb_qc = function(
    infiltr_out,
    metadata,
    sample_groups,
    save_plots = TRUE,
    outdir = "default"
){

  # prep dataframes for printing figures
  message("[", as.POSIXct(lubridate::now()), "] ... Prepping dataframes for figures and report")
  plot_tables = get_plot_tables(
    infiltr_out,
    metadata,
    sample_groups
  )

  mcp_relative_infil = plot_tables[["mcp_relative_infil"]]
  cb_all = plot_tables[["cb_all"]]
  cb_sign = plot_tables[["cb_sign"]]
  quantiseq = plot_tables[["quantiseq"]]

  mcp_relative_infil$model = "MCP-counter"
  cb_all$model = "CIBERSORT all"
  cb_sign$model = "CIBERSORT sign."
  quantiseq$model = "quanTIseq"


  # CIBERSORT p-values distributions
  pp_value = ggplot(infiltr_out[["cb_all"]], aes(x = `P-value`)) +
    geom_histogram(aes(y = after_stat(density)), colour = "black", fill = "white", bins = 100) +
    geom_density(alpha = 0.1, color = "dodgerblue") +
    geom_vline(aes(xintercept = 0.05), color = "firebrick", linetype = "dashed", linewidth = 0.5) +
    labs("CIBERSORT p-values distribution") +
    theme_bw() +
    theme(panel.grid.minor = element_blank())


  # CIBERSORT samples number before and after filtering
  cb_counts = cbind(
    table(cb_all$group)/22,
    table(cb_sign$group)/22
  ) %>% as.data.frame()
  cb_counts$group = rownames(cb_counts)
  colnames(cb_counts) = c("CIBERSORT all", "CIBERSORT sign.", "group")
  cb_counts = data.table::melt(cb_counts)

  p_cb_counts = ggplot(cb_counts, aes(x = group, y = value, fill = variable)) +
    geom_bar(stat = "identity",
             colour = "black",
             position = position_dodge(),
             width = 0.5) +
    geom_text(aes(label = paste0("n=", value)),
              position = position_dodge(width = 0.5),
              color = "black",
              vjust = -0.25,
              size = 3.5) +
    scale_fill_manual(values = c("grey75", "orange")) +
    labs(x = "", y = "sample counts", fill = "") +
    theme_bw() +
    theme(legend.position = "none",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.background = element_rect(fill = "#d6eaf8"),
          strip.text = element_text(size = 10))


  # CIBERSORT all VS filtered
  cb_combined = rbind(cb_all, cb_sign)

  cb_all_vs_sign = ggplot(cb_combined, aes(x = cell_type, y = value, fill = model)) +
    geom_violin(alpha = 0.5, scale = "width", trim = TRUE) +
    stat_summary(fun = "median", colour = "red", geom = "crossbar", size = 0.15,
                 position = position_dodge(0.9), width = 0.5) +
    facet_wrap(~ group, nrow = 3, strip.position = "right") +
    scale_fill_manual(values = c("grey75", "orange")) +
    scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 10)) +
    labs(x = "", y = "cell fractions", fill = "") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          legend.position = "bottom",
          panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          strip.background = element_rect(fill = "#d6eaf8"),
          strip.text = element_text(size = 10))


  # plot composite output figure
  cb_composite = cowplot::plot_grid(
    cowplot::plot_grid(
      p_cb_counts, pp_value,
      nrow = 1,
      align = "h",
      rel_widths = c(0.5, 0.5)
    ),
    cb_all_vs_sign,
    nrow = 2,
    rel_heights = c(2, 5)
  )

  plot(cb_composite)

}
