### Plot CIBERSORT QC ----------------------------------------------------------
#'
#' @author Andrea Del Cortona <andrea.delcortona@gmail.com>
#' @note print CIBERSORT QC
#' @note 2024-12-09
#'
#' @param infiltr_out infiltr_out output object
#' @param sample_groups The columns of the metadata matrix to be used for the
#'   comparisons between sample groups. E.g.: a column indicating which sample
#'   was treated and which one not.
#' @param my_palette A vector of color to be passed to the ggplot functions. If
#'   provided by the user, it must be the same length of number of factors in
#'   sample_groups.
#'   Default: "default", it uses standard ggplot2 palette.
#' @param save_plots save plots in pdf and png format
#'   Default: TRUE
#' @param outdir Output directory of the plots
#'   Default: "default", it prints in the current working directory
plot_cb_qc = function(
    infiltr_out,
    my_palette,
    save_plots = TRUE,
    outdir = "default"
){

  #### General ####

  # get number of metadata groups
  n_fact = metadata[[sample_groups]] %>%
    as.factor() %>%
    levels() %>%
    length()

  # check my palette
  if(length(my_palette) == 1){
    if(my_palette == "default"){
      my_palette = scales::hue_pal()(n_fact)
    }
  }

  # prep dataframes for printing figures
  message("[", as.POSIXct(lubridate::now()), "] ... Prepping dataframes for figures and report")
  plot_tables = get_plot_tables(
    infiltr_out
  )

  mcp_relative_infil = plot_tables[["mcp_relative_infil"]]
  cb_all = plot_tables[["cb_all"]]
  cb_sign = plot_tables[["cb_sign"]]
  quantiseq = plot_tables[["quantiseq"]]

  mcp_relative_infil$model = "MCP-counter"
  cb_all$model = "CIBERSORT all"
  cb_sign$model = "CIBERSORT sign."
  quantiseq$model = "quanTIseq"


  #### CIBERSORT QC ####

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

  cb_qll_vs_sign = ggplot(cb_combined, aes(x = cell_type, y = value, fill = model)) +
    geom_boxplot(outlier.shape = 1, outlier.size = 0.5) +
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
    #pp_value,
    cb_qll_vs_sign,
    nrow = 2,
    rel_heights = c(2, 5)
  )











  #### CIBERSORT QC ####


  quantiseq = infiltr_out[["quantiseq"]]

  # transfor and plot per cell abundances
  quantiseq$Sample = rownames(quantiseq)
  quantiseq$group = metadata[[sample_groups]]
  quantiseq = data.table::melt(quantiseq)
  colnames(quantiseq) = c("Sample", "group", "cell_type", "value")

  quantiseq$cell_type = factor(
    quantiseq$cell_type,
    levels = c(
      "B.cells", "Macrophages.M1", "Macrophages.M2", "Monocytes", "Neutrophils",
      "NK.cells", "T.cells.CD4", "T.cells.CD8", "Tregs", "Dendritic.cells", "Other"
    )
  )


}
