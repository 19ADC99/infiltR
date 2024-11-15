### Plot infiltR figures -------------------------------------------------------
#'
#' @author Andrea Del Cortona <andrea.delcortona@gmail.com>
#' @note wrapper for plotting infiltR summary figures
#' @note 2024-11-13
#'
#' @param infiltr_out infiltR output object
#' @param metadata A metadata matrix with samples as rows
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
plot_infiltR <- function(
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



  #### MCP-counter ####

  mcp_absolute_infil = infiltr_out[["mcp_counter"]]

  # t-test
  tests_abs = mcp_absolute_infil %>%
    rstatix::t_test(all_infiltrating ~ group)

  p_mcp_abs_all = ggplot(mcp_absolute_infil, aes(group, all_infiltrating, fill = group)) +
    geom_boxplot() +
    ggpubr::stat_compare_means(method = "anova", label.x = 0.75, label.y = 1) +
    scale_fill_manual(values = my_palette) +
    scale_y_log10() +
    xlab("") +
    ylab("cell fractions") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  for(k in 1:(nrow(tests_abs)-1)){
    p_mcp_abs_all = p_mcp_abs_all +
      ggplot2::annotate(
        "text",
        x = k+1,
        y = max(mcp_absolute_infil$all_infiltrating)*1.2,
        label = c(tests_abs$p.adj.signif)[[k]]
      )
  }


  # transfor and plot per cell abundances
  mcp_absolute_infil$Sample = rownames(mcp_absolute_infil)
  mcp_absolute_infil$group = metadata[[sample_groups]]
  mcp_absolute_infil = mcp_absolute_infil %>%
    dplyr::select(-c("all_infiltrating"))
  mcp_absolute_infil = data.table::melt(mcp_absolute_infil)
  colnames(mcp_absolute_infil) = c("group", "Sample", "cell_type", "value")

  mcp_absolute_infil$cell_type = factor(
    mcp_absolute_infil$cell_type,
    levels = c(
      "T cells", "CD8 T cells", "Cytotoxic lymphocytes", "B lineage",
      "NK cells", "Monocytic lineage", "Myeloid dendritic cells",
      "Neutrophils", "Endothelial cells", "Fibroblasts"
    )
  )

  ## retrieve sign annotation and sign positions
  # t-test
  tests_abscell = mcp_absolute_infil %>%
    dplyr::group_by(cell_type) %>%
    rstatix::t_test(value ~ group)
  tests_abscell = tests_abscell %>% rstatix::add_xy_position(x = "cell_type")

  my_pvalues = get_pvalues_positions(mcp_absolute_infil, tests_abscell)
  my_y_positions = my_pvalues[["my_y_positions"]]
  my_x_min = my_pvalues[["my_x_min"]]
  my_x_max = my_pvalues[["my_x_max"]]
  my_annotation = my_pvalues[["my_annotation"]]

  p_mcp_abs_cells = ggplot(mcp_absolute_infil, aes(cell_type, value, fill = group)) +
    geom_boxplot(outlier.shape = 1, outlier.size = 0.5) +
    scale_fill_manual(values = my_palette) +
    scale_y_log10() +
    xlab("") +
    ylab("cell fractions") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  for(k in 1:length(my_annotation)){
    p_mcp_abs_cells = p_mcp_abs_cells +
      ggplot2::annotate(
        "text",
        x = mean(c(my_x_min[[k]], my_x_max[[k]])),
        y = my_y_positions[[k]]*1.05,
        label = my_annotation[[k]]
      ) +
      ggplot2::annotate(
        "segment",
        x = my_x_min[[k]],
        xend = my_x_max[[k]],
        y = my_y_positions[[k]],
        yend = my_y_positions[[k]],
      )
  }
  p_mcp_abs_cells


  # # transfor and plot per relative cell abundances
  # mcp_relative_infil = infiltr_out[["mcp_counter_norm"]]
  # mcp_relative_infil$Sample = rownames(mcp_relative_infil)
  # mcp_relative_infil$group = metadata[[sample_groups]]
  # mcp_relative_infil = data.table::melt(mcp_relative_infil)
  # colnames(mcp_relative_infil) = c("group", "Sample", "cell_type", "value")
  #
  # mcp_relative_infil$cell_type = factor(
  #   mcp_relative_infil$cell_type,
  #   levels = c(
  #     "T cells", "CD8 T cells", "Cytotoxic lymphocytes", "B lineage",
  #     "NK cells", "Monocytic lineage", "Myeloid dendritic cells",
  #     "Neutrophils", "Endothelial cells", "Fibroblasts"
  #   )
  # )
  #
  # # t-test
  # tests_relcell = mcp_relative_infil %>%
  #   dplyr::group_by(cell_type) %>%
  #   rstatix::t_test(value ~ group)
  # tests_relcell = tests_relcell %>% rstatix::add_xy_position(x = "cell_type")
  #
  # my_pvalues = get_pvalues_positions(mcp_relative_infil, tests_relcell)
  # my_y_positions = my_pvalues[["my_y_positions"]]
  # my_x_min = my_pvalues[["my_x_min"]]
  # my_x_max = my_pvalues[["my_x_max"]]
  # my_annotation = my_pvalues[["my_annotation"]]
  #
  # p_mcp_rel_cells = ggplot(mcp_relative_infil, aes(cell_type, value, fill = group)) +
  #   geom_boxplot(outlier.shape = 1, outlier.size = 0.5) +
  #   scale_fill_manual(values = my_palette) +
  #   scale_y_log10() +
  #   xlab("") +
  #   ylab("cell fractions") +
  #   theme_bw() +
  #   theme(axis.text.x = element_text(angle = 45, hjust = 1))
  #
  # for(k in 1:length(my_annotation)){
  #   p_mcp_rel_cells = p_mcp_rel_cells +
  #     ggplot2::annotate(
  #       "text",
  #       x = mean(c(my_x_min[[k]], my_x_max[[k]])),
  #       y = my_y_positions[[k]]*1.05,
  #       label = my_annotation[[k]]
  #     ) +
  #     ggplot2::annotate(
  #       "segment",
  #       x = my_x_min[[k]],
  #       xend = my_x_max[[k]],
  #       y = my_y_positions[[k]],
  #       yend = my_y_positions[[k]],
  #     )
  # }




}
