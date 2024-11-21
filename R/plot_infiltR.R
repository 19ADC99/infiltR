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
plot_infiltR = function(
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
  tests_abs = tests_abs %>% rstatix::add_xy_position()

  my_pvalues = get_pvalues_positions(mcp_absolute_infil$all_infiltrating, tests_abs)
  my_y_positions = my_pvalues[["my_y_positions"]]
  my_x_min = my_pvalues[["my_x_min"]]
  my_x_max = my_pvalues[["my_x_max"]]
  my_annotation = my_pvalues[["my_annotation"]]

  p_mcp_abs_all = ggplot(mcp_absolute_infil, aes(group, all_infiltrating, fill = group)) +
    geom_boxplot(outlier.shape = 1, outlier.size = 0.5) +
    scale_fill_manual(values = my_palette) +
    scale_y_log10() +
    xlab("") +
    ylab("cell fractions") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  for(k in 1:(nrow(tests_abs))){
    p_mcp_abs_all = p_mcp_abs_all +
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
  p_mcp_abs_all

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


  # transfor and plot per relative cell abundances
  mcp_relative_infil = infiltr_out[["mcp_counter_norm"]]
  mcp_relative_infil$Sample = rownames(mcp_relative_infil)
  mcp_relative_infil$group = metadata[[sample_groups]]
  mcp_relative_infil = data.table::melt(mcp_relative_infil)
  colnames(mcp_relative_infil) = c("group", "Sample", "cell_type", "value")

  mcp_relative_infil$cell_type = factor(
    mcp_relative_infil$cell_type,
    levels = c(
      "T cells", "CD8 T cells", "Cytotoxic lymphocytes", "B lineage",
      "NK cells", "Monocytic lineage", "Myeloid dendritic cells",
      "Neutrophils", "Endothelial cells", "Fibroblasts"
    )
  )

  # t-test
  tests_relcell = mcp_relative_infil %>%
    dplyr::group_by(cell_type) %>%
    rstatix::t_test(value ~ group)
  tests_relcell = tests_relcell %>% rstatix::add_xy_position(x = "cell_type")

  my_pvalues = get_pvalues_positions(mcp_relative_infil, tests_relcell)
  my_y_positions = my_pvalues[["my_y_positions"]]
  my_x_min = my_pvalues[["my_x_min"]]
  my_x_max = my_pvalues[["my_x_max"]]
  my_annotation = my_pvalues[["my_annotation"]]

  p_mcp_rel_cells = ggplot(mcp_relative_infil, aes(cell_type, value, fill = group)) +
    geom_boxplot(outlier.shape = 1, outlier.size = 0.5) +
    scale_fill_manual(values = my_palette) +
    #scale_y_log10() +
    xlab("") +
    ylab("cell fractions") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  for(k in 1:length(my_annotation)){
    p_mcp_rel_cells = p_mcp_rel_cells +
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
  p_mcp_rel_cells



  #### CIBERSORT ####

  cb_all = infiltr_out[["cb_all"]]

  pp_value = ggplot(cb_all, aes(x = `P-value`)) +
    geom_histogram(aes(y = after_stat(density)), colour = "black", fill = "white", bins = 100) +
    geom_density(alpha = 0.1, color = "dodgerblue") +
    geom_vline(aes(xintercept = 0.05), color = "firebrick", linetype = "dashed", linewidth = 0.5) +
    theme_bw() +
    theme(panel.grid.minor = element_blank())
  pp_value

  # transfor and plot per cell abundances
  cb_all$Sample = rownames(cb_all)
  cb_all$group = metadata[[sample_groups]]
  cb_all = cb_all %>%
    dplyr::select(-c("P-value", "Correlation", "RMSE"))
  cb_all = data.table::melt(cb_all)
  colnames(cb_all) = c("Sample", "group", "cell_type", "value")

  cb_all$cell_type = factor(
    cb_all$cell_type,
    levels = c(
      "B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
      "T cells CD4 naive", "T cells CD4 memory resting", "T cells CD4 memory activated",
      "T cells follicular helper", "T cells regulatory (Tregs)", "T cells gamma delta",
      "NK cells resting", "NK cells activated", "Monocytes", "Macrophages M0",
      "Macrophages M1", "Macrophages M2", "Dendritic cells resting",
      "Dendritic cells activated", "Mast cells resting", "Mast cells activated",
      "Eosinophils", "Neutrophils"
    )
  )

  ## retrieve sign annotation and sign positions
  # t-test
  tests_cb_all = cb_all %>%
    dplyr::group_by(cell_type) %>%
    rstatix::t_test(value ~ group)
  tests_cb_all = tests_cb_all %>% rstatix::add_xy_position(x = "cell_type")

  my_pvalues = get_pvalues_positions(cb_all, tests_cb_all)
  my_y_positions = my_pvalues[["my_y_positions"]]
  my_x_min = my_pvalues[["my_x_min"]]
  my_x_max = my_pvalues[["my_x_max"]]
  my_annotation = my_pvalues[["my_annotation"]]

  p_cb_all = ggplot(cb_all, aes(cell_type, value, fill = group)) +
    geom_boxplot(outlier.shape = 1, outlier.size = 0.5) +
    scale_fill_manual(values = my_palette) +
    coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    xlab("") +
    ylab("cell fractions") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  for(k in 1:length(my_annotation)){
    p_cb_all = p_cb_all +
      ggplot2::annotate(
        "text",
        x = mean(c(my_x_min[[k]], my_x_max[[k]])),
        y = my_y_positions[[k]] + 0.0025,
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
  p_cb_all



  cb_sign = infiltr_out[["cb_sign"]]

  # transfor and plot per cell abundances
  cb_sign$Sample = rownames(cb_sign)
  cb_sign$group = metadata %>%
    dplyr::filter(sample %in% rownames(cb_sign)) %>%
    dplyr::select(dplyr::all_of(sample_groups)) %>%
    BiocGenerics::unlist()
  cb_sign = cb_sign %>%
    dplyr::select(-c("P-value", "Correlation", "RMSE"))
  cb_sign = data.table::melt(cb_sign)
  colnames(cb_sign) = c("Sample", "group", "cell_type", "value")

  cb_sign$cell_type = factor(
    cb_sign$cell_type,
    levels = c(
      "B cells naive", "B cells memory", "Plasma cells", "T cells CD8",
      "T cells CD4 naive", "T cells CD4 memory resting", "T cells CD4 memory activated",
      "T cells follicular helper", "T cells regulatory (Tregs)", "T cells gamma delta",
      "NK cells resting", "NK cells activated", "Monocytes", "Macrophages M0",
      "Macrophages M1", "Macrophages M2", "Dendritic cells resting",
      "Dendritic cells activated", "Mast cells resting", "Mast cells activated",
      "Eosinophils", "Neutrophils"
    )
  )

  ## retrieve sign annotation and sign positions
  # t-test
  tests_cb_sign = cb_sign %>%
    dplyr::group_by(cell_type) %>%
    rstatix::t_test(value ~ group)
  tests_cb_sign = tests_cb_sign %>% rstatix::add_xy_position(x = "cell_type")

  my_pvalues = get_pvalues_positions(cb_sign, tests_cb_sign)
  my_y_positions = my_pvalues[["my_y_positions"]]
  my_x_min = my_pvalues[["my_x_min"]]
  my_x_max = my_pvalues[["my_x_max"]]
  my_annotation = my_pvalues[["my_annotation"]]

  p_cb_sign = ggplot(cb_sign, aes(cell_type, value, fill = group)) +
    geom_boxplot(outlier.shape = 1, outlier.size = 0.5) +
    scale_fill_manual(values = my_palette) +
    #coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    xlab("") +
    ylab("cell fractions") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  for(k in 1:length(my_annotation)){
    p_cb_sign = p_cb_sign +
      ggplot2::annotate(
        "text",
        x = mean(c(my_x_min[[k]], my_x_max[[k]])),
        y = my_y_positions[[k]] + 0.0025,
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
  p_cb_sign


  #### quanTIseq ####

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

  ## retrieve sign annotation and sign positions
  # t-test
  tests_quanti = quantiseq %>%
    dplyr::group_by(cell_type) %>%
    rstatix::t_test(value ~ group)
  tests_quanti = tests_quanti %>% rstatix::add_xy_position(x = "cell_type")

  my_pvalues = get_pvalues_positions(quantiseq, tests_quanti)
  my_y_positions = my_pvalues[["my_y_positions"]]
  my_x_min = my_pvalues[["my_x_min"]]
  my_x_max = my_pvalues[["my_x_max"]]
  my_annotation = my_pvalues[["my_annotation"]]

  p_quanti = ggplot(quantiseq, aes(cell_type, value, fill = group)) +
    geom_boxplot(outlier.shape = 1, outlier.size = 0.5) +
    scale_fill_manual(values = my_palette) +
    #coord_cartesian(ylim = c(0, 1), expand = FALSE) +
    xlab("") +
    ylab("cell fractions") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  for(k in 1:length(my_annotation)){
    p_quanti = p_quanti +
      ggplot2::annotate(
        "text",
        x = mean(c(my_x_min[[k]], my_x_max[[k]])),
        y = my_y_positions[[k]] + 0.0025,
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
  p_quanti




}
