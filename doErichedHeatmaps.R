suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(magick)
})


doEnrichedHeatmaps <- function(
    mats_groups, titles_groups, out_dir,
    file_name = "two_groups_side_by_side.pdf",
    group_names = NULL,
    order_rows_by_rowmean = TRUE,
    palette = c("white","blue"),
    width_per_group = 2,
    height_per_panel = 2.0,
    meta_height_mm = 25,
    meta_gap_mm = 1,
    show_signal_legend = TRUE,
    legend_mode = c("capped","quantile"),  
    legend_q_hi = 0.95,                   
    legend_quantiles = c(0, .1, .25, .5, .75, .9, .95),  
    clip_high_q = NA_real_                 
){
  legend_mode <- match.arg(legend_mode)
  stopifnot(length(mats_groups) == length(titles_groups))
  G <- length(mats_groups)
  if (is.null(group_names)) group_names <- paste0("Group_", seq_len(G))
  stopifnot(length(group_names) == G)
  all_mats <- unlist(mats_groups, recursive = FALSE)
  stopifnot(all(vapply(all_mats, is.matrix, logical(1))))
  ncols <- vapply(all_mats, ncol, integer(1))
  stopifnot(length(unique(ncols)) == 1L)
  p <- ncols[1]
  if (isTRUE(order_rows_by_rowmean)) {
    mats_groups <- lapply(mats_groups, function(gl) lapply(gl, function(M)
      M[order(rowMeans(M, na.rm = TRUE), decreasing = TRUE), , drop = FALSE]))
    all_mats <- unlist(mats_groups, recursive = FALSE)
  }
  all_vals <- unlist(lapply(all_mats, as.numeric), use.names = FALSE)
  all_vals <- all_vals[is.finite(all_vals)]
  stopifnot(length(all_vals) > 0)
  vals_for_scale <- all_vals
  if (!is.na(clip_high_q)) {
    hi_clip <- as.numeric(quantile(all_vals, probs = clip_high_q, na.rm = TRUE))
    vals_for_scale[vals_for_scale > hi_clip] <- hi_clip
  }
  qs <- quantile(vals_for_scale,
                 probs = c(0,.1,.2,.3,.4,.5,.6,.7,.8,.9,.95,.98,.995,1),
                 na.rm = TRUE)
  col_fun <- colorRamp2(qs, colorRampPalette(palette)(length(qs)))
  make_ht <- function(M, title, top_anno = NULL){
    Heatmap(M,
            name = paste0(title, "_signal"),
            col = col_fun,
            cluster_rows = FALSE, cluster_columns = FALSE,
            show_row_names = FALSE, show_column_names = FALSE,
            row_title = title, row_title_gp = gpar(fontsize = 8, fontface = "bold"),
            row_title_rot = 90, use_raster = TRUE, raster_quality = 2,
            raster_device = "CairoPNG", show_heatmap_legend = FALSE,
            top_annotation = top_anno)
  }
  build_column <- function(mats, titles, gname){
    means_list <- lapply(mats, function(M) colMeans(M, na.rm = TRUE))
    stopifnot(all(vapply(means_list, length, integer(1)) == p))
    k <- length(mats)
    line_cols <- if (k <= 3) c("#E41A1C","#377EB8","#4DAF4A")[seq_len(k)] else grDevices::hcl.colors(k, "Dark 3")
    yr <- range(unlist(means_list), na.rm = TRUE)
    yr_diff <- diff(yr)
    top_meta <- HeatmapAnnotation(
      Top = function(index) {
        normalize01 <- function(v) if (yr_diff == 0) rep(0.5, length(v)) else (v - yr[1]) / yr_diff
      Top = function(index) {
        pushViewport(viewport())
        grid.text(gname, x = unit(0, "npc"), y = unit(1, "npc") - unit(1, "mm"),
                  just = c("left","top"), gp = gpar(fontsize = 10, fontface = "bold"))
        pushViewport(viewport(y = unit(5, "mm"),
                              height = unit(1, "npc") - unit(16, "mm"),
                              xscale = c(0,1), yscale = c(0,1)))
        x <- seq_along(index) / length(index)
        for (i in seq_along(means_list)) {
          grid.lines(x = x, y = normalize01(means_list[[i]][index]),
                     gp = gpar(col = line_cols[i], lwd = 1.2), default.units = "native")
        }
        grid.yaxis(gp = gpar(cex = 0.6))
        popViewport(); popViewport()
      },
      height = unit(meta_height_mm, "mm") + unit(5, "mm"),
      show_annotation_name = FALSE
    )
    first_ht <- make_ht(mats[[1]], titles[1], top_anno = top_meta)
    col_ht <- first_ht
    if (length(mats) > 1) {
      for (i in 2:length(mats)) col_ht <- col_ht %v% make_ht(mats[[i]], titles[i])
    }
    list(
      ht = col_ht,
      line_cols = line_cols,
      line_labels = paste0(gname, ": ", titles)
    )
  }
  cols_info <- lapply(seq_len(G), function(g)
    build_column(mats_groups[[g]], titles_groups[[g]], group_names[g]))
  cols <- lapply(cols_info, `[[`, "ht")
  line_cols_all   <- unlist(lapply(cols_info, `[[`, "line_cols"), use.names = FALSE)
  line_labels_all <- unlist(lapply(cols_info, `[[`, "line_labels"), use.names = FALSE)
  if (legend_mode == "capped") {
    lo <- min(all_vals, na.rm = TRUE)
    hi <- as.numeric(quantile(all_vals, probs = legend_q_hi, na.rm = TRUE))
    at_vals <- pretty(c(lo, hi), n = 5)
    at_vals <- at_vals[at_vals >= lo & at_vals <= hi]
    sig_legend <- Legend(
      title = "Signal",
      col_fun = col_fun,
      at = at_vals,
      labels = format(at_vals, digits = 3, trim = TRUE),
      direction = "vertical",
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    )
  } else { 
    at_q <- as.numeric(quantile(all_vals, probs = legend_quantiles, na.rm = TRUE))
    sig_legend <- Legend(
      title = "Signal (quantiles)",
      col_fun = col_fun,
      at = at_q,
      labels = paste0(legend_quantiles * 100, "%"),
      direction = "vertical",
      title_gp = gpar(fontsize = 9, fontface = "bold"),
      labels_gp = gpar(fontsize = 8)
    )
  }
  meta_legend <- Legend(
    title = "Metaplot", at = line_labels_all, type = "lines",
    legend_gp = gpar(lwd = 2, col = line_cols_all),
    labels_gp = gpar(fontsize = 8)
  )
  if (show_signal_legend) {
    right_legends <- packLegend(sig_legend, meta_legend, direction = "vertical", gap = unit(4, "mm"))
  } else {
    right_legends <- packLegend(meta_legend, direction = "vertical", gap = unit(4, "mm"))
  }
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_file <- file.path(out_dir, file_name)
  grDevices::cairo_pdf(out_file,
                       width  = max(2.5, G * width_per_group) + 1.6,
                       height = meta_height_mm/25.4 + (max(vapply(mats_groups, length, integer(1))) * height_per_panel),
                       family = "Helvetica")
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = G + 1,
                                             widths = unit.c(rep(unit(1, "null"), G), unit(16, "lines")))))
  for (g in seq_len(G)) {
    pushViewport(viewport(layout.pos.row = 1, layout.pos.col = g))
    draw(cols[[g]], newpage = FALSE)
    popViewport()
  }
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = G + 1))
  draw(right_legends, x = unit(0.5, "npc"), y = unit(0.5, "npc"))
  popViewport(); popViewport()
  dev.off()
  message("Wrote: ", out_file)
  invisible(list(file = out_file, colors = col_fun))
}





