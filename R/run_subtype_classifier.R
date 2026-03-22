#' Classify new patients into CTS using reference data and user healthy controls
#'
#' Users must provide both new case samples and new healthy control samples.
#' Healthy controls are used during batch correction to preserve healthy-vs-sepsis
#' biology while adjusting for dataset-specific technical effects.
#'
#' @param new_case_data Matrix or data frame of new case samples; genes in rows,
#'   samples in columns.
#' @param new_healthy_data Matrix or data frame of new healthy control samples;
#'   genes in rows, samples in columns.
#' @param gene_list Optional character vector of classifier genes.
#' @param make_heatmap Logical; whether to draw a heatmap of predicted case samples.
#' @param ntrees Number of trees for random forest.
#' @param save_plots Logical; whether to save heatmap and silhouette plots as PDF files.
#' @param output_dir Directory in which PDF files should be saved.
#' @param heatmap_file File name for the heatmap PDF.
#' @param silhouette_file File name for the silhouette PDF.
#'
#' @return A list containing predictions, corrected matrices, silhouette metrics,
#'   model objects, and output file paths.
#' @export
run_subtype_classifier <- function(new_case_data,
                                   new_healthy_data,
                                   gene_list = NULL,
                                   make_heatmap = TRUE,
                                   ntrees = 500,
                                   save_plots = FALSE,
                                   output_dir = ".",
                                   heatmap_file = "CTS_heatmap.pdf",
                                   silhouette_file = "CTS_silhouette.pdf") {

  data("exp_core_g", package = "ctsSubtypeR", envir = environment())
  data("core_samples", package = "ctsSubtypeR", envir = environment())

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }

  heatmap_path <- file.path(output_dir, heatmap_file)
  silhouette_path <- file.path(output_dir, silhouette_file)

  norm_res <- normalize_with_healthy_combat(
    exp_core_g = exp_core_g,
    core_samples = core_samples,
    new_case_data = new_case_data,
    new_healthy_data = new_healthy_data,
    gene_list = gene_list
  )

  exp_core_combat <- norm_res$corrected$reference
  exp_case_combat <- norm_res$corrected$new_case

  core_samples <- core_samples[colnames(exp_core_combat), , drop = FALSE]
  cts <- droplevels(factor(as.character(core_samples$CTS)))
  train_levels <- levels(cts)

  train_df <- as.data.frame(t(exp_core_combat))
  train_df$CTS <- cts

  rf <- randomForest::randomForest(
    CTS ~ .,
    data = train_df,
    ntree = as.integer(ntrees),
    proximity = TRUE
  )

  pred_case <- stats::predict(rf, newdata = as.data.frame(t(exp_case_combat)))
  pred_case <- factor(as.character(pred_case), levels = train_levels)

  pred_df <- data.frame(
    sample_id = colnames(exp_case_combat),
    CTS = pred_case
  )

  if (isTRUE(make_heatmap)) {
    ordered_pred_df <- pred_df[order(pred_df$CTS), , drop = FALSE]
    rownames(ordered_pred_df) <- ordered_pred_df$sample_id
    exp_case_ordered <- exp_case_combat[, ordered_pred_df$sample_id, drop = FALSE]

    ann_colors <- list(
      CTS = c(
        "1" = "royalblue",
        "2" = "#B2DF8A",
        "3" = "orange",
        "4" = "grey70"
      )
    )

    if (isTRUE(save_plots)) {
      grDevices::pdf(heatmap_path, width = 9, height = 11)
      on.exit(grDevices::dev.off(), add = TRUE)
    }

    pheatmap::pheatmap(
      exp_case_ordered,
      color = grDevices::colorRampPalette(
        rev(RColorBrewer::brewer.pal(11, "PuOr"))
      )(50),
      clustering_distance_rows = "correlation",
      cluster_rows = TRUE,
      cluster_cols = FALSE,
      show_rownames = TRUE,
      show_colnames = FALSE,
      annotation_col = ordered_pred_df[, "CTS", drop = FALSE],
      annotation_colors = ann_colors,
      scale = "row",
      treeheight_row = 0,
      border_color = NA
    )

    if (isTRUE(save_plots)) {
      grDevices::dev.off()
    }
  }

  sil <- NULL
  silhouette_saved <- FALSE

  keep <- !is.na(pred_case)
  if (sum(keep) > 2 && length(unique(pred_case[keep])) > 1) {
    rf_case <- randomForest::randomForest(
      x = as.data.frame(t(exp_case_combat[, keep, drop = FALSE])),
      y = droplevels(pred_case[keep]),
      ntree = as.integer(ntrees),
      proximity = TRUE
    )

    sil <- cluster::silhouette(
      as.integer(droplevels(pred_case[keep])),
      stats::as.dist(1 - rf_case$proximity)
    )

    if (isTRUE(save_plots)) {
      grDevices::pdf(silhouette_path, width = 8, height = 6)
      plot(sil,
           main = "Silhouette plot of predicted CTS classes",
           col = c("royalblue", "#B2DF8A", "orange", "grey70"))
      grDevices::dev.off()
      silhouette_saved <- TRUE
    }
  }

  list(
    predictions = pred_df,
    expression_corrected = norm_res$corrected,
    silhouette = sil,
    rf_model = rf,
    genes_used = norm_res$genes_used,
    gene_overlap = length(norm_res$genes_used),
    batch = norm_res$batch,
    bio_group = norm_res$bio_group,
    plot_files = list(
      heatmap = if (isTRUE(save_plots) && isTRUE(make_heatmap)) heatmap_path else NULL,
      silhouette = if (silhouette_saved) silhouette_path else NULL
    )
  )
}
