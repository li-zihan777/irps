#' irps S3 class and methods
#'
#' @description
#' The `irps` S3 object stores all results from an IRPS computation run,
#' including ssGSEA scores, IRPS scores, and metadata.
#'
#' @details
#' An `irps` object is a list with the following components:
#' \itemize{
#'   \item `ssgsea_scores`: matrix of ssGSEA enrichment scores (gene sets × samples)
#'   \item `irps_scores`: data.frame with per-sample IRPS, Favorable, Unfavorable scores
#'   \item `gene_sets`: list of gene sets used (named character vectors of feature IDs)
#'   \item `favorable_sets`: character vector of favorable gene set names
#'   \item `unfavorable_sets`: character vector of unfavorable gene set names
#'   \item `platform`: character, either "proteomics" or "transcriptomics"
#'   \item `n_features_used`: integer, total number of features detected in expression matrix
#'   \item `params`: list of parameters used for the computation
#' }
#'
#' @name irps-class
NULL

#' Create an irps object
#' @param ssgsea_scores Matrix of ssGSEA scores (gene sets × samples).
#' @param irps_scores Data frame of per-sample IRPS scores.
#' @param gene_sets Named list of gene sets used.
#' @param favorable_sets Names of favorable gene sets.
#' @param unfavorable_sets Names of unfavorable gene sets.
#' @param platform Character string: "proteomics" or "transcriptomics".
#' @param n_features_used Integer, number of features detected.
#' @param params List of additional parameters.
#' @return An S3 object of class `irps`.
#' @keywords internal
new_irps <- function(ssgsea_scores, irps_scores, gene_sets,
                     favorable_sets, unfavorable_sets,
                     platform, n_features_used, params = list()) {
  obj <- list(
    ssgsea_scores     = ssgsea_scores,
    irps_scores       = irps_scores,
    gene_sets         = gene_sets,
    favorable_sets     = favorable_sets,
    unfavorable_sets   = unfavorable_sets,
    platform          = platform,
    n_features_used   = n_features_used,
    params            = params
  )
  class(obj) <- "irps"
  obj
}

#' @export
print.irps <- function(x, ...) {
  cat("IRPS Object\n")
  cat(strrep("-", 50), "\n")
  cat(sprintf("Platform:          %s\n", x$platform))
  cat(sprintf("Features used:     %d\n", x$n_features_used))
  cat(sprintf("Gene sets:         %d (total)\n", length(x$gene_sets)))
  cat(sprintf("  Favorable:       %d\n", length(x$favorable_sets)))
  cat(sprintf("  Unfavorable:     %d\n", length(x$unfavorable_sets)))
  cat(sprintf("Samples:           %d\n", nrow(x$irps_scores)))
  cat(strrep("-", 50), "\n")

  if (!is.null(x$irps_scores$Subtype)) {
    by_group <- x$irps_scores %>%
      group_by(Subtype) %>%
      summarise(
        n = n(),
        mean_IRPS = mean(.data$IRPS, na.rm = TRUE),
        sd_IRPS = sd(.data$IRPS, na.rm = TRUE),
        .groups = "drop"
      )
    cat("IRPS by group:\n")
    for (i in seq_len(nrow(by_group))) {
      cat(sprintf("  %s (n=%d): %.4f +/- %.4f\n",
                  by_group$Subtype[i], by_group$n[i],
                  by_group$mean_IRPS[i], by_group$sd_IRPS[i]))
    }
  } else {
    cat("IRPS summary:\n")
    cat(sprintf("  Mean:    %.4f\n", mean(x$irps_scores$IRPS, na.rm = TRUE)))
    cat(sprintf("  SD:      %.4f\n", sd(x$irps_scores$IRPS, na.rm = TRUE)))
    cat(sprintf("  Median:  %.4f\n", stats::median(x$irps_scores$IRPS, na.rm = TRUE)))
    cat(sprintf("  Range:   [%.4f, %.4f]\n",
                min(x$irps_scores$IRPS, na.rm = TRUE),
                max(x$irps_scores$IRPS, na.rm = TRUE)))
  }
  invisible(x)
}

#' @export
summary.irps <- function(object, ...) {
  irps_vec <- object$irps_scores$IRPS
  list(
    platform        = object$platform,
    n_features      = object$n_features_used,
    n_gene_sets     = length(object$gene_sets),
    n_favorable     = length(object$favorable_sets),
    n_unfavorable   = length(object$unfavorable_sets),
    n_samples       = length(irps_vec),
    irps_mean       = mean(irps_vec, na.rm = TRUE),
    irps_sd         = sd(irps_vec, na.rm = TRUE),
    irps_median     = stats::median(irps_vec, na.rm = TRUE),
    irps_min        = min(irps_vec, na.rm = TRUE),
    irps_max        = max(irps_vec, na.rm = TRUE),
    gene_set_sizes  = vapply(object$gene_sets, length, integer(1))
  )
}

#' Plot an irps object
#'
#' @description
#' S3 plot method for `irps` objects. Generates publication-quality
#' visualizations by dispatching to the appropriate plotting function.
#'
#' @param x An `irps` object.
#' @param y Ignored (required by generic).
#' @param type Character, type of plot:
#'   \describe{
#'     \item{"boxplot"}{IRPS score boxplot by group (default).}
#'     \item{"heatmap"}{ssGSEA score heatmap with sample annotations.}
#'     \item{"radar_favorable"}{Radar chart of 5 favorable gene sets.}
#'     \item{"radar_all"}{Radar chart of all 8 immunotherapy-relevant sets.}
#'     \item{"portrait"}{6-dimension immune portrait faceted boxplot panel.}
#'   }
#' @param ... Additional arguments passed to the specific plot function.
#'
#' @return A ggplot object or (for radar charts) the data matrix invisibly.
#' @export
#' @method plot irps
plot.irps <- function(x, y, type = c("boxplot", "heatmap",
  "radar_favorable", "radar_all", "portrait"), ...) {
  type <- match.arg(type)
  switch(type,
    boxplot         = plot_irps_boxplot(x, ...),
    heatmap         = plot_irps_heatmap(x, ...),
    radar_favorable = plot_irps_radar(x, which = "favorable", ...),
    radar_all       = plot_irps_radar(x, which = "all", ...),
    portrait        = {
      port <- decompose_portrait(x)
      plot_portrait_panel(port, ...)
    }
  )
}
