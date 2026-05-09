#' Validate custom immune gene sets against MSigDB Hallmarks
#'
#' @description
#' Computes ssGSEA scores using MSigDB Hallmark gene sets and correlates
#' them with the custom immune gene set scores from an `irps` object.
#' This provides external validation that the custom gene sets capture
#' biologically meaningful immune signals.
#'
#' @param irps_obj An `irps` object returned by \code{\link{compute_irps}}.
#' @param expr_mat A numeric expression matrix (features × samples). Must
#'   be the same matrix used for \code{\link{compute_irps}}.
#' @param gene_map A data.frame mapping gene symbols to feature IDs.
#' @param gene_col,feature_col Column names in `gene_map` for gene symbols
#'   and feature IDs.
#' @param species Character, species identifier for msigdbr. Default: "Homo sapiens".
#' @param category Character, MSigDB category. Default: "H" (Hallmark).
#' @param immune_hallmarks Character vector of Hallmark gene set names to
#'   use. If `NULL`, a default set of 12 immune-related Hallmarks is used.
#' @param min_set_size Integer, minimum number of mapped features to retain
#'   a Hallmark set. Default: 3.
#'
#' @return A list with components:
#' \itemize{
#'   \item `correlation_matrix`: Spearman correlation matrix (custom × Hallmark).
#'   \item `hallmark_scores`: Matrix of Hallmark ssGSEA scores.
#'   \item `hallmark_stats`: Data frame of Hallmark scores by group (if Subtype
#'     is available).
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' result <- compute_irps(expr, gene_map, platform = "proteomics")
#' val <- validate_hallmark(result, expr, gene_map)
#'
#' # View strongest correlations
#' cor_mat <- val$correlation_matrix
#' which(cor_mat == max(cor_mat, na.rm = TRUE), arr.ind = TRUE)
#' }
validate_hallmark <- function(irps_obj,
                               expr_mat,
                               gene_map,
                               gene_col = "Gene",
                               feature_col = "FeatureID",
                               species = "Homo sapiens",
                               category = "H",
                               immune_hallmarks = NULL,
                               min_set_size = 3L) {
  if (!inherits(irps_obj, "irps")) {
    stop("'irps_obj' must be an irps object from compute_irps()")
  }
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop("Package 'msigdbr' is required for Hallmark validation. ",
         "Install with: install.packages('msigdbr')")
  }

  # ---- Default immune-related Hallmarks ----
  if (is.null(immune_hallmarks)) {
    immune_hallmarks <- c(
      "HALLMARK_INTERFERON_GAMMA_RESPONSE",
      "HALLMARK_INTERFERON_ALPHA_RESPONSE",
      "HALLMARK_INFLAMMATORY_RESPONSE",
      "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
      "HALLMARK_IL2_STAT5_SIGNALING",
      "HALLMARK_IL6_JAK_STAT3_SIGNALING",
      "HALLMARK_ALLOGRAFT_REJECTION",
      "HALLMARK_COMPLEMENT",
      "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
      "HALLMARK_APOPTOSIS",
      "HALLMARK_KRAS_SIGNALING_UP",
      "HALLMARK_P53_PATHWAY"
    )
  }

  # ---- Load MSigDB ----
  message("Loading MSigDB Hallmark gene sets from msigdbr...")
  hallmarks <- msigdbr::msigdbr(species = species, category = category)
  hallmark_sets <- split(hallmarks$gene_symbol, hallmarks$gs_name)

  available <- intersect(immune_hallmarks, names(hallmark_sets))
  if (length(available) == 0) {
    stop("None of the requested Hallmark gene sets found in msigdbr")
  }
  message(sprintf("  %d/%d Hallmarks available in msigdbr",
                  length(available), length(immune_hallmarks)))

  h_sets <- lapply(available, function(nm) hallmark_sets[[nm]])
  names(h_sets) <- gsub("HALLMARK_", "", available)

  # ---- Map Hallmark genes to features ----
  h_mapped <- map_genes_to_features(
    h_sets, gene_map, gene_col, feature_col, min_set_size
  )

  all_fids <- rownames(expr_mat)
  h_mapped <- lapply(h_mapped, function(fids) fids[fids %in% all_fids])
  h_mapped <- h_mapped[lengths(h_mapped) >= min_set_size]

  message(sprintf("  %d Hallmark sets mappable to expression features",
                  length(h_mapped)))

  if (length(h_mapped) < 2) {
    stop("Too few Hallmark sets with sufficient feature mapping")
  }

  # ---- Run Hallmark ssGSEA ----
  message("Running Hallmark ssGSEA...")
  h_gsc <- lapply(names(h_mapped), function(nm) {
    GSEABase::GeneSet(h_mapped[[nm]], setName = nm)
  })
  h_gsc <- GSEABase::GeneSetCollection(h_gsc)

  param <- GSVA::ssgseaParam(
    as.matrix(expr_mat), h_gsc,
    minSize = min_set_size, maxSize = 500, normalize = TRUE
  )
  hallmark_scores <- GSVA::gsva(param)

  # ---- Correlation: custom × Hallmark ----
  common_samples <- intersect(colnames(irps_obj$ssgsea_scores),
                              colnames(hallmark_scores))
  message(sprintf("  %d common samples for correlation", length(common_samples)))

  cor_mat <- stats::cor(
    t(irps_obj$ssgsea_scores[, common_samples]),
    t(hallmark_scores[, common_samples]),
    method = "spearman", use = "pairwise.complete.obs"
  )

  # ---- Hallmark scores by group ----
  hallmark_stats <- NULL
  if ("Subtype" %in% colnames(irps_obj$irps_scores)) {
    hallmark_long <- as.data.frame(hallmark_scores) %>%
      tibble::rownames_to_column("Signature") %>%
      tidyr::pivot_longer(-"Signature", names_to = "SampleID", values_to = "Score") %>%
      dplyr::left_join(
        irps_obj$irps_scores[, c("SampleID", "Subtype")],
        by = "SampleID"
      ) %>%
      dplyr::filter(!is.na(.data$Subtype))

    hallmark_stats <- hallmark_long %>%
      dplyr::group_by(.data$Signature, .data$Subtype) %>%
      dplyr::summarise(
        Mean = mean(.data$Score, na.rm = TRUE),
        SD   = sd(.data$Score, na.rm = TRUE),
        .groups = "drop"
      )
  }

  list(
    correlation_matrix = cor_mat,
    hallmark_scores    = hallmark_scores,
    hallmark_stats     = hallmark_stats
  )
}
