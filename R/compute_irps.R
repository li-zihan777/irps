#' Compute the Immune Response Potential Score (IRPS)
#'
#' @description
#' Core function of the irps package. Runs single-sample GSEA (ssGSEA)
#' on an expression matrix using 18 curated immune gene sets, then
#' computes the IRPS composite score: the mean of 5 favorable gene set
#' scores minus the mean of 3 unfavorable gene set scores.
#'
#' @param expr_mat A numeric matrix of expression values. Rows = features
#'   (proteins or genes), columns = samples. Row names should be feature IDs
#'   matching those in `gene_map`.
#' @param gene_map A data.frame mapping gene symbols to feature IDs in
#'   `expr_mat`. Must contain two columns as specified by `gene_col` and
#'   `feature_col`.
#' @param gene_sets A named list of gene symbol vectors. Default uses the
#'   18 built-in immune gene sets via \code{\link{get_immune_gene_sets}}.
#' @param favorable_sets Character vector naming the gene sets that
#'   contribute positively to IRPS. Default uses the 5 built-in favorable sets.
#' @param unfavorable_sets Character vector naming the gene sets that
#'   contribute negatively to IRPS. Default uses the 3 built-in unfavorable sets.
#' @param gene_col Character, column name in `gene_map` with gene symbols.
#' @param feature_col Character, column name in `gene_map` with feature IDs.
#' @param platform Character, "proteomics" or "transcriptomics".
#' @param groups Optional named character vector mapping sample IDs to
#'   group labels. Group information is attached to the output for
#'   downstream visualization and statistical testing.
#' @param min_set_size Integer, minimum number of mapped features required
#'   to retain a gene set. Default: 3.
#' @param ssgsea_min_size,ssgsea_max_size Integer, passed to
#'   \code{\link[GSVA]{ssgseaParam}} for minSize/maxSize.
#' @param ssgsea_normalize Logical, normalize ssGSEA scores. Default: TRUE.
#' @param verbose Logical, print progress messages. Default: TRUE.
#'
#' @return An S3 object of class `irps`. Use \code{print()},
#'   \code{summary()}, and \code{plot()} to explore results.
#'
#' @export
#'
#' @references
#' The IRPS methodology and gene set curation are described in:
#' (citation to the accompanying thesis/publication).
#'
#' @seealso \code{\link{get_immune_gene_sets}}, \code{\link{decompose_portrait}},
#'   \code{\link{classify_phenotype}}, \code{\link{plot.irps}}
#'
#' @examples
#' \dontrun{
#' # Simulated proteomics data
#' expr <- matrix(rnorm(1000 * 30), nrow = 1000, ncol = 30)
#' rownames(expr) <- paste0("P", seq_len(1000))
#' colnames(expr) <- paste0("Sample", seq_len(30))
#'
#' # Gene-to-protein mapping table
#' gene_map <- data.frame(
#'   Gene = c("B2M", "HLA-A", "STAT1", "CXCL10", "CD47",
#'            "TGFB1", "COL1A1", "S100A8"),
#'   Protein = paste0("P", 1:8),
#'   stringsAsFactors = FALSE
#' )
#' rownames(expr)[1:8] <- gene_map$Protein
#'
#' # Compute IRPS
#' result <- compute_irps(expr, gene_map, platform = "proteomics")
#' print(result)
#' plot(result, type = "boxplot")
#' }
compute_irps <- function(expr_mat,
                          gene_map,
                          gene_sets = get_immune_gene_sets(),
                          favorable_sets = .irps_favorable_sets,
                          unfavorable_sets = .irps_unfavorable_sets,
                          gene_col = "Gene",
                          feature_col = "FeatureID",
                          platform = c("proteomics", "transcriptomics"),
                          groups = NULL,
                          min_set_size = 3L,
                          ssgsea_min_size = 3L,
                          ssgsea_max_size = 500L,
                          ssgsea_normalize = TRUE,
                          verbose = TRUE) {

  platform <- match.arg(platform)

  # ---- Input validation ----
  if (!is.matrix(expr_mat) && !is.data.frame(expr_mat)) {
    stop("'expr_mat' must be a matrix or data.frame")
  }
  if (!is.data.frame(gene_map)) {
    stop("'gene_map' must be a data.frame")
  }
  if (is.null(rownames(expr_mat))) {
    stop("'expr_mat' must have rownames (feature IDs)")
  }
  if (is.null(colnames(expr_mat))) {
    stop("'expr_mat' must have colnames (sample IDs)")
  }
  if (!gene_col %in% colnames(gene_map)) {
    stop(sprintf("Column '%s' not found in gene_map", gene_col))
  }
  if (!feature_col %in% colnames(gene_map)) {
    stop(sprintf("Column '%s' not found in gene_map", feature_col))
  }

  # Convert to matrix if data.frame
  expr_mat <- as.matrix(expr_mat)

  # ---- Map gene sets to features ----
  if (verbose) cat("Mapping gene sets to expression features...\n")
  mapped_sets <- map_genes_to_features(
    gene_sets     = gene_sets,
    gene_map      = gene_map,
    gene_col      = gene_col,
    feature_col   = feature_col,
    min_set_size  = min_set_size
  )

  # Filter mapped sets to those present in expr_mat
  all_rownames <- rownames(expr_mat)
  mapped_sets <- lapply(mapped_sets, function(fids) {
    fids[fids %in% all_rownames]
  })
  mapped_sets <- mapped_sets[lengths(mapped_sets) >= min_set_size]

  n_genes_total <- sum(lengths(mapped_sets))
  if (verbose) {
    cat(sprintf("  %d gene sets retained, %d total feature hits\n",
                length(mapped_sets), n_genes_total))
  }

  if (length(mapped_sets) < 2) {
    stop("Too few gene sets with sufficient mapped features. ",
         "Only ", length(mapped_sets), " set(s) retained; need >= 2.")
  }

  # ---- Determine available favorable/unfavorable sets ----
  available_fav <- intersect(favorable_sets, names(mapped_sets))
  available_unf <- intersect(unfavorable_sets, names(mapped_sets))

  if (verbose) {
    cat(sprintf("Favorable sets: %d/%d available\n",
                length(available_fav), length(favorable_sets)))
    cat(sprintf("Unfavorable sets: %d/%d available\n",
                length(available_unf), length(unfavorable_sets)))
  }

  if (length(available_fav) == 0) {
    stop("No favorable gene sets available. IRPS cannot be computed.")
  }
  if (length(available_unf) == 0) {
    warning("No unfavorable gene sets available. IRPS will use only favorable scores.")
  }

  # ---- Build GSEABase GeneSetCollection ----
  if (verbose) cat("Building GSEABase GeneSetCollection...\n")
  gsc_list <- lapply(names(mapped_sets), function(nm) {
    GSEABase::GeneSet(mapped_sets[[nm]], setName = nm)
  })
  gsc <- GSEABase::GeneSetCollection(gsc_list)

  # ---- Run ssGSEA ----
  if (verbose) cat("Running ssGSEA...\n")
  param <- GSVA::ssgseaParam(
    expr_mat,
    gsc,
    minSize    = ssgsea_min_size,
    maxSize    = ssgsea_max_size,
    normalize  = ssgsea_normalize,
    verbose    = verbose
  )
  ssgsea_scores <- GSVA::gsva(param)
  if (verbose) cat(sprintf("ssGSEA complete: %d sets x %d samples\n",
                           nrow(ssgsea_scores), ncol(ssgsea_scores)))

  # ---- Compute IRPS ----
  if (verbose) cat("Computing IRPS...\n")
  score_t <- as.data.frame(t(ssgsea_scores))
  score_t$SampleID <- rownames(score_t)

  irps_df <- data.frame(
    SampleID       = colnames(ssgsea_scores),
    IRPS_Favorable = colMeans(ssgsea_scores[available_fav, , drop = FALSE]),
    IRPS_Unfavorable = colMeans(ssgsea_scores[available_unf, , drop = FALSE]),
    stringsAsFactors = FALSE
  )
  irps_df$IRPS <- irps_df$IRPS_Favorable - irps_df$IRPS_Unfavorable

  # ---- Attach groups if provided ----
  if (!is.null(groups)) {
    if (is.character(groups) && length(groups) == ncol(expr_mat)) {
      if (is.null(names(groups))) {
        names(groups) <- colnames(expr_mat)
      }
      irps_df$Subtype <- groups[irps_df$SampleID]
    } else {
      warning("'groups' must be a named character vector matching samples. Ignoring.")
    }
  }

  n_features_used <- length(unique(unlist(mapped_sets)))
  if (verbose) {
    cat(sprintf("\nIRPS computation complete. Mean IRPS = %.4f (SD = %.4f)\n",
                mean(irps_df$IRPS, na.rm = TRUE),
                sd(irps_df$IRPS, na.rm = TRUE)))
  }

  # ---- Build and return irps object ----
  new_irps(
    ssgsea_scores     = ssgsea_scores,
    irps_scores       = irps_df,
    gene_sets         = mapped_sets,
    favorable_sets     = available_fav,
    unfavorable_sets   = available_unf,
    platform          = platform,
    n_features_used   = n_features_used,
    params = list(
      min_set_size    = min_set_size,
      ssgsea_min_size = ssgsea_min_size,
      ssgsea_max_size = ssgsea_max_size,
      normalize       = ssgsea_normalize
    )
  )
}
