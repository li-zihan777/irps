#' Classify samples into immune phenotypes
#'
#' @description
#' Classifies each sample into one of three immune phenotypes based on
#' the ssGSEA scores computed by \code{\link{compute_irps}}. The
#' classification is adapted from the Chen & Mellman (Nature, 2017)
#' framework: Immune-Inflamed, Immune-Excluded, and Immune-Desert.
#'
#' @details
#' Because body-fluid proteomics cannot detect canonical T-cell markers
#' (CD3, CD4, CD8, FoxP3), this function uses surrogate axes:
#' \itemize{
#'   \item \strong{Immune Infiltration axis}: Mean of Macrophage, B_cell,
#'     NK_cell, Complement, Monocyte, and TLR_Signaling ssGSEA scores.
#'   \item \strong{Immune Exclusion axis}: Mean of Stromal_CAF,
#'     Immunosuppression, and Adhesion_Migration ssGSEA scores.
#' }
#'
#' Classification rule (based on median thresholds):
#' \itemize{
#'   \item \strong{Immune-Inflamed}: High infiltration AND low exclusion.
#'   \item \strong{Immune-Excluded}: High exclusion (regardless of infiltration).
#'   \item \strong{Immune-Desert}: Low infiltration AND not high exclusion.
#' }
#'
#' @param irps_obj An `irps` object returned by \code{\link{compute_irps}}.
#' @param infiltrate_sigs Character vector of gene set names for the
#'   infiltration axis. If `NULL`, defaults are used.
#' @param exclude_sigs Character vector of gene set names for the
#'   exclusion axis. If `NULL`, defaults are used.
#' @param infil_threshold,excl_threshold Numeric. Thresholds for
#'   infiltration and exclusion. If `NULL` (default), the median of each
#'   axis across all samples is used.
#'
#' @return A data.frame with columns: SampleID, Immune_Infiltration,
#'   Immune_Exclusion, Immune_Phenotype (and Subtype if available).
#'
#' @export
#'
#' @references
#' Chen DS, Mellman I. Elements of cancer immunity and the cancer-immune
#' set point. Nature. 2017;541(7637):321-330.
#'
#' @seealso \code{\link{compute_irps}}, \code{\link{plot_phenotype_bar}}
#'
#' @examples
#' \dontrun{
#' result <- compute_irps(expr, gene_map, platform = "proteomics")
#' pheno <- classify_phenotype(result)
#' head(pheno)
#' table(pheno$Immune_Phenotype)
#' }
classify_phenotype <- function(irps_obj,
                                infiltrate_sigs = NULL,
                                exclude_sigs = NULL,
                                infil_threshold = NULL,
                                excl_threshold = NULL) {
  if (!inherits(irps_obj, "irps")) {
    stop("'irps_obj' must be an irps object from compute_irps()")
  }

  ssgsea <- irps_obj$ssgsea_scores
  ssgsea_t <- as.data.frame(t(ssgsea))
  ssgsea_t$SampleID <- rownames(ssgsea_t)

  # Default infiltration axis
  if (is.null(infiltrate_sigs)) {
    infiltrate_sigs <- c("04_Macrophage", "08_B_cell", "09_NK_cell",
                          "17_Complement", "07_Monocyte", "18_TLR_Signaling")
  }
  infiltrate_sigs <- intersect(infiltrate_sigs, colnames(ssgsea_t))

  # Default exclusion axis
  if (is.null(exclude_sigs)) {
    exclude_sigs <- c("14_Stromal_CAF", "12_Immunosuppression",
                       "16_Adhesion_Migration")
  }
  exclude_sigs <- intersect(exclude_sigs, colnames(ssgsea_t))

  if (length(infiltrate_sigs) == 0) {
    stop("No infiltration gene sets available in the irps object")
  }
  if (length(exclude_sigs) == 0) {
    stop("No exclusion gene sets available in the irps object")
  }

  df <- ssgsea_t
  df$Immune_Infiltration <- rowMeans(
    df[, infiltrate_sigs, drop = FALSE], na.rm = TRUE)
  df$Immune_Exclusion <- rowMeans(
    df[, exclude_sigs, drop = FALSE], na.rm = TRUE)

  infil_med <- infil_threshold %||% stats::median(df$Immune_Infiltration, na.rm = TRUE)
  excl_med  <- excl_threshold  %||% stats::median(df$Immune_Exclusion, na.rm = TRUE)

  df$Immune_Phenotype <- ifelse(
    df$Immune_Infiltration > infil_med & df$Immune_Exclusion < excl_med,
    "Immune-Inflamed",
    ifelse(df$Immune_Exclusion > excl_med,
           "Immune-Excluded", "Immune-Desert")
  )

  # Attach Subtype if present
  if ("Subtype" %in% colnames(irps_obj$irps_scores)) {
    df$Subtype <- irps_obj$irps_scores$Subtype[
      match(df$SampleID, irps_obj$irps_scores$SampleID)]
  }

  df[, c("SampleID", "Immune_Infiltration", "Immune_Exclusion",
         "Immune_Phenotype",
         if ("Subtype" %in% colnames(df)) "Subtype" else NULL)]
}
