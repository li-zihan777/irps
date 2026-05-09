#' Decompose IRPS into biologically interpretable dimensions
#'
#' @description
#' Decomposes the IRPS composite score into 6 interpretable dimensions
#' for detailed characterization of the immune microenvironment. Each
#' dimension combines direct protein expression data with ssGSEA scores
#' where applicable. All dimensions are oriented such that higher values
#' indicate a more favorable immune microenvironment for immunotherapy.
#'
#' @details
#' The 6 dimensions are:
#' \describe{
#'   \item{D1_MHC_I_Integrity}{Direct z-score mean of MHC-I pathway proteins
#'     (HLA-A, HLA-C, B2M, TAP1, TAP2). Higher = better antigen presentation.}
#'   \item{D2_IFNgamma_Signal}{50% direct IFN-gamma proteins (CXCL10, STAT1)
#'     + 50% ssGSEA IFNgamma_Response gene set.}
#'   \item{D3_Cytotoxic_Potential}{ssGSEA Cytotoxic_Surrogate score.
#'     Note: canonical cytotoxic markers (GZMB, PRF1) are often undetectable
#'     in body-fluid proteomics; this uses detected surrogate proteins.}
#'   \item{D4_Tcell_Recruitment_Surrogate}{Mean of ssGSEA Chemokine_Inflammatory
#'     and Adhesion_Migration scores. Note: CD8A/B not detected in most
#'     body-fluid proteomics; chemokine/adhesion serves as surrogate.}
#'   \item{D5_Immune_Suppression}{Inverted composite: 40% direct immunosuppression
#'     proteins (TGFB1, COL1A1, S100A9, S100A8) + 60% ssGSEA Immunosuppression
#'     and Stromal_CAF. Higher = less suppression (better).}
#'   \item{D6_Immune_Infiltration}{Mean of ssGSEA Macrophage, Monocyte, and
#'     Complement scores. Note: CD3/CD4/CD8 not detected in body-fluid
#'     proteomics; myeloid/complement used as infiltration surrogates.}
#' }
#'
#' @param irps_obj An `irps` object returned by \code{\link{compute_irps}}.
#' @param expr_mat Optional numeric matrix of expression values (features ×
#'   samples) used for direct protein-level dimension components. If `NULL`,
#'   only ssGSEA-based dimensions are computed.
#' @param gene_map Optional data.frame mapping gene symbols to feature IDs.
#'   Required if `expr_mat` is provided.
#' @param gene_col,feature_col Column names in `gene_map` for gene symbols
#'   and feature IDs.
#'
#' @return A data.frame with columns: SampleID, Subtype (if available),
#'   D1-D6 dimension scores, and IRPS.
#'
#' @export
#'
#' @seealso \code{\link{compute_irps}}, \code{\link{plot_portrait_panel}}
#'
#' @examples
#' \dontrun{
#' result <- compute_irps(expr, gene_map, platform = "proteomics")
#' portrait <- decompose_portrait(result, expr_mat = expr, gene_map = gene_map)
#' head(portrait)
#' plot_portrait_panel(portrait)
#' }
decompose_portrait <- function(irps_obj,
                                expr_mat = NULL,
                                gene_map = NULL,
                                gene_col = "Gene",
                                feature_col = "FeatureID") {
  if (!inherits(irps_obj, "irps")) {
    stop("'irps_obj' must be an irps object from compute_irps()")
  }

  ssgsea_mat <- irps_obj$ssgsea_scores
  sample_ids <- colnames(ssgsea_mat)
  n_samples <- length(sample_ids)

  portrait <- data.frame(
    SampleID = sample_ids,
    stringsAsFactors = FALSE
  )

  has_direct <- !is.null(expr_mat) && !is.null(gene_map)

  # ---- D1: MHC-I Integrity (direct protein z-score) ----
  if (has_direct) {
    mhc1_genes <- c("HLA-A", "HLA-C", "B2M", "TAP1", "TAP2")
    mhc1_fids <- map_genes_to_features(
      list(MHC_I = mhc1_genes), gene_map, gene_col, feature_col, min_set_size = 1
    )
    if (length(mhc1_fids) > 0) {
      mhc1_fids <- intersect(mhc1_fids[[1]], rownames(expr_mat))
      if (length(mhc1_fids) >= 2) {
        mhc1_sub <- expr_mat[mhc1_fids, , drop = FALSE]
        mhc1_z <- t(scale(t(mhc1_sub)))
        portrait$D1_MHC_I_Integrity <- colMeans(mhc1_z, na.rm = TRUE)
      } else {
        portrait$D1_MHC_I_Integrity <- NA_real_
      }
    } else {
      portrait$D1_MHC_I_Integrity <- NA_real_
    }
  } else {
    portrait$D1_MHC_I_Integrity <- NA_real_
  }

  # ---- D2: IFN-gamma Signal (50% direct + 50% ssGSEA) ----
  ifng_direct <- rep(0, n_samples)
  if (has_direct) {
    ifng_genes <- c("CXCL10", "STAT1")
    ifng_fids <- map_genes_to_features(
      list(IFNG = ifng_genes), gene_map, gene_col, feature_col, min_set_size = 1
    )
    if (length(ifng_fids) > 0) {
      ifng_fids <- intersect(ifng_fids[[1]], rownames(expr_mat))
      if (length(ifng_fids) > 0) {
        ifng_sub <- expr_mat[ifng_fids, , drop = FALSE]
        if (length(ifng_fids) > 1) {
          ifng_z <- t(scale(t(ifng_sub)))
          ifng_direct <- colMeans(ifng_z, na.rm = TRUE)
        } else {
          ifng_direct <- as.numeric(scale(as.numeric(ifng_sub)))
        }
      }
    }
  }
  ifng_ssgsea <- rep(0, n_samples)
  if ("03_IFNgamma_Response" %in% rownames(ssgsea_mat)) {
    ifng_ssgsea <- as.numeric(ssgsea_mat["03_IFNgamma_Response", ])
  }
  portrait$D2_IFNgamma_Signal <- ifng_direct * 0.5 + ifng_ssgsea * 0.5

  # ---- D3: Cytotoxic Potential (ssGSEA surrogate) ----
  if ("15_Cytotoxic_Surrogate" %in% rownames(ssgsea_mat)) {
    portrait$D3_Cytotoxic_Potential <- as.numeric(
      ssgsea_mat["15_Cytotoxic_Surrogate", ])
  } else {
    portrait$D3_Cytotoxic_Potential <- NA_real_
  }

  # ---- D4: T Cell Recruitment Surrogate ----
  chemo <- rep(0, n_samples)
  if ("13_Chemokine_Inflammatory" %in% rownames(ssgsea_mat)) {
    chemo <- as.numeric(ssgsea_mat["13_Chemokine_Inflammatory", ])
  }
  adhesion <- rep(0, n_samples)
  if ("16_Adhesion_Migration" %in% rownames(ssgsea_mat)) {
    adhesion <- as.numeric(ssgsea_mat["16_Adhesion_Migration", ])
  }
  portrait$D4_Tcell_Recruitment_Surrogate <- (chemo + adhesion) / 2

  # ---- D5: Immune Suppression (inverted, higher = better) ----
  imm_direct <- rep(0, n_samples)
  if (has_direct) {
    supp_genes <- c("TGFB1", "COL1A1", "S100A9", "S100A8")
    supp_fids <- map_genes_to_features(
      list(Supp = supp_genes), gene_map, gene_col, feature_col, min_set_size = 1
    )
    if (length(supp_fids) > 0) {
      supp_fids <- intersect(supp_fids[[1]], rownames(expr_mat))
      if (length(supp_fids) > 0) {
        supp_sub <- expr_mat[supp_fids, , drop = FALSE]
        if (length(supp_fids) > 1) {
          supp_z <- t(scale(t(supp_sub)))
          imm_direct <- colMeans(supp_z, na.rm = TRUE)
        } else {
          imm_direct <- as.numeric(scale(as.numeric(supp_sub)))
        }
      }
    }
  }
  supp_ssgsea <- rep(0, n_samples)
  if ("12_Immunosuppression" %in% rownames(ssgsea_mat)) {
    supp_ssgsea <- supp_ssgsea + as.numeric(ssgsea_mat["12_Immunosuppression", ])
  }
  if ("14_Stromal_CAF" %in% rownames(ssgsea_mat)) {
    supp_ssgsea <- supp_ssgsea + as.numeric(ssgsea_mat["14_Stromal_CAF", ])
  }
  portrait$D5_Immune_Suppression <- -(imm_direct * 0.4 + supp_ssgsea * 0.6)

  # ---- D6: Immune Infiltration (surrogate) ----
  infil_sigs <- intersect(
    c("04_Macrophage", "07_Monocyte", "17_Complement"),
    rownames(ssgsea_mat)
  )
  if (length(infil_sigs) > 0) {
    portrait$D6_Immune_Infiltration <- colMeans(
      ssgsea_mat[infil_sigs, , drop = FALSE], na.rm = TRUE
    )
  } else {
    portrait$D6_Immune_Infiltration <- NA_real_
  }

  # ---- Attach Subtype and IRPS ----
  if ("Subtype" %in% colnames(irps_obj$irps_scores)) {
    portrait$Subtype <- irps_obj$irps_scores$Subtype[
      match(portrait$SampleID, irps_obj$irps_scores$SampleID)]
  }
  portrait$IRPS <- irps_obj$irps_scores$IRPS[
    match(portrait$SampleID, irps_obj$irps_scores$SampleID)]

  portrait
}
