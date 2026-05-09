#' Retrieve the built-in 18 immune gene sets
#'
#' @description
#' Returns a named list of 18 curated immune gene sets used for
#' ssGSEA-based immune microenvironment profiling. Each set contains
#' HGNC gene symbols selected from established immunology literature.
#'
#' The 5 sets contributing positively to IRPS (favorable):
#' \enumerate{
#'   \item Antigen Presentation (01_Antigen_Presentation)
#'   \item MHC-I Core (02_MHC_I_Core)
#'   \item IFN-gamma Response (03_IFNgamma_Response)
#'   \item Cytotoxic Surrogate (15_Cytotoxic_Surrogate)
#'   \item Chemokine/Inflammatory (13_Chemokine_Inflammatory)
#' }
#'
#' The 3 sets contributing negatively to IRPS (unfavorable):
#' \enumerate{
#'   \item Immunosuppression (12_Immunosuppression)
#'   \item Stromal/CAF (14_Stromal_CAF)
#'   \item MDSC/Neutrophil (05_MDSC_Neutrophil)
#' }
#'
#' @return A named list of 18 character vectors, each containing gene symbols.
#'
#' @export
#'
#' @examples
#' gs <- get_immune_gene_sets()
#' names(gs)
#' lengths(gs)
#' head(gs[["02_MHC_I_Core"]])
get_immune_gene_sets <- function() {
  .irps_gene_sets
}

#' Map gene symbols to feature IDs in an expression matrix
#'
#' @description
#' Given a list of gene sets (gene symbols) and a mapping table, this
#' function maps each gene to its corresponding feature ID in the
#' expression matrix. Only features present in the expression matrix
#' are retained, and gene sets with too few mapped features are
#' optionally dropped.
#'
#' @param gene_sets A named list of gene symbol vectors (as returned by
#'   \code{\link{get_immune_gene_sets}}).
#' @param gene_map A data.frame with at least two columns: one containing
#'   gene symbols and one containing feature IDs present in the expression
#'   matrix. Default column names are "Gene" and "FeatureID".
#' @param gene_col Character, column name in `gene_map` containing gene symbols.
#' @param feature_col Character, column name in `gene_map` containing feature IDs.
#' @param min_set_size Integer, minimum number of mapped features required
#'   to retain a gene set. Default: 3.
#'
#' @return A named list of character vectors, each containing feature IDs
#'   that match to the expression matrix.
#'
#' @export
#'
#' @examples
#' gs <- get_immune_gene_sets()
#' # Example gene-to-protein mapping
#' mapping <- data.frame(
#'   Gene = c("B2M", "HLA-A", "HLA-B", "CD8A"),
#'   Protein = c("P61769", "P04439", "P01889", "P01732"),
#'   stringsAsFactors = FALSE
#' )
#' mapped <- map_genes_to_features(gs, mapping,
#'   gene_col = "Gene", feature_col = "Protein")
#' lengths(mapped)
map_genes_to_features <- function(gene_sets, gene_map,
                                   gene_col = "Gene",
                                   feature_col = "FeatureID",
                                   min_set_size = 3L) {
  if (!is.list(gene_sets)) {
    stop("'gene_sets' must be a named list of gene symbol vectors")
  }
  if (!is.data.frame(gene_map)) {
    stop("'gene_map' must be a data.frame")
  }
  if (!gene_col %in% colnames(gene_map)) {
    stop(sprintf("Column '%s' not found in gene_map", gene_col))
  }
  if (!feature_col %in% colnames(gene_map)) {
    stop(sprintf("Column '%s' not found in gene_map", feature_col))
  }

  gene_map <- gene_map[!is.na(gene_map[[gene_col]]) &
                       !is.na(gene_map[[feature_col]]), ]
  gene_map <- gene_map[!duplicated(gene_map[[gene_col]]), ]

  mapped <- lapply(gene_sets, function(genes) {
    idx <- match(genes, gene_map[[gene_col]])
    fids <- gene_map[[feature_col]][idx]
    fids[!is.na(fids)]
  })

  sizes <- lengths(mapped)
  dropped <- names(mapped)[sizes < min_set_size]
  if (length(dropped) > 0) {
    message(sprintf("Dropping %d gene set(s) with < %d features: %s",
                    length(dropped), min_set_size,
                    paste(dropped, collapse = ", ")))
  }

  mapped <- mapped[sizes >= min_set_size]

  if (length(mapped) < 1) {
    stop("Too few gene sets with sufficient mapped features; check your gene_map")
  }

  mapped
}

# ---- Internal: 18 built-in gene sets ----
.irps_gene_sets <- list(
  "01_Antigen_Presentation" = c(
    "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F",
    "B2M", "TAP1", "TAP2", "TAPBP", "CALR", "CANX",
    "ERAP1", "ERAP2", "PSMB8", "PSMB9", "PSMB10",
    "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
    "HLA-DQA1", "HLA-DQB1", "CD74"
  ),
  "02_MHC_I_Core" = c(
    "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP"
  ),
  "03_IFNgamma_Response" = c(
    "STAT1", "STAT2", "STAT3", "IRF3", "IRF7", "IRF9",
    "ISG15", "GBP1", "GBP2", "GBP4", "GBP5",
    "JAK1", "JAK2", "IFNGR1", "IFNGR2", "CXCL10", "NLRC5"
  ),
  "04_Macrophage" = c(
    "CD14", "FCGR3A", "ITGAM", "LST1", "AIF1",
    "TYROBP", "FCER1G", "CTSS", "LYZ", "SPI1",
    "C1QA", "C1QB", "C1QC", "APOE", "CD68", "CSF1R"
  ),
  "05_MDSC_Neutrophil" = c(
    "S100A8", "S100A9", "S100A12", "ARG1", "CXCR2",
    "MPO", "ELANE", "FCGR3B", "CEACAM8"
  ),
  "06_Neutrophil" = c(
    "FCGR3B", "CXCR2", "MPO", "ELANE", "CEACAM8",
    "MMP8", "MMP9", "OLFM4"
  ),
  "07_Monocyte" = c(
    "CD14", "CD33", "FCGR3A", "ITGAM", "S100A8", "S100A9",
    "CCR2", "CSF1R", "CX3CR1"
  ),
  "08_B_cell" = c(
    "CD37", "MZB1", "JCHAIN", "SDC1", "CD79A", "CD79B", "MS4A1"
  ),
  "09_NK_cell" = c(
    "TYROBP", "FCER1G", "NCAM1", "KLRD1", "KLRK1",
    "CD247", "NCR1", "NCR3"
  ),
  "10_Eosinophil" = c(
    "RNASE2", "RNASE3", "EPX", "PRG2", "CLC", "CCR3", "SIGLEC8"
  ),
  "11_DC_cell" = c(
    "ITGAX", "WDFY4", "FCER1A", "LAMP3", "CCR7",
    "CD1C", "BATF3", "CLEC10A", "CLEC9A", "XCR1"
  ),
  "12_Immunosuppression" = c(
    "CD47", "VSIR", "TGFB1", "TGFB2", "TGFB3",
    "ARG1", "NT5E", "S100A8", "S100A9",
    "COL1A1", "COL1A2", "COL3A1", "FN1", "TNC"
  ),
  "13_Chemokine_Inflammatory" = c(
    "CXCL1", "CXCL2", "CXCL5", "CXCL8", "CXCL10", "CXCL17",
    "IL1B", "TNF", "IL6"
  ),
  "14_Stromal_CAF" = c(
    "COL1A1", "COL1A2", "COL3A1", "FN1", "TNC",
    "ACTA2", "TAGLN", "FAP", "PDGFRB",
    "TGFB1", "TGFB2", "TGFB3", "VEGFC"
  ),
  "15_Cytotoxic_Surrogate" = c(
    "SERPINB9", "TNFSF10", "TNFRSF10A", "TNFRSF10B",
    "FCER1G", "TYROBP", "CTSS", "LYZ", "CST7"
  ),
  "16_Adhesion_Migration" = c(
    "ICAM1", "ITGAL", "ITGAM", "ITGB2", "SELL", "CD44", "CD47"
  ),
  "17_Complement" = c(
    "C1QA", "C1QB", "C1QC", "C2", "C3", "C4A", "C4B", "C5",
    "CFB", "CFH", "CFI"
  ),
  "18_TLR_Signaling" = c(
    "TLR1", "TLR2", "TLR3", "TLR4", "TLR6", "TLR8",
    "MYD88", "CD14", "LYZ", "IRAK1", "IRAK4", "TICAM1"
  )
)

# IRPS configuration
.irps_favorable_sets <- c(
  "01_Antigen_Presentation", "02_MHC_I_Core",
  "03_IFNgamma_Response", "15_Cytotoxic_Surrogate",
  "13_Chemokine_Inflammatory"
)

.irps_unfavorable_sets <- c(
  "12_Immunosuppression", "14_Stromal_CAF", "05_MDSC_Neutrophil"
)
