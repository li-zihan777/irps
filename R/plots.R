#' Visualization functions for IRPS results
#'
#' @description
#' A comprehensive set of visualization functions for exploring
#' IRPS results, including boxplots, heatmaps, radar charts,
#' and dimension portrait panels.
#'
#' @name irps_plots
NULL

# ---- Default subtype color palette ----
.default_colors <- c(
  "Subtype1" = "#4472C4",
  "Subtype2" = "#ED7D31",
  "Subtype3" = "#E05C6B"
)

# ---- Internal: run statistical tests for plotting ----
.run_pairwise_stats <- function(df, value_col, group_col) {
  groups <- unique(df[[group_col]])
  if (length(groups) < 2) return(data.frame())

  kw <- stats::kruskal.test(
    as.formula(paste0(value_col, " ~ ", group_col)), data = df)

  pairs <- utils::combn(groups, 2, simplify = FALSE)
  pw <- lapply(pairs, function(p) {
    d <- df[df[[group_col]] %in% p, ]
    wt <- stats::wilcox.test(
      as.formula(paste0(value_col, " ~ ", group_col)),
      data = d, exact = FALSE)
    data.frame(group1 = p[1], group2 = p[2],
               p.value = wt$p.value, stringsAsFactors = FALSE)
  })
  pw <- dplyr::bind_rows(pw)
  if (nrow(pw) > 0) {
    pw$p.adj <- stats::p.adjust(pw$p.value, method = "BH")
    pw$sig <- dplyr::case_when(
      pw$p.adj < 0.001 ~ "***",
      pw$p.adj < 0.01  ~ "**",
      pw$p.adj < 0.05  ~ "*",
      TRUE             ~ "ns"
    )
  }
  pw
}

# ---- Internal: resolve group column ----
.resolve_group <- function(irps_obj) {
  if ("Subtype" %in% colnames(irps_obj$irps_scores)) {
    "Subtype"
  } else {
    NULL
  }
}

#' Plot IRPS boxplot
#'
#' @param irps_obj An `irps` object.
#' @param group_col Character, column in `irps_scores` for grouping.
#'   If `NULL`, auto-detects from the object.
#' @param colors Named character vector of group colors.
#' @param title Optional plot title.
#' @param ... Additional arguments passed to ggplot2.
#'
#' @return A ggplot object.
#' @export
plot_irps_boxplot <- function(irps_obj,
                               group_col = NULL,
                               colors = NULL,
                               title = "Immune Response Potential Score (IRPS)",
                               ...) {
  if (!inherits(irps_obj, "irps")) {
    stop("'irps_obj' must be an irps object")
  }

  group_col <- group_col %||% .resolve_group(irps_obj)
  df <- irps_obj$irps_scores

  if (is.null(colors)) colors <- .default_colors

  if (!is.null(group_col) && group_col %in% colnames(df)) {
    p <- ggplot(df, aes(x = .data[[group_col]], y = .data$IRPS,
                         fill = .data[[group_col]])) +
      geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA, linewidth = 0.6) +
      geom_jitter(aes(color = .data[[group_col]]), width = 0.12, size = 2, alpha = 0.7) +
      scale_fill_manual(values = colors) +
      scale_color_manual(values = colors) +
      labs(title = title, x = NULL, y = "IRPS (Favorable - Unfavorable)") +
      theme_classic(base_size = 13) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"),
            legend.position = "none")

    # Add KW test annotation
    kw <- stats::kruskal.test(
      as.formula(paste0("IRPS ~ ", group_col)), data = df)
    kw_label <- sprintf("KW p = %.2e", kw$p.value)
    y_max <- max(df$IRPS, na.rm = TRUE)
    y_min <- min(df$IRPS, na.rm = TRUE)
    p <- p + annotate("text", x = 2, y = y_max + (y_max - y_min) * 0.15,
                      label = kw_label, size = 3.5, hjust = 0.5, fontface = "italic")
  } else {
    p <- ggplot(df, aes(x = "", y = .data$IRPS)) +
      geom_boxplot(width = 0.5, alpha = 0.6, outlier.shape = NA, linewidth = 0.6,
                   fill = "steelblue") +
      geom_jitter(width = 0.12, size = 2, alpha = 0.7, color = "steelblue") +
      labs(title = title, x = NULL, y = "IRPS (Favorable - Unfavorable)") +
      theme_classic(base_size = 13) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  }
  p
}

#' Plot ssGSEA immune heatmap
#'
#' @param irps_obj An `irps` object.
#' @param group_col Character, column in `irps_scores` for annotation.
#' @param colors Named character vector of group colors.
#' @param show_rownames Logical, show gene set names. Default: TRUE.
#' @param main Title for the heatmap.
#' @param ... Additional arguments passed to \code{\link[pheatmap]{pheatmap}}.
#'
#' @return A pheatmap object (invisibly).
#' @export
plot_irps_heatmap <- function(irps_obj,
                               group_col = NULL,
                               colors = NULL,
                               show_rownames = TRUE,
                               main = "ssGSEA Immune Scores",
                               ...) {
  if (!inherits(irps_obj, "irps")) {
    stop("'irps_obj' must be an irps object")
  }

  group_col <- group_col %||% .resolve_group(irps_obj)
  if (is.null(colors)) colors <- .default_colors

  # Order samples by group then IRPS
  df <- irps_obj$irps_scores
  if (!is.null(group_col) && group_col %in% colnames(df)) {
    sample_order <- df %>%
      dplyr::arrange(.data[[group_col]], dplyr::desc(.data$IRPS)) %>%
      dplyr::pull("SampleID")
    ann_col <- data.frame(
      row.names = df$SampleID,
      Group = df[[group_col]]
    )
    colnames(ann_col) <- group_col
    ann_colors <- list(Group = colors)
    names(ann_colors) <- group_col
    gaps <- cumsum(table(ann_col[[group_col]])[unique(ann_col[[group_col]])])
  } else {
    sample_order <- df$SampleID
    ann_col <- NULL
    ann_colors <- NULL
    gaps <- NULL
  }

  # Match to available samples
  sample_order <- intersect(sample_order, colnames(irps_obj$ssgsea_scores))
  mat <- irps_obj$ssgsea_scores[, sample_order, drop = FALSE]

  pheatmap::pheatmap(
    mat,
    color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(100),
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    show_colnames = FALSE,
    show_rownames = show_rownames,
    annotation_col = ann_col,
    annotation_colors = ann_colors,
    gaps_col = gaps,
    main = main,
    fontsize = 10,
    border_color = NA,
    silent = TRUE,
    ...
  )
}

#' Plot IRPS radar chart
#'
#' @param irps_obj An `irps` object.
#' @param which Character, "favorable" (5 gene sets), "unfavorable" (3),
#'   or "all" (8 gene sets combined; unfavorable sets are inverted).
#' @param colors Named character vector of group colors.
#' @param title Optional title.
#' @param ... Additional arguments (currently unused).
#'
#' @return The radar chart is drawn as a side effect. Returns the
#'   (invisible) data matrix used for plotting.
#' @export
plot_irps_radar <- function(irps_obj,
                             which = c("favorable", "unfavorable", "all"),
                             colors = NULL,
                             title = NULL,
                             ...) {
  which <- match.arg(which)
  if (!inherits(irps_obj, "irps")) {
    stop("'irps_obj' must be an irps object")
  }

  if (is.null(colors)) colors <- .default_colors

  ssgsea <- irps_obj$ssgsea_scores
  irps_df <- irps_obj$irps_scores

  if (which == "favorable") {
    sigs <- irps_obj$favorable_sets
    if (is.null(title)) title <- "IRPS: Favorable Gene Sets"
  } else if (which == "unfavorable") {
    sigs <- irps_obj$unfavorable_sets
    if (is.null(title)) title <- "IRPS: Unfavorable Gene Sets"
  } else {
    sigs <- c(irps_obj$favorable_sets, irps_obj$unfavorable_sets)
    if (is.null(title)) title <- "IRPS: All 8 Immunotherapy-Relevant Sets"
  }
  sigs <- intersect(sigs, rownames(ssgsea))

  if (length(sigs) < 3) {
    stop("Need at least 3 gene sets for radar chart; only ", length(sigs), " available")
  }

  # Aggregate by group
  if ("Subtype" %in% colnames(irps_df)) {
    radar_data <- as.data.frame(t(ssgsea[sigs, , drop = FALSE]))
    radar_data$SampleID <- rownames(radar_data)
    radar_data <- dplyr::left_join(radar_data,
      irps_df[, c("SampleID", "Subtype")], by = "SampleID") %>%
      dplyr::filter(!is.na(.data$Subtype)) %>%
      dplyr::group_by(.data$Subtype) %>%
      dplyr::summarise(
        dplyr::across(dplyr::all_of(sigs), \(x) mean(x, na.rm = TRUE)),
        .groups = "drop"
      )
    radar_mat <- as.data.frame(radar_data[, sigs])
    rownames(radar_mat) <- radar_data$Subtype
  } else {
    # No groups: use overall means
    radar_mat <- as.data.frame(t(colMeans(ssgsea[sigs, , drop = FALSE])))
    rownames(radar_mat) <- "All"
  }

  # For "all", invert unfavorable sets so direction is consistent
  if (which == "all") {
    for (uf in intersect(irps_obj$unfavorable_sets, sigs)) {
      radar_mat[, uf] <- -radar_mat[, uf]
    }
    colnames(radar_mat)[colnames(radar_mat) %in% irps_obj$unfavorable_sets] <-
      paste0(colnames(radar_mat)[colnames(radar_mat) %in% irps_obj$unfavorable_sets], " (inv)")
  }

  # Clean up names
  colnames(radar_mat) <- gsub("^\\d+_", "", colnames(radar_mat))
  colnames(radar_mat) <- strtrim(colnames(radar_mat), 24)

  # Build fmsb format (max, min, data rows)
  rmin <- apply(radar_mat, 2, min, na.rm = TRUE) - 0.12
  rmax <- apply(radar_mat, 2, max, na.rm = TRUE) + 0.12
  radar_fmsb <- suppressWarnings(as.data.frame(rbind(rmax, rmin, radar_mat)))

  group_names <- rownames(radar_mat)
  group_colors <- colors[names(colors) %in% group_names]
  if (length(group_colors) == 0) group_colors <- unname(colors)[seq_len(min(length(colors), nrow(radar_mat)))]

  # Draw radar chart using base graphics (compatible with fmsb package if installed)
  # If fmsb is not installed, provide a message
  if (!requireNamespace("fmsb", quietly = TRUE)) {
    message("Package 'fmsb' is suggested for radar charts. ",
            "Install with: install.packages('fmsb')\n",
            "Returning the data matrix instead.")
    return(invisible(radar_fmsb))
  }

  graphics::par(mar = c(2, 4, 3, 2))
  fmsb::radarchart(
    radar_fmsb,
    axistype = 1,
    pcol = group_colors,
    pfcol = scales::alpha(group_colors, 0.2),
    plwd = 2.8,
    cglcol = "grey60", cglty = 1, cglwd = 0.8,
    axislabcol = "grey40",
    vlcex = 0.72, calcex = 0.6,
    title = title
  )
  graphics::legend(
    x = 1.3, y = 1.08,
    legend = group_names,
    col = group_colors,
    lwd = 2.8, bty = "n", cex = 0.8
  )

  invisible(radar_fmsb)
}

#' Plot 6-dimension immune portrait panel
#'
#' @param portrait_df A data.frame returned by \code{\link{decompose_portrait}}.
#' @param group_col Character, column name for grouping. Default: "Subtype".
#' @param colors Named character vector of group colors.
#' @param title Optional plot title.
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot object.
#' @export
plot_portrait_panel <- function(portrait_df,
                                 group_col = "Subtype",
                                 colors = NULL,
                                 title = "IRPS Decomposition: Immune Microenvironment Portrait",
                                 ...) {
  if (!is.data.frame(portrait_df)) {
    stop("'portrait_df' must be a data.frame from decompose_portrait()")
  }

  dim_cols <- c("D1_MHC_I_Integrity", "D2_IFNgamma_Signal",
                "D3_Cytotoxic_Potential", "D4_Tcell_Recruitment_Surrogate",
                "D5_Immune_Suppression", "D6_Immune_Infiltration")
  dim_cols <- intersect(dim_cols, colnames(portrait_df))

  if (length(dim_cols) == 0) {
    stop("No dimension columns found in portrait_df")
  }

  if (!is.null(group_col) && group_col %in% colnames(portrait_df)) {
    has_group <- TRUE
  } else {
    has_group <- FALSE
  }

  if (is.null(colors)) colors <- .default_colors

  # Pivot to long format
  dim_labels <- c(
    "D1_MHC_I_Integrity"             = "MHC-I Integrity\n(direct)",
    "D2_IFNgamma_Signal"             = "IFN-gamma\nResponse",
    "D3_Cytotoxic_Potential"         = "Cytotoxic\nPotential (surr.)",
    "D4_Tcell_Recruitment_Surrogate" = "T Cell Recruitment\n(surrogate)",
    "D5_Immune_Suppression"          = "Immune\nSuppression (inv)",
    "D6_Immune_Infiltration"         = "Immune\nInfiltration (surr.)"
  )

  portrait_long <- portrait_df
  portrait_long$Group <- if (has_group) {
    portrait_df[[group_col]]
  } else {
    "All"
  }

  portrait_long <- portrait_long %>%
    tidyr::pivot_longer(
      dplyr::all_of(dim_cols),
      names_to = "Dimension",
      values_to = "Score"
    ) %>%
    dplyr::filter(!is.na(.data$Score))

  portrait_long$Dimension <- factor(
    portrait_long$Dimension,
    levels = dim_cols,
    labels = dim_labels[dim_cols]
  )

  # Compute KW p-value per dimension
  if (has_group && length(unique(portrait_long$Group)) > 1) {
    dim_kw <- portrait_long %>%
      dplyr::group_by(.data$Dimension) %>%
      dplyr::summarise(
        kw_p = tryCatch(
          stats::kruskal.test(Score ~ Group, data = dplyr::cur_data())$p.value,
          error = function(e) NA_real_
        ),
        y_max = max(.data$Score, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      dplyr::mutate(kw_label = sprintf("KW p = %.2e", .data$kw_p))
  }

  p <- ggplot(portrait_long, aes(x = .data$Group, y = .data$Score,
                                  fill = .data$Group)) +
    geom_boxplot(width = 0.55, alpha = 0.6, outlier.shape = NA, linewidth = 0.4) +
    geom_jitter(aes(color = .data$Group), width = 0.1, size = 1.2, alpha = 0.4) +
    facet_wrap(~ Dimension, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    labs(x = NULL, y = "Dimension Score", title = title) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(size = 10, face = "bold")
    )

  if (has_group && exists("dim_kw")) {
    p <- p + geom_text(
      data = dim_kw,
      aes(x = length(unique(portrait_long$Group)) / 2 + 0.5,
          y = y_max * 1.15, label = kw_label),
      inherit.aes = FALSE, size = 3.2, fontface = "italic", color = "grey30"
    )
  }

  p
}

#' Plot immune phenotype distribution bar chart
#'
#' @param irps_obj An `irps` object, or a data.frame from
#'   \code{\link{classify_phenotype}}.
#' @param group_col Character, column for grouping. Default: "Subtype".
#' @param colors Named character vector of phenotype colors.
#' @param title Optional plot title.
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot object.
#' @export
plot_phenotype_bar <- function(irps_obj,
                                group_col = NULL,
                                colors = NULL,
                                title = "Immune Phenotype Distribution",
                                ...) {
  if (is.data.frame(irps_obj)) {
    pheno_df <- irps_obj
  } else if (inherits(irps_obj, "irps")) {
    pheno_df <- classify_phenotype(irps_obj)
  } else {
    stop("Input must be an 'irps' object or a data.frame from classify_phenotype()")
  }

  if (is.null(colors)) {
    colors <- c(
      "Immune-Inflamed" = "#E05C6B",
      "Immune-Excluded" = "#ED7D31",
      "Immune-Desert"   = "#92C5DE"
    )
  }

  group_col <- group_col %||% "Subtype"
  if (!group_col %in% colnames(pheno_df)) {
    stop(sprintf("Column '%s' not found in phenotype data", group_col))
  }

  prop_df <- pheno_df %>%
    dplyr::count(.data[[group_col]], .data$Immune_Phenotype) %>%
    dplyr::group_by(.data[[group_col]]) %>%
    dplyr::mutate(prop = n / sum(n)) %>%
    dplyr::ungroup()

  pheno_levels <- c("Immune-Inflamed", "Immune-Excluded", "Immune-Desert")
  prop_df$Immune_Phenotype <- factor(
    prop_df$Immune_Phenotype, levels = pheno_levels)

  ggplot(prop_df, aes(x = .data[[group_col]], y = prop,
                       fill = .data$Immune_Phenotype)) +
    geom_bar(stat = "identity", position = "fill", width = 0.6,
             color = "white", linewidth = 0.3) +
    scale_fill_manual(values = colors, name = "Immune Phenotype") +
    scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
    labs(x = NULL, y = "Proportion", title = title) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(face = "bold", size = 12)
    )
}

#' Plot IRPS vs key biomarker scatter
#'
#' @param irps_obj An `irps` object.
#' @param expr_mat Expression matrix used for compute_irps().
#' @param gene_map Data frame mapping gene symbols to feature IDs.
#' @param gene Character, gene symbol of the biomarker to plot (e.g., "B2M").
#' @param gene_col,feature_col Column names in `gene_map`.
#' @param colors Named character vector of group colors.
#' @param ... Additional arguments (currently unused).
#'
#' @return A ggplot object.
#' @export
plot_irps_scatter <- function(irps_obj, expr_mat, gene_map, gene,
                               gene_col = "Gene", feature_col = "FeatureID",
                               colors = NULL, ...) {
  if (!inherits(irps_obj, "irps")) {
    stop("'irps_obj' must be an irps object")
  }

  if (is.null(colors)) colors <- .default_colors

  # Resolve gene to feature ID
  fid <- gene_map[[feature_col]][gene_map[[gene_col]] == gene]
  if (length(fid) == 0) {
    stop(sprintf("Gene '%s' not found in gene_map", gene))
  }
  fid <- fid[1]

  if (!fid %in% rownames(expr_mat)) {
    stop(sprintf("Feature '%s' (for gene '%s') not found in expr_mat", fid, gene))
  }

  df <- data.frame(
    SampleID   = colnames(expr_mat),
    Expression = as.numeric(expr_mat[fid, ]),
    IRPS       = irps_obj$irps_scores$IRPS[
      match(colnames(expr_mat), irps_obj$irps_scores$SampleID)],
    stringsAsFactors = FALSE
  )

  if ("Subtype" %in% colnames(irps_obj$irps_scores)) {
    df$Subtype <- irps_obj$irps_scores$Subtype[
      match(df$SampleID, irps_obj$irps_scores$SampleID)]
  }

  df <- df[complete.cases(df[, c("Expression", "IRPS")]), ]

  rho <- stats::cor(df$IRPS, df$Expression, method = "spearman")
  p_val <- stats::cor.test(df$IRPS, df$Expression, method = "spearman")$p.value

  p <- ggplot(df, aes(x = .data$Expression, y = .data$IRPS))

  if ("Subtype" %in% colnames(df)) {
    p <- p + geom_point(aes(color = .data$Subtype), size = 2.5, alpha = 0.75) +
      scale_color_manual(values = colors)
  } else {
    p <- p + geom_point(size = 2.5, alpha = 0.75, color = "steelblue")
  }

  p <- p +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE,
                color = "grey40", linewidth = 0.8, linetype = "dashed") +
    annotate("text",
             x = min(df$Expression, na.rm = TRUE),
             y = max(df$IRPS, na.rm = TRUE),
             label = sprintf("Spearman rho = %.3f\np = %.2e", rho, p_val),
             hjust = 0, vjust = 1, size = 3.5, fontface = "italic") +
    labs(
      title = sprintf("IRPS vs %s Expression", gene),
      x = sprintf("%s Expression", gene),
      y = "IRPS"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  p
}
