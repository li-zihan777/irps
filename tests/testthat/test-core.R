library(testthat)
library(irps)

test_that("get_immune_gene_sets returns all 18 gene sets", {
  gs <- get_immune_gene_sets()
  expect_type(gs, "list")
  expect_equal(length(gs), 18)
  expect_true(all(grepl("^\\d{2}_", names(gs))))
  expect_true(all(lengths(gs) >= 5))
})

test_that("map_genes_to_features works correctly", {
  gs <- get_immune_gene_sets()
  gene_map <- data.frame(
    Gene      = c("B2M", "HLA-A", "HLA-B", "STAT1", "CXCL10", "CD8A", "CD47"),
    FeatureID = c("P001", "P002", "P003", "P004", "P005", "P006", "P007"),
    stringsAsFactors = FALSE
  )

  mapped <- map_genes_to_features(gs, gene_map)

  expect_type(mapped, "list")
  expect_true(all(lengths(mapped) > 0))
  # B2M should map to P001 in 02_MHC_I_Core
  expect_true("P001" %in% mapped[["02_MHC_I_Core"]])
})

test_that("map_genes_to_features drops sets below min_set_size", {
  gs <- get_immune_gene_sets()

  # With only 1 gene mapped, all 18 sets will have < 3 features
  gene_map <- data.frame(
    Gene      = c("B2M"),
    FeatureID = c("P001"),
    stringsAsFactors = FALSE
  )

  expect_error(
    map_genes_to_features(gs, gene_map, min_set_size = 3),
    "Too few gene sets"
  )

  # With 7 genes mapped, some sets will be retained and some dropped
  gene_map2 <- data.frame(
    Gene      = c("B2M", "HLA-A", "HLA-B", "STAT1", "CXCL10", "CD8A", "CD47"),
    FeatureID = c("P001", "P002", "P003", "P004", "P005", "P006", "P007"),
    stringsAsFactors = FALSE
  )
  expect_message(
    map_genes_to_features(gs, gene_map2, min_set_size = 3),
    "Dropping"
  )
})

test_that("map_genes_to_features errors on invalid input", {
  expect_error(
    map_genes_to_features("not_a_list", data.frame(Gene = "A", Protein = "P1")),
    "must be a named list"
  )
  expect_error(
    map_genes_to_features(list(A = "B2M"), "not_a_df"),
    "must be a data.frame"
  )
  expect_error(
    map_genes_to_features(list(A = "B2M"), data.frame(Gene = "A", Protein = "P1"),
      gene_col = "Missing"),
    "not found in gene_map"
  )
})

test_that("compute_irps runs end-to-end on simulated data", {
  set.seed(2024)
  n_prot <- 300
  n_samp <- 30

  expr <- matrix(rnorm(n_prot * n_samp, mean = 10, sd = 2),
                 nrow = n_prot, ncol = n_samp)
  rownames(expr) <- paste0("P", seq_len(n_prot))
  colnames(expr) <- paste0("S", seq_len(n_samp))

  gs <- get_immune_gene_sets()
  all_genes <- unique(unlist(gs))

  gene_map <- data.frame(
    Gene      = all_genes[1:min(120, length(all_genes))],
    FeatureID = paste0("P", seq_len(min(120, length(all_genes)))),
    stringsAsFactors = FALSE
  )

  result <- compute_irps(expr, gene_map, platform = "proteomics", verbose = FALSE)

  expect_s3_class(result, "irps")
  expect_true(is.matrix(result$ssgsea_scores))
  expect_true(is.data.frame(result$irps_scores))
  expect_true("IRPS" %in% colnames(result$irps_scores))
  expect_true("IRPS_Favorable" %in% colnames(result$irps_scores))
  expect_true("IRPS_Unfavorable" %in% colnames(result$irps_scores))
  expect_equal(nrow(result$irps_scores), n_samp)
  expect_true(length(result$favorable_sets) > 0)
  expect_true(length(result$unfavorable_sets) >= 0)
})

test_that("compute_irps with groups", {
  set.seed(2024)
  expr <- matrix(rnorm(200 * 20, mean = 10, sd = 2),
                 nrow = 200, ncol = 20)
  rownames(expr) <- paste0("P", seq_len(200))
  colnames(expr) <- paste0("S", seq_len(20))

  gs <- get_immune_gene_sets()
  all_genes <- unique(unlist(gs))
  gene_map <- data.frame(
    Gene      = all_genes[1:100],
    FeatureID = paste0("P", seq_len(100)),
    stringsAsFactors = FALSE
  )

  groups <- setNames(rep(c("High", "Low"), each = 10), colnames(expr))

  result <- compute_irps(expr, gene_map, groups = groups,
                         platform = "proteomics", verbose = FALSE)

  expect_s3_class(result, "irps")
  expect_true("Subtype" %in% colnames(result$irps_scores))
  expect_equal(unique(result$irps_scores$Subtype), c("High", "Low"))
})

test_that("decompose_portrait returns correct structure", {
  set.seed(2024)
  expr <- matrix(rnorm(200 * 20, mean = 10, sd = 2),
                 nrow = 200, ncol = 20)
  rownames(expr) <- paste0("P", seq_len(200))
  colnames(expr) <- paste0("S", seq_len(20))

  gs <- get_immune_gene_sets()
  all_genes <- unique(unlist(gs))
  gene_map <- data.frame(
    Gene      = all_genes[1:100],
    FeatureID = paste0("P", seq_len(100)),
    stringsAsFactors = FALSE
  )

  result <- compute_irps(expr, gene_map, verbose = FALSE)
  port <- decompose_portrait(result, expr_mat = expr, gene_map = gene_map)

  expect_true(is.data.frame(port))
  expected_dims <- c("D1_MHC_I_Integrity", "D2_IFNgamma_Signal",
                     "D3_Cytotoxic_Potential", "D4_Tcell_Recruitment_Surrogate",
                     "D5_Immune_Suppression", "D6_Immune_Infiltration")
  for (d in expected_dims) {
    expect_true(d %in% colnames(port))
  }
})

test_that("classify_phenotype returns correct structure", {
  set.seed(2024)
  expr <- matrix(rnorm(200 * 20, mean = 10, sd = 2),
                 nrow = 200, ncol = 20)
  rownames(expr) <- paste0("P", seq_len(200))
  colnames(expr) <- paste0("S", seq_len(20))

  gs <- get_immune_gene_sets()
  all_genes <- unique(unlist(gs))
  gene_map <- data.frame(
    Gene      = all_genes[1:100],
    FeatureID = paste0("P", seq_len(100)),
    stringsAsFactors = FALSE
  )

  result <- compute_irps(expr, gene_map, verbose = FALSE)
  pheno <- classify_phenotype(result)

  expect_true(is.data.frame(pheno))
  expect_true("Immune_Phenotype" %in% colnames(pheno))
  expect_true(all(pheno$Immune_Phenotype %in%
    c("Immune-Inflamed", "Immune-Excluded", "Immune-Desert")))
})

test_that("print, summary, and plot S3 methods work", {
  set.seed(2024)
  expr <- matrix(rnorm(200 * 20, mean = 10, sd = 2),
                 nrow = 200, ncol = 20)
  rownames(expr) <- paste0("P", seq_len(200))
  colnames(expr) <- paste0("S", seq_len(20))

  gs <- get_immune_gene_sets()
  all_genes <- unique(unlist(gs))
  gene_map <- data.frame(
    Gene      = all_genes[1:100],
    FeatureID = paste0("P", seq_len(100)),
    stringsAsFactors = FALSE
  )

  result <- compute_irps(expr, gene_map, verbose = FALSE)

  expect_output(print(result), "IRPS Object")
  expect_type(summary(result), "list")
  expect_true("irps_mean" %in% names(summary(result)))

  # plot.irps with various types
  expect_s3_class(plot(result, type = "boxplot"), "ggplot")
  expect_s3_class(plot_irps_boxplot(result), "ggplot")

  pheno <- classify_phenotype(result)
  expect_s3_class(plot_phenotype_bar(pheno, group_col = "Immune_Phenotype"), "ggplot")
})

test_that("gene set content is biologically valid", {
  gs <- get_immune_gene_sets()

  # MHC-I core set
  expect_true(all(c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP")
                  %in% gs[["02_MHC_I_Core"]]))

  # Immunosuppression set includes key players
  expect_true(all(c("TGFB1", "CD47", "ARG1") %in% gs[["12_Immunosuppression"]]))

  # Chemokine set
  expect_true("CXCL10" %in% gs[["13_Chemokine_Inflammatory"]])
})
