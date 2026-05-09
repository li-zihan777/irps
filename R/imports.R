# Package-level imports. These prevent R CMD check notes about
# undefined global variables or functions.

# ---- Standard R packages ----
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics legend par
#' @importFrom stats as.formula complete.cases cor cor.test kruskal.test
#'   median p.adjust sd setNames wilcox.test
#' @importFrom utils combn head
NULL

# ---- dplyr pipe and verbs ----
#' @importFrom dplyr %>%
#' @importFrom dplyr .data
#' @importFrom dplyr across all_of arrange bind_rows desc
#'   everything filter group_by inner_join left_join mutate n pull
#'   rename select summarise
NULL

# ---- ggplot2 ----
#' @importFrom ggplot2 aes annotate element_blank element_rect element_text
#'   facet_wrap geom_bar geom_boxplot geom_jitter geom_point geom_smooth
#'   geom_text ggplot labs scale_color_manual scale_fill_manual
#'   scale_y_continuous theme theme_classic xlab ylab
NULL

# ---- tidyr ----
#' @importFrom tidyr pivot_longer
NULL

# ---- tibble ----
#' @importFrom tibble rownames_to_column column_to_rownames
NULL
