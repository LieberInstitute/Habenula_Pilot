# Modified version of the violin plot function made by Louise Huuki.
# This is for an rse object rather than a sce object. This was modified
# by Louise and Bukola.
# Sep 30, 2022 - Bukola Ajanaku & Louise Huuki

my_plotExpression <- function(
    sce, genes, assay = "logcounts", ct = "Region", fill_colors = NULL,
    title = NULL
) {
  cat_df <- as.data.frame(colData(sce))[, ct, drop = FALSE]
  expression_long <- reshape2::melt(as.matrix(assays(sce)[[assay]][genes, , drop = FALSE]))
  
  cat <- cat_df[expression_long$Var2, ]
  expression_long <- cbind(expression_long, cat)
  
  expression_violin <- ggplot(
    data = expression_long, aes(x = cat, y = value, fill = cat)
  ) +
    geom_violin(scale = "width") +
    facet_wrap(
      ~Var1, ncol = 5, scales = "free_y"
    ) +
    labs(
      y = paste0("Expression (", assay, ")"),
      title = title
    ) +
    theme_bw(base_size = 35) +
    theme(
      legend.position = "None", axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1),
      strip.text.x = element_text(face = "italic")
    ) +
    stat_summary(
      fun = median,
      # fun.min = median,
      # fun.max = median,
      geom = "crossbar",
      width = 0.3
    )
  
  if (!is.null(fill_colors)) expression_violin <- expression_violin + scale_fill_manual(values = fill_colors)
  
  return(expression_violin)
}
