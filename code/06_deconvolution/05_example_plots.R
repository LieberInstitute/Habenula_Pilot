
library("tidyverse")
library("sessioninfo")
library("here")

## prep dirs ##
plot_dir <- here("plots", "06_deconvolution", "05_example_plots")
if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

set.seed(3042023)

example_colors2 <- c(Hb = "#B33D90", Thal = "#E69965", Other = "#59ABB1")

## est prop
set.seed(512)
example_prop <- matrix(data = rep(c(6,3), 4), nrow = 4, byrow = TRUE) - matrix(rnorm(8), nrow = 4)
example_prop[1,1] <- example_prop[1,1] -1
example_prop <- cbind(example_prop, 10 - rowSums(example_prop))/10
rowSums(example_prop)

colnames(example_prop) <- c("Hb", "Thal", "Other")
rownames(example_prop) <- paste0("sample_", 1:nrow(example_prop))
example_prop


example_prop_long <- reshape2::melt(example_prop) |>
    rename(Sample = Var1, cell_type = Var2, prop = value) |>
    mutate(cell_type = factor(cell_type, levels = c("Hb", "Thal", "Other")))

example_prop_plot <- DeconvoBuddies::plot_composition_bar(example_prop_long,
                                                          sample_col = "Sample",
                                                          x_col = "Sample",
                                                          add_text = FALSE) +
    scale_fill_manual(values = example_colors2) +
    theme_void()

ggsave(example_prop_plot, filename = here(plot_dir, "example_composition.png"), height = 3, width = 3)


