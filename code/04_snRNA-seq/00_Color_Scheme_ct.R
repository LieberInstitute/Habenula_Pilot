library("here")
library("jaffelab")
library("RColorBrewer")


cellTypecolors_9 <- 
  c(ExcNeuron = "#274ea6",
    InhibThalMed = "#6ed94f",
    OPC = "#ff67e3",
    Astro = "#296f00",
    Microglia = "#cb0029",
    ThalMedPVT = "#009770",
    Oligo = "#ff815f",
    ExcThalMed = "#d794b4",
    LHbN = "#ffb848")

preview_colors <- function(cell_colors) {
  par(las = 2) # make label text perpendicular to axis
  par(mar = c(5, 8, 4, 2)) # increase y-axis margin.
  barplot(rep(1, length(cell_colors)),
          col = cell_colors,
          horiz = TRUE,
          axes = FALSE,
          names.arg = names(cell_colors)
  )
}

png(here("plots", "cell_colors", "cellTypecolors_9.png"), height = 800)
  preview_colors(cellTypecolors_9)
dev.off()

# saving colors
save(cellTypecolors_9, file = here("processed-data", "cell_type_colors.Rdata"))
