
library("tidyverse")
library("ComplexHeatmap")
library("here")

output_path <- here("processed-data", "13_MAGMA","MAGMA_output")

datasets <- c("scz2022", "mdd2019edinburgh", "panic2019", "sud2020op")
names(datasets) <- datasets

magma_out_fn <- map(datasets, ~list.files(path = here(output_path, .x), pattern = "combo.gsa.out", full.names = TRUE))

magma_out <- map_dfr(magma_out_fn[2:4], ~read.table(.x, header = TRUE))

magma_out$dataset <- rep(datasets[2:4], each = length(unique(magma_out$VARIABLE)))

head(magma_out)

magma_out |> filter(VARIABLE == "Hb")

## dirs
plot_dir <- here("plots", "13_MAGMA")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

magama_p_matrix <- magma_out |>
  select(VARIABLE, dataset, P) |>
  pivot_wider(names_from = "dataset", values_from = "P") |>
  column_to_rownames("VARIABLE") |>
  as.matrix()

roworder <- c("Astrocyte","Endo","Excit.Thal","Inhib.Thal","Hb", "LHb", "MHb","Microglia", "OPC","Oligo")

magama_p_matrix <- magama_p_matrix[roworder,]

#              scz2022 mdd2019edinburgh panic2019 sud2020op
# Astrocyte  0.9935300         0.855180   0.44391 0.1784500
# Endo       0.0735120         0.027865   0.29779 0.0022801
# Excit.Thal 0.3703800         0.484070   0.58360 0.6848200
# Inhib.Thal 0.5479500         0.013134   0.26306 0.0969340
# LHb        0.5526100         0.904860   0.42161 0.6809000
# MHb        0.1795000         0.858670   0.34027 0.6494700
# Microglia  0.0050028         0.104720   0.31936 0.4231900
# OPC        0.9881500         0.747370   0.67581 0.6124800
# Oligo      0.0339570         0.347820   0.95005 0.7527400  

mypal = c("white",grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(50))

magama_p_matrix_log <- -log10(magama_p_matrix)
magama_p_matrix_anno <- round(magama_p_matrix_log,3)
magama_p_matrix_anno[magama_p_matrix > 0.05] <- ""
  
pdf(here(plot_dir, "MAGMA_pval_heatmap.pdf"))
Heatmap(magama_p_matrix_log,
        col = mypal,
        name = "-log10(p-val)",
        rect_gp = grid::gpar(col = "black", lwd = 1),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid::grid.text(magama_p_matrix_anno[i, j], x, y, gp = grid::gpar(fontsize = 10))
        })
dev.off()


