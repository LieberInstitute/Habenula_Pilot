
library("tidyverse")
library("ComplexHeatmap")
library("here")

output_path <- here("processed-data", "13_MAGMA","MAGMA_output")

datasets <- c("scz2022", "mdd2019edinburgh", "panic2019", "sud2020op")
names(datasets) <- datasets

magma_out_fn <- map_chr(datasets, ~list.files(path = here(output_path, .x), pattern = ".gsa.out", full.names = TRUE))

magama_out <- map_dfr(magma_out_fn, ~read.table(.x, header = TRUE))

magama_out$dataset <- rep(datasets, each = 9)

## dirs
plot_dir <- here("plots", "09_bulk_DE", "07_DE_plots")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

magama_out |>
  select(VARIABLE, dataset, P) |>
  pivot_wider(names_from = "dataset", values_from = "P")

# VARIABLE   scz2022 mdd2019edinburgh panic2019 sud2020op
# <chr>        <dbl>            <dbl>     <dbl>     <dbl>
# 1 Astrocyte  0.994             0.855      0.444   0.178  
# 2 Endo       0.0735            0.0279     0.298   0.00228
# 3 Excit.Thal 0.370             0.484      0.584   0.685  
# 4 Inhib.Thal 0.548             0.0131     0.263   0.0969 
# 5 LHb        0.553             0.905      0.422   0.681  
# 6 MHb        0.180             0.859      0.340   0.649  
# 7 Microglia  0.00500           0.105      0.319   0.423  
# 8 OPC        0.988             0.747      0.676   0.612  
# 9 Oligo      0.0340            0.348      0.950   0.753  
