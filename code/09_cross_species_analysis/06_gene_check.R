
library("here")
library("sessioninfo")
library("tidyverse")

load(here("processed-data", 
          "09_cross_species_analysis",
          "Hashikawa_homolog_modeling_results.Rdata"), verbose = TRUE)

mouse_modeling_results$all$enrichment$gene <- tolower(mouse_modeling_results$all$enrichment$gene)
mouse_modeling_results$neuron$enrichment$gene <- tolower(mouse_modeling_results$neuron$enrichment$gene)

search <- list(`Mhb.1` = c("TAC3", "CCK"),
            `Mhb.3` =  c("BHLHE22", "EBF3"),
            `Mhb.2` = c("CHAT"))

map2(search, names(search), function(genes, ct){
  
  hsap_modeling_results$all$enrichment |>
    filter(gene %in% genes) |>
    select(ends_with(ct), gene)
  
})

map2(search, names(search), function(genes, ct){
  
  mouse_modeling_results$neuron$enrichment |>
    filter(gene %in% tolower(genes)) 
})

head(mouse_modeling_results$all$enrichment)
