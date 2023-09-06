
library("spatialLIBD")
library("here")
library("sessioninfo")
library("purrr")

load(here("processed-data", 
     "09_cross_species_analysis",
     "Hashikawa_homolog_modeling_results.Rdata"), verbose = TRUE)
# hsap_modeling_results
# mouse_modeling_results

## human habenula t-stats
registration_t_stats <- hsap_modeling_results$enrichment[, grep("^t_stat", colnames(hsap_modeling_results$enrichment))]
colnames(registration_t_stats) <- gsub("^t_stat_", "", colnames(registration_t_stats))
rownames(registration_t_stats) <- hsap_modeling_results$enrichment$ensembl

cor_species <- map(mouse_modeling_results, ~layer_stat_cor(
  modeling_results = .x,
  stats = registration_t_stats,
  model_type = "enrichment",
  top_n = 100
))

## Subset our neuron to just Habenula subtypes
rownames(cor_species$neuron)[grep("Hb",rownames(cor_species$neuron))]
cor_species$neuron <- cor_species$neuron[grep("Hb",rownames(cor_species$neuron)),]

walk2(cor_species, names(cor_species), function(cor, n){
  pdf(here("plots", "09_cross_species_analysis",paste0("Hashikawa_mouse-",n,".pdf")))
  layer_stat_cor_plot(cor, max = max(cor))
  dev.off()
})


## Reproducibility information
print("Reproducibility information:")
options(width = 120)
session_info()

