
library("spatialLIBD")
library("here")
library("sessioninfo")
library("purrr")
library("ComplexHeatmap")

load(here("processed-data", 
     "09_cross_species_analysis",
     "Hashikawa_homolog_modeling_results.Rdata"), verbose = TRUE)
# hsap_modeling_results
# mouse_modeling_results

## human habenula t-stats
registration_t_stats <- map(hsap_modeling_results, function(mod){
  t_stats <- mod$enrichment[, grep("^t_stat", colnames(mod$enrichment))]
  colnames(t_stats) <- gsub("^t_stat_", "", colnames(t_stats))
  rownames(t_stats) <- mod$enrichment$ensembl ## ensembl is JAX.geneID from match(model_results$ensembl, rownames(stats))
  return(t_stats)
})

map(registration_t_stats, ~.x[1:5,1:5])

cor_species <- map2(mouse_modeling_results, registration_t_stats, 
                    ~layer_stat_cor(modeling_results = .x,
                                    stats = .y,
                                    model_type = "enrichment",
                                    top_n = 100
                    ))

cor_species[["allXneuron"]] <- layer_stat_cor(modeling_results = mouse_modeling_results$neuron,
                                              stats = registration_t_stats$all,
                                              model_type = "enrichment",
                                              top_n = 100)

map(cor_species, dim)

## Subset our neuron to just Habenula subtypes

walk2(cor_species, names(cor_species), function(cor, n){
  pdf(here("plots", "09_cross_species_analysis",paste0("Hashikawa_mouse-",n,".pdf")))
  layer_stat_cor_plot(cor, max = max(cor))
  dev.off()
})

#### Complex Heatmap versions ####
source(file = here("code", "99_paper_figs", "source_colors.R"))
sn_colors

cell_types <- rownames(cor_species$all)

ct_color_bar <- columnAnnotation(
  " " = cell_types,
  col = list(" " = sn_colors[cell_types]),
  show_legend = FALSE
)

theSeq <- seq(min(cor_species$all), max(cor_species$all), by = 0.01)
my.col <-
  grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "PRGn"))(length(theSeq))

sort(colnames(cor_species$all))

hsap_ct_order <- c("Astrocyte","Endo","Microglia","Oligo" ,"OPC",
                   "Excit.Thal", "Inhib.Thal", 
                   "MHb.1", "MHb.2", "MHb.3",
                   "LHb.1", "LHb.2", "LHb.3" , "LHb.4", "LHb.5", "LHb.6","LHb.7")

mm_ct_order <- c("Astrocyte1",  "Astrocyte2",  "Endothelial", "Epen", "Microglia",
                 "Mural","Oligo1","Oligo2", "Oligo3","OPC1", "OPC2", "OPC3",
                 "Neuron1","Neuron2","Neuron3","Neuron4","Neuron5","Neuron6","Neuron7","Neuron8")

names(cor_species)

combined_cor <- rbind(t(cor_species$all[hsap_ct_order,mm_ct_order]),
                      t(cor_species$allXneuron[hsap_ct_order,]))

pdf(here("plots", "09_cross_species_analysis","Hashikawa_mouse_complex-All.pdf"))
## combined and ordered 
Heatmap(combined_cor,
        name = "Cor",
        col = my.col,
        row_split = c(rep("All", ncol(cor_species$all)), rep("Neuron", ncol(cor_species$allXneuron))),
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        bottom_annotation = ct_color_bar,
        rect_gp = gpar(col = "black", lwd = 1))
## All vs All
Heatmap(t(cor_species$all),
        name = "Cor",
        col = my.col,
        bottom_annotation = ct_color_bar,
        rect_gp = gpar(col = "black", lwd = 1))
## All vs. Neuron
Heatmap(t(cor_species$allXneuron),
        name = "Cor",
        col = my.col,
        bottom_annotation = ct_color_bar,
        rect_gp = gpar(col = "black", lwd = 1))

dev.off()

## Neuron vs. Neuron
cell_types_neuron <- rownames(cor_species$neuron)

ct_color_bar_neuron <- columnAnnotation(
  " " = cell_types_neuron ,
  col = list(" " = sn_colors[cell_types_neuron ]),
  show_legend = FALSE
)

pdf(here("plots", "09_cross_species_analysis","Hashikawa_mouse_complex-Neuron.pdf"))
## unordred
Heatmap(t(cor_species$neuron[sort(rownames(cor_species$neuron)),]),
        name = "Cor",
        col = my.col,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        bottom_annotation = ct_color_bar_neuron,
        rect_gp = gpar(col = "black", lwd = 1))
## ordered
Heatmap(t(cor_species$neuron),
        name = "Cor",
        col = my.col,
        bottom_annotation = ct_color_bar_neuron,
        rect_gp = gpar(col = "black", lwd = 1))

dev.off()

## Reproducibility information
print("Reproducibility information:")
options(width = 120)
session_info()

