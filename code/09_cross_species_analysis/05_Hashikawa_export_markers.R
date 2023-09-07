
library("spatialLIBD")
library("here")
library("sessioninfo")
library("tidyverse")
library("xlsx")

load(here("processed-data", 
          "09_cross_species_analysis",
          "Hashikawa_homolog_modeling_results.Rdata"), verbose = TRUE)
# hsap_modeling_results

## export in excel for easy access

enrich_mouse <- map(mouse_modeling_results,"enrichment")
enrich_hsap <- map(hsap_modeling_results, "enrichment")

names(enrich_mouse) <- paste0("mouse_",names(enrich_mouse))
names(enrich_hsap) <- paste0("hsap_",names(enrich_hsap))

enrich <- c(enrich_mouse, enrich_hsap)
map(enrich, dim)



output_dir <- here("processed-data", "09_cross_species_analysis","Hashikawa_cross_species_enrichment")

export_enrich <- function(enrich_df = enrich$mouse_all, 
                          tag = "mouse_all"){
  message("export: ", tag)
  
  cell_types <- colnames(enrich_df)[grep("t_stat", colnames(enrich_df))]
  cell_types <- gsub("t_stat_", "", cell_types)
  names(cell_types) <- cell_types
  
  ct_enrich <- map(cell_types, function(ct) {
    enrich_df |> 
      select(gene, JAX_geneID = ensembl, ends_with(ct)) |>
      rownames_to_column("Ensembl") |>
      rename_all(~stringr::str_replace(., paste0("_",ct),"")) |>
      mutate(cell_type = ct) |>
      arrange(fdr) |>
      mutate(rank = row_number()) |>
      filter(rank <= 100)
  })
  
  walk2(ct_enrich, names(ct_enrich), ~write.xlsx(.x,
                                          file = here(output_dir, paste0("Hashikawa_enrichment-",tag,".xlsx")),
                                          sheetName= .y, 
                                          col.names=TRUE, 
                                          row.names=FALSE, 
                                          append=TRUE))
  }



walk2(enrich, names(enrich), ~export_enrich(.x, .y))



