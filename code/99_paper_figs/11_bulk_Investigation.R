## May 15, 2023 - Bukola Ajanaku
# Probing deconvoluted data to find trends alongside certain phenotypes.
# qrsh -l mem_free=20G,h_vmem=20G

library(SummarizedExperiment)
library(here)
library(tidyverse)
library(ggplot2)
library(jaffelab)

# loading deconvo data
load(file = here("processed-data", "99_paper_figs", "sce_objects", "prop_long.RDATA"),
     verbose = TRUE)

head(prop_long, 10)
    # Sample cellType     prop BrNum PrimaryDx factor_CT  Hb_sum Br_Order prop_perc
    # <chr>  <chr>       <dbl> <chr> <chr>     <fct>       <dbl> <fct>        <dbl>
    #   1 R18424 Astrocyte  0      Br55… Control   Astrocyte 0       Br5572        0   
    # 2 R18424 Endo       0      Br55… Control   Endo      0       Br5572        0   
    # 3 R18424 Microglia  0      Br55… Control   Microglia 0       Br5572        0   
    # 4 R18424 Oligo      0      Br55… Control   Oligo     0       Br5572        0   
    # 5 R18424 OPC        0      Br55… Control   OPC       0       Br5572        0   
    # 6 R18424 Inhib.Thal 1      Br55… Control   Inhib.Th… 0       Br5572      100   
    # 7 R18424 Excit.Thal 0      Br55… Control   Excit.Th… 0       Br5572        0   
    # 8 R18424 MHb        0      Br55… Control   MHb       0       Br5572        0   
    # 9 R18424 LHb        0      Br55… Control   LHb       0       Br5572        0   
    # 10 R18358 Astrocyte  0.0335 Br14… Schizo    Astrocyte 0.00897 Br1427        3.35

# loading final rse object 
load(here("processed-data", "02_bulk_qc", "count_data_bukola", 
          "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"), verbose = TRUE)
# rse_gene

# creating plotting directory
plot_dir <- here("plots", "99_paper_figs", "11_bulk_Investigation")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# adding Thalamus proportions and then dividing Habenula over Thalamus
sum_Prop <- prop_long |>
  filter(cellType %in% c("Excit.Thal", "Inhib.Thal")) |>
  group_by(BrNum) |>
  summarize(Thal_sum = sum(prop)) 

prop_long <- left_join(prop_long, sum_Prop) |>
  arrange(Br_Order) |>
  mutate(Hb_over_Thal = (Hb_sum / Thal_sum) )
    # %>% print(width = Inf)

# dropping confusing prop_perc as this was the percentages of just Hb_sum
prop_long$prop_perc <- NULL

# making data frame for ease of merge with rse ColData
comp_invest <- as.data.frame(prop_long[,c("Sample", "Hb_sum", "Thal_sum", "Hb_over_Thal",
                                          "BrNum")])
comp_invest$RNum = comp_invest$Sample


new_rse_colData <- as_tibble(colData(rse_gene))
new_rse_colData <- new_rse_colData |>
                    left_join(comp_invest[match(unique(comp_invest$RNum), comp_invest$RNum),],
                              copy = TRUE
                              )

new_rse_colData <- DataFrame(new_rse_colData)

# setting new colData 
colData(rse_gene) <- new_rse_colData

############## GRABBING PC VALUES ##############################################
# grabbing PC values 
pca <- prcomp(t(assays(rse_gene)$logcounts))
  
  ## % of the variance explained by each PC
pca_vars <- getPcaVars(pca)
pca_vars_labs<- paste0("PC", seq(along = pca_vars), ": ",
                       pca_vars, "% Var Expl")

  ## Joining PC and sample info
pca_data<-cbind(pca$x, colData(rse_gene))
  
  ## Add samples' phenotypes
pca_data<-as.data.frame(pca_data)
  
##### PLOTTING
pdf(file = here(plot_dir, "PC1_Investigation.pdf")) 

pos = position_jitter(width = 0.1, height = 0.1, seed = 1)

  ggplot(pca_data, 
         aes(x = PC1,
         y = Hb_sum,
         color = PrimaryDx),
         ) + 
    geom_jitter(position = pos) +
    ggtitle("Hb_sum_VS_Primary_Dx") +
    geom_text(position = pos,
              aes(color = "grey"),
              label = pca_data$BrNum)
  
dev.off()






# Grabbing metrics of interest
phenoMets <- c("AgeDeath", "PrimaryDx", "Flowcell", "mitoRate",
               "RIN")
# mainMets <c("Hb_sum", "Thal_sum", "Hb_over_Thal")










