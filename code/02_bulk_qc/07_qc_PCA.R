## Louise Huuki-Myers, Feb 2025
## Check quality metrics vs. PCs (reviews)

library("SummarizedExperiment")
library("tidyverse")
library("scran")
library("jaffelab")
library("sessioninfo")
library("here")
library("ggrepel")
library("GGally")
## plot dir

plot_dir <- here("plots", "02_bulk_qc", "07_qc_PCA")
if(!dir.exists(plot_dir)){
    dir.create(plot_dir)
}

# loading final rse object
load(here("processed-data", "02_bulk_qc", "count_data_bukola",
          "rse_gene_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata"), verbose = TRUE)

# loading deconvo data
load(file = here("processed-data", "99_paper_figs", "sce_objects", "prop_long.RDATA"),
     verbose = TRUE)

sum_Prop <- prop_long |>
    filter(cellType %in% c("Excit.Thal", "Inhib.Thal")) |>
    group_by(BrNum) |>
    summarize(Thal_sum = sum(prop))

## grab QC values
pd <- as.data.frame(colData(rse_gene)) |>
    left_join(sum_Prop) |>
    mutate(PrimaryDx =  gsub("Schizo", "SCZD", PrimaryDx))

qc_variables <- c("numReads", "numMapped", "numUnmapped", "overallMapRate", "concordMapRate", "totalMapped", "mitoMapped", "mitoRate", "rRNA_rate", "totalAssignedGene","Thal_sum")

pd_qc_long <- pd |>
    select(RNum, BrNum, PrimaryDx, all_of(qc_variables)) |>
    pivot_longer(!c(RNum,BrNum, PrimaryDx), names_to = "qc_var", values_to = "qc_value") |>
    mutate(qc_var = factor(qc_var, levels = qc_variables))

#### Calc PCs #####
set.seed(250218)
# grabbing PC values
pca <- prcomp(t(assays(rse_gene)$logcounts))

## % of the variance explained by each PC
pca_vars <- getPcaVars(pca)
pca_vars_labs<- paste0("PC", seq(along = pca_vars), ": ",
                       pca_vars, "% Var Expl")

colnames(pca$x) <- pca_vars_labs

## PC1 vs PC2

pc1v2 <- pca$x |>
    as_tibble()|>
    add_column(PrimaryDx = pd$PrimaryDx,
               BrNum = pd$BrNum) |>
    ggplot(aes(x = `PC1: 16.7% Var Expl`, y = `PC2: 8.61% Var Expl`)) +
    geom_point(aes(color = PrimaryDx)) +
    geom_text_repel(aes(label = ifelse(BrNum %in% c("Br5572", "Br5459"), BrNum, ""))) +
    theme_bw() +
    coord_equal()

ggsave(pc1v2, filename= here(plot_dir, "pc1v2.png"), height = 5)


pca_long <- reshape2::melt(pca$x[,1:10]) |>
    rename(RNum = Var1, PCA = Var2, pca_value = value) |>
    left_join(pd_qc_long, relationship = "many-to-many")

pca_qc_cor <- pca_long |>
    group_by(PCA, qc_var) |>
    summarize(cor = cor(pca_value, qc_value))

pca_qc_cor |> arrange(-abs(cor))

#### Plot PCAs vs. QCs ####

pca_vs_qc <- pca_long |>
    filter(grepl("[1|2|3|4|5]:", PCA)) |>
    ggplot(aes(x = pca_value, y = qc_value, color = PrimaryDx)) +
    geom_point(size = 1.5) +
    geom_text_repel(aes(label = BrNum), size = 1.5) +
    facet_grid(qc_var~PCA, scales = "free") +
    theme_bw()

ggsave(pca_vs_qc, filename = here(plot_dir, "pca_vs_qc.png"), width = 10, height = 13)


pca_vs_qc_cor  <-  pca_qc_cor |>
    ggplot(aes(x = qc_var, y = PCA, fill = cor)) +
    geom_tile() +
    geom_text(aes(label = round(cor, 2)), size = 2) +
    theme_bw() +
    scale_fill_gradient2(limits = c(-1, 1)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(pca_vs_qc_cor, filename = here(plot_dir , "pca_vs_qc_cor.png"))

#### GGPAIR ####
identical(rownames(pd), rownames(pca$x))

pca_top5 <- cbind(pd[,c("BrNum", "PrimaryDx", "Flowcell")], pca$x[,1:5])

ggpair_pca <- ggpairs(pca_top5,
                      columns = colnames(pca$x)[1:5],
                      aes(colour = PrimaryDx))

ggsave(ggpair_pca, filename = here(plot_dir, "ggpair_pca.png"))
