
library("SummarizedExperiment")
library("here")
library("sessioninfo")
library("recount")
library("jaffelab")
library("tidyverse")
library("GGally")

## dirs
data_dir <- here("processed-data", "03_bulk_pca", "02_multiregion_PCA")
if(!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

plot_dir <- here("plots", "03_bulk_pca", "02_multiregion_PCA")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

#### load data ####
## retreiv colData final rse 
load(here( "processed-data","rse_objects","rse_gene_Habenula_Pilot.rda"),verbose = TRUE)
pd <- colData(rse_gene)

## unfiltered rse_gene
load(here("preprocessed_data", "rse_gene_Roche_Habenula_PairedEnd_n73.Rdata")) # gene info
rse_gene_Hb <- rse_gene
dim(rse_gene_Hb)
# [1] 58037    73

rse_gene_Hb <- rse_gene_Hb[,rownames(pd)]
colData(rse_gene_Hb) <- pd

## control only
rse_gene_Hb <- rse_gene_Hb[,rse_gene_Hb$PrimaryDx == "Control"]
dim(rse_gene_Hb)
# [1] 58037    33

## other region data 
load("/dcs04/lieber/lcolladotor/dbDev_LIBD001/subsets/for_Louise_230628_nonHabenula/rse_gene_ribozero_nonHabenula_n784.rda", verbose = TRUE)
dim(rse_gene)
# [1] 58037   784

table(rse_gene$Region)
# Amygdala      BLA       CA  Caudate     dACC       DG    DLPFC    HIPPO      MeA     mPFC     sACC 
#      140       54       11       82       55       41      121       26       55       57      142 

table(rse_gene$Dataset)
# Astellas_DG   BrainSeq_Phase2_DLPFC   BrainSeq_Phase2_HIPPO BrainSeq_Phase3_Caudate     BrainSeq_Phase4and5 
#          30                      66                      12                      82                      32 
# psychENCODE_Mood         PTSD_BrainOmics                 VA_PTSD 
#              268                      75                     219 


colnames(colData(rse_gene))
# [1] "SAMPLE_ID"         "RNum"              "RIN"               "Region"            "Dataset"          
# [6] "BrNum"             "Dx"                "Age"               "Sex"               "Race"             
# [11] "Protocol"          "numReads"          "numMapped"         "numUnmapped"       "mitoMapped"       
# [16] "totalMapped"       "overallMapRate"    "concordMapRate"    "mitoRate"          "rRNA_rate"        
# [21] "totalAssignedGene" "bamFile" 

colnames(colData(rse_gene_Hb))

colnames(colData(rse_gene))[!colnames(colData(rse_gene)) %in% colnames(colData(rse_gene_Hb))]
# [1] "SAMPLE_ID" "Region"    "Dataset"   "Dx"        "Age"       "Protocol" 

## Add missing cols
rse_gene_Hb$SAMPLE_ID <- rse_gene_Hb$RNum
rse_gene_Hb$Region <- "Hb"
rse_gene_Hb$Dataset <- "Habenula_Pilot"
rse_gene_Hb$Protocol <- "RiboZeroGold"

colnames(colData(rse_gene_Hb))[which("AgeDeath" == colnames(colData(rse_gene_Hb)))] <- 'Age'
colnames(colData(rse_gene_Hb))[which("PrimaryDx" == colnames(colData(rse_gene_Hb)))] <- 'Dx'

## Add missing cols to Hb data
common_cols <- intersect(colnames(colData(rse_gene_Hb)), colnames(colData(rse_gene)))

colData(rse_gene_Hb) <- colData(rse_gene_Hb)[, common_cols]
colData(rse_gene) <- colData(rse_gene)[, common_cols]

## compare rowData
rowData(rse_gene_Hb)
rowData(rse_gene)

rownames(rse_gene) <- rowData(rse_gene)$gencodeID
## rowData is in same order
identical(rownames(rse_gene), rownames(rse_gene_Hb))

#### combine rse_gene 
(rse_gene <- cbind(rse_gene, rse_gene_Hb))
# class: RangedSummarizedExperiment 
# dim: 58037 817 
# metadata(0):
#   assays(1): counts
# rownames(58037): ENSG00000223972.5 ENSG00000227232.5 ... ENSG00000210195.2 ENSG00000210196.2
# rowData names(10): Length gencodeID ... gencodeTx meanExprs
# colnames(817): R10713_PilotRepeat R11449 ... R18422 R18423
# colData names(22): BrNum RNum ... Dataset Protocol

table(rse_gene$Region)

#### Run PCA ####
## get expression
gene_rpkm <- getRPKM(rse_gene, "Length")
gene_rpkm_filter <- gene_rpkm[rowMeans(gene_rpkm) > 0.1, ]

geneExprs_filter <- log2(gene_rpkm_filter + 1)
dim(geneExprs_filter)
# [1] 29507   817
pca <- prcomp(t(geneExprs_filter))
pca_vars <- getPcaVars(pca)
pca_vars_lab <- paste0(
  "PC", seq(along = pca_vars), ": ",
  pca_vars, "% Var Expl"
)

dim(pca$x)

pca_tab <- as.data.frame(colData(rse_gene)) |> 
  cbind(pca$x[, 1:10]) |> 
  mutate(Region = factor(Region, levels = c("Amygdala", "BLA","CA","MeA", "DG", "HIPPO", "dACC","sACC","DLPFC","mPFC","Caudate","Hb")))

colnames(pca_tab)
## create long table
pca_vars_lab_tb <- tibble(PC = ss(pca_vars_lab,":"), var_expl = pca_vars_lab) |>
  mutate(var_expl = fct_reorder(var_expl, as.integer(gsub("PC","", PC))))

levels(pca_vars_lab_tb$var_expl)[1:10]

pca_long <- pca_tab |> 
  pivot_longer(!c(1:22), names_to = "PC", values_to = "PC_val") |>
  left_join(pca_vars_lab_tb)

## save data 
save(pca_tab, pca_long, file = here(data_dir, "Multi_region_PCs.Rdata"))
# load(here(data_dir, "Multi_region_PCs.Rdata"), verbose = TRUE)
# pca_vars_lab <- unique(pca_long$var_expl)

#### PCA plots ####
"#ff8032"

region_colors <- c(Amygdala = "#ff9ccb",
                   BLA = "#c10040",
                   CA = "#ff7168",
                   MeA = "#b9008b",
                   DG = "#00960e",
                   HIPPO = "#99C71A",
                   dACC = "#0094fc",
                   sACC = "#014abf",
                   DLPFC = "#c495ff",
                   mPFC = "#8330b6",
                   Caudate = "#65717B",
                   # Hb = "#F2CA18",
                   # Hb = "#F9F50D"
                   Hb = "#FA9A09"
                   )

## ggpairs for pca ##
gg_pca <- ggpairs(pca_tab,
                  mapping = aes(color = Region),
                  columns = paste0("PC", 1:6)) +
  scale_fill_manual(values = region_colors)+
  scale_color_manual(values = region_colors)

ggsave(gg_pca, filename = here(plot_dir, "ggpairs_pca.png"), height = 11, width = 11)


pc_test <- pca_tab |>
  ggplot(aes(x = PC1, y = PC2, )) +
  geom_point(shape = 21, aes(fill = Region, color = Region == "Hb")) +
  theme_bw() +
  scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "#00000000"), guide = "none") +
  scale_fill_manual(values = region_colors) +
  labs(x = pca_vars_lab[[1]], y = pca_vars_lab[[2]]) 

ggsave(pc_test, filename = here(plot_dir, "Bulk_PC1vPC2_Region.png"))


## boxplots
pca_boxplots <- pca_long |>
  ggplot(aes(x = Region, y = PC_val, fill = Region)) +
  geom_boxplot()  +
  scale_fill_manual(values = region_colors) +
  facet_wrap(~var_expl, ncol = 2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None")
  
ggsave(pca_boxplots, filename = here(plot_dir, "Bulk_PCA_boxplots.png"))

pca_boxplots_free <- pca_long |>
  ggplot(aes(x = Region, y = PC_val, fill = Region)) +
  geom_boxplot()  +
  scale_fill_manual(values = region_colors) +
  facet_wrap(~var_expl, ncol = 2, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None")

ggsave(pca_boxplots_free, filename = here(plot_dir, "Bulk_PCA_boxplots_free.png"))

pc6_boxplot <- pca_long |>
  filter(PC == "PC6") |>
  ggplot(aes(x = Region, y = PC_val)) +
  geom_boxplot(aes(fill = Region), alpha = 0.5, outlier.shape = NA)  +
  geom_jitter(aes(fill = Region),shape = 21, width = 0.25) +
  scale_fill_manual(values = region_colors) +
  # facet_wrap(~var_expl, ncol = 2, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "None") +
  labs(y = pca_vars_lab[[6]])

ggsave(pc6_boxplot, filename = here(plot_dir, "Bulk_PC6_boxplot.png"), height = 4, width = 5)
ggsave(pc6_boxplot, filename = here(plot_dir, "Bulk_PC6_boxplot.pdf"), height = 4, width = 5)

pc_1v6 <- pca_tab  |>
  ggplot(aes(x = PC1, y = PC6)) +
  geom_point(shape = 21, aes(fill = Region, color = Region == "Hb")) +
  theme_bw() +
  scale_color_manual(values = c(`TRUE` = "black", `FALSE` = "#00000000"), guide = "none") +
  scale_fill_manual(values = region_colors) +
  labs(x = pca_vars_lab[[1]], y = pca_vars_lab[[6]]) 

ggsave(pc_1v6, filename = here(plot_dir, "Bulk_PC1vPC6_Region.png"), width = 5, height = 4)
ggsave(pc_1v6, filename = here(plot_dir, "Bulk_PC1vPC6_Region.pdf"), width = 5, height = 4)

## combine scatter + boxplot
library(patchwork)

pc1_patch <- pc_1v6 + theme(legend.position = "None") + 
  pc6_boxplot + theme(axis.title.y=element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y=element_blank())
 
ggsave(pc1_patch, filename = here(plot_dir, "Bulk_PC6_Region_combo.png"), height = 4)
ggsave(pc1_patch, filename = here(plot_dir, "Bulk_PC6_Region_combo.pdf"), height = 4)

# slurmjobs::job_single(name = "02_multiregion_PCA", memory = "10G", cores = 1, create_shell = TRUE, command = "Rscript 02_multiregion_PCA.R")

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()

