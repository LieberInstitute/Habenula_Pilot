library(dendextend)
library(dynamicTreeCut)
library(SingleCellExperiment)
library(scRNAseq)
library(batchelor)
library(scater)
library(scran)
library(uwot)
library(DropletUtils)
library(jaffelab)
library(Rtsne)
library(here)
library(utils)


load(here("processed-data","09_snRNA-seq_re-processed","matt_markers","markers.rda"), verbose = TRUE)
load(here("processed-data","09_snRNA-seq_re-processed","05_collapsedClustering.Rda"))


plotExpressionCustom <- function(sce, features, features_name, anno_name = "cellType",
                                 point_alpha=0.2, point_size=0.7, ncol=2, xlab = NULL,
                                 exprs_values = "logcounts", scales = "fixed"){
  scater::plotExpression(sce,
                         exprs_values = exprs_values,
                         features = features,
                         x = anno_name,
                         colour_by = anno_name,
                         ncol = ncol,
                         xlab = xlab,
                         point_alpha = point_alpha,
                         point_size = point_size,
                         add_legend = F,
                         scales = scales) +
  stat_summary(fun = median,
               fun.min = median,
               fun.max = median,
               geom = "crossbar",
               width = 0.3) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        strip.text = element_text(face = "italic")) +
  ggtitle(label=paste0(features_name, " markers"))
}


rownames(sce.all.hb)<-rowData(sce.all.hb)$gene_name
pdf(here("plots","09_snRNA-seq_re-processed", "regionSpecific_Hab-n7_marker-logExprs_collapsedClusters_JMS2022.pdf"), height=6, width=8)
for(i in 1:length(markers.mathys.custom)){
  print(
    plotExpressionCustom(sce = sce.all.hb,
                         features = markers.mathys.custom[[i]],
                         features_name = names(markers.mathys.custom)[[i]],
                         anno_name = "collapsedCluster")
  )
}
dev.off()



print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
