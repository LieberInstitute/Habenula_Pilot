library(here)
library(tidyverse)
library(SummarizedExperiment)
library(WGCNA)
library(jaffelab)
library(sessioninfo)

rse_path = here(
    'processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda'
)
tom_dir = here('processed-data', '19_wgcna', 'TOM_files')
out_path = here('processed-data', '19_wgcna', 'modules.rds')
deg_covariates = c(
    'PrimaryDx', 'AgeDeath', 'Flowcell', 'mitoRate', 'rRNA_rate', 'RIN',
    'totalAssignedGene', 'abs_ERCCsumLogErr', 'qSV1', 'qSV2', 'qSV3', 'qSV4',
    'qSV5', 'qSV6', 'qSV7', 'qSV8', 'tot.Hb', 'tot.Thal'
)

dir.create(tom_dir, recursive = TRUE, showWarnings = FALSE)

#   Load gene expression data and filter genes by minimum mean RPKM
rse_gene = get(load(rse_path))
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, length_var = 'Length')
rse_gene = rse_gene[rowMeans(assays(rse_gene)$rpkm) > 0.25,]

mod_deg <- model.matrix(
    as.formula(paste('~', paste(deg_covariates, collapse = " + "))),
    data = colData(rse_gene)
)

#   Log transform, regress out covariates, and transpose
exp_mat = log2(assays(rse_gene)$rpkm + 1) |>
    cleaningY(mod_deg, P = 3) |>
    t()

powers = c(c(1:10), seq(from = 12, to = 20, by = 2))
thres = pickSoftThreshold(exp_mat, powerVector = powers, verbose = 5)

net = blockwiseModules(
    exp_mat, power = thres$powerEstimate, TOMType = "unsigned",
    minModuleSize = 30, reassignThreshold = 0, mergeCutHeight = 0.25,
	numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE,
    verbose = 3, maxBlockSize = 30000,
    saveTOMFileBase = file.path(tom_dir, 'gene')
)
save(net, file = out_path)

session_info()
