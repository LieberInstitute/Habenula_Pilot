library(here)
library(tidyverse)
library(SummarizedExperiment)
library(jaffelab)
library(clusterProfiler)
library(sessioninfo)

net_path = here('processed-data', '19_wgcna', 'modules.rds')
rse_path = here(
    'processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda'
)
plot_dir = here('plots', '19_wgcna')
protected_deg_covariates = c('PrimaryDx', 'AgeDeath')

set.seed(0)
dir.create(plot_dir, showWarnings = FALSE)

#   Load WGCNA network and gene expression data, then filter genes by minimum
#   mean RPKM
net = readRDS(net_path)
rse_gene = get(load(rse_path))
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, length_var = 'Length')
rse_gene = rse_gene[rowMeans(assays(rse_gene)$rpkm) > 0.25,]

################################################################################
#   Plot module weights by diagnosis
################################################################################

me_df = net$MEs |>
    rownames_to_column('RNum') |>
    as_tibble() |>
    select(-ME0) |>
    left_join(
        colData(rse_gene) |>
            as_tibble(),
        by = 'RNum'
    ) |>
    mutate(
        PrimaryDx = factor(
            ifelse(PrimaryDx == "Schizo", "SCZD", "Control"),
            levels = c("Control", "SCZD")
        )
    )

plot_list = list()
p_val_list = list()
for (i in seq_len(length(grep('^ME[0-9]+$', colnames(me_df))))) {
    #   Get p-value of linear relationship with diagnosis
    lin_mod = lm(
        as.formula(
            paste0(
                'ME', i, ' ~ ', 
                paste(protected_deg_covariates, collapse = " + "), 
                ' - 1'
            )
        ),
        data = me_df
    )
    p_val_list[[i]] = summary(lin_mod)$coef['PrimaryDxSCZD', 4]
    
    plot_list[[i]] = ggplot(
            me_df,
            aes(
                x = PrimaryDx, y = !!sym(paste0('ME', i)), color = PrimaryDx
            )
        ) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter() +
        geom_text(
            label = paste('\np =', signif(p_val_list[[i]], 2), ''),
            x = Inf, y = Inf, hjust = 1, vjust = 1, color = 'black', size = 8
        ) +
        guides(color = "none") +
        labs(x = "Diagnosis", y = paste("Module", i)) +
        theme_bw(base_size = 20)
}

pdf(file.path(plot_dir, 'modules_by_dx.pdf'))
print(plot_list)
dev.off()

################################################################################
#   Gene ontology for genes within diagnosis-associated modules
################################################################################

top_modules = tibble(
        module_num = as.integer(seq_len(length(p_val_list))),
        p_val = p.adjust(unlist(p_val_list), "fdr")
    ) |>
    filter(p_val < 0.05)

message("Significant modules by p-value for linear relationship with diagnosis:")
print(top_modules)

go_list = list()
univ = rowData(rse_gene)$EntrezID[!is.na(rowData(rse_gene)$EntrezID)]

for (module_num in top_modules$module_num) {
    genes = rowData(rse_gene)$EntrezID[net$colors == module_num]
    genes = genes[!is.na(genes)]

    go_list = append(
        go_list,
        enrichGO(
            genes, univ = univ, OrgDb = "org.Hs.eg.db", ont = "ALL",
            readable = TRUE, pvalueCutoff = 1, qvalueCutoff = 0.05
        )
    )
}

#   Which modules have any enriched GO terms?
which_enriched = which(sapply(go_list, function(x) nrow(x@result)) >= 1)

#   Create a 'compareClusterResult' object of modules with enriched GO terms
go_list = go_list[which_enriched]
names(go_list) = top_modules$module_num[which_enriched]
go_list = merge_result(go_list)

pdf(file.path(plot_dir, 'GO_module_enrichment.pdf'))
dotplot(go_list, showCategory = 20)
dev.off()

session_info()
