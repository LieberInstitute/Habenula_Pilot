library(here)
library(tidyverse)
library(SummarizedExperiment)
library(jaffelab)
library(clusterProfiler)
library(rrvgo)
library(sessioninfo)

net_path = here('processed-data', '19_wgcna', 'modules.rds')
rse_path = here(
    'processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda'
)
plot_dir = here('plots', '19_wgcna')
protected_deg_covariates = c('PrimaryDx', 'AgeDeath')
deg_covariates = c(
    'PrimaryDx', 'AgeDeath', 'Flowcell', 'mitoRate', 'rRNA_rate', 'RIN',
    'totalAssignedGene', 'abs_ERCCsumLogErr', 'qSV1', 'qSV2', 'qSV3', 'qSV4',
    'qSV5', 'qSV6', 'qSV7', 'qSV8', 'tot.Hb', 'tot.Thal'
)
fdr_cutoff = 0.1

set.seed(0)
dir.create(plot_dir, showWarnings = FALSE)

################################################################################
#   Preprocess net and expression data
################################################################################

#   Load WGCNA network and gene expression data, then filter genes by minimum
#   mean RPKM
net = readRDS(net_path)
rse_gene = get(load(rse_path))
assays(rse_gene)$rpkm = recount::getRPKM(rse_gene, length_var = 'Length')
rse_gene = rse_gene[rowMeans(assays(rse_gene)$rpkm) > 0.25,]

mod_deg <- model.matrix(
    as.formula(paste('~', paste(deg_covariates, collapse = " + "))),
    data = colData(rse_gene)
)

#   Log transform expression, regress out covariates, and transpose
exp_mat = log2(assays(rse_gene)$rpkm + 1) |>
    cleaningY(mod_deg, P = 3) |>
    t()

################################################################################
#   Plot module weights by diagnosis (not manuscript ready)
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
    filter(p_val < fdr_cutoff)

message("Significant modules by p-value for linear relationship with diagnosis:")
print(top_modules)

univ = as.character(
    rowData(rse_gene)$EntrezID[!is.na(rowData(rse_gene)$EntrezID)]
)
for (ont_type in c("BP", "MF", "CC")) {
    go_list = list()

    for (module_num in top_modules$module_num) {
        genes = rowData(rse_gene)$EntrezID[net$colors == module_num]
        genes = genes[!is.na(genes)]

        go_list = append(
            go_list,
            enrichGO(
                genes, univ = univ, OrgDb = "org.Hs.eg.db", ont = ont_type,
                readable = TRUE, pvalueCutoff = 1, qvalueCutoff = fdr_cutoff
            )
        )
    }

    #   Which modules have any enriched GO terms?
    which_enriched = which(sapply(go_list, function(x) nrow(x@result)) >= 1)
    if (ont_type == "MF") {
        sig_modules = top_modules$module_num[which_enriched]
    }

    #   Create a 'compareClusterResult' object of modules with enriched GO terms
    go_list = go_list[which_enriched]
    names(go_list) = top_modules$module_num[which_enriched]
    go_list = merge_result(go_list)

    #   Dot plots
    pdf(
        file.path(
            plot_dir, sprintf('GO_module_enrichment_%s.pdf', ont_type)
        ),
        height = 14
    )
    print(dotplot(go_list, showCategory = 20))
    dev.off()
    
    #   Tree maps using semantic similarity of GO terms
    go_df = as.data.frame(go_list) |>
        as_tibble()
    
    sim_matrix = calculateSimMatrix(
        go_df$ID,
        orgdb = "org.Hs.eg.db",
        ont = ont_type,
        method = "Rel"
    )

    scores = -log10(go_df$qvalue)
    names(scores) = go_df$ID
    reduced_terms = reduceSimMatrix(
        sim_matrix,
        scores,
        threshold = 0.7,
        orgdb="org.Hs.eg.db"
    )

    pdf(file.path(plot_dir, sprintf('GO_treemap_%s.pdf', ont_type)))
    print(treemapPlot(reduced_terms))
    dev.off()
}

################################################################################
#   Plot module weights by diagnosis (manuscript ready)
################################################################################

#   For modules associated with both diagnosis and GO MF results, plot
#   manuscript-ready boxplots
label_df = tibble(
        module_num = factor(
            paste("Module", sig_modules),
            levels = paste("Module", sort(sig_modules))
        ),
        label = unlist(p_val_list[sig_modules]),
    ) |>
    mutate(label = paste('\np =', signif(label, 2), ''))

p = me_df |>
    select(RNum, PrimaryDx, all_of(paste0("ME", sig_modules))) |>
    pivot_longer(
        all_of(paste0("ME", sig_modules)), 
        names_to = "module_num", values_to = "weight"
    ) |>
    mutate(
        module_num = factor(
            str_replace(module_num, '^ME', 'Module '),
            levels = paste("Module", sort(sig_modules))
        )
    ) |>
    ggplot(aes(x = PrimaryDx, y = weight, color = PrimaryDx)) +
        geom_boxplot(outlier.shape = NA) +
        facet_wrap(~module_num) +
        geom_jitter() +
        geom_text(
            data = label_df,
            aes(label = label),
            x = Inf, y = Inf, hjust = 1, vjust = 0.5, color = 'black', size = 7
        ) +
        guides(color = "none") +
        labs(x = "Diagnosis", y = "Weight") +
        theme_bw(base_size = 20)
pdf(file.path(plot_dir, 'clean_boxplots_MF_modules.pdf'), height = 10)
print(p)
dev.off()

################################################################################
#   Heatmap of top-correlated genes with select module eigengenes
################################################################################

#   Correlate select MEs with gene expression to calculate "module membership".
#   Use gene symbol and format for plotting later
cor_df = cor(net$MEs[, paste0("ME", sig_modules)], exp_mat) |>
    abs() |>
    t() |>
    as.data.frame() |>
    rownames_to_column('gene_id') |>
    as_tibble() |>
    pivot_longer(
        all_of(paste0("ME", sig_modules)),
        names_to = "module", values_to = "cor_val"
    ) |>
    mutate(
        gene_id = rowData(rse_gene[gene_id,])$Symbol,
        module = factor(
            str_replace(module, '^ME', ''),
            levels = as.character(sort(sig_modules))
        )
    )

#   Grab genes with highest module membership for each module
top_genes = cor_df |>
    group_by(module) |>
    arrange(desc(cor_val)) |>
    slice_head(n = 5) |>
    ungroup()

#   Plot heatmap
p = cor_df |>
    filter(gene_id %in% top_genes$gene_id) |>
    mutate(
        gene_id = factor(
            gene_id,
            levels = unlist(
                lapply(
                    levels(module),
                    function(x) {
                        top_genes |>
                            filter(module == x) |>
                            pull(gene_id)
                    }
                )
            )
        )
    ) |>
    ggplot(aes(x = gene_id, y = module, fill = cor_val)) +
        geom_tile() +
        scale_fill_viridis_c() +
        theme_bw(base_size = 25) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + 
        labs(x = "Gene Symbol", y = "Module", fill = "Module\nMembership")
pdf(file.path(plot_dir, 'top_genes_heatmap.pdf'), height = 4.5, width = 10)
print(p)
dev.off()

session_info()
