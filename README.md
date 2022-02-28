# Habenula_Bulk

RNA-seq data from homogenate Habenula

## JHPCE

Location: `/dcl02/lieber/ajaffe/Roche_Habenula`

##Analysis summary 
In the exploration of the habenula bulk RNA-seq data two primary analysis have been completed. The first question we asked is how does habenula bulk data differ from other regions. This was performed in the script de_analysis_gene_region_habenula.R and the results of the gene ontology enrichment analysis are in tables/geneSet_output.csv. tables/de_stats_allExprs_region.csv and tables/de_stats_fdr10_sorted_region.csv are the results of the differential expression analysis by region. A look at some of the metrics showed that mitochondrial mapping rate was confounded in the hippocampus and as a result it was added to the model. The next analysis was a case vs control differential expression analysis. This analysis is limited by the small sample size and will be preliminary results for future grants. This was performed in the script de_analysis_gene_habenula.R and the results of the gene ontology enrichment analysis are in tables/geneSet_output_dx.csv. tables/de_stats_allExprs_dx.csv and tables/de_stats_fdr10_sorted_dx.csv are the results of the differential expression analysis by region.

## Ongoing work
Currently work on deconvultion of the habenula bulk data. We have also discussed doing a power analysis on the scz vs control differential expression to better estimate how many samples are needed in the grant. 

## Figure Index
https://docs.google.com/presentation/d/1F7X1H-g7t4aWQrx2hUxQyrUcZLjR2zvQ6VvsHerX7Ng/edit#slide=id.p

* check_confounders.pdf- boxplot of potentially confounded variables across brain regions.  
* DE_boxplots_byDx.pdf- significant gene case control boxplots.  
* DE_boxplots_byRegion_clean.pdf- habenula vs other significant genes with cleaning effect boxplots.  
* DE_boxplots_byRegion.pdf- habenula vs other significant genes boxplots. 
* de_by_region_habenulaDown.pdf- boxplot with jitter by region for down expressed genes in region analysis. 
* de_by_region_habenulaUp.pdf-boxplot with jitter by region for up expressed genes in region analysis. 
* hist_pval_dx.pdf- distribution of pvalues by diagnosis.  
* hist_pval_regions.pdf - distribution pvalues by region. 
* PCA_plots_gene_Exprs_all_dx.pdf - pca scatter plot by diagnosis for all regions. 
* PCA_plots_gene_Exprs_all_regions.pdf - pca scatter by regions.  
* PCA_plots_gene_Exprs_control.pdf - pca scatter for only control brains. 
* PCA_plots_gene_Exprs_dx.pdf - pca scatter for case control. 
* variable_metrics_dx.pdf -boxplot of potentially confounded variables by diagnosis  
* volcano_plot_all.pdf - needs fixing (corrupted file). 
* volcano_plot_sig.pdf - volcano plot of sig genes in region analysis.  
* voom_gene_dx.pdf - voom plot for model by diagnosis.  
* voom_gene_region.pdf- vooom plot for model by region.
