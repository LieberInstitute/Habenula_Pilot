* Analysis summary In the exploration of the habenula bulk RNA-seq data two primary analysis have been completed. 
* The first question we asked is how does habenula bulk data differ from other regions. 
* This was performed in the script `de_analysis_gene_region_habenula.R` and the results of the gene ontology enrichment analysis are in `tables/geneSet_output.csv`. 
* `tables/de_stats_allExprs_region.csv` and `tables/de_stats_fdr10_sorted_region.csv` are the results of the differential expression analysis by region.
* A look at some of the metrics showed that mitochondrial mapping rate was confounded in the hippocampus and as a result it was added to the model. 
* The next analysis was a case vs control differential expression analysis.
* This analysis is limited by the small sample size and will be preliminary results for future grants. 
* This was performed in the script `de_analysis_gene_habenula.R` and the results of the gene ontology enrichment analysis are in `tables/geneSet_output_dx.csv`.
* `tables/de_stats_allExprs_dx.csv` and `tables/de_stats_fdr10_sorted_dx.csv` are the results of the differential expression analysis by region.