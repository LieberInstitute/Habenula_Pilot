Original by *Renee Garcia-Flores*

#### 12_GSEA ####

|   |       Files     |
|---| --------------- |
| 1 | 01_GeneSetErichment.R |
| 2 | 02_GSE-plot.R | 
| 3 | 03_ExtractGenes.R | 
| 4 | gene_set_enrichment_plot_complex.R |
| 5 | marker_gene_set_enrichment.R |

#### Description ####

This folder has all the files to do Gene Set Enrichment Analysis

- **01_GeneSetErichment.R**       
GSEA with single cell data (three sets: genes for broad annotation, final annotation and top 25 gene markers) and differentially expressed genes from bulk data. Analysis made with gene_set_enrichment() from *spatialLIBD*.
For the top 25 gene markers' GSEA, functions from marker_gene_set_enrichment.R are needed.

- **02_GSE-plot.R**        
Script to plot heatmap with results. This script uses functions from gene_set_enrichment_plot_complex.R 
      
- **03_ExtractGenes.R**        
Script to obtain the list of genes that intersect in both data sets. 
