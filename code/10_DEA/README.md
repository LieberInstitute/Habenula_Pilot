Original by *Renee Garcia-Flores*

#### 10_DEA ####

|   |       Files     |
|---| --------------- |
| 1 | 02_DataExploration.R |
| 2 | 03_VariancePartition.R | 
| 3 | 04_DEA.R | 
| 4 | 05_DEAplotting.R |
| 5 | 06_GOenrichment.R |
| 6 | 07_GOplot.R |
| 7 | 08_Comparisons.R |

File 01_build_final_rse.R was move to code/02_bulk_qc as 05_build_final_rse.R

#### Description ####

This folder has all the files for DEA

- **02_DataExploration.R**       
This script has the plots to explore our rse object and select the variables for the DEA analysis. The outputs are boxplots (sample_variables - qc_metrics), correlation plots (AgeDeath - qc_metrics) and a heatmap that shows the correlation between PCs of the logcounts and all our metrics.

- **03_VariancePartition.R**        
This has the variance partition analysis, including getVarianceExplained() and canCorPairs() to explore correlation between variables and the variance explained by them.

- **04_DEA.R**        
Actual DEA with limma-voom pipeline. The formula we are using is *~ PrimaryDx + AgeDeath + Flowcell + mitoRate + rRNA_rate + RIN + totalAssignedGene + abs_ERCCsumLogErr + qSV1 + qSV2 + qSV3 + qSV4 + qSV5 + qSV6 + qSV7 + qSV8 + tot.Hb + tot.Thal*.
This includes DEa for genes, exons, jx and tx.

- **05_DEAplotting.R**  
Script to plot volcano plots for genes, exons, jx and tx.

- **06_GOenrichment.R**  
GO and KEGG enrichment using compareCluster() for genes, exons, jx and tx. 

- **07_GOplot.R**  
Plots GO and KEGG results. These plots are just for exploration, they are not publication ready

- **08_Comparisons.R**  
Script to obtain the EnsemblIDs from the differentially expressed genes, exons, jx and tx, as well as the multiple intersectiosn between each DEA


