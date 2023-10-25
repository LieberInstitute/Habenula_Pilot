# Habenula_Pilot

#### RNA-seq data from homogenate Habenula ####

### JHPCE

Location: `/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula`

### Analysis Summary

This project aims to answer two main questions: 1) What is the molecular taxonomy of the human habenula? and 2) Are there any differentially expressed genes between cases (SCZD) and control? We have been able to clean, cluster, annotate, and analyze the snRNA-seq data in order to produce the first ever molecular taxonomy of the human habenula. The results from this data have been leveraged to deconvolute our bulk RNA-seq samples (n = 69). We have dropped any brain samples with no Habenula, resulting in one drop (n = 68) and have decided to refer to our bulk RNA seq sample as Habenula-enriched samples as we perform our differential expression analyses. 

### Ongoing Work

We are currently working on the differential expression of genes between cases and controls as well as the trans-special analysis to determine if our Habenula markers are trans-specially conserved.

### Relevant Slides 
- **[Bulk QC](https://docs.google.com/presentation/d/1U4b2dCk3FI9uLYuBxrTyLwvCFAO2FxxPn5vkZCZyP2Y/edit?usp=sharing)**   
- **[Bulk PCA + snRNA QC + snRNA PCA](https://docs.google.com/presentation/d/1HWIBMhRQHeI8TPioMPTPJgeArAQ6AzHN3o4H_0l709E/edit?usp=sharing)**      
- **snRNA Full Analysis (2 Versions)** :
Clustering (and comparing my cluster schemes to Erik’s previous cluster schemes on this data set of which he had qc’ed differently), Annotations (the cluster method and resolution I settled on along with the markers I used to ultimately annotate each cluster), Mean Ratios (used to find the data driven gene markers), Summative Progress Report Heatmap
  - [Version 1](https://docs.google.com/presentation/d/1ua0hYzAk84n81v1r3w9_OqtHC_dxhN90xiaqlHtGRE8/edit?usp=sharing) - More biologically focused and shorter
  - [Version 2](https://docs.google.com/presentation/d/12kh6N5ALssipqBmgmU5y0pVn9z9AqzdZdskzhC6vWI0/edit?usp=sharing) - More computationally focused and longer


### Code Organization
Please note that each folder has an internal README.md file for clarity.

[01_bulk_speaqeasy](https://github.com/LieberInstitute/Habenula_Bulk/tree/master/code/01_bulk_speaqeasy) - bulk fastq files     
[07_cellranger](https://github.com/LieberInstitute/Habenula_Bulk/tree/master/code/07_cellranger) - sn transcriptomics data (only folder majorly out of order*)    
[02_bulk_qc](https://github.com/LieberInstitute/Habenula_Bulk/tree/master/code/02_bulk_qc) - bulk qc information     
[03_bulk_pca](https://github.com/LieberInstitute/Habenula_Bulk/tree/master/code/03_bulk_pca) - ran pca on filtered bulk data and investigated trends    
[04_snRNA-seq](https://github.com/LieberInstitute/Habenula_Bulk/tree/master/code/04_snRNA-seq) - full qc, pca, harmonization, clustering, and annotation journey  
[05_explore_sce](https://github.com/LieberInstitute/Habenula_Bulk/tree/master/code/05_explore_sce) - process of correcting annotations, finalizing identities of clusters, and collecting gene marker information.    
[06_deconvolution](https://github.com/LieberInstitute/Habenula_Bulk/tree/master/code/06_deconvolution) - process of bulk deconvolution  
[09_trans_special_analysis](https://github.com/LieberInstitute/Habenula_Bulk/tree/master/code/09_trans_special_analysis) - journey of performing trans-special analyses on snRNA habenula cluster data.   
[10_DEA](https://github.com/LieberInstitute/Habenula_Bulk/tree/master/code/10_DEA) - journey of differential expression analysis  
[70_GPR151](https://github.com/LieberInstitute/Habenula_Bulk/tree/master/code/70_GPR151) - side quest in checking GPR151 expression in our Habenula clusters  
[99_paper_figs](https://github.com/LieberInstitute/Habenula_Bulk/tree/master/code/99_paper_figs) - code used for all snRNA-seq paper figures including relevant bulk deconvolution plots 

<br />

> **Note on rse and sce objects**       
Official sce and rse objects used can be found in the processed-data folder on JHPCE. RSE objects has two renditions of official objects. The first official version is post-qc before dropping sampels and the second official version is post-dropping of sample based on deconvolution data (one sample had absolutely no Habenual).
> 




