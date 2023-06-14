**June 12th, 2023**     
*Bukola Ajanaku* 

#### 04_snRNA-seq ####

|   |       Files*     |
| --- | --------------- |
| 1 | 00_Color_Scheme_ct.R |
| 2 | 00_Gene_Marker_List.R | 
| 3 | 01_get_droplet_scores.R |
| 4 | 02_build_basic_sce.R | 
| 5 | 03_qc_metter.R | 
| 6 | 04_GLM_PCA.R | 
| 7 | 05_GLM_Harmony.R | 
| 8 | 06_Clustering.R | 
| 9 | 07_Gene_Marking.R |
| 10 | 07b_Marking_Clusters.R | 
| 11 | 07c_Pseudobulking.R | 
| 12 | 07d_Pseudo_Heatmap.R | 
| 13 |  08_Hierarchical_Clustering.R | 
| 14 | 09_Clustered_QC.R | 

**\* only listing main R files**


#### Description ####

This is the full snRNA-seq analysis journey that I went through during most of my time with LIBD. This was my take on the snRNA-seq data, similar to how Josh Stolz and Erik Nelson have done (in the Legacy folder). Please note that a great deal of this work was exploratory. All finialized pathways are in the 99_paper_figs folder.

- **00_Color_Scheme_ct.R**     
Unofficial unspecified color palette used for exploratory work.

- **00_Gene_Marker_List.R**    
Listing out the categories of cell types and their associated gene markers that are used to annotate the clusters of nuclei that were later made.

- **01_get_droplet_scores.R**      
Calculated droplet scores per cellranger (7th folder, only folder majorly out of order due to folder organization) per Sample. Saved output as RDATA files (in "processed-data/04_snRNA-seq/01_get_droplet_scores"). Running on cluster created .sh files and log folders. Also plotted droplet score plots with elbows for false discovery rate drop threshold in "plots/04_snRNA-seq/01_get_droplet_scores". 

- **02_build_basic_sce.R**     
Created sce object of all post-empty droplet sce objects of each sample (7 separate snRNA-seq Samples combined into main sce object). This is our pre-QC but post-emptyDroplets sce object saved as "processed-data/04_snRNA-seq/sce_objects/sce_hb_preQC.Rdata".

- **03_qc_metter.R**       
Using only the gene data for the rse object, I calculated mitochondrial percent, library size, and detected features per Sample. I dropped high mito_percent (using 3 mads approach), and low library size and detected features. Doublet scores were also calculated used the top 1000 HVGs and plotted. No Sample had a median of 5 or more doublets. 
The rest of this code is mainly leveraging the cell-type annotaion work creted by Erik Nelson. I call his annotated work "annoData". Because he cleaned the sce object slightly differently from how I cleaned it, there re nuclei involved in my sce object that are not found in his and vice versa. I plot this against different qc metrics for our dropped nuclei (violin plots of drops by mito_percent, lib size, detected features) and all of these plots can be found in "plots/04_snRNA-seq/03_qc_metter_plots".

- **04_GLM_PCA.R**
Loading post-qc sce object ("processed-data/04_snRNA-seq/sce_objects//sce_hb_postQC.Rdata"), I add lab information regarding run and then perform feature selection to later grab the top 2000 most variaible genes and run GLM PCA. I plot my GLM data by different metrics/groupings (Erik's annotations (annoData), Run, Sample, mito percent, and etc). I then use these similar metrics to plot with for my TSNEs and UMAPs. All plots cann be found in "plots/04_snRNA-seq/04_GLM_PCA_plots".  
**PLEASE NOTE!!!**  The post-GLM sce data was saved as 
"processed-data/04_snRNA-seq/sce_objects/sce_uncorrected_PCA.Rdata". This object was not harmonized. I later harmonize this data by sample. 


- **05_GLM_Harmony.R ** 
Corrected sce by Harmony and plotted by different metrics (similar to those plotted in 04_GLM_PCA.R) and was saved in "plots/04_snRNA-seq/05_GLM_Harmony_plots". Post-harmony GLM sce was saved as "processed-data/04_snRNA-seq/sce_objects/sce_harmony_by_Samp.Rdata"). 

- **06_Clustering.R**
This was exploratory code looking at different clustering methods for our data. Louvain was the example established by Stephanie Hicks btw whereas WalkTrap is one noted in the OSCA handbook and was done with our data set before through Erik Nelson. I compared my various clusters with Erik's annoData because he did a pretty good job overall with clustering and annotations. Louvain was not as clean as WalkTrap (also, plotting Louvain vs WalkTrap on TSNEs and comparing to TSNEs with Erik's annotations shows that Louvain focuses on splitting between thalamic neurons too much and fails to tease out difference amongst different cell types. Wasn't the resolution we needed). However, we have left room for improvement by using different nearest neighbors for each clustering method. All plots were kept in "plots/04_snRNA-seq/06_Clustering". The sce object for this code was kept in "processed-data/04_snRNA-seq/sce_objects/sce_post_clustering.Rdata".

- **07_Gene_Marking.R**
At some point of coding non-linearly, I ended up saving the post-clustered sce object with logcounts. I didn't need to make a new sce save but I didn't know that at the time, I'm glad I know that now. Using code developed by Louise Huuki-Myers (now part of her DeconvoBuddies package), I created violin plots for all three WalkTrap nearest neighbors that I selected (10, 20, 50), against our literature-based gene marker list. Also, note that columns in my sce's colData labeled "wT_20_Erik" for instance is the WalkTrap 20 method and I subscripted with Erik because it was the clustering method that Erik used not because they are Erik's notations. Remember, I was runnning many different clustering methods throughout this exploratory work. All plots were saved in "plots/04_snRNA-seq/07_Gene_Marking".

- **07b_Marking_Clusters.R** 
I annotated all clusters in Excel. I loaded them up inn this script and added them to the colData of my sce object. Please ignore the fastHplus lines of code. We ended up not using it at all.

- **07c_Pseudobulking.R**
Simple pseudobulking by sample. This is in preparation for the complex heatmap in later steps! Ran on cluster, creating .sh file and log. This sce was saved as "processed-data/04_snRNA-seq/sce_objects/sce_pseduobulked_wT.Rdata".

- **07d_Pseudo_Heatmap.R**  
Grabbed pseudo_bulked sce object and adding my annotated information by cluster number. I then create my Complex Heatmap between my clusters and my manual annotations. Some cluster identities were refined based on the heatmaps results. All plots were saved in "plots/04_snRNA-seq/07d_Pseudo_Heatmaps".

- **08_Hierarchical_Clustering.R** 
We tried hierarchical clustering but chose against it because of a host of reasons. The plots ended up in "plots/04_snRNA-seq/08_Hierarchical_Clustering".

- **09_Clustered_QC.R** 
Now that we've settled on walkTrap 10 as our clustering method and have felt settled with the identities of these clusters, this code plots the QC metrics of our clusters. I've also taken into account highlighting problematic clusters (cluster's we weren't fully sure of) as well as looking at the clusters at collapsed and granular resolutions. All plots were made in "plots/04_snRNA-seq/09_Clustered_QC". The sce object with all additional columns highlighting the problematic clusters and the collapsed resolution was saved as "processed-data/04_snRNA-seq/sce_objects/sce_post_09_clustered_qc.Rdata".
