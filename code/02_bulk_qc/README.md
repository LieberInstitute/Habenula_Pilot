June 12th, 2022 - Bukola Ajanaku

Folder: 02_bulk_qc
Files: 
  (1/3) 01_qc_data_Bukola.R
  (2/3) 02_build_objects_Bukola.R
  (3/3) 03_qc_plotter_Bukola.R

Description:

This folder is the bulk RNA-seq quality control aspect of the project.

(1/3) 01_qc_data_Bukola.R:
Working on all (gene, exon, junction, and transcipt) forms of our bulk RNA-seq data. Added phenotypic information to our rse object and boxplotted to visualize said phenotypes against covariates of interest. We had an issue of mislabeling for some samples (this was established from before my time with LIBD) so used a paralleled naming convention to correct those we would correct while drop whatever sample(s) we could not correctly identify. Then saved. Purpose of this script was to identify if there were any metrics that were driving trends in our bulk data. Couldn't find any.

(2/3) 02_build_objects_Bukola.R **
Here we are building the official rse objects by computing the qc metrics for our post-brain swap rse object. I made sure no NAs were kept and normalized our data for the porportions of 0s in each data type. I used addPerCellQC() to create the qc stats, filtered out genes with overall low expression values, added ensemble IDs to the rse object and saved! No bulk samples were dropped in this step.

ALL QUALITY CONTROLLED RSE OBJECTS ARE IN "~processed-data/02_bulk_qc/count_data_bukola/rse_tx_filt_Roche_Habenula_qcAndAnnotated_n69.Rdata" 

(3/3) 03_qc_plotter_Bukola.R
Made a series of plots to further investigate if any metrics where driving trends in the rse data post qc and filtration.





