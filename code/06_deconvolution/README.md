June 13th, 2023 - Bukola Ajanaku

Folder: 06_deconvolution
Files:

1) 01_Bulk_Annotations.R
Created different bulk resolutions of our finalized sn annotation resolution levels. Saved as sce_first_bulkTypes.RDATA in "processed-data/06_deconvolution/sce_objects".

2) 02a_bulkTypeSepHb_Exp_Data.R, 
3) 02b_bulkTypeBroadHb_Exp_Data.R, and
4) 02c_bulkTypeAllCollapse_Exp_Data.R
Ran marker stats code, made hockeystick plots, and saved R objects for the different types of bulk resolutions. We liked the bulkTypeSepHb version the best. PLEASE NOTE, MARKER STATS FOR THESE THREE RESOLUTIONS OF OUR DATA (VERY USEFUL FOR HEATMAPPING), ARE IN THEIR SELF-TITLED AND RESPECTIVE PLOTTING FOLDERS!

5) 03_run_Bisque.R
Another example of nonlinearity. Wherever I dropped the Hb ambig clusters, I dropped the dirty OPC samples that were making our large OPC cluster noisy. The clean OPC should only have 1202 nuclei. Here, I loaded the filtered rse object and final sce object with the bulkTypeSepHb. Grabbed marker starts, plotted marker expression, ran bisque and saved. 

The post OPC-cleaning sce marker stats are here:
"processed-data/06_deconvolution/run_Bisque/marker_stats_top_25_genes.csv"

The post OPC-cleaning bulk deconvolution cell-type proportion stats are here:
"processed-data/06_deconvolution/run_Bisque/est_prop_split_Hb_annotations.RDATA"

6) 04_explore_Bisque.R
Plotting composition results from the bulk deconvolution. Not much to see here unless your playing around with aesthetics! 
