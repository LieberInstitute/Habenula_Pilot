June 13, 2023 - Bukola Ajanaku

Folder: 05_explore_sce
Files:
(1/3) 01_Manually_Anno_sce.R

Description:
Now that we have our final sce object sofar ("processed-data/04_snRNA-seq/sce_objects/sce_post_09_clustered_qc.Rdata") and have decided to move forward with the WalkTrap 10, here we are exploring the annotated clusters to futher confirm their identities, find human Habenula gene markers, and create probes that will help deconvolute our bulk data.

1) 01_Manually_Anno_sce.R
Listing out the sce object that will be used in this folder as well as the relevant annotations and where they are kept. For our official sn annotations (snAnno), we collapsed everything but the LHb and MHb clusters. We also had a more arbitrary cluster that we called Hb but we aren't entirely sure if its L or M. Or if it was for sure Habenula actually. This cluster is later named Excit.Neuron and is dropped as it was an ambiguous cluster.

2) 02_Mean_Ratio_Explore.R
Using the DeconvoBuddies package by Louise Huuki-Myers, I was able to plot the top marker genes for the ultra fine resolution of clusters as well as our official snAnno resolution. I decided to investigate snAnno further rather than the split clusters. Marker stats for snAnno were housed here "processed-data/05_explore_sce/mean_ratio_data_from_02_mean_ratio_explore_file.Rdata". Not really important as we later change some things for the cluster.

3) FOLDER: 03_Plot_Investigate_snAnno
Created a folder for clean up and to focus on snAnno exploration. Made many changes and created final Annotations in this plot. Also, plot for progress report heatmap is found in this folder.