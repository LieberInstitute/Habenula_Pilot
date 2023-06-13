June 13th, 2023 - Bukola Ajanaku

Folder: 05_explore_sce/03_Plot_Investigate_snAnno
Files:
(1/9) 00_Info.R
(2/9) 01_snAnno_orig_Violin.R
(3/9) 02_sn_Anno_orig_Heatmap.R
(4/9) 03_New_Pseudobulk.R

Description:
This is a folder filled with scripts for plotting and corrected my annotations. Our finalized annotations were made in this folder.

1) 00_Info.R
A compiled list of more literature-based markers to confirm annotation identities.

2) 01_snAnno_orig_Violin.R
Violin plotted my original sce annotations with the original gene markers.

3) 02_sn_Anno_orig_Heatmap.R
Heatmapped my original pseudobulked sce annotations with the original gene markers.

4) 03_New_Pseudobulk.R
Relabeled LHb.6 as an endothelial cluster (Endo). Combined MHb.3 with MHb.2 using the mean expression plots. Peusdobulked and saved as "processed-data/04_snRNA-seq/sce_objects/sce_new_pseudobulk_with_snAnno_version2and3.Rdata".

5) 03b_New_Pseudobulk_Heatmaps.R
Using my new snAnno annotations, I heatmapped the pseudobulked data against the new markers list.

6) 04_newsnAnno.R
Ignore the repeats please. I used our new annotations to violin plot against an extensive list of marker genes.

7) 05_Updates_Annotations_meanExpression.R
I made and saved our new and finalized sn Annotation mean expression stats.

8) 06_Hockey_Stick_Plotter.R
I made hockey stick plots for the official annotations that indicated the top 25 marker genes.

9) 07_Plots_for_Progress_Report.R
Ignore the repeats again. I also added some renames to the clusters to shift them over by whatever was renamed or collapsed. I dropped the confusing Hb/Excit.Neuron Cluster and heatmapped our clusters based on literature and data-driven markers. Made for the progress report.