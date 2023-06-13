January 13th, 2023 - Bukola Ajanaku

Folder: 09_trans_special_analysis
Files:
(1/3) 00_creating_Wallace.R

Description:
This was supposed to be a transpecial analysis between our human Habenula data and the mouse Habenula data from the Wallace et al. 2019 paper.

1) 00_creating_Wallace.R
Loaded up the Wallace data (GSE146983). Tried my best to piece together the information posted on 
the site https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/2VFWF6 and combine for later analysis. Later found that no cell type information was shared.

2) 01_Wallace_sce_homolog_creator.R
Created sce objects between our human data and the mice data, subsetting for only homologous genes.

3) 02_Wallace_sce_calculations.R
Playing around with different locations in search of the cell-type data. None to be found. Cell-types needed in order to run spatialLIBD analysis.