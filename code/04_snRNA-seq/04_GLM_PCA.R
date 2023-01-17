# January 17, 2023 - Bukola Ajanaku
# Creating pca plots for filtered (removed empty droplets), quality controlled 
# (dropped high mito, low library size, low detected features, and any genes 
# with 0 counts across samples) sce object.
# Based on:
# 1) https://github.com/LieberInstitute/DLPFC_snRNAseq/blob/main/code/03_build_sce/04_GLM_PCA.R
# 2) https://www.stephaniehicks.com/biocdemo/articles/Demo.html
# qrsh -l mem_free=50G,h_vmem=50G

library("SingleCellExperiment")
library("scran")
library("scater")
library("scry")
library("here")
library("sessioninfo")
library("mbkmeans")

# loading post qc sce object 
load(here("processed-data", "04_snRNA-seq", "sce_objects", 
           "sce_hb_postQC.Rdata"))
sce <- sce_hb_postQC
rm(sce_hb_postQC)

# Normalization?
mbk <- mbkmeans(sce, whichAssay = "counts", reduceMethod = NA,
                clusters=10, batch_size = 500)
sce$mbk10 <- paste0("mbk", mbk$Clusters)
table(mbk$Clusters)

# Deviance featuring selection
sce <- devianceFeatureSelection(sce,
        assay = "counts", fam = "binomial", sorted = F,
        batch = as.factor(sce$Sample))

# Checking outputs for deviance ft selection function
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "binomial_deviance.pdf"))
  plot(sort(rowData(sce)$binomial_deviance, decreasing = T),
       type = "l", xlab = "ranked genes",
       ylab = "binomial deviance", main = "Feature Selection with Deviance"
  )
  abline(v = 2000, lty = 2, col = "red")
dev.off()

# Taking the GLM_PCA approach
sce <- nullResiduals(sce,
                     assay = "counts", fam = "binomial", # default params
                     type = "deviance")

# Running PCA 
sce <- scater::runPCA(sce, ncomponents = 50,
                      exprs_values = "binomial_deviance_residuals",
                      scale = TRUE, name = "GLMPCA_approx",
                      BSPARAM = BiocSingular::RandomParam())

# Plotting using GLM-PCA functions
pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "GLM_pca_plot.pdf"))
  plotReducedDim(sce, dimred = "GLMPCA_approx", colour_by = "Sample")
  plotReducedDim(sce, dimred = "GLMPCA_approx", colour_by = "mbk10")
dev.off()

# Normal PCA
# had to re-load sce 
sce <- computeSumFactors(sce, cluster = mbk$Clusters, min.mean = 0.1)
sce <- logNormCounts(sce)

sce <- scater::runPCA(sce, ncomponents = 50,
                      ntop = 1000,
                      scale = TRUE,
                      BSPARAM = BiocSingular::RandomParam())

sce$mbk10 <- paste0("mbk", mbk$Clusters)

pdf(here("plots", "04_snRNA-seq", "04_GLM_PCA_plots", "normal_pca_plot.pdf"))
  plotPCA(sce, colour_by = "Sample")
  plotPCA(sce, colour_by = "mbk10")
dev.off()



# Plotting PCA (personal function)
# pc_to_pc <- function (pcx, pcy, pc_df, colorbylist, dataType) {
#   # pcx and pcy are the pcas of interest
#   # pc_df is the output object from the PCA creator function
#   # colorbylist is the list of metrics to color plots buy (changes number of 
#   # plots per saved file)
#   # dataType options are: "exon" "tx" "jx" "gene"
#   # sampsdropped is one character string of brains that we dropped 
#   # (numdrop and sampsdropped are either both NA OR both not)
#   
#   # unlisting returned object from pca_creator function
#   pc_data = pc_df[[1]]
#   pc_variables = pc_df[[2]] 
#   
#   # Use pos to ensure jitter and text_repel share coordinates (prevents mislabeling).  
#   pos <- position_jitter(seed = 2)
#   
#   # prepping for saving plots into list
#   plot_list = list()
#   c = 1
#   
#   # plotting by coloring scheme
#   for(i in colorbylist){
#     
#     forhighlight = pc_data[pc_data$BrNum %in% c("Br1676","Br5459", "Br6323"),]
#     nothighlight = pc_data[! pc_data$BrNum %in% c("Br1676","Br5459", "Br6323"),]
#     
#     if (i == "PrimaryDx"){
#       plot = ggplot(pc_data, aes_string(x = pcx, y = pcy)) +
#         scale_color_manual(values=wes_palette(n=2, name="GrandBudapest1"))+
#         geom_jitter(aes_string(colour = i), data = nothighlight, alpha = 0.8, position = pos, size = 6) +
#         geom_jitter(aes_string(x = pcx, y = pcy), colour = "darkblue", 
#                     data = forhighlight, alpha = 0.8, position = pos, size = 6) +
#         geom_text_repel(aes_string(label = "BrNum"), position = pos, 
#                         color = "lightgrey", max.overlaps = 5) +
#         labs(x = x_titler, y = y_titler,
#              color = rename_vars[rename_vars$orig_var_name ==  i,]$var_plot_title,
#              caption = "Br1676, Br5459, and Br6323") +
#         theme_bw(base_size = 10) + 
#         theme(legend.position= "bottom", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
#               axis.text.x = element_text(vjust = 0.7), text = element_text(size=15),
#               axis.title = element_text(size=15),
#               plot.caption = element_text(color = "darkblue", face = "italic")) 
#       
#     } else if (i == "Flowcell"){
#       
#       plot = ggplot(pc_data, aes_string(x = pcx, y = pcy)) +
#         scale_color_manual(values=wes_palette(n=2, name="IsleofDogs1"))+
#         geom_jitter(aes_string(colour = i), data = nothighlight, alpha = 0.8, position = pos, size = 6) +
#         geom_jitter(aes_string(x = pcx, y = pcy), colour = "darkblue", 
#                     data = forhighlight, alpha = 0.8, position = pos, size = 6) +
#         geom_text_repel(aes_string(label = "BrNum"), position = pos, 
#                         color = "lightgrey", max.overlaps = 5) +
#         labs(x = x_titler, y = y_titler,
#              color = rename_vars[rename_vars$orig_var_name ==  i,]$var_plot_title,
#              caption = "Br1676, Br5459, and Br6323") +
#         theme_bw(base_size = 10) + 
#         theme(legend.position= "bottom", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
#               axis.text.x = element_text(vjust = 0.7), text = element_text(size=15),
#               axis.title = element_text(size=15),
#               plot.caption = element_text(color = "darkblue", face = "italic")) 
#       
#     } else {
#       
#       plot = ggplot(pc_data, aes_string(x = pcx, y = pcy)) +
#         scale_color_gradientn(colors = wes_palette("Darjeeling1", type = "continuous")) +
#         geom_jitter(aes_string(colour = i), data = nothighlight, alpha = 0.7,
#                     position = pos, size = 5) +
#         geom_jitter(aes_string(x = pcx, y = pcy), colour = "darkblue", 
#                     data = forhighlight, alpha = 0.8, position = pos, size = 6) +
#         geom_text_repel(aes_string(label = "BrNum"), position = pos, 
#                         color = "lightgrey", max.overlaps = 6) +
#         labs(x = x_titler, y = y_titler,
#              color = rename_vars[rename_vars$orig_var_name ==  i,]$var_plot_title,
#              caption = "Br1676, Br5459, and Br6323") +
#         theme_bw(base_size = 10) + 
#         theme(legend.position= "bottom", plot.margin=unit (c (1.5,2,1,2), 'cm'), 
#               axis.text.x = element_text(vjust = 0.7), text = element_text(size=15),
#               axis.title = element_text(size=15),
#               plot.caption = element_text(color = "darkblue", face = "italic")) 
#       
#     }
#     
#     plot_list[[c]] = plot 
#     c = c + 1
#   }
#   # Saving plots is easier in the function
# 
#     
#     pdf(file = here("plots", "03_bulk_pca", "pc_plots_bukola", type, orgbybrain, nameFILE))
#     print(plot_list)
#     dev.off()
#   }
#   



















## 