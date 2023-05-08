## May 2, 2023 - Bukola Ajanaku
# Plotting cell-type expression pre and post drop per sample 
# qrsh -l mem_free=50G,h_vmem=50G

# loading relevant libraries
library(SummarizedExperiment)
library(here)
library(SingleCellExperiment)
library(DeconvoBuddies)
library(tidyverse)
library(tibble)

# loading sce object with dropped ambig cluster
load(file = here("processed-data", "99_paper_figs", "sce_objects", 
                 "sce_final_preHbdrop.RDATA"))
sce <- sce_final_preHbdrop


# creating plot_dir
plot_dir <- here("plots", "99_paper_figs", "08_sce_Plot_Expression")
if(!dir.exists(plot_dir)){
  dir.create(plot_dir)
}

# changing OPC_noisy class into a cluster of it's own.
  # adding rownames of colData as a row for easier subsetting
OPC_noisy_Samps = c("Br5555", "Br1204", "Br1092")

# grabbing barcodes for  OPC
onlyOPC <- sce[, which(sce$final_Annotations == "OPC")]

# grabbing noisy OPC
onlyOPC <- onlyOPC[, which(onlyOPC$Barcode %in% OPC_noisy_Samps)]

# grabbing barcodes            
RowNos <- rownames(colData(onlyOPC))

# adding Nos to RowNos 
sce[, rownames(sce) %in% RowNos]$final_Annotations <- "OPC_noisy"

# check
table(sce$Sample, sce$final_Annotations)







# dropping OPC_noisy
table(sce$OPC_clean)
    # No   Yes 
    # 594 16437 
sce <- sce[, which(sce$OPC_clean == "Yes")]
# check again
table(sce$OPC_clean)
    # Yes 
    # 16437 

# grabbing relevant columns
pd <- as_tibble(colData(sce)[,c("Sample", 
                                    "final_Annotations", "NeuN", "Run")])
pd$total_nuclei <- NA
pd$ct_nuclei <- NA

# creating function that collects nuclei data
for(i in unique(pd$Sample)){

  pd[pd$Sample == i, ]$total_nuclei <- nrow(pd[pd$Sample == i, ])

  for(p in unique(pd$final_Annotations)){
    
    if(nrow(pd[pd$Sample == i & pd$final_Annotations == p, ]) == 0){
      tot_ct = 0
    } else{
      tot_ct = nrow(pd[pd$Sample == i & pd$final_Annotations == p, ])
      pd[pd$Sample == i & pd$final_Annotations == p, ]$ct_nuclei <- tot_ct
    }
  }
}

pd$prop <- ( pd$ct_nuclei / pd$total_nuclei )

# testing to make sure it adds up to %100
  # test <- pd[pd$Sample == i,]
  # tester <- unique(test$prop)
  # sum(tester)
   # [1] 99.90193


# creating sce composition plot
pdf(file = here(plot_dir, "sce_Comp_Expression.pdf"))
  
ggplot(pd, aes(fill = final_Annotations, y = prop, x = Sample)) +
    geom_col(stat = "identity")

dev.off()

## create composition bar plots (using plot_composition_bar)
pdf(here(plot_dir, "sce_Comp_Express_Bar.pdf"), width = 21, height = 12)

plot_composition_bar(prop_long = pd, sample_col = "Sample",
                     x_col = "Sample", ct_col = "final_Annotations")

dev.off()



# testing summary 

summary(pd[pd$Sample == i,])


#  