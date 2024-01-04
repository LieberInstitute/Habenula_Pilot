
library("SingleCellExperiment")
library("here")
library("tidyverse")

#### Load SCE ####
load(here("processed-data", "sce_objects", "sce_Habenula_Pilot.Rdata"), verbose = TRUE)

dim(sce)
# [1] 33848 16437

summary(as.numeric(table(sce$Sample)))
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 778    1604    2885    2348    3080    3405 

## Total UMIs
summary(sce$sum)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 605    8671   14995   18124   23624  170016

## non zero-expression per nuc
zero_sums <- colSums(counts(sce) != 0)
summary(zero_sums)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 481    3262    4734    4787    6054   13178

table(sce$Sample)

pd <- as.data.frame(colData(sce)) |> as_tibble()

pd |>
  group_by(Sample) |> 
  summarize(mean_umi = mean(sum),
            median_umi = median(sum)) |>
  summary()

# Sample             mean_umi       median_umi   
# Length:7           Min.   : 9419   Min.   : 8382  
# Class :character   1st Qu.:13774   1st Qu.:11842  
# Mode  :character   Median :18131   Median :16453  
#                    Mean   :22207   Mean   :20200  
#                    3rd Qu.:29240   3rd Qu.:27548  
#                    Max.   :41869   Max.   :37784  




