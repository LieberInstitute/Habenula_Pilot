
library("tidyverse")
# library("readxl")
library("GGally")
library("here")
library("sessioninfo")

## set up 
plot_dir <- here("plots", "14_RNAscope", "02_MHb_HALO")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load data
fn <- list(job1862 = here("processed-data", "14_RNAscope", "HALO_data", 
                          "MHbExperiment_Br5422_middleslice_20x_multistitchimage_maxIP_unmixed.nd2_job1862MHbExperiment_Br5422_middleslice_20x_multistitchimage_maxIP_unmixed_second_pass_object_results_FOV.csv"),
           job1864 = here("processed-data", "14_RNAscope", "HALO_data", 
                          "MHbExperiment_Br5422_middleslice_20x_multistitchimage_maxIP_unmixed_section.nd2_job1864MHbExperiment_Br5422_middleslice_20x_multistitchimage_maxIP_unmixed_second_pass_object_results.csv")
)

map(fn, file.exists)
halo <- map2(fn, names(fn), ~read_csv(.x) |> mutate(job = .y))

# all.equal(halo$job1862 |> select(-job), halo$job1864|> select(-job))

map(halo, dim)
# $job1862
# [1] 18084    55
# 
# $job1864
# [1] 11994    55
identical(colnames(halo$job1862), colnames(halo$job1864))
# [1] "Image Location"                  "Analysis Region"                 "Algorithm Name"                 
# [4] "Object Id"                       "XMin"                            "XMax"                           
# [7] "YMin"                            "YMax"                            "ALL THREE GENES CELLS"          
# [10] "NO POU4F1 ALL THREE GENES"       "NO POU4F1 CCK & EBF3 CELLS"      "NO POU4F1 CHAT & CCK LHb1 Cells"
# [13] "NO POU4F1 CHAT & EBF3 CELLS"     "NO POU4F1 CCK LHb1 CELLS"        "NO POU4F1 EBF3LHb4 CELLS"       
# [16] "NO FLOURESCENCE CELLS"           "CCk & EBF3 CELLS"                "CCK LHb1 CELLS"                 
# [19] "EBF3 LHb4 CELLS"                 "JUST POU4F1 CELLS"               "CHAT & CCK LHb1 CELLS"          
# [22] "CHAT LHb5 CELLS"                 "CHAT & EBF3 CELLS"               "NO POU4F1 CHAT LHb5 CELLS"      
# [25] "20x DAPI Positive"               "20x DAPI Positive Nucleus"       "20x DAPI Nucleus Intensity"     
# [28] "20x DAPI Positive Cytoplasm"     "20x DAPI Cytoplasm Intensity"    "25xSil Opal 520 Copies"         
# [31] "25xSil Opal 520 Area (µm²)"      "25xSil Opal 520 Classification"  "25xSil Opal 520 Cell Intensity" 
# [34] "25xSil Opal 520 Avg Intensity"   "20x Opal 570 Copies"             "20x Opal 570 Area (µm²)"        
# [37] "20x Opal 570 Classification"     "20x Opal 570 Cell Intensity"     "20x Opal 570 Avg Intensity"     
# [40] "20x Opal 620 Copies"             "20x Opal 620 Area (µm²)"         "20x Opal 620 Classification"    
# [43] "20x Opal 620 Cell Intensity"     "20x Opal 620 Avg Intensity"      "20x Opal 690 Copies"            
# [46] "20x Opal 690 Area (µm²)"         "20x Opal 690 Classification"     "20x Opal 690 Cell Intensity"    
# [49] "20x Opal 690 Avg Intensity"      "Cell Area (µm²)"                 "Cytoplasm Area (µm²)"           
# [52] "Nucleus Area (µm²)"              "Nucleus Perimeter (µm)"          "Nucleus Roundness"    

halo <- do.call("rbind", halo)

colnames(halo) <- gsub(" \\(µm.*\\)", "", colnames(halo))

#### create halo long ####
experiment <- tibble(probe = factor(c(690, 620, 570, 520)),
       marker = c("CCK", "EBF3","CHAT","POU4F1"),
       cluster = c("MHb.1", "Mhb.3", "Mhb.2", "Hb"))
# probe marker cluster
# <dbl> <chr>  <chr>  
# 1   690 CCK    MHb.1  
# 2   620 EBF3   Mhb.3  
# 3   570 CHAT   Mhb.2  
# 4   520 POU4F1 Hb 

halo_copies_long <- halo |>
  select(job, `Object Id`, XMin, XMax, YMin, YMax, ends_with("Copies")) |>
  pivot_longer(!c(job, `Object Id`, XMin, XMax, YMin, YMax), names_to = "probe", values_to = "copies") |>
  mutate(probe = factor(gsub("^2.*Opal (\\d+) Copies", "\\1", probe))) |>
  left_join(experiment) |>
  mutate(probe2 = paste0(probe, " ", marker, " (", cluster,")"))
  # mutate(probe = extract_numeric(probe))

# get quantile values for non-zero values
copy_cutoff <- halo_copies_long |> 
  filter(copies != 0) |> 
  group_by(probe, job) |> 
  summarize(q25 = quantile(copies,probs = 0.25),
            q50 = quantile(copies,probs = 0.5),
            q75 = quantile(copies,probs = 0.75),
            q95 = quantile(copies,probs = 0.95))

## different by job 
# probe job       q25   q50   q75   q95
# <fct> <chr>   <dbl> <dbl> <dbl> <dbl>
#   1 520   job1862     1     2     5  24  
# 2 520   job1864     1     2     6  24.1
# 3 570   job1862     1     2     5  27  
# 4 570   job1864     1     2     5  25  
# 5 620   job1862     1     1     2  13  
# 6 620   job1864     1     1     3  18  
# 7 690   job1862     1     2     7  30.3
# 8 690   job1864     1     2     7  32.1

## use to classify cells
halo_copies_long2 <- halo_copies_long |>
  left_join(copy_cutoff) |>
  mutate(copy_quant = case_when(copies > q75 ~ "q75",
                                copies > q50 ~ "q50",
                                copies > q25 ~ "q25",
                                copies > 0 ~ "q0",
                                copies == 0 ~ "None",
                                TRUE ~ NA
  )) |>
  mutate(copy_quant = as.ordered(copy_quant))

halo_copies_long2 |>
  count(probe2, job, copy_quant) |> 
  pivot_wider(names_from = "copy_quant", values_from = "n")

# probe2           job      None    q0   q25   q50   q75
# <chr>            <chr>   <int> <int> <int> <int> <int>
#   1 520 POU4F1 (Hb)  job1862 17632   195    63    92   102
# 2 520 POU4F1 (Hb)  job1864 11335   303    76   137   143
# 3 570 CHAT (Mhb.2) job1862 16527   685   229   290   353
# 4 570 CHAT (Mhb.2) job1864 10264   741   250   338   401
# 5 620 EBF3 (Mhb.3) job1862 16343  1026    NA   298   417
# 6 620 EBF3 (Mhb.3) job1864 10372   842    NA   441   339
# 7 690 CCK (MHb.1)  job1862 17710   161    51    78    84
# 8 690 CCK (MHb.1)  job1864 11475   235    62   102   120

halo |> filter(`Object Id` == 0)

#### plot with cell_quant 
cell_count_quant <- halo_copies_long2 |> 
  ggplot() +
  geom_rect(aes(
    xmin = XMin, xmax = XMax,
    ymin = YMin, ymax = YMax,
    fill = copy_quant
  )) +
  scale_fill_manual(values = copy_quant_colors, "Copy Quantile") +
  coord_equal()+
  theme_bw() +
  facet_grid(job~probe2)

ggsave(cell_count_quant, filename = here(plot_dir, paste0("MHb_cell_count_quant_facet.png")), height = 7, width = 9)
ggsave(cell_count_quant, filename = here(plot_dir, paste0("MHb_cell_count_quant_facet.pdf")), height = 7, width = 9)


#### distribution

copy_quant_colors <- c(None = "grey75", q0 = "#FECC5C", q25 = "#FD8D3C", q50 = "#F03B20", q75 = "#BD0026")

copies_density <- halo_copies_long2 |>
  # filter(job == "job1862") |>
  # filter(copies != 0) |>
  ggplot(aes(x = copies, fill = copy_quant)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  scale_fill_manual(values = copy_quant_colors, "Copy Quantile") +
  coord_cartesian(ylim=c(0, 1250), xlim = c(-1, 20)) +
  facet_grid(job~probe2)

ggsave(copies_density, filename = here(plot_dir, "MHbcopies_denisty.png"), height = 5)

#### hex plots ####

hex_copies_median <- ggplot(halo_copies_long2) +
  stat_summary_hex(aes(x = XMax, y = YMax, z = copies),
                   fun = median, bins = 100
  ) +
  # scale_fill_continuous(type = "viridis") + ## top value of 100 for visualization 
  scale_fill_continuous(type = "viridis", limits = c(1,20), "Median Copies\n(max 20)") + ## top value of 100 for visualization
  coord_equal() +
  theme_bw() +
  facet_grid(job~probe2)

ggsave(hex_copies_median, filename = here(plot_dir, paste0("hex_copies_median_facet.png")), height = 9, width = 12)



