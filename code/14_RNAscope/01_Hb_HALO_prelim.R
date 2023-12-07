
library("tidyverse")
library("readxl")
library("GGally")

library("here")
library("sessioninfo")

## set up 
plot_dir <- here("plots", "14_RNAscope", "01_Hb_HALO_prelim")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load data
fn <- here("processed-data", "14_RNAscope", "HALO_data", "Br5422_middleslice_20x_zstack_tiledlargeimage_multipointstitch_unmixed_maxIP_redo_tutorial.nd2_job1805middleslice_secondpass_CpI_median_non0_cell_intensity_object_results.xlsx")
file.exists(fn)
halo_prelim <- read_excel(fn)

dim(halo_prelim)
# [1] 18621    54
colnames(halo_prelim)
# [1] "Image Location"                        "Analysis Region"                      
# [3] "Algorithm Name"                        "Object Id"                            
# [5] "XMin"                                  "XMax"                                 
# [7] "YMin"                                  "YMax"                                 
# [9] "ALL THREE GENES CELLS"                 "NO POU4F1 ALL THREE GENES"            
# [11] "NO POU4F1 TLE2 & ONECUT2 CELLS"        "NO POU4F1 ONECUT2 & SEMA3D LHb1 Cells"
# [13] "NO POU4F1 SEMA3D&TLE2 CELLS"           "NO POU4F1 ONECUT2 LHb1 CELLS"         
# [15] "NO POU4F1 TLE2 LHb4 CELLS"             "NO FLOURESCENCE CELLS"                
# [17] "TLE2 & ONECUT2 CELLS"                  "ONECUT2 LHb1 CELLS"                   
# [19] "TLE2 LHb4 CELLS"                       "JUST POU4F1 CELLS"                    
# [21] "ONECUT2 SEMA3D LHb1 CELLS"             "SEMA3D LHb5 CELLS"                    
# [23] "SEMA3D & TLE2 CELLS"                   "NO POU4F1 SEMA3D LHb5 CELLS"          
# [25] "DAPI Positive"                         "DAPI Positive Nucleus"                
# [27] "DAPI Nucleus Intensity"                "DAPI Positive Cytoplasm"              
# [29] "DAPI Cytoplasm Intensity"              "520 POU4F1 Copies"                    
# [31] "520 POU4F1 Area (µm²)"                 "520 POU4F1 Classification"            
# [33] "520 POU4F1 Cell Intensity"             "520 POU4F1 Avg Intensity"             
# [35] "570 SEMA3D Copies"                     "570 SEMA3D Area (µm²)"                
# [37] "570 SEMA3D Classification"             "570 SEMA3D Cell Intensity"            
# [39] "570 SEMA3D Avg Intensity"              "620 TLE2 Copies"                      
# [41] "620 TLE2 Area (µm²)"                   "620 TLE2 Classification"              
# [43] "620 TLE2 Cell Intensity"               "620 TLE2 Avg Intensity"               
# [45] "690 ONECUT2 Copies"                    "690 ONECUT2 Area (µm²)"               
# [47] "690 ONECUT2 Classification"            "690 ONECUT2 Cell Intensity"           
# [49] "690 ONECUT2 Avg Intensity"             "Cell Area (µm²)"                      
# [51] "Cytoplasm Area (µm²)"                  "Nucleus Area (µm²)"                   
# [53] "Nucleus Perimeter (µm)"                "Nucleus Roundness"

## simplify colnames
colnames(halo_prelim) <- gsub(" \\(µm.*\\)", "", colnames(halo_prelim))


#### gpairs for copies ####
copies_tab <- halo_prelim |> 
  select(ends_with("Copies")) |>
  mutate(sum = rowSums(across(where(is.numeric))))

copies_tab |> count(sum == 0)
# `sum == 0`     n
# <lgl>      <int>
# 1 FALSE      14210
# 2 TRUE        4411

gg_copies <- copies_tab |> select(-sum) |> ggpairs(aes(alpha = 0.5))
ggsave(gg_copies, filename = here(plot_dir, "gg_copies.png"), height = 12, width = 12)

## intensity
halo_prelim |> 
  select(ends_with("Avg Intensity"))

gg_ave_intensity <- halo_prelim |> 
  select(ends_with("Avg Intensity")) |>
  ggpairs()
ggsave(gg_ave_intensity, filename = here(plot_dir, "gg_ave_intensity.png"), height = 9, width = 9)

gg_cell_intensity <- halo_prelim |> 
  select(ends_with("Cell Intensity")) |>
  ggpairs()
ggsave(gg_cell_intensity, filename = here(plot_dir, "gg_cell_intensity.png"), height = 9, width = 9)


probes <- c("520 POU4F1", "570 SEMA3D","620 TLE2","690 ONECUT2")

walk(probes, ~{
  halo_prelim |> 
    select(starts_with(.x)) |>
    ggpairs() |>
    ggsave(filename = here(plot_dir, paste0("gg_",gsub("^\\d+ ", "", .x),"_intensity.png")), height = 9, width = 9)
  })

#### Hex plots ####

## just one sample Br5422_middle
halo_prelim |> count(`Image Location`)

hex_nuc_count <- ggplot(halo_prelim) +
  geom_hex(aes(x = XMax, y = YMax), bins = 100) +
  scale_fill_continuous(type = "viridis") + 
  coord_equal() +
  theme_bw()

ggsave(hex_nuc_count, filename = here(plot_dir, "hex_nuc_nount.png"))



summary(halo_prelim$`Nucleus Area`)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 10.00   31.24   51.56   50.11   65.21  264.86 

halo_prelim |> arrange(`Nucleus Area`)

mean_nuc_area_hex <- ggplot(halo_prelim) +
  stat_summary_hex(aes(x = XMax, y = YMax, z = `Nucleus Area`),
                   fun = mean, bins = 100
  ) +
  scale_fill_continuous(type = "viridis") + ## top value of 100 for visualization 
  # scale_fill_continuous(type = "viridis", limits = c(5,100), "Mean Nuc Area\n(max 100)") + ## top value of 100 for visualization 
  coord_equal() +
  theme_bw()

ggsave(mean_nuc_area_hex, filename = here(plot_dir, "mean_nuc_area_hex.png"))

### gene expression

mean_nuc_area_hex <- ggplot(halo_prelim) +
  stat_summary_hex(aes(x = XMax, y = YMax, z = `Nucleus Area`),
                   fun = mean, bins = 100
  ) +
  scale_fill_continuous(type = "viridis") + ## top value of 100 for visualization 
  # scale_fill_continuous(type = "viridis", limits = c(5,100), "Mean Nuc Area\n(max 100)") + ## top value of 100 for visualization 
  coord_equal() +
  theme_bw()

ggsave(mean_nuc_area_hex, filename = here(plot_dir, "mean_nuc_area_hex.png"))


probes <- c("520 POU4F1", "570 SEMA3D","620 TLE2","690 ONECUT2")
names(probes) <- gsub("^\\d+ ", "", probes)
walk(probes, ~{
  
  p <- gsub("^\\d+ ", "", .x)
  
  hex <- ggplot(halo_prelim) +
    stat_summary_hex(aes(x = XMax, y = YMax, z = .data[[paste(.x, "Copies")]]),
                     fun = median, bins = 100
    ) +
    # scale_fill_continuous(type = "viridis") + ## top value of 100 for visualization 
    scale_fill_continuous(type = "viridis", limits = c(1,100), "Median Copies\n(max 100)") + ## top value of 100 for visualization
    coord_equal() +
    theme_bw() +
    labs(title = p)
  
    ggsave(hex, filename = here(plot_dir, paste0("hex_", p,"_copies_median.png")), height = 9, width = 9)
  
  })

## facet 
halo_copies_long <- halo_prelim |>
  select(`Object Id`, XMin, XMax, YMin, YMax, ends_with("Copies")) |>
  pivot_longer(!c(`Object Id`, XMin, XMax,YMin, YMax), names_to = "probe", values_to = "copies")

hex_copies_median <- ggplot(halo_copies_long) +
  stat_summary_hex(aes(x = XMax, y = YMax, z = copies),
                   fun = median, bins = 100
  ) +
  # scale_fill_continuous(type = "viridis") + ## top value of 100 for visualization 
  scale_fill_continuous(type = "viridis", limits = c(1,100), "Median Copies\n(max 100)") + ## top value of 100 for visualization
  coord_equal() +
  theme_bw() +
  facet_wrap(~probe)

ggsave(hex_copies_median, filename = here(plot_dir, paste0("hex_copies_median_facet.png")), height = 9, width = 9)


hex_copies_mean <- ggplot(halo_copies_long) +
  stat_summary_hex(aes(x = XMax, y = YMax, z = copies),
                   fun = mean, bins = 100
  ) +
  # scale_fill_continuous(type = "viridis") + ## top value of 100 for visualization 
  scale_fill_continuous(type = "viridis", limits = c(1,100), "Median Copies\n(max 100)") + ## top value of 100 for visualization
  coord_equal() +
  theme_bw() +
  facet_wrap(~probe)

ggsave(hex_copies_mean, filename = here(plot_dir, paste0("hex_copies_mean_facet.png")), height = 9, width = 9)


hex_copies_any <- ggplot(halo_copies_long) +
  stat_summary_hex(aes(x = XMax, y = YMax, z = copies > 0),
                   fun = any, bins = 100
  ) +
  # scale_fill_continuous(type = "viridis") + ## top value of 100 for visualization 
  # scale_fill_continuous(type = "viridis", limits = c(1,100), "Median Copies\n(max 100)") + ## top value of 100 for visualization
  coord_equal() +
  theme_bw() +
  facet_wrap(~probe)


ggsave(hex_copies_any, filename = here(plot_dir, paste0("hex_copies_any_facet.png")), height = 9, width = 9)

#### distribution

copies_density <- halo_copies_long |>
  filter(copies != 0) |>
  ggplot(aes(x = copies)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  xlim(0,50) +
  facet_wrap(~probe, ncol =1)

ggsave(copies_density, filename = here(plot_dir, "copies_denisty.png"))


halo_prelim |> count(POUF1_present = `520 POU4F1 Copies` > 0)
# POUF1_present     n
# <lgl>         <int>
# 1 FALSE         14555
# 2 TRUE           4066


cell_copies_facet <- halo_copies_long |>
  filter(copies > 0) |>
  ggplot() +
  geom_rect(aes(
    xmin = XMin, xmax = XMax,
    ymin = YMin, ymax = YMax,
    fill = copies
  )) +
  scale_fill_continuous(type = "viridis", limits = c(1,50), "Copies\n(max 50)")+
  coord_equal()+
  theme_bw() +
  facet_wrap(~probe)

ggsave(cell_copies_facet, filename = here(plot_dir, paste0("cell_copies_facet.png")), height = 9, width = 9)
ggsave(cell_copies_facet, filename = here(plot_dir, paste0("cell_copies_facet.pdf")), height = 9, width = 9)


## with POU4F1
halo_POU4F1_long <- halo_prelim |>
  select(`Object Id`, XMin, XMax, YMin, YMax, ends_with("Copies")) |>
  rename(POU4F1 = `520 POU4F1 Copies`) |>
  pivot_longer(!c(`Object Id`, XMin, XMax, YMin, YMax, POU4F1), names_to = "probe", values_to = "copies")

halo_POU4F1_long |>
  group_by(`Object Id`) |>
  arrange(-copies) |>
  slice(1) |>
  ungroup() |>
  count(POU4F1, probe)

halo_POU4F1_long |>
  filter(copies > 0) |>
  group_by(probe) |>
  count(POU4F1 > 0)

# probe              `POU4F1 > 0`     n
# <chr>              <lgl>        <int>
# 1 570 SEMA3D Copies FALSE          7958
# 2 570 SEMA3D Copies  TRUE          2814
# 3 620 TLE2 Copies    FALSE         6358
# 4 620 TLE2 Copies    TRUE          2275
# 5 690 ONECUT2 Copies FALSE         2489
# 6 690 ONECUT2 Copies TRUE           957


cell_max_filter <- halo_POU4F1_long |>
  group_by(`Object Id`) |>
  arrange(-copies) |>
  slice(1) |>
  filter(copies >= 5) |>
  ggplot() +
  geom_rect(aes(
    xmin = XMin, xmax = XMax,
    ymin = YMin, ymax = YMax,
    fill = probe
  )) +
  coord_equal()+
  theme_bw() 

ggsave(cell_max_filter, filename = here(plot_dir, paste0("cell_max_filter.png")), height = 9, width = 9)
ggsave(cell_max_filter, filename = here(plot_dir, paste0("cell_max_filter.pdf")), height = 9, width = 9)



#### Cell Type Annotations ####

# Which of those copy count bins are most useful for spatially visualizing our "cell types" of interest

halo_prelim |> count(`520 POU4F1 Classification`)
halo_prelim |> count(`570 SEMA3D Classification`)
halo_prelim |> count(`620 TLE2 Classification`)
halo_prelim |> count(`690 ONECUT2 Classification`)



cell_slide <- halo_prelim |>
  ggplot() +
  geom_rect(aes(
    xmin = XMin, xmax = XMax,
    ymin = YMin, ymax = YMax,
    fill = as.factor(`520 POU4F1 Classification`)
  )) +
  scale_fill_manual(values = c(`0` = "grey75", `1` = "#FECC5C", `2` = "#FD8D3C", `3` = "#F03B20", `4` = "#BD0026"), "class") +
  # scale_fill_manual(values = c(`0` = "grey75", `1` = "#D7191C", `2` = "#FDAE61", `3` = "#ABDDA4", `4` = "#2B83BA")) +
  coord_equal()+
  theme_bw() 

ggsave(cell_slide, filename = here(plot_dir, "cell_slide.png"))

## facet
halo_class_long <- halo_prelim |>
  select(`Object Id`, XMin, XMax,YMin, YMax, ends_with("Classification")) |>
  pivot_longer(!c(`Object Id`, XMin, XMax,YMin, YMax), names_to = "probe", values_to = "class") |>
  mutate(class = as.factor(class))

cell_class_facet <- ggplot(halo_class_long) +
  geom_rect(aes(
    xmin = XMin, xmax = XMax,
    ymin = YMin, ymax = YMax,
    fill = class
  )) +
  scale_fill_manual(values = c(`0` = "grey75", `1` = "#FECC5C", `2` = "#FD8D3C", `3` = "#F03B20", `4` = "#BD0026"), "class") +
  # scale_fill_manual(values = c(`0` = "grey75", `1` = "#D7191C", `2` = "#FDAE61", `3` = "#ABDDA4", `4` = "#2B83BA")) +
  coord_equal()+
  theme_bw() +
  facet_wrap(~probe)

ggsave(cell_class_facet, filename = here(plot_dir, paste0("cell_class_facet.png")), height = 9, width = 9)
ggsave(cell_class_facet, filename = here(plot_dir, paste0("cell_class_facet.pdf")), height = 9, width = 9)

## filter
cell_class_facet_filter <- halo_class_long |>
  filter(as.integer(class) > 2) |>
  ggplot() +
  geom_rect(aes(
    xmin = XMin, xmax = XMax,
    ymin = YMin, ymax = YMax,
    fill = class
  )) +
  scale_fill_manual(values = c(`0` = "grey75", `1` = "#FECC5C", `2` = "#FD8D3C", `3` = "#F03B20", `4` = "#BD0026"), "class") +
  # scale_fill_manual(values = c(`0` = "grey75", `1` = "#D7191C", `2` = "#FDAE61", `3` = "#ABDDA4", `4` = "#2B83BA")) +
  coord_equal()+
  theme_bw() +
  facet_wrap(~probe)

ggsave(cell_class_facet_filter, filename = here(plot_dir, paste0("cell_class_facet_filter.png")), height = 9, width = 9)
ggsave(cell_class_facet_filter, filename = here(plot_dir, paste0("cell_class_facet_filter.pdf")), height = 9, width = 9)



