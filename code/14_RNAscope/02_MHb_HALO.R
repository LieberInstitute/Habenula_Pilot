
library("tidyverse")
library("GGally")
library("here")
library("sessioninfo")
library("patchwork")

## set up
plot_dir <- here("plots", "14_RNAscope", "02_MHb_HALO")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

list.files("~/OneDrive - Johns Hopkins/habenulaPilot-paper/HALO_data/Medial_exp/", recursive = TRUE)
# [1] "Br5422/MHbExperiment_Br5422_middleslice_20x_multistitchimage_maxIP_unmixed.nd2_Sample1871MHbExperiment_Br5422_middle_secondpass_object_results.csv"
# [2] "Josh Ege Redo/MHb_ExperimentJ_Br8112_20x_EgeRedo_Sample1883_object_RESULTS.csv"

## load data
## local
datadir <- "~/OneDrive - Johns Hopkins/habenulaPilot-paper/HALO_data/Medial_exp/"

## jhpce
# datadir <- here("processed-data", "14_RNAscope", "HALO_data")

fn <- list(Br5422 = "Br5422/MHbExperiment_Br5422_middleslice_20x_multistitchimage_maxIP_unmixed.nd2_Sample1871MHbExperiment_Br5422_middle_secondpass_object_results.csv",
           Br8112 = "Josh Ege Redo/MHb_ExperimentJ_Br8112_20x_EgeRedo_Sample1883_object_RESULTS.csv"
)

fn <- map(fn, ~paste0(datadir, .x))

map(fn, file.exists)
halo <- map2(fn, names(fn), ~read_csv(.x) |> mutate(Sample = .y))

map(halo, dim)
# $Br5422
# [1] 18418    55
#
# $Br8112
# [1] 3098   55

## subset to common colnames
common_colnames <- intersect(colnames(halo$Br5422), colnames(halo$Br8112))
halo <- map(halo, ~.x[,common_colnames])

halo <- do.call("rbind", halo)

colnames(halo) <- gsub(" \\(µm.*\\)", "", colnames(halo))

#### create halo long ####
experiment <- rbind(
    tibble(probe = factor(c(690, 620, 570, 520)),
       marker = c("CCK", "EBF3","CHAT","POU4F1"),
       cluster = c("MHb.1", "Mhb.3", "Mhb.2", "Hb"),
       Sample = "Br5422"),
    tibble(probe = factor(c(690, 620, 570, 520)),
           marker = c("CCK", "BHLHE22","CHAT","CHRNB4"), ## * alt probes used in Josh's run
           cluster = c("MHb.1", "Mhb.3", "Mhb.2", "Hb"),
           Sample = "Br8112")) |>
    mutate(probe2 = paste0(probe, " ", marker, " (", cluster,")"))

# probe marker  cluster Sample probe2
# <fct> <chr>   <chr>   <chr>  <chr>
# 1 690   CCK     MHb.1   Br5422 690 CCK (MHb.1)
# 2 620   EBF3    Mhb.3   Br5422 620 EBF3 (Mhb.3)
# 3 570   CHAT    Mhb.2   Br5422 570 CHAT (Mhb.2)
# 4 520   POU4F1  Hb      Br5422 520 POU4F1 (Hb)
# 5 690   CCK     MHb.1   Br8112 690 CCK (MHb.1)
# 6 620   BHLHE22 Mhb.3   Br8112 620 BHLHE22 (Mhb.3) *
# 7 570   CHAT    Mhb.2   Br8112 570 CHAT (Mhb.2)
# 8 520   CHRNB4  Hb      Br8112 520 CHRNB4 (Hb) *

write.csv(experiment, file = here("processed-data", "14_RNAscope", "HALO_data", "Medial_exp", "Probes_MHb.csv"))

halo_copies_long <- halo |>
    select(Sample, `Object Id`, XMin, XMax, YMin, YMax, ends_with("Copies")) |>
    pivot_longer(!c(Sample, `Object Id`, XMin, XMax, YMin, YMax), names_to = "probe", values_to = "copies") |>
    mutate(probe = factor(gsub("^2.*Opal (\\d+) Copies", "\\1", probe))) |>
    left_join(experiment) |>
    filter((Sample != "Br8112" | probe !=  620)) ## this channel failed for this Sample

halo_copies_long |> count(Sample, probe2)
#   Sample probe2               n
# <chr>  <chr>            <int>
# 1 Br5422 520 POU4F1 (Hb)  18418
# 2 Br5422 570 CHAT (Mhb.2) 18418
# 3 Br5422 620 EBF3 (Mhb.3) 18418
# 4 Br5422 690 CCK (MHb.1)  18418
# 5 Br8112 520 CHRNB4 (Hb)   3098
# 6 Br8112 570 CHAT (Mhb.2)  3098
# 7 Br8112 690 CCK (MHb.1)   3098

#### ggpairs ####
halo_copies_wide <- halo_copies_long |>
    select(Sample, `Object Id`, copies, probe2) |>
    pivot_wider(id_cols = c(Sample, `Object Id`), names_from = probe2, values_from = copies)

# copies_tab |> count(sum == 0)
# copies_tab |> count(Sample, sum == 0)

halo_copies_wide |>
    filter(Sample == "Br5422") |>
    ggpairs(3:6) |> ## not all markers exist in both samples
    ggsave(filename = here(plot_dir, "MHb_gg_copies_Br5422.png"), height = 12, width = 12)

halo_copies_wide |>
    filter(Sample == "Br8112") |>
    ggpairs(c(4,6,7)) |> ## not all markers exist in both samples
    ggsave(filename = here(plot_dir, "MHb_gg_copies_Br8112.png"), height = 12, width = 12)
## expect  620 (BHLHE22) to be all 0s - filter earlier

#### quantile grouping ####
# get quantile values for non-zero values
copy_cutoff <- halo_copies_long |>
  filter(copies != 0) |>
  group_by(probe, Sample) |>
  summarize(q25 = quantile(copies,probs = 0.25),
            q50 = quantile(copies,probs = 0.5),
            q75 = quantile(copies,probs = 0.75),
            q95 = quantile(copies,probs = 0.95))


## use to classify cells
halo_copies_long_quant  <- halo_copies_long |>
  left_join(copy_cutoff) |>
  mutate(copy_quant = case_when(copies > q75 ~ "q75",
                                copies > q50 ~ "q50",
                                copies > q25 ~ "q25",
                                copies > 0 ~ "q0",
                                copies == 0 ~ "None",
                                TRUE ~ NA
  )) |>
  mutate(copy_quant = as.ordered(copy_quant))

halo_copies_long_quant  |>
  count(probe2, Sample, copy_quant) |>
  pivot_wider(names_from = "copy_quant", values_from = "n")

#### plot with cell_quant
copy_quant_colors <- c(None = "grey75", q0 = "#FECC5C", q25 = "#FD8D3C", q50 = "#F03B20", q75 = "#BD0026")

cell_count_quant <- halo_copies_long_quant  |>
  ggplot() +
  geom_rect(aes(
    xmin = XMin, xmax = XMax,
    ymin = YMin, ymax = YMax,
    fill = copy_quant
  )) +
  scale_fill_manual(values = copy_quant_colors, "Copy Quantile") +
  coord_equal()+
  theme_bw() +
  facet_grid(Sample~probe2)

ggsave(cell_count_quant, filename = here(plot_dir, paste0("MHb_cell_count_quant_facet.png")), height = 7, width = 9)
ggsave(cell_count_quant, filename = here(plot_dir, paste0("MHb_cell_count_quant_facet.pdf")), height = 7, width = 9)


#### distribution
copies_density <- halo_copies_long_quant  |>
  # filter(Sample == "Sample1862") |>
  # filter(copies != 0) |>
  ggplot(aes(x = copies, fill = copy_quant)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  scale_fill_manual(values = copy_quant_colors, "Copy Quantile") +
  coord_cartesian(ylim=c(0, 1250), xlim = c(-1, 20)) +
  facet_grid(Sample~probe2)

ggsave(copies_density, filename = here(plot_dir, "MHbcopies_denisty.png"), height = 5, width = 9)

#### hex plots ####
hex_copies_median <- ggplot(halo_copies_long) +
    stat_summary_hex(aes(x = XMax, y = YMax, z = copies),
                     fun = median, bins = 100
    ) +
    scale_fill_gradientn(
        name = "Median Copies",
        colors = rev(viridisLite::rocket(21)),
        na.value = "#CCCCCC50",
        limits = c(1,200)
    )+
    coord_equal() +
    theme_bw() +
    facet_grid(Sample~probe2)

ggsave(hex_copies_median, filename = here(plot_dir, paste0("MHb_hex_copies_median_facet.png")), height = 6, width = 9)

hex_copies_max <- ggplot(halo_copies_long) +
    stat_summary_hex(aes(x = XMax, y = YMax, z = copies),
                     fun = max, bins = 100
    ) +
    # scale_fill_continuous(type = "viridis") + ## top value of 100 for visualization
    scale_fill_gradientn(
        name = "Max Copies\n(1:200)",
        colors = rev(viridisLite::rocket(21)),
        na.value = "#CCCCCC50",
        limits = c(1,200) ## cap color scale don't filter data
    )+ coord_equal() +
    theme_bw() +
    facet_grid(Sample~probe2)

ggsave(hex_copies_max, filename = here(plot_dir, paste0("MHb_hex_copies_max_facet.png")), height = 6, width = 9)
ggsave(hex_copies_max, filename = here(plot_dir, paste0("MHb_hex_copies_max_facet.pdf")), height = 6, width = 9)


#### Bin by 100 ####
halo_copies_rank <- halo_copies_long |>
    group_by(Sample, probe) |>
    arrange(-copies) |>
    mutate(copies_rank = rank(-copies)) |>
    mutate(rank_cut = cut(copies_rank, c(0, 100, 200, 300, 400)) ,
           top400 = !is.na(rank_cut))

halo_copies_rank |> count(rank_cut)
halo_copies_rank |> count(top400)

halo_copies_rank |>
    filter(rank_cut == "(0,100]") |>
    arrange(copies) |>
    slice(1)

# Sample `Object Id`  XMin  XMax  YMin  YMax probe copies marker cluster probe2           copies_rank rank_cut top400
# <chr>        <dbl> <dbl> <dbl> <dbl> <dbl> <fct>  <dbl> <chr>  <chr>   <chr>                  <dbl> <fct>    <lgl>
# 1 Br5422        1247  6775  6802  2872  2902 520        5 POU4F1 Hb      520 POU4F1 (Hb)         87   (0,100]  TRUE
# 2 Br5422        4712  3963  4004  8686  8716 570       20 CHAT   Mhb.2   570 CHAT (Mhb.2)        97   (0,100]  TRUE
# 3 Br5422         605  6714  6743  1407  1446 620        3 EBF3   Mhb.3   620 EBF3 (Mhb.3)        89   (0,100]  TRUE #lowest
# 4 Br5422        1239  5817  5857  3183  3217 690       20 CCK    MHb.1   690 CCK (MHb.1)         97   (0,100]  TRUE
# 5 Br8112         599  2145  2184  2658  2691 520       12 CHRNB4 Hb      520 CHRNB4 (Hb)         92   (0,100]  TRUE
# 6 Br8112        1179  1656  1678  4524  4547 570       32 CHAT   Mhb.2   570 CHAT (Mhb.2)        99   (0,100]  TRUE
# 7 Br8112        1930   900   923  6667  6689 690        8 CCK    MHb.1   690 CCK (MHb.1)         94.5 (0,100]  TRUE

halo_copies_rank |>
    filter(rank_cut == "(300,400]") |>
    arrange(copies) |>
    slice(1)

#### confusion matrix ###
halo_copies_cat <- halo_copies_rank |>
    filter(!is.na(rank_cut)) |>
    filter(rank_cut == "(0,100]") |>
    ungroup() |>
    mutate(cat = paste0(probe, "_", marker, "_", rank_cut)) |>
    select(Sample, `Object Id`, cat)

halo_copies_cat2 <- halo_copies_cat |>
    left_join(halo_copies_cat, by = join_by(Sample, `Object Id`), relationship = "many-to-many")

confusion_top100 <- halo_copies_cat2 |>
    count(Sample, cat.x, cat.y) |>
    mutate(Sample2 = Sample) |> ## hack to make facet work
    group_by(Sample2) |>
    group_map(
        ~ggplot(.x, aes(cat.x, cat.y, fill = n)) +
            geom_tile() +
            geom_text(aes(label = n), color = "white") +
            facet_wrap(~Sample) +
            # scale_fill_gradient(name = "count", trans = "log") +
            theme_bw() +
            theme(axis.text.x = element_text(angle = 90),
                  axis.title.x=element_blank(),
                  axis.title.y=element_blank())
    )

ggsave(confusion_top100[[1]] + confusion_top100[[2]], filename = here(plot_dir, "MHb_confusion_top100.png"), height = 5, width = 11)


