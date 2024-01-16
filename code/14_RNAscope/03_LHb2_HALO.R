
library("tidyverse")
library("readxl")
library("GGally")
library("here")
library("sessioninfo")

## set up
plot_dir <- here("plots", "14_RNAscope", "03_LHb2_HALO")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load data
## local
list.files("~/OneDrive - Johns Hopkins/habenulaPilot-paper/HALO_data/Lateral_exp2/", recursive = TRUE)
# [1] "Br6462/LHbExperiment_Br6462_secondpass.xls"
# [2] "Br8112/LHbExperiment2_Br8112_thirdtoleft_20x_Sample1881__object_RESULTS.csv"

## load data
## local
datadir <- "~/OneDrive - Johns Hopkins/habenulaPilot-paper/HALO_data/Lateral_exp2/"
## JHPCE
# datadir <- here("processed-data", "14_RNAscope", "HALO_data", "Lateral_exp1")

halo <- list(
    Br6462 = read_excel(paste0(datadir, "Br6462/LHbExperiment_Br6462_secondpass.xls")) |> mutate(Sample = "Br6462"),
    Br8112 = read_csv(paste0(datadir, "Br8112/LHbExperiment2_Br8112_thirdtoleft_20x_job1881__object_RESULTS.csv")) |> mutate(Sample = "Br8112")
    )

map_int(halo, ncol)
map_int(halo, nrow)
# Br6462 Br8112
# 16383  15395

# ## flip y for Br8112
halo$Br8112$YMax <- max(halo$Br8112$YMax) - halo$Br8112$YMax
halo$Br8112$YMin <- max(halo$Br8112$YMax) - halo$Br8112$YMin

# halo$Br8112 |>
#     ggplot() +
#     geom_rect(aes(
#         xmin = XMin, xmax = XMax,
#         ymin = YMin, ymax = YMax
#     ))

## subset to common colnames
common_colnames <- intersect(colnames(halo$Br6462), colnames(halo$Br8112))
halo <- map(halo, ~.x[,common_colnames])
## rbind tables
halo <- do.call("rbind", halo)
colnames(halo) <- gsub(" \\(µm.*\\)", "", colnames(halo))

#### create halo long ####
experiment <- tibble(probe = factor(c(690, 620, 570, 520)),
       marker = c("ESRP1", "MCOLN3","CRH","POU4F1"),
       cluster = c("LHb.6", "LHb.3", "LHb.2", "Hb")) |>
    mutate(probe2 = paste0(probe, " ", marker, " (", cluster,")"))
# probe marker cluster probe2
# <fct> <chr>  <chr>   <chr>
# 1 690   ESRP1  LHb.6   690 ESRP1 (LHb.6)
# 2 620   MCOLN3 LHb.3   620 MCOLN3 (LHb.3)
# 3 570   CRH    LHb.2   570 CRH (LHb.2)
# 4 520   POU4F1 Hb      520 POU4F1 (Hb)

# write.csv(experiment, file = here("processed-data", "14_RNAscope", "HALO_data", "Lateral_exp2", "Probes_LHb2.csv"))

#### create halo long ####
halo_copies_long <- halo |>
  select(Sample, `Object Id`, XMin, XMax, YMin, YMax, ends_with("Copies")) |>
  pivot_longer(!c(Sample, `Object Id`, XMin, XMax, YMin, YMax), names_to = "probe", values_to = "copies") |>
  mutate(probe = factor(gsub("^2.*Opal (\\d+) Copies", "\\1", probe))) |>
  left_join(experiment)

#### ggpairs ####
halo_copies_wide <- halo_copies_long |>
    select(Sample, `Object Id`, copies, probe2) |>
    pivot_wider(id_cols = c(Sample, `Object Id`), names_from = probe2, values_from = copies)

# copies_tab |> count(sum == 0)
# copies_tab |> count(Sample, sum == 0)

gg_copies <- halo_copies_wide |> ggpairs(3:6, mapping=ggplot2::aes(colour = Sample))
ggsave(gg_copies, filename = here(plot_dir, "LHb2_gg_copies.png"), height = 12, width = 12)


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
halo_copies_long_quant <- halo_copies_long |>
  left_join(copy_cutoff) |>
  mutate(copy_quant = case_when(copies > q75 ~ "q75",
                                copies > q50 ~ "q50",
                                copies > q25 ~ "q25",
                                copies > 0 ~ "q0",
                                copies == 0 ~ "None",
                                TRUE ~ NA
  )) |>
  mutate(copy_quant = as.ordered(copy_quant))

halo_copies_long_quant |>
  count(probe2, Sample, copy_quant) |>
  pivot_wider(names_from = "copy_quant", values_from = "n")


#### plot with cell_quant
copy_quant_colors <- c(None = "grey75", q0 = "#FECC5C", q25 = "#FD8D3C", q50 = "#F03B20", q75 = "#BD0026")

cell_count_quant <- halo_copies_long_quant |>
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

ggsave(cell_count_quant, filename = here(plot_dir, paste0("LHb2_cell_count_quant_facet.png")), height = 4, width = 9)
ggsave(cell_count_quant, filename = here(plot_dir, paste0("LHb2_cell_count_quant_facet.pdf")), height = 4, width = 9)


#### distribution
copies_density <- halo_copies_long_quant |>
  # filter(Sample == "Sample1862") |>
  # filter(copies != 0) |>
  ggplot(aes(x = copies, fill = copy_quant)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  scale_fill_manual(values = copy_quant_colors, "Copy Quantile") +
  coord_cartesian(ylim=c(0, 2000), xlim = c(-1, 50)) +
  facet_grid(Sample~probe2)

ggsave(copies_density, filename = here(plot_dir, "LHb2_copies_denisty.png"), height = 5)

#### hex plots ####
hex_copies_median <- ggplot(halo_copies_long_quant) +
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

ggsave(hex_copies_median, filename = here(plot_dir, paste0("LHb2_hex_copies_median_facet.png")), height = 6, width = 9)

## max hex for supp
hex_copies_max <- halo_copies_long |>
    mutate(copies = ifelse(copies > 200, 200, copies)) |> # cap data at 200 counts for visualization
    ggplot() +
    stat_summary_hex(aes(x = XMax, y = YMax, z = copies),
                     fun = max, bins = 100
    ) +
    # scale_fill_continuous(type = "viridis") + ## top value of 100 for visualization
    scale_fill_gradientn(
        name = "Max Copies\n(capped at 200)",
        colors = rev(viridisLite::rocket(21)),
        na.value = "#CCCCCC50",
        limits = c(1,NA)
    )+ coord_equal() +
    theme_bw() +
    facet_grid(Sample~probe2) +
    theme(legend.position = "bottom")

ggsave(hex_copies_max, filename = here(plot_dir, paste0("LHb2_hex_copies_max_facet.png")), height = 5, width = 7)
ggsave(hex_copies_max, filename = here(plot_dir, paste0("LHb2_hex_copies_max_facet.pdf")), height = 5, width = 7)


## max quant cell_max_quant <- halo_copies_long_quant |>
cell_max_quant <- halo_copies_long_quant |>
  group_by(`Object Id`) |>
  filter(probe != "520",
         copies != 0) |>
  arrange(-copies) |>
  slice(1)

cell_max_quant |> group_by(probe2, Sample) |> count()
cell_max_quant |> group_by(probe2, Sample) |> count(copy_quant)

cell_max_quant |>
  filter(copy_quant >= "q50") |>
  group_by(probe2, Sample) |>
  count()

cell_max_quant_plot <- cell_max_quant |>
  filter(copy_quant >= "q50") |>
  ggplot() +
  geom_rect(aes(
    xmin = XMin, xmax = XMax,
    ymin = YMin, ymax = YMax,
    fill = probe2
  )) +
  coord_equal()+
  theme_bw() +
  facet_wrap(~Sample)

ggsave(cell_max_quant_plot, filename = here(plot_dir, paste0("LHb2_cell_max_quant.png")), height = 9, width = 9)
ggsave(cell_max_quant_plot, filename = here(plot_dir, paste0("LHb2_cell_max_quant.pdf")), height = 9, width = 9)

#### Bin by 100 ####
halo_copies_rank <- halo_copies_long |>
  group_by(Sample, probe) |>
  arrange(-copies) |>
  mutate(copies_rank = rank(-copies)) |>
  mutate(rank_cut = cut(copies_rank, c(0, 100, 200, 300, 400)))

halo_copies_rank |> count(rank_cut)

halo_copies_rank |>
  filter(rank_cut == "(0,100]") |>
  arrange(copies) |>
  slice(1)

# Sample `Object Id`  XMin  XMax  YMin  YMax probe copies marker cluster probe2             copies_rank rank_cut
# <chr>        <dbl> <dbl> <dbl> <dbl> <dbl> <fct>  <dbl> <chr>  <chr>   <chr>                    <dbl> <fct>
#     1 Br6462        5170  6918  6967  6068  6106 520       26 POU4F1 Hb      520 POU4F1 (Hb)           93   (0,100]
# 2 Br6462        5641  5785  5815  6533  6557 570       25 CRH    LHb.2   570 CRH (LHb.2)           97.5 (0,100]
# 3 Br6462        5658  6133  6179  6647  6692 620       24 MCOLN3 LHb.3   620 MCOLN3 (LHb.3)        93.5 (0,100]
# 4 Br6462         281 13751 13782  1079  1107 690        4 ESRP1  LHb.6   690 ESRP1 (LHb.6)         99.5 (0,100]
# 5 Br8112        5334  6961  6999  7838  7830 520      109 POU4F1 Hb      520 POU4F1 (Hb)           99.5 (0,100]
# 6 Br8112        3997  2859  2892  8397  8381 570       15 CRH    LHb.2   570 CRH (LHb.2)           97   (0,100]
# 7 Br8112       11864  2325  2359  2840  2823 620       21 MCOLN3 LHb.3   620 MCOLN3 (LHb.3)        96.5 (0,100]
# 8 Br8112         185  2723  2751 12638 12630 690        7 ESRP1  LHb.6   690 ESRP1 (LHb.6)         78   (0,100]

halo_copies_rank |>
  filter(rank_cut == "(300,400]") |>
  arrange(copies) |>
  slice(1)

#### top100 confusion matrix ####
halo_copies_cat <- halo_copies_rank |>
    filter(!is.na(rank_cut)) |>
    filter(rank_cut == "(0,100]") |>
    ungroup() |>
    mutate(cat = paste0(probe, "_", marker, "_", rank_cut)) |>
    select(Sample, `Object Id`, cat)

halo_copies_cat2 <- halo_copies_cat |>
    left_join(halo_copies_cat, by = join_by(Sample, `Object Id`), relationship = "many-to-many")

confusion_top100 <- halo_copies_cat2 |>
    mutate(cat.x = factor(cat.x),
           cat.y = factor(cat.y),
           Sample = factor(Sample)) |>
    group_by(cat.x, cat.y, Sample, .drop = FALSE) |>
    summarise(n=n()) |>
    ggplot(aes(cat.x, cat.y, fill = n)) +
    geom_tile() +
    geom_text(aes(label = n), color = "white") +
    facet_wrap(~Sample) +
    # scale_fill_gradient(name = "count", trans = "log") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())

ggsave(confusion_top100, filename = here(plot_dir, "LHb2_confusion_top100.png"), height = 4, width = 8)
ggsave(confusion_top100, filename = here(plot_dir, "LHb2_confusion_top100.pdf"), height = 4, width = 8)


#### cell plots ####
copy_cut_colors <- c(`(300,400]` = "#FECC5C", `(200,300]` = "#FD8D3C", `(100,200]` = "#F03B20", `(0,100]` = "#BD0026")

halo_copies_rank_cut <- halo_copies_rank |>
  filter(!is.na(rank_cut)) |>
  ggplot() +
  geom_rect(aes(
    xmin = XMin, xmax = XMax,
    ymin = YMin, ymax = YMax,
    fill = rank_cut
  )) +
  scale_fill_manual(values = copy_cut_colors, "Copy Quantile") +
  coord_equal()+
  theme_bw() +
  facet_grid(Sample~probe2)

ggsave(halo_copies_rank_cut, filename = here(plot_dir, paste0("LHb2_cell_count_rank_cut_facet.png")), height = 4, width = 9)
ggsave(halo_copies_rank_cut, filename = here(plot_dir, paste0("LHb2_cell_count_rank_cut_facet.pdf")), height = 4, width = 9)

rank_cut_density <- halo_copies_rank |>
  ggplot(aes(x = copies, fill = rank_cut)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  scale_fill_manual(values = copy_cut_colors) +
  coord_cartesian(ylim=c(0, 100)) +
  facet_grid(Sample~probe2, scales = "free_x")

ggsave(rank_cut_density, filename = here(plot_dir, "LHb2_rank_cut_denisty.png"), height = 5)

## top 100
cell_rank_top100 <- halo_copies_rank |>
  filter(rank_cut == "(0,100]", probe != 520) |>
  ggplot() +
  geom_rect(aes(
    xmin = XMin, xmax = XMax,
    ymin = YMin, ymax = YMax,
    fill = probe2
  )) +
  coord_equal()+
  theme_bw() +
  facet_wrap(~Sample)

ggsave(cell_rank_top100, filename = here(plot_dir, "LHb2_rank_top100.png"))
ggsave(cell_rank_top100, filename = here(plot_dir, "LHb2_rank_top100.pdf"))

#### shadow plots ####
## if nuclei has >1 top100 ID, pick marker w/ max copies
halo_copies_rank_ID  <- halo_copies_rank |>
    filter(rank_cut == "(0,100]",
           probe != 520) |>
    group_by(Sample, `Object Id`) |>
    arrange(-copies) |>
    slice(1)

halo_copies_rank_ID |> count() |> filter(n>1)


halo_copies_rank_cut_shadow <- halo_copies_rank |>
    filter(probe == 520) |>
    ggplot() +
    geom_rect(aes(
        xmin = XMin, xmax = XMax,
        ymin = YMin, ymax = YMax,
        fill = copies > 10
    )) +
    geom_point(data = halo_copies_rank_ID,
               aes(x = XMax,
                   y = YMax,
                   color = probe2
               ), size = 0.7) +
    scale_fill_manual(values = c(`FALSE`="#CCCCCC80", `TRUE` = "black"), "POU4F1 Copy >10") +
    scale_color_manual(values = c("570 CRH (LHb.2)" = "#0096FF", ## cell type colors
                                  "620 MCOLN3 (LHb.3)" = "#89CFF0",
                                  "690 ESRP1 (LHb.6)" = "#008080"), "Top100 Nuclei") +
    coord_equal()+
    theme_void() +
    facet_wrap(~Sample)

ggsave(halo_copies_rank_cut_shadow, filename = here(plot_dir, paste0("LHb2_cell_count_rank_cut_facet_shadow.pdf")), height = 3, width = 9)

## for main fig
halo_copies_rank_cut_shadow_Br8112 <- halo_copies_rank |>
    filter(probe == 520,
           Sample == "Br8112") |>
    ggplot() +
    geom_rect(aes(
        xmin = XMin, xmax = XMax,
        ymin = YMin, ymax = YMax,
        fill = copies > 10
    )) +
    geom_point(data = halo_copies_rank_ID |>
               filter(Sample == "Br8112"),
               aes(x = XMax,
                   y = YMax,
                   color = probe2
               ), size = 0.7) +
    scale_fill_manual(values = c(`FALSE`="#CCCCCC80", `TRUE` = "black"), "POU4F1 Copy >10") +
    scale_color_manual(values = c("570 CRH (LHb.2)" = "#0096FF", ## cell type colors
                                  "620 MCOLN3 (LHb.3)" = "#89CFF0",
                                  "690 ESRP1 (LHb.6)" = "#008080"), "Top100 Nuclei") +
    coord_equal()+
    theme_void() +
    facet_wrap(~Sample)

ggsave(halo_copies_rank_cut_shadow_Br8112, filename = here(plot_dir, paste0("LHb2_cell_count_rank_cut_facet_shadow_Br8112.pdf")), height = 4, width = 4)

## just POU4F1 inset
adj = 15
halo_copies_rank_cut_shadowIN_Br8112 <- halo_copies_rank |>
    filter(probe == 520,
           Sample == "Br8112") |>
    ggplot() +
    geom_rect(aes(
        xmin = XMin-adj, xmax = XMax+adj,
        ymin = YMin-(adj*2), ymax = YMax+(adj*2),
        fill = copies > 10
    )) +
    scale_fill_manual(values = c(`FALSE`="#CCCCCC80", `TRUE` = "black")) +
    coord_equal()+
    theme_void() +
    theme(legend.position = "None")

ggsave(halo_copies_rank_cut_shadowIN_Br8112, filename = here(plot_dir, paste0("LHb2_cell_count_rank_cut_facet_shadowIN_Br8112.pdf")), height = 1, width = 1)


#### Export top objects ####

halo_copies_rank |> group_by(probe, Sample) |> filter(copies_rank <= 10) |> arrange(probe,copies_rank) |> write_csv(file = here("processed-data", "14_RNAscope", "HALO_data", "Lateral_exp2", "LHb_top10_nuclei.csv"))


