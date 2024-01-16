
library("tidyverse")
library("readxl")
library("GGally")
library("here")
library("sessioninfo")

## set up
plot_dir <- here("plots", "14_RNAscope", "04_LHb1_HALO")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load data
## local
list.files("~/OneDrive - Johns Hopkins/habenulaPilot-paper/HALO_data/Lateral_exp1/", recursive = TRUE)
# [1] "Br5422/LHb_Experiment1_Br5422_middleslice_job1879_object_RESULTS.csv"
# [2] "Br6462/LHb_Experiment1_Br6462_rightslice_20x_job1886_object_RESULTS.csv"
# [3] "Br8112/LHbExperiment1_Br8112_thirdtoleft_20x_job1889_object_RESULTS.csv"
## TODO fix to here notation/JHPCE paths
fn <- list(Br5422 = "~/OneDrive - Johns Hopkins/habenulaPilot-paper/HALO_data/Lateral_exp1/Br5422/LHb_Experiment1_Br5422_middleslice_job1879_object_RESULTS.csv",
           Br6462 = "~/OneDrive - Johns Hopkins/habenulaPilot-paper/HALO_data/Lateral_exp1/Br6462/LHb_Experiment1_Br6462_rightslice_20x_job1886_object_RESULTS.csv",
           Br8112 = "~/OneDrive - Johns Hopkins/habenulaPilot-paper/HALO_data/Lateral_exp1/Br8112/LHbExperiment1_Br8112_thirdtoleft_20x_job1889_object_RESULTS.csv"
           )

map(fn, file.exists)
# halo <- map2(fn, names(fn), ~read_excel(.x) |> mutate(Sample = .y))
halo <- map2(fn, names(fn), ~read_csv(.x) |> mutate(Sample = .y))

map_int(halo, ncol)
map_int(halo, nrow)

## flip y for Br8112
halo$Br8112$YMax <- max(halo$Br8112$YMax) - halo$Br8112$YMax
halo$Br8112$YMin <- max(halo$Br8112$YMax) - halo$Br8112$YMin

## test
# halo$Br8112 |>
#     ggplot() +
#     geom_rect(aes(
#         xmin = XMin, xmax = XMax,
#         ymin = YMin, ymax = YMax
#     ))

## rbind tables
halo <- do.call("rbind", halo)
colnames(halo) <- gsub(" \\(µm.*\\)", "", colnames(halo))

experiment <- tibble(probe = factor(c(690, 620, 570, 520)),
                     marker = c("ONECUT2", "TLE2","SEMA3D","POU4F1"),
                     cluster = c("LHb.1", "LHb.4", "LHb.5/1", "Hb")) |>
    mutate(probe2 = paste0(probe, " ", marker, " (", cluster,")"))
# A tibble: 4 × 4
# probe marker  cluster probe2
# <fct> <chr>   <chr>   <chr>
# 1 690   ONECUT2 LHb.1   690 ONECUT2 (LHb.1)
# 2 620   TLE2    LHb.4   620 TLE2 (LHb.4)
# 3 570   SEMA3D  LHb.5/1 570 SEMA3D (LHb.5/1)
# 4 520   POU4F1  Hb      520 POU4F1 (Hb)

# write.csv(experiment, file = here("processed-data", "14_RNAscope", "HALO_data", "Lateral_exp2", "Probes_LHb1.csv"))

#### create halo long ####
halo_copies_long <- halo |>
  select(Sample, `Object Id`, XMin, XMax, YMin, YMax, ends_with("Copies")) |>
  pivot_longer(!c(Sample, `Object Id`, XMin, XMax, YMin, YMax), names_to = "probe", values_to = "copies") |>
  mutate(probe = factor(gsub("^2.*Opal (\\d+) Copies", "\\1", probe))) |>
  left_join(experiment)
  # mutate(probe = extract_numeric(probe))

#### ggpairs ####
halo_copies_wide <- halo_copies_long |>
    select(Sample, `Object Id`, copies, probe2) |>
    pivot_wider(id_cols = c(Sample, `Object Id`), names_from = probe2, values_from = copies)

# copies_tab |> count(sum == 0)
# copies_tab |> count(Sample, sum == 0)

gg_copies <- halo_copies_wide |> ggpairs(3:6, mapping=ggplot2::aes(colour = Sample))
ggsave(gg_copies, filename = here(plot_dir, "LHb1_gg_copies.png"), height = 12, width = 12)


#### quantile grouping ####
# get quantile values for non-zero values
copy_quant_cutoff <- halo_copies_long |>
  filter(copies != 0) |>
  group_by(probe, Sample) |>
  summarize(q25 = quantile(copies,probs = 0.25),
            q50 = quantile(copies,probs = 0.5),
            q75 = quantile(copies,probs = 0.75),
            q95 = quantile(copies,probs = 0.95))

# probe Sample   q25   q50   q75   q95
# <fct> <chr>  <dbl> <dbl> <dbl> <dbl>
# 1 520   Br5422     1     2     9  59
# 2 520   Br6462     1     2     6 129
# 3 520   Br8112     1     2     4  22

## use to classify cells
halo_copies_long_quant <- halo_copies_long |>
  left_join(copy_quant_cutoff) |>
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

# probe2               Sample     None    q0   q25   q50   q75
# <chr>                <chr>  <int> <int> <int> <int> <int>
#     1 520 POU4F1 (Hb)      Br5422 12325  1573   480   664   857
# 2 570 SEMA3D (LHb.5/1) Br5422  9847  2762  1258   654  1378
# 3 620 TLE2 (LHb.4)     Br5422 12234  1829   675   342   819
# 4 690 ONECUT2 (LHb.1)  Br5422 14414   902    NA   220   363


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

ggsave(cell_count_quant, filename = here(plot_dir, paste0("LHb1_cell_count_quant_facet.png")), height = 6, width = 9)
ggsave(cell_count_quant, filename = here(plot_dir, paste0("LHb1_cell_count_quant_facet.pdf")), height = 6, width = 9)


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

ggsave(copies_density, filename = here(plot_dir, "LHb1_copies_denisty.png"), height = 6, width = 9)

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

ggsave(hex_copies_median, filename = here(plot_dir, paste0("LHb1_hex_copies_median_facet.png")), height = 6, width = 9)

## TODO add limit to legend title
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
    ) + coord_equal() +
    theme_bw() +
    facet_grid(Sample~probe2) +
    theme(legend.position = "bottom")

ggsave(hex_copies_max, filename = here(plot_dir, paste0("LHb1_hex_copies_max_facet.png")), height = 6, width = 9)
ggsave(hex_copies_max, filename = here(plot_dir, paste0("LHb1_hex_copies_max_facet.pdf")), height = 7, width = 7)

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

ggsave(cell_max_quant_plot, filename = here(plot_dir, paste0("LHb1_cell_max_quant.png")), height = 9, width = 9)
ggsave(cell_max_quant_plot, filename = here(plot_dir, paste0("LHb1_cell_max_quant.pdf")), height = 9, width = 9)

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

# Sample    `Object Id`  XMin  XMax  YMin  YMax probe copies marker  cluster probe2          copies_rank rank_cut
# <chr>        <dbl> <dbl> <dbl> <dbl> <dbl> <fct>  <dbl> <chr>   <chr>   <chr>                 <dbl> <fct>
#     1 Br5422        6769  1925  1961  8481  8513 520       79 POU4F1  Hb      520 POU4F1 (Hb)       100   (0,100]
# 2 Br5422        1589  7922  7945  3041  3070 570       15 SEMA3D  LHb.5/1 570 SEMA3D (LH…        99.5 (0,100]
# 3 Br5422        1053  6364  6396  2262  2283 620       25 TLE2    LHb.4   620 TLE2 (LHb.…        96.5 (0,100]
# 4 Br5422         349  8477  8497   675   699 690       10 ONECUT2 LHb.1   690 ONECUT2 (L…        99   (0,100]

halo_copies_rank |>
  filter(rank_cut == "(300,400]") |>
  arrange(copies) |>
  slice(1)

#### confusion matrix ###

halo_copies_rank_wide <- halo_copies_rank |>
    ungroup() |>
    select(Sample, `Object Id`, marker, rank_cut) |>
    pivot_wider(id_cols = c(Sample, `Object Id`), names_from = marker, values_from = rank_cut)

confusion_SEMA3D_ONECUT2 <- halo_copies_rank_wide |>
    group_by(Sample) |>
    count(SEMA3D, ONECUT2) |>
    ggplot(aes(SEMA3D, ONECUT2, fill = n)) +
    geom_tile() +
    geom_text(aes(label = n)) +
    facet_wrap(~Sample) + scale_fill_gradient(name = "count", trans = "log")

ggsave(confusion_SEMA3D_ONECUT2, filename = here(plot_dir, "confusion_SEMA3D_ONECUT2.png"), height = 5, width = 9)

## better confusion
halo_copies_cat <- halo_copies_rank |>
    filter(!is.na(rank_cut)) |>
    filter(rank_cut == "(0,100]") |>
    ungroup() |>
    mutate(cat = paste0(probe, "_", marker, "_", rank_cut)) |>
    select(Sample, `Object Id`, cat)

halo_copies_cat2 <- halo_copies_cat |>
    left_join(halo_copies_cat, by = join_by(Sample, `Object Id`), relationship = "many-to-many")


halo_copies_cat2 |>
    mutate(cat.x = factor(cat.x),
           cat.y = factor(cat.y),
           Sample = factor(Sample)) |>
    group_by(cat.x, cat.y, Sample, .drop = FALSE) |>
    summarise(n=n()) |>
    arrange(n)

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

ggsave(confusion_top100, filename = here(plot_dir, "LHb1_confusion_top100.png"), height = 4, width = 9)
ggsave(confusion_top100, filename = here(plot_dir, "LHb1_confusion_top100.pdf"), height = 4, width = 9)

## copies scatter
halo_copies2_cat <- halo_copies_rank |>
    filter(!is.na(rank_cut)) |>
    ungroup() |>
    select(Sample, `Object Id`, probe2, copies) |>
    left_join(halo_copies_cat)

halo_copies2_cat |> count(cat)

copies_vs_top100_jitter <- halo_copies2_cat |>
    ggplot(aes(x = probe2, y = copies, color = cat)) +
    geom_jitter(alpha = .5) +
    theme_bw() +
    theme(legend.position = "bottom")

ggsave(copies_vs_top100_jitter, filename = here(plot_dir, "LHb1_copies_vs_top100_jitter.png"), width = 8)

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
  facet_wrap(Sample~probe2)

ggsave(halo_copies_rank_cut, filename = here(plot_dir, paste0("LHb1_cell_count_rank_cut_facet.png")), height = 5, width = 9)
ggsave(halo_copies_rank_cut, filename = here(plot_dir, paste0("LHb1_cell_count_rank_cut_facet.pdf")), height = 5, width = 9)


rank_cut_density <- halo_copies_rank |>
  # filter(Sample == "Sample1862") |>
  # filter(copies != 0) |>
  ggplot(aes(x = copies, fill = rank_cut)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  scale_fill_manual(values = copy_cut_colors) +
  coord_cartesian(ylim=c(0, 100)) +
  facet_grid(Sample~probe2, scales = "free_x")

ggsave(rank_cut_density, filename = here(plot_dir, "LHb1_rank_cut_denisty.png"), height = 5, width = 9)


## top 100
cell_rank_top100 <- halo_copies_rank |>
  filter(rank_cut == "(0,100]") |>
  ggplot() +
  geom_rect(aes(
    xmin = XMin, xmax = XMax,
    ymin = YMin, ymax = YMax,
    fill = probe2
  )) +
  coord_equal()+
  theme_bw() +
  facet_wrap(~Sample)

ggsave(cell_rank_top100, filename = here(plot_dir, "LHb1_rank_top100.png"), height = 5, width = 9)
ggsave(cell_rank_top100, filename = here(plot_dir, "LHb1_rank_top100.pdf"), height = 5, width = 9)

cell_rank_top100_point <-halo_copies_rank |>
    filter(rank_cut == "(0,100]")  |>
    ggplot() +
    geom_point(aes(x = XMax,
                   y = YMax,
                   color = probe2
               ), size = 0.7) +
    coord_equal()+
    theme_bw()

ggsave(cell_rank_top100_point, filename = here(plot_dir, "LHb1_rank_top100_point.png"))

cell_rank_top200 <- halo_copies_rank |>
    filter(rank_cut %in% c("(0,100]", "(100,200]")) |>
    ggplot() +
    geom_rect(aes(
        xmin = XMin, xmax = XMax,
        ymin = YMin, ymax = YMax,
        fill = probe2
    )) +
    coord_equal()+
    theme_bw() +
    facet_wrap(~Sample)

ggsave(cell_rank_top200, filename = here(plot_dir, "LHb1_rank_top200.png"))
ggsave(cell_rank_top200, filename = here(plot_dir, "LHb1_rank_top200.pdf"))


color_official_markers = c(
    "LHb.1" = c("#0085af"),
    "LHb.4" = c("#6F8FAF"),
    "LHb.5" = c("#40E0D0")
)

#### Shadow plots ####

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
    scale_color_manual(values = c("690 ONECUT2 (LHb.1)" = "#0085af", ## cell type colors
                                 "620 TLE2 (LHb.4)" = "#6F8FAF",
                                 "570 SEMA3D (LHb.5/1)" = "#40E0D0"), "Top100 Nuclei") +
    # scale_color_manual(values = c("690 ONECUT2 (LHb.1)" = "red",
    #                              "620 TLE2 (LHb.4)" = "blue",
    #                              "570 SEMA3D (LHb.5/1)" ="orange"), "Top100 Nuclei") +
    coord_equal()+
    theme_void() +
    facet_wrap(~Sample)

ggsave(halo_copies_rank_cut_shadow, filename = here(plot_dir, paste0("LHb1_cell_count_rank_cut_facet_shadow.pdf")), height = 3, width = 9)

# neon purple "#B026FF"

# halo_copies_rank_cut_shadow2 <- halo_copies_rank |>
#     filter(probe == 520) |>
#     ggplot() +
#     geom_rect(aes(
#         xmin = XMin, xmax = XMax,
#         ymin = YMin, ymax = YMax,
#         fill = copies > 10
#     )) +
#     geom_point(data = halo_copies_rank |>
#                    filter(rank_cut == "(0,100]",
#                           probe != 520),
#                aes(x = XMax,
#                    y = YMax,
#                    fill = probe2
#                ), shape = 21,
#                color = "black",
#                size = 1.2) +
#     # scale_fill_manual(values = c(`FALSE`="#CCCCCC80",
#     #                              `TRUE` = "magenta",
#     #                              "690 ONECUT2 (LHb.1)" = c("#0085af"),
#     #                              "620 TLE2 (LHb.4)" = c("#6F8FAF"),
#     #                              "570 SEMA3D (LHb.5/1)" = c("#40E0D0")), "Top100 Nuclei") +
#     scale_fill_manual(values = c(`FALSE`="#CCCCCC80",
#                                  `TRUE` = "magenta",
#                                  "690 ONECUT2 (LHb.1)" = c("#0085af"),
#                                  "620 TLE2 (LHb.4)" = c("#6F8FAF"),
#                                  "570 SEMA3D (LHb.5/1)" = c("#40E0D0")), "Top100 Nuclei") +
#     coord_equal()+
#     theme_void() +
#     facet_wrap(~Sample)
#
# # ggsave(halo_copies_rank_cut_shadow2, filename = here(plot_dir, paste0("LHb1_cell_count_rank_cut_facet_shadow.png")), height = 6, width = 9)
# ggsave(halo_copies_rank_cut_shadow2, filename = here(plot_dir, paste0("LHb1_cell_count_rank_cut_facet_shadow2.pdf")), height = 6, width = 10)

## main figure plot
halo_copies_rank_cut_shadow_Br5422 <- halo_copies_rank |>
    filter(probe == 520,
           Sample == "Br5422") |>
    ggplot() +
    geom_rect(aes(
        xmin = XMin, xmax = XMax,
        ymin = YMin, ymax = YMax,
        fill = copies > 10
    )) +
    geom_point(data = halo_copies_rank_ID |>
                   filter(Sample == "Br5422"),
               aes(x = XMax,
                   y = YMax,
                   color = probe2
               ), size = 0.7) +
    scale_fill_manual(values = c(`FALSE`="#CCCCCC80", `TRUE` = "black"), "POU4F1 Copy >10") +
    scale_color_manual(values = c("690 ONECUT2 (LHb.1)" = "#0085af", ## cell type colors
                                  "620 TLE2 (LHb.4)" = "#6F8FAF",
                                  "570 SEMA3D (LHb.5/1)" = "#40E0D0"), "Top100 Nuclei") +
    coord_equal()+
    theme_void() +
    facet_wrap(~Sample)

ggsave(halo_copies_rank_cut_shadow_Br5422, filename = here(plot_dir, paste0("LHb1_cell_count_rank_cut_facet_shadow_Br5422.pdf")), height = 4, width = 4)

halo_copies_rank |>
    filter(rank_cut == "(0,100]",
           probe != 520,
           Sample == "Br5422") |> count(probe)

#### Export top objects ####

halo_copies_rank |> group_by(probe, Sample) |> filter(copies_rank <= 10) |> arrange(probe,copies_rank) |> write_csv(file = here("processed-data", "14_RNAscope", "HALO_data", "Lateral_exp1", "LHb1_top10_nuclei.csv"))


