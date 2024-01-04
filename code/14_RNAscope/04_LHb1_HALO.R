
library("tidyverse")
library("readxl")
library("GGally")
library("here")
library("sessioninfo")

## set up 
plot_dir <- here("plots", "14_RNAscope", "03_LHb2_HALO")
if(!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)

## load data
fn <- list(Br6462 = here("processed-data", "14_RNAscope", "HALO_data", "Lateral_exp2",
                         "LHbExperiment_Br6462_secondpass.xls")
            )

map(fn, file.exists)
halo <- map2(fn, names(fn), ~read_excel(.x) |> mutate(job = .y))

# all.equal(halo$job1862 |> select(-job), halo$job1864|> select(-job))

map(halo, dim)
# $Br6462
# [1] 16383    55

## rbind tables
halo <- do.call("rbind", halo)
colnames(halo) <- gsub(" \\(µm.*\\)", "", colnames(halo))

#### create halo long ####
experiment <- tibble(probe = factor(c(690, 620, 570, 520)),
       marker = c("CCK", "EBF3","CHAT","POU4F1"),
       cluster = c("LHb.1", "LHb.3", "LHb.2", "Hb"))
# probe marker cluster
# <dbl> <chr>  <chr>  
# 1   690 CCK    LHb.1  
# 2   620 EBF3   LHb.3  
# 3   570 CHAT   LHb.2  
# 4   520 POU4F1 Hb 

# write.csv(experiment, file = here("processed-data", "14_RNAscope", "HALO_data", "Lateral_exp2", "Probes_LHb2.csv"))

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
#.  probe job      q25   q50   q75   q95
# 1 520   Br6462     1     3    10    25
# 2 570   Br6462     1     1     3    18
# 3 620   Br6462     1     2     4    16
# 4 690   Br6462     1     1     2     4

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

# probe2           job     None    q0   q25   q50   q75
# <chr>            <chr>  <int> <int> <int> <int> <int>
#   1 520 POU4F1 (Hb)  Br6462 14212   637   460   570   504
# 2 570 CHAT (LHb.2) Br6462 13370  1718    NA   692   603
# 3 620 EBF3 (LHb.3) Br6462 12618  1683   671   506   905
# 4 690 CCK (LHb.1)  Br6462 14803  1032    NA   310   238


#### plot with cell_quant 
copy_quant_colors <- c(None = "grey75", q0 = "#FECC5C", q25 = "#FD8D3C", q50 = "#F03B20", q75 = "#BD0026")

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
  facet_wrap(job~probe2)

ggsave(cell_count_quant, filename = here(plot_dir, paste0("LHb2_cell_count_quant_facet.png")), height = 7, width = 9)
ggsave(cell_count_quant, filename = here(plot_dir, paste0("LHb2_cell_count_quant_facet.pdf")), height = 7, width = 9)


#### distribution
copies_density <- halo_copies_long2 |>
  # filter(job == "job1862") |>
  # filter(copies != 0) |>
  ggplot(aes(x = copies, fill = copy_quant)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  scale_fill_manual(values = copy_quant_colors, "Copy Quantile") +
  coord_cartesian(ylim=c(0, 2000), xlim = c(-1, 50)) +
  facet_grid(job~probe2)

ggsave(copies_density, filename = here(plot_dir, "LHb2_copies_denisty.png"), height = 5)

#### hex plots ####

hex_copies_median <- ggplot(halo_copies_long2) +
  stat_summary_hex(aes(x = XMax, y = YMax, z = copies),
                   fun = median, bins = 100
  ) +
  # scale_fill_continuous(type = "viridis") + ## top value of 100 for visualization 
  scale_fill_continuous(type = "viridis", limits = c(1,50), "Median Copies\n(max 50)") + ## top value of 100 for visualization
  coord_equal() +
  theme_bw() +
  facet_grid(job~probe2)

ggsave(hex_copies_median, filename = here(plot_dir, paste0("LHb2_hex_copies_median_facet.png")), height = 4, width = 12)


## max quant cell_max_quant <- halo_copies_long2 |>
cell_max_quant <- halo_copies_long2 |>
  group_by(`Object Id`) |>
  filter(probe != "520", 
         copies != 0) |>
  arrange(-copies) |>
  slice(1)

cell_max_quant |> group_by(probe2, job) |> count()
cell_max_quant |> group_by(probe2, job) |> count(copy_quant)

cell_max_quant |>  
  filter(copy_quant >= "q50") |>
  group_by(probe2, job) |> 
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
  facet_wrap(~job)

ggsave(cell_max_quant_plot, filename = here(plot_dir, paste0("LHb2_cell_max_quant.png")), height = 9, width = 9)
ggsave(cell_max_quant_plot, filename = here(plot_dir, paste0("LHb2_cell_max_quant.pdf")), height = 9, width = 9)

#### Examine the relationship of counts in larger field of view ####

copies <- halo |> 
  select(ends_with("Copies")) |>
  filter(rowSums(across(where(is.numeric))) > 0)

gg_copies <- copies |>
  ggpairs(aes(alpha = 0.5)) 

ggsave(gg_copies, filename = here(plot_dir, "LHb2_gg_copies.png"), height = 12, width = 12)


#### Bin by 100 ####

halo_copies_rank <- halo_copies_long |>
  group_by(probe) |>
  arrange(-copies) |>
  mutate(copies_rank = rank(-copies)) |>
  mutate(rank_cut = cut(copies_rank, c(0, 100, 200, 300, 400)))
  
halo_copies_rank |> count(rank_cut)

halo_copies_rank |> 
  filter(rank_cut == "(0,100]") |>
  arrange(copies) |>
  slice(1)

# job    `Object Id`  XMin  XMax  YMin  YMax probe copies marker cluster probe2           copies_rank rank_cut
# <chr>        <dbl> <dbl> <dbl> <dbl> <dbl> <fct>  <dbl> <chr>  <chr>   <chr>                  <dbl> <fct>   
#   1 Br6462        5170  6918  6967  6068  6106 520       26 POU4F1 Hb      520 POU4F1 (Hb)         93   (0,100] 
# 2 Br6462        5641  5785  5815  6533  6557 570       25 CHAT   LHb.2   570 CHAT (LHb.2)        97.5 (0,100] 
# 3 Br6462        5658  6133  6179  6647  6692 620       24 EBF3   LHb.3   620 EBF3 (LHb.3)        93.5 (0,100] 
# 4 Br6462         281 13751 13782  1079  1107 690        4 CCK    LHb.1   690 CCK (LHb.1)         99.5 (0,100] 

halo_copies_rank |> 
  filter(rank_cut == "(300,400]") |>
  arrange(copies) |>
  slice(1)

# job    `Object Id`  XMin  XMax  YMin  YMax probe copies marker cluster probe2           copies_rank rank_cut 
# <chr>        <dbl> <dbl> <dbl> <dbl> <dbl> <fct>  <dbl> <chr>  <chr>   <chr>                  <dbl> <fct>    
# 1 Br6462        5201  7032  7061  6168  6194 520       13 POU4F1 Hb      520 POU4F1 (Hb)         397  (300,400]
# 2 Br6462        5597  5343  5379  6653  6683 570        6 CHAT   LHb.2   570 CHAT (LHb.2)        379  (300,400]
# 3 Br6462         354 12362 12397  1016  1055 620       10 EBF3   LHb.3   620 EBF3 (LHb.3)        394. (300,400]
# 4 Br6462          14 12680 12717   241   276 690        2 CCK    LHb.1   690 CCK (LHb.1)         394. (300,400]

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
  facet_wrap(job~probe2)

ggsave(halo_copies_rank_cut, filename = here(plot_dir, paste0("LHb2_cell_count_rank_cut_facet.png")), height = 7, width = 9)
ggsave(halo_copies_rank_cut, filename = here(plot_dir, paste0("LHb2_cell_count_rank_cut_facet.pdf")), height = 7, width = 9)


rank_cut_density <- halo_copies_rank |>
  # filter(job == "job1862") |>
  # filter(copies != 0) |>
  ggplot(aes(x = copies, fill = rank_cut)) +
  geom_histogram(binwidth = 1) +
  theme_bw() +
  scale_fill_manual(values = copy_cut_colors) +
  coord_cartesian(ylim=c(0, 100)) +
  facet_wrap(~probe2, scales = "free_x")

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
  facet_wrap(~job)

ggsave(cell_rank_top100, filename = here(plot_dir, "LHb2_rank_top100.png"))
ggsave(cell_rank_top100, filename = here(plot_dir, "LHb2_rank_top100.pdf"))
