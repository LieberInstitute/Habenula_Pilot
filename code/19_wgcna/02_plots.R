library(here)
library(tidyverse)
library(SummarizedExperiment)
library(WGCNA)
library(jaffelab)
library(sessioninfo)

net_path = here('processed-data', '19_wgcna', 'modules.rds')
rse_path = here(
    'processed-data', 'rse_objects', 'rse_gene_Habenula_Pilot.rda'
)
plot_dir = here('plots', '19_wgcna')

set.seed(0)
dir.create(plot_dir, showWarnings = FALSE)

net = readRDS(net_path)
rse_gene = get(load(rse_path))

me_df = net$MEs |>
    rownames_to_column('RNum') |>
    as_tibble() |>
    left_join(
        colData(rse_gene) |>
            as_tibble() |>
            select(RNum, PrimaryDx),
        by = 'RNum'
    ) |>
    mutate(
        PrimaryDx = factor(
            ifelse(PrimaryDx == "Schizo", "SCZD", "Control"),
            levels = c("Control", "SCZD")
        )
    )

plot_list = list()
for (i in seq_len(ncol(me_df) - 2)) {
    #   Get p-value of linear relationship with diagnosis
    lin_mod = lm(
        as.formula(paste0('ME', i - 1, ' ~ PrimaryDx')), data = me_df
    )
    p_val = signif(summary(lin_mod)$coef['PrimaryDxSCZD', 4], 2)
    
    plot_list[[i]] = ggplot(
            me_df,
            aes(
                x = PrimaryDx, y = !!sym(paste0('ME', i - 1)), color = PrimaryDx
            )
        ) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter() +
        geom_text(
            label = paste('\np =', p_val, ''),
            x = Inf, y = Inf, hjust = 1, vjust = 1, color = 'black', size = 8
        ) +
        guides(color = "none") +
        labs(x = "Diagnosis", y = paste("Module", i - 1)) +
        theme_bw(base_size = 20)
}

pdf(file.path(plot_dir, 'modules_by_dx.pdf'))
print(plot_list)
dev.off()

session_info()
