library(ggplot2)
library(patchwork)
library("rlang")
#removed   cutoff = 0.1 from below functionf
broad_boxplot <- function(data, xvar, yvar, fillvar, colorvar , pallet, title, subtitle){
  ignore <- NULL
  
  data <- data %>% group_by(cell_type)
  
  bp_broad <- data %>% 
    ggplot(aes(x = !!sym(xvar), y = !!sym(yvar), fill = !!sym(fillvar), color = !!sym(colorvar))) +
    geom_boxplot()+
    scale_fill_manual(values = pallet)+
    scale_color_manual(values = pallet)+
    labs(title = title,
         subtitle = subtitle)+
    theme(legend.position = "right")
  
#  bp_little <- data %>% filter(!big) %>%
#    ggplot(aes(x = !!sym(xvar), y = !!sym(yvar), fill = !!sym(fillvar), color = !!sym(colorvar))) +
#    geom_boxplot()+
#    scale_fill_manual(values = pallet) +
#   scale_color_manual(values = pallet)
  
  return(bp_broad)
}


#both_regions_bl_boxplot <- function(data_sacc, data_amyg, yvar,
#                                    subtitle_sacc = "sACC samples",
#                                    subtitle_amyg = "Amygdala samples",
#                                    title,
#                                    filename){
  
#  blb_sacc <- big_little_boxplot(data_sacc, xvar = "cell_type", yvar = yvar,
#                                 fillvar =  "PrimaryDx",
#                                 title = title,
#                                 subtitle = subtitle_sacc)
#  
#  blb_amyg <- big_little_boxplot(data_amyg, xvar = "cell_type", yvar = yvar,
#                                 fillvar =  "PrimaryDx",
#                                 title = NULL,
#                                 subtitle = subtitle_amyg)
#  
#  ggsave(plot = blb_sacc / blb_amyg, filename = filename, width = 10)
#}
