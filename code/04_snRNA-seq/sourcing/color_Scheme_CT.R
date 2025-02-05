# 00_Color_Scheme_ct.R 
# 2/28/23 - Bukola Ajanaku (advised and mentored by Louise Huuki-Myers)
## This creates a custom color scheme based on colors grabbed from 
## https://medialab.github.io/iwanthue/ . This operates as a function that returns 
## a list of HEX colors that can change lengths based on number of colors needed.
## Max amount for now will be 20.


library("here")
library("jaffelab")
library("RColorBrewer")

# updated color pallete (for 37 max)
color_pallete <- c( "#c16600",
                    "#6d36bb",
                    "#869343",
                    "#d933b2",
                    "#01bf6a",
                    "#ff7df4",
                    "#c0b200",
                    "#997aff",
                    "#9ad771",
                    "#951581",
                    "#00782c",
                    "#ec97ff",
                    "#d85c3e",
                    "#0157bb",
                    "#ee740d",
                    "#5c9fff",
                    "#a31907",
                    "#03c5e5",
                    "#fd3558",
                    "#6088c8",
                    "#ff4b8c",
                    "#006647",
                    "#bb0062",
                    "#3690c3",
                    "#793b76",
                    "#f9ba6f",
                    "#4bc0ff",
                    "#933300",
                    "#c2b1ff",
                    "#845c00",
                    "#eab3fa",
                    "#ff9d78",
                    "#3e5187",
                    "#db9b84",
                    "#9972c0",
                    "#ab7297",
                    "#bd555f")

# Previous color pallete (for 20 max)
# color_pallete <- c("#5b679a",
#   "#c2cf1a",
#   "#d779ff",
#   "#009421",
#   "#78319a",
#   "#00af80",
#   "#ff4b7f",
#   "#0091e1",
#   "#d45a00",
#   "#919fff",
#   "#ffb050",
#   "#faa0ff",
#   "#525b00",
#   "#ff64b8",
#   "#806200",
#   "#ffa5c7",
#   "#91332d",
#   "#ff886b",
#   "#8b3160",
#   "#b8766a")

grabColors <- function(num, start = 1){
  # num is the total number of colors needed
  # start defaults at the first color but can be re-written to start elsewhere
  return(color_pallete[start:(start + num - 1)])
}


## How Louise wanted me to do this:

# cellTypecolors_9 <- 
#   c(ExcNeuron = "#274ea6",
#     InhibThalMed = "#6ed94f",
#     OPC = "#ff67e3",
#     Astro = "#296f00",
#     Microglia = "#cb0029",
#     ThalMedPVT = "#009770",
#     Oligo = "#ff815f",
#     ExcThalMed = "#d794b4",
#     LHbN = "#ffb848")

# preview_colors <- function(cell_colors) {
#   par(las = 2) # make label text perpendicular to axis
#   par(mar = c(5, 8, 4, 2)) # increase y-axis margin.
#   barplot(rep(1, length(cell_colors)),
#           col = cell_colors,
#           horiz = TRUE,
#           axes = FALSE,
#           names.arg = names(cell_colors)
#   )
# }
# 
# png(here("plots", "cell_colors", "cellTypecolors_9.png"), height = 800)
#   preview_colors(cellTypecolors_9)
# dev.off()

# # saving colors
# save(cellTypecolors_9, file = here("processed-data", "cell_type_colors.Rdata"))
