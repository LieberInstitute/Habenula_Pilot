suppressMessages({
  library(GGally)
  library(data.table)
  library(stringi)
  library(ggtext)
  library(glue)
})

xggPairData <- NULL
xggColumnLabels <- NULL
# xggpairs will update these globals with main metrics for all the pairs
# and their cor and CAT values for the pairs plotted. used by gatherPairData()

## Function to gather pair data and augment gatheredData
## usage:
## allggData <- gatherPairData(xggPairData, allggData)
gatherPairData <- function(xggPairData, gatheredData) {
  if (is.null(gatheredData)) {
    gatheredData <- xggPairData[0, ]  # Create empty dataframe with same columns
  }
  # Iterate over each row in xggPairData
  for (i in 1:nrow(xggPairData)) {
    row <- xggPairData[i, ]

    # Convert labels to character to ensure proper sorting
    label1 <- as.character(row$x_label)
    label2 <- as.character(row$y_label)

    # Sort the labels lexicographically
    sorted_labels <- sort(c(label1, label2))
    new_x_label <- sorted_labels[1]
    new_y_label <- sorted_labels[2]

    # Check if the sorted pair already exists in corCatData
    existing_indices <- which(
      gatheredData$x_label == new_x_label &
      gatheredData$y_label == new_y_label
    )

    if (length(existing_indices) > 0) {
      # Extract existing cor_value and cat_perc
      existing_cor_value <- gatheredData$cor_value[existing_indices]
      existing_cat_perc <- gatheredData$cat_perc[existing_indices]

      # Compare with the current row's cor_value and cat_perc
      if (!all(existing_cor_value == row$cor_value) ||
          !all(existing_cat_perc == row$cat_perc)) {
        stop(paste(
          "Duplicate pair with differing values for",
          new_x_label, "and", new_y_label
        ))
      }
      # If identical, do nothing (retain existing entry)
    } else {
      # Create a new row with sorted labels and corresponding values
      new_row <- row
      new_row$x_label <- new_x_label
      new_row$y_label <- new_y_label
      # Set x and y to NA as they may vary across plots
      new_row$x <- NA
      new_row$y <- NA
      # Append the new row to corCatData
      gatheredData <- rbind(gatheredData, new_row)
    }
  }
  return(gatheredData)  ## return the updated gatheredData
}



## convert a named list of DE results to wide data for xggpairs
## each element of the list is a data.frame with topTable output
## from limma. The ensemblID is extracted from the rownames if not already present
deLst_ggwide <- function(deres_lst) {
  deres_lst <- lapply(deres_lst, function(x) {
    if (is.null(x$ensemblID)) { ## remove versioning from ensemblIDs
      x$ensemblID  <- sub(".*(ENS[TG]\\d+)\\.\\d+(.*)$", "\\1\\2", rownames(x), perl = TRUE)
    }
    return(x)
  })
  names(deres_lst) <- NULL
  combined_data <- rbindlist(deres_lst, idcol = "source")
  wide_data <- dcast(combined_data,
                     ensemblID ~ source,
                     value.var = c("t", "adj.P.Val", "P.Value", "logFC"))
  setDF(wide_data) ## ggpairs uses data.frame internally anyway
  rownames(wide_data) <- wide_data$ensemblID
  return(wide_data)
}

## adapted from FFPE::CATplot()
ggCATplot <- function(vec1, vec2, maxrank = 3000, x_col=0, y_col=0, fdr_overlap=NULL,
                pval_overlap=NULL, no.p=FALSE, xlab=NULL, ylab=NULL,
                standalone=FALSE) {
  # Sort vectors and get names
  vec1 <- names(sort(vec1))
  vec2 <- names(sort(vec2))
  # Adjust maxrank if needed
  maxrank <- min(maxrank, length(vec1), length(vec2))
  abl_slope=1/min(length(vec1), length(vec2))
  ranks <- 1:maxrank
  concordance <- sapply(ranks, function(i) length(intersect(vec1[1:i], vec2[1:i])) / i)
  pldata <- data.table(rank = ranks, concordance = concordance)
  # Calculate y-axis breaks
  max_y <- max(pldata$concordance)
  showdec <- "%.1f"
  #cat(">>> max_y:", max_y,"\n")
  #cat("    ybreaks before: ")
  #print(y_breaks)
  y_breaks <- seq(0, 1, by = 0.2)
  ## y_breaks_label is y_breaks without last element
  y_breaks_label <- y_breaks #y_breaks[-length(y_breaks)]
  line_color <- "gray45"
  if (max_y > 0.5) {
    line_color <- "#208030"
  } else if (max_y>0.2) { ## between 0.2 and 0.5
    #y_breaks <- seq(0, max_y, by = 0.1)  # Default
    line_color <- "#536857"
  }
  max_y <- fifelse(standalone, 1.05, 1.16)
  min_y <- fifelse(standalone, -0.07, -0.14)
  final_perc <- round(pldata$concordance[maxrank] * 100, 1)

  # Define consistent x-axis limits
  x_left_extension <- fifelse(standalone, maxrank*0.06, maxrank*0.25)
  #x_right_extension <- maxrank * 0.05
  x_right_extension <- 0
  x_limits <- c(-x_left_extension, maxrank + x_right_extension)
  y_limits <- c(min_y, max_y)
  early_rank <- fifelse(maxrank>=400, 350, 0)
  early_perc <- 0
  if (early_rank>0) {
    yovl = 0.5
    early_perc <- round(pldata$concordance[early_rank] * 100, 1)
    early_hjust <- fifelse(maxrank>2000, 0.5, 0.8)
    if (early_perc < 60 && early_perc > 30) {
      yovl = 0.8
    } else if (early_perc > 60) {
      yovl = 0.35
    }
  }
  ### ignore xggPairData when used standalone
  if (x_col>0 & y_col>0) { # Update global data frame
    x_index <- as.numeric(gsub("t_", "", x_col))
    y_index <- as.numeric(gsub("t_", "", y_col))
    row_index <- which((xggPairData$x == x_index & xggPairData$y == y_index) |
                        (xggPairData$x == y_index & xggPairData$y == x_index))
    xggPairData$cat_perc[row_index] <<- final_perc
    xggPairData$cat_350[row_index] <<- early_perc
  }

  perclabel=glue("**{final_perc}**%")
  borderColor  <- fifelse(standalone, "white", "#F0F0F0")
  gp <- ggplot(pldata) +
    geom_line(aes(x = rank, y = concordance), linewidth = 1, color = line_color)
    #geom_vline(xintercept = early_rank, linetype = "dashed", color = "#8080FF", linewidth=0.2) +
  if (early_rank>0) {
    gp <-  gp + annotate("segment", x = early_rank, xend = early_rank, y = -0.01, yend = 0.98,
                             linetype = "dashed", color = "#8080FF", linewidth = 0.2) +
            annotate("text", x = early_rank, y = y_limits[2], label = glue("{early_perc}%"),
                     color = "#4040D0",  hjust = early_hjust, vjust = 1.4, size = 3.2)
  }
  gp <- gp + geom_abline(intercept = 0, slope = abl_slope, linetype = 'dotted', color = "gray50", linewidth=0.4)
    # Add the elbow point to the plot
    #geom_vline(xintercept = pldata$rank[stable_start], linetype = "dashed", color = "blue") +
    #annotate("point", x = elbow$x, y = elbow$y, color = "red", size = 1.2) +
  if (!standalone) {
    gp <- gp + annotate("rect", xmin = x_limits[1], xmax = -x_left_extension/4,
            ymin = min_y, ymax = max_y, fill = "white", color = "white")
  }
  gp  <- gp + scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = x_limits) +
    scale_y_continuous(breaks = y_breaks, labels = y_breaks, limits = c(min_y, max_y)) +
    theme_minimal() +  theme(axis.ticks.y = element_blank(), panel.grid.minor = element_blank(),
                panel.grid.major = element_line(color="#F4F4F4"),
                axis.text.y = element_blank(), # axis.text.x = element_blank(),
                panel.border = element_rect(color = borderColor, fill=NA, linewidth=1)) +
    annotate("text", x = -x_left_extension, y = y_breaks_label,
             label = sprintf(showdec, y_breaks_label), hjust = -0.2, size = 3, color = "gray60") +
    annotate("rect", xmin = 0, xmax = x_limits[2], ymin = min_y, ymax = -0.04, fill = "white", color = "white") +
    geom_richtext(data = data.frame(x = maxrank, y = y_limits[2]),
                  aes(x = x, y = y, label = perclabel), fill = NA, label.color = NA,
                  hjust = 0.95, vjust = 0.91, size = 3.6, color="gray30")
  if (early_rank>0) {
    gp <- gp+ annotate("text", x = early_rank, y = min_y,
           label = sprintf("%d",early_rank),  hjust = 0.7, vjust=-0.2, size = 2.8, color = "#4040D0")
  }
  gp <- gp + annotate("text", x = maxrank, y = y_limits[1],,
             label = sprintf("n = %d",maxrank), hjust = 1.06, vjust=-0.3, size = 2.9, color = "gray45")
  if (!is.null(fdr_overlap)) {
    gp <- gp + annotate("text", x = maxrank/2, y = yovl,
            label = sprintf("%d \u2229 FDR<.05", fdr_overlap),
            hjust = 0.45, vjust = -0.12, size = 3.2, color = "gray40")
  }
  if (!no.p & !is.null(pval_overlap)) {
    gp <- gp + annotate("text", x = maxrank/2, y = yovl,
            label = sprintf("%d \u2229 p<.001", pval_overlap),
            hjust = 0.45, vjust = 1.12, size = 2.6, color = "gray40")
  }
  gp <- gp + coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE)
  if (standalone & !is.null(xlab)) gp <- gp + xlab(xlab)
  if (standalone & !is.null(ylab)) gp <- gp + ylab(ylab)

  return(gp)
}

ggScatter <- function(data, mapping, fdr = 0.05) {
  x_col <- quo_name(mapping$x)
  y_col <- quo_name(mapping$y)

  #data[is.na(data[[x_col]]), x_col] <- 0
  #data[is.na(data[[y_col]]), y_col] <- 0
  # Remove NA values instead of setting them to 0
  data <- data[complete.cases(data[c(x_col, y_col)]), ]
  # Calculate color based on adjusted p-values
  data$color <- ifelse(data[[gsub("t_", "adj.P.Val_", x_col)]] < fdr &
                         data[[gsub("t_", "adj.P.Val_", y_col)]] >= fdr, "#E04044",
                       ifelse(data[[gsub("t_", "adj.P.Val_", x_col)]] >= fdr &
                                data[[gsub("t_", "adj.P.Val_", y_col)]] < fdr, "#1150E0",
                              ifelse(data[[gsub("t_", "adj.P.Val_", x_col)]] < fdr &
                                       data[[gsub("t_", "adj.P.Val_", y_col)]] < fdr, "#D060C0",
                                     "gray50")))
  data$color <- fifelse(is.na(data$color), "gray20", data$color)
  # Calculate correlation
  #cor_value <- round(cor(data[[x_col]], data[[y_col]], use = "pairwise.complete.obs"), 2)
  cor_test <- cor.test(data[[x_col]], data[[y_col]], method = "spearman", use = "pairwise.complete.obs")
  cor_value <- round(cor_test$estimate, 2)
  p_val <- cor_test$p.value  # Update global data frame
  x_index <- as.numeric(gsub("t_", "", x_col))
  y_index <- as.numeric(gsub("t_", "", y_col))
  row_index <- which((xggPairData$x == x_index & xggPairData$y == y_index) |
                       (xggPairData$x == y_index & xggPairData$y == x_index))
  xggPairData$cor_value[row_index] <<- cor_value

  cor_label <- glue("R: **{cor_value}**, p: **{signif(p_val, 2)}**")

  p_x_col = gsub("t_", "P.Value_", x_col)
  p_y_col = gsub("t_", "P.Value_", y_col)
  adj_p_x_col = gsub("t_", "adj.P.Val_", x_col)
  adj_p_y_col = gsub("t_", "adj.P.Val_", y_col)

  # Compute replication statistics based on pi1, using each dataset as the
  # discovery set. Among significant genes in one set, what fraction in the
  # other are expected to be DE?
  p = data[[p_y_col]][data[[adj_p_x_col]] < fdr]
  if (length(p) > 0) {
    #   Compute pi1, adding a dummy p value of 1 if there's an error
    pi_discovery_x = tryCatch(
      1 - qvalue(p)$pi0,
      error = function(e) {
        1 - qvalue(c(1, p))$pi0
      }
    )
  } else {
    pi_discovery_x = 0
  }

  p = data[[p_x_col]][data[[adj_p_y_col]] < fdr]
  if (length(p) > 0) {
    #   Compute pi1, adding a dummy p value of 1 if there's an error
    pi_discovery_y = tryCatch(
      1 - qvalue(p)$pi0,
      error = function(e) {
        1 - qvalue(c(1, p))$pi0
      }
    )
  } else {
    pi_discovery_y = 0
  }

  rep_label = sprintf(
    "Rep. x->y: **%s%%**; y->x: **%s%%**",
    signif(100 * pi_discovery_x, 3),
    signif(100 * pi_discovery_y, 3)
  )

  # Calculate plot limits
  x_range <- range(data[[x_col]], na.rm = TRUE)
  y_range <- range(data[[y_col]], na.rm = TRUE)
  # Define extension factors
  ext_f_left <- 0.22
  ext_f_right <- 0.1
  ext_f_bot <- 0.15
  ext_f_top <- 0.17

  yax <- x_range[1] - 0.05*diff(x_range) # shifted y-axis position
  xay <- y_range[1] - 0.08*diff(y_range) # shifted x-axis position
  # Extend axis ranges
  x_limits <- c(x_range[1] - ext_f_left * diff(x_range),
                x_range[2] + ext_f_right * diff(x_range))
  y_limits <- c(y_range[1] - ext_f_bot * diff(y_range),
                y_range[2] + ext_f_top * diff(y_range))
  # Calculate axis breaks
  x_breaks <- pretty(x_range, n = 4)
  y_breaks <- pretty(y_range, n = 4)
  yuadj <- 0.02*diff(y_range)
  xuadj <- 0.02*diff(y_range)
  # Remove min and max breaks
  x_breaks <- x_breaks[x_breaks > x_range[1]-xuadj & x_breaks < x_range[2]]
  y_breaks <- y_breaks[y_breaks > y_range[1]-yuadj & y_breaks < y_range[2]-ext_f_top*diff(y_range)]

  ggplot(data = data, mapping = mapping) +
    geom_rect(aes(xmin = x_limits[1], xmax = x_range[2],
                  ymin = y_range[2], ymax = y_limits[2]),
              fill = "white", color = NA) +  # Top margin
    geom_rect(aes(xmin = x_limits[1], xmax = yax+xuadj,
                  ymin = y_limits[1], ymax = y_range[2]),
              fill = "white", color = NA) +  # Left margin
    geom_rect(aes(xmin = x_limits[1], xmax = x_range[2],
                  ymin = y_limits[1]-1, ymax = y_range[1]-yuadj),
              fill = "white", color = NA) +  # Bottom margin

    geom_point(alpha = 0.5, size = 0.6, aes(color = color)) +
    scale_color_identity() +
    scale_x_continuous(breaks = x_breaks,
                       expand = expansion(mult = c(ext_f_left, ext_f_right))) +
    scale_y_continuous(breaks = y_breaks,
                       expand = expansion(mult = c(ext_f_bot, ext_f_top))) +
    theme_minimal() +
    theme(axis.text.x = element_blank(),axis.text.y = element_blank(),
          axis.title = element_blank(), panel.grid.minor = element_line(color = "#F4F4F4")) +
    geom_text(data = data.frame(x = x_breaks, y = y_limits[1]),
              aes(x = x, y = y, label = sprintf("%d", x)),
              vjust = 0, hjust = 0.5, size = 2.8, color = "gray60") +
    # Add y-axis labels with white background
    geom_text(data = data.frame(x = yax, y = y_breaks),
              aes(x = x, y = y, label = sprintf("%d", y)),
              hjust = 1, vjust = 0.5, size = 2.8, color = "gray60") +
    # Add correlation text at the top, shifted right to align with x_range[1]
    geom_richtext(data = data.frame(x = x_range[1], y = y_range[2]),
                  aes(x = x, y = y, label = cor_label), hjust = 0.2, vjust = 0.2,
                  fill = NA, label.color = NA, size = 3.4, color = "gray30" ) +
    # Add replication text at the bottom
    geom_richtext(data = data.frame(x = x_range[1], y = y_range[1]),
                  aes(x = x, y = y, label = rep_label), hjust = 0.2, vjust = 0.7,
                  fill = NA, label.color = NA, size = 2, color = "gray30" ) +
    coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE)
}

ggVolcano <- function(data, source_name, x_col, no.p=FALSE) {
  logFC_col <- paste0("logFC_", source_name)
  P.Value_col <- paste0("P.Value_", source_name)
  adj.P.Val_col <- paste0("adj.P.Val_", source_name)

  # Remove rows with NAs for data.frame
  data <- data[!is.na(data[[logFC_col]]) & !is.na(data[[P.Value_col]]), ]

  significant <- data[[adj.P.Val_col]] < 0.05

  n_significant <- sum(significant)
  n_p_low <- sum(data[[P.Value_col]] < 0.001)

  # Get the actual label from column_labels
  title <- xggColumnLabels[as.numeric(source_name)]

  y_range <- range(-log10(data[[P.Value_col]]), na.rm = TRUE)
  # Calculate x-axis limits with padding
  x_range <- range(data[[logFC_col]], na.rm = TRUE)
  ext_factor_left <- 0.18
  ext_factor_right <- 0.05
  ext_factor_ybot <- 0.40 # bottom
  ext_factor_ytop <- ifelse(no.p, 0.22, 0.40) # top
  x_pad_right <- diff(x_range) * ext_factor_right  # 5% padding on the right
  x_pad_left <- diff(x_range) * ext_factor_left  # 18% padding on left side
  x_limits <- c(x_range[1] - x_pad_left, x_range[2] + x_pad_right)
  y_limits <- c(y_range[1] - ext_factor_ybot * diff(y_range), y_range[2] + ext_factor_ytop * diff(y_range))
  # Calculate y-axis limits
  max_y <- max(-log10(data[[P.Value_col]]), na.rm = TRUE)

  # Calculate y-axis breaks
  y_breaks <- pretty(c(0, max_y), n = 3)
  y_breaks <- y_breaks[y_breaks < max_y]

  x_breaks <- pretty(x_range, n = 3)  # Reduce number of breaks
  x_breaks <- x_breaks[x_breaks < x_limits[2] & x_breaks > x_limits[1]]

  # Calculate position for y-axis labels
  yax <- x_range[1] - 0.05 * diff(x_limits)
  xay <- y_range[1] - 0.05 * diff(y_limits)
  n_label = glue("FDR<.05: **{n_significant}**")
  n_plabel = glue("p<.001: {n_p_low}")
  gp <- ggplot(data, aes(x = .data[[logFC_col]], y = -log10(.data[[P.Value_col]]))) +
    geom_hline(yintercept = -log10(0.001), linetype = "dashed", linewidth = 0.3, color = "gray20") +
    geom_rect(aes(xmin = x_limits[1], xmax = yax,
                  ymin = y_limits[1], ymax = y_range[2]),
              fill = "#F4F4F4", color = NA) +  # Left margin
    geom_point(alpha = 0.5, size = 0.6, color = ifelse(significant, "#F06010", "gray40")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", linewidth=.3, color = "gray20") +
    annotate("rect", xmin = x_limits[1], xmax = x_limits[2], ymin = y_range[2]+0.04*diff(y_range),
             ymax = y_limits[2],      fill = "#F4F4F4", color = "#F4F4F4") +  # Top panel
    scale_x_continuous(expand = expansion(mult = c(ext_factor_left, ext_factor_right))) +
    scale_y_continuous(breaks = y_breaks,
                       expand = expansion(mult = c(ext_factor_ybot, ext_factor_ytop))) +
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_blank(),
          axis.title = element_text(size = 8),
          axis.text = element_text(size = 6),
          panel.background = element_rect(fill = "#F4F4F4", color = NA),
          plot.background = element_rect(fill = "#F4F4F4", color = NA),
          panel.grid.major = element_line(color = "#D8D8D8"),
          panel.grid.minor = element_line(color = "#E2E2E2"),
          axis.ticks = element_line(linewidth = 0.5) ) +
    geom_rect(aes(xmin = x_limits[1], xmax = x_limits[2],
                  ymin = y_limits[1], ymax = y_limits[1] + 0.16 * diff(y_range)),
              fill = "lightgray", color = NA) +  # Bottom panel
    annotate("text", x = mean(x_limits), y = y_limits[1], label = title, ## y = xay
             hjust = 0.5, vjust = -0.24, size = 3.2, color = "black") + # bottom title
    # Values label
    geom_richtext(data = data.frame(x = mean(x_limits), y = y_limits[2]),
                  aes(x = x, y = y, label = n_label), hjust = 0.5, vjust = 0.94,
                  fill = NA, label.color = NA, size = 3.2, color = "gray30" )
    if (!no.p) {
      gp <- gp + annotate("text", x = mean(x_limits), y = y_range[2], label = n_plabel,
                          hjust = 0.5, vjust = -0.5, size = 3, color = "gray30" )
    }
    # Add y and x-axis labels
    gp <-  gp +  geom_text(data = data.frame(x = yax, y = y_breaks),
              aes(x = x, y = y, label = sprintf("%d",y)),
              hjust = 1, vjust = 0.5, size = 2.8, color = "gray50") +
      geom_text(data = data.frame(x = x_breaks, y = xay),
                aes(x = x, y = y, label = sprintf("%.1f", x)),
                hjust = 0.5, vjust = 1, size = 2.8, color = "gray50") +
    coord_cartesian(xlim = x_limits, ylim = y_limits, expand = FALSE)
}

##xggpairs <- function(wide_data, column_labels, CAT.top=3000, no.p=FALSE) {
xggpairs <- function(deres_lst, CAT.top=3000, no.p=FALSE) {
  ## push column_labels into the global var to be used by ggVolcano()
  column_labels <- names(deres_lst)
  xggColumnLabels <<- column_labels
  wide_data <- deLst_ggwide(deres_lst)
  # Get t_ columns
  t_cols <- grep("t_", names(wide_data), value = TRUE)
  n_cols <- length(t_cols)

  # Create pairs for lower triangle only
  pairs <- expand.grid(x = 1:n_cols, y = 1:n_cols)
  pairs <- pairs[pairs$x > pairs$y, ]

  # Create and pre-populate the global data frame
  xggPairData <<- data.frame(
    x = pairs$x,
    y = pairs$y,
    x_label = column_labels[pairs$x],
    y_label = column_labels[pairs$y],
    cor_value = NA_real_,
    cat_perc = NA_real_,
    cat_350 = NA_real_,
    stringsAsFactors = FALSE
  )

  ggpairs(wide_data,
          columns = grep("t_", names(wide_data)),
          columnLabels = column_labels,
          showStrips = TRUE,
          axisLabels = "none",
          lower = list(continuous = ggScatter),
          upper = list(continuous = function(data, mapping, ...) {
            x_col <- quo_name(mapping$x)
            y_col <- quo_name(mapping$y)
            # Get P.Value columns and ensemblID
            pval_x_col <- gsub("t_", "P.Value_", x_col)
            pval_y_col <- gsub("t_", "P.Value_", y_col)
            adj_pval_x_col <- gsub("t_", "adj.P.Val_", x_col)
            adj_pval_y_col <- gsub("t_", "adj.P.Val_", y_col)

            # Calculate overlap counts
            fdr_overlap <- sum(wide_data[[adj_pval_x_col]] < 0.05 & wide_data[[adj_pval_y_col]] < 0.05, na.rm = TRUE)
            if (no.p) {
              pval_overlap <- 0
            } else {
              pval_overlap <- sum(wide_data[[pval_x_col]] < 0.001 & wide_data[[pval_y_col]] < 0.001, na.rm = TRUE)
            }
            ensembl_ids <- wide_data$ensemblID
            # Remove NAs
            valid_indices <- !is.na(wide_data[[pval_x_col]]) & !is.na(wide_data[[pval_y_col]])
            pval_x <- wide_data[[pval_x_col]][valid_indices]
            pval_y <- wide_data[[pval_y_col]][valid_indices]
            ensembl_ids <- ensembl_ids[valid_indices]
            names(pval_x) <- ensembl_ids
            names(pval_y) <- ensembl_ids
            # Sort named vectors
            vec1 <- sort(pval_x)
            vec2 <- sort(pval_y)
            ggCATplot(vec1, vec2, CAT.top, x_col, y_col, fdr_overlap, pval_overlap, no.p=no.p)
          }),
          diag = list(continuous = function(data, mapping, ...) {
            # Extract x column name from mapping
            x_col <- quo_name(mapping$x)

            src <- gsub("t_", "", x_col)
            ggVolcano(wide_data, src, x_col, no.p=no.p)
          })
  )
}
