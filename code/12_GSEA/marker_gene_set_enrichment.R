#' Evaluate the enrichment for a list of gene sets vs. top genes from mean ratio stats
#'
#' @param gene_list A named `list` object (could be a `data.frame`) where each
#' element of the list is a character vector of Ensembl gene IDs.
#' @param n_marker_gene A `numeric(1)` specifying the number of marker genes to select from `rank_ratio`
#' @param marker_stats Output table from `mean_ratio2()`
#'
#' @returnA table in long format with the enrichment results using
#' [stats::fisher.test()].
#' @export
#'
#' @examples#' asd_sfari <- utils::read.csv(
#'     system.file(
#'         "extdata",
#'         "SFARI-Gene_genes_01-03-2020release_02-04-2020export.csv",
#'         package = "spatialLIBD"
#'     ),
#'     as.is = TRUE
#' )
#'
#' ## Format them appropriately
#' asd_sfari_geneList <- list(
#'     Gene_SFARI_all = asd_sfari$ensembl.id,
#'     Gene_SFARI_high = asd_sfari$ensembl.id[asd_sfari$gene.score < 3],
#'     Gene_SFARI_syndromic = asd_sfari$ensembl.id[asd_sfari$syndromic == 1]
#' )
#'
#'## load marker stats
#' load(here::here("processed-data", "06_deconvolution", "run_Bisque", "marker_stats_top_25_genes.Rdata"), verbose = TRUE)
#'
#' ## Compute the gene set enrichment results
#' asd_sfari_enrichment <- marker_gene_set_enrichment(
#'     gene_list = asd_sfari_geneList,
#'     marker_stats = marker_stats
#' )
#'
#' ## Explore the results
#' asd_sfari_enrichment
#' 
marker_gene_set_enrichment <- function(gene_list, 
                                       n_marker_gene = 25,
                                       marker_stats){
  
  geneList_present <- lapply(gene_list, function(x) {
    x <- x[!is.na(x)]
    x[x %in% marker_stats$gene]
  })
  
  ## warn about low power for small geneLists
  geneList_length <- sapply(geneList_present, length)
  
  min_genes <- 25
  if (any(geneList_length < min_genes)) {
    warning(
      "Gene list with n < ",
      min_genes,
      " may have insufficent power for enrichment analysis: ",
      paste(names(geneList_length)[geneList_length < 200], collapse = " ,")
    )
  }
  
  
  enrichTab <-
    do.call(rbind, lapply(unique(marker_stats$cellType.target), function(i) {
      
      marker_stats_ct <- marker_stats |> filter(cellType.target == ct)
      
      
      tabList <- lapply(geneList_present, function(g) {
        table(
          Set = factor(marker_stats_ct$gene %in% g, c(FALSE, TRUE)),
          Marker = factor(marker_stats_ct$rank_ratio <= n_marker_genes, c(FALSE, TRUE))
        )
      })
      
      enrichList <- lapply(tabList, fisher.test, alternative = "greater")
      
      o <- data.frame(
        OR = vapply(enrichList, "[[", numeric(1), "estimate"),
        Pval = vapply(enrichList, "[[", numeric(1), "p.value"),
        test = ct,
        NumSig = vapply(tabList, function(x) {
          x[2, 2]
        }, integer(1)),
        SetSize = vapply(geneList_present, length, integer(1)),
        stringsAsFactors = FALSE
      )
      o$ID <- gsub(".odds ratio", "", rownames(o))
      rownames(o) <- NULL
      return(o)
      
    }))
}
