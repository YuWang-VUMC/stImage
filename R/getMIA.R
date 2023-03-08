#' getMIA function of Multimodal intersection analysis (MIA) proposed in
#' Moncada et al., 2020 (https://doi.org/10.1038/s41587-019-0392-8).
#' @inheritParams Seurat
#' @param singlecellobject A \code{Seurat} object of scRNAseq data.
#' @param sc_marker data frame of marker genes of cell types from scRNAseq data
#' @param st_marker data frame of marker genes of spatial clusters from
#' spatial transcriptomics data
#'
#' @return
#' @export
#'
#' @examples
getMIA <- function(singlecellobject,
                   sc_marker = NULL,
                   st_marker = NULL
                   ) {
  mia_list <- list()
  i <- 0
  assay <- singlecellobject@active.assay
  N <- nrow(singlecellobject@assays[[assay]])
  sc_num <- table(sc_marker$cluster)
  st_num <- table(st_marker$cluster)
  for (sc in 1:length(sc_num)) {
    for (st in 1:length(st_num)) {
      i <- i + 1
      m <- as.numeric(sc_num[sc])
      k <- as.numeric(st_num[st])
      q <- length(intersect(sc_marker[sc_marker$cluster == names(sc_num[sc]),
                                      "gene"],
                            st_marker[st_marker$cluster == names(st_num[st]),
                                      "gene"]))
      #Test for over-representation (enrichment)
      pvalue <- phyper(q-1, m, N-m, k, lower.tail = FALSE, log.p = FALSE)

      mia <- c(names(sc_num[sc]), names(st_num[st]), N, m, q, k, pvalue)
      mia_list[[i]] <- mia
    }
  }
  mia_df <- as.data.frame(do.call(rbind, mia_list))
  colnames(mia_df) <- c("sc_label", "st_label", "N", "m", "q", "k", "pvalue")
  mia_df$padj <- p.adjust(mia_df$pvalue, method = "BH")
  return(mia_df)
}


