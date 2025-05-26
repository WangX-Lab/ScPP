#' Title
#'
#' @param bulk_data Log2-normalized bulk expression data with genes in row and samples in column.
#' @param features Feature data of bulk samples, column1 are sample names (colname is "Sample") and column2 are feature (colname is "Feature") labels of each sample.
#' @param ref_group A character to indicate which feature is the control group.
#' @param Log2FC_cutoff Absolute cutoff value of fold change, default is 0.585.
#'
#' @return A gene list of feature markers with two binary groups.
#' @export
#'
#' @examples
marker_Binary <- function(bulk_data, features, ref_group, Log2FC_cutoff = 0.585){
  library(dplyr)
  
  if (missing(ref_group))
    stop("'ref_group' is missing or incorrect.")

  if (missing(bulk_data) || !class(bulk_data) %in% c("matrix", "data.frame"))
		stop("'bulk_data' is missing or incorrect.")

  if (missing(features) || !class(features) %in% c("matrix", "data.frame"))
		stop("'features' is missing or incorrect.") 


  ref = features$Sample[features$Feature == ref_group]
  tes = features$Sample[features$Feature != ref_group]
  
  ref_pos = which(colnames(bulk_data) %in% ref)
  tes_pos = which(colnames(bulk_data) %in% tes)
  
  pvalues <- apply(bulk_data, 1, function(x) {
    t.test(as.numeric(x)[tes_pos], as.numeric(x)[ref_pos])$p.value
  })
  
  log2FCs <- rowMeans(bulk_data[, tes_pos]) - rowMeans(bulk_data[, ref_pos])
  
  res <- data.frame(pvalue = pvalues, log2FC = log2FCs)
  res <- res[order(res$pvalue), ]
  res$fdr <- p.adjust(res$pvalue, method = "fdr")
  
  geneList <- list(
    gene_pos = res %>% filter(pvalue < 0.05, log2FC > Log2FC_cutoff ) %>% rownames(.),
    gene_neg = res %>% filter(pvalue < 0.05, log2FC < -Log2FC_cutoff ) %>% rownames(.)
  )
  
  if(length(geneList[[1]]) > 0 & length(geneList[[2]] > 0)){return(geneList)}

  else if (length(geneList[[1]]) == 0){warning("There is no genes positively correlated with the given feature in this bulk dataset.");geneList = list(gene_pos = geneList[[2]]);return(geneList)}
  else if (length(geneList[[2]]) == 0){warning("There is no genes negatively correlated with the given feature in this bulk dataset.");geneList = list(gene_neg = geneList[[1]]);return(geneList)}
}