#' Title
#'
#' @param bulk_data Log2-normalized bulk expression data with genes in row and samples in column.
#' @param features Feature data of bulk samples, such as TMB or CNA values of each sample.
#' @param estimate_cutoff Absolute cutoff value of correlation coefficient, default is 0.2.
#' @param method Method uses for cor.test, default is "spearman", another choice is "pearson".
#'
#' @return A gene list of feature markers.
#' @export
#'
#' @examples
marker_Continuous <- function(bulk_data, features, method = "spearman", estimate_cutoff = 0.2){
  library(dplyr)

  if (missing(bulk_data) || !class(bulk_data) %in% c("matrix", "data.frame"))
		stop("'bulk_data' is missing or incorrect.")

  if (missing(features) || !class(features) %in% c("integer", "numeric"))
		stop("'features' is missing or incorrect.") 


  CorrelationTest = apply(bulk_data,1,function(x){
    pvalue = cor.test(as.numeric(x),log2(as.numeric(features)+1), method = method)$p.value
    estimate = cor.test(as.numeric(x),log2(as.numeric(features)+1), method = method)$estimate
    res = cbind(pvalue,estimate)
    return(res)
  })
  
  resc = t(CorrelationTest)
  colnames(resc) = c("pvalue","estimate")
  resc = as.data.frame(resc[order(as.numeric(resc[,1]),decreasing = FALSE),])
  resc$fdr <- p.adjust(resc$pvalue,method = "fdr")
  
  geneList <- list(
    gene_pos = resc %>% filter(fdr < 0.05, estimate > estimate_cutoff) %>% rownames(.),
    gene_neg = resc %>% filter(fdr < 0.05, estimate < -estimate_cutoff)  %>% rownames(.)
  )
  
  if(length(geneList[[1]]) > 0 & length(geneList[[2]] > 0)){return(geneList)}

  else if (length(geneList[[1]]) == 0){warning("There is no genes positively correlated with the given feature in this bulk dataset.");geneList = list(gene_pos = geneList[[2]]);return(geneList)}
  else if (length(geneList[[2]]) == 0){warning("There is no genes negatively correlated with the given feature in this bulk dataset.");geneList = list(gene_neg = geneList[[1]]);return(geneList)}
}