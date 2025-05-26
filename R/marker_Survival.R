#' Title
#'
#' @param bulk_data Log2-normalized bulk expression data with genes in row and samples in column.
#' @param survival_data Survival data with time in column1 and status in column2. Rownames are sample name.
#'
#' @return A gene list.
#' @export
#'
#' @examples
marker_Survival <- function(bulk_data,survival_data){
  library(survival)
  library(dplyr)
  
  if (missing(bulk_data) || !class(bulk_data) %in% c("matrix", "data.frame"))
		stop("'bulk_data' is missing or incorrect.")

  if (missing(survival_data) || !class(survival_data) %in% c("matrix", "data.frame"))
		stop("'survival_data' is missing or incorrect.")
 

  SurvivalData <- data.frame(cbind(survival_data,t(bulk_data)))
  colnames(SurvivalData) = make.names(colnames(SurvivalData))
  var <- make.names(rownames(bulk_data))
  
  Model_Formula <- sapply(var, function(x) as.formula(paste("Surv(time, status) ~", x)))
  
  Model_all <- lapply(Model_Formula, function(x) coxph(x, data = SurvivalData))
  
  res <- lapply(seq_along(Model_all), function(i) {
    coef_summary <- summary(Model_all[[i]])$coefficients
    data.frame(
      variable = var[i],
      pvalue = coef_summary[,5],
      coef = coef_summary[,2]
    )
  }) %>% bind_rows()
  
  res <- res[order(res$pvalue), ]
  res$fdr <- p.adjust(res$pvalue, method = "fdr")

  geneList <- list(
    gene_pos = res %>% filter(fdr < 0.05, coef > 1) %>% pull(variable), #correalted with worse survival
    gene_neg = res %>% filter(fdr < 0.05, coef < 1) %>% pull(variable)  #correlated with better survival
  )
  if(length(geneList[[1]]) > 0 & length(geneList[[2]] > 0)){return(geneList)}

  else if (length(geneList[[1]]) == 0){warning("There is no genes negatively correlated with patients' prognosis in this bulk dataset.");geneList = list(gene_pos = geneList[[2]]);return(geneList)}
  else if (length(geneList[[2]]) == 0){warning("There is no genes positively correlated with patients' prognosis in this bulk dataset.");geneList = list(gene_neg = geneList[[1]]);return(geneList)}
}
