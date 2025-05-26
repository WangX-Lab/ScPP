#' Title
#'
#' @param sc_dataset A seurat object of single cell RNA sequencing data.
#' @param geneList A gene list correlated with interested features.
#' @param probs Cutoff value of ScPP, default is 0.2.
#'
#' @return Meta data with information of ScPP selected cells.
#' @export
#'
#' @examples
ScPP = function(sc_dataset, geneList, probs = 0.2){
  if(length(geneList) != 2){
    stop("This gene list do not have enough information correlated with interested feature.")
  }

  if (missing(sc_dataset) || class(sc_dataset) != "Seurat")
    stop("'sc_dataset' is missing or not a seurat object.")

  library(AUCell)
  seurat_version <- packageVersion("Seurat")
  if(!is.na(seurat_version) && seurat_version >= "5.0.0"){
    cellrankings = AUCell_buildRankings(sc_dataset@assays$RNA$data,plotStats = FALSE)
  }else{
    cellrankings = AUCell_buildRankings(sc_dataset@assays$RNA@data,plotStats = FALSE)
  }
  cellAUC = AUCell_calcAUC(geneList,cellrankings)

  metadata = as.data.frame(sc_dataset@meta.data)
  metadata$AUCup <- as.numeric(getAUC(cellAUC)["gene_pos", ])
  metadata$AUCdown <- as.numeric(getAUC(cellAUC)["gene_neg", ])
  
  downcells1 = rownames(metadata)[which(metadata$AUCup <= quantile(metadata$AUCup,probs = probs))]
  upcells1 = rownames(metadata)[which(metadata$AUCup >= quantile(metadata$AUCup,probs = (1-probs)))]
  downcells2 = rownames(metadata)[which(metadata$AUCdown >= quantile(metadata$AUCdown,probs = (1-probs)))]
  upcells2 = rownames(metadata)[which(metadata$AUCdown <= quantile(metadata$AUCdown,probs = probs))]
  
  ScPP_neg = intersect(downcells1,downcells2)
  ScPP_pos = intersect(upcells1,upcells2)
  
  metadata$ScPP <- ifelse(rownames(metadata) %in% ScPP_pos, "Phenotype+", "Background")
  metadata$ScPP <- ifelse(rownames(metadata) %in% ScPP_neg, "Phenotype-", metadata$ScPP)
  
  sc$ScPP = metadata$ScPP
  Idents(sc) = "ScPP"
  
  markers <- FindMarkers(sc, ident.1 = "Phenotype+", ident.2 = "Phenotype-")
  
  genes_pos <- rownames(markers[which(markers$avg_log2FC > 1 & markers$p_val_adj < 0.05),])
  if (length(genes_pos) == 0) {
    message("There are no genes significantly upregulated in Phenotype+ compared to Phenotype-.")
  }
  
  genes_neg <- rownames(markers[which(markers$avg_log2FC < -1 & markers$p_val_adj < 0.05),])
  if (length(genes_neg) == 0) {
    message("There are no genes significantly upregulated in Phenotype- compared to Phenotype+.")}
  
  res <- list(metadata = metadata,
              Genes_pos = genes_pos,
              Genes_neg = genes_neg)

  return(res)
}

