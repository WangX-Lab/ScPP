#' Title
#'
#' @param counts Count data of scRNA sequencing.
#' @param project Project name.
#' @param normalization.method Normalization method.
#' @param scale.factor 
#' @param selection.method 
#' @param resolution 
#' @param dims_Neighbors 
#' @param dims_TSNE 
#' @param dims_UMAP 
#'
#' @return A seurat object.
#' @export
#'
#' @examples
sc_Preprocess <- function(counts, project = "sc_preprocess",
                          normalization.method = "LogNormalize", scale.factor = 10000,
                          selection.method = "vst", resolution = 0.1,
                          dims_Neighbors = 1:20, dims_TSNE = 1:20, dims_UMAP = 1:20){
  library(Seurat)
  sc_dat <- CreateSeuratObject(counts = counts, project = project)
  sc_dat <- NormalizeData(object = sc_dat, normalization.method = normalization.method, scale.factor = scale.factor)
  sc_dat <- FindVariableFeatures(object = sc_dat, selection.method = selection.method)
  sc_dat <- ScaleData(object = sc_dat, features = rownames(sc_dat))
  sc_dat <- RunPCA(object = sc_dat, features = VariableFeatures(sc_dat))
  sc_dat <- FindNeighbors(object = sc_dat, dims = dims_Neighbors)
  sc_dat <- FindClusters(object = sc_dat, resolution = resolution)
  sc_dat <- RunTSNE(object = sc_dat, dims = dims_TSNE)
  sc_dat <- RunUMAP(object = sc_dat, dims = dims_UMAP)
  
  return(sc_dat)
}