#' Batch normalizes different datasets using Harmony
#' Returns a seurat object
#' @param matrix1 First count matrix dataset
#' @param matrix2 Second count matrix dataset
#'
#' @return combined seurat object
#' @export

batchNormalizeSeurat <- function(matrix1, matrix2) {
  sobject1 = CreateSeuratObject(matrix1,project = "one")
  sobject2 = CreateSeuratObject(matrix2,project = "two")

  combined_sobject <- merge(sobject1, y = c(sobject2), add.cell.ids = c("one","two"), project = "project")

  combined_sobject <- NormalizeData(combined_sobject) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose = FALSE)
  combined_sobject <- RunUMAP(combined_sobject, dims = 1:30)

  combined_sobject <- RunHarmony(combined_sobject, group.by.vars = "orig.ident")
  combined_sobject <- RunUMAP(combined_sobject, reduction = "harmony", dims = 1:30)
  combined_sobject <- FindNeighbors(combined_sobject, reduction = "harmony", dims = 1:30) %>% FindClusters()

  return(combined_sobject)
}
