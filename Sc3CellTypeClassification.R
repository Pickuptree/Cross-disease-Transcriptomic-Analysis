#' Runs sc3 clustering
#' Saves an Rdata file with sc3 object
#' Saves and returns data frame with cluster annotations
#'
#' @param seuratobject Seurat object with target cells
#' @param genes List of genes
#' @param k Number of clusters
#' @param morethan5K Boolean for whether there are more than 5K cells
#'
#' @return Dataframe with sc3 cluster annotations
#' @export

classifyWithSc3 <- function(seuratobject,genes,k,morethan5K) {

  sce <- as.SingleCellExperiment(seuratobject)

  rowData(sce)$feature_symbol <- genes
  counts(sce) <- as.matrix(counts(sce))
  logcounts(sce) <- as.matrix(logcounts(sce))

  sc3obj <- SC3::sc3(sce, ks = k, biology = TRUE, n_cores = 8)

  save(sc3obj,file="sc3OUTPUT.Rdata")

  if (morethan5k){
    sc3obj <- sc3_run_svm(sc3obj, ks = k)
    save(sce,file="sc3OUTPUT.Rdata")
  }
  df <- (colData(sce)[ , grep("sc3_", colnames(colData(sce)))])
  write.table(data.frame(df),'sc3clusters.tsv',quote = F,sep = '\t')

  return(df)
}
