#' Creates Monocle object
#' Data must be a matrix
#' Column names of data mustbe same as row names of meta
#' Row names of data must be same as row names of genes
#'
#' @param data Data matrix
#' @param meta metadata
#' @param genes gene list
#' @param dims dimensions to reduce to
#'
#' @return None
#' @export

createMonocle <- function(data,meta,genes,dims){
  cds <- new_cell_data_set(data,
                           cell_metadata = meta,
                           gene_metadata <- genes)
  cds <- preprocess_cds(cds, num_dim = dims)
  cds <- align_cds(cds, alignment_group = "batch")
  cds <- reduce_dimension(cds,reduction_method="UMAP")
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
}

