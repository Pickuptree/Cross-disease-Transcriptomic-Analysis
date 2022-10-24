#' Get marker genes from monocle subset, grouped by condition
#'
#' @param cds Monocle object or subset
#' @param cell_group Numerical value of which cell groupo to return marker genes of
#' @param path File path to save to
#'
#' @return None
#' @export


getMarkerGenes <- function(cds,cell_group,path){

  marker_test_res <- top_markers(cds, group_cells_by="condition",
                                 reference_cells=1000,
                                 genes_to_test_per_group = 100,
                                 cores=8)
  marker_test_res = marker_test_res[order(marker_test_res$marker_test_q_value),]


  marker_test_res1 <- filter(marker_test_res, cell_group==cell_group)
  write.table(data.frame(marker_test_res1),path,quote = F,sep = '\t')
}
