#' Uses dataframe with row names as official gene symbol and other symbol combined with "-"
#' Removes repeats
#' Returns tsv of only official gene symbols
#'
#' @param data dataframe with row names to be split
#'
#' @return None
#' @export

splitGeneNames <- function(data){
  genes <- rownames(data)
  genes <- str_split(genes, '-',simplify=TRUE)[,2]
  genes <- unique(genes)
  write.table(data.frame(genes),"officialgenesymbols.tsv",quote = F,sep = '\t')
}

