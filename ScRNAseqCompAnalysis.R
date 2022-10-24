#' This function does comparison analysis on two separate scRNA-seq datasets.
#' The rows are genes in both the count matrix and metadata
#' Two tsv report files will be produced for the following analysis:
#' Upregulated genes, downregulated genes
#'
#' @param data1 Path to first tsv file with count matrix
#' @param metadata1 Path to first tsv file with metadata
#' @param data2 Path to second tsv file with count matrix
#' @param metadata2 Path to second tsv file with metadata
#' @param subclass String with subclass of cell to be studied
#' @param ref Dataset (one or two) that will be used as the reference
#'
#' @return None
#' @export

scRNAseqComparisonAnalysis = function(data1,metadata1,data2,metadata2,subclass,ref){
  data1 = read.table(data1,header = T,row.names = 1,sep = '\t',check.names = F)
  metadata1 = read.table(metadata1,header = T,row.names = 1,sep = '\t',check.names = F)

  data2 = read.table(data2,header = T,row.names = 1,sep = '\t',check.names = F)
  metadata2=read.table(metadata2,header = T,row.names = 1,sep = '\t',check.names = F)

  #Create count matrix
  subset1 = filter(metadata1,Subclass==subclass)
  col1 = as.vector(row.names(subset1))
  filtereddata1 = data1[,col1]

  subset2 = filter(metadata2,Subclass==subclass)
  col2 = as.vector(row.names(subset2))
  filtereddata2 = data2[,col2]

  combined = merge(filtereddata1 ,filtereddata2,by = "row.names")

  matrix <- as.matrix(combined[ , -1])
  rownames(matrix) <- combined[ , 1]


  #Create metadata
  region = c(rep("one", as.integer(dim(subset1)[1])),rep("two",as.integer(dim(subset2)[1])))

  colnames = colnames(combined)
  colnames = colnames[-1]

  meta <- data.frame(condition = region,
                         row.names = colnames)
  meta$condition <- as.factor(meta$condition)

  #Finish count data

  data=matrix[,rownames(meta)]
  row_sub = apply(data, 1, function(row) sum(row !=0 )>=0.25*length(row))
  data=data[row_sub,]
  datamatrix <- as.matrix(data)
  rownames(datamatrix) <- data[ , 1]
  dataplusone <- as.matrix(data)+1

  #Deseq analysis
  dds <- DESeqDataSetFromMatrix(countData = dataplusone, colData = meta,
                                    design = ~condition)
  cerevsolf=relevel(dds$condition,ref = ref)
  cerevsolfdeseq=DESeq(dds)

  res <- results(cerevsolfdeseq)

  upregulated=res[res$padj<0.05&res$log2FoldChange>1,]
  upregulatedpadj = upregulated[order(upregulated$padj),]
  write.table(data.frame(upregulatedpadj),'upregulated.tsv',quote = F,sep = '\t')

  downregulated=res[res$padj<0.05&res$log2FoldChange<(-1),]
  downregulatedpadj = downregulated[order(downregulated$padj),]
  write.table(data.frame(downregulatedpadj),'downregulated.tsv',quote = F,sep = '\t')
}
