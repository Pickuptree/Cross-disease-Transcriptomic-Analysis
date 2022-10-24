#' Uses two lists of results from DESeq2
#' Returns two lists, one of common upregulated DEGs between the two lists, the other of downregulated DEGs
#' Uses q-value < 0.05 and |log2FC| < 1
#'
#' @param Adata dataset 1
#' @param Bdata dataset 2
#'
#' @return None
#' @export

findCommonDEGS <- function(Adata,Bdata){
  Aupreb=Adata[Adata$padj<0.05&Adata$log2FoldChange>1,]
  Aupreg = Aupreg[order(Aupreg$padj),]

  Adownreb=Adata[Adata$padj<0.05&Adata$log2FoldChange<(1),]
  Adownreg = Adownreg[order(Adownreg$padj),]

  Bupreb=Bdata[Bdata$padj<0.05&Bdata$log2FoldChange>1,]
  Bupreg = Bupreg[order(Bupreg$padj),]

  Bdownreb=Bdata[Bdata$padj<0.05&Bdata$log2FoldChange<(1),]
  Bdownreg = Bdownreg[order(Bdownreg$padj),]

  Upreg <- merge(Aupreg,Bupreg,by="row.names")
  Downreg <- merge(Adownreg,Bdownreg,by="row.names")

  write.table(data.frame(Upreg),'Common.Upreg.tsv',quote = F,sep = '\t')
  write.table(data.frame(Downreg),'Common.Downreg.tsv',quote = F,sep = '\t')
}

