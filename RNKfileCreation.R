#' Writes DEG object to rnk file
#'
#' @param DEG.object Differentially expressed genes
#' @param path File path to save to
#'
#' @return None
#' @export

getMonocleMarkerGenes <- function(DEG.object,path){
  DEG.object.rnk <- DEG.object %>% select(gene,avg_log2FC)
  DEG.object.rnk$gene <- str_to_upper(DEG.object.rnk$gene)
  rownames(DEG.object.rnk) <- NULL
  write.table(DEG.object.rnk,file=path,quote=F,sep="\t",row.names=F,col.names = F)
}

# all.markers <- FindAllMarkers(mouse.gbm,
#                               logfc.threshold = 0.1,
#                               test.use = "wilcox",
#                               min.pct = 0.05,
#                               base=2)
# path="markers/"
#
# for(i in 0:(length(unique(all.markers$cluster))-1)){
#   thiscluster <- all.markers %>% filter(cluster==i) %>% select(gene,p_val,avg_log2FC,pct.1,pct.2,p_val_adj)
#   write.csv(x = thiscluster,file = paste0(path,"cluster",i,".markers.csv"),row.names = FALSE)
#   rnk.me.baby(DEG.object=thiscluster,path=paste0(path,"cluster",i,".markers.rnk"))
# }
