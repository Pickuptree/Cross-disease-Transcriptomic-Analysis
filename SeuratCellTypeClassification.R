#' Clusters cells based on gene expression using Seurat
#' Returns a plot with clusters in color
#' Column names of countMatrix should be same as row names of metadata
#'
#' @param countMatrix Matrix with read counts
#' @param metadata Metadata with coordinates
#'
#' @return None
#' @export

classifyWithSeurat <- function(countMatrix,metadata) {
  #Seurat
  Seurat_Object=CreateSeuratObject(countMatrix,meta.data = metadata)
  Seurat_Object=NormalizeData(Seurat_Object)
  Seurat_Object=FindVariableFeatures(Seurat_Object)
  Seurat_Object=ScaleData(Seurat_Object)
  Seurat_Object=RunPCA(Seurat_Object,npcs = 50)
  Seurat_Object=RunUMAP(Seurat_Object,dims = 1:20)
  Seurat_Object=FindNeighbors(Seurat_Object,dims = 1:20)
  Seurat_Object=FindClusters(Seurat_Object,resolution = 0.5,algorithm = 2)
  Seurat_Object.markers=FindAllMarkers(Seurat_Object,only.pos = T)

  getPalette = colorRampPalette(brewer.pal(8, "Dark2"))
  color=getPalette(length(unique(Seurat_Object@meta.data$seurat_clusters)))
  DimPlot(Seurat_Object,cols = color,label = T,pt.size = 1)

  # plotting
  Cell_spacial_information=data.frame(cell=rownames(Seurat_Object@meta.data),
                                      x=Seurat_Object@meta.data$x,
                                      y=Seurat_Object@meta.data$x,
                                      Celltype=Seurat_Object@meta.data$seurat_clusters)


  includeberg = ggplot(Cell_spacial_information,aes(x=x,y=y,colour=Celltype))+geom_point()+
    xlim(min(Cell_spacial_information$x),max(Cell_spacial_information$x))+
    ylim(min(Cell_spacial_information$y),max(Cell_spacial_information$y)) + scale_color_manual(values = color)+ theme_bw() +
    theme(
      # get rid of panel grids
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # Change plot and panel background
      plot.background=element_rect(fill = "black"),
      panel.background = element_rect(fill = 'black'),
      # Change legend
      #legend.position = c(0.6, 0.07),
      #legend.direction = "horizontal",
      legend.background = element_rect(fill = "black", color = NA),
      legend.key = element_rect(color = "gray", fill = "black"),
      legend.title = element_text(color = "white"),
      legend.text = element_text(color = "white")
    )

}
