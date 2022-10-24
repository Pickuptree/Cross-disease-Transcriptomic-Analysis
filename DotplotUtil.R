#' Creates dotplot
#'
#' @param data Data file
#' @param features String of features
#' @param color_scale "viridis", “magma”, “plasma”, “inferno”, “civids”, “mako”, rocket”, and "turbo"
#' @param title1 title of first plot
#' @param title2 title of second plot
#'
#' @return None
#' @export

DotPlot <- function(data,features,color_scale,title1,title2) {
  DotPlot(data, features = features, dot.min=0.001, scale.by = "size")+
    #labs(title="TAM markers")+
    theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1))+
    #theme(axis.text.y=element_blank())+
    #coord_flip()+
    theme(axis.title=element_blank())+
    scale_colour_viridis(option=color_scale)+
    geom_point(aes(size=pct.exp), shape = 21, colour="black", stroke=0.5)+
    guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"),title =title1,order=1))+
    #guides(size=guide_legend(override.aes=list(shape=21, colour="black", fill="white"),title ="",order=1))+
    guides(colour=guide_colourbar(title=title2,frame.colour = "black",frame.linewidth=1,frame.linetype = 1,order=2))
    #guides(colour=guide_colourbar(title="",frame.colour = "black",frame.linewidth=1,frame.linetype = 1,order=2))
  ggsave("dotplot.pdf",width=5.6,height=3.6)
}
