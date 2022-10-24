#' Composition of sc3 clusters across 4 states
#' Saves and returns matrix with percent composition
#'
#' @param file path to clusters
#' @param state matrix with state annotation
#' @param state1 string with first state name
#' @param state2 string with second state name
#' @param state3 string with third state name
#' @param state4 string with fourth state name
#' @param k number of clusters
#' @param colname column name of matrix with clusters
#'
#' @return matrix
#' @export

findPercentComposition <- function(file, state,state1,state2,state3,state4,k,colname) {

  file <- read.table(file,header = T,row.names = 1,sep = '\t',check.names = F)
  sc3clustersperc = c("1","2","3","4")
  sc3clusters <- as.data.frame(file)
  state <- as.data.frame(state)
  sc3clusters <- cbind(sc3clusters,state)
  sc3clusters1 <- filter(sc3clusters,state ==state1)
  sc3clusters2 <- filter(sc3clusters,state ==state2)
  sc3clusters3 <- filter(sc3clusters,state ==state3)
  sc3clusters4 <- filter(sc3clusters,state ==state4)

  for (i in 1:k){
    total <- dim(dplyr::filter(sc3clusters,colname==i))[1]
    one <- dim(dplyr::filter(sc3clusters1,colname==i))[1]
    two <- dim(dplyr::filter(sc3clusters2,colname==i))[1]
    three <- dim(dplyr::filter(sc3clusters3,colname==i))[1]
    four <- dim(dplyr::filter(sc3clusters4,colname==i))[1]
    sc3clustersperc <- rbind(sc3clustersperc,c(one/total,two/total,three/total,four/total))
  }

  sc3clustersperc <- rbind(sc3clustersperc,c(dim(sc3clusters1)[1],
                                             dim(sc3clusters2)[1],
                                             dim(sc3clusters3)[1],
                                             dim(sc3clusters4)[1]))

  rownames(sc3clustersperc) <- c(0:k,"cellcount")
  write.table(data.frame(sc3clustersperc),'clustersperc.tsv',quote = F,sep = '\t')
  return(sc3clustersperc)
}
