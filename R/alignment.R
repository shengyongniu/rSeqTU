#' Generate GTF file for training.
#' It will generate a GTF (mocked.GTF) for training in the working directory.
#' @param genomelength The length of genome
#' @import Rsubread Rsamtools QuasR e1071 seqinr grid gridBase ggplot2 reshape2
#' @import plyr caret mlbench Gviz GenomicRanges
#' @return GTF file for training
#' @export

generate_GTF_training <- function(genomelength){
## Gerneration of mocked GTF file for data traning
baseGTFforward<-data.frame(NC="NC123",source="Refseq",feature="CDS",start=c(1:genomelength), end=c(1:genomelength), annotation1=".", strand="+", annotation2=".",note="positionID")
head(baseGTFforward)
baseGTFreverse<-data.frame(NC="NC123",source="Refseq",feature="CDS",start=c(1:genomelength), end=c(1:genomelength), annotation1=".", strand="-", annotation2=".",note="positionID")
head(baseGTFreverse)
baseGTF <- rbind(baseGTFforward, baseGTFreverse)
write.table(baseGTF, "mocked.GTF", row.names = F, col.names = F, sep="\t", quote = F)
}


