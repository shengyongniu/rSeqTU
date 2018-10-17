#' Generation of .NA file for data training
#' @import Rsubread Rsamtools QuasR e1071 seqinr grid gridBase ggplot2 reshape2
#' @import plyr caret mlbench Gviz GenomicRanges
#' @param bamfile The result of alignment in bam format.
#' @param ref_genome_name The file name of reference genome such as NC_000913.3
#' @param output_prefix The prefix of output file name.
#' @return .NA file for data training
#' @export

gen_NA <- function(bamfile, ref_genome_name, output_prefix){

# load bam file
#bamfile <- "Ecoli_alignment.BAM"
#bf <- BamFile(bamfile)
sortBam(file = bamfile, paste(output_prefix, ".sorted.bam", sep=""))
indexBam(file = paste(output_prefix, ".sorted.bam", sep=""))
bamfile <- paste(output_prefix, ".sorted.bam", sep="")
bf <- BamFile(bamfile)
seqlengths(bf)
param <- ScanBamParam(which=GRanges(ref_genome_name,
                                    IRanges(start=1, end=seqlengths(bf))))

p_param <- PileupParam(min_mapq = 15, min_base_quality = 10)
pileup(bf, scanBamParam=param, pileupParam = p_param)
res <- pileup(bf, scanBamParam=param, pileupParam=p_param)
plot(res$nucleotide, res$count ,pch=".", log="y", ylab="count (log scale)")

df <- dcast(res, seqnames+pos+strand ~ nucleotide, value.var="count", fun.aggregate=sum)
df_pos <- df[which(df$strand == "+"), ]
pos_hit <- rowSums(df_pos[,4:7])
df_pos <- cbind(df_pos, pos_hit)
df_pos <- df_pos[,c("pos","pos_hit")]
new_df <- data.frame(c(1:seqlengths(bf)))
colnames(new_df) <- "pos"
zz <- join(new_df, df_pos, type = "left", by = "pos")
zz[is.na(zz)] <- 0
df_pos <- zz

df_neg <- df[which(df$strand == "-"), ]
neg_hit <- rowSums(df_neg[,4:7])
df_neg <- cbind(df_neg, neg_hit)
df_neg <- df_neg[,c("pos","neg_hit")]
new_df <- data.frame(c(1:seqlengths(bf)))
colnames(new_df) <- "pos"
zz <- join(new_df, df_neg, type = "left", by = "pos")
zz[is.na(zz)] <- 0
df_neg <- zz
df_na <- join(df_pos, df_neg, type = "left", by = "pos")
df_na <- df_na[,2:3]
write.table(df_na, file=paste(output_prefix, ".NA", sep=""), col.names = FALSE, row.names = FALSE, sep="\t")
}
