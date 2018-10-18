#' Generate positive and negative TU training datasets
#' @import Rsubread Rsamtools QuasR e1071 seqinr grid gridBase ggplot2 reshape2
#' @import plyr caret mlbench Gviz GenomicRanges
#' @param file_RNAseqSignals The file of .NA file
#' @param fileNamePrefix The prefix of output file name
#' @param file_gff The .gff file of reference genome
#' @param genomeFile Reference genome sequence file
#' @return positive and negative TU training datasets
#' @export

gen_cTU_data <- function(file_RNAseqSignals, fileNamePrefix, file_gff, genomeFile){
source("seqTU_functions.R")
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)

# set some global variables
MeanQuantile = 0.1
FilterMatrixQuantile = 0.3
lambda = 0
MinimumMeanExpressionValue = 1
MinimumMedianExpressionValue = 1
MaximumGapProportionValue = 0.1

#Read EBCDs of each RNAseq data set (first column:Forward strand; second column:Reverse strand)
cat("Reading the plot file ... \n")
SelectedRNAseqData<-read.table(file_RNAseqSignals, head=F)

#Read Gene Annotation
cat("Reading gene annotations from the gff file ... \n")
x <- read.table(file_gff, sep="\t", quote="")
x1<-x[x$V3=="gene",]
tag=unlist(strsplit(as.character(x1[1,]$V9),";"))[grep("locus_tag=",unlist(strsplit(as.character(x1[1,]$V9),";")))]
tag=unlist(strsplit(tag,"="))[2]
x2<-within(x1[1,], y<-data.frame(do.call('rbind',strsplit(as.character(x1[1,]$V9), ';', fixed=TRUE))))
x3<-within(x2, y<-data.frame(do.call('rbind',strsplit(as.character(x2$y$X2), '=', fixed=TRUE))))
Ngene = data.frame(tag, x3$V4, x3$V5, x3$V7)
for (i in 2:nrow(x1)){
  tag=unlist(strsplit(as.character(x1[i,]$V9),";"))[grep("locus_tag=",unlist(strsplit(as.character(x1[i,]$V9),";")))]
  tag=unlist(strsplit(tag,"="))[2]
  x2<-within(x1[i,], y<-data.frame(do.call('rbind',strsplit(as.character(x1[i,]$V9), ';', fixed=TRUE))))
  x3<-within(x2, y<-data.frame(do.call('rbind',strsplit(as.character(x2$y$X2), '=', fixed=TRUE))))
  x4<-data.frame(tag, x3$V4, x3$V5, x3$V7)
  Ngene = rbind(Ngene,x4)
}

colnames(Ngene)=c("ID", "Start", "End", "Strand") # Set Column Names
Ngene=Ngene[order(Ngene[,2]),] #Order the list by start positions of each gene
head(Ngene)


#Select all genes on forward strand
AllForwardGenes=Ngene[Ngene[,4]=="+",]
#Select all genes on reverse strand
AllReverseGenes=Ngene[Ngene[,4]=="-",]
dim(AllForwardGenes)
dim(AllReverseGenes)
cat("Done with reading gene annotations from the gff file ... \n")
cat("In this genome, there are\n 1)",nrow(AllForwardGenes),"forward genes and\n 2)",nrow(AllReverseGenes),"reverse genes.\n")


#Read bacterial Genome Sequence
cat("Reading genome sequence for the fna file ... \n")
seq=read.fasta(genomeFile,seqonly=TRUE,as.string=TRUE)
seq=unlist(strsplit(seq[[1]], split="")) #Extract string
#length(seq) #Check genome length
cat("Done with reading genome sequence from the fna file ... \n")
cat("The input genome is",length(seq),"bp long.\n")

cat("Reading intergenic regions ... \n")
AllForwardIntergenicRegions=data.frame()
AllForwardIntergenicRegions=GenerateIGList(AllForwardGenes) #All intergenic regions on forward strand
AllReverseIntergenicRegions=data.frame()
AllReverseIntergenicRegions=GenerateIGList(AllReverseGenes) #All intergenic regions on reverse strand
cat("Done with reading intergenic regions ... \n")
cat("In this genome, there are\n 1)",nrow(AllForwardIntergenicRegions),"forward intergenic regions and\n 2)",nrow(AllReverseIntergenicRegions),"reverse intergenic regions.\n")


#Get mean expression level, mean gaps and lenght for each gene on forward strand;
AllForwardGenesExpressionMean=NULL;
AllForwardGenesExpressionMedian=NULL;
AllForwardGenesExpressionSD=NULL;
AllForwardGenesGapsProportion=NULL;
AllForwardGenesGapsLongest=NULL;
AllForwardGenesLength=NULL;
for(i in 1:nrow(AllForwardGenes)){
  AllForwardGenesExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllForwardGenes, i));
  AllForwardGenesExpressionMedian[i]=median(GetRegionExp(SelectedRNAseqData,AllForwardGenes, i));
  AllForwardGenesExpressionSD[i]=sd(GetRegionExp(SelectedRNAseqData,AllForwardGenes, i));
  AllForwardGenesGapsProportion[i]=mean(GetRegionGap(SelectedRNAseqData,AllForwardGenes, i));
  AllForwardGenesGapsLongest[i]=GetTheLengthOfLongestGap(GetRegionGap(SelectedRNAseqData,AllForwardGenes, i));
  AllForwardGenesLength[i]=AllForwardGenes[i,3]-AllForwardGenes[i,2]+1;
}

#Get mean expression level, mean gaps and lenght for each gene on reverse strand;
AllReverseGenesExpressionMean=NULL;
AllReverseGenesExpressionMedian=NULL;
AllReverseGenesExpressionSD=NULL;
AllReverseGenesGapsProportion=NULL;
AllReverseGenesGapsLongest=NULL;
AllReverseGenesLength=NULL;
for(i in 1:nrow(AllReverseGenes)){
  AllReverseGenesExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllReverseGenes, i));
  AllReverseGenesExpressionMedian[i]=median(GetRegionExp(SelectedRNAseqData,AllReverseGenes, i));
  AllReverseGenesExpressionSD[i]=sd(GetRegionExp(SelectedRNAseqData,AllReverseGenes, i));
  AllReverseGenesGapsProportion[i]=mean(GetRegionGap(SelectedRNAseqData,AllReverseGenes, i));
  AllReverseGenesGapsLongest[i]=GetTheLengthOfLongestGap(GetRegionGap(SelectedRNAseqData,AllReverseGenes, i));
  AllReverseGenesLength[i]=AllReverseGenes[i,3]-AllReverseGenes[i,2]+1;
}

#Get mean expression level, mean gaps and lenght for each intergenic region on forward strand;
AllForwardIntergenicRegionsExpressionMean=NULL;
AllForwardIntergenicRegionsExpressionMedian=NULL;
AllForwardIntergenicRegionsExpressionSD=NULL;
AllForwardIntergenicRegionsGapsProportion=NULL;
AllForwardIntergenicRegionsLength=NULL;
for(i in 1:(nrow(AllForwardIntergenicRegions))){
  AllForwardIntergenicRegionsExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllForwardIntergenicRegions, i));
  AllForwardIntergenicRegionsExpressionMedian[i]=median(GetRegionExp(SelectedRNAseqData,AllForwardIntergenicRegions, i));
  AllForwardIntergenicRegionsExpressionSD[i]=sd(GetRegionExp(SelectedRNAseqData,AllForwardIntergenicRegions, i));
  AllForwardIntergenicRegionsGapsProportion[i]=mean(GetRegionGap(SelectedRNAseqData,AllForwardIntergenicRegions, i));
  AllForwardIntergenicRegionsLength[i]=AllForwardIntergenicRegions[i,3]-AllForwardIntergenicRegions[i,2]+1;
}

#Get mean expression level, mean gaps and lenght for each intergenic region on reverse strand;
AllReverseIntergenicRegionsExpressionMean=NULL;
AllReverseIntergenicRegionsExpressionMedian=NULL;
AllReverseIntergenicRegionsExpressionSD=NULL;
AllReverseIntergenicRegionsGapsProportion=NULL;
AllReverseIntergenicRegionsLength=NULL;
for(i in 1:(nrow(AllReverseIntergenicRegions))){
  AllReverseIntergenicRegionsExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllReverseIntergenicRegions, i));
  AllReverseIntergenicRegionsExpressionMedian[i]=median(GetRegionExp(SelectedRNAseqData,AllReverseIntergenicRegions, i));
  AllReverseIntergenicRegionsExpressionSD[i]=sd(GetRegionExp(SelectedRNAseqData,AllReverseIntergenicRegions, i));
  AllReverseIntergenicRegionsGapsProportion[i]=mean(GetRegionGap(SelectedRNAseqData,AllReverseIntergenicRegions, i));
  AllReverseIntergenicRegionsLength[i]=AllReverseIntergenicRegions[i,3]-AllReverseIntergenicRegions[i,2]+1;
}

#Get length proportions of all forward intergenic regions
ForwardIntergenicRegionProportion=NULL
ForwardIntergenicRegionDeviate=NULL
for(i in 1:(nrow(AllForwardGenes)-1)){
  FullLength=AllForwardGenes[i:(i+1),][2,3]-AllForwardGenes[i:(i+1),][1,2]+1
  GeneLengthDifferent = abs(AllForwardGenes[i:(i+1),][2,3]+AllForwardGenes[i:(i+1),][1,2]-AllForwardGenes[i:(i+1),][2,2]-AllForwardGenes[i:(i+1),][1,3])
  ForwardIntergenicRegionProportion[i]=(AllForwardGenes[i:(i+1),][2,2]-AllForwardGenes[i:(i+1),][1,3]-1)/FullLength
  ForwardIntergenicRegionDeviate[i]=GeneLengthDifferent/FullLength
}
cat ("11.Proportion of intergenic region's length and sum length of flanking two genes on forward strand:",summary(ForwardIntergenicRegionProportion), "\n", sep = "\t")

summary(ForwardIntergenicRegionProportion)
summary(ForwardIntergenicRegionDeviate)


#Get length proportions of all reverse intergenic regions
ReverseIntergenicRegionProportion=NULL
ReverseIntergenicRegionDeviate=NULL
for(i in 1:(nrow(AllReverseGenes)-1)){
  FullLength=AllReverseGenes[i:(i+1),][2,3]-AllReverseGenes[i:(i+1),][1,2]+1
  GeneLengthDifferent = abs(AllReverseGenes[i:(i+1),][2,3]+AllReverseGenes[i:(i+1),][1,2]-AllReverseGenes[i:(i+1),][2,2]-AllReverseGenes[i:(i+1),][1,3])
  ReverseIntergenicRegionProportion[i]=(AllReverseGenes[i:(i+1),][2,2]-AllReverseGenes[i:(i+1),][1,3]-1)/FullLength
  ReverseIntergenicRegionDeviate[i]=GeneLengthDifferent/FullLength
}
cat ("12.Proportion of intergenic region's length and sum length of flanking two genes on reverse strand:",summary(ReverseIntergenicRegionProportion), "\n", sep= "\t")
summary(ReverseIntergenicRegionProportion)
summary(ReverseIntergenicRegionDeviate)

ForwardIntergenicRegionProportionNew=ForwardIntergenicRegionProportion
ReverseIntergenicRegionProportionNew=ReverseIntergenicRegionProportion

#Generate Forward simulated TU matrix
SimulatedForwardTUMatrixExpressionMean=data.frame()
SimulatedForwardTUMatrixExpressionSD=data.frame()
SimulatedForwardTUMatrixGapProportion=data.frame()
SimulatedForwardTUMatrixGapLongest=data.frame()
SimulatedForwardTUGeneFoldChange=NULL
num=1
for(i in 1:(nrow(AllForwardGenes))){
  if ((AllForwardGenesExpressionMean[i] >= MinimumMeanExpressionValue) & (AllForwardGenesExpressionMedian[i] >= MinimumMedianExpressionValue) & (AllForwardGenesGapsProportion < MaximumGapProportionValue) & (AllForwardGenesGapsLongest[i] < 50)){
    OneForwardSimulatedTU=GenerateOneSimulatedTU(AllForwardGenes, ForwardIntergenicRegionProportionNew, ForwardIntergenicRegionDeviate, i)
    for(j in 1:4){
      SimulatedForwardTUMatrixExpressionMean[num,j]=mean(GetRegionExp(SelectedRNAseqData, OneForwardSimulatedTU,j))
      SimulatedForwardTUMatrixExpressionSD[num,j]=sd(GetRegionExp(SelectedRNAseqData, OneForwardSimulatedTU,j))
      SimulatedForwardTUMatrixGapProportion[num,j]=mean(GetRegionGap(SelectedRNAseqData, OneForwardSimulatedTU,j))
      SimulatedForwardTUMatrixGapLongest[num,j]=GetTheLengthOfLongestGap(GetRegionGap(SelectedRNAseqData, OneForwardSimulatedTU,j))

    }
    SimulatedForwardTUGeneFoldChange[num] = FoldChange(mean(GetRegionExp(SelectedRNAseqData, OneForwardSimulatedTU,1)),mean(GetRegionExp(SelectedRNAseqData, OneForwardSimulatedTU,3))) # 1 is the 5' gene; 3 is the 3' gene.
    num = num +1
  }
}


#Generate Reverse simulated TU matrix
SimulatedReverseTUMatrixExpressionMean=data.frame()
SimulatedReverseTUMatrixExpressionSD=data.frame()
SimulatedReverseTUMatrixGapProportion=data.frame()
SimulatedReverseTUMatrixGapLongest=data.frame()
SimulatedReverseTUGeneFoldChange=NULL
num=1
for(i in 1:(nrow(AllReverseGenes))){
  if ((AllReverseGenesExpressionMean[i] >= MinimumMeanExpressionValue) & (AllReverseGenesExpressionMedian[i] >= MinimumMedianExpressionValue) & (AllReverseGenesGapsProportion < MaximumGapProportionValue) & (AllReverseGenesGapsLongest[i] < 50)){
    OneReverseSimulatedTU=GenerateOneSimulatedTU(AllReverseGenes, ReverseIntergenicRegionProportionNew, ReverseIntergenicRegionDeviate, i)
    for(j in 1:4){
      SimulatedReverseTUMatrixExpressionMean[num,j]=mean(GetRegionExp(SelectedRNAseqData, OneReverseSimulatedTU,j))
      SimulatedReverseTUMatrixExpressionSD[num,j]=sd(GetRegionExp(SelectedRNAseqData, OneReverseSimulatedTU,j))
      SimulatedReverseTUMatrixGapProportion[num,j]=mean(GetRegionGap(SelectedRNAseqData, OneReverseSimulatedTU,j))
      SimulatedReverseTUMatrixGapLongest[num,j]=GetTheLengthOfLongestGap(GetRegionGap(SelectedRNAseqData, OneReverseSimulatedTU,j))
    }
    SimulatedReverseTUGeneFoldChange[num] = FoldChange(mean(GetRegionExp(SelectedRNAseqData, OneReverseSimulatedTU,1)),mean(GetRegionExp(SelectedRNAseqData, OneReverseSimulatedTU,3))) # 1 is the 5' gene; 3 is the 3' gene.
    num = num +1
  }
}


OriginalSimulatedTwoStrandsTUMatrix=cbind(
  rbind(SimulatedForwardTUMatrixExpressionMean,SimulatedReverseTUMatrixExpressionMean),
  rbind(SimulatedForwardTUMatrixExpressionSD,SimulatedReverseTUMatrixExpressionSD),
  rbind(SimulatedForwardTUMatrixGapProportion,SimulatedReverseTUMatrixGapProportion),
  rbind(SimulatedForwardTUMatrixGapLongest,SimulatedReverseTUMatrixGapLongest),
  c(SimulatedForwardTUGeneFoldChange,SimulatedReverseTUGeneFoldChange)
)
colnames(OriginalSimulatedTwoStrandsTUMatrix)=c(
  "ExpressionMean1","ExpressionMean2","ExpressionMean3","ExpressionMean4",
  "ExpressionSD1","ExpressionSD2","ExpressionSD3","ExpressionSD4",
  "GapProportion1","GapProportion2","GapProportion3","GapProportion4",
  "GapLongest1","GapLongest2","GapLongest3","GapLongest4",
  "GeneFoldChange"
)


#colnames(OriginalSimulatedTwoStrandsTUMatrix)=c("GapPercentage", "FoldChange", "MeanExpressionIntergenicRegion", "MeanExpressionFiveEnd", "")


#geting real cases from the genome
TargetForwardTUMatrixExpressionMean=data.frame()
TargetForwardTUMatrixExpressionSD=data.frame()
TargetForwardTUMatrixGapProportion=data.frame()
TargetForwardTUMatrixGapLongest=data.frame()
TargetForwardTUGeneFoldChange=NULL
num=1
for(i in 1:(nrow(AllForwardGenes)-1)){
  OneForwardTargetTU=GenerateOneTargetTU(AllForwardGenes,i)
  for(j in 1:4){
    TargetForwardTUMatrixExpressionMean[num,j]=mean(GetRegionExp(SelectedRNAseqData, OneForwardTargetTU,j))
    TargetForwardTUMatrixExpressionSD[num,j]=sd(GetRegionExp(SelectedRNAseqData, OneForwardTargetTU,j))
    TargetForwardTUMatrixGapProportion[num,j]=mean(GetRegionGap(SelectedRNAseqData, OneForwardTargetTU,j))
    TargetForwardTUMatrixGapLongest[num,j]=GetTheLengthOfLongestGap(GetRegionGap(SelectedRNAseqData, OneForwardTargetTU,j))
  }
  TargetForwardTUGeneFoldChange[num] = FoldChange(mean(GetRegionExp(SelectedRNAseqData, OneForwardTargetTU,1)),mean(GetRegionExp(SelectedRNAseqData, OneForwardTargetTU,3))) # 1 is the 5' gene; 3 is the 3' gene.
  num = num +1
}
#which(is.na(TargetForwardTUMatrix[,1])) #check NaN elements
#which(is.na(TargetForwardTUMatrix[,2])) #check NaN elements
TargetReverseTUMatrixExpressionMean=data.frame()
TargetReverseTUMatrixExpressionSD=data.frame()
TargetReverseTUMatrixGapProportion=data.frame()
TargetReverseTUMatrixGapLongest=data.frame()
TargetReverseTUGeneFoldChange=NULL
num=1
for(i in 1:(nrow(AllReverseGenes)-1)){
  OneReverseTargetTU=GenerateOneTargetTU(AllReverseGenes,i)
  for(j in 1:4){
    TargetReverseTUMatrixExpressionMean[num,j]=mean(GetRegionExp(SelectedRNAseqData, OneReverseTargetTU,j))
    TargetReverseTUMatrixExpressionSD[num,j]=sd(GetRegionExp(SelectedRNAseqData, OneReverseTargetTU,j))
    TargetReverseTUMatrixGapProportion[num,j]=mean(GetRegionGap(SelectedRNAseqData, OneReverseTargetTU,j))
    TargetReverseTUMatrixGapLongest[num,j]=GetTheLengthOfLongestGap(GetRegionGap(SelectedRNAseqData, OneReverseTargetTU,j))
  }
  TargetReverseTUGeneFoldChange[num] = FoldChange(mean(GetRegionExp(SelectedRNAseqData, OneReverseTargetTU,1)),mean(GetRegionExp(SelectedRNAseqData, OneReverseTargetTU,3))) # 1 is the 5' gene; 3 is the 3' gene.
  num = num +1
}



#Output TargetSVM format
TargetPositiveTUMatrix=cbind(
  TargetForwardTUMatrixExpressionMean,
  TargetForwardTUMatrixExpressionSD,
  TargetForwardTUMatrixGapProportion,
  TargetForwardTUMatrixGapLongest,
  TargetForwardTUGeneFoldChange
)
TargetNegativeTUMatrix=cbind(
  TargetReverseTUMatrixExpressionMean,
  TargetReverseTUMatrixExpressionSD,
  TargetReverseTUMatrixGapProportion,
  TargetReverseTUMatrixGapLongest,
  TargetReverseTUGeneFoldChange
)

TargetTwoStrandsTUMatrix=cbind(
  rbind(TargetForwardTUMatrixExpressionMean,TargetReverseTUMatrixExpressionMean),
  rbind(TargetForwardTUMatrixExpressionSD,TargetReverseTUMatrixExpressionSD),
  rbind(TargetForwardTUMatrixGapProportion,TargetReverseTUMatrixGapProportion),
  rbind(TargetForwardTUMatrixGapLongest,TargetReverseTUMatrixGapLongest),
  c(TargetForwardTUGeneFoldChange,TargetReverseTUGeneFoldChange)
)
colnames(TargetPositiveTUMatrix)=colnames(TargetNegativeTUMatrix)=colnames(TargetTwoStrandsTUMatrix)=c(
  "ExpressionMean1","ExpressionMean2","ExpressionMean3","ExpressionMean4",
  "ExpressionSD1","ExpressionSD2","ExpressionSD3","ExpressionSD4",
  "GapProportion1","GapProportion2","GapProportion3","GapProportion4",
  "GapLongest1","GapLongest2","GapLongest3","GapLongest4",
  "GeneFoldChange"
)


#only keep fold chang < 10 in pseudo TU
#SimulatedPositiveTUMatrix=SimulatedPositiveTUMatrix[which(SimulatedPositiveTUMatrix[,2] < 10),]
SimulatedPositiveTUMatrix = OriginalSimulatedTwoStrandsTUMatrix
#SimulatedPositiveTUMatrix=rbind(SimulatedPositiveTUMatrix,TargetTwoStrandsTUMatrix[Temp454PositiveIndexInTraining,])
dim(SimulatedPositiveTUMatrix)
head(SimulatedPositiveTUMatrix)

#SimulatedNegativeTUMatrix = TargetTwoStrandsTUMatrix[Temp454NegativeIndex,]
SimulatedNegativeTUMatrix = TargetTwoStrandsTUMatrix
SimulatedNegativeTUMatrix = SimulatedNegativeTUMatrix[which(SimulatedNegativeTUMatrix[,"GeneFoldChange"] >= 7),]   # fold changes > 7
#SimulatedNegativeTUMatrix = SimulatedNegativeTUMatrix[which(SimulatedNegativeTUMatrix[,"ExpressionSD4" ]> 100),]
#SimulatedNegativeTUMatrix = SimulatedNegativeTUMatrix[which(SimulatedNegativeTUMatrix[,"GapLongest2"] > 100),]
dim(SimulatedNegativeTUMatrix)
head(SimulatedNegativeTUMatrix)



#Output SVM format for Training dataset
iter = c(1:ncol(TargetTwoStrandsTUMatrix))

sink(paste("SVM.TwoStrands.Training.",fileNamePrefix,sep=""))
#SimulatedTwoStrandsTUMatrix=rbind(SimulatedForwardTUMatrix, SimulatedReverseTUMatrix)
for(i in 1:nrow(SimulatedPositiveTUMatrix)){
  cat("+1")
  #   for (j in 1:ncol(SimulatedPositiveTUMatrix)){
  for (j in iter){
    cat("\t",j,":",SimulatedPositiveTUMatrix[i,j],sep="")
  }
  cat("\n")
}
for(i in 1:nrow(SimulatedNegativeTUMatrix)){
  cat("-1")
  #   for (j in 1:ncol(SimulatedNegativeTUMatrix)){
  for (j in iter){
    cat("\t",j,":",SimulatedNegativeTUMatrix[i,j],sep="")
  }
  cat("\n")
}
sink()


#Output SVM format for Target dataset
iter = c(1:ncol(TargetTwoStrandsTUMatrix))
sink(paste("SVM.TwoStrands.Target.",fileNamePrefix,sep=""))
for(i in 1:nrow(TargetPositiveTUMatrix)){
  cat("+1")
  for (j in iter){
    cat("\t",j,":",TargetPositiveTUMatrix[i,j],sep="")
  }
  cat("\n")
}
for(i in 1:nrow(TargetNegativeTUMatrix)){
  cat("-1")
  for (j in iter){
    cat("\t",j,":",TargetNegativeTUMatrix[i,j],sep="")
  }
  cat("\n")
}
sink()

write.table(TargetPositiveTUMatrix, file="TargetPositiveTUMatrix.txt", sep="\t", row.names = FALSE)
write.table(TargetNegativeTUMatrix, file="TargetNegativeTUMatrix.txt", sep="\t", row.names = FALSE)
write.table(SimulatedNegativeTUMatrix, file="SimulatedNegativeTUMatrix.txt", sep="\t", row.names = FALSE)
write.table(SimulatedPositiveTUMatrix, file="SimulatedPositiveTUMatrix.txt", sep="\t", row.names = FALSE)
}
