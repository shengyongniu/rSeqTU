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
  #Function: Generate intergenic regions according to a given gene list
  GenerateIGList<-function(InputGeneList){
    IGList=data.frame()
    for(i in 1:(nrow(InputGeneList)-1)){
      Name=paste(InputGeneList[i:(i+1),][1,1], InputGeneList[i:(i+1),][2,1])
      Start=InputGeneList[i:(i+1),][1,3]+1
      End=InputGeneList[i:(i+1),][2,2]-1
      Strand=InputGeneList[i:(i+1),][1,4]
      IGList[i,1]=Name
      IGList[i,2]=Start
      IGList[i,3]=End
      IGList[i,4]=Strand
    }
    colnames(IGList)=c("ID", "Start", "End", "Strand") # Set Column Names
    return(IGList)
  } #END Function: Generate intergenic regions according to a given gene list


  #Function: Get all expression values of a region
  GetRegionExp<-function(SelectedRNAseqData, InputGeneList, i){
    if(InputGeneList[i,4] == "+"){
      SubExp=SelectedRNAseqData[InputGeneList[i,2]:InputGeneList[i,3],1]
    }else{
      SubExp=SelectedRNAseqData[InputGeneList[i,2]:InputGeneList[i,3],2]
    }
    return(SubExp)
  }#END Function: Get all expression values of a region

  #Function: Get all expression values of a region in a TU final table
  GetRegionExp_FinalTable<-function(SelectedRNAseqData, start, end, strand){
    if(strand == "+"){
      SubExp=SelectedRNAseqData[start:end, 1]
    }else{
      SubExp=SelectedRNAseqData[start:end, 2]
    }
    return(SubExp)
  }#END Function: Get all expression values of a region


  #Function: Get all gap values of a region
  GetRegionGap<-function(SelectedRNAseqData, InputGeneList, i){
    if(InputGeneList[i,4] == "+"){
      SubGap=(SelectedRNAseqData[InputGeneList[i,2]:InputGeneList[i,3],1]==0)*1
    }else{
      SubGap=(SelectedRNAseqData[InputGeneList[i,2]:InputGeneList[i,3],2]==0)*1
    }
    return(SubGap)
  }#END Function: Get all gap values of a region

  #Function: Get the length of the longest gap (continuous "1")
  GetTheLengthOfLongestGap <- function(x){
    maxGapLength=0;
    continuousGap=0;
    for(i in 1:length(x)){
      if(x[i]==1){
        continuousGap = continuousGap + 1;
      }else{
        if(maxGapLength < continuousGap){
          maxGapLength = continuousGap
        }
        continuousGap = 0;
      }
    }
    return(maxGapLength);
  }#END: Function: Get the length of the longest gap (continuous "1")



  getIntergenicInfo<-function(SelectedRNAseqData){
    AllForwardIntergenicRegionsExpressionMean=NULL;
    AllForwardIntergenicRegionsExpressionMedian=NULL;
    AllForwardIntergenicRegionsGapsProportion=NULL;
    AllForwardIntergenicRegionsLength=NULL;
    for(i in 1:(nrow(AllForwardIntergenicRegions))){
      AllForwardIntergenicRegionsExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllForwardIntergenicRegions, i));
      AllForwardIntergenicRegionsExpressionMedian[i]=median(GetRegionExp(SelectedRNAseqData,AllForwardIntergenicRegions, i));
      AllForwardIntergenicRegionsGapsProportion[i]=mean(GetRegionGap(SelectedRNAseqData,AllForwardIntergenicRegions, i));
      AllForwardIntergenicRegionsLength[i]=AllForwardIntergenicRegions[i,3]-AllForwardIntergenicRegions[i,2]+1;
    }
    #Get mean expression level, mean gaps and lenght for each intergenic region on reverse strand
    AllReverseIntergenicRegionsExpressionMean=NULL;
    AllReverseIntergenicRegionsExpressionMedian=NULL;
    AllReverseIntergenicRegionsGapsProportion=NULL;
    AllReverseIntergenicRegionsLength=NULL;
    for(i in 1:(nrow(AllReverseIntergenicRegions))){
      AllReverseIntergenicRegionsExpressionMean[i]=mean(GetRegionExp(SelectedRNAseqData,AllReverseIntergenicRegions, i));
      AllReverseIntergenicRegionsExpressionMedian[i]=median(GetRegionExp(SelectedRNAseqData,AllReverseIntergenicRegions, i));
      AllReverseIntergenicRegionsGapsProportion[i]=mean(GetRegionGap(SelectedRNAseqData,AllReverseIntergenicRegions, i));
      AllReverseIntergenicRegionsLength[i]=AllReverseIntergenicRegions[i,3]-AllReverseIntergenicRegions[i,2]+1;
    }
    MyList<- list("a"=AllForwardIntergenicRegionsExpressionMean, "b"=AllReverseIntergenicRegionsExpressionMean, "c"=AllForwardIntergenicRegionsExpressionMedian, "d"=AllReverseIntergenicRegionsExpressionMedian);
    return(MyList);
  }#end function getIntergenicInfo<-function(SelectedRNAseqData){


  #Function: Generate one simulated TU
  GenerateOneSimulatedTU<-function(InputGeneList, IntergenicRegionProportion, IntergenicRegionDeviate, i){
    SimIGLength=as.integer((InputGeneList[i,3]-InputGeneList[i,2]+1)*sample(IntergenicRegionProportion,1)) #Need a real probability
    if(SimIGLength %% 2 == 1){SimIGLength=SimIGLength+1}
    MiddlePoint=as.integer(median(InputGeneList[i,2]:InputGeneList[i,3]))
    #get the maximal AT-rich region in the perturbation
    GCvalue = 1
    for (j in 1:as.integer((InputGeneList[i,3]-InputGeneList[i,2]+1)*median(IntergenicRegionDeviate)/2)){
      SubSeq=seq[(MiddlePoint-(0.5*SimIGLength)+1-j):(SimIGLength+MiddlePoint-(0.5*SimIGLength)-j)]
      if( GC(SubSeq)<=GCvalue & length(SubSeq)>1){
        GCvalue = GC(SubSeq)
        MiddlePointNew=MiddlePoint-j
      }
    }
    for (j in 1:as.integer((InputGeneList[i,3]-InputGeneList[i,2]+1)*median(IntergenicRegionDeviate)/2)){
      SubSeq=seq[(MiddlePoint-(0.5*SimIGLength)+1+j):(SimIGLength+MiddlePoint-(0.5*SimIGLength)+j)]
      if( GC(SubSeq)<=GCvalue & length(SubSeq)>1){
        GCvalue = GC(SubSeq)
        MiddlePointNew=MiddlePoint+j
      }
    }
    MiddlePoint = MiddlePointNew
    LeftSubGeneStart=InputGeneList[i,2]
    LeftSubGeneEnd=MiddlePoint-(0.5*SimIGLength)
    IGSubGeneStart=MiddlePoint-(0.5*SimIGLength)+1
    IGSubGeneEnd=SimIGLength+MiddlePoint-(0.5*SimIGLength)
    RightSubGeneStart=SimIGLength+MiddlePoint-(0.5*SimIGLength)+1
    RightSubGeneEnd=InputGeneList[i,3]
    SubGeneList=NULL
    SubGeneList=c("LeftSubGene", LeftSubGeneStart, LeftSubGeneEnd, as.character(InputGeneList[i,4]));
    SubGeneList=rbind(SubGeneList, c("IGSubGene", IGSubGeneStart, IGSubGeneEnd, as.character(InputGeneList[i,4])));
    SubGeneList=rbind(SubGeneList, c("RightSubGene", RightSubGeneStart, RightSubGeneEnd, as.character(InputGeneList[i,4])));
    SubGeneList=rbind(SubGeneList, c("FullSubGene", InputGeneList[i,2], InputGeneList[i,3], as.character(InputGeneList[i,4])));
    row.names(SubGeneList)=NULL
    SubGeneList=as.data.frame(SubGeneList)
    SubGeneList[,2]=as.numeric(as.character(SubGeneList[,2]))
    SubGeneList[,3]=as.numeric(as.character(SubGeneList[,3]))
    return(SubGeneList)
  }#END Function: Generate one simulated TU

  #Function: Get gap percentage of a Intergenic region over a full TU
  GetGapPercentageofIntergenicRegionOverATU<-function(NAExp, InputGeneList, i){
    if(InputGeneList[i,4] == "+"){
      if(sum((NAExp[InputGeneList[4,2]:InputGeneList[4,3],1]==0)*1)==0){
        SubGaps=0
      }else{
        SubGaps=sum((NAExp[InputGeneList[i,2]:InputGeneList[i,3],1]==0)*1)/sum((NAExp[InputGeneList[4,2]:InputGeneList[4,3],1]==0)*1)
      }
    }else{
      if(sum((NAExp[InputGeneList[4,2]:InputGeneList[4,3],2]==0)*1)==0){
        SubGaps=0
      }else{
        SubGaps=sum((NAExp[InputGeneList[i,2]:InputGeneList[i,3],2]==0)*1)/sum((NAExp[InputGeneList[4,2]:InputGeneList[4,3],2]==0)*1)
      }
    }
    return(SubGaps)
  }#END Function: Get gap percentage of a Intergenic region over a full TU



  #Function: Calculate foldchanges
  FoldChange<-function(num1, num2){
    MAX=NULL
    MAX=max(num1,num2)
    MIN=NULL
    if(min(num1,num2) <= 0){
      MIN=0.00001
    }else{
      MIN=min(num1, num2)
    }
    Change=MAX/MIN
    return(Change)
  }#END Function: Calculate foldchanges

  #Function: Get New Start and End positions for removing the position bias
  NewStartEnd<-function(Start, End){
    Range=c(Start:End)
    NewStartInFactor=round(length(Range)/100*5,0)
    NewEndInFactor=length(Range)-NewStartInFactor
    NewStart=Range[NewStartInFactor]
    NewEnd=Range[NewEndInFactor]
    Positions=c(NewStart, NewEnd)
    return(Positions)
  }#END Function: Get New Start and End positions for removing the position bias



  #Function: Get Expression of a region with its position bias removed
  GetRegionExpRemovePositionBias<-function(NAExp, InputGeneList, i){
    if(InputGeneList[i,4] == "+"){
      NewRange=NewStartEnd(InputGeneList[i,2], InputGeneList[i,3])
      NewRange=c(NewRange[1]:NewRange[2])
      SubExp=NAExp[NewRange,1]
    }else{
      SubExp=NAExp[NewRange,2]
    }
    return(SubExp)
  }#END Function: Get Expression of a region with its position bias removed

  #Generate Target TU SVM format
  GenerateOneTargetTU<-function(InputGeneList, i){
    OneTUList=data.frame()
    OneTUList=c("5'Gene", InputGeneList[i:(i+1),][1,2], InputGeneList[i:(i+1),][1,3], as.character(InputGeneList[i:(i+1),][1,4]));
    OneTUList=rbind(OneTUList, c("IntergenicRegion", (InputGeneList[i:(i+1),][1,3]+1), (InputGeneList[i:(i+1),][2,2]-1),as.character( InputGeneList[i:(i+1),][1,4])));
    OneTUList=rbind(OneTUList, c("3'Gene", as.numeric(as.character(InputGeneList[i:(i+1),][2,2])), InputGeneList[i:(i+1),][2,3], as.character(InputGeneList[i:(i+1),][2,4])));
    OneTUList=rbind(OneTUList, c("FullTU", InputGeneList[i:(i+1),][1,2], InputGeneList[i:(i+1),][2,3], as.character(InputGeneList[i:(i+1),][1,4])));
    row.names(OneTUList)=NULL
    OneTUList=as.data.frame(OneTUList)
    OneTUList[,2]=as.numeric(as.character(OneTUList[,2]))
    OneTUList[,3]=as.numeric(as.character(OneTUList[,3]))
    colnames(OneTUList)=c("ID", "Start", "End", "Strand") # Set Column Names
    return(OneTUList)
  }

  GenerateFinalTUTable<-function(TargetTUsPrediction, SelectedGenes, PTT){
    #TargetTUsPrediction=ForwardPredictResult
    #SelectedGenes=AllForwardGenes
    #PTT=ptt
    colnames(TargetTUsPrediction)[1] = "TUresult"
    TargetTUName1Name2StartEnd=data.frame()
    j=1
    for(i in 1:(nrow(SelectedGenes)-1)){
      TargetTUName1Name2StartEnd=rbind(TargetTUName1Name2StartEnd, cbind(SelectedGenes[i,1], SelectedGenes[(i+1),1],GenerateOneTargetTU(SelectedGenes,i)[4,2:3]))
      j=j+1
    }
    TargetTUsPrediction = (cbind(TargetTUsPrediction, TargetTUName1Name2StartEnd))
    #   head(TargetTUsPrediction)

    colnames(TargetTUsPrediction)[2:3]=c("GeneName1", "GeneName2")
    GeneProduct1=NULL
    GeneProduct1=as.character(PTT$Product[match(TargetTUsPrediction$GeneName1, PTT$Synonym)])
    GeneProduct1[which(is.na(GeneProduct1))]=as.character(TargetTUsPrediction$GeneName1[which(is.na(GeneProduct1))])
    #   GeneProduct1[1669:1675]
    GeneProduct2=NULL
    GeneProduct2=as.character(PTT$Product[match(TargetTUsPrediction$GeneName2, PTT$Synonym)])
    GeneProduct2[which(is.na(GeneProduct2))]=as.character(TargetTUsPrediction$GeneName2[which(is.na(GeneProduct2))])
    TargetTUsPrediction=data.frame(TargetTUsPrediction, GeneProduct1=GeneProduct1, GeneProduct2=GeneProduct2)

    #   head(TargetTUsPrediction)

    FinalTUAll=data.frame()
    FinalTUOne=data.frame()
    FinalTUOneNameStart=as.character(TargetTUsPrediction[1,"GeneName1"])
    FinalTUOneProductStart=as.character(TargetTUsPrediction[1,"GeneProduct1"])
    FinalTUOneNameEnd=as.character("")
    FinalTUOneProductEnd=as.character("")
    FinalTUOnePosStart=TargetTUsPrediction[1,4]
    FinalTUOnePosEnd=TargetTUsPrediction[1,5]
    FinalTUOne=data.frame(FinalTUOneNameStart,FinalTUOneNameEnd,FinalTUOnePosStart,FinalTUOnePosEnd,FinalTUOneProductStart,FinalTUOneProductEnd)
    flag=TargetTUsPrediction[1,1]
    for(i in c(1:nrow(TargetTUsPrediction))){
      if(TargetTUsPrediction[i,1]==1){  #TU result = 1
        FinalTUOne[,2]=paste(FinalTUOne[,2],TargetTUsPrediction[i,"GeneName2"], sep=";")
        FinalTUOne[,6]=paste(FinalTUOne[,6],TargetTUsPrediction[i,"GeneProduct2"], sep=";")
        FinalTUOne[,4]=TargetTUsPrediction[i,5]
        flag=1
      }else{                                   #TU result = -1
        FinalTUAll=rbind(FinalTUAll,FinalTUOne)
        FinalTUOne=NULL
        FinalTUOneNameStart=as.character(TargetTUsPrediction[i,"GeneName2"])
        FinalTUOneProductStart=as.character(TargetTUsPrediction[i,"GeneProduct2"])
        FinalTUOneNameEnd=as.character("")
        FinalTUOneProductEnd=as.character("")
        if(i==nrow(TargetTUsPrediction)){
          FinalTUOnePosStart=SelectedGenes[nrow(TargetTUsPrediction)+1,2]
        }else{
          FinalTUOnePosStart=TargetTUsPrediction[i+1,4]
        }
        FinalTUOnePosEnd=TargetTUsPrediction[i,5]
        FinalTUOne=data.frame(FinalTUOneNameStart,FinalTUOneNameEnd,FinalTUOnePosStart,FinalTUOnePosEnd,FinalTUOneProductStart,FinalTUOneProductEnd)
        flag=-1
      }
    }#for

    FinalTUAll=rbind(FinalTUAll,FinalTUOne)
    FinalTUAll$allNames=paste(FinalTUAll[,"FinalTUOneNameStart"],FinalTUAll[,"FinalTUOneNameEnd"],sep="")
    FinalTUAll$allProducts=paste(FinalTUAll[,"FinalTUOneProductStart"],FinalTUAll[,"FinalTUOneProductEnd"],sep="")
    FinalTUAll<-FinalTUAll[,c("FinalTUOnePosStart","FinalTUOnePosEnd","allNames","allProducts")]
    return(FinalTUAll)
  } #end of function

  # Generate Target TU SVM format
  GenerateOneTargetTU<-function(InputGeneList, i){
    OneTUList=data.frame()
    OneTUList=c("5'Gene", InputGeneList[i:(i+1),][1,2], InputGeneList[i:(i+1),][1,3], as.character(InputGeneList[i:(i+1),][1,4]));
    OneTUList=rbind(OneTUList, c("IntergenicRegion", (InputGeneList[i:(i+1),][1,3]+1), (InputGeneList[i:(i+1),][2,2]-1),as.character( InputGeneList[i:(i+1),][1,4])));
    OneTUList=rbind(OneTUList, c("3'Gene", as.numeric(as.character(InputGeneList[i:(i+1),][2,2])), InputGeneList[i:(i+1),][2,3], as.character(InputGeneList[i:(i+1),][2,4])));
    OneTUList=rbind(OneTUList, c("FullTU", InputGeneList[i:(i+1),][1,2], InputGeneList[i:(i+1),][2,3], as.character(InputGeneList[i:(i+1),][1,4])));
    row.names(OneTUList)=NULL
    OneTUList=as.data.frame(OneTUList)
    OneTUList[,2]=as.numeric(as.character(OneTUList[,2]))
    OneTUList[,3]=as.numeric(as.character(OneTUList[,3]))
    colnames(OneTUList)=c("ID", "Start", "End", "Strand") # Set Column Names
    return(OneTUList)
  }

  # Generate the final TU tables
  GenerateFinalTUTable<-function(TargetTUsPrediction, SelectedGenes, TargetTUsPredictionFreq){
    colnames(TargetTUsPrediction) = "TUresult"
    TargetTUName1Name2StartEnd=data.frame()
    j=1
    for(i in 1:(nrow(SelectedGenes)-1)){
      TargetTUName1Name2StartEnd=rbind(TargetTUName1Name2StartEnd, cbind(SelectedGenes[i,1], SelectedGenes[(i+1),1],GenerateOneTargetTU(SelectedGenes,i)[4,2:3]))
      j=j+1
    }
    TargetTUsPrediction = (cbind(TargetTUsPrediction, TargetTUName1Name2StartEnd))
    FinalTUAll=data.frame()
    FinalTUOne=data.frame()
    FinalTUOneNameStart=as.character(TargetTUsPrediction[1,2])
    FinalTUOneNameEnd=as.character("")
    FinalTUOnePosStart=TargetTUsPrediction[1,4]
    FinalTUOnePosEnd=TargetTUsPrediction[1,5]
    FinalTUOneFreq=0
    FinalTUOne=data.frame(FinalTUOneNameStart,FinalTUOneNameEnd,FinalTUOnePosStart,FinalTUOnePosEnd,FinalTUOneFreq)
    flag=TargetTUsPrediction[1,1]
    geneNum=0
    for(i in c(1:nrow(TargetTUsPrediction))){
      if(TargetTUsPrediction[i,1]==1){  #TU result = 1
        FinalTUOne[,2]=paste(FinalTUOne[,2],TargetTUsPrediction[i,3], sep=" ")
        FinalTUOne[,4]=TargetTUsPrediction[i,5]
        FinalTUOne[,5] = FinalTUOne[,5] + TargetTUsPredictionFreq[i,1]
        flag=1
        geneNum = geneNum+1
        if (i==nrow(TargetTUsPrediction)){
          FinalTUOne[,5] = (FinalTUOne[,5])/(geneNum)
        }
      }else{                                   #TU result = -1
        FinalTUOne[,5] = (FinalTUOne[,5] + TargetTUsPredictionFreq[i,1])/(geneNum+1)
        FinalTUAll=rbind(FinalTUAll,FinalTUOne)
        FinalTUOne=NULL
        flag=-1
        geneNum = 1
        FinalTUOneFreq = TargetTUsPredictionFreq[i,1]
        FinalTUOneNameStart=as.character(TargetTUsPrediction[i+1,2])
        FinalTUOneNameEnd=as.character("")
        FinalTUOnePosStart=TargetTUsPrediction[i+1,4]
        FinalTUOnePosEnd=TargetTUsPrediction[i,5]
        FinalTUOne=data.frame(FinalTUOneNameStart,FinalTUOneNameEnd,FinalTUOnePosStart,FinalTUOnePosEnd,FinalTUOneFreq)
      }
    }
    if(FinalTUOne[,3] == TRUE && FinalTUOne[,4] == TRUE){
      FinalTUAll=rbind(FinalTUAll,FinalTUOne)
    }
    return(FinalTUAll)
  }

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
  if ((AllForwardGenesExpressionMean[i] >= MinimumMeanExpressionValue) & (AllForwardGenesExpressionMedian[i] >= MinimumMedianExpressionValue) & (AllForwardGenesGapsProportion[i] < MaximumGapProportionValue) & (AllForwardGenesGapsLongest[i] < 50)){
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
  if ((AllReverseGenesExpressionMean[i] >= MinimumMeanExpressionValue) & (AllReverseGenesExpressionMedian[i] >= MinimumMedianExpressionValue) & 
      (AllReverseGenesGapsProportion[i] < MaximumGapProportionValue) & (AllReverseGenesGapsLongest[i] < 50)){
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
