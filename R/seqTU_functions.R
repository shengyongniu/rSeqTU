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



