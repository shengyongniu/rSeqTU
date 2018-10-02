# rSeqTU
An R package for prediction of Transcription Unit in Bacteria

# rSeqTU - a machine learning based R package for identification of bacterial transcriptional units
These R scripts will be wrapped up as an R package. Users will be able to install this package from Bioconductor.

## Introduction
A transcriptional unit (TU) is composed of one or multiple consecutive genes on the same strand of a
bacterial genome. The genes within a TU are transcribed into a single mRNA to respond to specific
growth conditions, and the TU are regulated by one promoter. To delineate the transcriptional regulatory
networks, it is a crucial first step to identify TUs within a bacterial genome accurately. To allow users to
efficiently perform TU identification on their machine and provide a more accurate prediction, we
develop an R package, named rSeqTU. rSeqTU R package can automatically select essential TU features
through a random forest algorithm in a machine learning framework. Besides, rSeqTU performs nearly
98% accuracy in most of our testing TU cases, e.g., public RNA-Seq datasets of E. coli. Users will be able to 
install the rSeqTU package from Bioconductor and input their customized RNA-Seq dataset to conduct
de-multiplexing, quality controlling, reads alignment, random-forest-based feature selection, prediction
model training, and TU prediction. Moreover, rSeqTU presents results in graphical visualizations and
interactive tables for customized downstream analysis. rSeqTU also output read count matrix of both
genes and TUs for further differentially expression analysis.

## Quick Start
library(Rsubread)
library(Rsamtools)
library(QuasR)
library(e1071)
library(seqinr)
library(grid)
library(gridBase)
library(ggplot2)
library(reshape2)
library(plyr)
library(caret)
library(mlbench)
library(Gviz)
library(GenomicRanges)


setwd("your_working_dir")
sampleFile <- "sampleFile.txt"
genomeFile <- "GCF_000005845.2_ASM584v2_genomic.fna"
genome_gff <- "GCF_000005845.2_ASM584v2_genomic.gff"
proj <- qAlign(sampleFile, genomeFile, paired="fr", clObj = makeCluster(detectCores()))

bam_M9.files <- "M9Enrich_S1_L001_R1_001_1241d6627572d.bam"
bam_Ri.files <- "RiEnrich_S3_L001_R1_001_1241d2758dc5e.bam"

asSam(bam_M9.files, destination=sub("\\.bam", "", bam_M9.files))
asSam(bam_Ri.files, destination=sub("\\.bam", "", bam_Ri.files))
sam_M9.file <- "M9Enrich_S1_L001_R1_001_1241d6627572d.sam"
sam_Ri.file <- "RiEnrich_S3_L001_R1_001_1241d2758dc5e.sam"


qQCReport(proj, pdfFilename="qc_report10.pdf")
alignmentStats(proj)


# Extract quality scores
qs <- qualityScores(filename="SRR400619.fastq",nreads=100)
# Check dimension of qs
dim(qs)
# Check first few elements of qs with head
head(qs)
boxplot(qs)

# FeatureCounts
# non-strand speciifc (duplicate column one and col two // and strand specific
fc_M9 <- featureCounts(bam_M9.files, annot.ext = "mocked.GTF", isGTFAnnotationFile=TRUE,
                       GTF.featureType  = "CDS", GTF.attrType="position_ID", countMultiMappingReads = TRUE, strandSpecific = 0)
fc_Ri <- featureCounts(bam.files, annot.ext = "mocked.GTF", isGTFAnnotationFile=TRUE,
                       GTF.featureType  = "CDS", GTF.attrType="position_ID", countMultiMappingReads = TRUE, strandSpecific = 0)

# See what slots are stored in fc
names(fc)
# Take a look at the featurecounts stats
fc$stat
# Take a look at the dimensions to see the number of genes
dim(fc$counts)
# Take a look at the first 6 lines
head(fc$counts)

# Take a look at the first 6 lines
head(fc$counts)
head(fc$annotation)

save(fc, "fc.RData")

generate_GTF_training(4641628)
gen_NA("Ecoli_alignment.BAM", "NC_000913.3", "Ecoli")
gen_cTU_data("Ri.NA", "test", genome_gff, genomeFile)
TU_SVM("SimulatedPositiveTUMatrix.txt", "SimulatedNegativeTUMatrix.txt", "TargetPositiveTUMatrix.txt", "TargetNegativeTUMatrix.txt",
      "Ri.NA", genome_gff, "Ri", "NC_000913")


