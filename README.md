# Introduction
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

# Result Demo
<img src="IGV_TU_demo.png" width="1200" />

# Quick Start

### Install package
```R
library(devtools)
install_github("s18692001/rSeqTU")
library(rSeqTU)
```

### Alignment
```R
library(QuasR)
library(Rsamtools)
# Set up your working directory
setwd("working_dir")

# Set up parameters of file paths for alignment
sampleFile <- "sampleFile.txt"
genomeFile <- "GCF_000005845.2_ASM584v2_genomic.fna"
genome_gff <- "GCF_000005845.2_ASM584v2_genomic.gff"
proj <- qAlign(sampleFile, genomeFile, paired="fr", clObj = makeCluster(detectCores()))
```

```
    create 1 genomic alignment(s)
will start in ..9s..8s..7s..6s..5s..4s..3s..2s..1s
Testing the compute nodes...OK
Loading QuasR on the compute nodes...OK
Available cores:
nodeNames
  syniu 
     16 
Performing genomic alignments for 1 samples. See progress in the log file:
local/QuasR_log_1ccc01b031f84.txt
Genomic alignments have been created successfully
```

### Quality Check
```R
# get Quality Check report and statistics
qQCReport(proj, pdfFilename="qc_report_test.pdf")
alignmentStats(proj)
```
<img src="QC/QC1.png" width="800" />
<img src="QC/QC2.png" width="800" />
<img src="QC/QC3.png" width="800" />
<img src="QC/QC4.png" width="800" />
<img src="QC/QC5.png" width="800" />
<img src="QC/QC6.png" width="800" />
<img src="QC/QC7.png" width="800" />
<img src="QC/QC8.png" width="800" />


### .NA file generation
```R
# Generate .NA file for constructing cTU
gen_NA("M9Enrich_S2_L001_R1_001_1ccc063fe8630.bam", "NC_000913.3", "M9")
gen_cTU_data("M9.NA", "test", genome_gff, genomeFile)
```

### SVM modeling and prediction
```R
# Train model and generate prediction result in .bedgraph format
TU_SVM("SimulatedPositiveTUMatrix.txt", "SimulatedNegativeTUMatrix.txt", "TargetPositiveTUMatrix.txt", "TargetNegativeTUMatrix.txt","M9.NA", genome_gff, "M9", "NC_000913")
```

### Visualization in IGV (Import the .bedgraph, reference genome and annotation files)
<img src="IGV_TU_demo.png" width="1200" />






