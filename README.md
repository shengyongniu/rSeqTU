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
Please check vignette.

## Demo
[[https://github.com/s18692001/rSeqTU/blob/master/IGV_TU_demo.png]]
