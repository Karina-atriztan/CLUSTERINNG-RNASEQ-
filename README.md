# RNA-seq-data-analysis
Here I describe the significant steps for fungal RNA-seq data Analysis

# Objetive
Main of these repository is to concatenate tested software and pipileines for fungal  RNA-seq Data analysis.

# Biological phenomena analyzed 
I desribe the followed strategies transcriptional comparison of the response of the fungus Trichoderma atroviride to injury and predation in a course time experiment

# Sequencing Technology used and libraries obtained
42 libraries were sequencing througth Illumina TruSeq 1X75 single-end.


# Quality analysis and reads mapping

Quality of the RNA-seq data were analyzed by FastQC Version 0.11.6. Around 10 millions of high-quality read per library were obtained.
Cleaned reads were mapped to the new genome reference of T. atroviride IMI206040 (Atriztán et al., in prep) using HISAT2 version 2.1.0. Specifications are desribed in repository 1

# Mapped reads counting and Differential expressed Analysis (DE).
For mapping reads counting to each gene HTseq version 0.14.1. Specifications are described in repository 1

To performed DE we use the Bioconductor package EdgeR and Limma.
Libraries Preprocessing:
1. Filter genes with more than five reads per libraries and present in at least two different libraries: 
(counts = counts[rowSums(cpm(counts) >= 2) >=5,])
2. 
