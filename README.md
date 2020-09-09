# RNA-seq-data-analysis
Here I describe the significant steps for fungal RNA-seq data Analysis

## Objetive
Main of these repository is to concatenate tested software and pipileines for fungal  RNA-seq Data analysis.

## Biological phenomena analyzed 
I desribe the followed strategies transcriptional comparison of the response of the fungus Trichoderma atroviride to injury and predation in a course time experiment

##### Sequencing Technology used and libraries obtained
42 libraries were sequencing througth Illumina TruSeq 1X75 single-end.


##### Quality analysis and reads mapping

Quality of the RNA-seq data were analyzed by FastQC Version 0.11.6. Around 10 millions of high-quality read per library were obtained.
Cleaned reads were mapped to the new genome reference of T. atroviride IMI206040 (Atriztán et al., in prep) using HISAT2 version 2.1.0. Specifications are desribed in repository 1

##### Mapped reads counting and Differential expressed Analysis (DE).
For mapping reads counting to each gene HTseq version 0.14.1. Specifications are described in repository 1

To performed DE we use the Bioconductor package EdgeR and Limma.
Libraries Preprocessing:
1. Filter genes with more than five reads per libraries and present in at least two different libraries and data was converted to log2CMP: 
```
counts = counts[rowSums(cpm(counts) >= 2) >=5,]
logcounts= cpm(counts,log=TRUE)
```
2. Using filtered data, we contructed a DGElist, group the data per replicates using the "grp" function:
```
grp = c("C0","C0","C0","C0","C0",
        "I30","I30","I30",
      "I90","I90","I90",  
      "I4","I4","I4",
        "I8","I8","I8",
        "D30","D30","D30",
      "D90","D90","D90",
        "D4","D4","D4",
        "D8", "D8", "D8",
      "OF30","OF30","OF30",
      "OF90","OF90","OF90",
      "OF4","OF4","OF4",
      "OF8", "OF8", "OF8")
      
 dge = DGEList(counts=counts, group=grp)
```
3. We analized the libraries dispersion using a MDSplot:
```
plotMDS(grp, col=c(rep("blue",3), rep(")))
```
4. Normalization was performed used  the  "Normalization" function:
```
Normalization=calcNormFactors(dge)
```
5. Due that we will to perform a multifactorial analysis we need to calculate the  GLMTagwise disperssion
```
dge = estimateGLMTagwiseDisp(dge)
```

6. We used the GLM approach to compare treatments, to do this, we contructed a matrix (one column per group):
```
design = model.matrix(~0+group, data=dge$samples) 

colnames(design) = levels(dge$samples$group)
design
```
6. Now, we can compare any of the treatment groups using the contrast argument of the glmQLFTest

