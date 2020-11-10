# RNA-seq-data-analysis
Here I describe the significant steps for T. atroviride RNA-seq data analysis using in the paper Atriztán-Hernández K. and Herrera-Estrella A. in press.

## Objetive
Main of these repository is to concatenate tested software and pipileines for T. atroviride  RNA-seq Data analysis.

## Biological phenomena analyzed 
I desribe the followed strategies transcriptional comparison of the response of the fungus Trichoderma atroviride to mechanical injury and predation under a course time experiment

##### Sequencing Technology used and libraries obtained
42 libraries were sequencing througth Illumina TruSeq 1X75 single-end.


##### Quality analysis and reads mapping

Quality of the RNA-seq data were analyzed by FastQC Version 0.11.6. Around 10 millions of high-quality read per library were obtained.
Cleaned reads were mapped to the new genome reference of T. atroviride IMI206040 (Atriztán-Hernández et al., in prep) using HISAT2 version 2.1.0. Output files were converted to BAM files for visualization on desktop app IGV://software.broadinstitute.org/software/igv/home). We use the code as next:

```
hisat2 -q  -x HISAT_mapping/TAIMI -U trimmed/library1.fq.gz --dta --dta-cufflinks  -S HISAT_mapping/SAM/library1.sam
samtools view -bS HISAT_mapping/SAM/library1.sam > HISAT_mapping/SAM/library1.bam
samtools sort HISAT_mapping/SAM/library1.bam -o  HISAT_mapping/SAM/library1.sorted.bam

```

##### Mapped reads counting and Differential expressed Analysis (DE).
For mapping reads counting to each gene HTseq version 0.14.1. using the next code:
```
htseq-count -i ID --nonunique all -m intersection-strict  --additional-attr Parent  -t mRNA   -f sam  HISAT_mapping/SAM/library1.sam  Trichoderma_atrovirideIMI206040.gff3  > results_library1.txt

```


## Functional annotation of T. atroviride proteome
To perform GO term enrichment and selection of specific genes, we decide to "de novo" annotate the T. atroviride proteome. Followint The InterProScan/5.41 pipiline we run the next code:

```
interproscan.sh  --goterms --pathways --iprlookup -i Trichoderma_atroviride_IMI206040.proteins.fa -b output/IprScan
```

# DEG analysis
To performed DE we use the Bioconductor package EdgeR and Limma.
#### Libraries Preprocessing:

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
7. Now, we can compare any of the treatment groups using the contrast argument of the glmQLFTest
example:
```
qlf.Treatment1vscontrol =glmQLFTest(fit, contrast=my.contrasts[,"Treatment1vscontrol"])
topTags(qlf.Treatment1vscontrol)
dim(qlf.Treatment1vscontrol)

```
8. Now we filter the significat DEG using a DFR value >0.005 and save the corresponding table.

```
deTabTreatment1vscontrol = deTab = topTags(qlf.Treatment1vscontrol, n=Inf)$table
deTabTreatment1vscontrol[c(2),]
deGenesTreatment1vscontrol = rownames(deTabTreatment1vscontrol)[deTabTreatment1vscontrol$FDR < 0.005]
length(deGenesTreatment1vscontrol)
 # To obtain a plot with up and down regulated genes 
plotSmear(dge[,c("columns where are the traetments and control")], de.tags=deGenesTreatment1vscontrol, cex=0.7,xlab="Average logCPM", ylab="logFC")
abline(h=c(-1,1), col='blue')        
#To save the DEG
write.table(deTabTreatment1vscontrol[deGenesTreatment1vscontrol,], file=paste(outpathcount2, "DEGTreatment1vscontrol.txt", sep=""), row.names=TRUE, col.names=NA, quote=FALSE, sep="\t")


```






