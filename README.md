# RNA-seq-data-analysis
Here I describe the significant steps for T. atroviride RNA-seq data analysis using in the paper Atriztán-Hernández K. and Herrera-Estrella A. 2021.

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
module load FastQC/0.11.2
fastqc /path-to-.fastq.gz* --outdir /path-to-outputdirectory

#Trimmomatic
module load Trimmomatic/0.32
java -jar $TRIMMOMATIC SE -phred33 “path.to-.fastq.gz” “outputdirectory” ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

#HISAT2
module load  hisat2/2.1.0
#create index
hisat2-build -p 8 -f “genome.fasta”   “Output-directory”

hisat2 -q  -x “path to index” -U “path to trimed-.fq.gz” --dta --dta-cufflinks  -S “path-to output-.sam”
![image](https://user-images.githubusercontent.com/62311305/127053116-b0b7f745-433e-4953-9967-5f2bdee03d12.png)


```

##### Mapped reads counting and Differential expressed Analysis (DE).
For mapping reads counting to each gene HTseq version 0.14.1. using the next code:
```
module load HTSeq/0.12.4
htseq-count -i ID --nonunique all -m intersection-strict  --additional-attr Parent  -t mRNA   -f sam  “path-to-.sam-files”  “genome-annoatation-.gff3”  > results_.txt



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
grp = c("libraries by licates")
      
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

## K-means and herarchical clusering 

To analaze if there are spcific clusters for time or condtions we performed a herarchical clustering and extract clusters using the K-means method

### Bioconductor Packages used

BiocManager::install("ggplot2")
BiocManager::install("ggdendro")
BiocManager::install("reshape2")
BiocManager::install("factoextra")
BiocManager::install("purrr")
BiocManager::install("dendextend")

```
data.exprs <- read.table(file = "all.txt", header = TRUE, sep="\t", comment.char="", row.names = 1)
```

This step  is for to scalate the data expresion because genes with high expression can affect the posterior analysis

```
data.exprs  <- as.matrix(data.exprs)
data.exprs  <- as.matrix(data.exprs)
scaledata <- t(scale(t(data.exprs))) # Centers and scales data.
scaledata <- scaledata[complete.cases(scaledata),]

```

Once that data is scalated we performed the clusterization of rows and columns using the herarchical clustering + Pearson and Spearman correlation by cluster

```
hr <- hclust(as.dist(1-cor(t(scaledata), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
hc <- hclust(as.dist(1-cor(scaledata, method="spearman")), method="complete") # Clusters columns by Spearman correlation.
```

## here we use a combinmation of herarchical clusterin and k-means.

Here I’m asking for  (k=65): THis, because we before observed that 65 is the best option 
to analyze the clusters in our sample, it depends of size ans complexity of the sample

```

hclustk4 = cutree(hr, k=65)
tiff("dendrogram-k65.tiff", width = 6, height = 8, res= 100, units = "in") #Para exportar la imagen en alta calidad

plot(TreeR,
     leaflab = "none",
     main = "Gene Clustering",
     ylab = "Height")
colored_bars(hclustk4, TreeR, sort_by_labels_order = T,
             y_shift=-0.1, rowLabels = c("k=65"),cex.rowLabels=1.2)
dev.off()

```

Set the minimum module size
minModuleSize = 50;


### NOW WE GO TO  PERFORMED A CLUSTERING USING THE WGCNA PACKAGE TO EXTRACT THE FORMED CLUSTERS TO FUTURE ANALYSIS

Module identification using dynamic tree cut

```
BiocManager::install("WGCNA")
library(WGCNA)
````
WE WILL TO USE THE SAME DENDROGRAM BEFORE CONSTRUCTED AND WILL  CUTTED USING CUTREEDYNAMIC
```
dynamicMods = cutreeDynamicTree(dendro = hr, maxTreeHeight = 0.5, deepSplit = TRUE, minModuleSize = 25);
```

The following command gives the module labels and the size of each module. Lable 0 is reserved for unassigned genes
table(dynamicMods)

```

dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
write.table(table(dynamicColors), file="keycolors65K.txt", sep="\t")
plotDendroAndColors(hr, dynamicColors, "k=65", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and module colors")

````
Discard the unassigned genes, and focus on the rest
restGenes= (dynamicColors != "grey")

Extract modules 

```
gene.names =row.names(data.exprs)
n=5660
SubGeneNames =gene.names [1:n]
module_colors= setdiff(unique(dynamicColors), "grey")
for (color in module_colors){
  module= SubGeneNames[which(dynamicColors==color)]
  write.table(module, paste("module_",color, ".txt",sep=""), sep="\t", row.names=FALSE, col.names=FALSE,quote=FALSE)
}

```

If your K number produces clusters with high correlation (say above 0.85) then consider reducing the number of clusters.

```
cor(kClustcentroids)

```

To calculate the scores for a single cluster, in this case 2 we’ll extract the core data for cluster 2,then subset the scaled data by cluster =2. Then, we’ll calculate the ‘score’ by correlating each gene#with the cluster core. We can then plot the results for each gene with the core overlayed:
Subset the cores molten dataframe so we can plot the core

```
core2 <- Kmolten[Kmolten$cluster=="11",]
```
get cluster 2 or other of of interest
```
K2 <- (scaledata[kClusters==11,])
str(K2)
```
Save the table of the genes for the core
```
write.table(K2, file="ID-cluster33.txt", sep="\t")

```

Calculate the correlation with the core

```
corscore <- function(x){cor(x,core2$value)}
score <- apply(K2, 1, corscore)

```

Get the data frame into long format for plotting
```
K2molten <- melt(K2)
colnames(K2molten) <- c('gene','sample','value')
```

Add the score

```
K2molten <- merge(K2molten,score, by.x='gene',by.y='row.names', all.x=T)
colnames(K2molten) <- c('gene','sample','value','score')
```

Order the dataframe by score, to do this first create an ordering factor
```
K2molten$order_factor <- 1:length(K2molten$gene)
```
Order the dataframe by score

```
K2molten <- K2molten[order(K2molten$gene),]
```

set the order by setting the factors
```
K2molten$order_factor <- factor(K2molten$order_factor , levels = K2molten$order_factor)
```

