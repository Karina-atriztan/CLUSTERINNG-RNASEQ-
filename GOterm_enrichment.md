# GO term enrichment using Top GO version 2.38.1 in Bioconductor
```
library(topGO)
```
To perform GO enrichment we needed:
1. Differential expressed genes (myInteresingGenes)
2. GO term annotation of T. atrovirideIMI206040 (geneID2GO)
3. List of GO and genes taht will use as Universe (geneNames)

We selected our interesting genes by K-means method, and herarchical clustering.

1. We call our DEG list
```
DEGlist1= read.table ("DEGlist1.txt",header=T, row.names = 1, sep="\t")

```
2. Now we call the GO annotation of the organism using the function "readMappings"

```
geneID2GO <- readMappings("GOannotation_T.atroviride.txt",sep="\t", IDsep = ";")
```

3. To obtain a gene universe we used the list of the  point 2

```
geneNames <- names (geneID2GO) 
```

4. Now we selected the names if the Interesting genes
```
myInterestingGenes1= rownames(DEGlist1)
```

5. Next step is to search the interesting genes in the gene Universe
```
geneList1 <- factor(as.integer(geneNames %in% myInterestingGenes1))
names(geneList1) <- geneNames
```

6. We reated the GOdata object that contain lists and tables above created
```
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList1, annot =annFUN.gene2GO, gene2GO = geneID2GO) 
Interesting-genes=allGenes; annot(anotacion de temrinos GO-genes), gene2GO:universo de genes (los mismosque la anotación)

resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(GOdata, classic = resultFis,orderBy = "weight", ranksOf = "classic", topNodes = 50)

