rm(list=ls())
library(org.Mm.eg.db) # Genome wide annotation for Mouse, primarily based on mapping using Entrez Gene identifiers. ref: https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html
library(topGO)

# load all genes from DESeq2 result data
geneList <- read.csv(file = "../../../projects/nfd/results/DE_genes/MllAF9_vs_NF.csv", row.names = 1, as.is = TRUE)
dim(geneList)


# convert geneList to topGO input format
geneUni <- geneList$padj
names(geneUni) <- rownames(geneList) 


# GO enrichment analysis
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneUni, geneSel = function(p) p < 0.01, description = "Mouse GO Enrichment Analysis", ID = "Ensembl", mapping = "org.Mm.eg.db", annot = annFUN.org)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
GenTable(GOdata, classicFisher = resultFisher, topNodes = 10)
