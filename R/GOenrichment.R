# maybe use both up and down regulated genes with all genes as geneUni
# maybe use clusterProfiler instead of topGO: https://www.r-bloggers.com/go-analysis-using-clusterprofiler/ & https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html


rm(list=ls())
library(org.Mm.eg.db) # Genome wide annotation for Mouse, primarily based on mapping using Entrez Gene identifiers. ref: https://bioconductor.org/packages/release/data/annotation/html/org.Mm.eg.db.html
library(topGO)

# load all genes from DESeq2 result data
#geneList <- read.csv2(file = "../results/DEgenes/codingGenes/postHocFiltering/trt_vs_ctrl.csv", row.names = 1, as.is = TRUE)
#geneList <- read.csv2(file = "../results/DEgenes/codingGenes/postHocFiltering/relap_vs_trt.csv", row.names = 1, as.is = TRUE)
geneList <- read.csv2(file = "../results/DEgenes/codingGenes/postHocFiltering/relap_vs_ctrl.csv", row.names = 1, as.is = TRUE)

#geneList <- read.csv2(file = "../results/DEgenes/codingGenes/lfcThreshold/trt_vs_ctrl.csv", row.names = 1, as.is = TRUE)
#geneList <- read.csv2(file = "../results/DEgenes/codingGenes/lfcThreshold/relap_vs_trt.csv", row.names = 1, as.is = TRUE)
#geneList <- read.csv2(file = "../results/DEgenes/codingGenes/lfcThreshold/relap_vs_ctrl.csv", row.names = 1, as.is = TRUE)

dim(geneList)
head(geneList)
# remove NAs for gene list
geneList <- geneList[!is.na(geneList$padj),]

# up-regulated gene universe
#geneList[geneList$log2FoldChange < 1,"padj"] <- 1
# down-regulated gene universe
geneList[geneList$log2FoldChange > -1,"padj"] <- 1




# convert geneList to topGO input format
geneUni <- geneList$padj
names(geneUni) <- rownames(geneList) 
length(geneUni)



# GO enrichment analysis
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneUni, geneSel = function(p) p < 0.01, description = "Mouse GO Enrichment Analysis", ID = "Ensembl", mapping = "org.Mm.eg.db", annot = annFUN.org)
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")


tableResult <- GenTable(GOdata, classicFisher = resultFisher, topNodes = 10000)
tableResult$fdr <- p.adjust(p = tableResult$classicFisher, method = "fdr")



# full gene universe
#write.csv2(x = tableResult, "../results/GOenrichment/codingGenes/postHocFiltering/trt_vs_ctrl_up.csv")
#write.csv2(x = tableResult, "../results/GOenrichment/codingGenes/postHocFiltering/relap_vs_trt_up.csv")
#write.csv2(x = tableResult, "../results/GOenrichment/codingGenes/postHocFiltering/relap_vs_ctrl_up.csv")

#write.csv2(x = tableResult, "../results/GOenrichment/codingGenes/postHocFiltering/trt_vs_ctrl_down.csv")
#write.csv2(x = tableResult, "../results/GOenrichment/codingGenes/postHocFiltering/relap_vs_trt_down.csv")
write.csv2(x = tableResult, "../results/GOenrichment/codingGenes/postHocFiltering/relap_vs_ctrl_down.csv")



#write.csv2(x = tableResult, "../results/GOenrichment/codingGenes/lfcThreshold/fullGeneUniverse/trt_vs_ctrl_up.csv")
#write.csv2(x = tableResult, "../results/GOenrichment/codingGenes/lfcThreshold/fullGeneUniverse/relap_vs_trt_up.csv")
#write.csv2(x = tableResult, "../results/GOenrichment/codingGenes/lfcThreshold/fullGeneUniverse/relap_vs_ctrl_up.csv")

#write.csv2(x = tableResult, "../results/GOenrichment/codingGenes/lfcThreshold/fullGeneUniverse/trt_vs_ctrl_down.csv")
#write.csv2(x = tableResult, "../results/GOenrichment/codingGenes/lfcThreshold/fullGeneUniverse/relap_vs_trt_down.csv")
#write.csv2(x = tableResult, "../results/GOenrichment/codingGenes/lfcThreshold/fullGeneUniverse/relap_vs_ctrl_down.csv")









#GOdata <- new("topGOdata", ontology = "CC", allGenes = geneUni, geneSel = function(p) p < 0.01, description = "Mouse GO Enrichment Analysis", ID = "Ensembl", mapping = "org.Mm.eg.db", annot = annFUN.org)
#resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#GenTable(GOdata, classicFisher = resultFisher, topNodes = 1000)


#GOdata <- new("topGOdata", ontology = "MF", allGenes = geneUni, geneSel = function(p) p < 0.01, description = "Mouse GO Enrichment Analysis", ID = "Ensembl", mapping = "org.Mm.eg.db", annot = annFUN.org)
#resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#GenTable(GOdata, classicFisher = resultFisher, topNodes = 1000)
