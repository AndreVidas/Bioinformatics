rm(list=ls())
setwd("~/Documents/projects/mieDitte/scripts/")

###################################################
### DEseq 2
###################################################
library(DESeq2)
library(biomaRt)

#### load data

# experimental data
# load mouse data
data <- read.table("../data/featureCounts/annotated_combined.counts", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
dim(data)
head(data)
colnames(data)


# filter mouse data so that only coding genes are used
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
ensmusg_geneName <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "gene_biotype"), filters='biotype', values=c('protein_coding'), mart = mouse)
data <- data[match(ensmusg_geneName$ensembl_gene_id, rownames(data), nomatch = FALSE),]
dim(data)




geneNames <- data.frame(ensmusg = rownames(data), data$symbol, stringsAsFactors = FALSE)

# remove gene symbol from data.frame
data <- data[,1:12]


###################################################
### DEA all samples
###################################################
# make deseq2 object
print(head(data))
condition <- c("ctrl","relap","relap","relap","ctrl","ctrl","ctrl", "trt", "trt", "trt", "trt", "relap")


coldata <- data.frame(condition = condition)
print(coldata)
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design= ~ condition)


dds <- DESeq(dds)
resultsNames(dds)

#### trt vs ctrl
# post hoc filtering version
resLFCshrinked_exp <- lfcShrink(dds, contrast=c("condition","trt","ctrl"))

# MA plot
pdf(file = "../results/MA_plots/codingGenes/postHocFiltering/trt_vs_ctrl.pdf")
plotMA(resLFCshrinked_exp, alpha = 0.05, ylim=c(-4,6))
dev.off()

# insert gene names into results
resLFCshrinked_exp$gene <- geneNames[match(rownames(resLFCshrinked_exp), geneNames$ensmusg),"data.symbol"]

# save shrinked results, ordered by adj. p-value
write.csv2(x = resLFCshrinked_exp[order(resLFCshrinked_exp$padj),], file = "../results/DEgenes/codingGenes/postHocFiltering/trt_vs_ctrl.csv")



# lfcThreshold version
resLFCshrinked_exp <- lfcShrink(dds, contrast=c("condition","trt","ctrl"), lfcThreshold = 1)

# MA plot
pdf(file = "../results/MA_plots/codingGenes/lfcThreshold/trt_vs_ctrl.pdf")
plotMA(resLFCshrinked_exp, alpha = 0.05, ylim=c(-4,6))
dev.off()

# insert gene names into results
resLFCshrinked_exp$gene <- geneNames[match(rownames(resLFCshrinked_exp), geneNames$ensmusg),"data.symbol"]

# save shrinked results, ordered by adj. p-value
write.csv2(x = resLFCshrinked_exp[order(resLFCshrinked_exp$padj),], file = "../results/DEgenes/codingGenes/lfcThreshold/trt_vs_ctrl.csv")






#### relap vs trt
# post hoc filtering version
resLFCshrinked_exp <- lfcShrink(dds, contrast=c("condition","relap","trt"))

# MA plot
pdf(file = "../results/MA_plots/codingGenes/postHocFiltering/relap_vs_trt.pdf")
plotMA(resLFCshrinked_exp, alpha = 0.05, ylim=c(-6,4))
dev.off()

# insert gene names into results
resLFCshrinked_exp$gene <- geneNames[match(rownames(resLFCshrinked_exp), geneNames$ensmusg),"data.symbol"]

# save shrinked results, ordered by adj. p-value
write.csv2(x = resLFCshrinked_exp[order(resLFCshrinked_exp$padj),], file = "../results/DEgenes/codingGenes/postHocFiltering/relap_vs_trt.csv")



# lfcThreshold version
resLFCshrinked_exp <- lfcShrink(dds, contrast=c("condition","relap","trt"), lfcThreshold = 1)

# MA plot
pdf(file = "../results/MA_plots/codingGenes/lfcThreshold/relap_vs_trt.pdf")
plotMA(resLFCshrinked_exp, alpha = 0.05, ylim=c(-6,4))
dev.off()

# insert gene names into results
resLFCshrinked_exp$gene <- geneNames[match(rownames(resLFCshrinked_exp), geneNames$ensmusg),"data.symbol"]

# save shrinked results, ordered by adj. p-value
write.csv2(x = resLFCshrinked_exp[order(resLFCshrinked_exp$padj),], file = "../results/DEgenes/codingGenes/lfcThreshold/relap_vs_trt.csv")







#### relap vs ctrl
# post hoc filtering version
resLFCshrinked_exp <- lfcShrink(dds, contrast=c("condition","relap","ctrl"))

# MA plot
pdf(file = "../results/MA_plots/codingGenes/postHocFiltering/relap_vs_ctrl.pdf")
plotMA(resLFCshrinked_exp, alpha = 0.05, ylim=c(-7,4))
dev.off()

# insert gene names into results
resLFCshrinked_exp$gene <- geneNames[match(rownames(resLFCshrinked_exp), geneNames$ensmusg),"data.symbol"]

# save shrinked results, ordered by adj. p-value
write.csv2(x = resLFCshrinked_exp[order(resLFCshrinked_exp$padj),], file = "../results/DEgenes/codingGenes/postHocFiltering/relap_vs_ctrl.csv")



# lfcThreshold version
resLFCshrinked_exp <- lfcShrink(dds, contrast=c("condition","relap","ctrl"), lfcThreshold = 1)

# MA plot
pdf(file = "../results/MA_plots/codingGenes/lfcThreshold/relap_vs_ctrl.pdf")
plotMA(resLFCshrinked_exp, alpha = 0.05, ylim=c(-7,4))
dev.off()

# insert gene names into results
resLFCshrinked_exp$gene <- geneNames[match(rownames(resLFCshrinked_exp), geneNames$ensmusg),"data.symbol"]

# save shrinked results, ordered by adj. p-value
write.csv2(x = resLFCshrinked_exp[order(resLFCshrinked_exp$padj),], file = "../results/DEgenes/codingGenes/lfcThreshold/relap_vs_ctrl.csv")







sessionInfo()

