# ref: http://www.sthda.com/english/wiki/ggplot2-texts-add-text-annotations-to-a-graph-in-r-software
rm(list=ls())
library(ggplot2)
library(ggrepel)

#### load data
genes <- read.csv2(file = "../results/DEgenes/codingGenes/postHocFiltering/trt_vs_ctrl.csv", row.names = 1)
#genes <- read.csv2(file = "../results/DEgenes/codingGenes/postHocFiltering/relap_vs_trt.csv", row.names = 1)
#genes <- read.csv2(file = "../results/DEgenes/codingGenes/postHocFiltering/relap_vs_ctrl.csv", row.names = 1)

#genes <- read.csv2(file = "../results/DEgenes/codingGenes/lfcThreshold/trt_vs_ctrl.csv", row.names = 1)
#genes <- read.csv2(file = "../results/DEgenes/codingGenes/lfcThreshold/relap_vs_trt.csv", row.names = 1)
#genes <- read.csv2(file = "../results/DEgenes/codingGenes/lfcThreshold/relap_vs_ctrl.csv", row.names = 1)

dim(genes)
genes <- genes[!is.na(genes$padj),]
dim(genes)
nrow(genes[genes$padj < 0.05 & (genes$log2FoldChange < -1 | genes$log2FoldChange > 1),])


#genes$Gene <- rownames(genes) 
head(genes)
alphaLevel <- 0.05 
geneNamesThres <- 1e-100
logFC_thres <- 1 # when changing this, remember to change the FC legends in the plot
top_up_down <- 0
pAdjMethod <- "FDR"


# remove na
dim(genes)
genes <- genes[!is.na(genes$padj),]
dim(genes)

genes$Significance <- NA
genes[genes$padj > alphaLevel, "Significance"] <- "Not signicant"
genes[genes$padj < alphaLevel & ( genes$log2FoldChange <= logFC_thres & genes$log2FoldChange >= -logFC_thres), "Significance"] <- paste0("Sig., abs(log2(FC)) < ", logFC_thres)
#genes[genes$padj < alphaLevel & ( genes$log2FoldChange > logFC_thres | genes$log2FoldChange < -logFC_thres), "Significance"] <-  paste0("Sig., abs(log2(FC)) > ", logFC_thres)
genes[genes$padj < alphaLevel & ( genes$log2FoldChange > logFC_thres), "Significance"] <-  paste0("Sig., log2(FC) > ", logFC_thres)
genes[genes$padj < alphaLevel & ( genes$log2FoldChange < -logFC_thres), "Significance"] <-  paste0("Sig., log2(FC) < -", logFC_thres)
head(genes, n = 45)
dim(genes)


# find 10 most up and down regulated genes
up_top <- head(genes[order(-genes$log2FoldChange),], n = top_up_down)
down_top <- head(genes[order(genes$log2FoldChange),], n = top_up_down)

pdf(file = "../results/volcano/codingGenes/noLabels/trt_vs_ctrl.pdf", width = 4, height = 4)
#pdf(file = "../results/volcano/codingGenes/noLabels/relap_vs_trt.pdf", width = 4, height = 4)
#pdf(file = "../results/volcano/codingGenes/noLabels/relap_vs_ctrl.pdf")
ggplot(genes, aes(x = log2FoldChange, y = -log10(padj))) + xlab(expression(log[2]("FC"))) +
    ylab(bquote(-log[10](.(pAdjMethod)))) + geom_point(aes(color = Significance)) +
    #scale_color_manual(values = c("Sig., abs(log2(FC)) > 1" = "red", "Not signicant" = "grey", "Sig., abs(log2(FC)) < 1" = "blue")) + # remember to change threshold valuse of logFC here for legend
    scale_color_manual(values = c("Sig., log2(FC) > 1" = "red", "Sig., log2(FC) < -1" = "blue", "Not signicant" = "grey", "Sig., abs(log2(FC)) < 1" = "black")) + # remember to change threshold valuse of logFC here for legend
    theme_bw(base_size = 12) + theme(legend.position = "bottom") +
    geom_hline(yintercept = -log10(alphaLevel), colour="red", linetype="dashed") + geom_vline(xintercept = logFC_thres, colour="red", linetype="dashed") + geom_vline(xintercept = -logFC_thres, colour="red", linetype="dashed") +
    geom_text_repel(
        #data = subset(genes, padj < geneNamesThres),
        data = rbind(up_top, down_top),
        aes(label = gene),
        size = 5,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
    )
dev.off()


sessionInfo()

