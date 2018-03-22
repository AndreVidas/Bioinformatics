rm(list=ls())


###################################################
### DEseq 2
###################################################
library(DESeq2)

#### load data

# experimental data
sampleTable_exp <- data.frame(sampleName = c("Galaxy92-[htseq-count_on_sg26_4d_ctrl_R1]", "Galaxy94-[htseq-count_on_sg26_4d_ctrl_R2]", "Galaxy96-[htseq-count_on_sg26_4d_ctrl_R3]", "Galaxy98-[htseq-count_on_sg26_4d_dox_R1]", "Galaxy100-[htseq-count_on_sg26_4d_dox_R2]", "Galaxy102-[htseq-count_on_sg26_4d_dox_R3]"), fileName = c("Galaxy92-[htseq-count_on_sg26_4d_ctrl_R1].tabular", "Galaxy94-[htseq-count_on_sg26_4d_ctrl_R2].tabular", "Galaxy96-[htseq-count_on_sg26_4d_ctrl_R3].tabular", "Galaxy98-[htseq-count_on_sg26_4d_dox_R1].tabular", "Galaxy100-[htseq-count_on_sg26_4d_dox_R2].tabular", "Galaxy102-[htseq-count_on_sg26_4d_dox_R3].tabular"), condition = c("ctrl", "ctrl", "ctrl", "trt", "trt", "trt"))
ddsHTSeq_exp <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_exp,
                                          directory = "dataChunqin/",
                                          design= ~ condition)
head(counts(ddsHTSeq_exp))


# negative control
sampleTable_NC <- data.frame(sampleName = c("Galaxy79-[htseq-count_on_sgNCI_4d_ctrl_R1]", "Galaxy81-[htseq-count_on_sgNCI_4d_ctrl_R2]", "Galaxy83-[htseq-count_on_sgNCI_4d_ctrl_R3]", "Galaxy85-[htseq-count_on_sgNCI_4d_dox_R1]", "Galaxy87-[htseq-count_on_sgNCI_4d_dox_R2]", "Galaxy89-[htseq-count_on_sgNCI_4d_dox_R3]"), fileName = c("Galaxy79-[htseq-count_on_sgNCI_4d_ctrl_R1].tabular", "Galaxy81-[htseq-count_on_sgNCI_4d_ctrl_R2].tabular", "Galaxy83-[htseq-count_on_sgNCI_4d_ctrl_R3].tabular", "Galaxy85-[htseq-count_on_sgNCI_4d_dox_R1].tabular", "Galaxy87-[htseq-count_on_sgNCI_4d_dox_R2].tabular", "Galaxy89-[htseq-count_on_sgNCI_4d_dox_R3].tabular"), condition = c("ctrl", "ctrl", "ctrl", "trt", "trt", "trt"))
ddsHTSeq_NC <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_NC,
                                       directory = "dataChunqin/",
                                       design= ~ condition)
head(counts(ddsHTSeq_NC))




# DEA exp. data
ddsHTSeq_exp <- DESeq(ddsHTSeq_exp)
res_exp <- results(ddsHTSeq_exp, contrast=c("condition","trt","ctrl"))
resultsNames(ddsHTSeq_exp)
resLFCshrinked_exp <- lfcShrink(ddsHTSeq_exp, coef=2)

# MA plot
pdf(file = "dataChunqin/results/MA_plots/expData.pdf")
plotMA(resLFCshrinked_exp, ylim=c(-2.6,3.9))
dev.off()


# remove NA's from results
resLFCshrinked_exp <- na.exclude(resLFCshrinked_exp)


# up-regulated exp
up_exp <- resLFCshrinked_exp[resLFCshrinked_exp$padj < 0.05 & resLFCshrinked_exp$log2FoldChange > 1,]
dim(na.exclude(up_exp[order(up_exp$padj),]))
up_exp[order(up_exp$padj),]
up_exp[order(-up_exp$log2FoldChange),]

# down-regulated exp
down_exp <- resLFCshrinked_exp[resLFCshrinked_exp$padj < 0.05 & resLFCshrinked_exp$log2FoldChange < -1,]
dim(na.exclude(down_exp[order(down_exp$padj),]))
down_exp[order(down_exp$padj),]
down_exp[order(down_exp$log2FoldChange),]



# DEA negative control
ddsHTSeq_NC <- DESeq(ddsHTSeq_NC)
res_NC <- results(ddsHTSeq_NC, contrast=c("condition","trt","ctrl"))
resultsNames(ddsHTSeq_NC)
resLFCshrinked_NC <- lfcShrink(ddsHTSeq_NC, coef=2)

# MA plot
pdf(file = "dataChunqin/results/MA_plots/NC.pdf")
plotMA(resLFCshrinked_NC, ylim=c(-0.6,0.7))
dev.off()

# negative control MA plot for comparison exp. data MA plot
pdf(file = "dataChunqin/results/MA_plots/NC_copmparison.pdf")
plotMA(resLFCshrinked_NC, ylim=c(-2.6,3.9))
dev.off()


# remove NA's from results
resLFCshrinked_NC <- na.exclude(resLFCshrinked_NC)


# up-regulated exp
up_NC <- resLFCshrinked_NC[resLFCshrinked_NC$padj < 0.05 & resLFCshrinked_NC$log2FoldChange > 1,]
dim(na.exclude(up_NC[order(up_NC$padj),]))
up_NC[order(up_NC$padj),]
up_NC[order(-up_NC$log2FoldChange),]

# down-regulated exp
down_NC <- resLFCshrinked_NC[resLFCshrinked_NC$padj < 0.05 & resLFCshrinked_NC$log2FoldChange < -1,]
dim(na.exclude(down_NC[order(down_NC$padj),]))
down_NC[order(down_NC$padj),]
down_NC[order(down_NC$log2FoldChange),]





# remove DE genes in exp. data that intersects with negative control
geneIntersection_up <- intersect(rownames(up_exp), rownames(up_NC))
geneIntersection_down <- intersect(rownames(down_exp), rownames(down_NC))


dim(up_exp)
up_exp <- up_exp[!rownames(up_exp) %in% geneIntersection_up,]
dim(up_exp)

dim(down_exp)
down_exp <- down_exp[!rownames(down_exp) %in% geneIntersection_down,]
dim(down_exp)


#### save DE lists
# result list without thresholds
write.csv(resLFCshrinked_exp, file = "dataChunqin/results/DE_gene_list/exp.csv", quote = FALSE)

# p adj. < 0.05 & logFC threshold (-)1
write.csv(up_exp, file = "dataChunqin/results/DE_gene_list/up_exp_p_0d05_logFC_1.csv", quote = FALSE)
write.csv(down_exp, file = "dataChunqin/results/DE_gene_list/down_exp_p_0d05_logFC_1.csv", quote = FALSE)


# p adj. < 0.05 & logFC threshold (-)2
# up-regulated exp
up_exp <- resLFCshrinked_exp[resLFCshrinked_exp$padj < 0.05 & resLFCshrinked_exp$log2FoldChange > 2,]
dim(na.exclude(up_exp[order(up_exp$padj),]))
up_exp[order(up_exp$padj),]
up_exp[order(-up_exp$log2FoldChange),]
# down-regulated exp
down_exp <- resLFCshrinked_exp[resLFCshrinked_exp$padj < 0.05 & resLFCshrinked_exp$log2FoldChange < -2,]
dim(na.exclude(down_exp[order(down_exp$padj),]))
down_exp[order(down_exp$padj),]
down_exp[order(down_exp$log2FoldChange),]


write.csv(up_exp, file = "dataChunqin/results/DE_gene_list/up_exp_p_0d05_logFC_2.csv", quote = FALSE)
write.csv(down_exp, file = "dataChunqin/results/DE_gene_list/down_exp_p_0d05_logFC_2.csv", quote = FALSE)



# p adj. < 0.05 & logFC threshold (-)1.5
# up-regulated exp
up_exp <- resLFCshrinked_exp[resLFCshrinked_exp$padj < 0.05 & resLFCshrinked_exp$log2FoldChange > 1.5,]
dim(na.exclude(up_exp[order(up_exp$padj),]))
up_exp[order(up_exp$padj),]
up_exp[order(-up_exp$log2FoldChange),]
# down-regulated exp
down_exp <- resLFCshrinked_exp[resLFCshrinked_exp$padj < 0.05 & resLFCshrinked_exp$log2FoldChange < -1.5,]
dim(na.exclude(down_exp[order(down_exp$padj),]))
down_exp[order(down_exp$padj),]
down_exp[order(down_exp$log2FoldChange),]

write.csv(up_exp, file = "dataChunqin/results/DE_gene_list/up_exp_p_0d05_logFC_1d5.csv", quote = FALSE)
write.csv(down_exp, file = "dataChunqin/results/DE_gene_list/down_exp_p_0d05_logFC_1d5.csv", quote = FALSE)







sessionInfo()
