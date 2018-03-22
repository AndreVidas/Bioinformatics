rm(list=ls())
library(DESeq2)
library(ComplexHeatmap)
library(circlize)


#### load count data
sampleTable_exp <- data.frame(sampleName = c("sg26_4d_ctrl_R1", "sg26_4d_ctrl_R2", "sg26_4d_ctrl_R3", "sg26_4d_dox_R1", "sg26_4d_dox_R2", "sg26_4d_dox_R3"), fileName = c("Galaxy92-[htseq-count_on_sg26_4d_ctrl_R1].tabular", "Galaxy94-[htseq-count_on_sg26_4d_ctrl_R2].tabular", "Galaxy96-[htseq-count_on_sg26_4d_ctrl_R3].tabular", "Galaxy98-[htseq-count_on_sg26_4d_dox_R1].tabular", "Galaxy100-[htseq-count_on_sg26_4d_dox_R2].tabular", "Galaxy102-[htseq-count_on_sg26_4d_dox_R3].tabular"), condition = c("ctrl", "ctrl", "ctrl", "trt", "trt", "trt"))
data_exp <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_exp,
                                       directory = "dataChunqin/",
                                       design= ~ condition)
data_exp <- as.matrix(counts(data_exp))
head(data_exp)



#### load DE gene lists
down_exp <- read.csv("dataChunqin/results/DE_gene_list/down_exp_p_0d05_logFC_2.csv", row.names = 1)
up_exp <- read.csv("dataChunqin/results/DE_gene_list/up_exp_p_0d05_logFC_2.csv", row.names = 1)
head(down_exp)
head(up_exp)

# merge up and down together
DE_genes <- rbind(down_exp, up_exp)



# find expression values from DE genes
DE_genes_expression <- data_exp[rownames(data_exp) %in% rownames(DE_genes),]

DE_genes_expression <- log2(DE_genes_expression)
is.na(DE_genes_expression) <- sapply(DE_genes_expression, is.infinite)



png(filename = "dataChunqin/results/heatmap/heatmap.png", height = 800, width = (22*nrow(DE_genes_expression)))
Heatmap(t(DE_genes_expression), name = "log2(RNA expression)", na_col = "black",
        bottom_annotation_height = unit(3, "cm"),  col = colorRamp2(c(0, max(DE_genes_expression, na.rm = TRUE)), c("blue", "yellow")), show_row_dend = FALSE, show_column_dend = FALSE)
dev.off()


sessionInfo()
