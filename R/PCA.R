rm(list=ls())
library(DESeq2)
library(ggplot2)

#### load data
# exp. data
sampleTable_exp <- data.frame(sampleName = c("sg26_4d_ctrl_R1", "sg26_4d_ctrl_R2", "sg26_4d_ctrl_R3", "sg26_4d_dox_R1", "sg26_4d_dox_R2", "sg26_4d_dox_R3"), fileName = c("Galaxy92-[htseq-count_on_sg26_4d_ctrl_R1].tabular", "Galaxy94-[htseq-count_on_sg26_4d_ctrl_R2].tabular", "Galaxy96-[htseq-count_on_sg26_4d_ctrl_R3].tabular", "Galaxy98-[htseq-count_on_sg26_4d_dox_R1].tabular", "Galaxy100-[htseq-count_on_sg26_4d_dox_R2].tabular", "Galaxy102-[htseq-count_on_sg26_4d_dox_R3].tabular"), condition = c("ctrl", "ctrl", "ctrl", "trt", "trt", "trt"))
data_exp <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_exp,
                                           directory = "dataChunqin/",
                                           design= ~ condition)
data_exp <- counts(data_exp)
head(data_exp)


# negative control
sampleTable_NC <- data.frame(sampleName = c("sgNCI_4d_ctrl_R1", "sgNCI_4d_ctrl_R2", "sgNCI_4d_ctrl_R3", "sgNCI_4d_dox_R1", "sgNCI_4d_dox_R2", "sgNCI_4d_dox_R3"), fileName = c("Galaxy79-[htseq-count_on_sgNCI_4d_ctrl_R1].tabular", "Galaxy81-[htseq-count_on_sgNCI_4d_ctrl_R2].tabular", "Galaxy83-[htseq-count_on_sgNCI_4d_ctrl_R3].tabular", "Galaxy85-[htseq-count_on_sgNCI_4d_dox_R1].tabular", "Galaxy87-[htseq-count_on_sgNCI_4d_dox_R2].tabular", "Galaxy89-[htseq-count_on_sgNCI_4d_dox_R3].tabular"), condition = c("ctrl", "ctrl", "ctrl", "trt", "trt", "trt"))
data_NC <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_NC,
                                          directory = "dataChunqin/",
                                          design= ~ condition)


data_NC <- counts(data_NC)
head(data_NC)


# make PCA
pca_data_exp <- prcomp(t(data_exp))
pca_data_NC <- prcomp(t(data_NC))


# calculate the variances in percentages covered by components.
pca_data_exp_perc <- round(100*pca_data_exp$sdev^2/sum(pca_data_exp$sdev^2),1)
pca_data_NC_perc <- round(100*pca_data_NC$sdev^2/sum(pca_data_NC$sdev^2),1)


# create a data frame with principal component 1 (PC1), PC2, Conditions and sample names
df_pca_data_exp <- data.frame(PC1 = pca_data_exp$x[,1], PC2 = pca_data_exp$x[,2], sample = colnames(data_exp), condition = rep(c("ctrl","trt"),each=3))
df_pca_data_NC <- data.frame(PC1 = pca_data_NC$x[,1], PC2 = pca_data_NC$x[,2], sample = colnames(data_NC), condition = rep(c("ctrl","trt"),each=3))

# plot exp data
pdf(file = "dataChunqin/results/PCA/PCA_exp.pdf")
ggplot(df_pca_data_exp, aes(PC1,PC2, color = condition))+
    geom_point(size=4)+
    labs(x=paste0("PC1 (",pca_data_exp_perc[1],"%)"), y=paste0("PC2 (",pca_data_exp_perc[2],"%)"))+
    geom_text(aes(label=sample),hjust=0.5, vjust=-1) + ylim(-50000, 75000)
dev.off()

# plot negative control
pdf(file = "dataChunqin/results/PCA/PCA_NC.pdf")
ggplot(df_pca_data_NC, aes(PC1,PC2, color = condition))+
    geom_point(size=4)+
    labs(x=paste0("PC1 (",pca_data_NC_perc[1],"%)"), y=paste0("PC2 (",pca_data_NC_perc[2],"%)"))+
    geom_text(aes(label=sample),hjust=0.5, vjust=-1) + ylim(-50000, 50000)
dev.off()


sessionInfo()
