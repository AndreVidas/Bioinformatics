rm(list=ls())
setwd("/Users/gdn417/Documents/projects/mieDitte/scripts")
library(DESeq2)
library(ggplot2)
library(Rtsne)
library(ggrepel)

#### load data
# exp. data
data_mouse <- read.table("../data/featureCounts/annotated_combined.counts", header = TRUE, stringsAsFactors = FALSE, row.names = 1)
colnames(data_mouse)
head(data_mouse)


# remove symbol column in data
data_mouse <- data_mouse[,1:12]

colnames(data_mouse)


mouse_condition <- c("ctrl","relap","relap","relap","ctrl","ctrl","ctrl", "trt", "trt", "trt", "trt", "relap")


coldata <- data.frame(condition = mouse_condition)
print(coldata)
dds <- DESeqDataSetFromMatrix(countData = data_mouse,
                              colData = coldata,
                              design= ~ condition)


# normalize data_exp
dds <- estimateSizeFactors(dds)
data_exp <- as.matrix(counts(dds, normalized=TRUE))
head(data_exp)


# make PCA
pca_data_exp <- prcomp(t(data_exp), center = TRUE)


# calculate the variances in percentages covered by components.
pca_data_exp_perc <- round(100*pca_data_exp$sdev^2/sum(pca_data_exp$sdev^2),1)


# create a data frame with principal component 1 (PC1), PC2, Conditions and sample names
df_pca_data_exp <- data.frame(PC1 = pca_data_exp$x[,1], PC2 = pca_data_exp$x[,2], sample = colnames(data_exp), condition = coldata)


# plot exp data
pdf(file = "../results/dimReduction/PCA_exp.pdf", width = 12, height = 12)
ggplot(df_pca_data_exp, aes(PC1,PC2, color = condition))+
    geom_point(size=4)+
    labs(x=paste0("PC1 (",pca_data_exp_perc[1],"%)"), y=paste0("PC2 (",pca_data_exp_perc[2],"%)"))+
    geom_text_repel(
        aes(label = sample),
        size = 5, # USER: gene label size
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines"))
dev.off()



######################################################################################################
#### T-SNE
######################################################################################################

#### make t-sne with WTCas9s
data.tsne <- Rtsne(t(data_exp), theta=0.0, pca=TRUE,max_iter=1000, perplexity = 2)
# display object summary
data.tsne
# display embedding coordinates
data.tsne$Y

# create a data frame with principal component 1 (PC1), PC2, Conditions and sample names
mouse_condition <- c("ctrl","relap","relap","relap","ctrl","ctrl","ctrl", "trt", "trt", "trt", "trt", "relap")
mouse_number <- sapply(colnames(data_exp), FUN = function(x) unlist(strsplit(x, split = "Sample"))[2])
#mouse_recipient <- substr(sapply(colnames(data_exp), FUN = function(x) unlist(strsplit(x, split = "\\."))[4]),1,1)
data.tsne <- data.frame(dim1 = data.tsne$Y[,1], dim2 = data.tsne$Y[,2], sample = colnames(data_exp), condition = mouse_condition, mouseNumber = mouse_number)#, recipient = mouse_recipient)



pdf(file = "../results/dimReduction/tsne.pdf", width = 8, height = 8)
ggplot(data.tsne, aes(dim1,dim2, color = condition))+ ggtitle("t-SNE")+
    geom_point(size=4)+
    #labs(x=paste0("PC1 (",pca_data_exp_perc[1],"%)"), y=paste0("PC2 (",pca_data_exp_perc[2],"%)"))+
    #geom_text(aes(label=mouse_number),hjust=0.5, vjust=-1) + ylim(-50000, 75000)
    geom_text_repel(aes(label=paste0(condition,"_",mouse_number))) #+ ylim(-50000, 75000)
dev.off()


sessionInfo()