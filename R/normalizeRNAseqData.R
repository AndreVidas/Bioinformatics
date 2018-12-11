library(DESeq2)

# load data
data <- read.csv("HT_htseq.csv", header = TRUE, stringsAsFactors = FALSE)
colData <- colnames(data)[2:ncol(data)]
data <- data[,-1] # remove first column corresponding to rownames (included as a column since it contains duplicates not allowed in rownames)


# convert data to DESeq2 format
colData <- data.frame(condition = c("ctrl", "ctrl", "ctrl", "case", "case", "case"))
dds <- DESeqDataSetFromMatrix(countData = data, colData = colData, design= ~ condition)


# normalize data
dds <- estimateSizeFactors(dds)
normalized_data <- as.matrix(counts(dds, normalized=TRUE))

# save data
write.csv(x = normalized_data, file = "HT_htseq_normalized.csv")

