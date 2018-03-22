rm(list=ls())
library(DESeq2)
library(pathview)
library(clusterProfiler)
library(biomaRt)

# load gene universe
sampleTable_exp <- data.frame(sampleName = c("sg26_4d_ctrl_R1", "sg26_4d_ctrl_R2", "sg26_4d_ctrl_R3", "sg26_4d_dox_R1", "sg26_4d_dox_R2", "sg26_4d_dox_R3"), fileName = c("Galaxy92-[htseq-count_on_sg26_4d_ctrl_R1].tabular", "Galaxy94-[htseq-count_on_sg26_4d_ctrl_R2].tabular", "Galaxy96-[htseq-count_on_sg26_4d_ctrl_R3].tabular", "Galaxy98-[htseq-count_on_sg26_4d_dox_R1].tabular", "Galaxy100-[htseq-count_on_sg26_4d_dox_R2].tabular", "Galaxy102-[htseq-count_on_sg26_4d_dox_R3].tabular"), condition = c("ctrl", "ctrl", "ctrl", "trt", "trt", "trt"))
data_exp <- DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_exp,
                                       directory = "dataChunqin/",
                                       design= ~ condition)
data_exp <- as.matrix(counts(data_exp))
geneUniverse <- rownames(data_exp)


# load DE gene list
DE_genes_down <- read.csv("dataChunqin/results/DE_gene_list/down_exp_p_0d05_logFC_1.csv", row.names = 1, stringsAsFactors = FALSE)
DE_genes_up <- read.csv("dataChunqin/results/DE_gene_list/up_exp_p_0d05_logFC_1.csv", row.names = 1, stringsAsFactors = FALSE)
DE_genes <- rbind(DE_genes_down, DE_genes_up)
head(DE_genes)

# find entrez gene id for each gene
ensembl <- useDataset("hsapiens_gene_ensembl",mart = useMart("ensembl"))


DE_gene_entrez_id <- unique(getBM(attributes = c("hgnc_symbol", "entrezgene"),    
                    filters = "hgnc_symbol",
                    values = rownames(DE_genes),
                    mart = ensembl) )
head(DE_gene_entrez_id)

geneUniverse_entrez_id <- unique(getBM(attributes = c("hgnc_symbol", "entrezgene"),    
                                  filters = "hgnc_symbol",
                                  values = geneUniverse,
                                  mart = ensembl) )
head(geneUniverse_entrez_id)

pathwayEnrichment <- enrichKEGG(as.character(DE_gene_entrez_id$entrezgene), organism = "hsa", keyType = "kegg", pvalueCutoff = 0.05,
           pAdjustMethod = "BH", universe = as.character(geneUniverse_entrez_id$entrezgene), minGSSize = 10, maxGSSize = 500,
           qvalueCutoff = 0.2, use_internal_data = FALSE)

# pathway enrichment summary
pathwayEnrichment$ID
pathwayEnrichment$p.adjust


#### create pathway plots

# convert DE genes to right input format: numeric vector containing logFC with entrezID as names
DE_genes$hgnc_symbol <- toupper(rownames(DE_genes)) # create new column with upper case gene names so that they match biomarts hgnc_symbols
head(DE_genes, n = 100)
DE_genes$entrez_id <- DE_gene_entrez_id[match(DE_genes$hgnc_symbol, DE_gene_entrez_id$hgnc_symbol), "entrezgene"]
pathway_genes <- DE_genes$log2FoldChange
names(pathway_genes) <- DE_genes$entrez_id



pv.out <- pathview(gene.data = pathway_genes, pathway.id = "00100", species = "hsa", out.suffix = "chunqin", kegg.dir = "dataChunqin/results/pathwayEnrichment/", limit = list(gene=c(-round(max(abs(min(pathway_genes)),max(pathway_genes))),round(max(abs(min(pathway_genes)),max(pathway_genes)))), cpd=1), low = "yellow", high = "red", bins = list(gene=8,cpd=10))
pv.out <- pathview(gene.data = pathway_genes, pathway.id = "00900", species = "hsa", out.suffix = "chunqin", kegg.dir = "dataChunqin/results/pathwayEnrichment/", limit = list(gene=c(-round(max(abs(min(pathway_genes)),max(pathway_genes))),round(max(abs(min(pathway_genes)),max(pathway_genes)))), cpd=1), low = "yellow", high = "red", bins = list(gene=8,cpd=10))


