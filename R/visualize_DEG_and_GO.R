rm(list=ls())
library(ggplot2)
library(DESeq2)
library(gplots)



######################################################################################################
# heatmap
######################################################################################################

# load gene counts
data_exp <- read.table("../data/featureCounts/annotated_combined.counts", header = TRUE, stringsAsFactors = FALSE, row.names = 1)

# reorder columns
data_exp <- data_exp[,c(1,5,6,7,8,9,10,11,12,2,3,4,13)] 
colnames(data_exp)

# load data
#trt_vs_ctrl <- read.csv2("../results/DEgenes/codingGenes/postHocFiltering/trt_vs_ctrl.csv", row.names = 1, stringsAsFactors = FALSE)
#relap_vs_trt <- read.csv2("../results/DEgenes/codingGenes/postHocFiltering/relap_vs_trt.csv", row.names = 1, stringsAsFactors = FALSE)
#relap_ctrl <- read.csv2("../results/DEgenes/codingGenes/postHocFiltering/relap_vs_ctrl.csv", row.names = 1, stringsAsFactors = FALSE)



# genes of interest

# Genes related to T cell cytotoxicity/functionality:
t_cell_cytotoxicity <- c("Cd8a","Cd4","Ifng","Gzmb","Pdcd1","Tnf","Prf1","Fasl","Cd69","Cxcr3","Tigit","Lag3","Klrg1","Tbx21","Eomes")

# Genes related to DC activation/priming capacity (SERPINB14 = Chicken gene): 
DC_activation <- c("Ccr7","Cd40","Cd80","Cd86","Il12b","H2-K1","SERPINB14","Tap1","Tap2","Batf3")

# Genes related to immune suppressive mechanisms (Tregs, MDSCs):
immune_suppres_mech <- c("Foxp3","Il10","Ly6c1","Ly6g","Tgfb1","Cd274","Ido1","Arg1","Nos2","Ctla4")

# Genes related to ECM remodeling:
ECM_remodeling <- c("Spp1","Sparc","Vcan","Tnc","Col1a1")

# Chemokines: "Ccl21" & "Il8" do not exist in the data results. Ccl19 has two different ENSMUSG-ID's, so it's removed because of ambiguity
chemokines <- c("Ccl2","Ccl3","Ccl4","Ccr5","Ccl5","Ccl17","Ccr4","Ccr7","Ccl21","Ccl22","Ccl25","Cxcl1","Cxcr1","Cxcr2","Il8","Cxcl9","Cxcl10","Cxcr3","Cxcl11","Cxcl12","Cxcr4","Cxcl16")

# Genes related to NK cells and functionality:
NK_cells <- c("Klrk1","Klra1","Klra2","Klra3","Klra4","Klra5","Klra6","Klra7","Klra8","Klra9","Klra10","Kit")


# new heatmaps (13-08-2019)
#cytotoxic_effector_functions <- c("Gzmb", "Cd8a", "Ifng", "Eomes", "Ccl3", "Cxcr3", "Tbx21", "Ccl25", "Prf1", "Klrk1") 
cytotoxic_effector_functions <- c("Gzmb", "Cd8a", "Ifng", "Eomes", "Ccl3", "Cxcr3", "Tbx21", "Prf1", "Klrk1") 
immune_suppression_and_metastasis <- c("Cd4", "Ccl17", "Ctla4", "Il10", "Nos2", "Cd274", "Arg1", "Foxp3", "Tgfb1", "Mrc1", "Ly6g6c")
DC_activation_and_functions <- c("Cd86", "Cd40", "Batf3", "Cxcl9", "Cxcl10", "Cxcl11", "Ccl4")



dim(data_exp)
head(data_exp)
colnames(data_exp)

# normalize data
sizeFactors <- estimateSizeFactorsForMatrix(data_exp[,1:(ncol(data_exp)-1)], locfunc=median)
colSums(data_exp[,1:(ncol(data_exp)-1)])

for(i in 1:length(sizeFactors)){
    data_exp[,i] <- data_exp[,i] / sizeFactors[i]
}

colSums(data_exp[,1:(ncol(data_exp)-1)])



##### t_cell_cytotoxicity
DE_genes_exp <- data_exp[data_exp$symbol %in% t_cell_cytotoxicity,]
nrow(DE_genes_exp) == length(t_cell_cytotoxicity)

# set gene names as rownames
rownames(DE_genes_exp) <- DE_genes_exp$symbol

# remove gene symbol from data.frame
DE_genes_exp <- DE_genes_exp[,1:(ncol(DE_genes_exp)-1)]


# define colors
color.palette  <- colorRampPalette(c("green", "black", "red"))(n=600)


# plot heatmap
pdf(file = "../results/heatmap/t_cell_cytotoxicity_heatmap.pdf")
heatmap.2(as.matrix(DE_genes_exp[,1:12]), scale = "row", col = color.palette, trace = "none", margins=c(9,10), density.info="none", cexRow = 1.0, cexCol = 0.7, dendrogram = "none", Rowv = TRUE, Colv = FALSE,  main = "t_cell_cytotoxicity")
dev.off()



##### DC_activation
DE_genes_exp <- data_exp[data_exp$symbol %in% DC_activation,]
nrow(DE_genes_exp) == length(DC_activation)

# set gene names as rownames
rownames(DE_genes_exp) <- DE_genes_exp$symbol

# remove gene symbol from data.frame
DE_genes_exp <- DE_genes_exp[,1:(ncol(DE_genes_exp)-1)]


# define colors
color.palette  <- colorRampPalette(c("green", "black", "red"))(n=600)

# plot heatmap
pdf(file = "../results/heatmap/DC_activation_heatmap.pdf")
heatmap.2(as.matrix(DE_genes_exp[,1:12]), scale = "row", col = color.palette, trace = "none", margins=c(9,10), density.info="none", cexRow = 1.0, cexCol = 0.7, dendrogram = "none", Rowv = TRUE, Colv = FALSE, main = "DC_activation")
dev.off()



##### immune_suppres_mech
DE_genes_exp <- data_exp[data_exp$symbol %in% immune_suppres_mech,]
nrow(DE_genes_exp) == length(immune_suppres_mech)

# set gene names as rownames
rownames(DE_genes_exp) <- DE_genes_exp$symbol

# remove gene symbol from data.frame
DE_genes_exp <- DE_genes_exp[,1:(ncol(DE_genes_exp)-1)]


# define colors
color.palette  <- colorRampPalette(c("green", "black", "red"))(n=600)

# plot heatmap
pdf(file = "../results/heatmap/immune_suppres_mech_heatmap.pdf")
heatmap.2(as.matrix(DE_genes_exp), scale = "row", col = color.palette, trace = "none", margins=c(9,10), density.info="none", cexRow = 1.0, cexCol = 0.7, dendrogram = "none", Rowv = TRUE, Colv = FALSE, main = "immune_suppres_mech")
dev.off()



##### ECM_remodeling
DE_genes_exp <- data_exp[data_exp$symbol %in% ECM_remodeling,]
nrow(DE_genes_exp) == length(ECM_remodeling)

# set gene names as rownames
rownames(DE_genes_exp) <- DE_genes_exp$symbol

# remove gene symbol from data.frame
DE_genes_exp <- DE_genes_exp[,1:(ncol(DE_genes_exp)-1)]


# define colors
color.palette  <- colorRampPalette(c("green", "black", "red"))(n=600)

# plot heatmap
pdf(file = "../results/heatmap/ECM_remodeling_heatmap.pdf")
heatmap.2(as.matrix(DE_genes_exp), scale = "row", col = color.palette, trace = "none", margins=c(9,10), density.info="none", cexRow = 1.0, cexCol = 0.7, dendrogram = "none", Rowv = TRUE, Colv = FALSE, main = "ECM_remodeling")
dev.off()



##### Chemokines
DE_genes_exp <- data_exp[data_exp$symbol %in% chemokines,]
nrow(DE_genes_exp) == length(chemokines)

# set gene names as rownames
rownames(DE_genes_exp) <- DE_genes_exp$symbol

# remove gene symbol from data.frame
DE_genes_exp <- DE_genes_exp[,1:(ncol(DE_genes_exp)-1)]


# define colors
color.palette  <- colorRampPalette(c("green", "black", "red"))(n=600)

# plot heatmap
pdf(file = "../results/heatmap/chemokines_heatmap.pdf")
heatmap.2(as.matrix(DE_genes_exp), scale = "row", col = color.palette, trace = "none", margins=c(9,10), density.info="none", cexRow = 1.0, cexCol = 0.7, dendrogram = "none", Rowv = TRUE, Colv = FALSE, main = "chemokines")
dev.off()



##### NK_cells 
DE_genes_exp <- data_exp[data_exp$symbol %in% NK_cells,]
nrow(DE_genes_exp) == length(NK_cells)

# set gene names as rownames
rownames(DE_genes_exp) <- DE_genes_exp$symbol

# remove gene symbol from data.frame
DE_genes_exp <- DE_genes_exp[,1:(ncol(DE_genes_exp)-1)]


# define colors
color.palette  <- colorRampPalette(c("green", "black", "red"))(n=600)

# plot heatmap
pdf(file = "../results/heatmap/NK_cells_heatmap.pdf")
heatmap.2(as.matrix(DE_genes_exp), scale = "row", col = color.palette, trace = "none", margins=c(9,10), density.info="none", cexRow = 1.0, cexCol = 0.7, dendrogram = "none", Rowv = TRUE, Colv = FALSE, main = "NK_cells")
dev.off()


#### new heatmaps 13-08-2019

##### cytotoxic_effector_functions
DE_genes_exp <- data_exp[data_exp$symbol %in% cytotoxic_effector_functions,]
nrow(DE_genes_exp) == length(cytotoxic_effector_functions)

# set gene names as rownames
rownames(DE_genes_exp) <- DE_genes_exp$symbol

# remove gene symbol from data.frame
DE_genes_exp <- DE_genes_exp[,1:(ncol(DE_genes_exp)-1)]


# define colors
#color.palette  <- colorRampPalette(c("green", "black", "red"))(n=600)
color.palette  <- colorRampPalette(c("blue", "white", "red"))(n=600)

# plot heatmap
pdf(file = "../results/heatmap/13082019/colorBlindFriendly_italic/cytotoxic_effector_functions.pdf")
heatmap.2(as.matrix(DE_genes_exp), scale = "row", col = color.palette, trace = "none", margins=c(9,10), density.info="none", cexRow = 1.4, cexCol = 0.7, dendrogram = "none", Rowv = TRUE, Colv = FALSE, main = "cytotoxic_effector_functions", labRow=as.expression(lapply(rownames(DE_genes_exp), function(a) bquote(italic(.(a))))))
dev.off()




##### immune_suppression_and_metastasis
DE_genes_exp <- data_exp[data_exp$symbol %in% immune_suppression_and_metastasis,]
nrow(DE_genes_exp) == length(immune_suppression_and_metastasis)

# set gene names as rownames
rownames(DE_genes_exp) <- DE_genes_exp$symbol

# remove gene symbol from data.frame
DE_genes_exp <- DE_genes_exp[,1:(ncol(DE_genes_exp)-1)]


# define colors
#color.palette  <- colorRampPalette(c("green", "black", "red"))(n=600)
color.palette  <- colorRampPalette(c("blue", "white", "red"))(n=600)

# plot heatmap
pdf(file = "../results/heatmap/13082019/colorBlindFriendly_italic/immune_suppression_and_metastasis.pdf")
heatmap.2(as.matrix(DE_genes_exp), scale = "row", col = color.palette, trace = "none", margins=c(9,10), density.info="none", cexRow = 1.4, cexCol = 0.7, dendrogram = "none", Rowv = TRUE, Colv = FALSE, main = "immune_suppression_and_metastasis", labRow=as.expression(lapply(rownames(DE_genes_exp), function(a) bquote(italic(.(a))))))
dev.off()




##### DC_activation_and_functions
DE_genes_exp <- data_exp[data_exp$symbol %in% DC_activation_and_functions,]
nrow(DE_genes_exp) == length(DC_activation_and_functions)

# set gene names as rownames
rownames(DE_genes_exp) <- DE_genes_exp$symbol

# remove gene symbol from data.frame
DE_genes_exp <- DE_genes_exp[,1:(ncol(DE_genes_exp)-1)]


# define colors
#color.palette  <- colorRampPalette(c("green", "black", "red"))(n=600)
color.palette  <- colorRampPalette(c("blue", "white", "red"))(n=600)

# plot heatmap
pdf(file = "../results/heatmap/13082019/colorBlindFriendly_italic/DC_activation_and_functions.pdf")
heatmap.2(as.matrix(DE_genes_exp), scale = "row", col = color.palette, trace = "none", margins=c(9,10), density.info="none", cexRow = 1.4, cexCol = 0.7, dendrogram = "none", Rowv = TRUE, Colv = FALSE, main = "DC_activation_and_functions", labRow=as.expression(lapply(rownames(DE_genes_exp), function(a) bquote(italic(.(a))))))
dev.off()




######################################################################################################
# GO visualization
######################################################################################################


# load data
#trt_vs_ctrl_up <- read.csv2("../results/GOenrichment/codingGenes/postHocFiltering/trt_vs_ctrl_up.csv", row.names = 1, stringsAsFactors = FALSE)
#relap_vs_trt_up <- read.csv2("../results/GOenrichment/codingGenes/postHocFiltering/relap_vs_trt_up.csv", row.names = 1, stringsAsFactors = FALSE)
relap_vs_trt_down <- read.csv2("../results/GOenrichment/codingGenes/postHocFiltering/relap_vs_trt_down.csv", row.names = 1, stringsAsFactors = FALSE)

# change the go term description to upper case for the first letter
capFirst <- function(s) {
    paste(toupper(substring(s, 1, 1)), substring(s, 2), sep = "")
}

#trt_vs_ctrl_up$Term <- capFirst(trt_vs_ctrl_up$Term)
relap_vs_trt_down$Term <- capFirst(relap_vs_trt_down$Term)

# old
#data_tmp <- trt_vs_ctrl_up[trt_vs_ctrl_up$GO.ID %in% c("GO:0001906", "GO:0035458", "GO:0060333", "GO:0002250", "GO:0001913"),]
#data_tmp <- relap_vs_trt_up[relap_vs_trt_up$GO.ID %in% c("GO:0001915", "GO:0002419", "GO:0002686"),]
#data_tmp <- relap_vs_trt_down[relap_vs_trt_down$GO.ID %in% c("GO:0035456", "GO:0045087", "GO:0034341", "GO:0097529"),]

# new (for full gene universe)
#data_tmp <- trt_vs_ctrl_up[trt_vs_ctrl_up$GO.ID %in% c("GO:0045087", "GO:0001816", "GO:0002252", "GO:0035458", "GO:0034341", "GO:0002250"),]
#data_tmp <- relap_vs_trt_up[relap_vs_trt_up$GO.ID %in% c("GO:0001913", "GO:0002476", "GO:0002428"),]
#data_tmp <- relap_vs_trt_down[relap_vs_trt_down$GO.ID %in% c("GO:0034097", "GO:0002376", "GO:0045087", "GO:0034341", "GO:0030595", "GO:0034340"),]

# new (for full gene universe) 13-08-2019
#data_tmp <- trt_vs_ctrl_up[trt_vs_ctrl_up$GO.ID %in% c("GO:0045087","GO:0001816","GO:0002252","GO:0035458","GO:0034341","GO:0002250"),]
data_tmp <- relap_vs_trt_down[relap_vs_trt_down$GO.ID %in% c("GO:0034097","GO:0045087","GO:0034341","GO:0030595","GO:0001816","GO:0035458"),]








data_tmp[data_tmp$classicFisher == "< 1e-30",c("classicFisher", "fdr")] <- 1e-30
#data <- data.frame(Genes = rep(c("Annotated","Significant"), each=nrow(data_tmp)), GO.Term = rep(data_tmp$Term,2), Gene.Count = c(data_tmp$Annotated, data_tmp$Significant), p.value = as.numeric(data_tmp$classicFisher), fdr = as.numeric(data_tmp$fdr), stringsAsFactors = FALSE)
data <- data.frame(Genes = rep(c("Expected","Significant"), each=nrow(data_tmp)), GO.Term = rep(data_tmp$Term,2), Gene.Count = c(data_tmp$Expected, data_tmp$Significant), p.value = as.numeric(data_tmp$classicFisher), fdr = as.numeric(data_tmp$fdr), stringsAsFactors = FALSE)
data$asterix <- NA


# for p-value
#data[data$p.value < 0.05,"asterix"] <- "*"
#data[data$p.value < 0.01,"asterix"] <- "**"
#data[data$p.value < 0.001,"asterix"] <- "***"
#data[data$p.value < 0.0001,"asterix"] <- "****"
#data[data$p.value < 0.0001,"asterix"] <- "****"

# for fdr
data[data$fdr >= 0.05,"asterix"] <- "n.s."
data[data$fdr < 0.05,"asterix"] <- "*"
data[data$fdr < 0.01,"asterix"] <- "**"
data[data$fdr < 0.001,"asterix"] <- "***"
data[data$fdr < 0.0001,"asterix"] <- "****"
data[data$fdr < 0.0001,"asterix"] <- "****"


# remove half of the asterix rows (set them to "", so that they don't appear twice)
#data[(nrow(data)/2+1):nrow(data),"asterix"] <- ""
data[1:(nrow(data)/2),"asterix"] <- ""


# order categories
#data$GO.Term <- factor(data$GO.Term,levels = c("Response to interferon-gamma", "Innate immune response", "Cytokine production", "Cellular response to interferon-beta", "Immune effector process", "Adaptive immune response")[6:1]) # trt_vs_ctrl_up
data$GO.Term <- factor(data$GO.Term,levels = c("Response to interferon-gamma", "Innate immune response", "Cytokine production", "Cellular response to interferon-beta", "Leukocyte chemotaxis", "Response to cytokine")[6:1]) # relap_vs_trt_down

#### Stacked barplot with multiple groups
#p-values
#pdf(file = "../results/GOenrichment/codingGenes/lfcThreshold/fullGeneUniverse/13082019/plot/trt_vs_ctrl_up_pvalue.pdf", width = 10, height = 5)
#pdf(file = "../results/GOenrichment/codingGenes/lfcThreshold/fullGeneUniverse/13082019/plot/relap_vs_trt_up_pvalue.pdf", width = 10, height = 5)
#pdf(file = "../results/GOenrichment/codingGenes/lfcThreshold/fullGeneUniverse/13082019/plot/relap_vs_trt_down_pvalue.pdf", width = 10, height = 5)

#fdr
#pdf(file = "../results/GOenrichment/codingGenes/postHocFiltering/plot/exp_signi/trt_vs_ctrl_up_fdr.pdf", width = 10, height = 5)
#pdf(file = "../results/GOenrichment/codingGenes/postHocFiltering/plot/relap_vs_trt_up_fdr.pdf", width = 10, height = 5)
pdf(file = "../results/GOenrichment/codingGenes/postHocFiltering/plot/exp_signi/relap_vs_trt_down_fdr.pdf", width = 10, height = 5)
ggplot(data=data, aes(x=GO.Term, y=Gene.Count, fill=Genes, width=.6)) +
    geom_bar(stat="identity", position=position_dodge()) +
    geom_text(aes(label=asterix, angle = 90), size=7, position = position_nudge(y = 10)) +
    coord_flip() +
    theme(axis.title.x = element_text(size = 20), axis.text.x = element_text(size = 20), axis.text.y = element_text(size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black")) + xlab("") + ylab("Gene Count") +
    scale_fill_manual("Genes", values = c("Significant" = "blue", "Expected" = "gray"))
    #theme_bw()
dev.off()


