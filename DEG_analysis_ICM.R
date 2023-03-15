rm(list=ls())
setwd("E:/Rlastest/ICM_MI_transcriptomics analysis")
library(dplyr)
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library("FactoMineR")
library("factoextra") 
library(limma)

#read in count-data
#To get started with this analysis, downloaded the file available online from GEO and extract the relevant files.
ICM_data <- read.table("GSE48166_Cufflinks_FPKM.txt.gz",sep= "\t", fill=T, header = T,quote = '', comment.char="!")
ICM_count_matrix <- ICM_data[,c(2,11:42)]

sum(duplicated(ICM_data[,1]))
ICM_count_matrix[,1][duplicated(ICM_count_matrix[,1])]

ICM_count_matrix <- ICM_count_matrix %>% group_by(geneSymbol) %>% summarise(across(everything(),mean))
ICM_count_matrix <- as.data.frame(ICM_count_matrix)
geneSymbol <- ICM_count_matrix$geneSymbol
ICM_count_matrix <- ICM_count_matrix[,-1]
rownames(ICM_count_matrix) <- geneSymbol 
ICM_count_matrix <- apply(ICM_count_matrix,2,as.numeric)

###
ICM_log <- log2(ICM_count_matrix + 1)
rownames(ICM_log) <- geneSymbol
grouplist <- data.frame(group = rep(c("Ctrl","ICM"),c(16,16)))

exprSet <- ICM_log
rownames(exprSet) <- geneSymbol 
design <- model.matrix(~0+factor(grouplist$group))
colnames(design) <- levels(factor(grouplist$group))
rownames(design) <- colnames(ICM_count_matrix)
contrast.matrix <- makeContrasts("ICM-Ctrl",levels=design)

deg <- function(exprSet, design, contrast.matrix){
  fit <- lmFit(exprSet,design)
  fit <- contrasts.fit(fit, contrast.matrix)
  fit <- eBayes(fit)
  tempOutput <- topTable(fit,coef=1,n=Inf)
  nrDEG <- na.omit(tempOutput)
  return(nrDEG)
}


nrDGE <- deg(exprSet = exprSet,design=design,contrast.matrix = contrast.matrix)

########################
DEG <- nrDGE[,c(1,4)]
DEG$state <- "stable"
DEG$state[DEG$logFC > 1 & DEG$P.Value < 0.01]  <- "Up"
DEG$state[DEG$logFC < -1 & DEG$P.Value < 0.01]  <- "Down"
sum(DEG$state=="Up")
sum(DEG$state=="Down")

########################Volcano plot
library(ggpubr)
DEG$logPvalue <- -log10(DEG$P.Value)
DEG$gene_names <- rownames(DEG)
ggscatter(DEG,x="logFC",y="logPvalue",size = 1,color= "state") + ylab("-log10P.value")

ggscatter(DEG, x="logFC",y="logPvalue", color = "state",size = 1, repel = T,
          label = "gene_names",
          label.select = rownames(DEG)[DEG$state != "stable"] ,
          font.label = 7,
          palette = c("#00AFBB", "#E7B800", "#FC4E07") )
ggsave('ICM_volcano.png')



######################heatmap
library(pheatmap)
data <- ICM_log
de_genes <- na.omit(match(DEG$gene_names[DEG$state != "stable"],rownames(data)))
data <- data[de_genes,]
pheatmap(data,show_colnames = T ,show_rownames = T)
ac <- data.frame(group=rep(c("Control","ICM"),c(16,16)))
rownames(ac) <- colnames(data) 
pheatmap(data,show_colnames = T,show_rownames = T,
         annotation_col=ac,display_numbers = T)



