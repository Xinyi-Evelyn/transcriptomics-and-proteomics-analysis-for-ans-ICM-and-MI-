#remove all variables and functions from the current environment 
rm(list = ls())
setwd("E:/Rlastest/ICM_MI_transcriptomics analysis")
#load the required packages
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(GEOquery)


######step1 Load data with platform and annotation information####

#GSE113871 <- getGEO("GSE113871",destdir = ".",AnnotGPL = T, getGPL = T)
#GSE48166 <- getGEO("GSE48166",destdir = ".",AnnotGPL = T,getGPL = T)
#Since the original article has been retracted, the downloded data is invalid####
#We can download the data directly form the GEO website


gset_113871 <- read.table("GSE113871_Organoid_Infarction_vs_Control_deseq_124360_1337.txt.gz",head=T,sep="\t",
                          quote = "",fill = T,comment.char = "!")


######Step2 clean data 
# when we use the symbol as the rownames, there will be duplictae row names
# we can use the mean as the duplicate genes expression
est_113871 <- gset_113871[,c(1,5:10)]
#check theh duplicate rownames
est_113871[,1][duplicated(est_113871[,1])]

# take the mean of each group 
est_113871 <- est_113871 %>% group_by(Symbol) %>% summarise(across(everything(),mean))

# rename the rows and columns 
Rownames <- as.vector(as.matrix(est_113871[,1]))

est_113871 <- est_113871[,-1]
rownames(est_113871) <- Rownames

Colnames <- rep(c("Control","MI"),c(3,3))
Colnames[1:3]<- paste(Colnames[1:3],1:3)
Colnames[4:6]<- paste(Colnames[4:6],1:3)
colnames(est_113871) <- Colnames

est_113871 <- apply(est_113871,2,as.numeric)

########step3  test the each group distribution####
library("FactoMineR")
library("factoextra") 

data <- t(est_113871)
data <- as.data.frame(data)
data$group <- rep(c("control","MI"),c(3,3))
#remove the last column -- group information
data.pca <- PCA(data[,-ncol(data)], graph = FALSE)

#Visualize Principal Component Analysis
fviz_pca_ind(data.pca,geom.ind = "point",col.ind = data$group,
             palette = c("#00AFBB", "#E7B800"),addEllipses = TRUE,legend.title = "Groups")

ggsave('MI_PCA.png')

#Visualize heatmap
library(pheatmap)
data <- est_113871
rownames(data) <- Rownames 
top_1000 <- names(tail(sort(apply(data,1,sd)),1000))
dat <- data[top_1000,]

dat <- t(scale(t(dat)))

pheatmap(dat,show_colnames = T ,show_rownames = F)
ac <- data.frame(g=rep(c("control","MI"),c(3,3)))
rownames(ac) <- colnames(dat) 

pheatmap(dat,show_colnames =F,show_rownames = F,
         annotation_col=ac,filename = 'heatmap_top1000_sd.png')

######Step4 Since we take the mean value, here we need take the integer
library(DESeq2)
data <- ceiling(data)
coldata <- data.frame(condition =rep(c("control","MI"),c(3,3)), row.names = colnames(data))
dds <- DESeqDataSetFromMatrix(countData = data,
                              colData = coldata,
                              design= ~ condition)
dds <- DESeq(dds)
resultsNames(dds) # lists the coefficients
res <- results(dds, name="condition_MI_vs_control")
#res <- results(dds, contrast=c("condition", "MI", "control"))

degs_MI = as.data.frame(res)
degs_MI = na.omit(degs_MI)

####Step5 Volcano plot
library(ggpubr)

degs_MI$group <- "STABLE"
degs_MI$group[degs_MI$pvalue < 0.01 & degs_MI$log2FoldChange > 2 ] <- "UP"
degs_MI$group[degs_MI$pvalue < 0.01 & degs_MI$log2FoldChange < -2 ] <- "DOWN"
table(degs_MI$group)
df_gene <- degs_MI[,c(2,5,7)]
df_gene$pvalue <- -log10(df_gene$pvalue)
df_gene$names <- rownames(df_gene)


ggscatter(df_gene, x = "log2FoldChange", y = "pvalue",size=0.1,color = 'group') + theme_test()
ggsave("MI.png")

ggscatter(df_gene, x = "log2FoldChange", y ="pvalue", color = "group",size = 0.1, repel = T,
            label = "names",
            label.select = rownames(df_gene)[df_gene$group != "STABLE"] ,
            font.label = 7,
            palette = c("#00AFBB", "#E7B800", "#FC4E07") )
ggsave('volcano.png')

####heatmap

data <- est_113871
rownames(data) <- Rownames 
df_gene_data <- data[rownames(df_gene[df_gene$group != "STABLE",]),]
df_gene_data <- t(scale(t(df_gene_data)))

pheatmap(df_gene_data,show_colnames = T ,show_rownames = F)
ac <- data.frame(group=rep(c("Control","MI"),c(3,3)))
rownames(ac) <- colnames(df_gene_data) 
pheatmap(df_gene_data,show_colnames = T,show_rownames = T,
         annotation_col=ac,filename = 'heatmap_df_gene.png')
  
  









