#install.packages(c(“DESeq2”, “Biobase”, “gplots”, “ggplot2”, “RColorBrewer”, “pheatmap”, “matrixStats”, “lattice”, “gtools”, “dplyr”, “BiocManager”)
#BiocManager::install("tximportData")
#BiocManager::install("tximport")
#BiocManager::install("DESeq2")
#BiocManager::install("vsn")
library("DESeq2")
library("ggplot2")
library("ggplot")
library("pheatmap")
dir <- "E:/Users/Pedro/Downloads E/DOC/RNASEQ_R_PRACTICE/DADOS_1"
sampleinfo <- read.table(file.path(dir,"SampleInfo.txt"), header=TRUE)
coldata <- read.csv2("E:/Users/Pedro/Downloads E/DOC/RNASEQ_R_PRACTICE/DADOS_1/GSE60450_LactationGenewiseCounts.csv")
coldatab <- subset(coldata,select=-c(EntrezGeneID,Length))
dds <- DESeqDataSetFromMatrix(countData = coldatab,colData = sampleinfo,design = ~ Status)
dds <- DESeq(dds)
res <- results(dds)
#Resultados alternativos
##res <- results(dds, contrast=c("Status","virgin","lactate"))
##res <- results(dds, contrast=c("Status","virgin","pregnant"))
##res <- results(dds, contrast=c("Status","lactate","pregnant"))
#Graficos
##plotMA(res, ylim=c(-2,2))
##plotCounts(dds, gene=which.min(res$padj), intgroup="Status")
#ggplot2
##d <- plotCounts(dds, gene=which.min(res$padj), intgroup="Status",returnData=TRUE)
##ggplot(d, aes(x=Status, y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + scale_y_log10(breaks=c(25,100,400))
#Heatmap
##select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
##df <- as.data.frame(colData(dds)[,c("Status","CellType")])
##pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df)
#Contagem de outliers
##W <- res$stat
##maxCooks <- apply(assays(dds)[["cooks"]],1,max)
##idx <- !is.na(W)
##plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic",ylab="maximum Cook's distance per gene",ylim=c(0,5), cex=.4, col=heat.colors(10))
##m <- ncol(dds)
##p <- 3
##abline(h=qf(.99, p, m - p))
#PCA
##vsd <- vst(dds, blind=FALSE)
##plotPCA(vsd, intgroup=c("Status", "CellType"))