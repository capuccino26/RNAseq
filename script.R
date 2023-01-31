#install.packages(c(“DESeq2”, “Biobase”, “gplots”, “ggplot2”, “RColorBrewer”, “pheatmap”, “matrixStats”, “lattice”, “gtools”, “dplyr”, “BiocManager”)
#BiocManager::install("tximportData")
#BiocManager::install("tximport")
#BiocManager::install("DESeq2")
#BiocManager::install("vsn")
library("DESeq2")
library("ggplot2")
library("pheatmap")

#Settings for DESeq2
sampleinfo <- read.table(file.path("/starAligned/SampleInfo.txt"), header=TRUE)
coldata <- read.csv2("/starAligned/counts.csv")
coldatab <- subset(coldata,select=-c(EntrezGeneID,Length))

#DESeq2 Execution
dds <- DESeqDataSetFromMatrix(countData = coldatab,colData = sampleinfo,design = ~ Status)
dds <- DESeq(dds)
res <- results(dds)

#Exporting Results
write.table(as.data.frame(res), file="/starAligned/treated/RES_deseq.csv")
#print(res)

#Graphs
##Create file for data
pdf("/starAligned/treated/img.pdf")
###This data must end with "dev.off()"

#general
plotMA(res, ylim=c(-2,2))
plotCounts(dds, gene=which.min(res$padj), intgroup="Status")
##ggplot2
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="Status",returnData=TRUE)
ggplot(d, aes(x=Status, y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + scale_y_log10(breaks=c(25,100,400))
##Heatmap obsolete
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Status","CellType")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df)
##Pheatmap corrected
data <- read.csv2("/GSE60450_LactationGenewiseCounts.csv",sep=';',header=TRUE,row.names="EntrezGeneID")
data <- subset(data,select=-c(Length))
data_subset <- as.matrix(data[rowSums(data)>1000000,])
##DENDROGRAM
###install.packages("dendextend")
###library(dendextend)
my_hclust_gene <- hclust(dist(data_subset), method = "complete")
as.dendrogram(my_hclust_gene) %>% plot(horiz = TRUE)
##Outliers counting
W <- res$stat
maxCooks <- apply(assays(dds)[["cooks"]],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic",ylab="maximum Cook's distance per gene",ylim=c(0,5), cex=.4, col=heat.colors(10))
m <- ncol(dds)
p <- 3
abline(h=qf(.99, p, m - p))
##PCA
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("Status", "CellType"))
##Close file for data
dev.off()
