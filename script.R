# Packages
install.packages("DESeq2")
install.packages("ggplots")
install.packages("ggplot2")
install.packages("pheatmap")
install.packages("Biobase")
install.packages("RColorBrewer")
install.packages("matrixStats")
install.packages("lattice")
install.packages("gtools")
install.packages("dplyr")
install.packages("BiocManager")
BiocManager::install("tximportData")
BiocManager::install("tximport")
BiocManager::install("DESeq2")
BiocManager::install("vsn")
library("DESeq2")
library("ggplot2")
library("ggplot")
library("pheatmap")

# Files directory
dir <- "C:/Users/Seq/Data"

# Sampleinfo data file
sampleinfo <- read.table(file.path(dir,"SampleInfo.txt"), header=TRUE)

# Count data file
coldata <- read.csv2("C:/Users/Seq/Data/GSE60450_LactationGenewiseCounts.csv")

# Removal of columns for downstream analysis
coldatab <- subset(coldata,select=-c(EntrezGeneID,Length))

# DESeq2 run
dds <- DESeqDataSetFromMatrix(countData = coldatab,colData = sampleinfo,design = ~ Status)
dds <- DESeq(dds)
res <- results(dds)

# Plots
plotMA(res, ylim=c(-2,2))
plotCounts(dds, gene=which.min(res$padj), intgroup="Status")

# Plots through ggplot
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="Status",returnData=TRUE)
ggplot(d, aes(x=Status, y=count)) + geom_point(position=position_jitter(w=0.1,h=0)) + scale_y_log10(breaks=c(25,100,400))

# Heatmap
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Status","CellType")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,cluster_cols=FALSE, annotation_col=df)

# Outliers counting
W <- res$stat
maxCooks <- apply(assays(dds)[["cooks"]],1,max)
idx <- !is.na(W)
plot(rank(W[idx]), maxCooks[idx], xlab="rank of Wald statistic",ylab="maximum Cook's distance per gene",ylim=c(0,5), cex=.4, col=heat.colors(10))
m <- ncol(dds)
p <- 3
abline(h=qf(.99, p, m - p))

#PCA
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("Status", "CellType"))