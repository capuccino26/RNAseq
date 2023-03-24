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

#Merge table res/entrezID (Obsoleto, tabelas com tamanhos diferentes)
##resb<-results(dds,tidy=TRUE)
##dataid<-subset(coldata,select="EntrezGeneID") #For lactation Data
##names(dataid)<-("row")
##dataid<-merge(resb,dataid)
##colnames(dataid)[1]<-"EntrezGeneID"

#Merge table res/EntrezID
ids<-subset(coldata,select="EntrezGeneID")
data<-as.data.frame(res)
dataid<-cbind(ids,data)

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
##Heatmap DESeq2
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[,c("Status","CellType")])
ntd <- normTransform(dds)
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
library("vsn")
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("Status", "CellType"))
##vioplot
library("vioplot")
vioplot(coldata,col=c(rep("orange",3),rep("lightblue",3)),ylab="reads",las=2)
##dispersion
plotDispEsts(dds,ylim=c(1e-6,1e2))
##results histograms
hist(res$pvalue,breaks=20,col="grey")
hist(res$padj,breaks=20,col="grey")
##Shrinking
###Check the coeficient number using "resultsNames(dds)"
BiocManager::install("apeglm")
shrink_res <- lfcShrink(dds,res=res,coef=3)

##Volcano
###Convert to dataframe
data<-as.data.frame(res)
###Plot inicial
ggplot(data=data,aes(x=log2FoldChange,y=pvalue))+geom_point()
###Convert p to -10log
p<-ggplot(data=data, aes(x=log2FoldChange, y=-log10(pvalue))) + geom_point()
###Extra formating
####Theme Minimal
p2<-p+theme_minimal()
####Vertical Lines for log2FoldChange threshold
p2<-p2+geom_vline(xintercept=c(-0.6, 0.6), col="red")
####Horizontal line for p-values
p2<-p2+geom_hline(yintercept=-log10(0.05), col="red")
####The significantly differentially expressed genes are the ones found in the upper-left and upper-right corners.
####Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)
#####Add a column of NAs
data$diffexpressed <- "NO"
#####if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
data$diffexpressed[data$log2FoldChange > 0.6 & data$pvalue < 0.05] <- "UP"
#####if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
data$diffexpressed[data$log2FoldChange < -0.6 & data$pvalue < 0.05] <- "DOWN"
#####Re-plot but this time color the points with "diffexpressed"
p <- ggplot(data=data, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed)) + geom_point() + theme_minimal()
p2 <- p+geom_vline(xintercept=c(-0.6, 0.6), col="red")+geom_hline(yintercept=-log10(0.05), col="red")
####Change colors as desired
p3 <- p2 + scale_color_manual(values=c("blue", "black", "red"))
####Create a variable for colors:
mycolors<-c("blue", "black", "red")
names(mycolors)<-c("DOWN","UP","NO")
p3 <- p2 + scale_color_manual(values=mycolors)
####Naming the diff expressed genes:
#####Need the EntrezGeneID from the raw data:
dataraw<-read.csv2("/GSE60450_LactationGenewiseCounts.csv")
data$datalabel <- NA
data$datalabel[data$diffexpressed != "NO"] <- dataraw$EntrezGeneID[data$diffexpressed != "NO"]
ggplot(data=data, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=datalabel)) + geom_point() + theme_minimal() + geom_text()
####Organize the labels using the "ggrepel" package and the geom_text_repel() function
library(ggrepel)
####Plot adding up all layers we have seen so far
ggplot(data=data, aes(x=log2FoldChange, y=-log10(pvalue), col=diffexpressed, label=datalabel)) + geom_point() + theme_minimal() + geom_text_repel() + scale_color_manual(values=c("blue", "black", "red")) + geom_vline(xintercept=c(-0.6, 0.6), col="red") + geom_hline(yintercept=-log10(0.05), col="red")

#Simple comparative heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrixB <- as.matrix(sampleDists)
pheatmap(sampleDistMatrixB, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists, col = colors)

#Comparative heatmap
library("genefilter")
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
anno <- as.data.frame(colData(vsd)[, c("CellType","Status")])
pheatmap(mat, annotation_col = anno)

##Close file for data
dev.off()
