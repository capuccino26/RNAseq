#Open aligned data after SamIndex from Snakefile
aligned <- featureCounts(here::here("/starAligned/Aligned.sortedByCoordPICARD.SAMindex.bai"), isPairedEnd=TRUE)
##Optional: Open directly from STAR
aligned <- featureCounts(here::here("/starAligned/Aligned.sortedByCoord.bam"), isPairedEnd=TRUE)
##Optional: Using known genome
aligned <- -a /rawData/genome.gtf -o /starAligned/Aligned.sortedByCoord.bam

#Save and format Table as CSV
write.table(x=data.frame(aligned$annotation[,c("GeneID","Length")],aligned$counts,stringsAsFactors=False),file="/starAligned/counts.csv",quote=FALSE,row.names=FALSE)
##Open data in R for downstream analysis
coldata <- read.csv2("/starAligned/counts.csv")
rownames(coldata) <- coldata$GeneID
