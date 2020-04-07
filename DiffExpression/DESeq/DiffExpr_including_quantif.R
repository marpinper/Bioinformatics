
##Differential expression analysis with quantification 

##Using the .bam files provided and conditions: treated and untreated

##Importing libraries
x<-c("Rsamtools","readxl","GenomicFeatures","GenomicAlignments","DESeq2","gplots","RColorBrewer","IRanges","GenomicRanges","TxDb.Hsapiens.UCSC.hg38.knownGene")
lapply(x, require, character.only = TRUE)

setwd(" ")
dir<-"/path/to/bam/files/.bam"

sampletable<-as.data.frame(read.table("sample_table.txt", sep='\t',header=TRUE))
sampletable<-sampletable[1:13,1:9]
filenames <- file.path(dir, paste0(sampletable[,2],".bam"))

bamfiles <- BamFileList(filenames, yieldSize=2000000)

##Careful with the genome being used (has to be the same as the one used for mapping)
ebg <- exonsBy(TxDb.Hsapiens.UCSC.hg38.knownGene, by="gene")

##Build summarized  object form the quantified .bam files
se <- summarizeOverlaps(features=ebg, reads=bamfiles,
                        mode="Union",
                        singleEnd=TRUE,
                        ignore.strand=FALSE,
                        fragments=FALSE )


se <- se[ rowSums(assay(se)) >= 5, ]

colSums(assay(se))

colData(se)

rowData(se)

str(metadata(rowData(se)))

colData(se) <- DataFrame(sampletable)


##dds object
dds <- DESeqDataSet(se, design = ~  Batch) ##treated, untreated

dds$Batch <- relevel(dds$Batch, ref = "untreated")


###Plotting the counts

##We can see how genes with low counts (bottom left-hand corner) seem to be excessively 
#variable on the ordinary logarithmic scale, while the rlog transform compresses 
##differences for the low count genes for which the data provide little information about differential expression.

rld <- rlog(dds)

par( mfrow = c( 1, 2 ) )

dds <- estimateSizeFactors(dds)

plot( log2( 1 + counts(dds, normalized=TRUE)[ , 1:2] ),
      col=rgb(0,0,0,.2), pch=16, cex=0.3 )
plot( assay(rld)[ , 1:2],
      col=rgb(0,0,0,.2), pch=16, cex=0.3 )


###Sample distances:
##Heatmap of sample-to-sample distances using the rlog-transformed values.

sampleDists <- dist( t( assay(rld) ) )

sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <- paste( rld$Condition, sep="-" )
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
hc <- hclust(sampleDists)
heatmap.2( sampleDistMatrix, Rowv=as.dendrogram(hc),
           symm=TRUE, trace="none", col=colors,
           margins=c(2,10), labCol=FALSE )


###Principal Component Analysis

plotPCA(rld, intgroup = c("Batch"))



##differential expression analysis pipeline

dds <- DESeq(dds)
               

res <- results(dds, contrast=c("Batch", "treated", "untreated"))


##annotating genes

library("AnnotationDbi")
library("Homo.sapiens")

res$hgnc_symbol <- mapIds(Homo.sapiens,
                          keys=row.names(res),
                          column="SYMBOL",
                          keytype="ENTREZID",
                          multiVals="first")

##sort results
resOrdered <- res[order(res$pvalue),]

##filter significant results
resSig <- subset(resOrdered, padj < 0.1)


##divide by up or down regulation
up<-resSig[  (resSig$log2FoldChange)>0,]
upreg<-up[ order(- up$log2FoldChange),]

down<-resSig[  (resSig$log2FoldChange)<0,]
downreg<-down[ order(- down$log2FoldChange),]

write.csv(as.data.frame(upreg), file="upregulated_significative.csv")
write.csv(as.data.frame(downreg), file="downreg_significative.csv")


###Results 
#In MA each gene is represented with a dot. Genes with an adjusted p value below a threshold (here 0.1, the default) are shown in red.
#Gene clustering: The heatmap becomes more interesting if we do not look at absolute expression strength but rather at the amount by which each gene 
#deviates in a specific #sample from the gene’s average across all samples.

DESeq2::plotMA(res, ylim=c(-5,5))

##Heatmap of the most significant genes

library("pheatmap")
mat <- assay(rld)[ head(order(res$padj),30), ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(rld)[,"Condition"])
pheatmap(mat, annotation_col=df)






