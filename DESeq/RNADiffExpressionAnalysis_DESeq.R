library("DESeq2")
library("RColorBrewer")
library("pheatmap")

dir<-("/path/to/counts/")##*._ReadsPerGene.out.tab if its a STAR mapping output
setwd(dir)

ff <- list.files( path = "/path/to/counts/", 
                  pattern = "_ReadsPerGene.out.tab$", full.names = TRUE )
counts.files <- lapply( ff, read.table, skip = 4 )
counts <- as.data.frame( sapply( counts.files, function(x) x[ , 4 ] ) ) ##we need the 4th column

ff <- gsub( "[.]ReadsPerGene[.]out[.]tab", "", ff )
ff <- gsub( "[.]/counts/", "", ff )
colnames(counts) <- c("sample1a","sample1b","sample1c","sample2a","sample2b","sample2c","sample3a",
                      "sample3b","sample3c","sample4a","sample4b","sample4c","sample5a","sample5b","sample5c") ##5 conditions, 3 biological replicates
row.names(counts) <- counts.files[[1]]$V1
cts<- as.matrix(counts)


sampleFiles <- c("sample1a_ReadsPerGene.out.tab","sample1b_ReadsPerGene.out.tab","sample1c_ReadsPerGene.out.tab",
                 "sample2a_ReadsPerGene.out.tab","sample2b_ReadsPerGene.out.tab","sample2c_ReadsPerGene.out.tab",
                 "sample3a_ReadsPerGene.out.tab","sample3b_ReadsPerGene.out.tab","sample3c_ReadsPerGene.out.tab",
                 "sample4a_ReadsPerGene.out.tab","sample4b_ReadsPerGene.out.tab","sample4c_ReadsPerGene.out.tab",
                 "sample5a_ReadsPerGene.out.tab","sample5b_ReadsPerGene.out.tab","sample5c_ReadsPerGene.out.tab")

sampleTable<- data.frame(row.names=c("sample1a","sample1b","sample1c","sample2a","sample2b","sample2c","sample3a",
                      "sample3b","sample3c","sample4a","sample4b","sample4c","sample5a","sample5b","sample5c"), condition=as.factor(c(rep("sample1",3),
                                        rep("sample2", 3), rep("sample3", 3), rep("sample4", 3),rep("sample5", 3))))


dds_2 <- DESeqDataSetFromMatrix(countData = cts,colData = sampleTable,design = ~ condition)

##Pre-filtering
dds_2 <- dds_2[ rowSums(counts(dds_2)) > 1, ]


dds_2 <- DESeq(dds_2) ##internal normalization in the DESeq pipeline

Comp_1<-results(dds_2, contrast=c("condition","sample1","sample2"))
Comp_2<-results(dds_2, contrast=c("condition","sample4","sample3"))
##... al comparisons wanted

Comp_1_resSig <- Comp_1[which(Comp_1$padj <0.05),]
head(Comp_1_resSig[order(Comp_1_resSig$padj, decreasing = TRUE),])
Comp_2_resSig <- Comp_2[which(Comp_2$padj <0.05),]
head(Comp_2_resSig[order(Comp_2_resSig$padj, decreasing = TRUE),])
##...



nrow(Comp_1_resSig) 
nrow(Comp_2_resSig) 
##...

rld_2 <- rlog(dds_2, blind=F) ##rlog

sampleDists <- dist(t(assay(rld_2)))

sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld_2$condition, rld_2$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
