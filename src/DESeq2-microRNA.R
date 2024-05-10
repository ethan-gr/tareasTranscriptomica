# ==============================================================================
# IMPORTS

suppressMessages((library(DESeq2)))
suppressMessages(library(optparse))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(dplyr))

# ==============================================================================

# Create directory for outputs
name <- "fc-microRNAs"
output_dir <- paste("results/", name, "-DESeq2", sep="")
dir.create(output_dir)

# Select the feature counts object
fc <- readRDS(paste("results/", name, sep=""))

samples <- factor(c("E12","E12","E20","E20"))

# Veifying col data
coldata <- data.frame(samples)
rownames(coldata) <- c("E12-A","E12-B","E20-A","E20-B")

# Design matrix
design <- model.matrix(~0+samples)
colnames(design) <- levels(samples)

# Verify column order
colnames(fc$counts) <- c("E12-A","E12-B","E20-A","E20-B")

# ==============================================================================

# creating DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData=fc$counts, colData=coldata, design=design)

mcols(dds)$basepairs <- fc$annotation$Length

# pre-filtering data
nrow(dds)
keep = rowSums(counts(dds)>1) >= 3
dds <- dds[keep,]

vsd <- vst(dds, blind = FALSE, nsub = 500)

# ==============================================================================

# Differential Expression Analysis 
dds <- DESeq(dds)
res <- results(dds)

## the genes that were significantly differentialy expressed

up = (res$log2FoldChange > 1) & (res$padj < 1e-2)
down = (res$log2FoldChange < -1) & (res$padj < 1e-2)

# ==============================================================================

# volcano plot

png(file=paste(output_dir, "/", name, "-Volcano-plot.png", sep=""))

plot(res$log2FoldChange[!(up | down)], -log10(res$padj[!(up | down)]), pch=19, col="gray", cex=0.4,
     xlab = "log2 Expression fold change", ylab = "-log10 padj", main="Volcano plot", xlim=c(-10, 15), ylim=c(0, 350))

points(res$log2FoldChange[up], -log10(res$padj[up]), pch = 19, col = "red", cex = 0.4)

points(res$log2FoldChange[down], -log10(res$padj[down]), pch=19, col="blue", cex=0.4)

abline(h=10, col="black", lty=3)
abline(v=c(-1,1), col="black", lty=3)

dev.off()
# ==============================================================================

# Heatmap
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 20)
mat  <- assay(vsd)[ topVarGenes, ]
mat  <- mat - rowMeans(mat)
png(file=paste(output_dir, "/", name, "-heatmap.png", sep=""))
pheatmap(mat,cluster_rows=TRUE, show_rownames=FALSE , scale = "column")
