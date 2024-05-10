# ==============================================================================
# IMPORTS

library(edgeR)
library(pheatmap)

# ==============================================================================

# Create directory for outputs
name <- "fc-microRNAs"
output_dir <- paste("results/", name, "-edgeR", sep="")
dir.create(output_dir)

# Select the feature counts object
fc <- readRDS(paste("results/", name, ".Rds", sep=""))

# Define the samples
samples <- factor(c("E12", "E12", "E20", "E20"))
colnames(fc$counts) <- c("E12A", "E12B", "E20A", "E20B")

# Load data as DGEList
DGEList.fc <- DGEList(counts=fc$counts, group=samples, genes=fc$annotation[,c("GeneID","Length")])

# ==============================================================================

# Delete lowly expressed genes (genes with cpm < 1 in at least 2 replicates)
DGEList.fc <- DGEList.fc[rowSums(cpm(DGEList.fc) > 1 ) >= 2,]

# Make desgin matrix
design <- model.matrix(~0 + samples)
colnames(design) <- levels(samples)

# Calculate norm factors
DGEList.fc <- calcNormFactors(DGEList.fc)
DGEList.fc <- estimateDisp(DGEList.fc, design=design)

# ==============================================================================

# Generate MDS
colors <- c("red", "green", "red", "green")
png(file=paste(output_dir, "/", name, "-MDS.png", sep=""))
plotMDS(DGEList.fc, cex=0.8, col=colors)
dev.off()

# Generate PCA
PCA_log2CPM = prcomp(t(cpm(DGEList.fc, log=TRUE)),center=TRUE, scale=TRUE)

# Genereate scree plot
PCA_perc_var = round ( ((PCA_log2CPM$sdev^2 / sum(PCA_log2CPM$sdev^2))*100), 1)

png(file=paste(output_dir, "/", name, "-Scree.png", sep=""))
barplot(PCA_perc_var, names=colnames(PCA_log2CPM$x), 
        main="Scree plot microRNA hippocampus", xlab="Principal Components", ylab="Percent Variation")

dev.off()

# Compare PC1 & PC2
PCA_var_lab = c(paste("PC1 - ", PCA_perc_var[1], "% var.", sep=""), paste("PC2 - ", PCA_perc_var[2], "% var.", sep=""))

png(file=paste(output_dir, "/", name, "-PCA.png", sep=''))
plot(PCA_log2CPM$x[,"PC1"], PCA_log2CPM$x[,"PC2"], 
     main="PCA of microRNA from hippocampus", xlab=PCA_var_lab[1], ylab=PCA_var_lab[2], col=colors, pch=19)

legend(60,-35, legend = c("5xFAD","WT","5xFAD","WT"), pch=19, col=colors, cex=0.7)

dev.off()

# ==============================================================================

# Make model
fit <- glmFit(DGEList.fc, design)
contrasts <- makeContrasts(E12vsE20=E20-E12, levels=design)
lrt_hip <- glmLRT(fit, contrast=contrasts[, "E12vsE20"])
DEG_lrt <- as.data.frame(topTags(lrt_hip, n=length(DGEList.fc$counts[,1])))

# Get DEGs
up = (DEG_lrt$logFC > 1) & (DEG_lrt$FDR < 1e2)
write.table(DEG_lrt[up,], paste(output_dir,"/", name, "-DEG-Up.tsv", sep=""), sep="\t", quote=FALSE, row.names=FALSE)

down = (DEG_lrt$logFC < -1) & (DEG_lrt$FDR < 1e-2)
write.table(DEG_lrt[down,], paste(output_dir, "/", name, "-DEG-Down.tsv", sep=""), sep="\t", quote=FALSE, row.names=FALSE)

# ==============================================================================

# Volcano plot

png(file=paste(output_dir, "/", name, "-Volcano-plot.png", sep=""))

plot(DEG_lrt$logFC[!(up | down)], -log10(DEG_lrt$FDR[!(up | down)]), pch=19, col="gray", cex=0.4, 
     xlab="log2 Expression fold change", ylab="-log10 FDR", main="Volcano plot E12vsE20", xlim=c(-10,10),ylim=c(0,30))

points(DEG_lrt$logFC[up], -log10(DEG_lrt$FDR[up]), pch=19, col="red", cex=0.4)

points(DEG_lrt$logFC[down], -log10(DEG_lrt$FDR[down]), pch=19, col="blue", cex=0.4)

abline(h=5, col="black", lty=3)
abline(v=c(-1,1), col="black", lty=3)

dev.off()
# ==============================================================================

# Normalize

# Calculate FPKM and FPKM average
log2_fpkm_hip_RNA = rpkm(DGEList.fc, DGEList.fc$genes$Length, log=TRUE)
log2_fpkm_hip_RNA_average = rpkmByGroup(DGEList.fc, log=TRUE)

# Define a function to calculate TPM from FPKM
fpkm2tpm_log2 <- function(fpkm) { fpkm - log2(sum(2^fpkm)) + log2(1e6) }

# Apply function 
log2_tpm_hip_RNA = apply(log2_fpkm_hip_RNA, 2, fpkm2tpm_log2)
log2_tpm_hip_RNA_average = apply(log2_fpkm_hip_RNA_average, 2, fpkm2tpm_log2)

# Save normalized values
write.table(log2_tpm_hip_RNA, paste(output_dir, "/", name, "-TPM-log2.tsv", sep=""), sep="\t", quote=FALSE)
write.table(log2_tpm_hip_RNA_average, paste(output_dir, "/", name, "-TPM-log2-average-egdeR.tsv", sep=""), sep="\t", quote=FALSE)

# Heatmap
diff <- log2_tpm_hip_RNA[(up | down), ]
png(file=paste(output_dir, "/", name, "-heatmap.png", sep=""))
pheatmap(diff, show_rownames=FALSE, cluster_rows=TRUE, cluster_cols=TRUE, scale="column")
dev.off()
