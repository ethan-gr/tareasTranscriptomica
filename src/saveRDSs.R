# Preparamos el ambiente
library(Rsubread)

# Definimos los conjuntos de archivos bam
bfiles.RNAseq.bwt <- c("results/alignments/aln-bwt-5xFAD-A.bam",
                       "results/alignments/aln-bwt-5xFAD-B.bam",
                       "results/alignments/aln-bwt-WT-A.bam",
                       "results/alignments/aln-bwt-WT-B.bam"
)
bfiles.RNAseq.ht2 <- c("results/alignments/aln-ht2-5xFAD-A.bam",
                       "results/alignments/aln-ht2-5xFAD-B.bam",
                       "results/alignments/aln-ht2-WT-A.bam",
                       "results/alignments/aln-ht2-WT-B.bam"
)
bfiles.microRNAs <- c("results/alignments/aln-bwt-E12-A.bam",
                      "results/alignments/aln-bwt-E12-B.bam",
                      "results/alignments/aln-bwt-E20-A.bam",
                      "results/alignments/aln-bwt-E20-B.bam"
)


# counts para los RNAseq hechos con Hisat2
fc.RNAseq.ht2 <- featureCounts(
  files=bfiles.RNAseq.ht2,
  annot.ext="data/indexes/mm10-ensembl_99-genes.gtf",
  isGTFAnnotationFile=TRUE, 
  useMetaFeatures=TRUE, 
  largestOverlap=TRUE, 
  isPairedEnd=TRUE, 
  requireBothEndsMapped=TRUE, 
  nthreads=15
)

# counts para los RNAseq hechos con Bowtie
fc.RNAseq.bwt <- featureCounts(
  files=bfiles.RNAseq.bwt,
  annot.ext="data/indexes/mm10-ensembl_99-genes.gtf",
  isGTFAnnotationFile=TRUE, 
  useMetaFeatures=TRUE, 
  largestOverlap=TRUE, 
  isPairedEnd=TRUE, 
  requireBothEndsMapped=TRUE, 
  nthreads=15
)

# counts para los microRNAs
fc.microRNAs <- featureCounts(
  files=bfiles.microRNAs,
  annot.ext="data/indexes/mm10-ensembl_99-genes.gtf",
  isGTFAnnotationFile=TRUE, 
  useMetaFeatures=TRUE, 
  largestOverlap=TRUE,
  requireBothEndsMapped=TRUE, 
  nthreads=15
)

# counts para los RNAseq hechos con Hisat2
saveRDS(fc.RNAseq.ht2, file="results/counts/fc-RNAseq-ht2.Rds")
# counts para los RNAseq hechos con Bowtie
saveRDS(fc.RNAseq.bwt, file="results/counts/fc-RNAseq-bwt.Rds")
# counts para los microRNAs
saveRDS(fc.microRNAs, file="results/counts/fc-microRNAs.Rds")
