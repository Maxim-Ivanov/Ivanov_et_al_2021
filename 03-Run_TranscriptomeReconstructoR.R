# This is the pipeline for de novo annotation of A.thaliana transcriptome (without correction by TAIR10 or Araport11)

library(TranscriptomeReconstructoR)
packageVersion("TranscriptomeReconstructoR") # 0.9.0
library(tidyverse)
library(rtracklayer)

bam_dir <- "." # change to the folder containing Direct RNA-seq, CAGE-seq, PAT-seq and plaNET-seq BAM files produced by 01-Download_and_remap_published_data.sh
drs_bamfiles <- list.files(bam_dir, pattern = "^ont.*final.bam$", full.names = TRUE)
tss_bamfiles <- list.files(bam_dir, pattern = "^cage.*final.bam$", full.names = TRUE)
pas_bamfiles <- list.files(bam_dir, pattern = "^pat.*final.bam$", full.names = TRUE)
nascent_bamfiles <- list.files(bam_dir, pattern = "^planet.*final.bam$", full.names = TRUE)

# Load BAM files:
long_reads <- load_BAM_files(drs_bamfiles, mode = "long_read")
tss_data <- load_BAM_files(tss_bamfiles, mode = "tss") 
pas_data <- load_BAM_files(pas_bamfiles, mode = "pas")
nascent_data <- load_BAM_files(nascent_bamfiles, mode = "nascent", ngs_mode = "PE")

# Export the loaded data as tracks for genome browser:
write_grl_as_bed12(nanopore, "Direct_RNAseq.bed")
# (for faster access, execute "bgzip Direct_RNAseq.bed && tabix -p bed Direct_RNAseq.bed.gz" in the system shell)
tss_data %>% merge_GRanges() %>% save_GRanges_as_bedGraph("CAGEseq.bedgraph.gz")
pas_data %>% merge_GRanges() %>% save_GRanges_as_bedGraph("PATseq.bedgraph.gz")
save_GRanges_as_bedGraph(nascent_data, "plaNETseq.bedgraph.gz")
# (if using IGV browser, consider converting the bedGraph files to TDF for faster access)

# Call TSS and PAS coordinates:
tss <- call_TCs(tss_data)
pas <- call_TCs(pas_data)

# Export TSS and PAS as BED files for the genome browser:
export(tss, "TSS.bed", format = "BED")
export(pas, "PAS.bed", format = "BED")

# Extend Direct RNA-seq reads to the nearby TSS and PAS:
long_reads_2 <- extend_long_reads_to_TSS_and_PAS(long_reads, tss, pas)

# Adjust exon borders of Direct RNA-seq reads:
long_reads_3 <- adjust_exons_of_long_reads(long_reads_2)

# Detect alignment errors in Direct RNA-seq reads:
long_reads_4 <- detect_alignment_errors(long_reads_3)

# Call transcript and gene models from Direct RNA-seq reads:
out <- call_transcripts_and_genes(long_reads_4)
hml_genes <- out[[1]]
hml_tx <- out[[2]]
fusion_genes <- out[[3]]
fusion_tx <- out[[4]]
reads_free <- out[[5]]

# Call continuous intervals of nascent transcription:
trans <- call_transcribed_intervals(nascent_data)
transcribed <- trans[[1]]
gaps <- trans[[2]]

# Add nascent transcription data to the new annotation:
results <- process_nascent_intervals(hml_genes, transcribed, tss, pas, reads_free, gaps)
hml_genes_v2 <- results[[1]]
hml_genes_v2_RT <- results[[2]]
lncrna <- results[[3]]

# Export the new annotation as a set of BED files: 
export(hml_genes_v2, "Called_genes.bed", format = "BED")
export(hml_genes_v2_RT, "Called_genes_with_RT_tails.bed", format = "BED")
write_grl_as_bed12(hml_tx, "Called_transcripts.bed")
export(fusion_genes, "Fusion_genes.bed", format = "BED")
write_grl_as_bed12(fusion_tx, "Fusion_transcripts.bed")
export(lncrna, "Called_lncRNAs.bed", format = "BED")

# Also save all results as RDS file:
final_out <- list("hml_genes" = hml_genes, "hml_genes_v2" = hml_genes_v2, "hml_genes_v2_RT" = hml_genes_v2_RT, 
                  "hml_tx" = hml_tx, "fusion_genes" = fusion_genes, "fusion_tx" = fusion_tx, "lncrna" = lncrna)
saveRDS(final_out, "final_out.RDS")
