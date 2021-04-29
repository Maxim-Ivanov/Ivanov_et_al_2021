# This is the pipeline for de novo annotation of S.cerevisiae transcriptome (without correction by SacCer3 gene annotation);
# Observe that adjust_exons_of_long_reads() and detect_alignment_errors() calls were skipped, because yeast have only a few intronic genes;

library(TranscriptomeReconstructoR)
library(tidyverse)
library(rtracklayer)
library(devtools)

devtools::source_url("https://github.com/Maxim-Ivanov/Ivanov_et_al_2021/blob/main/07-Custom_functions.R?raw=TRUE")

bam_dir <- "." # change to the folder containing Direct RNA-seq, CAGE-seq, 3'READS and NET-seq BAM files produced by 04-Download_and_remap_published_Saccharomyces_datasets.sh
drs_bamfiles <- list.files(bam_dir, pattern = "^Garalde2018.*final.bam$", full.names = TRUE)
tss_bamfiles <- list.files(bam_dir, pattern = "^Lu2019.*final.bam$", full.names = TRUE)
pas_bamfiles <- list.files(bam_dir, pattern = "^Liu2017.*final.bam$", full.names = TRUE)
nascent_bamfiles <- list.files(bam_dir, pattern = "^Marquardt2014.*final.bam$", full.names = TRUE)

# Load BAM files:
long_reads <- load_BAM_files(drs_bamfiles, mode = "long_read")
tss_data <- load_BAM_files(tss_bamfiles, mode = "tss") 
nascent_data <- load_BAM_files(nascent_bamfiles, mode = "nascent", ngs_mode = "SE", skip_split_aln = TRUE)

# Load 3'READS data and filter for PAS-supporting reads:
pas_data <- lapply(pas_bamfiles, load_3pREADS_BAM_file)
names(pas_data) <- pas_bamfiles %>% basename() %>% str_replace(".bam$", "")
pas_data <- as(pas_data, "CompressedGRangesList")

# Export the loaded data as tracks for genome browser:
write_grl_as_bed12(nanopore, "Yeast_Direct_RNAseq.bed")
tss_data %>% merge_GRanges() %>% save_GRanges_as_bedGraph("Yeast_CAGEseq.bedgraph.gz")
pas_data %>% merge_GRanges() %>% save_GRanges_as_bedGraph("Yeast_3pREADS.bedgraph.gz")
save_GRanges_as_bedGraph(nascent_data, "Yeast_NETseq.bedgraph.gz")

# Call TSS and PAS coordinates:
tss <- call_TCs(tss_data, min_tpm = 0.25)
pas <- call_TCs(pas_data, min_tpm = 0.25)

# Export TSS and PAS as BED files for the genome browser:
export(tss, "TSS.bed", format = "BED")
export(pas, "PAS.bed", format = "BED")

# Extend Direct RNA-seq reads to the nearby TSS and PAS:
long_reads_2 <- extend_long_reads_to_TSS_and_PAS(long_reads, tss, pas, extend_along_guides = FALSE)

# Call transcript and gene models from Direct RNA-seq reads:
out <- call_transcripts_and_genes(long_reads_2)
hml_genes <- out[[1]]
hml_tx <- out[[2]]
fusion_genes <- out[[3]]
fusion_tx <- out[[4]]
reads_free <- out[[5]]

# Call continuous intervals of nascent transcription:
trans <- call_transcribed_intervals(nascent_data, min_signal = 5)
transcribed <- trans[[1]]
gaps <- trans[[2]]

# Add nascent transcription data to the new annotation:
results <- process_nascent_intervals(hml_genes, transcribed, tss, pas, reads_free, gaps)
hml_genes_v2 <- results[[1]]
hml_genes_v2_RT <- results[[2]]
lncrna <- results[[3]]

# Export the new annotation as a set of BED files: 
export(hml_genes_v2, "Called_yeast_genes.bed", format = "BED")
export(hml_genes_v2_RT, "Called_yeast_genes_with_RT_tails.bed", format = "BED")
write_grl_as_bed12(hml_tx, "Called_yeast_transcripts.bed")
write_grl_as_bed12(fusion_tx, "Fusion_yeast_transcripts.bed")
export(lncrna, "Called_yeast_lncRNAs.bed", format = "BED")

# Also save all results as RDS file:
final_out <- list("hml_genes" = hml_genes, "hml_genes_v2" = hml_genes_v2, "hml_genes_v2_RT" = hml_genes_v2_RT, 
                  "hml_tx" = hml_tx, "fusion_genes" = fusion_genes, "fusion_tx" = fusion_tx, "lncrna" = lncrna)
saveRDS(final_out, "final_out_yeast.RDS")
