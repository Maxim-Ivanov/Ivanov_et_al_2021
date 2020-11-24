# This R code is used to recover the strandedness of chrRNA-seq reads from Jia et al., 2020 (PMID 32541953);
# 
# When cDNA is prepared by Takara Bio (Clontech) SMARTer PCR cDNA Synthesis kit, the flanking sequences are identical (mirrored):
# AAGCAGTGGTATCAACGCAGAGTAC---RNA_MOLECULE---AAAAAAAAAAAAAAAGTACTCTGCGTTGATACCACTGCTT
# TTCGTCACCATAGTTGCGTCTCATG------------------TTTTTTTTTTTTTTTCATGAGACGCAACTATGGTGACGAA
# Thus, it is not possible to distinguish between the first strand and the second strand of cDNA;
# The full-length cDNA sequencing results in unstranded data;
# 
# However, Jia et al. used this kit with custom oligonucleotides;
# The SMART CDS Primer II A was replaced by NEB Universal miRNA Cloning Linker (which was ligated to 3' end of RNA molecule)
# First strand cDNA synthesis was primed by an oligo complementary to the miRNA Cloning Linker;
# Since the miRNA Cloning Linker and the SMARTer II A Oligonucleotide have different sequences, the flanking sequences are not mirrored anymore:
# AAGCAGTGGTATCAACGCAGAGTAC---RNA_MOLECULE---AAAAAAAAAAAAAAACTGTAGGCACCATCAATGTACTCTGCGTTGATACCACTGCTT
# TTCGTCACCATAGTTGCGTCTCATG------------------TTTTTTTTTTTTTTTGACATCCGTGGTAGTTACATGAGACGCAACTATGGTGACGAA
# Thus, it is possible to guess the strand orientation of the original RNA molecule from the full-length cDNA-seq read by detecting the flanking adapter sequences;

# (TranscriptomeReconstructoR is assumed to be installed)
library(GenomicAlignments)
library(tidyverse)

skip_duplicated_reads <- function(bam) {
  nms <- names(bam)
  idx <- which(duplicated(nms) | duplicated(nms, fromLast = TRUE))
  return(bam[-idx])
}

process_stranded_Clontech_BAM <- function(bamfile, neb_adapter = DNAString("CTGTAGGCACCATCAAT"), filter_by_term_subaln = TRUE, abs_threshold = 30, rel_threshold = 0.05) {
  message("Loading ", bamfile, "..."); flush.console()
  bam <- readGAlignments(bamfile, use.names = TRUE, param = ScanBamParam(what = "seq")) # SEQ is converted to Watson strand
  # (as a result, both forward and reverse reads from genes encoded on the Watson strand always appear as if they were having the NEB adapter downstream)
  # (similarly, genes encoded on the Crick strand always appear as having the NEB adapter upstream)
  message("\t", length(bam), " input reads;"); flush.console()
  # Skip chimeric alignments:
  bam <- skip_duplicated_reads(bam)
  message("\t", length(bam), " non-chimeric reads;"); flush.console()
  if (isTRUE(filter_by_term_subaln)) {
    # Skip reads with suspiciously short 5'-most or 3'-most exonic subalignments
    bam <- TranscriptomeReconstructoR:::filter_ga_by_terminal_subalignments(bam, abs_threshold = abs_threshold, rel_threshold = rel_threshold)
    message("\t", length(bam), " reads after filtering by length of terminal subalignments;"); flush.console()
  }
  # Detect NEB adapters:
  gr <- pmapToAlignments(granges(bam), bam) %>% unname()
  gr2 <- resize(gr, ifelse(width(gr) >= 80, width(gr) - 40, width(gr)), fix = "center")
  seq <- bam %>% as.data.frame() %>% .$seq %>% as("DNAStringSet")
  up <- subseq(seq, start = 1, end = start(gr2))
  down <- subseq(seq, start = end(gr2), end = width(seq))
  m1 <- elementNROWS(vmatchPattern(reverseComplement(neb_adapter), up, max.mismatch = 4)) == 1
  m2 <- elementNROWS(vmatchPattern(neb_adapter, down, max.mismatch = 4)) == 1
  rev <- m1 & !m2
  fw <- m2 & !m1
  message("\t", sum(fw), " reads on the Watson strand;"); flush.console()
  message("\t", sum(rev), " reads on the Crick strand;"); flush.console()
  # Extract stranded reads:
  bam_rev <- bam[rev]
  strand(bam_rev) <- "-"
  bam_fw <- bam[fw]
  strand(bam_fw) <- "+"
  bam2 <- c(bam_rev, bam_fw) %>% sort()
  message("\tStrand info recovered for ", round(length(bam2) / length(bam) * 100, 1), "% reads;"); flush.console()
  return(bam2)
}


# Restore the strandedness:
bamfile <- "chr_merged_final.bam"
chr_bam <- process_stranded_Clontech_BAM(bamfile)
# Convert GAlignments to GRanges and merge into coverage:
chr_cov <- chr_bam %>% granges() %>% TranscriptomeReconstructoR:::convert_GRanges_to_coverage()
# Export as bedGraph:
TranscriptomeReconstructoR::save_GRanges_as_bedGraph(chr_cov, "chrRNAseq.bedgraph.gz")
# (the bedGraph track is used to get screenshots on Figures 4C and S2)
