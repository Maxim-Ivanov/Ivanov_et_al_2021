# This code allows to recalculate the stats and reproduce the figures used in the manuscript;

library(tidyverse)
library(readxl)
library(GenomicFeatures)
library(TxDb.Athaliana.BioMart.plantsmart28)
txdb_tair <- TxDb.Athaliana.BioMart.plantsmart28
library(rtracklayer)

source("05-Custom_functions.R")

# Load the new annotation returned by 03-Run_TranscriptomeReconstructoR.R:
final_out <- readRDS("final_out.RDS")
hml_genes_v2 <- final_out$hml_genes_v2
hml_genes_v2_RT <- final_out$hml_genes_v2_RT
hml_tx <- final_out$hml_tx
lncrna <- final_out$lncrna

# Load TAIR10 genes:
genes_tair <- genes(txdb_tair, columns = c("gene_id", "tx_type"))
mcols(genes_tair)$tx_type <- mcols(genes_tair)$tx_type %>% unlist() %>% unname()
ebt_tair <- exonsBy(txdb_tair, by = "tx", use.names = TRUE)

# Load Araport11 genes:
# (downloaded from https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_GFF3_genes_transposons.201606.gff.gz)
txdb_ara <- makeTxDbFromGFF("Araport11_GFF3_genes_transposons.201606.gff.gz")
genes_ara <- genes(txdb_ara, columns = c("gene_id", "tx_type"))
mcols(genes_ara)$tx_type <- mcols(genes_ara)$tx_type %>% unlist() %>% unname()
ebt_ara <- exonsBy(txdb_ara, by = "tx", use.names = TRUE)
# Adjust seqinfo:
seqinfo(genes_ara, new2old = as.integer(c(1:5, 7, 6))) <- seqinfo(genes_tair)
seqinfo(ebt_ara, new2old = as.integer(c(1:5, 7, 6))) <- seqinfo(genes_tair)

# Load lncRNA annotation from Zhao et al. 2018 (PMID 30498193):
# (download from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6265284/bin/41467_2018_7500_MOESM6_ESM.xlsx and save as Tab-delimited TXT)
zhao <- suppressWarnings(read_excel("41467_2018_7500_MOESM6_ESM.xlsx", skip = 1)) %>% dplyr::select(`LncRNA ID`, Coordinates, Strand) %>% dplyr::rename(gene_id = `LncRNA ID`, coord = Coordinates, strand = Strand)
zhao <- zhao[-nrow(zhao), ]
zhao <- filter(zhao, strand %in% c("+", "-"))
chrom <- str_sub(zhao$coord, 4, 4)
coord <- str_sub(zhao$coord, start = 6)
spl <- str_split(coord, fixed(".."))
start <- lapply(spl, `[`, 1) %>% unlist() %>% as.integer()
end <- lapply(spl, function(x) { return(x[length(x)]) }) %>% unlist() %>% as.integer()
zhao <- GRanges(seqnames = chrom, ranges = IRanges(start, end = end), strand = zhao$strand, gene_id = zhao$gene_id)
seqinfo(zhao) <- seqinfo(txdb_tair)

# Load TSS-seq, 3'DRS-seq, plaNET-seq and pNET-seq data:
bg_dir <- "." # change to the folder with bedGraph files produced by 01-Download_and_remap_published_data.sh
tss_data <- read_stranded_bedGraph("tss_merged_fw_rev.bedgraph.gz", dir = bg_dir, seqinfo = seqinfo(txdb_tair))
drs_data <- read_stranded_bedGraph("drs_merged_fw_rev.bedgraph.gz", dir = bg_dir, seqinfo = seqinfo(txdb_tair))
planet_data <- read_stranded_bedGraph("planet_first_base.bedgraph.gz", dir = bg_dir, seqinfo = seqinfo(txdb_tair))
pnet_data <- read_stranded_bedGraph("pnet_whole_insert.bedgraph.gz", dir = bg_dir, seqinfo = seqinfo(txdb_tair))

# Load TIF-seq reads (hen2-2 mutant):
tif_data <- import("tif_hen2_merged_final.bam", format = "BED", seqinfo = seqinfo(txdb_tair))

# Load names of sppRNA-containing genes from Thomas et al. 2020 (PMID 32444691):
# (download Supplementary Data 1 from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7244574/bin/41467_2020_16390_MOESM4_ESM.xlsx)
spprna_names <- read_excel("41467_2020_16390_MOESM4_ESM.xlsx", skip = 1, col_names = FALSE, .name_repair = "minimal")[[1]]


### Analyze length of the readthrough tails --------------------------------------------------

analyze_nascent_flanks_of_genes(hml_genes_v2_RT)


### Compare borders of TAIR10 vs Araport11 genes ---------------------------------------------

t1 <- tibble("gene_id" = mcols(genes_tair)$gene_id, "idx_tair" = 1:length(genes_tair))
t2 <- tibble("gene_id" = mcols(genes_ara)$gene_id, "idx_ara" = 1:length(genes_ara))
tbl <- inner_join(t1, t2, by = "gene_id")
tair_par <- genes_tair[tbl$idx_tair]
ara_par <- genes_ara[tbl$idx_ara]

compare_gene_borders(ara_par, tair_par, title = "Difference of gene borders (Araport11 vs TAIR10)", width = 8, height = 5) # SFig. 1
plot_percent_overlap(ara_par, tair_par, title = "Percent overlap (Araport11 vs TAIR10)", width = 8, height = 5)


### Analyze overlaps between called genes and TAIR/Araport -----------------------------------

tbl1 <- classify_genes_by_overlap_with_ref_genes(hml_genes_v2, genes_tair) %>% as_tibble() %>% dplyr::select(type, Overlap) %>% mutate(ann = "TAIR10")
tbl2 <- classify_genes_by_overlap_with_ref_genes(hml_genes_v2, genes_ara) %>% as_tibble() %>% dplyr::select(type, Overlap) %>% mutate(ann = "Araport11")
tbl <- bind_rows(tbl1, tbl2)
table(tbl)

# Fig. 2A:
tbl$ann <- tbl$ann %>% as_factor() %>% relevel(ref = "TAIR10")
tbl$type <- tbl$type %>% factor(levels = c("HC", "MC", "LC"))
title <- "Overlaps between called and known genes"
p <- ggplot(tbl, aes(x = type, fill = Overlap)) + geom_bar(colour = "black") + xlab(NULL) + ylab("Number of genes") + ggtitle(title) +
  scale_fill_manual(values = c(No = "white", Unique = "grey40", Multiple = "grey80")) + theme_bw() + facet_grid(. ~ ann)
for (ext in c(".png", ".pdf")) {
  ggsave(paste0(title, ext), plot = p, width = 5, height = 5, units = "in")
}


### Find matched pairs of called/known genes ---------------------------------------------------

hml_par1 <- classify_genes_by_overlap_with_ref_genes(hml_genes_v2, genes_tair) %>% `[`(mcols(.)$Overlap == "Unique")
tair_par <- findOverlaps(hml_par1, genes_tair, select = "first") %>% genes_tair[.]

hml_par2 <- classify_genes_by_overlap_with_ref_genes(hml_genes_v2, genes_ara) %>% `[`(mcols(.)$Overlap == "Unique")
ara_par <- findOverlaps(hml_par2, genes_ara, select = "first") %>% genes_ara[.]


### Find matched trios of called/TAIR/Araport genes ------------------------------------------

t1 <- tibble(gene_id = mcols(hml_par1)$name, idx1 = 1:length(hml_par1))
t2 <- tibble(gene_id = mcols(hml_par2)$name, idx2 = 1:length(hml_par2))
tbl <- inner_join(t1, t2, by = "gene_id")
hml_trio <- hml_par1[tbl$idx1]
mcols(hml_trio)$tair_mate <- tair_par[tbl$idx1] %>% granges()
mcols(hml_trio)$ara_mate <- ara_par[tbl$idx2] %>% granges()


### Compare borders of called genes vs TAIR/Araport ------------------------------------------

grp1 <- mcols(hml_par1)$type %>% factor(levels = c("HC", "MC", "LC"))
# Fig. 2C (upper):
compare_gene_borders(hml_par1, tair_par, groups = grp1, title = "Difference of gene borders (called genes vs TAIR10)", width = 6, height = 4)
# Fig. 2B (left):
plot_percent_overlap(hml_par1, tair_par, groups = grp1, title = "Percent overlap (called genes vs TAIR10)", width = 4, height = 6)
compute_stats_on_percent_overlap(hml_par1, tair_par, threshold = 0.9)

grp2 <- mcols(hml_par2)$type %>% factor(levels = c("HC", "MC", "LC"))
# Fig. 2C (lower):
compare_gene_borders(hml_par2, ara_par, groups = grp2, title = "Difference of gene borders (called genes vs Araport11)", width = 6, height = 4)
# Fig. 2B (right):
plot_percent_overlap(hml_par2, ara_par, groups = grp2, title = "Percent overlap (called genes vs Araport11)", width = 4, height = 6)
compute_stats_on_percent_overlap(hml_par2, ara_par, threshold = 0.9)


### Metagene plot of independent TSS and PAS signal tracks ---------------------------------------

# Subset trios to HC genes:
hc_trio <- hml_trio[mcols(hml_trio)$type == "HC"]

# Generate windows around TSS and PAS:
tss_win_hc <- granges(hc_trio) %>% generate_windows_around_tss_or_pas(mode = "start", width = 50)
tss_win_tair <- mcols(hc_trio)$tair_mate %>% generate_windows_around_tss_or_pas(mode = "start", width = 50)
tss_win_ara <- mcols(hc_trio)$ara_mate %>% generate_windows_around_tss_or_pas(mode = "start", width = 50)

pas_win_hc <- granges(hc_trio) %>% generate_windows_around_tss_or_pas(mode = "end", width = 50)
pas_win_tair <- mcols(hc_trio)$tair_mate %>% generate_windows_around_tss_or_pas(mode = "end", width = 50)
pas_win_ara <- mcols(hc_trio)$ara_mate %>% generate_windows_around_tss_or_pas(mode = "end", width = 50)

# Draw metagenes around TSS:
tss_m1 <- metageneMatrix(tss_data, tss_win_hc, scaling = FALSE, skip.zeros = FALSE, skip.outliers = FALSE)
tss_m2 <- metageneMatrix(tss_data, tss_win_tair, scaling = FALSE, skip.zeros = FALSE, skip.outliers = FALSE)
tss_m3 <- metageneMatrix(tss_data, tss_win_ara, scaling = FALSE, skip.zeros = FALSE, skip.outliers = FALSE)
tss_matlist <- list("HC genes" = tss_m1, "TAIR10" = tss_m2, "Araport11" = tss_m3)
drawMetagenePlot(tss_matlist, x.axis = seq(-25, 24), vline = 0, linetype = "dotted", title = "TSS-seq around gene starts", 
                 xlabel = "Windows 50 bp centered at gene starts", ylabel = "TSS-seq signal", width = 5, height = 8, units = "in") # Fig. 2D

# Draw metagenes around PAS:
pas_m1 <- metageneMatrix(drs_data, pas_win_hc, scaling = FALSE, skip.zeros = FALSE, skip.outliers = FALSE)
pas_m2 <- metageneMatrix(drs_data, pas_win_tair, scaling = FALSE, skip.zeros = FALSE, skip.outliers = FALSE)
pas_m3 <- metageneMatrix(drs_data, pas_win_ara, scaling = FALSE, skip.zeros = FALSE, skip.outliers = FALSE)
pas_matlist <- list("HC genes" = pas_m1, "TAIR10" = pas_m2, "Araport11" = pas_m3)
drawMetagenePlot(pas_matlist, x.axis = seq(-25, 24), vline = 0, linetype = "dotted", title = "Helicos 3\' DRS-seq around gene ends", 
                 xlabel = "Windows 50 bp centered at gene ends", ylabel = "3\' DRS-seq signal", width = 5, height = 8, units = "in") # Fig. 2E


##### Analyze internal exons in re-discovered known genes (only matched pairs of genes are considered) ---------------------------------

# Fig. 3B:
out1_m <- find_novel_internal_exons(hml_tx, ebt_tair, hml_par1, tair_par, xlab_nbreaks = 8, title = "Difference of matched exon borders in matched gene pairs (called vs TAIR10)") 

#out2_m <- find_novel_internal_exons(hml_tx, ebt_ara, hml_par2, ara_par, xlab_nbreaks = 8, title = "Difference of matched exon borders in matched gene pairs (called vs Araport11)")

# Fig. 3A:
tbl <- tibble("TAIR10" = lapply(out1_m, length), "Araport11" = lapply(out2_m, length), type = names(out1_m)) %>% 
  pivot_longer(cols = c("TAIR10", "Araport11"), names_to = "ann") %>% arrange(ann)
tbl$ann <- tbl$ann %>% as_factor() %>% relevel(ref = "TAIR10")
tbl$type <- tbl$type %>% factor(levels = c("No overlap", "Exact match", "IR", "Alt donor", "Alt acceptor", "Other"))
title <- "Classification of internal exons in matched gene pairs"
p <- ggplot(tbl, aes(x = ann, y = value, fill = type)) + geom_bar(stat = "identity", width = 0.6, colour = "white") + xlab(NULL) + 
  ylab("Number of exons") + theme_bw() + theme(legend.title = element_blank()) + ggtitle(title)
for (ext in c(".png", ".pdf")) {
  ggsave(paste0(title, ext), plot = p, width = 4, height = 6, units = "in")
}


##### Metagene plot of plaNET-seq signal at donor sites of the novel exons ----------------------------------------------

out_tair <- classify_exons(hml_tx, ebt_tair)
lapply(out_tair, length)

# Fig. 3C (left):
ml1 <- vector("list", length(out_tair))
names(ml1) <- names(out_tair)
for (i in seq_along(out_tair)) {
  exons <- out_tair[[i]]
  windows <- resize(exons, 1, "end") %>% resize(50, "center")
  ml1[[i]] <- metageneMatrix(planet_data, windows, scaling = FALSE, skip.zeros = FALSE)
}
drawMetagenePlot(ml1, title = "plaNET-seq signal over donor sites (TAIR10)", xlabel = "50 bp windows centered on donor sites", vline = 0, x.axis = seq(-24, 25), width = 5, height = 6, units = "in")

out_ara <- classify_exons(hml_tx, ebt_ara)
lapply(out_ara, length)

# Fig. 3C (right):
ml2 <- vector("list", length(out_ara))
names(ml2) <- names(out_ara)
for (i in seq_along(out_ara)) {
  exons <- out_ara[[i]]
  windows <- resize(exons, 1, "end") %>% resize(50, "center")
  ml2[[i]] <- metageneMatrix(planet_data, windows, scaling = FALSE, skip.zeros = FALSE)
}
drawMetagenePlot(ml2, title = "plaNET-seq signal over donor sites (Araport)", xlabel = "50 bp windows centered on donor sites", vline = 0, x.axis = seq(-24, 25), width = 5, height = 6, units = "in")


### Annotate called genes by overlap with TAIR, Araport and Zhao 2018 --------------------------------------------------------------

# Find novel genes which do not overlaps with any known gene:
mcols(hml_genes_v2)$novel <- hml_genes_v2 %outside% genes_tair & hml_genes_v2 %outside% genes_ara & hml_genes_v2 %outside% zhao
hml_novel <- hml_genes_v2[mcols(hml_genes_v2)$novel]
length(hml_novel)
table(mcols(hml_novel)$type)
mcols(lncrna)$novel <- lncrna %outside% genes_tair & lncrna %outside% genes_ara & lncrna %outside% zhao
lncrna_novel <- lncrna[mcols(lncrna)$novel]
length(lncrna_novel)

# How many of them are antisense to known genes?
mcols(hml_novel)$as <- overlapsAny(hml_novel, genes_tair, ignore.strand = TRUE) | overlapsAny(hml_novel, genes_ara, ignore.strand = TRUE)
sum(mcols(hml_novel)$as)
table(mcols(hml_novel)$type[mcols(hml_novel)$as])
mcols(lncrna_novel)$as <- overlapsAny(lncrna_novel, genes_tair, ignore.strand = TRUE) | overlapsAny(lncrna_novel, genes_ara, ignore.strand = TRUE)
sum(mcols(lncrna_novel)$as)

# Antisense to TAIR10 only:
hml_as2 <- overlapsAny(hml_novel, genes_tair, ignore.strand = TRUE)
table(mcols(hml_novel)$type[hml_as2])
lnc_as2 <- overlapsAny(lncrna_novel, genes_tair, ignore.strand = TRUE)
sum(lnc_as2)

# Barplot of percent novel intergenic and antisense genes (Fig. 4A):
hc <- mcols(hml_novel)$type == "HC"
mc <- mcols(hml_novel)$type == "MC"
lc <- mcols(hml_novel)$type == "LC"
hml_as_tair <- overlapsAny(hml_novel, genes_tair, ignore.strand = TRUE)
lnc_as_tair <- overlapsAny(lncrna_novel, genes_tair, ignore.strand = TRUE)
hml_as_ara <- overlapsAny(hml_novel, genes_ara, ignore.strand = TRUE)
lnc_as_ara <- overlapsAny(lncrna_novel, genes_ara, ignore.strand = TRUE)

tbl <- tibble(total = c(sum(hc), sum(mc), sum(lc), length(lncrna_novel)), 
              as_tair = c(sum(hc & hml_as_tair), sum(mc & hml_as_tair), sum(lc & hml_as_tair), sum(lnc_as_tair)),
              as_ara = c(sum(hc & hml_as_ara), sum(mc & hml_as_ara), sum(lc & hml_as_ara), sum(lnc_as_ara)),
              type = c("HC", "MC", "LC", "tr.RNA"))
tbl2 <- tbl %>% mutate(ig_tair = total - as_tair, ig_ara = total - as_ara) %>% dplyr::select(-total) %>% 
  pivot_longer(cols = c(as_tair, as_ara, ig_tair, ig_ara)) %>% separate(name, into = c("as", "ann"), sep = "_")
tbl2$as <- ifelse(tbl2$as == "as", "Antisense", "Intergenic") %>% as_factor() %>% relevel(ref = "Antisense")
tbl2$type <- tbl2$type %>% factor(levels = c("HC", "MC", "LC", "lncRNA"))
tbl2$ann <- ifelse(tbl2$ann == "tair", "TAIR10", "Araport11") %>% as_factor() %>% relevel(ref = "TAIR10")
title <- "Novel genes and lncRNAs vs TAIR10 and Araport11"
p <- ggplot(tbl2, aes(x = type, y = value, fill = as)) + geom_bar(stat = "identity", colour = "white") + facet_wrap(. ~ ann) + 
  theme_bw() + theme(legend.title = element_blank()) + xlab(NULL) + ylab("Number of genes") + ggtitle(title)
for (ext in c(".png", ".pdf")) {
  ggsave(paste0(title, ext), plot = p, width = 5, height = 4, units = "in")
}


### Confirm novel genes by pNET-seq (Fig. 4B) --------------------------------------------------------------------------

metagene_matrix_with_flanks <- function(signal, win, extend_by = 100L, mlen = c(20, 100, 20), skip_outliers = 0.99) {
  target_w <- width(win) + 2 * extend_by
  win_ext <- suppressWarnings(resize(win, target_w, "center") %>% trim())
  good <- width(win_ext) == target_w
  win <- win[good]
  m1 <- metageneMatrix(signal, flank(win, 100), scaling = TRUE, matrix.length = mlen[[1]], skip.zeros = FALSE, skip.outliers = FALSE, skip.top.obs = TRUE)
  m2 <- metageneMatrix(signal, win, scaling = TRUE, matrix.length = mlen[[2]], skip.zeros = FALSE, skip.outliers = FALSE, skip.top.obs = TRUE)
  m3 <- metageneMatrix(signal, flank(win, 100, start = FALSE), scaling = TRUE, matrix.length = mlen[[3]], skip.zeros = FALSE, skip.outliers = FALSE, skip.top.obs = TRUE)
  mat <- cbind(m1, m2, m3)
  if (is.numeric(skip_outliers) && is.atomic(skip_outliers)) {
    gene_cov <- rowSums(mat, na.rm = TRUE)
    q <- quantile(gene_cov, skip_outliers)
    outliers <- gene_cov > q
    mat <- mat[!outliers, ]
    message("Skipped ", sum(outliers), " potential outliers;"); flush.console()
  }
  return(mat)
}

m1 <- metagene_matrix_with_flanks(pnet_cov, hml_novel)
m2 <- metagene_matrix_with_flanks(pnet_cov, lncrna_novel)
ml1 <- list("Novel HC MC LC genes" = m1, "Novel transient RNAs" = m2)
drawMetagenePlot(ml1, title = "Metagene pNET-seq novel genes", x.axis = seq(-19, 120), vline = c(0, 100), xlabel = "Scaled genes with 100 bp flanks", width = 7, height = 4, units = "in")


##### Plot TIF-seq last bases around TSS of sppRNA-containing genes --------------------------------------

# Subset HML trios to sppRNA-containing pairs:
hml_trio_spp <- hml_trio[names(mcols(hml_trio)$tair_mate) %in% spprna_names]

# Make windows around TSS:
w1 <- mcols(hml_trio_spp)$tair_mate %>% resize(1, "start") %>% resize(600, "center")
w2 <- mcols(hml_trio_spp)$ara_mate %>% resize(1, "start") %>% resize(600, "center")
w3 <- hml_trio_spp %>% granges() %>% resize(1, "start") %>% resize(600, "center")
all_w <- list("TAIR10" = w1, "Araport11" = w2, "Called genes" = w3)

# Resize TIF-seq reads to the last base and convert to coverage:
tif_hen2_data <- tif_hen2 %>% resize(1, "end") %>% convert_GRanges_to_coverage()

# Metagene plot of TIF-seq last base signal around TSS of sppRNA-containing genes ():
ml <- vector("list", length(all_w))
names(ml) <- names(all_w)
for (i in seq_along(all_w)) {
  ml[[i]] <- metageneMatrix(tif_hen2_data, all_w[[i]], scaling = TRUE, matrix.length = 120, skip.zeros = FALSE, skip.outliers = FALSE, skip.top.obs = TRUE)
}
drawMetagenePlot(ml, x.axis = seq(-59, 60), vline = 0, filename = "TIFseq hen2 over TSS of sppRNA genes", width = 7, height = 5, units = "in",
                 xlabel = "Windows 600 bp centered at TSS (scaled to 5 bp bins)", ylabel = "Average signal of 3\' bases of TIF-seq reads")
