# This file contains R functions which have to be sourced before executing the 04-Reproduce_figures.R code;

analyze_nascent_flanks_of_genes <- function(genes, mode = "down", min_width = 50) {
  stopifnot(!is.null(mcols(genes)$thick))
  stopifnot(!is.null(mcols(genes)$type))
  stopifnot(mode %in% c("up", "down"))
  thick <- mcols(genes)$thick
  if (mode == "down") {
    message("Downstream intervals of nascent transcription (RT tails):")
    w <- pgap(resize(thick, 0, "end"), resize(genes, 0, "end")) %>% width()
  } else {
    message("Upstream intervals of nascent transcription:")
    w <- pgap(resize(genes, 0, "start"), resize(thick, 0, "start")) %>% width()
  }
  hc <- mcols(genes)$type == "HC"
  mc <- mcols(genes)$type == "MC"
  lc <- mcols(genes)$type == "LC"
  wider <- w >= min_width
  message("- Fraction of intervals longer than ", min_width, " bp:")
  message("\tAll genes: ", sum(wider), " / ", length(genes), " = ", round(mean(wider) * 100, 1), "%")
  message("\tHC genes: ", sum(hc & wider), " / ", sum(hc), " = ", round(mean(wider[hc]) * 100, 1), "%")
  message("\tMC genes: ", sum(mc & wider), " / ", sum(mc), " = ", round(mean(wider[mc]) * 100, 1), "%")
  message("\tLC genes: ", sum(lc & wider), " / ", sum(lc), " = ", round(mean(wider[lc]) * 100, 1), "%")
  message("- Median interval length:")
  message("\tAll genes: ", median(w[wider]), " bp")
  message("\tHC genes: ", median(w[hc & wider]), " bp")
  message("\tMC genes: ", median(w[mc & wider]), " bp")
  message("\tLC genes: ", median(w[lc & wider]), " bp")
}

# -----------------------------------------------------------------------------------------------------------

skip_pairs_with_wrong_strands <- function(g1, g2) {
  stopifnot(length(g1) == length(g2))
  if (any(strand(g1) != strand(g2))) {
    bad <- which(strand(g1) != strand(g2))
    message("Skipped ", length(bad), " gene pairs with incompatible strandedness;")
    g1 <- g1[-bad]
    g2 <- g2[-bad]
  }
  return(list(g1, g2))
}

# -----------------------------------------------------------------------------------------------------------

compare_gene_borders <- function(g1, g2, groups = NULL, type = "freqpoly", title = "Difference of gene borders", 
                                 facet_labels = c("5\' borders", "3\' borders"), ylab = "Number of gene pairs", 
                                 bins = 50, size = 1, xlim = c(-500, 500), width = 7, height = 6, units = "in", xlab_nbreaks = NULL) {
  # borders of g1 are plotted relative to g2
  # Negative values: border of g1 is located upstream from the respective g2 border
  # Positive values: g1 starts or ends downstream from g2
  clean <- skip_pairs_with_wrong_strands(g1, g2)
  g1 <- clean[[1]]; g2 <- clean[[2]]
  t1 <- tibble("dist" = ifelse(strand(g1) == "+", start(g1) - start(g2), end(g2) - end(g1)), "type" = facet_labels[[1]])
  t2 <- tibble("dist" = ifelse(strand(g1) == "+", end(g1) - end(g2), start(g2) - start(g1)), "type" = facet_labels[[2]])
  tbl <- bind_rows(t1, t2)
  tbl$type <- tbl$type %>% as_factor() %>% relevel(ref = facet_labels[[1]])
  p <- ggplot(tbl, aes(x = dist, colour = rep(groups, 2)))
  if (type == "barplot") {
    p <- p + geom_bar(stat = "count", fill = "grey60")
  } else {
    p <- p + geom_freqpoly(bins = bins, size = size)
  }
  p <- p + xlim(xlim[[1]], xlim[[2]]) + ggtitle(title) + geom_vline(xintercept = 0, linetype = "dotted") +
    xlab("Difference of coordinates, bp") + ylab(ylab) + facet_grid(. ~ type) + theme_bw() + theme(legend.title = element_blank())
  if (!is.null(xlab_nbreaks) && is.numeric(xlab_nbreaks) && length(xlab_nbreaks) == 1) {
    p <- p + scale_x_continuous(breaks = scales::breaks_extended(n = xlab_nbreaks))
  }
  for (ext in c(".png", ".pdf")) {
    ggsave(paste0(title, ext), plot = p, width = width, height = height, units = units)
  }
}

# ------------------------------------------------------------------------------------------------------------

plot_percent_overlap <- function(g1, g2, groups = NULL, title = "Percent overlap", bins = 50, size = 1, width = 7, height = 7, units = "in") {
  clean <- skip_pairs_with_wrong_strands(g1, g2)
  g1 <- clean[[1]]; g2 <- clean[[2]]
  over <- width(pintersect(g1, g2)) / width(punion(g1, g2, fill.gap = TRUE))
  p <- ggplot(enframe(over), aes(x = value, colour = groups)) + geom_freqpoly(bins = bins, size = size) + geom_vline(xintercept = 1, linetype = "dotted") +
    xlab("Overlap = Intersection / Union") + ylab("Number of gene pairs") + ggtitle(title) + theme_bw()
  for (ext in c(".pdf", ".png")) {
    ggsave(paste0(title, ext), plot = p, width = width, height = height, units = units)
  }
}

# ------------------------------------------------------------------------------------------------------------

compute_stats_on_percent_overlap <- function(g1, g2, threshold) {
  stopifnot(!is.null(mcols(g1)$type))
  clean <- skip_pairs_with_wrong_strands(g1, g2)
  g1 <- clean[[1]]; g2 <- clean[[2]]
  overlap <- width(pintersect(g1, g2)) / width(punion(g1, g2, fill.gap = TRUE))
  valid <- overlap >= threshold
  hc <- mcols(g1)$type == "HC"
  mc <- mcols(g1)$type == "MC"
  lc <- mcols(g1)$type == "LC"
  message("All genes: ", sum(valid), " / ", length(g1), " = ", round(mean(valid) * 100, 1), "%")
  message("HC genes: ", sum(valid[hc]), " / ", sum(hc), " = ", round(mean(valid[hc]) * 100, 1), "%")
  message("MC genes: ", sum(valid[mc]), " / ", sum(mc), " = ", round(mean(valid[mc]) * 100, 1), "%")
  message("LC genes: ", sum(valid[lc]), " / ", sum(lc), " = ", round(mean(valid[lc]) * 100, 1), "%")
}

# ------------------------------------------------------------------------------------------------------------

classify_genes_by_overlap_with_ref_genes <- function(g, gref) {
  old_names <- names(g)
  names(g) <- 1:length(g)
  # Genes without overlap with any reference gene:
  out1 <- g[g %outside% gref]
  if (length(out1) > 0) {
    mcols(out1)$Overlap <- "No"
  }
  # Uniquely matched genes:
  g_one <- g[countOverlaps(g, gref) == 1]
  hits <- findOverlaps(g_one, gref)
  sh <- subjectHits(hits)
  dupl <- duplicated(sh) | duplicated(sh, fromLast = TRUE)
  hits_uniq <- hits[!dupl]
  idx <- queryHits(hits_uniq)
  out2 <- g_one[idx]
  if (length(out2) > 0) {
    mcols(out2)$Overlap <- "Unique"
  }
  # Other genes (multiple match on either side):
  out3 <- g[countOverlaps(g, gref) > 1]
  if (length(idx) > 0) {
    out3 <- c(out3, g_one[-idx])
  } else {
    out3 <- c(out3, g_one)
  }
  if (length(out3) > 0) {
    mcols(out3)$Overlap <- "Multiple"
  }
  out <- c(out1, out2, out3)
  mcols(out)$Overlap <- mcols(out)$Overlap %>% factor(levels = c("Multiple", "Unique", "No"))
  out <- out[order(as.integer(names(out)))]
  names(out) <- old_names
  return(out)
}

# -----------------------------------------------------------------------------------------------------

generate_windows_around_tss_or_pas <- function(gr, mode, width) {
  stopifnot(mode %in% c("start", "end"))
  stopifnot(is.numeric(width) && is.atomic(width))
  out <- suppressWarnings(gr %>% resize(1, mode) %>% resize(width, "center") %>% GenomicRanges::trim())
  bad <- width(out) < width
  if (any(bad)) {
    message("Skipped ", sum(bad), " trimmed windows;")
    out <- out[!bad]
  }
  return(out)
}

# -----------------------------------------------------------------------------------------------------

mark_internal_exons <- function(grl) {
  fw <- elementNROWS(grl) %>% lapply(function(x) { seq(1, x) }) %>% unlist() %>% unname()
  rev <- elementNROWS(grl) %>% lapply(function(x) { seq(x, 1) }) %>% unlist() %>% unname()
  gr <- unlist(grl, use.names = FALSE)
  mcols(gr)$internal <- fw != 1 & rev != 1
  out <- relist(gr, grl)
  return(out)
}


# -------------------------------------------------------------------------------------------------------

extract_gid_from_txid <- function(gr, mode) {
  stopifnot(!is.null(names(gr)))
  stopifnot(mode %in% c("called", "txdb"))
  if (mode == "called") {
    out <- names(gr) %>% str_split("_") %>% lapply(`[`, 1:3) %>% lapply(str_c, collapse = "_") %>% unlist() %>% unname()
  } else {
    out <- names(gr) %>% str_split("\\.") %>% lapply(`[`, 1) %>% unlist() %>% unname()
  }
  return(out)
}

# ------------------------------------------------------------------------------------------------------

filter_internal_exons_in_chosen_genes <- function(tx, mode, chosen_gid) {
  tx <- mark_internal_exons(tx)
  tx_unl <- unlist(tx)
  mcols(tx_unl)$gene_id <- extract_gid_from_txid(tx_unl, mode = mode)
  names(tx_unl) <- NULL
  out <- tx_unl[mcols(tx_unl)$internal] %>% unique() %>% `[`(mcols(.)$gene_id %in% chosen_gid)
  return(out)
}

# -------------------------------------------------------------------------------------------------------

find_novel_internal_exons <- function(tx1, tx2, g1 = NULL, g2 = NULL, max_diff = 5L, title = "Difference of matched exon borders", 
                                      width = 6, height = 4, xlab_nbreaks = NULL) {
  stopifnot(!is.null(g1) && !is.null(g1))
  # g1 and g2 are expected to be pairs of matched genes:
  stopifnot(length(g1) == length(g2))
  # Retain only internal exons in matched gene pairs:
  message("Number of unique internal exons in matched gene pairs:")
  ex1 <- tx1 %>% filter_internal_exons_in_chosen_genes(mode = "called", chosen_gid = mcols(g1)$name)
  message("\t", length(ex1), " exons in called genes;")
  ex2 <- tx2 %>% filter_internal_exons_in_chosen_genes(mode = "txdb", chosen_gid = mcols(g2)$gene_id)
  message("\t", length(ex2), " exons in known genes;")
  # Classify called exons by overlap with known exons:
  message("Classification of called exons:")
  # 1) Novel exons outside of known exons:
  p1 <- ex1[ex1 %outside% ex2]
  message("\t", length(p1), " are outside of known exons;")
  # 2) Re-discovered exons:
  ex1_over <- ex1[ex1 %over% ex2]
  hits <- findOverlaps(ex1_over, ex2)
  par1 <- ex1_over[queryHits(hits)]
  par2 <- ex2[subjectHits(hits)]
  diff_up <- pgap(resize(par1, 0, "start"), resize(par2, 0, "start")) %>% width()
  diff_down <- pgap(resize(par1, 0, "end"), resize(par2, 0, "end")) %>% width()
  # Consider only valid overlaps:
  valid <- diff_up <= max_diff & diff_down <= max_diff
  hits <- hits[valid]
  par1 <- par1[valid]
  par2 <- par2[valid]
  diff_up <- diff_up[valid]
  diff_down <- diff_down[valid]
  # In case of multiple valid overlaps, choose the best one:
  diff_total <- diff_up + diff_down
  best <- tapply(diff_total, queryHits(hits), find_uniq_min) %>% unlist() %>% unname()
  hits <- hits[best]
  idx <- queryHits(hits)
  p2 <- ex1_over[idx]
  message("\t", length(p2), " match to a known exon (with up to ", max_diff, " bp border offset);")
  par2 <- ex2[subjectHits(hits)]
  # Plot how exactly do they match:
  compare_gene_borders(p2, par2, type = "barplot", title = title, facet_labels = c("Acceptor site", "Donor site"),
                       ylab = "Number of exon pairs", bins = 2 * max_diff, size = 1, xlim = c(-max_diff, max_diff), 
                       width = width, height = height, xlab_nbreaks = xlab_nbreaks)
  # 3) Novel exons with overlap to known exons:
  if (length(idx) > 0) {
    other <- ex1_over[-idx]
  } else {
    other <- ex1_over
  }
  message("\t", length(other), " weakly overlap with known exons:")
  # Classify them by overlap with known donor and acceptor sites:
  accep_ext <- resize(ex2, 1, "start") %>% resize(2 * max_diff, "center")
  donor_ext <- resize(ex2, 1, "end") %>% resize(2 * max_diff, "center")
  p3_start <- resize(other, 1, "start")
  p3_end <- resize(other, 1, "end")
  a <- p3_start %over% accep_ext
  b <- p3_end %over% donor_ext
  message("\t\t", sum(a & b), " IR events;")
  message("\t\t", sum(a & !b), " alternative donor;")
  message("\t\t", sum(!a & b), " alternative acceptor;")
  message("\t\t", sum(!a & !b), " alternative donor and acceptor;")
  out <- c("No_overlap" = p1, "Exact_match" = p2, "IR" = other[a & b], "Alt_donor" = other[a & !b], "Alt_acceptor" = other[!a & b], "Other" = other[!a & !b])
  return(out)
}

# -------------------------------------------------------------------------------------------------------

find_uniq_max <- function(x) {
  out <- rep(FALSE, times = length(x))
  idx <- BiocGenerics::which.max(x)
  out[[idx]] <- TRUE
  return(out)
}

# ------------------------------------------------------------------------------------------------------------------------

find_uniq_min <- function(x) {
  out <- rep(FALSE, times = length(x))
  idx <- BiocGenerics::which.min(x)
  out[[idx]] <- TRUE
  return(out)
}

# -----------------------------------------------------------------------------------------------------------

classify_exons <- function(tx1, tx2, offset = 5) {
  ex1_int <- tx1 %>% mark_internal_exons() %>% unlist() %>% unname() %>% `[`(mcols(.)$internal) %>% unique()
  ex2_int <- tx2 %>% mark_internal_exons() %>% unlist() %>% unname() %>% `[`(mcols(.)$internal) %>% unique()
  ex2_all <- tx2 %>% unlist() %>% unname() %>% unique()
  a1 <- ex1_int %>% resize(1, "start")
  d1 <- ex1_int %>% resize(1, "end")
  a2 <- ex2_int %>% resize(1, "start") %>% resize(offset * 2, "center")
  d2 <- ex2_int %>% resize(1, "end") %>% resize(offset * 2, "center")
  over_acc <- a1 %over% a2
  over_donor <- d1 %over% d2
  out <- list("Rediscovered" = ex1_int[over_acc & over_donor], "Novel_exons" = ex1_int[ex1_int %outside% ex2_all])
  return(out)
}

# ----------------------------------------------------------------------------------------------------------

findOutOfBounds <- function(gr) {
  lookup <- as.list(seqlengths(seqinfo(gr)))
  lengths <- lookup[as.character(seqnames(gr))]
  out <- (start(gr) < 0) | (end(gr) > lengths)
  return(out)
}

calcMatrixLength <- function(intervals, matrix.length) {
  max_width <- max(width(intervals))
  min_width <- min(width(intervals))
  if (is.numeric(matrix.length)) {
    out <- matrix.length
  } else if (max_width == min_width) {
    out <- max_width
  } else if (matrix.length == "max") {
    out <- max_width
  } else if (matrix.length == "min") {
    out <- min_width
  } else if (matrix.length == "mean") {
    out <- round(mean(width(intervals)))
  } else if (matrix.length == "median") {
    out <- round(median(width(intervals)))
  } else {
    message("Not sure how to determine the number of bins. Using median interval length by default..."); flush.console()
    out <- round(median(width(intervals)))
  }
  message("Matrix length = ", out); flush.console()
  return(out)
}

expandOrShrink <- function(x, mlen) {
  if (length(x) == mlen) {
    return(as.numeric(x))
  } else {
    runLength(x) <- runLength(x) * mlen
    return(colMeans(matrix(x, ncol = mlen)))###
  }
}

trimOrFill <- function(x, mlen, anchor, na.as.zeros) {
  if (length(x) == mlen) {
    out <- x
  } else if (length(x) > mlen) {
    if (anchor == "start") {
      out <- x[1:mlen]
    } else if (anchor == "end") {
      out <- x[(length(x)-mlen+1):length(x)]
    }
  } else {
    if (isTRUE(na.as.zeros)) { vals <- 0 } else { vals <- NA }
    to_add <- Rle(values = vals, lengths = mlen - length(x))
    if (anchor == "start") {
      out <- c(x, to_add)
    } else if (anchor == "end") {
      out <- c(to_add, x)
    }
  }
  return(as.numeric(out))
}

calcMatrix <- function(signal, intervals, strand = NULL, max_w, min_w, scaling, mlen, anchor, na.as.zeros) {
  if (!is.null(strand)) {
    intervals <- intervals[strand(intervals) == strand]
    signal <- signal[strand(signal) == strand]
  }
  rlelist <- round(coverage(signal, weight = "score"), 8) # to avoid small negative and positive values instead of zeros
  rlelist_split <- rlelist[intervals]
  rlelist_split <- revElements(rlelist_split, strand(intervals) == "-")
  if ((max_w != min_w) || (max_w != mlen)) {
    if (isTRUE(scaling)) {
      if (length(rlelist_split) > 0) { message("Expanding and shrinking (this can be slow)..."); flush.console() }
      numlist <- lapply(rlelist_split, expandOrShrink, mlen = mlen)
    } else {
      if (length(rlelist_split) > 0) { message("Trimming and filling..."); flush.console() }
      numlist <- lapply(rlelist_split, trimOrFill, mlen = mlen, anchor = anchor, na.as.zeros = na.as.zeros)
    }
  } else {
    numlist <- as(rlelist_split, "NumericList")
  }
  mat <- do.call(rbind, numlist)
  return(mat)
}

average_bin <- function(x, shrink_method) {
  if (all(is.na(x))) {
    return(NA)
  } else if (shrink_method == "mean") {
    return(sum(x, na.rm = TRUE) / length(x))
  } else if (shrink_method == "median") {
    return(median(x, na.rm = TRUE))
  } else if (shrink_method == "mode") {
    freq <- table(x)
    return(as.numeric(names(freq)[which.max(freq)]))# the most frequently observed value; can be very slow
    #return(as.numeric(names(sort(-table(x)))[1]))
  }
}

shrink_row <- function(x, mask, shrink_method) {
  return(as.numeric(tapply(x, mask, average_bin, shrink_method = shrink_method, simplify = FALSE)))
}

shrink_to_bins <- function(mat, binsize, shrink_method) {
  mask <- rep(seq(1, ncol(mat) / binsize), each = binsize)
  out <- t(apply(mat, 1, shrink_row, mask = mask, shrink_method = shrink_method))
  return(out)
}

metageneMatrix <- function(signal, intervals, scaling = FALSE, matrix.length = NA, anchor = "start", na.as.zeros = FALSE, 
                           skip.zeros = TRUE, skip.outliers = 0.995, skip.top.obs = FALSE, n.top.obs = 3, 
                           equal.weights = FALSE, antisenseMode = FALSE, shrink = FALSE, binsize = 5, shrink_method = "mean") {
  library(GenomicRanges)
  out_of_bound <- findOutOfBounds(intervals)
  if (sum(out_of_bound) > 0) {
    message(sum(out_of_bound), " intervals were out-of-bounds;"); flush.console()
    intervals <- intervals[!out_of_bound]
  }
  names(intervals) <- 1:length(intervals) # enumerate intervals
  mlen <- calcMatrixLength(intervals, matrix.length)
  max_w <- max(width(intervals)); min_w <- min(width(intervals))
  if (any(strand(intervals)=="*")) {
    message(round(sum(strand(intervals)=="*") / length(intervals) * 100, 1), "% intervals are unstranded!"); flush.console()
  }
  if (any(strand(signal)=="*")) {
    message(round(sum(strand(signal)=="*") / length(signal) * 100, 1), "% signal values are unstranded!"); flush.console()
  }
  if (isTRUE(antisenseMode)) {
    stranded <- strand(signal) %in% c("+", "-")
    strand(signal)[stranded] <- ifelse(strand(signal)[stranded]=="+", "-", "+")
    message("Antisense mode: signal was flipped to the opposite strand;"); flush.console()
  }
  interval_strands <- list("+", "-", "*")
  signal_strands <- list(c("+", "*"), c("-", "*"), c("+", "-", "*"))
  matlist <- vector("list", 3)
  for (i in 1:3) {
    curr_signal <- signal[strand(signal) %in% signal_strands[[i]]]
    curr_intervals <- intervals[strand(intervals) %in% interval_strands[[i]]]
    if (length(curr_intervals) > 0) {
      curr_mat <- calcMatrix(signal = curr_signal, intervals = curr_intervals, max_w = max_w, min_w = min_w, 
                             scaling = scaling, mlen = mlen, anchor = anchor, na.as.zeros = na.as.zeros)
      rownames(curr_mat) <- names(curr_intervals) # to trace back the original interval (because matrix rows are permuted at this step)
      matlist[[i]] <- curr_mat
    }
  }
  mat <- do.call(rbind, matlist)
  mat <- mat[order(as.numeric(rownames(mat))), ] # restore the original order of intervals
  if (isTRUE(shrink) && is.numeric(binsize) && length(binsize) == 1 && shrink_method %in% c("mean", "median", "mode")) {
    if (ncol(mat) %% binsize == 0) {
      message("Averaging signal within ", binsize, " bp bins;")
      mat <- shrink_to_bins(mat, binsize, shrink_method) # when it is required to average signal within N bp bins but scaling == FALSE
      message("Now matrix length is ", ncol(mat), "!")
    } else {
      message("Check the binsize parameter!")
    }
  }
  gene_cov <- rowSums(mat, na.rm = TRUE)
  if (is.numeric(skip.outliers) & length(skip.outliers) == 1) {
    q <- quantile(gene_cov, skip.outliers)
    outliers <- gene_cov > q
    mat <- mat[!outliers, ]
    message("Skipped ", sum(outliers), " potential outliers;"); flush.console()
    gene_cov <- rowSums(mat, na.rm = TRUE) # recalculate gene_cov
  }
  if (isTRUE(skip.top.obs) && is.numeric(n.top.obs) && length(n.top.obs) == 1) {
    mat_sorted <- apply(mat, 2, sort, decreasing = TRUE, na.last = TRUE) # sort each column independently from other columns
    max_vals <- mat_sorted[n.top.obs, ] # find N'th top value
    for (j in 1:ncol(mat)) {
      to_skip <- mat[, j] >= max_vals[j]
      mat[to_skip, j] <- NA # in each column, change N top observations to NA
    }
    message("Skipped ", n.top.obs, " top observations in each bin;")
    gene_cov <- rowSums(mat, na.rm = TRUE)
  }
  if (isTRUE(skip.zeros)) {
    zeros <- gene_cov == 0
    if (sum(!zeros) < 100) {
      stop("Too little intervals with non-zero signal! Consider changing skip.zeros to FALSE!")
    }
    if (sum(zeros) > 0) {
      mat <- mat[!zeros, ]
      message("Skipped ", sum(zeros), " intervals with zero signal;"); flush.console()
    }
  }
  if (isTRUE(equal.weights)) {
    mat <- t(apply(mat, 1, function(x) { x / sum(x, na.rm=TRUE) }))
    message("All intervals were assigned equal weights;"); flush.console()
  }
  return(mat)
}

# -----------------------------------------------------------------------------------------------------

drawMetagenePlot <- function(mat_list, drawCI = TRUE, x.axis = FALSE, title = "Metagene_plot", filename = NULL, xlabel = "Intervals of interest", 
                             ylabel = "Average signal", vline = FALSE, hline = FALSE, linetype = "solid", width = NA, height = NA, units = "in", plotPDF = TRUE, 
                             alpha = 1, alpha_CI = 0.25, scale.plot = NA, ylim = NA, custom.colors = NA, out_dir = ".") {
  library(ggplot2)
  if (identical(x.axis, FALSE)) {
    x.axis <- seq_len(ncol(mat_list[[1]]))
  }
  stopifnot(length(x.axis) == ncol(mat_list[[1]]))
  long_df <- data.frame()
  for (i in seq(length(mat_list))) {
    curr_mat <- mat_list[[i]]
    curr_name <- names(mat_list)[[i]]
    if (!is.na(scale.plot) && is.numeric(scale.plot) && length(scale.plot)==length(mat_list)) {
      curr_mat <- curr_mat * scale.plot[[i]]
    }
    avg <- apply(curr_mat, 2, mean, na.rm = TRUE)
    message("Max value ", i, " = ", max(avg, na.rm = TRUE)); flush.console()
    df <- data.frame("pos" = x.axis, "avg" = avg, "group" = curr_name)
    if (isTRUE(drawCI)) {
      sem <- apply(curr_mat, 2, sd, na.rm = TRUE) / sqrt(nrow(curr_mat))
      lower <- avg - 1.96 * sem; upper <- avg + 1.96 * sem
      df <- cbind(df, data.frame("lower" = lower, "upper" = upper))
    }
    long_df <- rbind(long_df, df)
  }
  p <- ggplot(long_df, aes(x = pos, y = avg, color = group)) + geom_line(size = 1, alpha = alpha) + ggtitle(title) + xlab(xlabel) + ylab(ylabel) + theme_bw()
  if (isTRUE(drawCI)) {
    p <- p + geom_ribbon(aes(x = pos, ymin = lower, ymax = upper, fill = group), alpha = alpha_CI, linetype = "blank")
  }
  if (is.numeric(vline)) {
    p <- p + geom_vline(xintercept = vline, linetype = linetype)
  }
  if (is.numeric(hline)) {
    p <- p + geom_hline(yintercept = hline, linetype = linetype)
  }
  if (!all(is.na(ylim))) {
    p <- p + ylim(ylim[[1]], ylim[[2]])
  }
  if (!all(is.na(custom.colors))) {
    p <- p + scale_colour_manual(values = custom.colors) + scale_fill_manual(values = custom.colors)
  }
  if (is.null(filename)) {
    filename <- title
  }
  filename <- sub("\n", " ", filename)
  ggsave(filename = file.path(out_dir, paste0(filename, ".png")), plot = p, width = width, height = height, units = units)
  if (isTRUE(plotPDF)) {
    ggsave(filename = file.path(out_dir, paste0(filename, ".pdf")), plot = p, width = width, height = height, units = units)
  }
}

# -----------------------------------------------------------------------------------------------------------

read_stranded_bedGraph <- function(fname, dir = ".", seqinfo = NULL) {
  fpath <- file.path(dir, fname)
  gr <- GenomicRanges::trim(rtracklayer::import(fpath, format = "bedGraph", seqinfo = seqinfo))
  strand(gr) <- ifelse(score(gr) >= 0, "+", "-")
  score(gr) <- abs(score(gr))
  return(sort(gr))
}
