#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# build_bb_table.R
#
# Combines the per-pair outputs of COUNT_PAIR into:
#   bb_sites.tsv        one row per pair x site x alt allele
#   bb_pair_summary.tsv one row per pair
#
# Reads from the working directory:
#   <pair>.donor_lofreq.tsv, <pair>.recipient_lofreq.tsv, <pair>.counts.tsv
#
# Deliberately applies NO thresholds: every LoFreq-called SNP is kept with its
# raw read counts in both hosts. Thresholds / depth / masks live in
# R/as_bb_input.R. Base R only.
# ------------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop("missing argument ", flag)
  args[i + 1]
}
pairs_path  <- get_arg("--pairs")
mask_path   <- get_arg("--mask")
out_sites   <- get_arg("--out-sites")
out_summary <- get_arg("--out-summary")

lofreq_cols <- c("pair_id", "chrom", "pos", "ref", "alt", "qual", "af", "dp", "sb")

read_tsv_base <- function(path, header = TRUE, col_names = NULL) {
  if (file.info(path)$size == 0) {
    out <- as.data.frame(setNames(replicate(length(col_names), character(0), simplify = FALSE),
                                  col_names), stringsAsFactors = FALSE)
    return(out)
  }
  d <- read.delim(path, header = header, sep = "\t", quote = "",
                  colClasses = "character", comment.char = "")
  if (!is.null(col_names)) names(d) <- col_names
  d
}

read_lofreq <- function(pattern) {
  files <- list.files(".", pattern = pattern, full.names = TRUE)
  do.call(rbind, lapply(files, read_tsv_base, header = FALSE, col_names = lofreq_cols))
}

donor_calls <- read_lofreq("\\.donor_lofreq\\.tsv$")
recip_calls <- read_lofreq("\\.recipient_lofreq\\.tsv$")
counts <- do.call(rbind, lapply(list.files(".", pattern = "\\.counts\\.tsv$", full.names = TRUE),
                                read_tsv_base))

pairs <- read.delim(pairs_path, sep = "\t", colClasses = "character", check.names = FALSE)

# ---- candidate (pair, site, alt) rows: union of donor and recipient calls ----
key <- function(d) paste(d$pair_id, d$chrom, d$pos, d$ref, d$alt, sep = "|")
donor_calls$key <- key(donor_calls)
recip_calls$key <- key(recip_calls)

sites <- unique(rbind(donor_calls[, c("pair_id", "chrom", "pos", "ref", "alt")],
                      recip_calls[, c("pair_id", "chrom", "pos", "ref", "alt")]))
sites$key <- key(sites)
sites$called_in_donor     <- sites$key %in% donor_calls$key
sites$called_in_recipient <- sites$key %in% recip_calls$key

# LoFreq annotations (donor and recipient)
lf <- function(calls, prefix) {
  m <- match(sites$key, calls$key)
  data.frame(af   = as.numeric(calls$af[m]),
             qual = as.numeric(calls$qual[m]),
             sb   = as.numeric(calls$sb[m])) |>
    setNames(paste0(prefix, c("_lofreq_af", "_lofreq_qual", "_lofreq_sb")))
}
sites <- cbind(sites, lf(donor_calls, "donor"), lf(recip_calls, "recipient"))

# ---- attach re-counted alleles from bcftools mpileup ----
counts$ckey <- paste(counts$pair_id, counts$chrom, counts$pos, sep = "|")
m <- match(paste(sites$pair_id, sites$chrom, sites$pos, sep = "|"), counts$ckey)

split_int <- function(x) lapply(strsplit(ifelse(is.na(x) | x == ".", "", x), ","), as.integer)

allele_counts <- function(role) {
  alts <- strsplit(ifelse(is.na(m), "", counts$alts[m]), ",")
  AD   <- split_int(counts[[paste0(role, "_AD")]][m])
  ADF  <- split_int(counts[[paste0(role, "_ADF")]][m])
  ADR  <- split_int(counts[[paste0(role, "_ADR")]][m])
  idx  <- mapply(function(a, alt) match(alt, a), alts, sites$alt)   # position among ALTs
  pick <- function(v, i) if (length(v) == 0 || is.na(i) || i > length(v)) 0L else v[i]
  ref_n   <- vapply(AD, function(v) if (length(v)) v[1] else 0L, integer(1))
  alt_n   <- mapply(function(v, i) pick(v, i + 1L), AD, idx)
  alt_fwd <- mapply(function(v, i) pick(v, i + 1L), ADF, idx)
  alt_rev <- mapply(function(v, i) pick(v, i + 1L), ADR, idx)
  depth   <- vapply(AD, function(v) sum(v, na.rm = TRUE), integer(1))
  out <- data.frame(ref_n, alt_n, depth,
                    af = ifelse(depth > 0, alt_n / depth, NA_real_),
                    alt_fwd, alt_rev)
  setNames(out, paste0(role, c("_ref", "_alt", "_depth", "_af", "_alt_fwd", "_alt_rev")))
}
sites <- cbind(sites, allele_counts("donor"), allele_counts("recipient"))

# ---- PE/PPE (and anything else in the mask BED): flag, do not drop ----
bed <- read.delim(gzfile(mask_path), header = FALSE, colClasses = c("character", "integer", "integer"),
                  col.names = c("chrom", "start", "end"), comment.char = "#")[, 1:3]
pos_i <- as.integer(sites$pos)
sites$mask_pe_ppe <- vapply(seq_len(nrow(sites)), function(i) {
  any(bed$chrom == sites$chrom[i] & bed$start < pos_i[i] & pos_i[i] <= bed$end)
}, logical(1))

# ---- join pair metadata (donor/recipient ids, cohort, anything else supplied) ----
sites$key <- NULL
sites <- merge(pairs, sites, by = "pair_id", all.y = TRUE, sort = FALSE)
sites$pos <- as.integer(sites$pos)
sites <- sites[order(sites$pair_id, sites$chrom, sites$pos, sites$alt), ]

write.table(sites, out_sites, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

# ---- per-pair summary (every pair in the manifest appears, even with 0 sites) ----
summ <- do.call(rbind, lapply(pairs$pair_id, function(p) {
  d <- sites[sites$pair_id == p, ]
  data.frame(
    pair_id                  = p,
    n_sites                  = nrow(d),
    n_called_donor           = sum(d$called_in_donor),
    n_called_recipient       = sum(d$called_in_recipient),
    n_called_both            = sum(d$called_in_donor & d$called_in_recipient),
    n_called_donor_unmasked  = sum(d$called_in_donor & !d$mask_pe_ppe),
    median_donor_depth       = if (nrow(d)) median(d$donor_depth) else NA,
    median_recipient_depth   = if (nrow(d)) median(d$recipient_depth) else NA
  )
}))
write.table(summ, out_summary, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

message(sprintf("bb_sites.tsv: %d rows across %d pairs", nrow(sites), length(unique(sites$pair_id))))
