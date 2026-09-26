# ------------------------------------------------------------------------------
# as_bb_input.R
#
# Turns the pipeline output (bb_sites.tsv) into the data frame that
# code/simulation/R/estimate.R::fit_bb() expects: columns pair, nu_D, k_R, n_R.
#
# All analysis choices live here, as arguments, so they can be swept without
# re-running Nextflow:
#
#   min_donor_af / max_donor_af   which donor iSNVs count as "polymorphic"
#   T                             recipient calling threshold. Recipient counts
#                                 with k/n < T are set to 0 (censored), matching
#                                 the threshold correction in site_loglik(): pass
#                                 the SAME T to fit_bb().
#   min_depth                     required in BOTH donor and recipient; sites
#                                 failing it are dropped, not recorded as zero
#   exclude_pe_ppe                drop sites flagged by the PE/PPE mask
#   require_donor_call            TRUE = sites must be LoFreq-called in the donor
#                                 (the Sobel Leonard definition). Recipient-only
#                                 variants are de novo / not informative here.
#
# Usage:
#   source("code/simulation/R/estimate.R")
#   source("code/bottleneck_prep/R/as_bb_input.R")
#   sites <- read_bb_sites("path/to/results/bb_sites.tsv")
#   dat   <- as_bb_input(sites, T = 0.01)
#   fit   <- fit_bb(dat, T = 0.01)
#   summarise_fit(fit)
# ------------------------------------------------------------------------------

read_bb_sites <- function(path) {
  d <- read.delim(path, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  d$called_in_donor     <- as.logical(d$called_in_donor)
  d$called_in_recipient <- as.logical(d$called_in_recipient)
  d$mask_pe_ppe         <- as.logical(d$mask_pe_ppe)
  d
}

as_bb_input <- function(sites,
                        min_donor_af       = 0.01,
                        max_donor_af       = 0.99,
                        T                  = 0.01,
                        min_depth          = 20,
                        exclude_pe_ppe     = TRUE,
                        require_donor_call = TRUE) {
  d <- sites
  if (require_donor_call) d <- d[d$called_in_donor, ]
  if (exclude_pe_ppe)     d <- d[!d$mask_pe_ppe, ]
  d <- d[d$donor_depth >= min_depth & d$recipient_depth >= min_depth, ]
  d <- d[d$donor_af >= min_donor_af & d$donor_af <= max_donor_af, ]

  k_R <- d$recipient_alt
  k_R[d$recipient_alt / d$recipient_depth < T] <- 0L   # censor below threshold

  data.frame(
    pair = d$pair_id,
    pos  = d$pos,
    alt  = d$alt,
    nu_D = d$donor_af,          # from re-counted reads, not LoFreq's AF
    k_R  = as.integer(k_R),
    n_R  = as.integer(d$recipient_depth),
    stringsAsFactors = FALSE
  )
}

# Sites per pair surviving the filters: check this before fitting anything.
bb_site_counts <- function(dat) {
  as.data.frame(table(pair = dat$pair), responseName = "n_sites")
}
