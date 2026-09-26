# ------------------------------------------------------------------------------
# viralbottleneck.R
#
# Runs ViralBottleneck (Zheng, Johnson & Hughes 2025, Virus Evol 11:veaf071;
# erratum 2026, veag055) on the pipeline output bb_sites.tsv.
#
# Install once:
#   remotes::install_github("BowenArchaman/ViralBottleneck")
#
# ViralBottleneck expects one CSV per sample in the WORKING DIRECTORY, named
# <sample>.csv, with columns  Pos, segment, A, T, C, G, Non/Syn  (integer read
# counts), plus a 2-column Donor/Recipient table. The functions below write
# those files to a temporary folder, run the package there and return the
# result table, so you never have to setwd() yourself.
#
# Two modes:
#   pooled = FALSE  one estimate per pair (the package's native behaviour)
#   pooled = TRUE   one estimate across ALL pairs. The package sums site
#                   log-likelihoods within a pair, so stacking every pair's
#                   sites into one pseudo-pair (site ID = pair_id + position)
#                   gives exactly the pooled likelihood with a shared Nb.
#                   With Mtb's handful of sites per pair this is the estimate
#                   to look at; per-pair values are mainly diagnostics.
#
# Usage:
#   source("code/bottleneck_prep/R/as_bb_input.R")      # read_bb_sites()
#   source("code/bottleneck_prep/R/viralbottleneck.R")
#   sites <- read_bb_sites("code/bottleneck_prep/results/colangeli/bb_sites.tsv")
#   run_viralbottleneck(sites, variant_calling = 0.01, pooled = TRUE)
# ------------------------------------------------------------------------------

# Sites that go to ViralBottleneck: one row per pair x position.
# Pre-filters are limited to things the package cannot do itself (mask,
# LoFreq support). Depth and frequency thresholds are passed to the package.
vb_select_sites <- function(sites, exclude_pe_ppe = TRUE, require_donor_call = TRUE) {
  d <- sites
  if (require_donor_call) d <- d[d$called_in_donor, ]
  if (exclude_pe_ppe)     d <- d[!d$mask_pe_ppe, ]
  d <- d[!duplicated(d[, c("pair_id", "chrom", "pos")]), ]   # multi-allelic: one row per site
  d
}

vb_sample_table <- function(d, role, segment) {
  data.frame(
    Pos       = as.integer(d$pos),
    segment   = segment,
    A         = as.integer(d[[paste0(role, "_A")]]),
    T         = as.integer(d[[paste0(role, "_T")]]),
    C         = as.integer(d[[paste0(role, "_C")]]),
    G         = as.integer(d[[paste0(role, "_G")]]),
    `Non/Syn` = "Non",                      # not annotated; only used if NonSyn_or_Syn != "All"
    check.names = FALSE
  )
}

# Writes the per-sample CSVs into `dir` and returns the Donor/Recipient table.
write_viralbottleneck_input <- function(sites, dir, pooled = FALSE,
                                        exclude_pe_ppe = TRUE, require_donor_call = TRUE) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  d <- vb_select_sites(sites, exclude_pe_ppe, require_donor_call)
  if (!nrow(d)) stop("No sites left after filtering.")
  write_one <- function(x, name) {
    write.csv(x, file.path(dir, paste0(name, ".csv")), row.names = FALSE)
  }
  if (pooled) {
    # segment = pair_id keeps site IDs unique across pairs
    write_one(vb_sample_table(d, "donor", d$pair_id),     "pooled_donor")
    write_one(vb_sample_table(d, "recipient", d$pair_id), "pooled_recipient")
    return(data.frame(Donor = "pooled_donor", Recipient = "pooled_recipient"))
  }
  # Per pair. Files are named <pair_id>_donor / <pair_id>_recipient because a
  # sample can appear in several pairs with a different site set in each.
  # Note: ViralBottleneck builds IDs as "<Donor>-<Recipient>", so avoid "-" in pair_id.
  pairs <- unique(d$pair_id)
  for (p in pairs) {
    dp <- d[d$pair_id == p, ]
    write_one(vb_sample_table(dp, "donor", "1"),     paste0(p, "_donor"))
    write_one(vb_sample_table(dp, "recipient", "1"), paste0(p, "_recipient"))
  }
  data.frame(Donor = paste0(pairs, "_donor"), Recipient = paste0(pairs, "_recipient"))
}

# Runs one ViralBottleneck method and returns its table with a pair_id column.
#   variant_calling      donor minor-allele threshold AND recipient calling
#                        threshold T (the package uses one value for both)
#   donor/recipient_depth_threshold  minimum site depth in each host
#   Nbmax                upper end of the Nb grid; the exact method costs
#                        ~Nbmax^2 per site, so keep it modest for Mtb
run_viralbottleneck <- function(sites,
                                method                    = "Beta_binomial_Exact",
                                variant_calling           = 0.01,
                                donor_depth_threshold     = 20,
                                recipient_depth_threshold = 20,
                                Nbmin                     = 1,
                                Nbmax                     = 200,
                                pooled                    = FALSE,
                                exclude_pe_ppe            = TRUE,
                                require_donor_call        = TRUE,
                                keep_dir                  = NULL) {
  if (!requireNamespace("ViralBottleneck", quietly = TRUE)) {
    stop('Install first: remotes::install_github("BowenArchaman/ViralBottleneck")')
  }
  dir <- if (is.null(keep_dir)) tempfile("vb_") else keep_dir
  pairs_tbl <- write_viralbottleneck_input(sites, dir, pooled, exclude_pe_ppe, require_donor_call)

  old <- setwd(dir)
  on.exit(setwd(old), add = TRUE)

  ob  <- ViralBottleneck::CreateTransmissionObject(pairs_tbl)
  res <- ViralBottleneck::Bottleneck_size_Calculation(
    transmission_ob           = ob,
    method                    = method,
    donor_depth_threshold     = donor_depth_threshold,
    recipient_depth_threshold = recipient_depth_threshold,
    variant_calling           = variant_calling,
    Nbmin                     = Nbmin,
    Nbmax                     = Nbmax
  )
  res$pair_id <- if (pooled) "pooled" else sub("_donor$", "", res$donor)
  res$method  <- method
  res$variant_calling <- variant_calling
  res[, c("pair_id", setdiff(names(res), "pair_id"))]
}
