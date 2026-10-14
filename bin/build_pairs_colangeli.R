#!/usr/bin/env Rscript
# ---------------------------------------------------------------------------
# build_pairs_colangeli.R
#
# Builds the curated donor-recipient pairs file for the Colangeli 2020 (Vitoria,
# Brazil) household cohort, for use as --pairs input to bottleneck_prep.
#
# Inputs
#   --si    Supplementary Information PDF (Nat Commun 41467_2020_18699_MOESM1_ESM)
#   --sra   SRA run metadata CSV for PRJNA475130, with at least the columns
#           run_accession, library_name, collection_date. Produced by:
#             pysradb metadata PRJNA475130 --detailed --saveto colangeli_sra.csv
#           or the "Metadata" download from the SRA Run Selector.
#   --out   Output TSV path.
#
# What is automated and what is not
#   AUTOMATED  - parsing Supplementary Table 2 (index ID, contact ID, months
#                between diagnoses) out of the SI PDF text layer
#              - joining study IDs to SRA run accessions and collection dates
#              - deriving which pairs had index/secondary re-defined, from the
#                enrolment-ID morphology, and checking the count against the
#                seven stated in the paper's Methods
#   MANUAL     - which pairs have a community-unique RFLP pattern. The paper
#                marks these in bold in Supplementary Table 2 and bold is not
#                recoverable from the PDF text layer; the page was rendered and
#                read visually. The result is the RFLP_UNIQUE constant below.
#                If this script is ever pointed at a revised SI, that constant
#                must be re-checked by eye.
#
# Base R only, no package dependencies.
# ---------------------------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) {
    if (is.null(default)) stop("missing required argument: ", flag, call. = FALSE)
    return(default)
  }
  args[i + 1L]
}
si_pdf   <- get_arg("--si")
sra_csv  <- get_arg("--sra")
out_tsv  <- get_arg("--out")

for (f in c(si_pdf, sra_csv)) {
  if (!file.exists(f)) stop("input not found: ", f, call. = FALSE)
}
if (nchar(Sys.which("pdftotext")) == 0L) {
  stop("pdftotext not found on PATH (install poppler-utils / poppler)", call. = FALSE)
}

# --- Pairs whose RFLP pattern was unique in the Vitoria community ------------
# Read visually from the bold rows of Supplementary Table 2. Named by index case.
# The paper states 14 such pairs; that count is asserted below.
RFLP_UNIQUE <- c("VV3015", "SE3001-2", "SE3063", "VT3010", "VT3050", "CA3064",
                 "CA3077", "VT3046-3", "CA3078-3", "VT3067", "CA3010",
                 "SE3010-2", "VV3012-4", "CA3038")

N_PAIRS       <- 24L  # pairs retained after the paper's exclusion cascade
N_REDEFINED   <- 7L   # Methods: contacts whose TB predated the original index

# --- 1. Parse Supplementary Table 2 ------------------------------------------
si_txt <- system2("pdftotext", c(shQuote(si_pdf), "-"), stdout = TRUE, stderr = FALSE)
si_txt <- paste(si_txt, collapse = "\n")

start <- regexpr("Supplementary Table 2", si_txt, fixed = TRUE)
if (start < 0L) stop("could not locate 'Supplementary Table 2' in the SI text", call. = FALSE)
seg <- substr(si_txt, start, start + 6000L)

# Participant IDs look like VV3015, CA3043-1, SE3001-2: two letters, four
# digits, optional -N enrolment suffix. The table lists all index IDs, then all
# contact IDs, in row order.
id_pat <- "\\b(?:VV|CA|SE|VT)[0-9]{4}(?:-[0-9]+)?\\b"
ids <- regmatches(seg, gregexpr(id_pat, seg, perl = TRUE))[[1]]
if (length(ids) < 2L * N_PAIRS) {
  stop("found only ", length(ids), " participant IDs; expected at least ",
       2L * N_PAIRS, ". The SI layout may have changed.", call. = FALSE)
}
donor_study_id     <- ids[seq_len(N_PAIRS)]
recipient_study_id <- ids[N_PAIRS + seq_len(N_PAIRS)]

# Interval is printed as "years (months)", e.g. "<1 (11)" or "3-4 (45)".
# Capture the parenthesised months.
mon_pat <- "(?:<1|[0-9]-[0-9])\\s*\\(([0-9]+)\\)"
mon_hit <- regmatches(seg, gregexpr(mon_pat, seg, perl = TRUE))[[1]]
months  <- sub(".*\\(([0-9]+)\\).*", "\\1", mon_hit)
if (length(months) < N_PAIRS) {
  stop("found only ", length(months), " interval values; expected ", N_PAIRS,
       call. = FALSE)
}
months <- months[seq_len(N_PAIRS)]

stopifnot(
  length(donor_study_id)     == N_PAIRS,
  length(recipient_study_id) == N_PAIRS,
  !anyDuplicated(c(donor_study_id, recipient_study_id))
)

# --- 2. SRA metadata ---------------------------------------------------------
sra <- utils::read.csv(sra_csv, stringsAsFactors = FALSE, check.names = FALSE)
need <- c("run_accession", "library_name", "collection_date")
missing_cols <- setdiff(need, names(sra))
if (length(missing_cols)) {
  stop("SRA file is missing column(s): ", paste(missing_cols, collapse = ", "),
       "\nRe-export with `pysradb metadata ... --detailed`, which includes the ",
       "submitter-supplied fields.", call. = FALSE)
}
sra$library_name <- trimws(sra$library_name)

dups <- sra$library_name[duplicated(sra$library_name)]
if (length(dups)) {
  stop("library_name maps to more than one run for: ", paste(unique(dups), collapse = ", "),
       "\nRuns may be split across lanes; decide how to merge before pairing.",
       call. = FALSE)
}

run_of  <- setNames(sra$run_accession,   sra$library_name)
date_of <- setNames(sra$collection_date, sra$library_name)

wanted  <- c(donor_study_id, recipient_study_id)
unmatched <- setdiff(wanted, names(run_of))
if (length(unmatched)) {
  stop("no SRA run found for study ID(s): ", paste(unmatched, collapse = ", "),
       "\nCheck that library_name in the SRA export carries the study IDs.",
       call. = FALSE)
}
# Runs present in the project but not used: expected to include the pair the
# paper excluded after sequencing. Reported, not an error.
unused <- setdiff(names(run_of), wanted)

blank_to_na <- function(x) ifelse(is.na(x) | x == "", "NA", x)

# --- 3. Direction ------------------------------------------------------------
# In the seven pairs where the Methods re-defined index and secondary, the donor
# carries the -N contact enrolment suffix and the recipient has the bare stem.
has_suffix <- function(x) grepl("-[0-9]+$", x)
redefined  <- has_suffix(donor_study_id) & !has_suffix(recipient_study_id)
if (sum(redefined) != N_REDEFINED) {
  stop("derived ", sum(redefined), " re-defined pairs but the Methods state ",
       N_REDEFINED, ". Do not trust the direction column until this is resolved.",
       call. = FALSE)
}

if (!all(RFLP_UNIQUE %in% donor_study_id)) {
  stop("RFLP_UNIQUE names an index case not present in the parsed table: ",
       paste(setdiff(RFLP_UNIQUE, donor_study_id), collapse = ", "), call. = FALSE)
}

# --- 4. Assemble -------------------------------------------------------------
notes <- rep("IDs+interval from Suppl Table 2; direction IC->HHC by design", N_PAIRS)
notes[redefined] <- paste0(notes[redefined],
  "; donor enrolled as HHC but TB predated original index - re-defined per Methods")

d_date <- blank_to_na(unname(date_of[donor_study_id]))
r_date <- blank_to_na(unname(date_of[recipient_study_id]))

# Flag pairs where the recipient specimen predates the donor specimen. This is
# not itself a direction error -- the paper's rule uses case occurrence, not
# sampling -- but it marks the pairs with the least independent support.
gap <- suppressWarnings(as.integer(as.Date(r_date) - as.Date(d_date)))
odd <- !is.na(gap) & gap < 0L
notes[odd] <- paste0(notes[odd], "; CHECK: recipient sampled ", abs(gap[odd]),
                     " day(s) BEFORE donor, ", months[odd],
                     " months between diagnoses")

pairs <- data.frame(
  pair_id                   = sprintf("COL_%03d", seq_len(N_PAIRS)),
  cohort                    = "Colangeli2020",
  donor_sample_id           = unname(run_of[donor_study_id]),
  recipient_sample_id       = unname(run_of[recipient_study_id]),
  donor_study_id            = donor_study_id,
  recipient_study_id        = recipient_study_id,
  link_type                 = "household_index_to_contact",
  months_between_diagnoses  = months,
  donor_collection_date     = d_date,
  recipient_collection_date = r_date,
  rflp_unique_in_community  = ifelse(donor_study_id %in% RFLP_UNIQUE, "yes", "no"),
  direction_redefined       = ifelse(redefined, "yes", "no"),
  notes                     = notes,
  stringsAsFactors          = FALSE
)

stopifnot(
  nrow(pairs) == N_PAIRS,
  sum(pairs$rflp_unique_in_community == "yes") == length(RFLP_UNIQUE),
  sum(pairs$direction_redefined == "yes")      == N_REDEFINED,
  !anyDuplicated(c(pairs$donor_sample_id, pairs$recipient_sample_id))
)

dir.create(dirname(out_tsv), showWarnings = FALSE, recursive = TRUE)
utils::write.table(pairs, out_tsv, sep = "\t", quote = FALSE,
                   row.names = FALSE, na = "NA")

message(sprintf("wrote %s: %d pairs, %d unique-RFLP, %d direction re-defined",
                out_tsv, nrow(pairs),
                sum(pairs$rflp_unique_in_community == "yes"),
                sum(pairs$direction_redefined == "yes")))
if (length(unused)) {
  message("SRA runs in the project not used by any pair (expected: the pair ",
          "excluded after sequencing): ", paste(unused, collapse = ", "))
}
