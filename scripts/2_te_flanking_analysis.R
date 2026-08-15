# =====================================================================
# Cichlid TE Flanking-Region Intersection
# ---------------------------------------------------------------------
# For one species, extends each single-copy ortholog by +/-10 kb,
# intersects the flanking windows against the RepeatMasker TE annotation
# with GenomicRanges, computes per-ortholog TE statistics, and writes a
# master table plus a set of plots and complementary data tables.
#
# Run once per species; edit the USER INPUTS block below as needed
# =====================================================================

# ---- USER INPUTS: change these for each species run ----
te_file         <- "data/TE_Oreochromis_niloticus_AF_genome.fa.out"
orthogroup_file <- "data/OUTPUT_TABLE_orthogroup_gene_info.tsv"
species_name    <- "Oreochromis_niloticus_AF"
output_dir      <- "insights/Oreochromis_niloticus_AF_intersection_outputs"
gff3_file       <- "data/GENE_Oreochromis_niloticus_AF.annotation.gff3"  # optional, for contig lengths


NULL_SEED <- 42   # fixed so the sampling is reproducible between runs
NULL_N    <- NA   # NA = match the number of orthologs retained; or set an integer
# --------------------------------------------------------

# ---- Analysis constants ----
FLANK_BP   <- 10000   # window size per side
BIN_BP     <- 500     # bin width for the positional density plots
# Classes (top-level class = text before the first '/' for a given TE classification) are discovered
# automatically from each te_file input file -- every class found gets its own
# pct_<class> column and nothing is folded into an "Other" bucket. You never
# need to edit this for a new annotation file.
#
# CLASS_ORDER is ONLY a display-ordering preference: any of these classes that
# are present appear first (in this order), and any other classes found in the
# file are appended after them alphabetically. Listing a class here that isn't
# in the data does nothing; a class in the data that isn't listed here is still
# included automatically.
CLASS_ORDER <- c("DNA", "LTR", "LINE", "SINE", "RC", "Unknown", "Simple_repeat", "Low_complexity", "tRNA", "rRNA", "snRNA")


main <- function() {

  # Loud, visually distinct banner so it is always obvious which species
  species_banner(species_name)
  message("=== Cichlid TE flanking-region analysis: ", species_name, " ===")

  ## ---- 0. Dependencies & input validation -------------------------
  check_packages()
  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(IRanges)
    library(ggplot2)
    library(dplyr)
    library(data.table)
  })
  # Ensure required inputs are given and non-zero
  for (f in c(te_file, orthogroup_file)) {
    if (!file.exists(f)) stop("Required input file not found: ", f, call. = FALSE)
  }
  use_gff3 <- nzchar(gff3_file) && file.exists(gff3_file)
  if (!use_gff3) {
    message("NOTE: GFF3 not found/!set -- downstream windows will be capped at the ",
            "max TE coordinate per contig instead of true contig length.")
  }

  ## ---- 1. Parse TE annotation -------------------------------------
  # This parses the .out TE annotation files and returns a GRanges object that represents the entire TE annotation for the species
  message("Parsing TE annotation file...")
  te_gr <- parse_te_file(te_file)
  message(sprintf("Loaded %s TE annotations", format(length(te_gr), big.mark = ",")))

  ## ---- 1b. Drop sub-100 bp TE fragments ------------------
  # Sub-100 bp fragments are mostly noise / spurious annotations. GRanges already
  # knows each TE's length via width(), so drop anything narrower than 100 bp.
  keep_te   <- width(te_gr) >= 100L
  n_dropped <- sum(!keep_te)
  te_gr     <- te_gr[keep_te]
  message(sprintf("Applied 100 bp filter: dropped %s TEs (< 100 bp), %s remain",
                  format(n_dropped, big.mark = ","),
                  format(length(te_gr), big.mark = ",")))

  ## ---- 2. Contig lengths ------------------------------------------
  if (use_gff3) {
    message("Reading contig lengths from GFF3 (##sequence-region)...")
    contig_len <- parse_contig_lengths(gff3_file) # Returns a file containing all the contigs and the lengths of each one
    message(sprintf("Found lengths for %s contigs", format(length(contig_len), big.mark = ",")))
  } else { # If there is no gene annotation file provided, each contig length is estimated taking the maximum end coordinate of the rightmost TE annotation on each contig
    contig_len <- tapply(end(te_gr), as.character(seqnames(te_gr)), max)
  }

  ## ---- 2b. Total genome length (for EDIT 5: genome-wide breakdown) --
  # A trustworthy GFF3 gives one ##sequence-region per contig, so summing those
  # is the true genome length. Any species with a bad/empty GFF3 falls back to summing the per-contig max TE end, which
  # under-estimates slightly. genome_len_est flags when the fallback was used so
  # the limitation can be surfaced in messages and on the plot.
  if (use_gff3 && length(contig_len) > 1L) { # This just checks, if the gff file contains the typical ##sequence-region line per contig, each giving that contigs true length, just add them up to get the total genome length
    total_genome_len <- sum(as.numeric(contig_len), na.rm = TRUE)
    genome_len_est   <- FALSE
  } else {
    total_genome_len <- sum(tapply(end(te_gr), as.character(seqnames(te_gr)), max))
    genome_len_est   <- TRUE
    message("NOTE: total genome length ESTIMATED from max TE end per contig ",
            "(no usable ##sequence-region header) -- genome-wide percentages are approximate.")
  }
  message(sprintf("Total genome length%s: %s bp",
                  if (genome_len_est) " (estimated)" else "",
                  format(total_genome_len, big.mark = ",")))

  ## ---- 3. Parse ortholog table ------------------------------------
  message("Filtering orthologs for species: ", species_name, "...")
  orth <- parse_orthologs(orthogroup_file, species_name) # Returns a data table with only ortholog info for current species of interest
  message(sprintf("Found %s orthologs", format(nrow(orth), big.mark = ",")))

  ## ---- 3b. Remove orthologs whose gene BODY is >=90% TE (EDIT 2) ---
  # some "orthologs" are really transposon-derived
  # -- their own gene body is almost entirely TE. Intersect each gene body itself
  # (gene_start..gene_end) against the 100 bp-filtered TE set, and if >=90% of the
  # gene-body length is TE, drop it from ALL downstream analysis. Coverage is
  # measured against the GENE BODY, not the flank, and uses merged bp (no
  # double-counting). Discarded orthologs are logged to a table.
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  message("Screening gene bodies for TE-derived orthologs (>=90% TE)...")
  te_flag <- flag_te_orthologs(orth, te_gr, threshold = 0.90)
  if (nrow(te_flag$discarded)) {
    disc_tsv <- file.path(output_dir, "discarded_TE_orthologs.tsv")
    data.table::fwrite(te_flag$discarded, disc_tsv, sep = "\t")
    message(sprintf("Removed %s TE-derived orthologs (>=90%% gene body is TE); log: %s",
                    format(nrow(te_flag$discarded), big.mark = ","), disc_tsv))
  } else {
    message("No orthologs exceeded the 90% gene-body TE threshold.")
  }
  orth <- orth[!te_flag$is_te, , drop = FALSE]
  message(sprintf("%s orthologs retained after TE-ortholog removal",
                  format(nrow(orth), big.mark = ",")))

  ## ---- 4. Build flanking windows (+/-10kb) ------------------------
  message("Building flanking windows (+/-", FLANK_BP, "bp)...")
  fw <- build_flank_windows(orth, contig_len) # returns two things packaged into one list - the win_info dataframe storing window coordinates and sizes for both upstream and downstream flanks + a windows_gr GRanges objects with one row per valid window - storing coordinates, which ortholog it belongs to etc.
  message(sprintf("Built %s flanking windows over %s orthologs",
                  format(length(fw$windows_gr), big.mark = ","),
                  format(nrow(orth), big.mark = ",")))

  ## ---- 5. Intersection --------------------------------------------
  message("Running findOverlaps()...")
  hits <- intersect_windows(fw$windows_gr, te_gr) # a data table where each row is one window-TE overlap, with all the info you need about it. This is what gets passed to compute_stats().
  message(sprintf("Found %s window-TE overlaps", format(nrow(hits), big.mark = ",")))

  ## ---- 6. Per-ortholog statistics ---------------------------------
  message("Computing TE statistics per ortholog...")
  master <- compute_stats(orth, fw$win_info, hits) # Giving the original orthogroup info table for our species, the win_info table and the TE genomic ranges object - run a bunch of stats and create the master data table which contains a bunch of information about the flanking regions of each ortholog

  ## ---- 7. Write master table --------------------------------------
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  out_tsv <- file.path(output_dir, "master_table.tsv")
  message("Writing master table...")
  data.table::fwrite(master, out_tsv, sep = "\t")
  message(sprintf("Wrote %s rows to %s", format(nrow(master), big.mark = ","), out_tsv))

  ## ---- 8. Random baselines  ---------------
  # Two comparison sets, sampled at the same size as the ortholog set and written
  # as small TSVs. The cross-species density FIGURE is drawn by the aggregation
  # script, which needs all four species at once.
  n_null <- if (is.na(NULL_N)) nrow(orth) else as.integer(NULL_N)
  message(sprintf("Sampling %s random %s bp genome chunks...",
                  format(n_null, big.mark = ","), 2L * FLANK_BP))
  rc <- sample_random_chunks(contig_len, te_gr, n_null)
  data.table::fwrite(rc, file.path(output_dir, "null_random_chunks.tsv"), sep = "\t")
  message(sprintf("  median %%TE, random chunks: %.1f%%", median(rc$pct_te, na.rm = TRUE)))

  if (use_gff3) {
    message(sprintf("Sampling %s random genes from the GFF3...",
                    format(n_null, big.mark = ",")))
    rg <- sample_random_genes(gff3_file, contig_len, te_gr, n_null)
    data.table::fwrite(rg, file.path(output_dir, "null_random_genes.tsv"), sep = "\t")
    message(sprintf("  median %%TE, random gene flanks: %.1f%%",
                    median(rg$pct_te, na.rm = TRUE)))
  } else {
    message("NOTE: gff3_file not set -- skipping the random-GENE baseline ",
            "(random chunks still written).")
  }

  ## ---- 9. Plots ---------------------------------------------------
  message("Generating plots...")
  make_plots(master, hits, output_dir, species_name,
             te_gr, total_genome_len, genome_len_est)

  message("Done. Output written to: ", normalizePath(output_dir))
}


# =====================================================================
# Helper functions
# =====================================================================

# Print a loud, obvious banner (via cat, so it is visually distinct from the
# normal message() lines) announcing which species/run everything below belongs
# to, e.g.  ======== METRIACLIMA_ZEBRA_AF ========
species_banner <- function(txt) {
  bar <- paste(rep("=", 60), collapse = "")
  cat("\n", bar, "\n======== ", toupper(txt), " ========\n", bar, "\n\n", sep = "")
}

# Function to check if all packages installed, if not, tells you to install them
check_packages <- function() {
  needed <- list(
    GenomicRanges = "Bioconductor", IRanges = "Bioconductor",
    ggplot2 = "CRAN", dplyr = "CRAN", data.table = "CRAN")
  missing <- needed[!vapply(names(needed), requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    msg <- paste(sprintf("  - %s (install from %s)", names(missing), unlist(missing)),
                 collapse = "\n")
    stop("Missing required packages:\n", msg,
         "\nBioconductor: install.packages('BiocManager'); BiocManager::install(c('GenomicRanges','IRanges'))",
         "\nCRAN: install.packages(c('ggplot2','dplyr','data.table'))", call. = FALSE)
  }
}

# Parse a RepeatMasker .out file into a GRanges.
# Format: 15 whitespace-delimited cols (+ optional trailing '*' on
# overlapping rows). 2 header lines + 1 blank line are skipped.

# BASICALLY - the .out file has TWO header lines, not the usual one. Because of this strangeness, this script skips BOTH headers and the blank line before the actual content (thus why it's skip = 3) - opting to give it's own column names. We only keep the necessary columns, rest are discarded for our purposes.
parse_te_file <- function(path) {
  d <- data.table::fread(
    path, skip = 3, header = FALSE, sep = " ", quote = "",
    fill = TRUE, strip.white = TRUE, blank.lines.skip = TRUE,
    select = c(5, 6, 7, 9, 10, 11),
    col.names = c("chrom", "start", "end", "strand", "repeat_name", "class_family"))

  # RepeatMasker strand: '+' = plus, 'C' = minus (complement).
  strand <- ifelse(d$strand == "C", "-", "+")

  # TE class = top-level class, the text before the first '/'; if there is no
  # '/', the whole string (e.g. Simple_repeat, Unknown, Low_complexity).
  

  te_class  <- sub("/.*", "", d$class_family)

  # Take all parsed columns and represent them with a GRanges object
  GRanges(
    seqnames = d$chrom,
    ranges   = IRanges(start = d$start, end = d$end),
    strand   = strand,
    repeat_name = d$repeat_name,
    te_family   = d$class_family,
    te_class    = te_class)
}

# reads the GFF3 gene annotation file to get the length of each chromosome/contig

# Pull contig lengths from the ##sequence-region directives. 
# There is one at the start of each chrsmone/contig in the file. These are
# scattered through the (large) GFF3, so grep them out rather than
# reading the whole file into R (however, if grep fails, it falls bckack to reading the GFF3 in chunks of 5000 and pulling out those lines manually)
parse_contig_lengths <- function(path) {
  lines <- tryCatch( # try the following, but if it fails, don't crash
    readLines(pipe(paste("grep '^##sequence-region'", shQuote(path)))),
    warning = function(w) character(0), error = function(e) character(0)) # if failure - return an empty vector
  if (!length(lines)) {  # fallback: scan header lines directly (if vector is empty, meaning it failed)
    con <- file(path, "r"); on.exit(close(con)) # read the file, close upon completion
    lines <- character(0) # 
    repeat {
      chunk <- readLines(con, n = 5000)
      if (!length(chunk)) break
      sr <- grep("^##sequence-region", chunk, value = TRUE)
      lines <- c(lines, sr)
      if (any(!startsWith(chunk, "#")) && !length(sr)) break
    }
  }
  parts <- strsplit(trimws(lines), "\\s+")
  setNames(
    as.numeric(vapply(parts, `[`, character(1), 4)),
    vapply(parts, `[`, character(1), 2))
}

parse_orthologs <- function(path, species) {
  d <- data.table::fread(path, sep = "\t", header = TRUE, data.table = FALSE)
  d <- d[d$species == species, ] # Filter the entire OUTPUT_TABLE_orthogroup... info table down to only our current species of interest
  if (!nrow(d)) stop("No orthologs found for species '", species,
                     "'. Check species_name spelling.", call. = FALSE)
  
  # Create a new data frame with only columns that are needed
  data.frame(
    orthogroup = d$orthogroup, species = d$species,
    gene_name  = d$gene_name,  protein_id = d$protein_id,
    chromosome = d$chromosome,
    gene_start = as.integer(d$start), gene_end = as.integer(d$end),
    gene_strand = d$strand, stringsAsFactors = FALSE)
}

# flag orthologs whose OWN gene body is essentially a TE.
# For each ortholog, intersect its gene body (gene_start..gene_end on its
# chromosome) with the 100 bp-filtered TE set, compute the fraction of the
# gene-body length covered by TE using MERGED bp (reduced, so bp are never
# double-counted), and mark it for removal if that fraction is >= threshold.
# Returns a logical mask (is_te, one per input ortholog, same order) plus a
# data.frame of the discarded orthologs with the supporting numbers.
flag_te_orthologs <- function(orth, te_gr, threshold = 0.90) {
  n <- nrow(orth)
  gene_len <- orth$gene_end - orth$gene_start + 1L
  gene_gr  <- GRanges(
    seqnames = orth$chromosome,
    ranges   = IRanges(start = orth$gene_start, end = orth$gene_end),
    og_idx   = seq_len(n))

  te_bp <- numeric(n)
  ov <- findOverlaps(gene_gr, te_gr, ignore.strand = TRUE)
  if (length(ov)) {
    qh <- queryHits(ov); sh <- subjectHits(ov)
    # clip each TE to the gene body so nothing outside the body is counted
    clip_start <- pmax(start(gene_gr)[qh], start(te_gr)[sh])
    clip_end   <- pmin(end(gene_gr)[qh],   end(te_gr)[sh])
    bp <- reduced_bp_by_group(clip_start, clip_end, qh)  # merged bp per gene body
    te_bp[as.integer(names(bp))] <- bp
  }

  frac  <- te_bp / gene_len
  is_te <- frac >= threshold

  discarded <- data.frame(
    orthogroup        = orth$orthogroup[is_te],
    gene_name         = orth$gene_name[is_te],
    chromosome        = orth$chromosome[is_te],
    gene_start        = orth$gene_start[is_te],
    gene_end          = orth$gene_end[is_te],
    gene_length_bp    = gene_len[is_te],
    te_bp_covered     = te_bp[is_te],
    coverage_fraction = round(frac[is_te], 4),
    stringsAsFactors  = FALSE)

  list(is_te = is_te, discarded = discarded)
}

# Build upstream and downstream flank windows, capping at contig
# boundaries [1, contig_length]. Returns the windows as a GRanges
# (one entry per non-empty flank) plus a per-ortholog win_info frame.
build_flank_windows <- function(orth, contig_len) {
  n <- nrow(orth) # Store the number of orthologs
  clen <- unname(contig_len[orth$chromosome]) # for each ortholog, look up its contig's (note the chromosome name here is a bit of a misnomer, it can mean contig as well) length - so clen is a vector of same number of lengths as n, one per ortholog, in the same order. unname() strips the vector names: a lookup miss returns an NA name, and those NA names otherwise propagate into the GRanges below and trigger "missing values in 'row.names'".
  clen[is.na(clen)] <- .Machine$integer.max  # unknown contig: no upper cap

  # This gets the up_start, up_end (basically the entire 10 000 base pair sequence before the ortholog starts) and the down_start and down_end - each of these is a vector containing all of these positions for every single ortholog - so if there are 8000 orthologs total, each of these variables will have 8000 vectors for ever single one of them
  up_start   <- pmax(1L, orth$gene_start - FLANK_BP)
  up_end     <- orth$gene_start - 1L
  down_start <- orth$gene_end + 1L
  down_end   <- pmin(clen, orth$gene_end + FLANK_BP)

  # Check if the upstream window is valid — i.e. it actually has at least 1bp. Could be invalid if the gene starts at position 1 (no upstream room at all). Same thing done with downstream.
  up_valid   <- up_end   >= up_start
  down_valid <- down_end >= down_start

  # Calculate actual size in bp of each window. If invalid, set to 0.
  up_bp   <- ifelse(up_valid,   up_end   - up_start   + 1L, 0L)
  down_bp <- ifelse(down_valid, down_end - down_start + 1L, 0L)

  # Build a data frame storing all the window coordinates and sizes for every ortholog. og_idx is just the row number (1 to 8,334) so you can match back to the ortholog later. Invalid windows get NA for their coordinates.
  win_info <- data.frame(
    og_idx = seq_len(n),
    upstream_window_start = ifelse(up_valid, up_start, NA_integer_),
    upstream_window_end   = ifelse(up_valid, up_end,   NA_integer_),
    upstream_actual_bp    = as.integer(up_bp),
    downstream_window_start = ifelse(down_valid, down_start, NA_integer_),
    downstream_window_end   = ifelse(down_valid, down_end,   NA_integer_),
    downstream_actual_bp    = as.integer(down_bp),
    stringsAsFactors = FALSE)

  # Assemble the GRanges of valid windows.
  idx    <- c(which(up_valid),               which(down_valid)) # Get the indices of all valid upstream and downstream windows, combined into one vector
  
  # Get the actual start and end coordinates of all valid windows
  starts <- c(up_start[up_valid],            down_start[down_valid])
  ends   <- c(up_end[up_valid],              down_end[down_valid])
  
  # Label each window as upstream or downstream so its known which side of the gene it came from
  flank  <- c(rep("upstream", sum(up_valid)), rep("downstream", sum(down_valid)))
  
# Package all valid windows into a GRanges object. Each row is one window (upstream or downstream), with its chromosome, coordinates, which ortholog it belongs to (og_idx), which flank it is, and the gene's own start/end (needed later for the density plot).
  windows_gr <- GRanges(
    seqnames = orth$chromosome[idx],
    ranges   = IRanges(start = starts, end = ends),
    og_idx   = idx,
    flank    = flank,
    gene_start = orth$gene_start[idx],
    gene_end   = orth$gene_end[idx])

  list(windows_gr = windows_gr, win_info = win_info) # This function returns a list containing two named objects - the GRanges windows object and the actual data frame with all windows information
}

# Overlap windows vs TEs, clip each TE to its window, return a tidy
# data.table of clipped hits (one row per window-TE overlap).
intersect_windows <- function(windows_gr, te_gr) {
  ov <- findOverlaps(windows_gr, te_gr, ignore.strand = TRUE) # Finds every place where a flanking window overlaps a TE and returns a 'hits' object where each row is one window-TE overlap pair. 
  qh <- queryHits(ov); sh <- subjectHits(ov) # qh = the index of the window in each overlap. sh = the index of the TE in each overlap. So if window 5 overlaps TE 203, you'd have qh = 5, sh = 203 as one row.

  win_start <- start(windows_gr)[qh]; win_end <- end(windows_gr)[qh] #Get the start and end coordinates of the window involved in each overlap.
  te_start  <- start(te_gr)[sh];      te_end  <- end(te_gr)[sh] # Get the start and end coordinates of the TE involved in each overlap.

  data.table::data.table(
    og_idx     = windows_gr$og_idx[qh],             # Which ortholog this overlap belongs to.
    flank      = windows_gr$flank[qh],              # Whether it's upstream or downstream.
    te_id      = sh,                                # unique TE annotation id
    te_class   = te_gr$te_class[sh],                # The TE class (DNA, LTR, LINE etc.).
    clip_start = pmax(win_start, te_start),         # clip to window boundary (ie. a TE may extend beyong an overlap, so we take the overlap region only)
    clip_end   = pmin(win_end,   te_end),
    
    # These three things, the start and end of gene and the midpoint of the TE are used in TE density later
    gene_start = windows_gr$gene_start[qh],
    gene_end   = windows_gr$gene_end[qh])
    # te_mid     = (te_start + te_end) / 2)
}

# Sum of merged (reduced) interval widths per group -- collapses
# overlapping annotations so bp are not double-counted. All ranges in
# one group share a chromosome (grouped within an ortholog), so plain
# IRanges is sufficient.
reduced_bp_by_group <- function(start, end, group) {
  ir  <- IRanges(start = start, end = end)
  red <- reduce(split(ir, group))
  vapply(width(red), sum, numeric(1))
}

# Merged (non-double-counted) TE bp per group, for ANY set of windows carrying an
# og_idx. Same logic the ortholog flanks use, so the baselines are comparable:
# same 100 bp-filtered te_gr, same reduce() so overlapping TEs count once.
merged_te_bp_by_group <- function(gr, te_gr, n_groups) {
  out <- numeric(n_groups)
  ov  <- findOverlaps(gr, te_gr, ignore.strand = TRUE)
  qh  <- queryHits(ov); sh <- subjectHits(ov)
  if (!length(qh)) return(out)
  grp <- factor(gr$og_idx[qh], levels = seq_len(n_groups))
  cs  <- pmax(start(gr)[qh], start(te_gr)[sh])
  ce  <- pmin(end(gr)[qh],   end(te_gr)[sh])
  out[] <- as.numeric(reduced_bp_by_group(cs, ce, grp))
  out
}

# sample n random windows of `size` bp from anywhere in
# the genome, weighted by contig length so long contigs are hit proportionally.
# Windows may land on genes
sample_random_chunks <- function(contig_len, te_gr, n, size = 2L * FLANK_BP,
                                 seed = NULL_SEED) {
  set.seed(seed)
  cl <- contig_len[!is.na(contig_len) & contig_len >= size]
  if (!length(cl))
    stop("No contig is long enough to sample a ", size, " bp chunk.", call. = FALSE)
  room <- as.numeric(cl) - size + 1
  ctg  <- sample(names(cl), n, replace = TRUE, prob = room / sum(room))
  st   <- as.integer(floor(runif(n) * (as.numeric(cl[ctg]) - size + 1)) + 1)
  gr   <- GRanges(seqnames = ctg,
                  ranges   = IRanges(start = st, width = size),
                  og_idx   = seq_len(n))
  bp   <- merged_te_bp_by_group(gr, te_gr, n)
  data.frame(chunk_id = seq_len(n), chromosome = ctg,
             start = st, end = st + size - 1L,
             window_bp = size, te_bp_covered = bp,
             pct_te = 100 * bp / size, stringsAsFactors = FALSE)
}

# read every gene row out of a GFF3 (col 3 == "gene").
# Chunked readLines rather than awk/grep so this also works on Windows.
parse_gff_genes <- function(path) {
  con <- file(path, "r"); on.exit(close(con))
  keep <- list(); i <- 0L
  repeat {
    chunk <- readLines(con, n = 200000L, warn = FALSE)
    if (!length(chunk)) break
    hit <- grep("\tgene\t", chunk, value = TRUE, fixed = TRUE)
    if (length(hit)) { i <- i + 1L; keep[[i]] <- hit }
  }
  keep <- unlist(keep, use.names = FALSE)
  if (!length(keep)) stop("No gene rows found in ", path, call. = FALSE)
  f <- strsplit(keep, "\t", fixed = TRUE)
  data.frame(chromosome = vapply(f, `[`, character(1), 1L),
             gene_start = as.integer(vapply(f, `[`, character(1), 4L)),
             gene_end   = as.integer(vapply(f, `[`, character(1), 5L)),
             stringsAsFactors = FALSE)
}

# sample n genes at random from the FULL annotation (not
# just the single-copy orthologs), build the same +/-FLANK_BP flanks, and report
# %TE per gene. Tyler's baseline -- fairer than random chunks, because gene
# neighbourhoods are structurally unlike random genome.
sample_random_genes <- function(gff3_file, contig_len, te_gr, n, seed = NULL_SEED) {
  genes <- parse_gff_genes(gff3_file)
  set.seed(seed + 1L)
  g  <- genes[sample(seq_len(nrow(genes)), min(n, nrow(genes))), , drop = FALSE]
  ng <- nrow(g)

  clen <- contig_len[g$chromosome]
  clen[is.na(clen)] <- .Machine$integer.max     # unknown contig: no upper cap
  up_s <- pmax(1L, g$gene_start - FLANK_BP); up_e <- g$gene_start - 1L
  dn_s <- g$gene_end + 1L;                   dn_e <- pmin(clen, g$gene_end + FLANK_BP)
  uv <- up_e >= up_s; dv <- dn_e >= dn_s

  idxw <- c(which(uv), which(dv))
  gr <- GRanges(seqnames = g$chromosome[idxw],
                ranges   = IRanges(start = as.integer(c(up_s[uv], dn_s[dv])),
                                   end   = as.integer(c(up_e[uv], dn_e[dv]))),
                og_idx   = idxw)
  flank_bp <- ifelse(uv, up_e - up_s + 1L, 0L) + ifelse(dv, dn_e - dn_s + 1L, 0L)
  bp <- merged_te_bp_by_group(gr, te_gr, ng)
  data.frame(gene_idx = seq_len(ng), chromosome = g$chromosome,
             gene_start = g$gene_start, gene_end = g$gene_end,
             total_flank_bp = flank_bp, te_bp_covered = bp,
             pct_te = ifelse(flank_bp > 0, 100 * bp / flank_bp, NA_real_),
             stringsAsFactors = FALSE)
}

# This function takes the hits table and computes all the statistics for each ortholog, then assembles the master table.This function takes the hits table and computes all the statistics for each ortholog, then assembles the master table.
# It calculates total flank size for each ortholog, total BP covered by TEs per ortholog, number of distinct TE annotations per ortholog, BP covered per TE class per ortholog etc.
# It then returns a master data frame with all this information
compute_stats <- function(orth, win_info, hits) {
  n <- nrow(orth)

  total_flank_bp <- win_info$upstream_actual_bp + win_info$downstream_actual_bp

  # --- TE bp covered per ortholog (both flanks merged) ---
  te_bp <- numeric(n)
  if (nrow(hits)) {
    bp <- reduced_bp_by_group(hits$clip_start, hits$clip_end, hits$og_idx)
    te_bp[as.integer(names(bp))] <- bp
  }

  # --- distinct TE annotations per ortholog ---
  n_te <- integer(n)
  if (nrow(hits)) {
    cnt <- hits[, .(k = data.table::uniqueN(te_id)), by = og_idx]
    n_te[cnt$og_idx] <- cnt$k
  }

  # --- per-class bp covered per ortholog ---
  # Discover the classes present in this file. Preferred (CLASS_ORDER) ones
  # come first when present; any others are appended alphabetically. Every
  # class found gets a column -- nothing is dropped or folded.
  present <- if (nrow(hits)) sort(unique(hits$te_class)) else character(0)
  classes <- c(intersect(CLASS_ORDER, present), setdiff(present, CLASS_ORDER))
  message("  TE classes found (", length(classes), "): ",
          paste(classes, collapse = ", "))
  class_mat <- matrix(0, nrow = n, ncol = length(classes),
                      dimnames = list(NULL, classes))
  if (nrow(hits)) {
    key <- paste(hits$og_idx, hits$te_class, sep = "")
    bpc <- reduced_bp_by_group(hits$clip_start, hits$clip_end, key)
    kk  <- do.call(rbind, strsplit(names(bpc), "", fixed = TRUE))
    rows <- as.integer(kk[, 1])
    cols <- match(kk[, 2], classes)
    class_mat[cbind(rows, cols)] <- bpc
  }

  total_safe <- ifelse(total_flank_bp > 0, total_flank_bp, NA_real_)

  master <- data.frame(
    orthogroup = orth$orthogroup, species = orth$species,
    gene_name  = orth$gene_name,  protein_id = orth$protein_id,
    chromosome = orth$chromosome,
    gene_start = orth$gene_start, gene_end = orth$gene_end,
    gene_strand = orth$gene_strand,
    win_info[, -1],                                 # drop og_idx column
    total_flank_bp = as.integer(total_flank_bp),
    flank_complete = total_flank_bp == (2L * FLANK_BP),
    n_te_annotations = n_te,
    te_bp_covered = te_bp,
    pct_flank_te  = round(te_bp / total_safe * 100, 4),
    te_density_per_kb = round(n_te / (total_safe / 1000), 4),
    stringsAsFactors = FALSE)

  # Per-class percentage columns.
  for (cl in classes) {
    master[[paste0("pct_", cl)]] <- round(class_mat[, cl] / total_safe * 100, 4)
  }
  # Orthologs with no flank (total_flank_bp == 0) -> NA pct; set counts 0.
  master$pct_flank_te[is.na(master$pct_flank_te)] <- 0
  master
}

# Build a stable colour map for TE classes: fixed colours for the known
# classes, an auto-generated qualitative palette for any extras, plus grey
# for the non-TE remainder.
te_class_palette <- function(classes) {
  known <- c(DNA = "#E64B35", LTR = "#4DBBD5", LINE = "#00A087",
             SINE = "#3C5488", RC = "#F39B7F", Unknown = "#8491B4",
             Simple_repeat = "#91D1C2", Low_complexity = "#DC0000",
             tRNA = "#7E6148", rRNA = "#B09C85", snRNA = "#FDB462")
  extra <- setdiff(classes, names(known))
  if (length(extra))
    known <- c(known, setNames(hcl.colors(length(extra), "Dark 3"), extra))
  c(known[classes], `non-TE` = "grey85")
}

# per-gene, per-bin TE COVERAGE FRACTION along the flanks.
# For every window-TE overlap, convert its clipped coordinates into signed
# distance-from-gene, merge overlapping TEs within each (gene, side) so coverage
# can never exceed the bin, split the merged intervals across BIN_BP bins, and
# return one row per (gene, bin) with cov_frac in [0, 1]. x is signed: negative
# upstream, positive downstream. Genes/bins with no TE simply have no row (=0).
per_gene_bin_coverage <- function(hits, flank_bp = FLANK_BP, bin_bp = BIN_BP) {
  empty <- data.frame(og_idx = integer(0), x = numeric(0), cov_frac = numeric(0))
  if (!nrow(hits)) return(empty)

  # absolute distance interval [s, e] (1..flank_bp) for one side of one hit
  side_ivl <- function(h, side) {
    if (!nrow(h)) return(NULL)
    if (side == "upstream") {          # distance measured back from gene_start
      d_near <- h$gene_start - h$clip_end
      d_far  <- h$gene_start - h$clip_start
    } else {                           # distance measured forward from gene_end
      d_near <- h$clip_start - h$gene_end
      d_far  <- h$clip_end   - h$gene_end
    }
    s <- pmax(1L, as.integer(d_near)); e <- pmin(as.integer(flank_bp), as.integer(d_far))
    ok <- e >= s
    data.frame(og_idx = h$og_idx[ok], side = side, s = s[ok], e = e[ok],
               stringsAsFactors = FALSE)
  }
  ivl <- rbind(side_ivl(hits[hits$flank == "upstream", ],   "upstream"),
               side_ivl(hits[hits$flank == "downstream", ], "downstream"))
  if (is.null(ivl) || !nrow(ivl)) return(empty)

  # merge overlapping intervals within each (gene, side) so bp are not double-counted
  key <- paste(ivl$og_idx, ivl$side, sep = "|")
  red <- reduce(split(IRanges(ivl$s, ivl$e), key))
  merged <- data.frame(
    key   = rep(names(red), lengths(red)),
    start = unlist(start(red)), end = unlist(end(red)),
    stringsAsFactors = FALSE)
  kp <- do.call(rbind, strsplit(merged$key, "|", fixed = TRUE))
  merged$og_idx <- as.integer(kp[, 1]); merged$side <- kp[, 2]

  # distribute each merged interval across the bins it spans; per-bin covered bp
  bin_lo <- (merged$start - 1L) %/% bin_bp + 1L
  bin_hi <- (merged$end   - 1L) %/% bin_bp + 1L
  reps   <- bin_hi - bin_lo + 1L
  ri     <- rep(seq_len(nrow(merged)), reps)
  bink   <- unlist(mapply(seq.int, bin_lo, bin_hi, SIMPLIFY = FALSE), use.names = FALSE)
  b_lo   <- (bink - 1L) * bin_bp + 1L
  b_hi   <-  bink       * bin_bp
  cov_bp <- pmin(merged$end[ri], b_hi) - pmax(merged$start[ri], b_lo) + 1L

  agg <- aggregate(cov_bp,
                   by = list(og_idx = merged$og_idx[ri], side = merged$side[ri], bink = bink),
                   FUN = sum)
  abs_mid <- (agg$bink - 0.5) * bin_bp
  data.frame(
    og_idx   = agg$og_idx,
    x        = ifelse(agg$side == "upstream", -abs_mid, abs_mid),
    cov_frac = pmin(1, agg$x / bin_bp),
    stringsAsFactors = FALSE)
}

# whole-genome TE class composition as % OF GENOME LENGTH.
# For every class, merge (reduce) its intervals per contig so bp are not
# double-counted, sum across contigs, and divide by total genome length. Mirrors
# the flank composition (Plot #3) but with the genome as the denominator, so the
# two are directly comparable. Returns a named vector of % for the requested
# classes (classes absent from the genome come back as 0).
genome_wide_composition <- function(te_gr, total_genome_len, classes) {
  key <- paste(as.character(seqnames(te_gr)), te_gr$te_class, sep = "\r")
  bp  <- reduced_bp_by_group(start(te_gr), end(te_gr), key)
  cl  <- sub("^.*\r", "", names(bp))
  per_class <- tapply(bp, cl, sum)
  out <- setNames(numeric(length(classes)), classes)
  common <- intersect(names(per_class), classes)
  out[common] <- per_class[common] / total_genome_len * 100
  round(out, 4)
}

# ---------------------------------------------------------------------
# Plots - TE density plot, percentage TE vs not TE plot, TE class proportions plot, flank completeness plot, top20 orthologs covered by TE plot
# ---------------------------------------------------------------------
make_plots <- function(master, hits, output_dir, species,
                       te_gr = NULL, total_genome_len = NA_real_,
                       genome_len_est = FALSE) {
  d_density <- file.path(output_dir, "plots", "density")
  d_comp    <- file.path(output_dir, "plots", "composition")
  d_summ    <- file.path(output_dir, "plots", "summary")
  for (d in c(d_density, d_comp, d_summ))
    dir.create(d, showWarnings = FALSE, recursive = TRUE)

  base_theme <- theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
  n_orth <- nrow(master)

  ## --- 1. TE coverage along flank, by density tier (EDIT 3) --------
  # Shows the MEAN fraction of each 500 bp bin covered by TE as one clean line PER TIER (low / medium / high),
  # plus a bold black overall-mean line. NOTE: we deliberately do NOT draw one
  # line per gene -- with 8000+ orthologs and near-binary per-gene coverage (a
  # 500 bp bin is usually either inside a TE ~1 or not ~0) that produced an
  # unreadable mesh. A small faint random sample of individual genes is
  # kept as background texture (set N_SAMPLE <- 0 to remove it). Also
  # writes the binned overall-mean profile (density_profile.tsv) that the
  # cross-species aggregation script consumes.
  pgc <- per_gene_bin_coverage(hits)              # og_idx, x, cov_frac  (in [0,1])

  # Tier each gene by FIXED cutoffs on pct_flank_te (the % of the gene's WHOLE
  # +/-10 kb flank that is TE) -- NOT tertiles. Tertiles force an uninformative
  # 33/33/33 split every time; fixed thresholds make the tier percentages a real,
  # comparable-across-species result. Rule: low < 20%, medium 20-50%, high > 50%.
  # Adjust the two cutoffs below if the biology calls for it.
  TIER_LOW_MAX  <- 20   # pct_flank_te below this  = low
  TIER_HIGH_MIN <- 50   # pct_flank_te above this  = high  (between = medium)
  master_tier <- factor(
    ifelse(master$pct_flank_te <  TIER_LOW_MAX,  "low",
    ifelse(master$pct_flank_te <= TIER_HIGH_MIN, "medium", "high")),
    levels = c("low", "medium", "high"))
  tier_n   <- table(master_tier)
  tier_pct <- round(100 * as.numeric(tier_n) / n_orth, 1)
  names(tier_pct) <- names(tier_n)
  tier_lab <- sprintf("%.1f%% low, %.1f%% medium, %.1f%% high",
                      tier_pct["low"], tier_pct["medium"], tier_pct["high"])

  if (nrow(pgc)) {
    pgc$tier <- master_tier[pgc$og_idx]
    # Overall mean coverage per bin ACROSS ALL genes (bold black line). A gene
    # absent from a bin contributes 0, so divide the summed coverage by n_orth.
    avg <- aggregate(cov_frac ~ x, data = pgc, FUN = sum)
    avg$mean_cov_frac <- avg$cov_frac / n_orth

    # Mean coverage per bin WITHIN each tier -> one clean line per tier. Divide
    # each tier's summed per-bin coverage by that tier's gene count, so genes in
    # the tier with no TE in a bin correctly count as 0.
    tier_avg <- aggregate(cov_frac ~ x + tier, data = pgc, FUN = sum)
    tier_avg$mean_cov_frac <-
      tier_avg$cov_frac / as.numeric(tier_n[as.character(tier_avg$tier)])
    tier_avg$tier <- factor(tier_avg$tier, levels = c("low", "medium", "high"))

    # Faint background sample of individual genes -- texture only, NOT the focus.
    # Small + very transparent so it can never become the old hairball.
    N_SAMPLE <- 150
    set.seed(1)                                   # reproducible sample
    og_ids <- unique(pgc$og_idx)
    samp   <- pgc[pgc$og_idx %in% sample(og_ids, min(N_SAMPLE, length(og_ids))), ]

    tier_cols <- c(low = "#4DBBD5", medium = "#F39B7F", high = "#E64B35")
    p1 <- ggplot() +
      geom_line(data = samp, aes(x = x, y = cov_frac, group = og_idx, colour = tier),
                alpha = 0.05, linewidth = 0.15) +                  # faint texture
      geom_line(data = tier_avg, aes(x = x, y = mean_cov_frac, colour = tier),
                linewidth = 1.2) +                                 # per-tier means
      geom_line(data = avg, aes(x = x, y = mean_cov_frac),
                colour = "black", linewidth = 1.3) +               # overall mean
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
      annotate("text", x = 0, y = Inf, label = "gene", vjust = 1.5, size = 3,
               colour = "grey30") +
      scale_colour_manual(values = tier_cols, name = "Density tier") +
      scale_x_continuous(breaks = seq(-FLANK_BP, FLANK_BP, by = 2500),
                         labels = function(b) paste0(abs(b) / 1000, "kb")) +
      guides(colour = guide_legend(override.aes = list(alpha = 1, linewidth = 1))) +
      labs(title = paste0("TE coverage along flanks by density tier — ", species),
           subtitle = paste0("Mean fraction of each ", BIN_BP,
                             "bp bin covered by TE per tier; bold black = overall mean (n = ",
                             n_orth, " orthologs).\nTiers (pct_flank_te: low <",
                             TIER_LOW_MAX, "%, medium ", TIER_LOW_MAX, "-", TIER_HIGH_MIN,
                             "%, high >", TIER_HIGH_MIN, "%): ", tier_lab),
           x = "Distance from gene  (<- upstream | downstream ->)",
           y = "Mean fraction of bin covered by TE") +
      base_theme +
      theme(plot.subtitle = element_text(size = 9))   # keep long subtitle from clipping
    ggsave(file.path(d_density, "te_density_per_gene_lines.pdf"),
           p1 + ylim(0, 1), width = 12, height = 6)

    # binned mean profile for the cross-species combined plot (EDIT 3)
    data.table::fwrite(
      data.frame(species = species, x = avg$x,
                 mean_cov_frac = round(avg$mean_cov_frac, 6)),
      file.path(output_dir, "density_profile.tsv"), sep = "\t")
  } else {
    message("  No TE coverage in flanks -- skipping per-gene density plot.")
  }

  # tier percentages companion file (X% low / Y% medium / Z% high)
  data.table::fwrite(
    data.frame(tier = names(tier_pct), n_genes = as.integer(tier_n),
               pct_genes = as.numeric(tier_pct)),
    file.path(output_dir, "density_tier_summary.tsv"), sep = "\t")

  ## --- 2. % TE vs non-TE (distribution across orthologs) -----------
  med <- median(master$pct_flank_te, na.rm = TRUE)
  p2 <- ggplot(master, aes(x = species, y = pct_flank_te)) +
    geom_violin(fill = "#4C9F70", alpha = 0.5, colour = NA) +
    geom_boxplot(width = 0.12, outlier.size = 0.4, fill = "white") +
    geom_hline(yintercept = med, linetype = "dashed", colour = "firebrick") +
    annotate("text", x = 0.6, y = med, label = sprintf("median = %.1f%%", med),
             vjust = -0.6, hjust = 0, size = 3.2, colour = "firebrick") +
    labs(title = paste0("Percent of flanking region that is TE — ", species),
         subtitle = paste0("Distribution across ", n_orth, " orthologs"),
         x = NULL, y = "% of flank covered by TE") +
    base_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  ggsave(file.path(d_comp, "pct_te_vs_nonTE.pdf"), p2 + ylim(0, 100), width = 9, height = 5)

  ## --- 3. TE class proportions (mean composition) ------------------
  # Derive the class set from the master's pct_<class> columns (preserving
  # their order) so this adapts to whatever classes are present.
  pct_cols <- setdiff(grep("^pct_", names(master), value = TRUE), "pct_flank_te")
  plot_classes <- sub("^pct_", "", pct_cols)
  comp <- data.frame(
    te_class = factor(plot_classes, levels = plot_classes),
    mean_pct = vapply(pct_cols, function(c) mean(master[[c]], na.rm = TRUE), numeric(1)))
  # Add a non-TE remainder so the stacked bar sums to exactly 100%.
  # (Classes are measured independently and may mutually overlap, so the
  # remainder absorbs that small excess rather than the bar exceeding 100.)
  non_te <- max(0, 100 - sum(comp$mean_pct))
  comp_plot <- rbind(comp,
                     data.frame(te_class = "non-TE", mean_pct = non_te))
  comp_plot$te_class <- factor(comp_plot$te_class,
                               levels = c(plot_classes, "non-TE"))
  pal <- te_class_palette(plot_classes)
  p3 <- ggplot(comp_plot, aes(x = species_name, y = mean_pct, fill = te_class)) +
    geom_col(width = 0.5) +
    scale_fill_manual(values = pal) +
    labs(title = paste0("Mean TE class composition of flanks — ", species),
         subtitle = paste0("Averaged across ", n_orth, " orthologs (per 20kb flank)"),
         x = NULL, y = "% of flanking region", fill = "TE class") +
    base_theme + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
  ggsave(file.path(d_comp, "te_class_proportions.pdf"), p3 + ylim(0, 100), width = 9, height = 5)

  ## --- 4. Flank completeness ---------------------------------------
  n_complete <- sum(master$flank_complete)
  comp_df <- data.frame(
    status = c("complete (20kb)", "truncated (contig edge)"),
    count  = c(n_complete, n_orth - n_complete))
  p4 <- ggplot(comp_df, aes(x = status, y = count, fill = status)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_text(aes(label = count), vjust = -0.3, size = 3.5) +
    scale_fill_manual(values = c("complete (20kb)" = "#4C9F70",
                                 "truncated (contig edge)" = "#D98880")) +
    labs(title = paste0("Flank completeness — ", species),
         subtitle = "Orthologs with a full 20kb flank vs truncated near a contig edge",
         x = NULL, y = "Number of orthologs") +
    base_theme
  ggsave(file.path(d_summ, "flank_completeness.pdf"), p4, width = 8, height = 5)

  ## --- 5. (extra) Top 20 orthologs by TE coverage ------------------
  top <- master %>%
    arrange(desc(pct_flank_te)) %>%
    head(20) %>%
    mutate(label = ifelse(nzchar(gene_name) & gene_name != ".",
                          paste0(gene_name, " (", orthogroup, ")"), orthogroup),
           label = factor(label, levels = rev(label)))
  p5 <- ggplot(top, aes(x = label, y = pct_flank_te)) +
    geom_col(fill = "#3C5488") +
    coord_flip() +
    labs(title = paste0("Top 20 orthologs by flank TE coverage — ", species),
         x = NULL, y = "% of flank covered by TE") +
    base_theme
  ggsave(file.path(d_comp, "top20_orthologs_by_te.pdf"),
         p5 + ylim(0, 100),
         width = 10, height = 6)

  ## --- 6. Genome-wide vs flank TE breakdown ---------------
  # Same-style composition, but computed for the WHOLE genome (% of total genome
  # length) placed next to the flank-only composition (reusing the Plot #3 means).
  # Both sides use the 100 bp-filtered TE set, so the comparison shows whether
  # flanks are enriched/depleted for particular TE classes vs the genome background.
  if (!is.null(te_gr)) {
    genome_classes <- sort(unique(te_gr$te_class))
    uni <- union(plot_classes, genome_classes)
    all_classes <- c(intersect(CLASS_ORDER, uni), setdiff(uni, CLASS_ORDER))

    gw <- genome_wide_composition(te_gr, total_genome_len, all_classes)
    fl <- setNames(numeric(length(all_classes)), all_classes)
    fl_present <- intersect(plot_classes, all_classes)
    fl[fl_present] <- comp$mean_pct[match(fl_present, comp$te_class)]

    brk <- rbind(
      data.frame(source = "whole genome", te_class = all_classes, pct = as.numeric(gw)),
      data.frame(source = "flanks",       te_class = all_classes, pct = as.numeric(fl)))
    # non-TE remainder per bar so each column sums to 100 (classes may mutually
    # overlap; remainder absorbs the small excess rather than exceeding 100)
    brk <- rbind(brk, data.frame(
      source   = c("whole genome", "flanks"), te_class = "non-TE",
      pct      = c(max(0, 100 - sum(gw)), max(0, 100 - sum(fl)))))
    brk$te_class <- factor(brk$te_class, levels = c(all_classes, "non-TE"))
    brk$source   <- factor(brk$source, levels = c("whole genome", "flanks"))

    est_note <- if (genome_len_est)
      "  [genome length ESTIMATED from max TE end/contig -- approximate]" else ""
    p6 <- ggplot(brk, aes(x = source, y = pct, fill = te_class)) +
      geom_col(width = 0.6) +
      scale_fill_manual(values = te_class_palette(all_classes)) +
      labs(title = paste0("Genome-wide vs flank TE composition — ", species),
           subtitle = paste0("% of bp that is each TE class (100 bp filter both sides)",
                             est_note),
           x = NULL, y = "% of bp", fill = "TE class") +
      base_theme
    ggsave(file.path(d_comp, "genome_vs_flank_breakdown.pdf"),
           p6 + ylim(0, 100), width = 8, height = 6)

    data.table::fwrite(
      data.frame(te_class = all_classes,
                 pct_genome = as.numeric(gw), pct_flank = as.numeric(fl)),
      file.path(output_dir, "genome_wide_composition.tsv"), sep = "\t")
  }

  message(sprintf("  Saved plots under %s", file.path(output_dir, "plots")))
}


# =====================================================================
main()
