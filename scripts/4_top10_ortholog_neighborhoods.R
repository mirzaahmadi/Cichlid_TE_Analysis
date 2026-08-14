# =====================================================================
# Top-10 Cross-Species-Variable Orthologs — Genomic Neighborhood Plots
# ---------------------------------------------------------------------
# For the 10 single-copy orthologs whose flank TE coverage varies most
# across the four African cichlids (highest SD, from
# 6/insights/V3_cross_species_aggregation.R's cross_species_sd_top10.tsv),
# draws the +/-10 kb genomic neighborhood of the ortholog gene in all four
# species side by side, with every TE in the window drawn as a rectangle
# coloured by TE class.
#
# Run from inside the "8" folder.
# =====================================================================

# ---- USER INPUTS ----
top10_file     <- "data/cross_species_sd_top10.tsv"
gene_info_file <- "data/OUTPUT_TABLE_orthogroup_gene_info.tsv"

# One RepeatMasker .out per species. Named vector order = display order
# requested for the figures (Oreochromis, Metriaclima, Neolamprologus,
# Pundamilia); SPECIES_ORDER below is just names(te_files).
te_files <- c(
  Oreochromis_niloticus_AF    = "data/TE_Oreochromis_niloticus_AF_genome.fa.out",
  Metriaclima_zebra_AF        = "data/TE_Metriaclima_zebra_AF_genome.fa.out",
  Neolamprologus_brichardi_AF = "data/TE_Neolamprologus_brichardi_AF_genome.fa.out",
  Pundamilia_nyererei_AF      = "data/TE_Pundamilia_nyererei_AF_genome.fa.out")
SPECIES_ORDER <- names(te_files)

FLANK_BP   <- 10000   # window size per side of the gene, matches V3
output_dir <- "insights/top10_ortholog_neighborhoods"
# ----------------------


main <- function() {

  banner("TOP-10 ORTHOLOG NEIGHBORHOOD PLOTS")
  check_packages()
  suppressPackageStartupMessages({
    library(GenomicRanges)
    library(IRanges)
    library(data.table)
    library(dplyr)
    library(ggplot2)
  })

  for (f in c(top10_file, gene_info_file, te_files))
    if (!file.exists(f)) stop("Required input file not found: ", f, call. = FALSE)
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  ## ---- 1. Top-10 table (already ranked most->least variable) ------
  top10 <- data.table::fread(top10_file, data.table = FALSE)
  message(sprintf("Loaded top-10 table: %d orthogroups", nrow(top10)))

  ## ---- 2. Gene coordinates for these 10 orthogroups, our 4 species -
  gene_info_all <- data.table::fread(gene_info_file, data.table = FALSE)
  gene_info <- gene_info_all %>%
    dplyr::filter(orthogroup %in% top10$orthogroup, species %in% SPECIES_ORDER)
  # sanity check: every orthogroup should have exactly one gene row per species
  n_per_og <- table(gene_info$orthogroup)
  if (any(n_per_og != length(SPECIES_ORDER)))
    warning("Some top-10 orthogroups are missing a gene row for one or more species; ",
            "their neighborhood plot will simply omit that species' track.")

  ## ---- 3. Parse + filter each species' TE annotation (once each) --
  # Same parse_te_file() and >=100bp filter as V3_te_flanking_analysis.R,
  # so the TE set here matches what the rest of the pipeline analyzed.
  te_gr_list  <- list()
  contig_caps <- list()   
  for (s in SPECIES_ORDER) {
    message("Parsing TE annotation for ", s, "...")
    te_gr <- parse_te_file(te_files[[s]])
    te_gr <- te_gr[width(te_gr) >= 100L]              # sub-100bp fragment filter
    message(sprintf("  %s TEs retained (>=100bp)", format(length(te_gr), big.mark = ",")))
    te_gr_list[[s]]  <- te_gr
    contig_caps[[s]] <- tapply(end(te_gr), as.character(seqnames(te_gr)), max)
  }

  ## ---- 4. Reconstruct per-TE positions for the 10x4 gene windows --
  message("Building +/-", FLANK_BP, "bp windows and intersecting against TEs...")
  te_table   <- build_te_table(top10$orthogroup, gene_info, te_gr_list, contig_caps)
  gene_table <- build_gene_table(top10$orthogroup, gene_info)
  message(sprintf("  %s TE fragments found across all windows", format(nrow(te_table), big.mark = ",")))

  ## ---- 5. Write the reconstructed per-TE table ---------------------
  te_out_tsv <- file.path(output_dir, "top10_neighborhood_TEs.tsv")
  data.table::fwrite(te_table, te_out_tsv, sep = "\t")
  message("Wrote reconstructed per-TE table: ", te_out_tsv)

  ## ---- 6. Shared TE-class colour palette (consistent w/ rest of pipeline) --
  all_classes     <- sort(unique(te_table$te_class))
  ordered_classes <- c(intersect(CLASS_ORDER, all_classes), setdiff(all_classes, CLASS_ORDER))
  pal             <- te_class_palette(ordered_classes)

  ## ---- 7. One plot per orthogroup + a combined multi-page PDF -----
  message("Drawing neighborhood plots...")
  plots <- vector("list", nrow(top10))
  for (i in seq_len(nrow(top10))) {
    og <- top10$orthogroup[i]
    p  <- make_neighborhood_plot(og, top10[i, ], te_table, gene_table, ordered_classes, pal)
    plots[[i]] <- p
    fname <- sprintf("top10_neighborhood_%02d_%s.pdf", i, og)
    ggsave(file.path(output_dir, fname), p, width = 10, height = 8)
  }

  combined_pdf <- file.path(output_dir, "top10_neighborhoods_combined.pdf")
  pdf(combined_pdf, width = 10, height = 8)
  for (p in plots) print(p)
  invisible(dev.off())
  message("Wrote combined multi-page PDF: ", combined_pdf)

  message("Done. Output written to: ", normalizePath(output_dir))
}


# =====================================================================
# Helpers copied verbatim (logic-identical) from V3_te_flanking_analysis.R
# so TE parsing / filtering / colours stay consistent with the rest of the
# pipeline. Not sourced from that file, per the "standalone script" spec.
# =====================================================================

CLASS_ORDER <- c("DNA", "LTR", "LINE", "SINE", "RC", "Unknown", "Simple_repeat",
                  "Low_complexity", "tRNA", "rRNA", "snRNA")

# Parse a RepeatMasker .out file into a GRanges. 2 header lines + 1 blank
# line are skipped; RepeatMasker strand 'C' means minus.
parse_te_file <- function(path) {
  d <- data.table::fread(
    path, skip = 3, header = FALSE, sep = " ", quote = "",
    fill = TRUE, strip.white = TRUE, blank.lines.skip = TRUE,
    select = c(5, 6, 7, 9, 10, 11),
    col.names = c("chrom", "start", "end", "strand", "repeat_name", "class_family"))

  strand   <- ifelse(d$strand == "C", "-", "+")
  te_class <- sub("/.*", "", d$class_family)   # text before first '/', else whole string

  GRanges(
    seqnames = d$chrom,
    ranges   = IRanges(start = d$start, end = d$end),
    strand   = strand,
    repeat_name = d$repeat_name,
    te_family   = d$class_family,
    te_class    = te_class)
}

# Stable colour map: fixed colours for known classes, auto-generated
# qualitative palette for any extras, grey for the non-TE remainder.
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

banner <- function(txt) {
  bar <- paste(rep("=", 60), collapse = "")
  cat("\n", bar, "\n======== ", toupper(txt), " ========\n", bar, "\n\n", sep = "")
}

check_packages <- function() {
  needed <- list(GenomicRanges = "Bioconductor", IRanges = "Bioconductor",
                  data.table = "CRAN", dplyr = "CRAN", ggplot2 = "CRAN")
  missing <- needed[!vapply(names(needed), requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    msg <- paste(sprintf("  - %s (install from %s)", names(missing), unlist(missing)), collapse = "\n")
    stop("Missing required packages:\n", msg,
         "\nBioconductor: install.packages('BiocManager'); BiocManager::install(c('GenomicRanges','IRanges'))",
         "\nCRAN: install.packages(c('data.table','dplyr','ggplot2'))", call. = FALSE)
  }
}


# =====================================================================
# New logic for this script: coordinate re-centering, window
# reconstruction, and plotting.
# =====================================================================

# Map one genomic coordinate to a gene-relative, orientation-corrected
# coordinate: 0 = gene start (TSS), positive = downstream in the direction
# of transcription, negative = upstream. For '+' genes this is just
# (pos - gene_start); for '-' genes we flip around gene_end so that every
# panel reads 5'->3' left-to-right regardless of genomic strand.
# NOTE: `pos` is typically many TE coordinates for ONE gene (strand a
# scalar), while `strand`/`gene_start`/`gene_end` are the vectorized case in
# build_gene_table(). ifelse()'s output length follows its *condition* arg,
# not `yes`/`no` -- so a scalar `strand` against a long `pos` would silently
# truncate to length 1. rep_len() everything to a common length first avoids
# that trap in either calling pattern.
orient_pos <- function(pos, gene_start, gene_end, strand) {
  n <- max(length(pos), length(strand))
  pos        <- rep_len(pos, n)
  gene_start <- rep_len(gene_start, n)
  gene_end   <- rep_len(gene_end, n)
  strand     <- rep_len(strand, n)
  ifelse(strand == "-", gene_end - pos, pos - gene_start)
}

# Same idea for an interval [s, e]: transform both ends, then re-sort them
# (the '-' strand transform reverses order) so relative_start <= relative_end.
orient_interval <- function(s, e, gene_start, gene_end, strand) {
  a <- orient_pos(s, gene_start, gene_end, strand)
  b <- orient_pos(e, gene_start, gene_end, strand)
  data.frame(relative_start = pmin(a, b), relative_end = pmax(a, b))
}

# For each orthogroup x species: build the [-FLANK_BP, +FLANK_BP] window
# around the gene, intersect it against that species' (100bp-filtered) TE
# GRanges, clip hits to the window, and convert to oriented relative
# coordinates. Returns one row per (orthogroup, species, TE fragment).
build_te_table <- function(orthogroups, gene_info, te_gr_list, contig_caps) {
  rows <- vector("list", 0L)
  for (og in orthogroups) {
    for (s in names(te_gr_list)) {
      gi <- gene_info[gene_info$orthogroup == og & gene_info$species == s, ]
      if (!nrow(gi)) next               # missing species for this orthogroup -- skip
      gi <- gi[1, ]
      te_gr <- te_gr_list[[s]]

      # Window = gene body +/- FLANK_BP, clamped to [1, contig end]. No GFF3
      # is used, so (as in V3's no-GFF3 fallback) the upper cap is the max TE
      # coordinate seen on that contig; if the contig has no TEs at all the
      # cap is simply skipped (there's nothing there to clip against anyway).
      cap <- contig_caps[[s]][gi$chromosome]
      win_start <- max(1L, gi$start - FLANK_BP)
      win_end   <- gi$end + FLANK_BP
      if (!is.na(cap)) win_end <- min(win_end, as.integer(cap))
      if (win_end < win_start) win_end <- gi$end + FLANK_BP   # safety fallback

      win_gr <- GRanges(seqnames = gi$chromosome, ranges = IRanges(win_start, win_end))
      ov <- findOverlaps(win_gr, te_gr, ignore.strand = TRUE)
      if (!length(ov)) next
      sh <- subjectHits(ov)

      clip_s <- pmax(win_start, start(te_gr)[sh])
      clip_e <- pmin(win_end,   end(te_gr)[sh])
      rel    <- orient_interval(clip_s, clip_e, gi$start, gi$end, gi$strand)

      # Flip TE strand symbol too when the gene is '-', so the drawn TE
      # orientation stays meaningful in the re-oriented (5'->3') panel.
      te_strand  <- as.character(strand(te_gr)[sh])
      disp_strand <- if (gi$strand == "-") c(`+` = "-", `-` = "+", `*` = "*")[te_strand] else te_strand

      rows[[length(rows) + 1]] <- data.frame(
        orthogroup = og, species = s,
        te_class   = te_gr$te_class[sh],
        relative_start = rel$relative_start, relative_end = rel$relative_end,
        strand = unname(disp_strand),
        repeat_name = te_gr$repeat_name[sh],   # kept for optional labeling
        te_family   = te_gr$te_family[sh],     # kept for optional labeling
        stringsAsFactors = FALSE)
    }
  }
  dplyr::bind_rows(rows)
}

# Gene-body rows for plotting: same orient_interval() transform applied to
# the gene's own start/end, so the gene box uses exactly the same coordinate
# system as the TEs (always comes out as [0, gene_length] by construction).
build_gene_table <- function(orthogroups, gene_info) {
  gi  <- gene_info[gene_info$orthogroup %in% orthogroups, ]
  rel <- orient_interval(gi$start, gi$end, gi$start, gi$end, gi$strand)
  data.frame(orthogroup = gi$orthogroup, species = gi$species,
             gene_name = gi$gene_name, strand = gi$strand,
             relative_start = rel$relative_start, relative_end = rel$relative_end,
             stringsAsFactors = FALSE)
}

# One figure per orthogroup: 4 stacked species tracks (facet), gene as a
# thick arrowed box above the baseline, TEs as coloured boxes below it.
make_neighborhood_plot <- function(og, top10_row, te_table, gene_table, ordered_classes, pal) {
  te_sub   <- te_table[te_table$orthogroup == og, ]
  gene_sub <- gene_table[gene_table$orthogroup == og, ]
  te_sub$species   <- factor(te_sub$species,   levels = SPECIES_ORDER)
  gene_sub$species <- factor(gene_sub$species, levels = SPECIES_ORDER)
  te_sub$te_class  <- factor(te_sub$te_class,  levels = ordered_classes)

  gene_name  <- top10_row$gene_name[1]
  title_lab  <- if (nzchar(gene_name) && gene_name != ".") paste0(gene_name, " (", og, ")") else og
  sd_val     <- top10_row$SD[1]

  # Small arrowhead near the 3' end of each gene box, indicating the
  # (already re-oriented, so always left-to-right) direction of transcription.
  arrow_len   <- pmin(300, pmax(1, (gene_sub$relative_end - gene_sub$relative_start) * 0.15))
  gene_sub$arrow_x0 <- gene_sub$relative_end - arrow_len
  gene_sub$arrow_x1 <- gene_sub$relative_end

  p <- ggplot() +
    geom_hline(yintercept = 0, colour = "grey75", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_rect(data = te_sub,
              aes(xmin = relative_start, xmax = relative_end, ymin = -0.55, ymax = -0.15,
                  fill = te_class), colour = NA, alpha = 0.9) +
    geom_rect(data = gene_sub,
              aes(xmin = relative_start, xmax = relative_end, ymin = 0.15, ymax = 0.55),
              fill = "grey20", colour = "black", linewidth = 0.3) +
    geom_segment(data = gene_sub,
                 aes(x = arrow_x0, xend = arrow_x1, y = 0.35, yend = 0.35),
                 colour = "white", linewidth = 1,
                 arrow = arrow(length = unit(0.08, "inches"), type = "closed")) +
    geom_text(data = gene_sub,
              aes(x = (relative_start + relative_end) / 2, y = 0.72, label = gene_name),
              size = 2.8, fontface = "bold") +
    facet_wrap(vars(species), ncol = 1, strip.position = "right") +
    scale_fill_manual(values = pal, name = "TE class") +
    scale_x_continuous(name = "Distance from gene TSS (kb)",
                        labels = function(b) b / 1000) +
    coord_cartesian(ylim = c(-0.75, 0.9)) +
    labs(title = title_lab,
         subtitle = sprintf(
           "Cross-species SD of %% flank TE coverage = %.2f  |  +/-%dkb window, oriented 5'->3' (arrow = gene direction)",
           sd_val, FLANK_BP / 1000),
         y = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          plot.subtitle = element_text(size = 9),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          strip.text.y.right = element_text(angle = 0, hjust = 0))
  p
}


# =====================================================================
main()
