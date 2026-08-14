# =====================================================================
# Cichlid TE Flanking-Region — Cross-Species Aggregation
# ---------------------------------------------------------------------
# The main pipeline (V2_te_flanking_analysis.R) runs ONE species at a time
# and writes, per species: master_table.tsv and density_profile.tsv.
#
# This SEPARATE, self-contained script reads those already-written per-species
# files for all four species and produces the cross-species outputs:
#   - combined per-bin density line plot (one line per species + overall mean)
#   - orthogroup x species heatmap of pct_flank_te (absolute + row z-score,
#     rows hierarchically clustered; missing combinations shown grey)
#   - a cross-species SD table (median + SD across the 4 species) and the
#     top-10 most variable orthogroups
#
# RUN ORDER: run V2_te_flanking_analysis.R once for EACH of the four species
# first (so every output_dir below already contains master_table.tsv and
# density_profile.tsv), THEN run this script.
# =====================================================================

# ---- USER INPUTS: the four per-species output directories + names ----
# These must match the output_dir / species_name used in the main pipeline runs.
species_dirs <- c(
  Metriaclima_zebra_AF        = "insights/Metriaclima_zebra_AF_intersection_outputs",
  Pundamilia_nyererei_AF      = "insights/Pundamilia_nyererei_AF_intersection_outputs",
  Oreochromis_niloticus_AF    = "insights/Oreochromis_niloticus_AF_intersection_outputs",
  Neolamprologus_brichardi_AF = "insights/Neolamprologus_brichardi_AF_intersection_outputs")

# Where to write the combined outputs.
agg_output_dir <- "insights/cross_species_aggregation"

# which species' gene symbols to prefer when labelling
# orthogroups. Oreochromis is the best-annotated of the four (78% real symbols vs
# LOC placeholders); the others are used as a fallback where it has none.
GENE_NAME_PRIORITY <- "Oreochromis_niloticus_AF"
# ----------------------------------------------------------------------


main <- function() {
  
  # Loud banner (via cat) so it is obvious this is the aggregation stage.
  agg_banner("CROSS-SPECIES AGGREGATION")
  
  check_packages()
  suppressPackageStartupMessages({
    library(ggplot2)
    library(data.table)
  })
  
  species_names <- names(species_dirs)
  dir.create(agg_output_dir, showWarnings = FALSE, recursive = TRUE)
  
  ## ---- 0. Validate the per-species inputs are present -------------
  for (s in species_names) {
    mt <- file.path(species_dirs[[s]], "master_table.tsv")
    dp <- file.path(species_dirs[[s]], "density_profile.tsv")
    if (!file.exists(mt))
      stop("Missing master_table.tsv for ", s, " (expected ", mt,
           "). Run the main pipeline for this species first.", call. = FALSE)
    if (!file.exists(dp))
      stop("Missing density_profile.tsv for ", s, " (expected ", dp,
           "). Re-run the main pipeline for this species first.", call. = FALSE)
  }
  message("Found all per-species inputs for: ", paste(species_names, collapse = ", "))
  
  ## ---- 1. Combined density profile plot ------------------
  message("Building combined cross-species density plot...")
  combined_density_plot(species_dirs, agg_output_dir)
  
  ## ---- 2. Build the orthogroup x species matrix of pct_flank_te ---
  message("Assembling orthogroup x species matrix (pct_flank_te)...")
  mat <- build_pct_matrix(species_dirs)
  message(sprintf("  %s orthogroups x %s species",
                  format(nrow(mat), big.mark = ","), ncol(mat)))
  
  ## ---- 3. Heatmaps: absolute + row z-score ---------------
  message("Drawing clustered heatmaps (absolute + z-score)...")
  draw_heatmaps(mat, agg_output_dir)
  
  ## ---- 4. Cross-species SD table + top-10 ---------------
  message("Writing cross-species SD table and top-10 most variable...")
  gene_names <- build_gene_names(species_dirs)          # MODIFIED (MEETING 17)
  message(sprintf("  %s of %s orthogroups have a gene symbol",
                  format(sum(nzchar(gene_names)), big.mark = ","),
                  format(length(gene_names), big.mark = ",")))
  sd_table_and_top10(mat, agg_output_dir, gene_names)
  
  message("Done. Cross-species outputs written to: ", normalizePath(agg_output_dir))
}


# =====================================================================
# Helper functions
# =====================================================================

# Visually distinct banner, matching the main pipeline's species_banner style.
agg_banner <- function(txt) {
  bar <- paste(rep("=", 60), collapse = "")
  cat("\n", bar, "\n======== ", toupper(txt), " ========\n", bar, "\n\n", sep = "")
}

check_packages <- function() {
  needed <- c(ggplot2 = "CRAN", data.table = "CRAN")
  missing <- needed[!vapply(names(needed), requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing))
    stop("Missing required packages: ", paste(names(missing), collapse = ", "),
         "\nInstall with: install.packages(c('ggplot2','data.table'))", call. = FALSE)
}

# read each species' density_profile.tsv, overlay one coloured
# line per species (its mean per-bin TE coverage) plus one bold line = the mean
# across all species at each bin.
combined_density_plot <- function(species_dirs, out_dir) {
  prof <- data.table::rbindlist(lapply(names(species_dirs), function(s) {
    d <- data.table::fread(file.path(species_dirs[[s]], "density_profile.tsv"))
    d$species <- s
    d
  }))
  # overall mean across species at each x bin
  overall <- prof[, .(mean_cov_frac = mean(mean_cov_frac)), by = x]
  
  p <- ggplot(prof, aes(x = x, y = mean_cov_frac, colour = species)) +
    geom_line(linewidth = 0.8) +
    geom_line(data = overall, aes(x = x, y = mean_cov_frac),
              inherit.aes = FALSE, colour = "black", linewidth = 1.3) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    annotate("text", x = 0, y = Inf, label = "gene", vjust = 1.5, size = 3,
             colour = "grey30") +
    scale_x_continuous(breaks = seq(min(prof$x), max(prof$x), by = 2500),
                       labels = function(b) paste0(abs(b) / 1000, "kb")) +
    labs(title = "Mean TE coverage along flanks — all species",
         subtitle = "Coloured lines = per-species mean per-bin coverage; bold line = mean across species",
         x = "Distance from gene  (<- upstream | downstream ->)",
         y = "Mean fraction of bin covered by TE", colour = "Species") +
    theme_bw(base_size = 12) + theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(out_dir, "combined_density_all_species.pdf"),
         p + ylim(0, 1), width = 12, height = 6)   # fixed 0-1 to match the per-species density plots
}

# Read each master_table.tsv, keep orthogroup + pct_flank_te, and pivot to a
# matrix: rows = orthogroups, cols = species. Missing combinations are NA.
build_pct_matrix <- function(species_dirs) {
  long <- data.table::rbindlist(lapply(names(species_dirs), function(s) {
    d <- data.table::fread(file.path(species_dirs[[s]], "master_table.tsv"),
                           select = c("orthogroup", "pct_flank_te"))
    d$species <- s
    d
  }))
  wide <- as.data.frame(data.table::dcast(long, orthogroup ~ species,
                                          value.var = "pct_flank_te",
                                          fun.aggregate = mean))
  mat <- as.matrix(wide[, setdiff(names(wide), "orthogroup"), drop = FALSE])
  rownames(mat) <- wide$orthogroup
  # keep species column order stable = the order in species_dirs
  mat[, names(species_dirs), drop = FALSE]
}

# one display label per orthogroup, taken from the
# gene_name column ALREADY present in every master_table.tsv (no new input file).
# The four species can each carry a different symbol for the same orthogroup, so
# we pick one, in this preference order:
#   1. a real gene symbol (not a "LOC..." placeholder) from GENE_NAME_PRIORITY
#   2. a real gene symbol from any other species (species_dirs order)
#   3. any non-empty name, LOC placeholders included
#   4. "" -> callers fall back to the bare orthogroup ID
# Step 1 is pinned to a named species rather than "first in species_dirs" so that
# re-ordering species_dirs later does not silently change the labels.
build_gene_names <- function(species_dirs) {
  nm <- lapply(names(species_dirs), function(s) {
    d <- data.table::fread(file.path(species_dirs[[s]], "master_table.tsv"),
                           select = c("orthogroup", "gene_name"))
    d <- d[!duplicated(d$orthogroup), ]
    setNames(as.character(d$gene_name), as.character(d$orthogroup))
  })
  names(nm) <- names(species_dirs)
  ord <- c(intersect(GENE_NAME_PRIORITY, names(nm)),
           setdiff(names(nm), GENE_NAME_PRIORITY))
  nm  <- nm[ord]

  ok   <- function(x) !is.na(x) & nzchar(x) & x != "."
  real <- function(x) ok(x) & !grepl("^LOC", x)

  ogs  <- unique(unlist(lapply(nm, names), use.names = FALSE))
  cols <- vapply(nm, function(tab) unname(tab[ogs]), character(length(ogs)))
  dim(cols) <- c(length(ogs), length(nm))

  out <- rep("", length(ogs))
  for (j in seq_len(ncol(cols))) {          # pass 1: real symbols, priority first
    v <- cols[, j]; take <- !nzchar(out) & real(v)
    out[take] <- v[take]
  }
  for (j in seq_len(ncol(cols))) {          # pass 2: LOC placeholders as fallback
    v <- cols[, j]; take <- !nzchar(out) & ok(v)
    out[take] <- v[take]
  }
  setNames(out, ogs)
}

# "OG0001234" -> "tvp23b (OG0001234)" when a name is known,
# else the bare ID. Matches the main pipeline's top-20 plot labelling.
og_label <- function(og, gene_names) {
  og <- as.character(og)
  g  <- unname(gene_names[og])
  g[is.na(g)] <- ""
  ifelse(nzchar(g), paste0(g, " (", og, ")"), og)
}

# two heatmaps sharing the same biological species column order (Oreo, Metria,
# Neo, Punda). (a) ABSOLUTE: rows clustered hierarchically (NA-safe via row-mean
# imputation used ONLY for the clustering distance) with the clustering tree drawn
# on the row side. (b) Z-SCORE: rows ordered by Oreochromis TE% (high->low), no
# tree. Missing cells are grey. Colours unchanged (magma).
draw_heatmaps <- function(mat, out_dir) {
  # ---- reorder columns into biological species order: Oreo, Metria, Neo, Punda ----
  BIO_ORDER <- c("Oreochromis_niloticus_AF", "Metriaclima_zebra_AF",
                 "Neolamprologus_brichardi_AF", "Pundamilia_nyererei_AF")
  mat <- mat[, c(intersect(BIO_ORDER, colnames(mat)),
                 setdiff(colnames(mat), BIO_ORDER)), drop = FALSE]
  sp_levels <- colnames(mat)
  n <- nrow(mat); k <- ncol(mat)

  # ---- hierarchical row clustering (impute NA just for the distance calc) ----
  mat_imp <- mat
  rmean <- rowMeans(mat_imp, na.rm = TRUE)
  na_idx <- which(is.na(mat_imp), arr.ind = TRUE)
  if (nrow(na_idx)) mat_imp[na_idx] <- rmean[na_idx[, 1]]
  mat_imp[is.na(mat_imp)] <- mean(mat_imp, na.rm = TRUE)  # any all-NA row
  hc  <- hclust(dist(mat_imp))
  ord <- rownames(mat)[hc$order]

  # ---- hclust -> elbow dendrogram segments in (leaf-position, height) space ----
  dendro_segments <- function(h) {
    m   <- length(h$order)
    pos <- numeric(m); pos[h$order] <- seq_len(m)
    node_x <- numeric(nrow(h$merge))
    getpos <- function(kk) if (kk < 0) pos[-kk] else node_x[kk]
    geth   <- function(kk) if (kk < 0) 0       else h$height[kk]
    out <- vector("list", nrow(h$merge))
    for (i in seq_len(nrow(h$merge))) {
      a <- h$merge[i, 1]; b <- h$merge[i, 2]; ht <- h$height[i]
      xa <- getpos(a); xb <- getpos(b)
      node_x[i] <- (xa + xb) / 2
      out[[i]] <- data.frame(
        x    = c(xa, xb, xa),   y    = c(geth(a), geth(b), ht),
        xend = c(xa, xb, xb),   yend = c(ht,      ht,      ht))
    }
    do.call(rbind, out)
  }

  # ---- long form for ggplot: numeric row (ypos) + species (xpos) positions ----
  mk_long <- function(m, value_name, row_order) {
    df <- as.data.frame(as.table(m), stringsAsFactors = FALSE)
    names(df) <- c("orthogroup", "species", value_name)
    df$ypos <- match(df$orthogroup, row_order)
    df$xpos <- match(df$species, sp_levels)
    df
  }
  tile_theme <- theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1),
          plot.margin = margin(6, 6, 14, 26))   # extra bottom/left so rotated species labels aren't clipped

  # ---- (a) absolute: raw pct_flank_te, capped colour scale, with row dendrogram ----
  abs_df <- mk_long(mat, "pct", ord)
  cap <- as.numeric(quantile(abs_df$pct, 0.98, na.rm = TRUE))  # cap so extremes don't wash out mid-range
  abs_df$pct_cap <- pmin(abs_df$pct, cap)

  # map the clustering tree into the empty strip left of the tiles (heights -> x, leaves -> y)
  seg    <- dendro_segments(hc)
  maxh   <- max(seg$y)
  leaf_x <- 0.35; DW <- 1.6                     # leaf edge just left of the tiles; tree width
  tx     <- function(hh) leaf_x - (hh / maxh) * DW
  segp   <- data.frame(x = tx(seg$y), y = seg$x, xend = tx(seg$yend), yend = seg$xend)
  root_x <- tx(maxh)

  p_abs <- ggplot(abs_df, aes(xpos, ypos, fill = pct_cap)) +
    geom_tile() +
    geom_segment(data = segp, aes(x = x, y = y, xend = xend, yend = yend),
                 inherit.aes = FALSE, linewidth = 0.15, colour = "grey30") +
    scale_fill_viridis_c(option = "magma", na.value = "grey75",
                         name = sprintf("pct_flank_te\n(capped @ %.0f)", cap)) +
    scale_x_continuous(breaks = seq_len(k), labels = sp_levels,
                       limits = c(root_x - 0.1, k + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    labs(title = "Flank TE coverage across species (absolute)",
         subtitle = paste0(nrow(mat), " orthogroups, rows hierarchically clustered; grey = missing"),
         x = NULL, y = "Orthogroups (clustered)") +
    tile_theme
  ggsave(file.path(out_dir, "heatmap_absolute.pdf"), p_abs, width = 7, height = 10)

  # ---- (b) row z-score: rows ordered by Oreochromis TE% (high->low), no dendrogram ----
  oreo  <- mat[, "Oreochromis_niloticus_AF"]
  ord_z <- rownames(mat)[order(oreo, decreasing = FALSE, na.last = FALSE)]  # highest at top
  rm2  <- rowMeans(mat, na.rm = TRUE)
  rsd  <- apply(mat, 1, sd, na.rm = TRUE)
  z    <- (mat - rm2) / rsd            # rm2/rsd recycle down columns -> per-row
  z[!is.finite(z)] <- NA               # rows with SD 0 (or a single value) -> NA
  z_df <- mk_long(z, "z", ord_z)
  zmax <- max(abs(z_df$z), na.rm = TRUE)
  p_z <- ggplot(z_df, aes(xpos, ypos, fill = z)) +
    geom_tile() +
    scale_fill_viridis_c(option = "magma", limits = c(-zmax, zmax),
                         na.value = "grey75", name = "row z-score") +
    scale_x_continuous(breaks = seq_len(k), labels = sp_levels,
                       limits = c(0.5, k + 0.5), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.5, n + 0.5), expand = c(0, 0)) +
    labs(title = "Flank TE coverage across species (row z-score)",
         subtitle = "Which species stands out WITHIN each orthogroup; rows ordered by Oreochromis TE% (high to low); grey = missing / SD 0",
         x = NULL, y = "Orthogroups (ordered by Oreochromis TE%)") +
    tile_theme
  ggsave(file.path(out_dir, "heatmap_zscore.pdf"), p_z, width = 7, height = 10)
}

# SD table (orthogroups x species pct_flank_te + median, mean, SD across
# the 4 species), ranked by SD; write it out and pull the top-10 most variable.
sd_table_and_top10 <- function(mat, out_dir, gene_names) {          
  sd_tbl <- data.frame(orthogroup = rownames(mat),
                       gene_name  = unname(gene_names[rownames(mat)]),  
                       mat, check.names = FALSE,
                       stringsAsFactors = FALSE)
  sd_tbl$median <- apply(mat, 1, median, na.rm = TRUE)
  sd_tbl$mean   <- rowMeans(mat, na.rm = TRUE)              
  sd_tbl$SD     <- apply(mat, 1, sd, na.rm = TRUE)
  sd_tbl <- sd_tbl[order(-sd_tbl$SD), ]
  
  data.table::fwrite(sd_tbl, file.path(out_dir, "cross_species_sd_table.tsv"), sep = "\t")
  
  top10 <- utils::head(sd_tbl[!is.na(sd_tbl$SD), ], 10)
  data.table::fwrite(top10, file.path(out_dir, "cross_species_sd_top10.tsv"), sep = "\t")
  
  # small focused figure of just the top-10 (separate from the main heatmap)
  long <- data.table::melt(data.table::as.data.table(top10),
                           id.vars = "orthogroup",
                           measure.vars = colnames(mat),
                           variable.name = "species", value.name = "pct")
  # MODIFIED (MEETING 17): label rows "gene (OG…)" instead of the bare ID
  long$label <- og_label(long$orthogroup, gene_names)
  long$label <- factor(long$label,
                       levels = rev(og_label(top10$orthogroup, gene_names)))
  # white text on dark (low-value) cells, black on light cells, so labels stay readable
  rng <- range(long$pct, na.rm = TRUE)
  long$txtcol <- ifelse(is.na(long$pct), "black",
                        ifelse((long$pct - rng[1]) / (rng[2] - rng[1]) < 0.55, "white", "black"))
  p <- ggplot(long, aes(species, label, fill = pct)) +
    geom_tile(colour = "grey90") +
    geom_text(aes(label = ifelse(is.na(pct), "NA", sprintf("%.1f", pct)), colour = txtcol),
              size = 3) +
    scale_colour_identity() +
    scale_fill_viridis_c(option = "magma", na.value = "grey75", name = "pct_flank_te") +
    labs(title = "Top 10 most cross-species-variable orthogroups (highest SD)",
         subtitle = "Flank TE coverage per species; ranked by SD across the 4 species",
         x = NULL, y = NULL) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1),
          plot.margin = margin(6, 6, 14, 20))   # extra bottom/left so 'Metriaclima...' isn't clipped
  ggsave(file.path(out_dir, "top10_variable_orthogroups.pdf"), p, width = 8, height = 6)
}


# =====================================================================
main()
