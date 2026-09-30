#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------
# plot_coverage.R — Per-genotype coverage depth plots for ViroConstrictor output
#
# Usage:
#   Rscript plot_coverage.R [results_dir] [output_dir] [min_coverage]
#
#   results_dir   — pipeline output dir (default: the current directory).
#                   Searched recursively for
#                   consensus/coverages/<sample>_coverage.tsv files.
#   output_dir    — where PNGs are written (default: <results_dir>/coverage_plots)
#   min_coverage  — depth threshold below which the line is drawn red (default 20)
#
# Produces one PNG per genotype (IA, IB, IIIA, ...), x = genome position,
# y = depth on a log10 scale, one line per sample. Line segments where depth
# falls below <min_coverage> are coloured red.
# ---------------------------------------------------------------------------

suppressWarnings(suppressMessages({
  library(ggplot2)
  library(patchwork)
}))

# --- Arguments --------------------------------------------------------------
args        <- commandArgs(trailingOnly = TRUE)
results_dir <- ifelse(length(args) >= 1, args[1], ".")
output_dir  <- ifelse(length(args) >= 2, args[2], file.path(results_dir, "coverage_plots"))
min_cov     <- ifelse(length(args) >= 3, as.numeric(args[3]), 20)

if (!dir.exists(results_dir)) stop(sprintf("results_dir not found: %s", results_dir))
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# --- Locate per-position coverage TSVs --------------------------------------
cov_files <- list.files(
  results_dir,
  pattern    = "_coverage\\.tsv$",
  recursive  = TRUE,
  full.names = TRUE
)
# Keep only the per-position coverage files under consensus/coverages/
cov_files <- cov_files[grepl("consensus/coverages/", cov_files)]

if (length(cov_files) == 0) stop("No consensus/coverages/*_coverage.tsv files found.")

message(sprintf("Found %d coverage file(s).", length(cov_files)))

# --- Read + annotate with sample and genotype -------------------------------
# Genotype is the suffix of the RefID directory (e.g. RefID~AB279732.1_IIIA -> IIIA)
read_one <- function(path) {
  df <- read.table(path, sep = "\t", header = FALSE,
                   col.names = c("position", "depth"),
                   colClasses = c("integer", "numeric"))
  refid    <- sub(".*/RefID~([^/]+)/consensus/.*", "\\1", path)   # AB279732.1_IIIA
  genotype <- sub(".*_", "", refid)                               # IIIA
  sample   <- sub("_coverage\\.tsv$", "", basename(path))
  df$sample   <- sample
  df$genotype <- genotype
  df$refid    <- refid
  df
}

all_cov <- do.call(rbind, lapply(cov_files, read_one))

# --- Amplicon coordinates per genotype (from primer BED files) --------------
# For each genotype, use the genotype-matched primers (prefix HAV<genotype>_)
# and define each amplicon's span as min(start) .. max(end) of its primers.
bed_files <- list.files(
  results_dir, pattern = "_primers\\.bed$", recursive = TRUE, full.names = TRUE
)
bed_files <- bed_files[grepl("/primers/", bed_files) &
                       !grepl("removedprimers", bed_files)]

amplicons_for <- function(genotype) {
  # Pick one BED for this genotype (coords are reference-based, identical across samples)
  bed <- bed_files[grepl(sprintf("RefID~[^/]*_%s/", genotype), bed_files)][1]
  if (is.na(bed)) return(NULL)
  b <- read.table(bed, sep = "\t", header = FALSE,
                  col.names = c("ref", "start", "end", "name", "score", "strand"),
                  stringsAsFactors = FALSE)
  # Use ALL mapped primers (any genotype prefix) keyed by amplicon number. Cross-
  # genotype primers map to nearly the same positions, so they fill in boundaries
  # where the genotype-matched primer failed to map (e.g. missing HAVIB_1_LEFT).
  keep <- grepl("^HAV[A-Za-z0-9]+_[0-9]+_(LEFT|RIGHT)", b$name, ignore.case = TRUE)
  b <- b[keep, ]
  if (nrow(b) == 0) return(NULL)
  b$amplicon <- as.integer(sub("^HAV[A-Za-z0-9]+_([0-9]+)_.*", "\\1", b$name))
  agg <- do.call(rbind, lapply(split(b, b$amplicon), function(g) {
    data.frame(amplicon = g$amplicon[1],
               start    = min(g$start),
               end      = max(g$end))
  }))
  agg <- agg[order(agg$amplicon), ]
  # Stagger overlapping neighbours: odd amplicons on the top row, even on the bottom
  agg$row <- ifelse(agg$amplicon %% 2 == 1, 2L, 1L)
  agg$mid <- (agg$start + agg$end) / 2
  agg
}

# For a log10 y-axis, floor display depth at 1 (log10(0) is undefined).
# Threshold logic still uses the raw depth value.
all_cov$depth_disp <- pmax(all_cov$depth, 1)
# Red overlay: display value only where raw depth < threshold, NA elsewhere so
# geom_line breaks the red line outside the low-coverage regions.
all_cov$depth_low  <- ifelse(all_cov$depth < min_cov, pmax(all_cov$depth, 1), NA_real_)

genotypes <- sort(unique(all_cov$genotype))
message(sprintf("Genotypes: %s", paste(genotypes, collapse = ", ")))

# --- Shared amplicon track: ALL genotypes, each in its own scheme -----------
# The same three-genotype track is drawn under every coverage plot. Each
# genotype occupies a labelled band; within a band, odd amplicons sit on the
# upper sub-row and even on the lower, so overlapping neighbours stay visible.
gt_order   <- c("IA", "IB", "IIIA")                       # vertical order (top -> bottom)
gt_present <- c(gt_order[gt_order %in% genotypes],
                setdiff(sort(genotypes), gt_order))       # append any others
n_gt       <- length(gt_present)

track_list <- list()
for (i in seq_along(gt_present)) {
  g <- gt_present[i]
  a <- amplicons_for(g)
  if (is.null(a)) next
  a$genotype <- g
  band_base  <- (n_gt - i) * 2                            # top genotype = highest band
  a$y        <- band_base + ifelse(a$amplicon %% 2 == 1, 1, 0)  # odd above even
  a$band_mid <- band_base + 0.5
  track_list[[g]] <- a
}
track_df <- do.call(rbind, track_list)

# Per-genotype colours and y-axis band labels
base_pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")
gt_cols  <- setNames(base_pal[seq_len(n_gt)], gt_present)
band_labels <- tapply(track_df$band_mid, track_df$genotype, function(x) x[1])
band_labels <- band_labels[!is.na(band_labels)]

# Shared x range so the identical track aligns across every plot
global_xmax <- max(max(all_cov$position), max(track_df$end))

p_track <- ggplot(track_df) +
  geom_segment(aes(x = start, xend = end, y = y, yend = y, colour = genotype),
               linewidth = 4, lineend = "butt") +
  geom_text(aes(x = mid, y = y, label = amplicon),
            colour = "white", fontface = "bold", size = 2.7) +
  scale_colour_manual(values = gt_cols, guide = "none") +
  scale_y_continuous(breaks = as.numeric(band_labels), labels = names(band_labels)) +
  coord_cartesian(xlim = c(0, global_xmax),
                  ylim = c(min(track_df$y) - 0.8, max(track_df$y) + 0.8)) +
  labs(x = "Genome position (bp)", y = "Amplicon scheme") +
  theme_bw(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.ticks.y       = element_blank()
  )

# --- One plot per genotype (shared track underneath) ------------------------
for (gt in genotypes) {
  d <- all_cov[all_cov$genotype == gt, ]
  n_samples <- length(unique(d$sample))

  # --- Top panel: coverage depth ---
  p_cov <- ggplot(d, aes(x = position, group = sample)) +
    geom_line(aes(y = depth_disp), colour = "grey25", linewidth = 0.4) +
    geom_line(aes(y = depth_low), colour = "red", linewidth = 0.5) +
    geom_hline(yintercept = min_cov, linetype = "dashed",
               colour = "red", alpha = 0.6) +
    annotate("text", x = global_xmax, y = min_cov,
             label = sprintf("%gx", min_cov), hjust = 1.1, vjust = -0.4,
             colour = "red", size = 3) +
    scale_y_log10(
      breaks = c(1, 10, 100, 1000, 10000, 100000),
      labels = scales::comma
    ) +
    coord_cartesian(xlim = c(0, global_xmax)) +
    labs(
      title    = sprintf("HAV genotype %s — sequencing depth", gt),
      subtitle = sprintf("%d sample(s) · reference %s · red = depth < %gx · tracks: all genotype schemes",
                         n_samples, unique(d$refid)[1], min_cov),
      y        = "Depth (log10)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title       = element_text(face = "bold"),
      panel.grid.minor = element_blank(),
      axis.title.x     = element_blank(),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank()
    )

  # --- Combine with aligned x-axes ---
  p <- p_cov / p_track + plot_layout(heights = c(4, 2.4))

  outfile <- file.path(output_dir, sprintf("coverage_%s.png", gt))
  ggsave(outfile, p, width = 11, height = 7.5, dpi = 150)
  message(sprintf("  wrote %s (%d sample(s))", outfile, n_samples))
}

message("Done.")
