#!/usr/bin/env Rscript
# Diagnose tree directory matching issues

library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  cat("Usage: Rscript diagnose_trees.R <out_base> <batch_fa> <nc_lineages.tsv>\n")
  quit(status = 1)
}

out_base <- args[1]
batch_fa <- args[2]
nc_file <- args[3]

trees_dir <- file.path(out_base, "trees")

cat("=== DIAGNOSTIC: Tree Directory Matching ===\n\n")

# 1. Check FASTA headers
cat("1. FASTA HEADERS from batch_fa:\n")
cat("   File: ", batch_fa, "\n", sep = "")

if (!file.exists(batch_fa)) {
  cat("   ERROR: File not found!\n\n")
} else {
  lines <- readLines(batch_fa)
  headers <- lines[grep("^>", lines)]
  headers <- sub("^>", "", headers)
  
  cat("   Found ", length(headers), " sequences:\n", sep = "")
  for (i in seq_along(headers[1:min(5, length(headers))])) {
    normalized <- word(headers[i], 1)
    cat("     [", i, "] Full: '", headers[i], "'\n", sep = "")
    cat("         Normalized: '", normalized, "'\n", sep = "")
  }
  if (length(headers) > 5) cat("     ... (", length(headers) - 5, " more)\n", sep = "")
  cat("\n")
}

# 2. Check NextClade seqNames
cat("2. NEXTCLADE SEQUENCE NAMES:\n")
cat("   File: ", nc_file, "\n", sep = "")

if (!file.exists(nc_file)) {
  cat("   ERROR: File not found!\n\n")
} else {
  nc <- read_tsv(nc_file, show_col_types = FALSE)
  if ("seqName" %in% names(nc)) {
    seq_names <- nc$seqName
    cat("   Found ", length(seq_names), " sequences:\n", sep = "")
    for (i in seq_along(seq_names[1:min(5, length(seq_names))])) {
      normalized <- word(seq_names[i], 1)
      cat("     [", i, "] seqName: '", seq_names[i], "'\n", sep = "")
      cat("         Normalized: '", normalized, "'\n", sep = "")
    }
    if (length(seq_names) > 5) cat("     ... (", length(seq_names) - 5, " more)\n", sep = "")
    cat("\n")
  } else {
    cat("   ERROR: No 'seqName' column in NextClade output!\n\n")
  }
}

# 3. Check tree directories
cat("3. ACTUAL TREE DIRECTORIES:\n")
cat("   Path: ", trees_dir, "\n", sep = "")

if (!dir.exists(trees_dir)) {
  cat("   ERROR: trees_dir not found!\n\n")
} else {
  tree_dirs <- list.dirs(trees_dir, recursive = FALSE, full.names = FALSE)
  cat("   Found ", length(tree_dirs), " directories:\n", sep = "")
  for (i in seq_along(tree_dirs[1:min(5, length(tree_dirs))])) {
    treefile_exists <- file.exists(file.path(trees_dir, tree_dirs[i], "tree.treefile"))
    status <- if_else(treefile_exists, "✓", "✗")
    cat("     [", i, "] '", tree_dirs[i], "' ", status, "\n", sep = "")
  }
  if (length(tree_dirs) > 5) cat("     ... (", length(tree_dirs) - 5, " more)\n", sep = "")
  cat("\n")
}

# 4. Compare: which names are present in each source?
if (length(headers) > 0 && length(seq_names) > 0 && length(tree_dirs) > 0) {
  headers_norm <- word(headers, 1)
  seq_names_norm <- word(seq_names, 1)
  
  cat("4. NAME MATCHING:\n\n")
  
  # Check each NextClade name
  match_results <- tibble(
    seq_name = seq_names,
    seq_norm = word(seq_names, 1),
    in_fasta = seq_norm %in% headers_norm,
    dir_exists = seq_norm %in% tree_dirs,
    treefile_exists = file.exists(file.path(trees_dir, seq_norm, "tree.treefile"))
  )
  
  print(match_results, n = Inf)
  
  cat("\n=== SUMMARY ===\n")
  cat("Sequences in FASTA: ", length(headers), "\n")
  cat("Sequences in NextClade: ", length(seq_names), "\n")
  cat("Tree directories created: ", length(tree_dirs), "\n")
  cat("Matches (FASTA ↔ NextClade): ", sum(match_results$in_fasta), " / ", nrow(match_results), "\n")
  cat("Tree dirs found: ", sum(match_results$dir_exists), " / ", nrow(match_results), "\n")
  cat("Treefiles present: ", sum(match_results$treefile_exists), " / ", nrow(match_results), "\n")
}
