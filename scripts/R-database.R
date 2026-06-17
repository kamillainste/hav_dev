library(DBI)
library(RSQLite)
library(tidyverse)
library(seqinr)

# --- CONNECT ---
db_file <- "data/local_datasets/eksterne_sek/eksterne.db"
con <- dbConnect(SQLite(), dbname = db_file)

# --- SCHEMA ---
table_exists <- dbExistsTable(con, "ext_seqs")

if (!table_exists) {
  dbExecute(con, "
  CREATE TABLE ext_seqs (
    sequence_id     TEXT PRIMARY KEY,
    genotype        TEXT,
    variant         TEXT,
    comment         TEXT,
    length_bp       INTEGER NOT NULL,
    collection_date TEXT,
    sequence        TEXT NOT NULL
  );
  ")
} else {
  cols <- dbGetQuery(con, "PRAGMA table_info(ext_seqs);")$name
  if (!"genotype" %in% cols) {
    dbExecute(con, "ALTER TABLE ext_seqs ADD COLUMN genotype TEXT;")
  }
}

dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_variant ON ext_seqs (variant);")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_date ON ext_seqs (collection_date);")
dbExecute(con, "CREATE INDEX IF NOT EXISTS idx_genotype ON ext_seqs (genotype);")

# --- INSERT FUNCTION ---
insert_seq <- function(sequence_id,
                       sequence,
                       genotype = NA,
                       variant = NA,
                       collection_date = NA,
                       comment = NA) {
  
  if (!grepl("^[ACGTUacgtu]+$", sequence)) {
    stop(paste("Invalid sequence:", sequence_id))
  }
  
  collection_date <- if (is.na(collection_date)) {
    NA_character_
  } else {
    as.character(as.Date(collection_date))
  }
  
  dbExecute(con,
            "INSERT OR REPLACE INTO ext_seqs
     (sequence_id, genotype, variant, comment, length_bp, collection_date, sequence)
     VALUES (?, ?, ?, ?, ?, ?, ?)",
            params = list(
              sequence_id,
              genotype,
              variant,
              comment,
              nchar(sequence),
              collection_date,
              toupper(sequence)
            )
  )
}

# --- FASTA + METADATA IMPORT ---
import_fasta_with_metadata <- function(
    fasta_file,
    metadata_tsv = "data/local_datasets/eksterne_sek/externe_metadata.tsv"
) {
  
  # metadata
  meta <- readr::read_tsv(
    metadata_tsv,
    col_types = cols(
      sequence_id = col_character(),
      genotype = col_character(),
      variant = col_character(),
      comment = col_character(),
      collection_date = col_character()
    )
  )
  
  if (any(duplicated(meta$sequence_id))) {
    stop("Duplicate sequence_id entries in metadata TSV")
  }
  
  # FASTA
  fasta_list <- read.fasta(
    file = fasta_file,
    seqtype = "DNA",
    as.string = FALSE,
    set.attributes = FALSE
  )
  
  raw_ids <- names(fasta_list)
  
  # strict validation: no whitespace allowed
  if (any(grepl("\\s", raw_ids))) {
    bad <- raw_ids[grepl("\\s", raw_ids)]
    stop(
      paste0(
        "Invalid FASTA headers: whitespace not allowed.\n",
        paste(bad, collapse = "\n")
      )
    )
  }
  
  clean_ids <- raw_ids
  
  # collision: within FASTA
  if (any(duplicated(clean_ids))) {
    dup <- clean_ids[duplicated(clean_ids)]
    stop(
      paste0(
        "Duplicate FASTA sequence IDs:\n",
        paste(unique(dup), collapse = "\n")
      )
    )
  }
  
  # collision: existing DB
  existing_ids <- dbGetQuery(con, "SELECT sequence_id FROM ext_seqs")$sequence_id
  overlap <- intersect(clean_ids, existing_ids)
  
  if (length(overlap) > 0) {
    stop(
      paste0(
        "Sequence IDs already exist in database:\n",
        paste(overlap, collapse = "\n")
      )
    )
  }
  
  # insert
  for (i in seq_along(fasta_list)) {
    
    seq_vec <- fasta_list[[i]]
    seq_str <- paste(seq_vec, collapse = "")
    
    id <- clean_ids[i]
    
    row <- meta %>% filter(sequence_id == id)
    
    if (nrow(row) == 0) {
      stop(paste("Missing metadata for:", id))
    }
    
    insert_seq(
      sequence_id = id,
      sequence = seq_str,
      genotype = row$genotype,
      variant = row$variant,
      collection_date = row$collection_date,
      comment = row$comment
    )
  }
}

# --- EXAMPLES ---
insert_seq("BRCA1", "ATGCGTACGTAGCTAGCTAGCTAGCTAGCTA", genotype = "G1", variant = "Ext1")
insert_seq("TP53",  "ATGGAGGAGCCGCAGTCAGATCCTAGCGTCG", genotype = "G2", variant = "Ext2")

import_fasta_with_metadata("data/local_datasets/eksterne_sek/Request RIV-HAV16-069.fa")

# --- QUERY ---
genes_df <- dbGetQuery(con, "SELECT * FROM ext_seqs;")
print(genes_df)

# --- DUMP ENTIRE DATABASE ---
export_database <- function(
    con,
    fasta_file = "ext_seqs",
    metadata_file = "ext_seqs_metadata"
) {
  
  df <- dbGetQuery(con, "SELECT * FROM ext_seqs")
  
  today <- format(Sys.Date(), "%Y-%m-%d")
  
  fasta_path <- paste0(fasta_file, "_", today, ".fasta")
  tsv_path   <- paste0(metadata_file, "_", today, ".tsv")
  
  # --- metadata export ---
  readr::write_tsv(df, tsv_path)
  
  # --- FASTA export ---
  fasta_conn <- file(fasta_path, open = "w")
  
  for (i in seq_len(nrow(df))) {
    
    header <- df$sequence_id[i]
    seq <- df$sequence[i]
    
    if (is.na(header) || is.na(seq)) {
      stop(paste("Invalid row in DB at index", i))
    }
    
    writeLines(paste0(">", header), fasta_conn)
    writeLines(seq, fasta_conn)
  }
  
  close(fasta_conn)
}

export_database(con)

dbDisconnect(con)