# ── Logging ───────────────────────────────────────────────────────────────────
# Snakemake only redirects output to `log:` for `shell:` directives. For
# `script:` rules the log file stays empty unless the script writes to it, which
# leaves failures visible only in the per-job cluster log. Redirect explicitly
# so show-failed-logs surfaces errors in the main slurm log.
if (length(snakemake@log) > 0) {
  .log_path <- snakemake@log[[1]]
  dir.create(dirname(.log_path), showWarnings = FALSE, recursive = TRUE)
  .log_con <- file(.log_path, open = "wt")
  sink(.log_con, type = "output", split = TRUE)
  sink(.log_con, type = "message")
}

library(tidyverse)

# ── Snakemake I/O ─────────────────────────────────────────────────────────────
master_metadata_path   <- snakemake@config$metadata$master_sheet
checkm2_report_dirs    <- snakemake@input$checkm2_reports
flye_info_files        <- snakemake@input$flye_info
gtdbtk_dir             <- snakemake@input$gtdbtk_summary
bakta_dirs             <- snakemake@input$bakta_dirs
antismash_overview_csv <- snakemake@input$antismash_overviews
asm_qc_summary_path    <- snakemake@input$asm_qc_summary
tree_membership_path   <- snakemake@input$tree_membership
antismash_fails_path   <- snakemake@input$antismash_fails
output_tsv             <- snakemake@output$tsv
output_mibig           <- snakemake@output$mibig
output_product         <- snakemake@output$product

# ── Constants ─────────────────────────────────────────────────────────────────
PATHOGEN_COLS <- c(
  "Enterococcus mundtii DSM 4838",
  "Enterobacter mori DSM 26271",
  "Staphylococcus epidermidis AU 24",
  "Escherichia coli DSM 498",
  "Pseudomonas putida DSM 6125",
  "Acinetobacter baylyi DSM 14961",
  "Enterococcus faecium DSM 20477",
  "Staphylococcus aureus DSM 20231",
  "Klebsiella oxytoca DSM 5175",
  "Acinetobacter baumannii DSM 300007",
  "Pseudomonas aeruginosa DSM 19880",
  "Enterobacter cloacae DSM 30054"
)

# ── Helpers ───────────────────────────────────────────────────────────────────
convert_eskape <- function(x) {
  score <- case_when(
    is.na(x)                          ~ NA_real_,
    x == "111"                        ~ 1,
    x == "000"                        ~ 0,
    x %in% c("001", "010", "100")    ~ 1 / 3,
    x %in% c("011", "110", "101")    ~ 2 / 3,
    TRUE                              ~ NA_real_
  )
  round(score, 0)
}

# ── Master metadata sheet ──────────────────────────────────────────────────────
message("Loading master metadata sheet...")
if (!is.null(master_metadata_path) && file.exists(master_metadata_path)) {
  master <- read_csv(master_metadata_path,
                     col_types = cols(.default = col_character()),
                     show_col_types = FALSE) %>%
    rename_with(~ str_replace_all(., "Actinetobacter", "Acinetobacter")) %>%
    filter(!is.na(SEQID)) %>%
    select(
      SEQID, ID,
      any_of(c("DNA conc. (ng/µl)",
               "Screening against Minimal Antibiotic Resistance Platform (Wright Lab)")),
      any_of(PATHOGEN_COLS)
    ) %>%
    mutate(across(any_of(PATHOGEN_COLS), convert_eskape)) %>%
    rename_with(~ paste("inhibit", .), any_of(PATHOGEN_COLS))
} else {
  message("Master metadata sheet not found — skipping bioactivity.")
  master <- tibble(
    SEQID = character(), ID = character()
  )
}

# ── Parse pipeline outputs ─────────────────────────────────────────────────────
message("Parsing assembly and QC metrics...")

# Flye — direct assembly_info.txt files (one per sample)
flye_raw <- map_dfr(flye_info_files, function(f) {
  read_delim(f, delim = "\t", col_types = cols(.default = "c"), show_col_types = FALSE) %>%
    rename_all(~ gsub(" ", "_", .)) %>%
    mutate(SEQID = str_remove(basename(f), "\\.assembly_info\\.txt$"))
})
flye <- if (nrow(flye_raw) > 0) {
  flye_raw %>%
    group_by(SEQID) %>%
    filter(length == max(as.numeric(length))) %>%
    ungroup() %>%
    select(SEQID, `Genome size` = length, circ., Coverage = cov.)
} else {
  tibble(SEQID = character(), `Genome size` = character(),
         circ. = character(), Coverage = character())
}

# CheckM2 — genome completeness and contamination
checkm2_raw <- map_dfr(checkm2_report_dirs, function(dir) {
  f <- file.path(dir, "quality_report.tsv")
  if (!file.exists(f)) return(NULL)
  read_tsv(f, col_types = cols(), show_col_types = FALSE) %>%
    mutate(SEQID = basename(dir))
})
checkm2 <- if (nrow(checkm2_raw) > 0) {
  checkm2_raw %>% select(SEQID, Completeness, Contamination)
} else {
  tibble(SEQID = character(), Completeness = numeric(), Contamination = numeric())
}

# GTDB-Tk — taxonomy
message("Parsing GTDB-Tk classification...")
files_gtdb <- list.files(gtdbtk_dir, pattern = "\\.summary\\.tsv$",
                         full.names = TRUE, recursive = TRUE)
gtdb_tk <- if (length(files_gtdb) > 0) {
  map_dfr(files_gtdb, ~ read_tsv(.x, col_types = cols(.default = "c"), show_col_types = FALSE)) %>%
    rename(SEQID = user_genome) %>%
    # gtdbtk was run with --extension fa.gz against files named
    # {sample}.flye.fa.gz, so it only strips ".fa.gz" and reports genome
    # IDs like "{sample}.flye". Strip the leftover ".flye" so SEQID lines
    # up with every other table (flye/checkm2/bakta/master), which all use
    # the bare "{sample}" as SEQID.
    mutate(SEQID = str_remove(SEQID, "\\.flye$")) %>%
    separate(classification, sep = ";",
             into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
             remove = FALSE) %>%
    mutate(
      simple_classification = paste(
        str_remove(Phylum, "p__"),
        str_remove(Genus,  "g__"),
        sep = " - g__"
      ),
      Species = case_when(
        !str_detect(Species, " ") ~ paste(str_remove(Genus, "g__"), "sp."),
        TRUE                      ~ str_remove(Species, "s__")
      )
    ) %>%
    select(SEQID, classification, simple_classification, Species, Phylum, Genus) %>%
    # list.files(recursive = TRUE) can match the same summary more than once
    # (e.g. a copy at the top level as well as under classify/), which stacks a
    # duplicate row per genome and doubles the whole output table. The
    # 2026-08-06 run produced 1098 rows for 549 genomes this way.
    distinct(SEQID, .keep_all = TRUE)
} else {
  tibble(SEQID = character(), classification = character(),
         simple_classification = character(),
         Species = character(), Phylum = character(), Genus = character())
}

# Bakta — rRNA and tRNA counts
message("Parsing Bakta annotations...")
bakta_raw <- map_dfr(bakta_dirs, function(dir) {
  seq_id <- basename(dir)
  f <- file.path(dir, paste0(seq_id, ".tsv"))
  if (!file.exists(f)) return(NULL)

  data <- read_delim(f, skip = 5, delim = "\t", show_col_types = FALSE)

  rRNA <- data %>%
    filter(str_detect(Product, "ribosomal RNA")) %>%
    count(Product) %>%
    pivot_wider(names_from = Product, values_from = n, values_fill = 0L)
  for (col in c("5S ribosomal RNA", "16S ribosomal RNA", "23S ribosomal RNA"))
    if (!col %in% colnames(rRNA)) rRNA[[col]] <- 0L

  tRNA <- data %>%
    filter(str_detect(Product, "tRNA-[a-zA-Z]{3}")) %>%
    mutate(Product = str_remove(Product, "\\([a-zA-Z]+\\)$")) %>%
    filter(!str_detect(Product, "(pseudo)")) %>%
    distinct(Product) %>%
    summarise(Number_of_tRNA = n())

  bind_cols(rRNA, tRNA) %>% mutate(SEQID = seq_id)
})
hq_bakta <- if (nrow(bakta_raw) > 0) {
  bakta_raw %>%
    mutate(`# 5S rRNA/16S rRNA / 23S rRNA` =
             paste(`5S ribosomal RNA`, `16S ribosomal RNA`, `23S ribosomal RNA`, sep = "/")) %>%
    select(SEQID, `5S ribosomal RNA`, `16S ribosomal RNA`, `23S ribosomal RNA`,
           `# 5S rRNA/16S rRNA / 23S rRNA`, Number_of_tRNA)
} else {
  tibble(SEQID = character(), `5S ribosomal RNA` = numeric(),
         `16S ribosomal RNA` = numeric(), `23S ribosomal RNA` = numeric(),
         `# 5S rRNA/16S rRNA / 23S rRNA` = character(), Number_of_tRNA = numeric())
}

# Assembly QC gate — the hard filter applied before tree inference
message("Parsing assembly QC verdicts...")
asm_qc <- if (!is.null(asm_qc_summary_path) && file.exists(asm_qc_summary_path)) {
  read_tsv(asm_qc_summary_path, col_types = cols(.default = col_character()),
           show_col_types = FALSE) %>%
    transmute(
      SEQID              = sample,
      Contigs            = as.numeric(contigs),
      Passed_assembly_QC = status == "PASS",
      Assembly_QC_reason = if_else(status == "PASS", NA_character_, reason)
    )
} else {
  tibble(SEQID = character(), Contigs = numeric(),
         Passed_assembly_QC = logical(), Assembly_QC_reason = character())
}

# Tree membership — bac120 marker recovery and final tree inclusion.
# Passing assembly QC is necessary but not sufficient: GTDB-Tk's align step
# drops genomes with too few amino acids in the marker alignment, so those
# genomes have no tip even though their annotation and BGC results are kept.
message("Parsing tree membership...")
tree_membership <- if (!is.null(tree_membership_path) &&
                       file.exists(tree_membership_path)) {
  read_tsv(tree_membership_path, col_types = cols(.default = col_character()),
           show_col_types = FALSE) %>%
    transmute(
      SEQID                 = sample,
      In_bac120_MSA         = in_bac120_msa == "TRUE",
      N_bac120_markers      = suppressWarnings(as.numeric(n_bac120_markers)),
      In_bac120_tree        = in_bac120_tree == "TRUE",
      Tree_exclusion_reason = if_else(tree_exclusion_reason == "-",
                                      NA_character_, tree_exclusion_reason)
    )
} else {
  tibble(SEQID = character(), In_bac120_MSA = logical(),
         N_bac120_markers = numeric(), In_bac120_tree = logical(),
         Tree_exclusion_reason = character())
}

# antiSMASH failures — genomes antiSMASH could not process. These are recorded
# rather than fatal, so they must be visible here or they would look simply like
# genomes with no BGCs.
message("Parsing antiSMASH failures...")
antismash_fails <- if (!is.null(antismash_fails_path) &&
                       file.exists(antismash_fails_path)) {
  read_tsv(antismash_fails_path, col_types = cols(.default = col_character()),
           show_col_types = FALSE) %>%
    transmute(SEQID = sample, AntiSMASH_failure = reason)
} else {
  tibble(SEQID = character(), AntiSMASH_failure = character())
}

# ── Merge pipeline data ────────────────────────────────────────────────────────
message("Compiling core assembly table...")
seq_data <- flye %>%
  left_join(checkm2,  by = "SEQID") %>%
  left_join(gtdb_tk,  by = "SEQID") %>%
  left_join(hq_bakta, by = "SEQID") %>%
  left_join(asm_qc,   by = "SEQID") %>%
  left_join(tree_membership, by = "SEQID") %>%
  left_join(antismash_fails, by = "SEQID") %>%
  mutate(
    AntiSMASH_OK = is.na(AntiSMASH_failure),
    GenomeQ = case_when(
      as.numeric(Completeness) >= 90 & as.numeric(Contamination) < 5 &
        as.numeric(`23S ribosomal RNA`) >= 1 & as.numeric(`16S ribosomal RNA`) >= 1 &
        as.numeric(`5S ribosomal RNA`) >= 1 & as.numeric(Number_of_tRNA) >= 18 ~ "High quality genome",
      as.numeric(Completeness) >= 50 & as.numeric(Contamination) < 10           ~ "Medium quality genome",
      as.numeric(Completeness) <  50 & as.numeric(Contamination) < 10           ~ "Low quality genome",
      TRUE                                                                        ~ "Unannotated genomes"
    )
  )

# ── Join master metadata ───────────────────────────────────────────────────────
message("Joining master metadata...")
results <- seq_data %>%
  left_join(master, by = "SEQID")

# ── Write main results ─────────────────────────────────────────────────────────
# One row per genome. Any join fanning out is a bug, not something to paper over,
# so say so rather than silently deduplicating.
if (anyDuplicated(results$SEQID)) {
  n_dup <- sum(duplicated(results$SEQID))
  warning("aggregated results contain ", n_dup, " duplicate SEQID rows; ",
          "a join has fanned out. Deduplicating, but this should be fixed.")
  results <- distinct(results, SEQID, .keep_all = TRUE)
}
message("Writing aggregated results TSV: ", nrow(results), " genomes.")
write_tsv(results, file = output_tsv)

# ── antiSMASH BGC processing ───────────────────────────────────────────────────
message("Processing antiSMASH biosynthetic gene clusters...")

overview_data_combined <- map_dfr(antismash_overview_csv, function(f) {
  if (!file.exists(f)) return(NULL)
  read_delim(f, delim = ";", col_types = cols(.default = col_character()),
             show_col_types = FALSE) %>%
    mutate(SEQID = str_remove(basename(f), "_overview\\.csv$"))
})

if (nrow(overview_data_combined) > 0) {
  overview_data_combined <- overview_data_combined %>%
    mutate(qualifiers_product = str_remove_all(qualifiers_product, '^c\\(|\\)|"'))

  BGC_types <- overview_data_combined %>%
    mutate(qualifiers_product_list = strsplit(qualifiers_product, ", ")) %>%
    unnest(qualifiers_product_list) %>%
    filter(!is.na(qualifiers_product_list), qualifiers_product_list != "") %>%
    count(SEQID, qualifiers_product_list, name = "number_product")

  categorize_cluster <- function(x) {
    if (is.na(x)) return("BGC not in MiBiG database")
    cats <- c(
      if (grepl("NRP",        x, ignore.case = TRUE)) "NRP",
      if (grepl("Polyketide", x, ignore.case = TRUE)) "Polyketide",
      if (grepl("RiPP",       x, ignore.case = TRUE)) "RiPP",
      if (grepl("Terpene",    x, ignore.case = TRUE)) "Terpene",
      if (grepl("Alkaloid",   x, ignore.case = TRUE)) "Alkaloid",
      if (grepl("Saccharide", x, ignore.case = TRUE)) "Saccharide"
    )
    if (length(cats) == 0 || grepl("Other", x, ignore.case = TRUE)) cats <- c(cats, "Other")
    paste(cats, collapse = "+")
  }

  mibig_summary <- overview_data_combined %>%
    filter(qualifiers_contig_edge != "True") %>%
    mutate(simplified_cluster_type = sapply(knowncblast_cluster_type, categorize_cluster)) %>%
    count(SEQID, simplified_cluster_type, name = "Number_mibig_type")

  write_tsv(mibig_summary, file = output_mibig)
  write_tsv(BGC_types,     file = output_product)
} else {
  write_tsv(
    tibble(SEQID = character(), simplified_cluster_type = character(), Number_mibig_type = numeric()),
    file = output_mibig
  )
  write_tsv(
    tibble(SEQID = character(), qualifiers_product_list = character(), number_product = numeric()),
    file = output_product
  )
}

message("Compilation complete!")
