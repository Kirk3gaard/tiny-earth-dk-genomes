library(tidyverse)
library(jsonlite)

# ── Snakemake I/O ─────────────────────────────────────────────────────────────
sample_id          <- snakemake@wildcards$sample
antismash_dir      <- snakemake@input[["data"]]
overview_path      <- snakemake@output[["overview_path"]]
cluster_blast_path <- snakemake@output[["cluster_blast_path"]]
known_cblast_path  <- snakemake@output[["known_cblast_path"]]
areas_path         <- snakemake@output[["areas_path"]]

# ── Load JSON ─────────────────────────────────────────────────────────────────
json_file <- file.path(antismash_dir, paste0(sample_id, ".json"))
if (!file.exists(json_file)) {
  stop("antiSMASH JSON not found: ", json_file)
}

message("[", sample_id, "] Loading antiSMASH JSON...")
records <- fromJSON(json_file, flatten = FALSE)$records

# Helper: collapse list-columns to character strings
collapse_lists <- function(df) {
  mutate(df, across(
    where(is.list),
    ~ map_chr(., function(x) {
      if (is.null(x) || length(x) == 0) NA_character_ else paste(x, collapse = ", ")
    })
  ))
}

# ── Features (regions & protoclusters) ────────────────────────────────────────
message("[", sample_id, "] Extracting features...")
features <- records[, c("id", "features")] |>
  as.data.frame() |>
  unnest(features) |>
  unnest_wider(qualifiers, names_sep = "_") |>
  collapse_lists()

has_bgc <- nrow(filter(features, type %in% c("region", "proclusters"))) > 0

if (!has_bgc) {
  msg <- paste0(
    sample_id, ": no regions or protoclusters detected by antiSMASH.\n",
    "Consider re-running via the antiSMASH web interface."
  )
  message(msg)
  writeLines(msg, file.path(dirname(overview_path), paste0(sample_id, "_no_bgc.txt")))
  empty_df <- tibble()
  write_delim(empty_df, file = overview_path,      delim = ";")
  write_delim(empty_df, file = cluster_blast_path, delim = ";")
  write_delim(empty_df, file = known_cblast_path,  delim = ";")
  write_delim(empty_df, file = areas_path,         delim = ";")
  quit(save = "no", status = 0)
}

features_bgc <- features |>
  filter(type %in% c("region", "proclusters")) |>
  select(
    id,
    region_number = qualifiers_region_number,
    type, location,
    qualifiers_contig_edge,
    qualifiers_product,
    qualifiers_candidate_cluster_numbers,
    qualifiers_rules
  )

# ── Areas ─────────────────────────────────────────────────────────────────────
message("[", sample_id, "] Extracting areas...")
areas <- records[, c("id", "areas")] |>
  as.data.frame() |>
  unnest(areas) |>
  unnest_wider(protoclusters, names_sep = "_") |>
  select(id, start, end, products) |>
  mutate(
    products = map_chr(products, str_c, collapse = ","),
    location = paste0("[", start, ":", end, "]")
  )

# ── ClusterBlast ──────────────────────────────────────────────────────────────
message("[", sample_id, "] Extracting clusterblast results...")
clusterblast <- records[, "modules", drop = FALSE]$modules$`antismash.modules.clusterblast`

general_blast <- as.data.frame(clusterblast$general) |>
  unnest(results) |>
  unnest(ranking) |>
  select(prefix, record_id, region_number, total_hits, ranking) |>
  unnest_wider(ranking, names_sep = "_") |>
  collapse_lists() |>
  group_by(region_number, record_id) |>
  mutate(across(starts_with("ranking_"), ~ map_chr(as.list(.), str_c, collapse = ","))) |>
  distinct(.keep_all = TRUE) |>
  ungroup() |>
  mutate(
    region_number = as.character(region_number),
    id = as.character(record_id)
  ) |>
  rename_with(~ str_replace(., "^ranking_", "generalcblast_"), starts_with("ranking_")) |>
  select(id, region_number, total_hits, everything(), -record_id)

known_blast_raw <- as.data.frame(clusterblast$knowncluster) |>
  unnest(results) |>
  unnest(ranking, keep_empty = TRUE) |>
  select(prefix, record_id, region_number, total_hits, ranking) |>
  unnest_wider(ranking, names_sep = "_") |>
  collapse_lists()

if ("ranking_accession" %in% colnames(known_blast_raw)) {
  message("[", sample_id, "] MiBiG hits found — building known clusterblast table...")
  known_blast <- known_blast_raw |>
    group_by(region_number, record_id) |>
    mutate(across(
      starts_with("ranking_"),
      ~ if (is.list(.)) map_chr(., str_c, collapse = ",") else as.character(.)
    )) |>
    distinct(region_number, record_id, .keep_all = TRUE) |>
    ungroup() |>
    mutate(
      region_number = as.character(region_number),
      id = as.character(record_id)
    ) |>
    rename_with(~ str_replace(., "^ranking_", "knowncblast_"), starts_with("ranking_")) |>
    select(id, region_number, total_hits,
           knowncblast_accession, knowncblast_description,
           knowncblast_cluster_type, knowncblast_similarity)
  overview <- left_join(features_bgc, known_blast, by = c("id", "region_number"))
  write_delim(known_blast_raw, file = known_cblast_path, delim = ";")
} else {
  message("[", sample_id, "] No MiBiG hits — writing feature overview without MiBiG columns.")
  overview <- left_join(
    features_bgc,
    mutate(known_blast_raw,
      region_number = as.character(region_number),
      id = as.character(record_id)
    ),
    by = c("id", "region_number")
  )
  write_delim(tibble(), file = known_cblast_path, delim = ";")
}

# ── Write outputs ─────────────────────────────────────────────────────────────
message("[", sample_id, "] Writing outputs...")
write_delim(overview,      file = overview_path,      delim = ";")
write_delim(general_blast, file = cluster_blast_path, delim = ";")
write_delim(areas,         file = areas_path,         delim = ";")
message("[", sample_id, "] Done.")
