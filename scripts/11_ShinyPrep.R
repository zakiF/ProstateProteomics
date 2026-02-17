#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(preprocessCore)
})

normalize_sample_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("-", "_", x)
  x <- gsub("^c(?=[0-9])", "c_", x, perl = TRUE)
  x <- gsub("^c([0-9]{4}_[0-9]+)$", "c_\\1", x)
  x
}

parse_proteins_with_sample_type <- function(markers_text) {
  proteins_list <- list()
  if (is.na(markers_text) || markers_text == "") return(proteins_list)

  has_tissue <- grepl("Tissue-based\\s*\\(", markers_text, ignore.case = TRUE)
  has_plasma <- grepl("Plasma-based\\s*\\(", markers_text, ignore.case = TRUE)
  has_urine <- grepl("Urine-based\\s*\\(", markers_text, ignore.case = TRUE)

  extract_group <- function(pattern, sample_type) {
    m <- regmatches(markers_text, gregexpr(pattern, markers_text, ignore.case = TRUE, perl = TRUE))[[1]]
    if (length(m) == 0) return(NULL)
    vals <- gsub(".*\\(([^)]+)\\).*", "\\1", m)
    vals <- unlist(strsplit(vals, ",\\s*"))
    vals <- trimws(gsub("^and\\s+", "", vals, ignore.case = TRUE))
    vals <- vals[vals != ""]
    lapply(vals, function(v) list(name = v, sample_type = sample_type))
  }

  if (has_tissue || has_plasma || has_urine) {
    if (has_tissue) proteins_list <- c(proteins_list, extract_group("Tissue-based\\s*\\(([^)]+)\\)", "Tissue"))
    if (has_plasma) proteins_list <- c(proteins_list, extract_group("Plasma-based\\s*\\(([^)]+)\\)", "Plasma"))
    if (has_urine) proteins_list <- c(proteins_list, extract_group("Urine-based\\s*\\(([^)]+)\\)", "Urine"))
  } else {
    vals <- unlist(strsplit(markers_text, ",\\s*"))
    vals <- trimws(gsub("^and\\s+", "", vals, ignore.case = TRUE))
    vals <- vals[vals != ""]
    proteins_list <- lapply(vals, function(v) list(name = v, sample_type = NA_character_))
  }
  proteins_list
}

root_dir <- if (basename(getwd()) == "scripts") normalizePath("..") else normalizePath(".")
shiny_data_dir <- file.path(root_dir, "scripts", "ShinyData")
if (!dir.exists(shiny_data_dir)) dir.create(shiny_data_dir, recursive = TRUE)

paths <- list(
  expr = file.path(root_dir, "output", "data", "allSamples_expression_matrix.csv"),
  olink_meta = file.path(root_dir, "output", "data", "allSamples_metadata.csv"),
  clin = file.path(root_dir, "output", "data", "mCRPC_metadata.csv"),
  stats = file.path(root_dir, "output", "06_mCRPC_analysis", "HR_values_mCRPC_specific_multi_pre-allProteins.csv"),
  dep_sig = file.path(root_dir, "output", "06_mCRPC_analysis", "HR_values_mCRPC_specific_multivariate_cell_surface.csv"),
  matched = file.path(root_dir, "output", "09_public", "Matched_Olink_SWATHMS_Proteins.csv"),
  public_raw = file.path(root_dir, "data", "public", "41467_2018_3573_MOESM4_ESM.xlsx"),
  public_meta = file.path(root_dir, "data", "public", "sample_metadata.csv"),
  biomarker = file.path(root_dir, "data", "updated_2024", "Protein biomarker pubs mCRPC_zaki_edit.xlsx"),
  cs = file.path(root_dir, "output", "data", "CS_human.txt"),
  drug_target = file.path(root_dir, "data", "DrugBank", "drugbank_all_target_polypeptide_ids.csv", "pharmacologically_active.csv"),
  drug_vocab = file.path(root_dir, "data", "DrugBank", "drugbank vocabulary.csv"),
  drug_manual = file.path(root_dir, "data", "Drug_review", "protein_drug_associations.tsv"),
  out = file.path(shiny_data_dir, "11_shiny_prep.rds"),
  out_rdata = file.path(shiny_data_dir, "11_shiny_prep.RData")
)

required <- c("expr", "olink_meta", "clin", "stats", "dep_sig", "matched", "public_raw", "public_meta")
missing_files <- required[!file.exists(unlist(paths[required]))]
if (length(missing_files) > 0) {
  stop("Missing required input files: ", paste(missing_files, collapse = ", "))
}

cat("Loading inputs...\n")
expr <- readr::read_csv(paths$expr, show_col_types = FALSE)
olink_meta <- readr::read_csv(paths$olink_meta, show_col_types = FALSE)
clin <- readr::read_csv(paths$clin, show_col_types = FALSE)
protein_stats_raw <- readr::read_csv(paths$stats, show_col_types = FALSE)
df_deps <- readr::read_csv(paths$dep_sig, show_col_types = FALSE)
matched_proteins <- readr::read_csv(paths$matched, show_col_types = FALSE)
public_raw <- readxl::read_excel(paths$public_raw, sheet = "Area - proteins")
public_meta <- readr::read_csv(paths$public_meta, show_col_types = FALSE)

colnames(expr)[1] <- "OlinkID"

protein_stats <- protein_stats_raw
if (!"Multi_Coefficient" %in% names(protein_stats) && "Coefficient" %in% names(protein_stats)) {
  protein_stats <- protein_stats %>% mutate(Multi_Coefficient = Coefficient)
}
if (!"Multi_HR" %in% names(protein_stats) && "HR" %in% names(protein_stats)) {
  protein_stats <- protein_stats %>% mutate(Multi_HR = HR)
}
if (!"Multi_FDR_num" %in% names(protein_stats) && "fdr_cont" %in% names(protein_stats)) {
  protein_stats <- protein_stats %>% mutate(Multi_FDR_num = fdr_cont)
}
if (!"Multi_Sig" %in% names(protein_stats) && "Sig_fdr_cont" %in% names(protein_stats)) {
  protein_stats <- protein_stats %>%
    mutate(Multi_Sig = ifelse(tolower(as.character(Sig_fdr_cont)) %in% c("yes", "sig", "true"), "Yes", "No"))
}
if (!"Uni_Coefficient" %in% names(protein_stats)) {
  protein_stats <- protein_stats %>% mutate(Uni_Coefficient = NA_real_)
}
for (cn in c("Multi_Coefficient", "Multi_HR", "Multi_FDR_num", "Multi_Sig")) {
  if (!cn %in% names(protein_stats)) protein_stats[[cn]] <- NA
}

protein_stats <- protein_stats %>%
  mutate(
    Uni_Coefficient = as.numeric(Uni_Coefficient),
    Multi_Coefficient = as.numeric(Multi_Coefficient),
    Multi_HR = as.numeric(Multi_HR),
    Multi_FDR_num = as.numeric(Multi_FDR_num)
  ) %>%
  select(OlinkID, Assay, Uni_Coefficient, Multi_Coefficient, Multi_HR, Multi_FDR_num, Multi_Sig) %>%
  distinct(OlinkID, .keep_all = TRUE)

df_hr_all <- protein_stats_raw %>%
  mutate(HR = as.numeric(HR), Lower_Bound = as.numeric(Lower_Bound), Upper_Bound = as.numeric(Upper_Bound)) %>%
  filter(!is.na(Assay), Assay != "") %>%
  distinct(OlinkID, Assay, .keep_all = TRUE)

dep_unique <- df_deps %>%
  arrange(Multi_FDR_num, desc(abs(Multi_Coefficient))) %>%
  distinct(Assay, .keep_all = TRUE)

assay_map <- bind_rows(
  protein_stats %>% select(OlinkID, Assay),
  matched_proteins %>% select(OlinkID, Assay)
) %>%
  filter(!is.na(Assay), !is.na(OlinkID)) %>%
  distinct(OlinkID, .keep_all = TRUE)

expr_long <- expr %>%
  pivot_longer(cols = -OlinkID, names_to = "sample_id", values_to = "NPX") %>%
  mutate(sample_id_norm = normalize_sample_id(sample_id)) %>%
  left_join(assay_map, by = "OlinkID")

if ("HCI_cID" %in% names(olink_meta)) {
  olink_meta <- olink_meta %>% mutate(sample_id_norm = normalize_sample_id(HCI_cID))
} else if ("sample_id" %in% names(olink_meta)) {
  olink_meta <- olink_meta %>% mutate(sample_id_norm = normalize_sample_id(sample_id))
} else {
  stop("Cannot find HCI/sample id in allSamples_metadata.csv")
}

npx_long <- expr_long %>%
  left_join(
    olink_meta %>% select(sample_id_norm, cohort, txt_stat_ordered) %>% distinct(),
    by = "sample_id_norm"
  )

# Treatment-naive mCRPC IDs (cohort C pre-treatment)
c_pre_ids <- olink_meta %>%
  filter(cohort == "C", grepl("^Pre", txt_stat_ordered)) %>%
  distinct(sample_id_norm) %>%
  pull(sample_id_norm)

clin_sample_col <- if ("HCI_cID" %in% names(clin)) {
  "HCI_cID"
} else if ("sample_id" %in% names(clin)) {
  "sample_id"
} else {
  stop("Cannot find sample ID column in mCRPC_metadata.csv (expected HCI_cID or sample_id)")
}

analysis_meta <- clin %>%
  mutate(
    sample_id = as.character(.data[[clin_sample_col]]),
    sample_id_norm = normalize_sample_id(sample_id),
    log10_psa = log10(as.numeric(psa) + 0.1),
    log10_alp = log10(as.numeric(alk_phos)),
    log10_ldh = log10(as.numeric(ldh)),
    OS = as.numeric(OS),
    death = as.numeric(death)
  ) %>%
  select(sample_id, sample_id_norm, OS, death, log10_psa, log10_alp, log10_ldh) %>%
  distinct(sample_id_norm, .keep_all = TRUE) %>%
  filter(sample_id_norm %in% c_pre_ids)

# Restrict NPX to treatment-naive mCRPC only
npx_long_c_pre <- npx_long %>%
  filter(cohort == "C", grepl("^Pre", txt_stat_ordered), sample_id_norm %in% analysis_meta$sample_id_norm)

public_long <- public_raw %>%
  pivot_longer(cols = -Protein, names_to = "Sample_Name", values_to = "Intensity") %>%
  left_join(public_meta, by = "Sample_Name")

public_avg <- public_long %>%
  group_by(Protein, Biological_ID, Group, Group_Number) %>%
  summarise(Intensity_Avg = mean(Intensity, na.rm = TRUE), .groups = "drop")

public_wide <- public_avg %>%
  pivot_wider(id_cols = Protein, names_from = Biological_ID, values_from = Intensity_Avg)
public_matrix <- public_wide %>% select(-Protein) %>% as.matrix()
public_log2 <- log2(public_matrix + 1)
public_norm <- preprocessCore::normalize.quantiles(public_log2)
colnames(public_norm) <- colnames(public_log2)
rownames(public_norm) <- public_wide$Protein

biom_unique_proteins <- character(0)
if (file.exists(paths$biomarker)) {
  protein_table <- readxl::read_excel(paths$biomarker)
  marker_col <- protein_table$Edited_name
  all_entries <- lapply(marker_col, parse_proteins_with_sample_type)
  all_proteins_raw <- lapply(all_entries, function(study) sapply(study, function(p) p$name)) %>% unlist()
  biom_unique_proteins <- sort(unique(all_proteins_raw))
}

cs_prot <- character(0)
if (file.exists(paths$cs)) {
  cs <- read.delim(paths$cs)
  if ("Approved.symbol" %in% names(cs)) cs_prot <- unique(cs$Approved.symbol)
}

drug_annotated_proteins <- character(0)
if (file.exists(paths$drug_target) && file.exists(paths$drug_vocab) && file.exists(paths$drug_manual)) {
  df_target <- read.csv(paths$drug_target)
  df_vocab <- read.csv(paths$drug_vocab)
  df_target_manual <- read.delim(paths$drug_manual, sep = "\t", stringsAsFactors = FALSE)

  df_target_sub <- df_target %>%
    select(Gene.Name, Drug.IDs) %>%
    mutate(Drug.IDs = str_split(Drug.IDs, ",\\s*")) %>%
    tidyr::unnest(Drug.IDs) %>%
    filter(!is.na(Gene.Name), !is.na(Drug.IDs), Drug.IDs != "") %>%
    distinct(Gene.Name, Drug.IDs) %>%
    group_by(Gene.Name) %>%
    summarise(DrugBank.ID = str_c(sort(unique(Drug.IDs)), collapse = "; "), .groups = "drop")

  df_vocab_sub <- df_vocab %>% select(DrugBank.ID, Common.name)

  df_target_ann <- df_target_sub %>%
    mutate(
      DrugBank.ID_list = str_split(DrugBank.ID, ";\\s*"),
      drugname = purrr::map_chr(DrugBank.ID_list, function(ids) {
        names <- df_vocab_sub %>%
          filter(DrugBank.ID %in% ids) %>%
          pull(Common.name) %>%
          unique()
        str_c(names, collapse = "; ")
      })
    ) %>%
    select(-DrugBank.ID_list)

  dt_manual <- df_target_manual %>% filter(Drug_Association == "Yes") %>% pull(Protein)
  drug_annotated_proteins <- sort(unique(c(
    df_target_ann %>% filter(!is.na(DrugBank.ID), DrugBank.ID != "") %>% pull(Gene.Name),
    dt_manual
  )))
}

prep <- list(
  created_at = Sys.time(),
  npx_long = npx_long,
  npx_long_c_pre = npx_long_c_pre,
  protein_stats = protein_stats,
  df_hr_all = df_hr_all,
  dep_unique = dep_unique,
  analysis_meta = analysis_meta,
  matched_proteins = matched_proteins,
  assay_map = assay_map,
  public_avg = public_avg,
  public_norm = public_norm,
  biom_unique_proteins = biom_unique_proteins,
  cs_prot = cs_prot,
  drug_annotated_proteins = drug_annotated_proteins,
  time_points = c(24, 36)
)

saveRDS(prep, paths$out)
prep_shiny <- prep
save(prep_shiny, file = paths$out_rdata, compress = "xz")

writeLines(
  c(
    "ShinyData bundle generated by scripts/11_ShinyPrep.R",
    paste0("Created: ", as.character(Sys.time())),
    "Files:",
    "- 11_shiny_prep.rds",
    "- 11_shiny_prep.RData (object name: prep_shiny)"
  ),
  con = file.path(shiny_data_dir, "README.txt")
)

cat("Saved prep data to:", paths$out, "\n")
cat("Saved prep .RData to:", paths$out_rdata, "\n")
cat("Samples in analysis_meta:", nrow(analysis_meta), "\n")
cat("NPX rows (all):", nrow(npx_long), "\n")
cat("NPX rows (cohort C pre):", nrow(npx_long_c_pre), "\n")
cat("Public proteins:", nrow(public_norm), "\n")
cat("Proteins with non-missing Multi_Coefficient:", sum(!is.na(protein_stats$Multi_Coefficient)), "\n")
