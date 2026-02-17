## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(here)
library(tidyverse)
library(readxl)
library(patchwork)
library(preprocessCore)
library(survival)
library(timeROC)
library(ComplexHeatmap)
library(circlize)
library(grid)

# Set default theme
theme_set(theme_bw())


## ----directories--------------------------------------------------------------
wd <- list()
wd$main <- here()
wd$data <- file.path(wd$main, "data")
wd$public <- file.path(wd$data, "public")
wd$output <- file.path(wd$main, "output")
wd$outData <- file.path(wd$output, "data")
wd$out06 <- file.path(wd$output, "06_mCRPC_analysis")
wd$out09 <- file.path(wd$output, "09_public")
wd$outCurr <- file.path(wd$output, "10_integration")

# Create output directory if it doesn't exist
if (!file.exists(wd$outCurr)) {
  dir.create(wd$outCurr, recursive = TRUE)
  cat("✓ Created output directory:", wd$outCurr, "\n")
} else {
  cat("✓ Output directory exists:", wd$outCurr, "\n")
}


## -----------------------------------------------------------------------------
# Public SWATH-MS data (same source used in script 09)
public_raw <- read_excel(
  file.path(wd$public, "41467_2018_3573_MOESM4_ESM.xlsx"),
  sheet = "Area - proteins"
)
public_meta <- readr::read_csv(
  file.path(wd$public, "sample_metadata.csv"),
  show_col_types = FALSE
)

# Average technical replicates to biological sample level
public_long <- public_raw %>%
  pivot_longer(cols = -Protein, names_to = "Sample_Name", values_to = "Intensity") %>%
  left_join(public_meta, by = "Sample_Name")
public_avg <- public_long %>%
  group_by(Protein, Biological_ID, Group, Group_Number) %>%
  summarise(Intensity_Avg = mean(Intensity, na.rm = TRUE), .groups = "drop")

# Log2 + quantile normalization
public_wide <- public_avg %>%
  pivot_wider(id_cols = Protein, names_from = Biological_ID, values_from = Intensity_Avg)
public_matrix <- public_wide %>% select(-Protein) %>% as.matrix()
public_log2 <- log2(public_matrix + 1)
public_norm <- preprocessCore::normalize.quantiles(public_log2)
colnames(public_norm) <- colnames(public_log2)
rownames(public_norm) <- public_wide$Protein

cat("Public data loaded:", nrow(public_norm), "proteins x", ncol(public_norm), "samples\n")


## -----------------------------------------------------------------------------
# Olink expression (all cohorts) and metadata
olink_expr <- readr::read_csv(
  file.path(wd$outData, "allSamples_expression_matrix.csv"),
  show_col_types = FALSE
)
colnames(olink_expr)[1] <- "OlinkID"

olink_meta <- readr::read_csv(
  file.path(wd$outData, "allSamples_metadata.csv"),
  show_col_types = FALSE
)

# Clinical metadata for survival variables (requested: use OS from this file)
meta2_olink <- read.csv(file.path(wd$outData, "proteomics_metadata.csv"))
meta2_olink <- meta2_olink %>%
  mutate(
    psa_num = readr::parse_number(as.character(psa), na = c("", "NA", "N/A", "n/a", "N/a")),
    alk_phos_num = readr::parse_number(as.character(alk_phos), na = c("", "NA", "N/A", "n/a", "N/a")),
    hemoGlo_num = readr::parse_number(as.character(hemoGlo), na = c("", "NA", "N/A", "n/a", "N/a"))
  ) %>%
  mutate(
    HCI_cID = paste0("c_", gsub("-", "_", HCI_cID)),
    death = ifelse(tolower(death_status) == "y", 1, 0),
    log10_psa = log10(psa_num + 0.1),
    log10_alp = log10(alk_phos_num),
    log10_ldh = log10(hemoGlo_num)
  ) %>%
  select(-psa_num, -alk_phos_num, -hemoGlo_num)

# Curated mCRPC cohort C clinical data (used for CPPS objective)
meta_mcrpc <- read.csv(file.path(wd$outData, "mCRPC_metadata.csv")) %>%
  mutate(
    log10_psa = log10(psa + 0.1),
    log10_alp = log10(alk_phos),
    log10_ldh = log10(ldh)
  )

# Olink long format
df_npx_long <- olink_expr %>%
  pivot_longer(-OlinkID, names_to = "HCI_cID", values_to = "NPX") %>%
  left_join(
    olink_meta %>% select(HCI_cID, cohort, txt_stat_ordered),
    by = "HCI_cID"
  )

# Cohort C pre-treatment subset for CPPS prognostic modeling
df_npx_long_c_pre <- df_npx_long %>%
  filter(cohort == "C", grepl("^Pre", txt_stat_ordered))

meta_in <- meta_mcrpc %>% rename(sample_id = HCI_cID)

cat("Olink data loaded:", n_distinct(df_npx_long$OlinkID), "proteins\n")


## -----------------------------------------------------------------------------
matched_proteins <- readr::read_csv(
  file.path(wd$out09, "Matched_Olink_SWATHMS_Proteins.csv"),
  show_col_types = FALSE
)

matched_olink_ids <- unique(matched_proteins$OlinkID)
matched_public_ids <- unique(matched_proteins$Protein)

# Public matrix restricted to matched proteins
public_norm_matched <- public_norm[rownames(public_norm) %in% matched_public_ids, , drop = FALSE]

cat("Matched proteins:", nrow(matched_proteins), "\n")


## -----------------------------------------------------------------------------
# Read the protein biomarker table
protein_table <- read_excel(file.path(wd$data, "updated_2024/Protein biomarker pubs mCRPC_zaki_edit.xlsx"))

# Use Edited_name column for protein names (cleaner version)
protein_markers_col <- protein_table$`Edited_name`
# Function to parse sample-type-specific proteins
# Handles formats like "Tissue-based (DST, COL3A1) Urine-based (COL1A1, COL3A1)"
parse_proteins_with_sample_type <- function(markers_text, study_index) {
  proteins_list <- list()

  if (is.na(markers_text) || markers_text == "") {
    return(proteins_list)
  }

  # Check if the text contains sample type markers
  has_tissue <- grepl("Tissue-based\\s*\\(", markers_text, ignore.case = TRUE)
  has_plasma <- grepl("Plasma-based\\s*\\(", markers_text, ignore.case = TRUE)
  has_urine <- grepl("Urine-based\\s*\\(", markers_text, ignore.case = TRUE)

  if (has_tissue || has_plasma || has_urine) {
    # Parse each sample type separately

    # Extract Tissue-based proteins
    if (has_tissue) {
      tissue_match <- regmatches(markers_text, gregexpr("Tissue-based\\s*\\(([^)]+)\\)", markers_text, ignore.case = TRUE, perl = TRUE))[[1]]
      if (length(tissue_match) > 0) {
        tissue_proteins <- gsub(".*\\(([^)]+)\\).*", "\\1", tissue_match)
        tissue_proteins <- unlist(strsplit(tissue_proteins, ",\\s*"))
        tissue_proteins <- trimws(tissue_proteins)
        tissue_proteins <- gsub("^and\\s+", "", tissue_proteins, ignore.case = TRUE)
        tissue_proteins <- trimws(tissue_proteins)
        tissue_proteins <- tissue_proteins[tissue_proteins != ""]

        for (p in tissue_proteins) {
          proteins_list[[length(proteins_list) + 1]] <- list(name = p, sample_type = "Tissue")
        }
      }
    }

    # Extract Plasma-based proteins
    if (has_plasma) {
      plasma_match <- regmatches(markers_text, gregexpr("Plasma-based\\s*\\(([^)]+)\\)", markers_text, ignore.case = TRUE, perl = TRUE))[[1]]
      if (length(plasma_match) > 0) {
        plasma_proteins <- gsub(".*\\(([^)]+)\\).*", "\\1", plasma_match)
        plasma_proteins <- unlist(strsplit(plasma_proteins, ",\\s*"))
        plasma_proteins <- trimws(plasma_proteins)
        plasma_proteins <- gsub("^and\\s+", "", plasma_proteins, ignore.case = TRUE)
        plasma_proteins <- trimws(plasma_proteins)
        plasma_proteins <- plasma_proteins[plasma_proteins != ""]

        for (p in plasma_proteins) {
          proteins_list[[length(proteins_list) + 1]] <- list(name = p, sample_type = "Plasma")
        }
      }
    }

    # Extract Urine-based proteins
    if (has_urine) {
      urine_match <- regmatches(markers_text, gregexpr("Urine-based\\s*\\(([^)]+)\\)", markers_text, ignore.case = TRUE, perl = TRUE))[[1]]
      if (length(urine_match) > 0) {
        urine_proteins <- gsub(".*\\(([^)]+)\\).*", "\\1", urine_match)
        urine_proteins <- unlist(strsplit(urine_proteins, ",\\s*"))
        urine_proteins <- trimws(urine_proteins)
        urine_proteins <- gsub("^and\\s+", "", urine_proteins, ignore.case = TRUE)
        urine_proteins <- trimws(urine_proteins)
        urine_proteins <- urine_proteins[urine_proteins != ""]

        for (p in urine_proteins) {
          proteins_list[[length(proteins_list) + 1]] <- list(name = p, sample_type = "Urine")
        }
      }
    }
  } else {
    # No sample type specified - parse normally
    proteins <- unlist(strsplit(markers_text, ",\\s*"))
    proteins <- trimws(proteins)
    proteins <- gsub("^and\\s+", "", proteins, ignore.case = TRUE)
    proteins <- trimws(proteins)
    proteins <- proteins[proteins != ""]

    for (p in proteins) {
      proteins_list[[length(proteins_list) + 1]] <- list(name = p, sample_type = NA)
    }
  }

  return(proteins_list)
}

# Extract all proteins with their sample types
all_protein_entries <- list()
for (i in 1:nrow(protein_table)) {
  markers <- protein_markers_col[i]
  proteins_in_study <- parse_proteins_with_sample_type(markers, i)
  all_protein_entries[[i]] <- proteins_in_study
}

# Create unique protein identifiers
# Just use protein names without sample type suffix
all_proteins_raw <- lapply(all_protein_entries, function(study) {
  sapply(study, function(p) p$name)
}) %>% unlist()


biom_unique_proteins <- sort(unique(all_proteins_raw))


## -----------------------------------------------------------------------------
# Downloaded the file from DrugBank
df_target <- read.csv(file.path(wd$data, "DrugBank/drugbank_all_target_polypeptide_ids.csv/pharmacologically_active.csv"))
df_vocab <- read.csv(file.path(wd$data, "DrugBank/drugbank vocabulary.csv"))


## -----------------------------------------------------------------------------
df_target_manual <- read.delim(file.path(wd$data, "Drug_review/protein_drug_associations.tsv"),
                               sep = "\t", 
                               stringsAsFactors = FALSE)


## -----------------------------------------------------------------------------
df_target_sub <- df_target %>% select(Gene.Name, Drug.IDs) %>% 
  mutate(Drug.IDs = str_split(Drug.IDs, ",\\s*")) %>%
  tidyr::unnest(Drug.IDs) %>%
  filter(!is.na(Gene.Name), !is.na(Drug.IDs), Drug.IDs != "") %>%
  distinct(Gene.Name, Drug.IDs) %>%
  group_by(Gene.Name) %>%
  summarise(
    DrugBank.ID = str_c(sort(unique(Drug.IDs)), collapse = "; "),
    .groups = "drop"
  )


## -----------------------------------------------------------------------------
df_vocab_sub <- df_vocab %>% select(DrugBank.ID, Common.name)
# Replace each DrugBank ID with the corresponding Common name
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


## -----------------------------------------------------------------------------
# Proteins that are drug target from manual curations
dt_manual <- df_target_manual %>% filter(Drug_Association == "Yes") %>% pull (Protein)

# Final unique list of proteins considered drug-annotated (DrugBank + manual)
drug_annotated_proteins <- sort(unique(c(
  df_target_ann %>% filter(!is.na(DrugBank.ID), DrugBank.ID != "") %>% pull(Gene.Name),
  dt_manual
)))

df_druggable_index <- tibble(
  Protein = drug_annotated_proteins,
  is_Dx = "Yes"
)

cat("Drug-annotated unique proteins:", length(drug_annotated_proteins), "\n")


## -----------------------------------------------------------------------------
# List of cell surface markers
cs <- read.delim(file.path(wd$outData, "CS_human.txt"))
cs_prot <- unique(cs$Approved.symbol)


## -----------------------------------------------------------------------------
df_DEPs <- read.csv(file.path(wd$out06, "HR_values_mCRPC_specific_multivariate_cell_surface.csv"))

# Section 03 all-protein multivariate HR results (2000+ proteins)
df_hr_all <- read.csv(file.path(wd$out06, "HR_values_mCRPC_specific_multi_pre-allProteins.csv"))
df_hr_all <- df_hr_all %>%
  filter(!is.na(Assay), Assay != "") %>%
  distinct(OlinkID, Assay, .keep_all = TRUE)

# Considered prognostic
prog_prot <- df_DEPs %>% filter(Multi_Sig == "Yes") %>% pull(Assay) %>% unique()


## -----------------------------------------------------------------------------
cat("=== Objective 2: CPPS from user-defined protein sets ===\n")

# De-duplicate Assay symbols once (pick one row for duplicated assays)
dep_unique <- df_DEPs %>%
  arrange(Multi_FDR_num, desc(abs(Multi_Coefficient))) %>%
  distinct(Assay, .keep_all = TRUE)

prepare_cpps_protein_table <- function(proteins, method = c("multivariate", "univariate", "equal_weight")) {
  method <- match.arg(method)

  pt <- tibble(Assay = proteins) %>%
    left_join(
      dep_unique %>% select(Assay, OlinkID, Uni_Coefficient, Multi_Coefficient),
      by = "Assay"
    )

  if (method == "multivariate") {
    pt <- pt %>% mutate(Coefficient = Multi_Coefficient)
  } else if (method == "univariate") {
    pt <- pt %>% mutate(Coefficient = Uni_Coefficient)
  } else {
    pt <- pt %>% mutate(Coefficient = 1)
  }

  # Allow proteins without published coefficients by assigning unit weight
  pt <- pt %>%
    mutate(Coefficient = ifelse(is.na(Coefficient), 1, Coefficient)) %>%
    filter(!is.na(OlinkID)) %>%
    distinct(OlinkID, .keep_all = TRUE)

  pt
}

build_cpps_dataset <- function(protein_table) {
  cpps_scores <- df_npx_long_c_pre %>%
    filter(OlinkID %in% protein_table$OlinkID) %>%
    select(HCI_cID, OlinkID, NPX) %>%
    left_join(protein_table %>% select(OlinkID, Coefficient), by = "OlinkID") %>%
    group_by(HCI_cID) %>%
    summarise(CPPS = sum(Coefficient * NPX, na.rm = TRUE), .groups = "drop")

  df_cpps <- meta_mcrpc %>%
    inner_join(cpps_scores, by = "HCI_cID") %>%
    transmute(
      SampleID = HCI_cID,
      CPPS = CPPS,
      log10_psa = log10_psa,
      log10_alp = log10_alp,
      log10_ldh = log10_ldh,
      time_to_death = OS,
      death = death
    ) %>%
    filter(!is.na(CPPS), !is.na(time_to_death), !is.na(death))

  cpps_cut <- median(df_cpps$CPPS, na.rm = TRUE)
  df_cpps %>%
    mutate(
      CPPS_risk = ifelse(CPPS >= cpps_cut, "High", "Low"),
      CPPS_risk = factor(CPPS_risk, levels = c("Low", "High"))
    )
}

run_cpps_model <- function(df_cpps, model_type = c("CPPS_only", "Clinical_plus_CPPS"), eval_times = c(24, 36)) {
  model_type <- match.arg(model_type)
  df_cpps_complete <- df_cpps %>%
    filter(!is.na(log10_psa), !is.na(log10_alp), !is.na(log10_ldh))

  if (model_type == "CPPS_only") {
    fit <- coxph(Surv(time_to_death, death) ~ CPPS, data = df_cpps_complete)
  } else {
    fit <- coxph(Surv(time_to_death, death) ~ log10_psa + log10_alp + log10_ldh + CPPS, data = df_cpps_complete)
  }

  fit_bin <- coxph(Surv(time_to_death, death) ~ CPPS_risk, data = df_cpps_complete)
  lp <- predict(fit, type = "lp")

  roc_obj <- timeROC(
    T = df_cpps_complete$time_to_death,
    delta = df_cpps_complete$death,
    marker = lp,
    cause = 1,
    times = eval_times,
    iid = TRUE
  )
  auc_vals <- if (length(roc_obj$AUC) == length(eval_times) + 1) roc_obj$AUC[2:(length(eval_times) + 1)] else roc_obj$AUC[1:length(eval_times)]

  fit_summary <- summary(fit)
  fit_bin_summary <- summary(fit_bin)

  tibble(
    Model = model_type,
    N = nrow(df_cpps_complete),
    Events = sum(df_cpps_complete$death),
    HR_CPPS_cont = fit_summary$conf.int[nrow(fit_summary$conf.int), "exp(coef)"],
    HR_CPPS_cont_L95 = fit_summary$conf.int[nrow(fit_summary$conf.int), "lower .95"],
    HR_CPPS_cont_U95 = fit_summary$conf.int[nrow(fit_summary$conf.int), "upper .95"],
    P_CPPS_cont = fit_summary$coefficients[nrow(fit_summary$coefficients), "Pr(>|z|)"],
    HR_High_vs_Low = fit_bin_summary$conf.int[1, "exp(coef)"],
    HR_High_vs_Low_L95 = fit_bin_summary$conf.int[1, "lower .95"],
    HR_High_vs_Low_U95 = fit_bin_summary$conf.int[1, "upper .95"],
    P_High_vs_Low = fit_bin_summary$coefficients[1, "Pr(>|z|)"],
    AUC_24m = round(auc_vals[1], 3),
    AUC_36m = round(auc_vals[2], 3)
  )
}

# Example run
example_proteins <- dep_unique %>%
  filter(Multi_Sig == "Yes") %>%
  arrange(desc(Multi_HR)) %>%
  slice_head(n = 10) %>%
  pull(Assay)

cpps_protein_table <- prepare_cpps_protein_table(example_proteins, method = "multivariate")
df_cpps_example <- build_cpps_dataset(cpps_protein_table)
cpps_perf_cpps <- run_cpps_model(df_cpps_example, model_type = "CPPS_only")
cpps_perf_combined <- run_cpps_model(df_cpps_example, model_type = "Clinical_plus_CPPS")
cpps_performance <- bind_rows(cpps_perf_cpps, cpps_perf_combined)

print(cpps_performance)


## -----------------------------------------------------------------------------
# Map OlinkID -> Assay (pick one mapping per assay)
assay_map <- bind_rows(
  dep_unique %>% select(OlinkID, Assay),
  matched_proteins %>% select(OlinkID, Assay)
) %>%
  filter(!is.na(Assay)) %>%
  distinct(OlinkID, .keep_all = TRUE)

fit_olink_hr <- function(proteins) {
  prot_tbl <- tibble(Assay = proteins) %>%
    left_join(assay_map, by = "Assay") %>%
    filter(!is.na(OlinkID)) %>%
    distinct(Assay, .keep_all = TRUE)

  hr_df <- df_npx_long %>%
    filter(cohort == "C", grepl("^Pre", txt_stat_ordered), OlinkID %in% prot_tbl$OlinkID) %>%
    left_join(prot_tbl %>% select(OlinkID, Assay), by = "OlinkID") %>%
    left_join(meta2_olink %>% select(HCI_cID, OS, death, log10_psa, log10_ldh, log10_alp), by = "HCI_cID") %>%
    group_by(Assay, OlinkID) %>%
    group_modify(~{
      dd <- .x %>% filter(!is.na(NPX), !is.na(OS), !is.na(death), !is.na(log10_psa), !is.na(log10_ldh), !is.na(log10_alp))
      if (nrow(dd) < 15 || length(unique(dd$death)) < 2) {
        return(tibble(HR = NA_real_, Lower95 = NA_real_, Upper95 = NA_real_, Pval = NA_real_, N = nrow(dd), Events = sum(dd$death)))
      }
      fit <- coxph(Surv(OS, death) ~ NPX + log10_psa + log10_ldh + log10_alp, data = dd)
      ss <- summary(fit)
      tibble(
        HR = ss$conf.int["NPX", "exp(coef)"],
        Lower95 = ss$conf.int["NPX", "lower .95"],
        Upper95 = ss$conf.int["NPX", "upper .95"],
        Pval = ss$coefficients["NPX", "Pr(>|z|)"],
        N = nrow(dd),
        Events = sum(dd$death)
      )
    }) %>%
    ungroup() %>%
    filter(!is.na(HR))

  hr_df
}

plot_olink_hr_forest <- function(hr_df, title = "Olink HR forest plot (cohort C pre-treatment)") {
  hr_df %>%
    mutate(Assay = forcats::fct_reorder(Assay, HR)) %>%
    ggplot(aes(x = HR, y = Assay)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
    geom_errorbarh(aes(xmin = Lower95, xmax = Upper95), height = 0.18, color = "#2C7FB8") +
    geom_point(size = 2.8, color = "#D95F0E") +
    scale_x_log10() +
    labs(title = title, x = "Hazard Ratio (log scale)", y = "Protein") +
    theme_bw()
}

# Example HR forest with top 15 significant proteins
proteins_for_hr <- dep_unique %>%
  filter(Multi_Sig == "Yes") %>%
  arrange(Multi_FDR_num) %>%
  slice_head(n = 15) %>%
  pull(Assay)

hr_olink <- fit_olink_hr(proteins_for_hr)
p_hr <- plot_olink_hr_forest(hr_olink)
print(p_hr)
ggsave(file.path(wd$outCurr, "Objective3_HR_forest_olink.pdf"), p_hr, width = 8, height = 6)


## -----------------------------------------------------------------------------
plot_expression_cross_platform <- function(proteins) {
  # Define color palettes
  pal.olink <- c("Local" = "#AEC7E8",
                 "mHSPC" = "#FFBB78",
                 "mCRPC" = "#98DF8A")

  # Public data group colors (Benign, Local_PC, mCRPC)
  pal.public <- c("Benign" = "grey50",
                  "Local_PC" = "#AEC7E8",
                  "mCRPC" = "#98DF8A")

  # Olink data (A/B/C, pre-treatment)
  olink_df <- df_npx_long %>%
    filter(cohort %in% c("A", "B", "C"), grepl("^Pre", txt_stat_ordered)) %>%
    left_join(assay_map, by = "OlinkID") %>%
    filter(Assay %in% proteins) %>%
    mutate(
      Group = recode(cohort, A = "Local", B = "mHSPC", C = "mCRPC"),
      Value = NPX
    ) %>%
    select(Assay, Group, Value)

  # Public tissue data (Benign, Local_PC, mCRPC)
  public_df <- as.data.frame(public_norm) %>%
    tibble::rownames_to_column("Protein") %>%
    pivot_longer(-Protein, names_to = "Biological_ID", values_to = "Value") %>%
    left_join(
      public_avg %>% select(Biological_ID, Group) %>% distinct(),
      by = "Biological_ID"
    ) %>%
    filter(Group %in% c("Benign", "Local_PC", "mCRPC")) %>%
    left_join(matched_proteins %>% select(Protein, Assay) %>% distinct(), by = "Protein") %>%
    filter(Assay %in% proteins) %>%
    select(Assay, Group, Value)

  blank_panel <- function(title_txt) {
    ggplot() +
      annotate("text", x = 1, y = 1, label = "No data", color = "grey45", size = 3.5) +
      xlim(0.5, 1.5) + ylim(0.5, 1.5) +
      labs(title = title_txt) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
        panel.border = element_rect(color = "grey80", fill = NA)
      )
  }

  row_plots <- lapply(proteins, function(p) {
    d_ol <- olink_df %>% filter(Assay == p)
    d_pb <- public_df %>% filter(Assay == p)

    p_ol <- if (nrow(d_ol) > 0) {
      d_ol %>%
        mutate(Group = factor(Group, levels = names(pal.olink))) %>%
        ggplot(aes(x = Group, y = Value, fill = Group)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.75) +
        geom_jitter(width = 0.18, size = 0.9, alpha = 0.5) +
        scale_fill_manual(values = pal.olink, drop = FALSE) +
        labs(title = paste0(p, " | Olink"), x = NULL, y = "NPX") +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          plot.title = element_text(size = 10, face = "bold")
        )
    } else {
      blank_panel(paste0(p, " | Olink"))
    }

    p_pb <- if (nrow(d_pb) > 0) {
      d_pb %>%
        mutate(Group = factor(Group, levels = names(pal.public))) %>%
        ggplot(aes(x = Group, y = Value, fill = Group)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.75) +
        geom_jitter(width = 0.18, size = 0.9, alpha = 0.5) +
        scale_fill_manual(values = pal.public, drop = FALSE) +
        labs(title = paste0(p, " | Public"), x = NULL, y = "log2 intensity") +
        theme_bw() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          plot.title = element_text(size = 10, face = "bold")
        )
    } else {
      blank_panel(paste0(p, " | Public"))
    }

    p_ol + p_pb + patchwork::plot_layout(ncol = 2, widths = c(1, 1))
  })

  patchwork::wrap_plots(row_plots, ncol = 1) +
    patchwork::plot_annotation(title = "Cross-platform expression by protein (Left: Olink, Right: Public)")
}

# Example expression plot with matched prognostic proteins
proteins_for_expr <- matched_proteins %>%
  filter(Assay %in% prog_prot) %>%
  distinct(Assay) %>%
  slice_head(n = 6) %>%
  pull(Assay)

p_expr <- plot_expression_cross_platform(proteins_for_expr)
print(p_expr)
ggsave(file.path(wd$outCurr, "Objective3_Expression_cross_platform.pdf"), p_expr, width = 14, height = 2.6 * length(proteins_for_expr))


## -----------------------------------------------------------------------------
build_integration_panels_gg <- function(proteins, cpps_method = c("multivariate", "univariate", "equal_weight")) {
  cpps_method <- match.arg(cpps_method)
  proteins <- unique(proteins)
  if (length(proteins) == 0) stop("No proteins provided.")

  # Define scaling function (min-max per row)
  pheatmap.scale <- function(x) {
    rmin <- apply(x, 1, min, na.rm = TRUE)
    rmax <- apply(x, 1, max, na.rm = TRUE)
    rrng <- rmax - rmin
    z <- (x - rmin) / rrng
    z[!is.finite(z)] <- NA_real_
    z
  }

  hr_all_tbl <- if (exists("df_hr_all", inherits = TRUE)) {
    get("df_hr_all", inherits = TRUE)
  } else {
    read.csv(file.path(wd$out06, "HR_values_mCRPC_specific_multi_pre-allProteins.csv")) %>%
      filter(!is.na(Assay), Assay != "") %>%
      distinct(OlinkID, Assay, .keep_all = TRUE)
  }

  # Order proteins by HR (high -> low), keep CPPS as top row
  hr_order_df <- tibble(Assay = proteins, input_idx = seq_along(proteins)) %>%
    mutate(Assay_key = toupper(Assay)) %>%
    left_join(
      hr_all_tbl %>%
        transmute(Assay_key = toupper(Assay), HR = as.numeric(HR)),
      by = "Assay_key"
    ) %>%
    arrange(is.na(HR), desc(HR), input_idx)
  proteins_ord <- hr_order_df$Assay
  rows_all <- c(proteins_ord, "CPPS")
  prot_levels <- c(rev(proteins_ord), "CPPS")

  ann_df <- tibble(
    Protein = rows_all,
    `In Olink` = ifelse(Protein == "CPPS", NA_integer_, as.integer(Protein %in% unique(assay_map$Assay))),
    `In Public` = ifelse(Protein == "CPPS", NA_integer_, as.integer(Protein %in% unique(matched_proteins$Assay))),
    `Biom` = ifelse(Protein == "CPPS", NA_integer_, as.integer(toupper(Protein) %in% toupper(unique(biom_unique_proteins)))),
    `CS` = ifelse(Protein == "CPPS", NA_integer_, as.integer(toupper(Protein) %in% toupper(unique(cs_prot)))),
    `Dx` = ifelse(Protein == "CPPS", NA_integer_, as.integer(toupper(Protein) %in% toupper(unique(drug_annotated_proteins)))),
    `Prog` = ifelse(Protein == "CPPS", NA_integer_, as.integer(Protein %in% unique(dep_unique$Assay[dep_unique$Multi_Sig == "Yes"])))
  ) %>%
    mutate(Protein = factor(Protein, levels = prot_levels))

  ann_cells <- ann_df %>%
    select(-Protein) %>%
    mutate(Protein = ann_df$Protein) %>%
    pivot_longer(-Protein, names_to = "Column", values_to = "Value") %>%
    mutate(
      Column = factor(Column, levels = c("In Public", "Biom", "CS", "Dx", "Prog")),
      fill_key = case_when(
        is.na(Value) ~ NA_character_,
        Value == 0 ~ "no",
        Column == "In Public" & Value == 1 ~ "plasma_yes",
        Value == 1 ~ "other_yes",
        TRUE ~ NA_character_
      )
    )

  name_cells <- ann_df %>%
    transmute(Protein, Column = factor("Protein", levels = c("Protein", "In Public", "Biom", "CS", "Dx", "Prog")), label = as.character(Protein))

  cpps_sum_public <- ann_df %>%
    filter(as.character(Protein) != "CPPS") %>%
    summarise(score = sum(`In Public`, na.rm = TRUE)) %>%
    pull(score)
  cpps_sum_biom <- ann_df %>%
    filter(as.character(Protein) != "CPPS") %>%
    summarise(score = sum(`Biom`, na.rm = TRUE)) %>%
    pull(score)
  cpps_sum_cs <- ann_df %>%
    filter(as.character(Protein) != "CPPS") %>%
    summarise(score = sum(`CS`, na.rm = TRUE)) %>%
    pull(score)
  cpps_sum_dx <- ann_df %>%
    filter(as.character(Protein) != "CPPS") %>%
    summarise(score = sum(`Dx`, na.rm = TRUE)) %>%
    pull(score)
  cpps_sum_prog <- ann_df %>%
    filter(as.character(Protein) != "CPPS") %>%
    summarise(score = sum(`Prog`, na.rm = TRUE)) %>%
    pull(score)

  cpps_score_cells <- tibble(
    Protein = factor(rep("CPPS", 5), levels = prot_levels),
    Column = factor(c("In Public", "Biom", "CS", "Dx", "Prog"), levels = c("Protein", "In Public", "Biom", "CS", "Dx", "Prog")),
    label = c(
      as.character(cpps_sum_public),
      as.character(cpps_sum_biom),
      as.character(cpps_sum_cs),
      as.character(cpps_sum_dx),
      as.character(cpps_sum_prog)
    )
  )

  p_ann <- ggplot() +
    geom_tile(
      data = ann_cells,
      aes(x = Column, y = Protein, fill = fill_key),
      color = "grey60", linewidth = 0.3, width = 0.95, height = 0.95
    ) +
    geom_tile(
      data = name_cells,
      aes(x = Column, y = Protein),
      fill = "white", color = "grey60", linewidth = 0.3, width = 0.95, height = 0.95
    ) +
    geom_text(
      data = name_cells,
      aes(x = Column, y = Protein, label = label),
      size = 3
    ) +
    geom_text(
      data = cpps_score_cells,
      aes(x = Column, y = Protein, label = label),
      size = 3.2,
      fontface = "bold"
    ) +
    scale_fill_manual(
      values = c(
        "no" = "grey99",
        "olink_yes" = "#D69526",
        "plasma_yes" = "#FBB4AE",
        "other_yes" = "grey40"
      ),
      na.value = "white",
      guide = "none"
    ) +
    scale_x_discrete(drop = FALSE, limits = c("Protein", "In Public", "Biom", "CS", "Dx", "Prog")) +
    labs(x = NULL, y = NULL, title = "Annotations") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      panel.grid = element_blank(),
      plot.title = element_text(size = 10, face = "bold")
    )

  olink_long <- df_npx_long %>%
    filter(cohort %in% c("A", "C"), grepl("^Pre", txt_stat_ordered)) %>%
    left_join(assay_map, by = "OlinkID") %>%
    filter(!is.na(Assay), Assay %in% proteins) %>%
    mutate(Group = recode(cohort, A = "Local", C = "mCRPC")) %>%
    transmute(Assay, Sample = HCI_cID, Group, Expr = NPX)

  # min-max scaling within each protein for Olink panel (row-wise matrix scaling)
  olink_mat <- olink_long %>%
    select(Assay, Sample, Expr) %>%
    distinct() %>%
    pivot_wider(names_from = Sample, values_from = Expr)
  olink_mat_num <- as.matrix(olink_mat %>% select(-Assay))
  rownames(olink_mat_num) <- olink_mat$Assay
  olink_z_mat <- pheatmap.scale(olink_mat_num)
  olink_z_long <- as.data.frame(olink_z_mat) %>%
    tibble::rownames_to_column("Assay") %>%
    pivot_longer(-Assay, names_to = "Sample", values_to = "Expr_scaled") %>%
    left_join(olink_long %>% select(Sample, Group) %>% distinct(), by = "Sample")

  olink_z_mean <- olink_z_long %>%
    group_by(Assay, Group) %>%
    summarise(MeanExpr = mean(Expr_scaled, na.rm = TRUE), .groups = "drop") %>%
    complete(Assay = proteins, Group = c("Local", "mCRPC"), fill = list(MeanExpr = NA_real_)) %>%
    mutate(Assay = as.character(Assay)) %>%
    bind_rows(tibble(Assay = "CPPS", Group = c("Local", "mCRPC"), MeanExpr = NA_real_)) %>%
    mutate(Assay = factor(Assay, levels = prot_levels), Group = factor(Group, levels = c("Local", "mCRPC")))

  p_olink_mean <- olink_z_mean %>%
    ggplot(aes(y = Assay)) +
    geom_segment(
      data = ~ dplyr::filter(.x, Group == "Local"),
      aes(x = 0, xend = MeanExpr, yend = Assay),
      color = "grey40",
      position = position_nudge(y = -0.16),
      linewidth = 0.55,
      na.rm = TRUE
    ) +
    geom_point(
      data = ~ dplyr::filter(.x, Group == "Local"),
      aes(x = MeanExpr),
      color = pal.olink["Local"],
      position = position_nudge(y = -0.16),
      size = 2.5,
      na.rm = TRUE
    ) +
    geom_segment(
      data = ~ dplyr::filter(.x, Group == "mCRPC"),
      aes(x = 0, xend = MeanExpr, yend = Assay),
      color = "grey40",
      position = position_nudge(y = 0.16),
      linewidth = 0.55,
      na.rm = TRUE
    ) +
    geom_point(
      data = ~ dplyr::filter(.x, Group == "mCRPC"),
      aes(x = MeanExpr),
      color = pal.olink["mCRPC"],
      position = position_nudge(y = 0.16),
      size = 2.5,
      na.rm = TRUE
    ) +
    #scale_x_continuous(limits = c(0, 1)) +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Mean min-max scaled", y = NULL, title = "Olink Expression") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      panel.grid = element_blank(),
      plot.title = element_text(size = 10, face = "bold")
    )

  # FC is kept on raw scale
  olink_fc <- olink_long %>%
    group_by(Assay, Group) %>%
    summarise(MeanExpr = mean(Expr, na.rm = TRUE), .groups = "drop") %>%
    complete(Assay = proteins, Group = c("Local", "mCRPC"), fill = list(MeanExpr = NA_real_)) %>%
    pivot_wider(names_from = Group, values_from = MeanExpr) %>%
    mutate(FC = mCRPC - Local) %>%
    transmute(Assay, Dataset = "Olink", log2FC = FC) %>%
    bind_rows(tibble(Assay = "CPPS", Dataset = "Olink", log2FC = NA_real_))

  p_olink_fc <- olink_fc %>%
    mutate(Assay = factor(Assay, levels = prot_levels), Dataset = factor(Dataset, levels = "Olink")) %>%
    ggplot(aes(x = Dataset, y = Assay, fill = log2FC)) +
    geom_tile(color = "white", linewidth = 0.3, width = 0.95, height = 0.95, na.rm = TRUE) +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0, oob = scales::squish, na.value = "white") +
    scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = "Olink FC") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 10, face = "bold", angle = 30, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank(),
      #panel.border = element_blank(),
      plot.title = element_text(size = 10, face = "bold")
    )

  public_long <- as.data.frame(public_norm) %>%
    tibble::rownames_to_column("Protein") %>%
    pivot_longer(-Protein, names_to = "Biological_ID", values_to = "Value") %>%
    left_join(public_avg %>% select(Biological_ID, Group) %>% distinct(), by = "Biological_ID") %>%
    filter(Group %in% c("Local_PC", "mCRPC")) %>%
    left_join(matched_proteins %>% select(Protein, Assay) %>% distinct(), by = "Protein") %>%
    filter(!is.na(Assay), Assay %in% proteins) %>%
    transmute(Assay, Sample = Biological_ID, Group, Expr = Value)

  # min-max scaling within each protein for Public panel (row-wise matrix scaling)
  public_mat <- public_long %>%
    select(Assay, Sample, Expr) %>%
    distinct() %>%
    pivot_wider(names_from = Sample, values_from = Expr)
  public_mat_num <- as.matrix(public_mat %>% select(-Assay))
  rownames(public_mat_num) <- public_mat$Assay
  public_z_mat <- pheatmap.scale(public_mat_num)
  public_z_long <- as.data.frame(public_z_mat) %>%
    tibble::rownames_to_column("Assay") %>%
    pivot_longer(-Assay, names_to = "Sample", values_to = "Expr_scaled") %>%
    left_join(public_long %>% select(Sample, Group) %>% distinct(), by = "Sample")

  public_z_mean <- public_z_long %>%
    group_by(Assay, Group) %>%
    summarise(MeanExpr = mean(Expr_scaled, na.rm = TRUE), .groups = "drop") %>%
    complete(Assay = proteins, Group = c("Local_PC", "mCRPC"), fill = list(MeanExpr = NA_real_)) %>%
    mutate(Assay = as.character(Assay)) %>%
    bind_rows(tibble(Assay = "CPPS", Group = c("Local_PC", "mCRPC"), MeanExpr = NA_real_)) %>%
    mutate(Assay = factor(Assay, levels = prot_levels), Group = factor(Group, levels = c("Local_PC", "mCRPC")))

  p_public_mean <- public_z_mean %>%
    ggplot(aes(y = Assay)) +
    geom_segment(
      data = ~ dplyr::filter(.x, Group == "Local_PC"),
      aes(x = 0, xend = MeanExpr, yend = Assay),
      color = "grey40",
      position = position_nudge(y = -0.16),
      linewidth = 0.55,
      na.rm = TRUE
    ) +
    geom_point(
      data = ~ dplyr::filter(.x, Group == "Local_PC"),
      aes(x = MeanExpr),
      color = pal.public["Local_PC"],
      position = position_nudge(y = -0.16),
      size = 2.5,
      na.rm = TRUE
    ) +
    geom_segment(
      data = ~ dplyr::filter(.x, Group == "mCRPC"),
      aes(x = 0, xend = MeanExpr, yend = Assay),
      color = "grey40",
      position = position_nudge(y = 0.16),
      linewidth = 0.55,
      na.rm = TRUE
    ) +
    geom_point(
      data = ~ dplyr::filter(.x, Group == "mCRPC"),
      aes(x = MeanExpr),
      color = pal.public["mCRPC"],
      position = position_nudge(y = 0.16),
      size = 2.5,
      na.rm = TRUE
    ) +
    # scale_x_continuous(limits = c(0, 1)) +
    scale_x_reverse() +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Mean min-max scaled", y = NULL, title = "Public Expression") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.position = "none",
      panel.grid = element_blank(),
      plot.title = element_text(size = 10, face = "bold")
    )

  # FC is kept on raw normalized scale
  public_fc <- public_long %>%
    group_by(Assay, Group) %>%
    summarise(MeanExpr = mean(Expr, na.rm = TRUE), .groups = "drop") %>%
    complete(Assay = proteins, Group = c("Local_PC", "mCRPC"), fill = list(MeanExpr = NA_real_)) %>%
    pivot_wider(names_from = Group, values_from = MeanExpr) %>%
    mutate(FC = mCRPC - Local_PC) %>%
    transmute(Assay, Dataset = "Public", log2FC = FC) %>%
    bind_rows(tibble(Assay = "CPPS", Dataset = "Public", log2FC = NA_real_))

  p_public_fc <- public_fc %>%
    mutate(Assay = factor(Assay, levels = prot_levels), Dataset = factor(Dataset, levels = "Public")) %>%
    ggplot(aes(x = Dataset, y = Assay, fill = log2FC)) +
    geom_tile(color = "white", linewidth = 0.3, width = 0.95, height = 0.95, na.rm = TRUE) +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0, oob = scales::squish, na.value = "white") +
    scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = "Public FC") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 10, face = "bold", angle = 30, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank(),
      #panel.border = element_blank(),
      plot.title = element_text(size = 10, face = "bold")
    )

  hr_df <- tibble(Assay = proteins) %>%
    mutate(Assay_key = toupper(Assay)) %>%
    left_join(
      hr_all_tbl %>%
        transmute(
          Assay_key = toupper(Assay),
          HR = as.numeric(HR),
          L95 = as.numeric(Lower_Bound),
          U95 = as.numeric(Upper_Bound)
        ),
      by = "Assay_key"
    )

  cpps_tbl <- prepare_cpps_protein_table(proteins, method = cpps_method)
  cpps_df <- build_cpps_dataset(cpps_tbl) %>%
    filter(!is.na(CPPS), !is.na(time_to_death), !is.na(death))

  cpps_hr <- NA_real_
  cpps_l95 <- NA_real_
  cpps_u95 <- NA_real_
  if (nrow(cpps_df) >= 15 && length(unique(cpps_df$death)) > 1) {
    cpps_cut_cox <- median(cpps_df$CPPS, na.rm = TRUE)
    cpps_df <- cpps_df %>%
      mutate(
        CPPS_risk = ifelse(CPPS >= cpps_cut_cox, "High", "Low"),
        CPPS_risk = factor(CPPS_risk, levels = c("Low", "High"))
      )
    fit_cpps <- coxph(Surv(time_to_death, death) ~ CPPS_risk, data = cpps_df)
    ss_cpps <- summary(fit_cpps)
    cpps_hr <- ss_cpps$conf.int["CPPS_riskHigh", "exp(coef)"]
    cpps_l95 <- ss_cpps$conf.int["CPPS_riskHigh", "lower .95"]
    cpps_u95 <- ss_cpps$conf.int["CPPS_riskHigh", "upper .95"]
  }

  hr_df <- bind_rows(
    hr_df,
    tibble(Assay = "CPPS", HR = cpps_hr, L95 = cpps_l95, U95 = cpps_u95)
  ) %>%
    mutate(
      Assay = factor(Assay, levels = prot_levels),
      log2HR = log2(HR),
      log2L95 = log2(L95),
      log2U95 = log2(U95)
    )

  hr_xmin <- suppressWarnings(min(hr_df$log2L95, na.rm = TRUE))
  hr_xmax <- suppressWarnings(max(hr_df$log2U95, na.rm = TRUE))
  if (!is.finite(hr_xmin)) hr_xmin <- -1
  if (!is.finite(hr_xmax)) hr_xmax <- 1
  hr_xlim <- c(min(-1, hr_xmin), max(1, hr_xmax))

  p_hr <- hr_df %>%
    ggplot(aes(y = Assay)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey55") +
    geom_segment(
      data = hr_df %>% filter(is.finite(log2L95), is.finite(log2U95)),
      aes(x = log2L95, xend = log2U95, y = Assay, yend = Assay),
      linewidth = 0.5
    ) +
    geom_point(
      data = hr_df %>% filter(is.finite(log2HR)),
      aes(x = log2HR),
      size = 2.6
    ) +
    coord_cartesian(xlim = hr_xlim) +
    labs(x = "log2(HR)", y = NULL, title = "HR (OS)") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 10, face = "bold")
    )

  p_ann + p_public_mean + p_public_fc + p_olink_fc + p_olink_mean + p_hr +
    patchwork::plot_layout(ncol = 6, widths = c(3, 2, 0.6, 0.6, 2, 2.5))
}

# Example usage
# proteins_for_complex <- c("SDC1", "EDF1", "CKAP4", "LMNB2", "LETM1", "COX5B", "HIP1R", "HMOX2", "DPY30", "WFDC2", "RCC1")
#proteins_for_complex <- c("FETUB", "MOCS2", "MSLN", "USP28", "SEZ6L2", "TOMM20", "CCN4", "CXADR", "PAXX", "MUC13")
proteins_for_complex <- c("DPY30", "NFU1", "LETM1", "GDF15", "TPR", "CST7", "FKBP5", "LMNB2", "SDC1", "KRT19")
p_multi <- build_integration_panels_gg(proteins_for_complex)
print(p_multi)
ggsave(file.path(wd$outCurr, "Objective3_Integration_Multipanel_ggplot.pdf"), 
       p_multi, width = 10.7, height = 0.4 * length(proteins_for_complex) + 3)

