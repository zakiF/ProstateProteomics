normalize_sample_id <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("-", "_", x)
  x <- gsub("^c(?=[0-9])", "c_", x, perl = TRUE)
  x <- gsub("^c([0-9]{4}_[0-9]+)$", "c_\\1", x)
  x
}

parse_proteins <- function(txt) {
  txt %>%
    strsplit("[,;[:space:]]+") %>%
    unlist() %>%
    trimws() %>%
    .[nzchar(.)] %>%
    unique()
}

extract_auc_values <- function(roc_obj, eval_times) {
  if (length(roc_obj$AUC) == length(eval_times) + 1) {
    return(as.numeric(roc_obj$AUC[2:(length(eval_times) + 1)]))
  }
  as.numeric(roc_obj$AUC[1:length(eval_times)])
}

prepare_cpps_protein_table <- function(proteins, protein_stats, method = c("multivariate", "univariate", "equal_weight")) {
  method <- match.arg(method)

  pt <- tibble(Assay = proteins) %>%
    left_join(
      protein_stats %>% select(Assay, OlinkID, Uni_Coefficient, Multi_Coefficient),
      by = "Assay"
    )

  if (method == "multivariate") {
    pt <- pt %>% mutate(Coefficient = Multi_Coefficient)
  } else if (method == "univariate") {
    pt <- pt %>% mutate(Coefficient = Uni_Coefficient)
  } else {
    pt <- pt %>% mutate(Coefficient = 1)
  }

  pt %>%
    mutate(Coefficient = ifelse(is.na(Coefficient), 1, Coefficient)) %>%
    filter(!is.na(OlinkID)) %>%
    distinct(OlinkID, .keep_all = TRUE)
}

build_cpps_dataset <- function(protein_table, prep) {
  npx_src <- if (!is.null(prep$npx_long_c_pre)) prep$npx_long_c_pre else prep$npx_long
  cpps_scores <- npx_src %>%
    filter(OlinkID %in% protein_table$OlinkID) %>%
    select(sample_id_norm, OlinkID, NPX) %>%
    left_join(protein_table %>% select(OlinkID, Coefficient), by = "OlinkID") %>%
    group_by(sample_id_norm) %>%
    summarise(CPPS = sum(Coefficient * NPX, na.rm = TRUE), .groups = "drop")

  df_cpps <- prep$analysis_meta %>%
    select(sample_id_norm, OS, death, log10_psa, log10_alp, log10_ldh) %>%
    inner_join(cpps_scores, by = "sample_id_norm") %>%
    transmute(
      SampleID = sample_id_norm,
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

build_integration_panels_gg <- function(proteins, prep, cpps_method = c("multivariate", "univariate", "equal_weight")) {
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

  pal.olink <- c("Local" = "#AEC7E8", "mCRPC" = "#98DF8A")
  pal.public <- c("Local_PC" = "#AEC7E8", "mCRPC" = "#98DF8A")

  hr_all_tbl <- if (!is.null(prep$df_hr_all)) prep$df_hr_all else prep$protein_stats
  if (!"HR" %in% names(hr_all_tbl) && "Multi_HR" %in% names(hr_all_tbl)) {
    hr_all_tbl <- hr_all_tbl %>% mutate(HR = Multi_HR)
  }
  if (!"Lower_Bound" %in% names(hr_all_tbl)) hr_all_tbl$Lower_Bound <- NA_real_
  if (!"Upper_Bound" %in% names(hr_all_tbl)) hr_all_tbl$Upper_Bound <- NA_real_

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

  dep_unique <- if (!is.null(prep$dep_unique)) prep$dep_unique else prep$protein_stats
  biom_unique_proteins <- if (!is.null(prep$biom_unique_proteins)) prep$biom_unique_proteins else character(0)
  cs_prot <- if (!is.null(prep$cs_prot)) prep$cs_prot else character(0)
  drug_annotated_proteins <- if (!is.null(prep$drug_annotated_proteins)) prep$drug_annotated_proteins else character(0)

  ann_df <- tibble(
    Protein = rows_all,
    `In Tissue` = ifelse(Protein == "CPPS", NA_integer_, as.integer(Protein %in% unique(prep$matched_proteins$Assay))),
    `Biomarker` = ifelse(Protein == "CPPS", NA_integer_, as.integer(toupper(Protein) %in% toupper(unique(biom_unique_proteins)))),
    `Surface Protein` = ifelse(Protein == "CPPS", NA_integer_, as.integer(toupper(Protein) %in% toupper(unique(cs_prot)))),
    `Druggable` = ifelse(Protein == "CPPS", NA_integer_, as.integer(toupper(Protein) %in% toupper(unique(drug_annotated_proteins))))
  ) %>%
    mutate(Protein = factor(Protein, levels = prot_levels))

  ann_cells <- ann_df %>%
    select(-Protein) %>%
    mutate(Protein = ann_df$Protein) %>%
    pivot_longer(-Protein, names_to = "Column", values_to = "Value") %>%
    mutate(
      Column = factor(Column, levels = c("In Tissue", "Biomarker", "Surface Protein", "Druggable")),
      fill_key = case_when(
        is.na(Value) ~ NA_character_,
        Value == 0 ~ "no",
        Column == "In Tissue" & Value == 1 ~ "plasma_yes",
        Value == 1 ~ "other_yes",
        TRUE ~ NA_character_
      )
    )

  name_cells <- ann_df %>%
    transmute(Protein, Column = factor("Protein", levels = c("Protein", "In Tissue", "Biomarker", "Surface Protein", "Druggable")), label = as.character(Protein))

  cpps_sum_tissue <- ann_df %>% filter(as.character(Protein) != "CPPS") %>% summarise(score = sum(`In Tissue`, na.rm = TRUE)) %>% pull(score)
  cpps_sum_biom <- ann_df %>% filter(as.character(Protein) != "CPPS") %>% summarise(score = sum(`Biomarker`, na.rm = TRUE)) %>% pull(score)
  cpps_sum_cs <- ann_df %>% filter(as.character(Protein) != "CPPS") %>% summarise(score = sum(`Surface Protein`, na.rm = TRUE)) %>% pull(score)
  cpps_sum_dx <- ann_df %>% filter(as.character(Protein) != "CPPS") %>% summarise(score = sum(`Druggable`, na.rm = TRUE)) %>% pull(score)

  cpps_score_cells <- tibble(
    Protein = factor(rep("CPPS", 4), levels = prot_levels),
    Column = factor(c("In Tissue", "Biomarker", "Surface Protein", "Druggable"), levels = c("Protein", "In Tissue", "Biomarker", "Surface Protein", "Druggable")),
    label = c(as.character(cpps_sum_tissue), as.character(cpps_sum_biom), as.character(cpps_sum_cs), as.character(cpps_sum_dx))
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
      size = 3.5
    ) +
    geom_text(
      data = cpps_score_cells,
      aes(x = Column, y = Protein, label = label),
      size = 3.8,
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
    scale_x_discrete(drop = FALSE, limits = c("Protein", "In Tissue", "Biomarker", "Surface Protein", "Druggable")) +
    labs(x = NULL, y = NULL, title = "Annotations") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 10, face = "bold")
    )

  olink_long <- prep$npx_long %>%
    filter(cohort %in% c("A", "C"), grepl("^Pre", txt_stat_ordered)) %>%
    filter(!is.na(Assay), Assay %in% proteins) %>%
    mutate(Group = recode(cohort, A = "Local", C = "mCRPC")) %>%
    transmute(Assay, Sample = sample_id_norm, Group, Expr = NPX)

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
    scale_y_discrete(drop = FALSE) +
    labs(x = "Scaled Exp", y = NULL, title = "Plasma Exp") +
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
    transmute(Assay, Dataset = "Plasma", log2FC = FC) %>%
    bind_rows(tibble(Assay = "CPPS", Dataset = "Plasma", log2FC = NA_real_))

  p_olink_fc <- olink_fc %>%
    mutate(Assay = factor(Assay, levels = prot_levels), Dataset = factor(Dataset, levels = "Plasma")) %>%
    ggplot(aes(x = Dataset, y = Assay, fill = log2FC)) +
    geom_tile(color = "white", linewidth = 0.3, width = 0.95, height = 0.95, na.rm = TRUE) +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0, oob = scales::squish, na.value = "white") +
    scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = "FC") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 10, face = "bold", angle = 30, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank(),
      plot.title = element_text(size = 10, face = "bold")
    )

  public_long <- as.data.frame(prep$public_norm) %>%
    tibble::rownames_to_column("Protein") %>%
    pivot_longer(-Protein, names_to = "Biological_ID", values_to = "Value") %>%
    left_join(prep$public_avg %>% select(Biological_ID, Group) %>% distinct(), by = "Biological_ID") %>%
    filter(Group %in% c("Local_PC", "mCRPC")) %>%
    left_join(prep$matched_proteins %>% select(Protein, Assay) %>% distinct(), by = "Protein") %>%
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
    scale_x_reverse() +
    scale_y_discrete(drop = FALSE) +
    labs(x = "Scaled Exp", y = NULL, title = "Tissue Exp") +
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
    transmute(Assay, Dataset = "Tissue", log2FC = FC) %>%
    bind_rows(tibble(Assay = "CPPS", Dataset = "Tissue", log2FC = NA_real_))

  p_public_fc <- public_fc %>%
    mutate(Assay = factor(Assay, levels = prot_levels), Dataset = factor(Dataset, levels = "Tissue")) %>%
    ggplot(aes(x = Dataset, y = Assay, fill = log2FC)) +
    geom_tile(color = "white", linewidth = 0.3, width = 0.95, height = 0.95, na.rm = TRUE) +
    scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0, oob = scales::squish, na.value = "white") +
    scale_y_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL, title = "FC") +
    theme_bw() +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 10, face = "bold", angle = 30, hjust = 1),
      legend.position = "none",
      panel.grid = element_blank(),
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

  cpps_tbl <- prepare_cpps_protein_table(proteins, prep$protein_stats, method = cpps_method)
  cpps_df <- build_cpps_dataset(cpps_tbl, prep) %>%
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
    fit_cpps <- survival::coxph(survival::Surv(time_to_death, death) ~ CPPS_risk, data = cpps_df)
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
build_analysis <- function(protein_symbols, prep) {
  stats_sub <- prep$protein_stats %>%
    filter(Assay %in% protein_symbols) %>%
    distinct(Assay, .keep_all = TRUE)

  shiny::validate(shiny::need(nrow(stats_sub) > 0, "No selected proteins with available coefficients."))

  cpps_tbl <- prepare_cpps_protein_table(stats_sub$Assay, prep$protein_stats, method = "multivariate")
  df <- build_cpps_dataset(cpps_tbl, prep)

  shiny::validate(shiny::need(nrow(df) > 10, "Too few valid samples after filtering."))

  df_complete <- df %>% filter(!is.na(log10_psa), !is.na(log10_alp), !is.na(log10_ldh))

  cox_combined <- survival::coxph(survival::Surv(time_to_death, death) ~ log10_psa + log10_alp + log10_ldh + CPPS, data = df_complete)

  # Constant clinical-only model (independent of selected proteins)
  clinical_ids <- prep$analysis_meta %>% distinct(sample_id_norm) %>% pull(sample_id_norm)
  clinical_df <- prep$analysis_meta %>%
    filter(sample_id_norm %in% clinical_ids) %>%
    transmute(
      time_to_death = OS,
      death = death,
      log10_psa = log10_psa,
      log10_alp = log10_alp,
      log10_ldh = log10_ldh
    ) %>%
    filter(!is.na(time_to_death), !is.na(death), !is.na(log10_psa), !is.na(log10_alp), !is.na(log10_ldh))
  cox_clinical <- survival::coxph(survival::Surv(time_to_death, death) ~ log10_psa + log10_alp + log10_ldh, data = clinical_df)

  roc_clinical <- timeROC::timeROC(
    T = clinical_df$time_to_death,
    delta = clinical_df$death,
    marker = predict(cox_clinical, type = "lp"),
    cause = 1,
    times = prep$time_points,
    iid = TRUE
  )

  roc_combined <- timeROC::timeROC(
    T = df_complete$time_to_death,
    delta = df_complete$death,
    marker = predict(cox_combined, type = "lp"),
    cause = 1,
    times = prep$time_points,
    iid = TRUE
  )

  auc_tbl <- tibble(
    time_months = prep$time_points,
    AUC_Clinical = extract_auc_values(roc_clinical, prep$time_points),
    AUC_Clinical_CPPS = extract_auc_values(roc_combined, prep$time_points)
  )

  km_fit <- survival::survfit(survival::Surv(time_to_death, death) ~ CPPS_risk, data = df)
  km_plot <- survminer::ggsurvplot(
    km_fit,
    data = df,
    pval = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.28,
    palette = c("grey", "#e41a1c"),
    legend.labs = c("Low Risk", "High Risk"),
    title = "CPPS Kaplan-Meier",
    xlab = "Time (months)",
    ylab = "Overall survival probability"
  )

  prot_km <- lapply(stats_sub$Assay, function(gene) {
    npx_src <- if (!is.null(prep$npx_long_c_pre)) prep$npx_long_c_pre else prep$npx_long
    pd <- npx_src %>%
      filter(Assay == gene) %>%
      select(sample_id_norm, NPX) %>%
      group_by(sample_id_norm) %>%
      summarise(NPX = mean(NPX, na.rm = TRUE), .groups = "drop") %>%
      left_join(prep$analysis_meta %>% select(sample_id_norm, OS, death), by = "sample_id_norm") %>%
      filter(!is.na(NPX), !is.na(OS), !is.na(death)) %>%
      mutate(group = ifelse(NPX >= median(NPX, na.rm = TRUE), "High", "Low"), group = factor(group, levels = c("Low", "High")))

    if (nrow(pd) < 8) return(NULL)
    fit <- survival::survfit(survival::Surv(OS, death) ~ group, data = pd)
    gp <- survminer::ggsurvplot(
      fit, data = pd, pval = TRUE, conf.int = FALSE, risk.table = FALSE,
      palette = c("grey", "#e41a1c"), legend = "none", title = gene,
      xlab = "Months", ylab = "OS"
    )
    gp$plot
  })
  prot_km <- Filter(Negate(is.null), prot_km)

  integration_panel <- build_integration_panels_gg(stats_sub$Assay, prep = prep, cpps_method = "multivariate")

  list(
    proteins = stats_sub$Assay,
    stats_sub = stats_sub,
    df = df,
    df_complete = df_complete,
    km_plot = km_plot,
    protein_km_plots = prot_km,
    auc_tbl = auc_tbl,
    roc_clinical = roc_clinical,
    roc_combined = roc_combined,
    summary_plot = integration_panel
  )
}
