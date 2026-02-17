build_complex_protein_heatmap <- function(proteins, cpps_method = c("multivariate", "univariate", "equal_weight")) {
  cpps_method <- match.arg(cpps_method)
  
  proteins <- unique(proteins)
  if (length(proteins) == 0) stop("No proteins provided.")
  
  # Allow running this chunk standalone by loading HR table on demand
  hr_all_tbl <- if (exists("df_hr_all", inherits = TRUE)) {
    get("df_hr_all", inherits = TRUE)
  } else {
    read.csv(file.path(wd$out06, "HR_values_mCRPC_specific_multi_pre-allProteins.csv")) %>%
      filter(!is.na(Assay), Assay != "") %>%
      distinct(OlinkID, Assay, .keep_all = TRUE)
  }
  
  # 1) Binary status columns
  olink_detected <- proteins %in% unique(assay_map$Assay)
  public_detected <- proteins %in% unique(matched_proteins$Assay)
  # Case-insensitive matching for annotation sources
  is_cell_surface <- toupper(proteins) %in% toupper(unique(cs_prot))
  is_druggable <- toupper(proteins) %in% toupper(unique(drug_annotated_proteins))
  is_prognostic <- proteins %in% unique(dep_unique$Assay[dep_unique$Multi_Sig == "Yes"])
  
  # 2) Olink log2FC (mCRPC - Local) using pre-treatment only
  olink_fc <- df_npx_long %>%
    filter(cohort %in% c("A", "C"), grepl("^Pre", txt_stat_ordered)) %>%
    left_join(assay_map, by = "OlinkID") %>%
    filter(!is.na(Assay), Assay %in% proteins) %>%
    group_by(Assay, cohort) %>%
    summarise(mu = mean(NPX, na.rm = TRUE), .groups = "drop") %>%
    mutate(cohort_lbl = recode(cohort, A = "Local", C = "mCRPC")) %>%
    select(-cohort) %>%
    pivot_wider(names_from = cohort_lbl, values_from = mu) %>%
    mutate(Olink_log2FC = mCRPC - Local) %>%
    select(Assay, Olink_log2FC) %>%
    group_by(Assay) %>%
    summarise(Olink_log2FC = mean(Olink_log2FC, na.rm = TRUE), .groups = "drop")
  
  # 3) Public log2FC (mCRPC - Local_PC), following script 09 style
  public_fc <- as.data.frame(public_norm) %>%
    tibble::rownames_to_column("Protein") %>%
    pivot_longer(-Protein, names_to = "Biological_ID", values_to = "Value") %>%
    left_join(public_avg %>% select(Biological_ID, Group) %>% distinct(), by = "Biological_ID") %>%
    filter(Group %in% c("Local_PC", "mCRPC")) %>%
    group_by(Protein, Group) %>%
    summarise(mu = mean(Value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = mu) %>%
    mutate(Public_log2FC = mCRPC - Local_PC) %>%
    select(Protein, Public_log2FC) %>%
    left_join(matched_proteins %>% select(Protein, Assay) %>% distinct(), by = "Protein") %>%
    filter(!is.na(Assay), Assay %in% proteins) %>%
    group_by(Assay) %>%
    summarise(Public_log2FC = mean(Public_log2FC, na.rm = TRUE), .groups = "drop")
  
  # 4) Protein-level HR from section 03 all-protein multivariate output
  hr_protein <- tibble(Assay = proteins) %>%
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
    ) %>%
    select(Assay, HR, L95, U95)
  
  # Keep CPPS first, order protein rows by HR (highest -> lowest), NAs at bottom
  protein_order <- tibble(Assay = proteins, input_idx = seq_along(proteins)) %>%
    left_join(hr_protein %>% select(Assay, HR), by = "Assay") %>%
    arrange(is.na(HR), desc(HR), input_idx) %>%
    pull(Assay)
  
  # 5) CPPS HR row (only for HR panel)
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
  
  # Build row index including CPPS (top row)
  rows_all <- c("CPPS", protein_order)
  n_rows <- length(rows_all)
  
  # Binary matrices (protein rows always 0/1; CPPS row stays NA/blank)
  mat_detect <- matrix(NA_real_, nrow = n_rows, ncol = 2,
                       dimnames = list(rows_all, c("Olink", "Public")))
  prot_idx <- match(proteins, rows_all)
  mat_detect[prot_idx, "Olink"] <- 0
  mat_detect[prot_idx, "Public"] <- 0
  mat_detect[prot_idx, "Olink"] <- as.numeric(olink_detected)
  mat_detect[prot_idx, "Public"] <- as.numeric(public_detected)
  
  mat_cs <- matrix(NA_real_, nrow = n_rows, ncol = 1, dimnames = list(rows_all, "CS"))
  mat_dx <- matrix(NA_real_, nrow = n_rows, ncol = 1, dimnames = list(rows_all, "Dx"))
  mat_prog <- matrix(NA_real_, nrow = n_rows, ncol = 1, dimnames = list(rows_all, "Prog"))
  mat_cs[prot_idx, 1] <- 0
  mat_dx[prot_idx, 1] <- 0
  mat_prog[prot_idx, 1] <- 0
  mat_cs[prot_idx, 1] <- as.numeric(is_cell_surface)
  mat_dx[prot_idx, 1] <- as.numeric(is_druggable)
  mat_prog[prot_idx, 1] <- as.numeric(is_prognostic)
  
  # FC matrix (CPPS row blank)
  mat_fc <- matrix(NA_real_, nrow = n_rows, ncol = 2,
                   dimnames = list(rows_all, c("Olink FC", "Public FC")))
  if (nrow(olink_fc) > 0) {
    mat_fc[olink_fc$Assay, "Olink FC"] <- olink_fc$Olink_log2FC
  }
  if (nrow(public_fc) > 0) {
    mat_fc[public_fc$Assay, "Public FC"] <- public_fc$Public_log2FC
  }
  
  # Row-wise z-score distributions for Local vs mCRPC in each platform
  olink_z_long <- df_npx_long %>%
    filter(cohort %in% c("A", "C"), grepl("^Pre", txt_stat_ordered)) %>%
    left_join(assay_map, by = "OlinkID") %>%
    filter(!is.na(Assay), Assay %in% proteins) %>%
    mutate(Group = recode(cohort, A = "Local", C = "mCRPC")) %>%
    group_by(Assay) %>%
    mutate(z = as.numeric(scale(NPX))) %>%
    ungroup()
  
  public_z_long <- as.data.frame(public_norm) %>%
    tibble::rownames_to_column("Protein") %>%
    pivot_longer(-Protein, names_to = "Biological_ID", values_to = "Value") %>%
    left_join(public_avg %>% select(Biological_ID, Group) %>% distinct(), by = "Biological_ID") %>%
    filter(Group %in% c("Local_PC", "mCRPC")) %>%
    left_join(matched_proteins %>% select(Protein, Assay) %>% distinct(), by = "Protein") %>%
    filter(!is.na(Assay), Assay %in% proteins) %>%
    group_by(Assay) %>%
    mutate(z = as.numeric(scale(Value))) %>%
    ungroup()
  
  olink_local_list <- lapply(rows_all, function(rw) {
    if (rw == "CPPS") return(numeric(0))
    olink_z_long %>% filter(Assay == rw, Group == "Local") %>% pull(z) %>% .[is.finite(.)]
  })
  names(olink_local_list) <- rows_all
  olink_mcrpc_list <- lapply(rows_all, function(rw) {
    if (rw == "CPPS") return(numeric(0))
    olink_z_long %>% filter(Assay == rw, Group == "mCRPC") %>% pull(z) %>% .[is.finite(.)]
  })
  names(olink_mcrpc_list) <- rows_all
  
  public_local_list <- lapply(rows_all, function(rw) {
    if (rw == "CPPS") return(numeric(0))
    public_z_long %>% filter(Assay == rw, Group == "Local_PC") %>% pull(z) %>% .[is.finite(.)]
  })
  names(public_local_list) <- rows_all
  public_mcrpc_list <- lapply(rows_all, function(rw) {
    if (rw == "CPPS") return(numeric(0))
    public_z_long %>% filter(Assay == rw, Group == "mCRPC") %>% pull(z) %>% .[is.finite(.)]
  })
  names(public_mcrpc_list) <- rows_all
  
  # Mean z-score vectors (one value per row, CPPS stays NA)
  olink_mu_local <- sapply(rows_all, function(rw) {
    if (rw == "CPPS") return(NA_real_)
    vv <- olink_local_list[[rw]]
    if (length(vv) == 0) return(NA_real_)
    mean(vv, na.rm = TRUE)
  })
  olink_mu_mcrpc <- sapply(rows_all, function(rw) {
    if (rw == "CPPS") return(NA_real_)
    vv <- olink_mcrpc_list[[rw]]
    if (length(vv) == 0) return(NA_real_)
    mean(vv, na.rm = TRUE)
  })
  public_mu_local <- sapply(rows_all, function(rw) {
    if (rw == "CPPS") return(NA_real_)
    vv <- public_local_list[[rw]]
    if (length(vv) == 0) return(NA_real_)
    mean(vv, na.rm = TRUE)
  })
  public_mu_mcrpc <- sapply(rows_all, function(rw) {
    if (rw == "CPPS") return(NA_real_)
    vv <- public_mcrpc_list[[rw]]
    if (length(vv) == 0) return(NA_real_)
    mean(vv, na.rm = TRUE)
  })
  
  z_cap <- suppressWarnings(max(abs(c(olink_mu_local, olink_mu_mcrpc, public_mu_local, public_mu_mcrpc)), na.rm = TRUE))
  if (!is.finite(z_cap) || z_cap <= 0) z_cap <- 1
  z_lim <- c(-z_cap, z_cap)
  
  # HR vectors (include CPPS as final row)
  hr_df_all <- tibble(
    RowID = rows_all,
    HR = NA_real_,
    L95 = NA_real_,
    U95 = NA_real_
  )
  hr_df_all <- hr_df_all %>%
    left_join(hr_protein %>% rename(RowID = Assay), by = "RowID", suffix = c("", ".new")) %>%
    mutate(
      HR = coalesce(HR.new, HR),
      L95 = coalesce(L95.new, L95),
      U95 = coalesce(U95.new, U95)
    ) %>%
    select(RowID, HR, L95, U95)
  hr_df_all[hr_df_all$RowID == "CPPS", c("HR", "L95", "U95")] <- list(cpps_hr, cpps_l95, cpps_u95)
  
  # Color scales
  col_bin <- c("0" = "grey85", "1" = "black")
  fc_cap <- suppressWarnings(max(abs(mat_fc), na.rm = TRUE))
  if (!is.finite(fc_cap) || fc_cap <= 0) fc_cap <- 1
  col_fc <- circlize::colorRamp2(c(-fc_cap, 0, fc_cap), c("#2166AC", "white", "#B2182B"))
  
  # Main heatmap blocks
  ht_detect <- Heatmap(
    mat_detect,
    name = "Detected",
    col = col_bin,
    na_col = "white",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    row_names_gp = gpar(fontsize = 10),
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 8),
    width = unit(c(7, 7), "mm"),
    rect_gp = gpar(col = "grey70"),
    show_heatmap_legend = FALSE
  )
  
  ht_cs <- Heatmap(
    mat_cs,
    name = "CS",
    col = col_bin,
    na_col = "white",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    show_row_names = FALSE,
    column_names_rot = 0,
    column_names_gp = gpar(fontsize = 8),
    width = unit(10, "mm"),
    rect_gp = gpar(col = "grey50"),
    show_heatmap_legend = FALSE
  )
  
  ht_dx <- Heatmap(
    mat_dx,
    name = "Dx",
    col = col_bin,
    na_col = "white",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    show_row_names = FALSE,
    column_names_rot = 0,
    column_names_gp = gpar(fontsize = 8),
    width = unit(10, "mm"),
    rect_gp = gpar(col = "grey50"),
    show_heatmap_legend = FALSE
  )
  
  ht_prog <- Heatmap(
    mat_prog,
    name = "Prog",
    col = col_bin,
    na_col = "white",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    row_names_side = "left",
    show_row_names = FALSE,
    column_names_rot = 0,
    column_names_gp = gpar(fontsize = 8),
    width = unit(12, "mm"),
    rect_gp = gpar(col = "grey70"),
    show_heatmap_legend = FALSE
  )
  
  # Explicit spacer to create a visible gap between column 3 and 4 blocks
  ht_spacer <- Heatmap(
    matrix(NA_real_, nrow = n_rows, ncol = 1, dimnames = list(rows_all, "")),
    col = c("0" = "white"),
    na_col = "white",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    width = unit(3.5, "mm"),
    rect_gp = gpar(col = NA),
    show_heatmap_legend = FALSE
  )
  
  # Thin separator to keep flag columns visually distinct
  ht_sep <- Heatmap(
    matrix(NA_real_, nrow = n_rows, ncol = 1, dimnames = list(rows_all, "")),
    col = c("0" = "white"),
    na_col = "white",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    show_column_names = FALSE,
    width = unit(0.8, "mm"),
    rect_gp = gpar(col = NA),
    show_heatmap_legend = FALSE
  )
  
  ht_fc_olink <- Heatmap(
    mat_fc[, "Olink FC", drop = FALSE],
    name = "Olink FC",
    col = col_fc,
    na_col = "white",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 8),
    width = unit(14, "mm"),
    rect_gp = gpar(col = "grey80"),
    show_heatmap_legend = FALSE
  )
  
  ht_fc_public <- Heatmap(
    mat_fc[, "Public FC", drop = FALSE],
    name = "Public FC",
    col = col_fc,
    na_col = "white",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = FALSE,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 8),
    width = unit(14, "mm"),
    rect_gp = gpar(col = "grey80")
  )
  
  # Separate mean-bar columns (no overlap with FC)
  olink_local_anno <- rowAnnotation(
    `Olink Local` = anno_barplot(
      olink_mu_local,
      baseline = 0,
      which = "row",
      ylim = z_lim,
      axis = FALSE,
      border = FALSE,
      bar_width = 0.72,
      gp = gpar(fill = pal.olink["Local"], col = NA),
      width = unit(8, "mm")
    ),
    annotation_name_side = "top",
    annotation_name_rot = 45
  )
  olink_mcrpc_anno <- rowAnnotation(
    `Olink mCRPC` = anno_barplot(
      olink_mu_mcrpc,
      baseline = 0,
      which = "row",
      ylim = z_lim,
      axis = FALSE,
      border = FALSE,
      bar_width = 0.72,
      gp = gpar(fill = pal.olink["mCRPC"], col = NA),
      width = unit(8, "mm")
    ),
    annotation_name_side = "top",
    annotation_name_rot = 45
  )
  public_local_anno <- rowAnnotation(
    `Public Local` = anno_barplot(
      public_mu_local,
      baseline = 0,
      which = "row",
      ylim = z_lim,
      axis = FALSE,
      border = FALSE,
      bar_width = 0.72,
      gp = gpar(fill = pal.public["Local_PC"], col = NA),
      width = unit(8, "mm")
    ),
    annotation_name_side = "top",
    annotation_name_rot = 45
  )
  public_mcrpc_anno <- rowAnnotation(
    `Public mCRPC` = anno_barplot(
      public_mu_mcrpc,
      baseline = 0,
      which = "row",
      ylim = z_lim,
      axis = FALSE,
      border = FALSE,
      bar_width = 0.72,
      gp = gpar(fill = pal.public["mCRPC"], col = NA),
      width = unit(8, "mm")
    ),
    annotation_name_side = "top",
    annotation_name_rot = 45
  )
  
  # Row-wise HR forest panel (log2 scale so symmetric around HR=1)
  hr_x <- log2(hr_df_all$HR)
  hr_l <- log2(hr_df_all$L95)
  hr_u <- log2(hr_df_all$U95)
  
  hr_anno <- rowAnnotation(
    `HR (OS)` = anno_empty(
      border = TRUE,
      width = unit(45, "mm"),
      zoom = TRUE
    ),
    annotation_name_side = "top",
    annotation_name_rot = 0
  )
  
  ht_all <- ht_detect + ht_spacer + ht_cs + ht_sep + ht_dx + ht_sep + ht_prog + ht_sep +
    olink_local_anno + olink_mcrpc_anno + ht_sep +
    ht_fc_olink + ht_fc_public + ht_sep +
    public_local_anno + public_mcrpc_anno + hr_anno
  
  draw(ht_all,
       heatmap_legend_side = "right",
       annotation_legend_side = "right",
       column_title = "Protein integration summary")
  
  # Draw CI and points in HR panel viewport after layout
  decorate_annotation("HR (OS)", {
    y_pos <- rev(seq_len(n_rows)) - 0.5
    ok <- which(is.finite(hr_x) & is.finite(hr_l) & is.finite(hr_u))
    if (length(ok) > 0) {
      x_min <- min(hr_l[ok], log2(0.5), na.rm = TRUE)
      x_max <- max(hr_u[ok], log2(4), na.rm = TRUE)
    } else {
      x_min <- log2(0.5)
      x_max <- log2(4)
    }
    pad <- 0.08 * (x_max - x_min + 1e-8)
    x_lims <- c(x_min - pad, x_max + pad)
    
    # x-axis ticks at clinically familiar HR values
    tick_hr <- c(0.5, 1, 2, 4)
    tick_x <- log2(tick_hr)
    
    # baseline at HR=1
    grid.lines(
      x = unit(rep((log2(1) - x_lims[1]) / diff(x_lims), 2), "npc"),
      y = unit(c(0, 1), "npc"),
      gp = gpar(col = "grey55", lty = 2)
    )
    
    # CIs + points
    for (i in ok) {
      x0 <- (hr_l[i] - x_lims[1]) / diff(x_lims)
      x1 <- (hr_u[i] - x_lims[1]) / diff(x_lims)
      xp <- (hr_x[i] - x_lims[1]) / diff(x_lims)
      yp <- (y_pos[i]) / n_rows
      grid.lines(x = unit(c(x0, x1), "npc"), y = unit(c(yp, yp), "npc"), gp = gpar(col = "black", lwd = 1))
      grid.points(x = unit(xp, "npc"), y = unit(yp, "npc"), pch = 16, size = unit(2.2, "mm"))
    }
    
    # Tick marks and labels
    for (j in seq_along(tick_x)) {
      if (tick_x[j] < x_lims[1] || tick_x[j] > x_lims[2]) next
      tx <- (tick_x[j] - x_lims[1]) / diff(x_lims)
      grid.lines(x = unit(c(tx, tx), "npc"), y = unit(c(0, 0.02), "npc"), gp = gpar(col = "black"))
      grid.text(as.character(tick_hr[j]), x = unit(tx, "npc"), y = unit(-0.02, "npc"), gp = gpar(fontsize = 6))
    }
    
    # Axis title at bottom of HR panel
    grid.text("Hazard Ratio", x = unit(0.5, "npc"), y = unit(-0.08, "npc"), gp = gpar(fontsize = 7))
  })
}
