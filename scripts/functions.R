theme_custom <- 
  theme_bw() +
  theme(panel.grid = element_blank())


# ------- Colours --------- #

alpha = 1*255

red <- grDevices::col2rgb('#FE1F04')
orange <- grDevices::col2rgb('#FF8C22')
yellow <- grDevices::col2rgb('#FFC700')
green <- grDevices::col2rgb('#27AE55')
teal <- grDevices::col2rgb('#077183')
turqoise <- grDevices::col2rgb('#00C7E1')
lightblue <- grDevices::col2rgb('#A2D9F5')
darkblue <- grDevices::col2rgb('#00559E')
purple <- grDevices::col2rgb('#6A27AE')
pink <- grDevices::col2rgb('#FF51B8')

red<- grDevices::rgb(red[1], red[2], red[3], alpha, maxColorValue = 255)
orange <- grDevices::rgb(orange[1], orange[2], orange[3], alpha, maxColorValue = 255)
yellow <- grDevices::rgb(yellow[1], yellow[2], yellow[3], alpha, maxColorValue = 255)
green <- grDevices::rgb(green[1], green[2], green[3], alpha, maxColorValue = 255)
teal <- grDevices::rgb(teal[1], teal[2], teal[3], alpha, maxColorValue = 255)
turqoise <- grDevices::rgb(turqoise[1], turqoise[2], turqoise[3], alpha, maxColorValue = 255)
lightblue <-  grDevices::rgb(lightblue[1], lightblue[2], lightblue[3], alpha, maxColorValue = 255)
darkblue <- grDevices::rgb(darkblue[1], darkblue[2], darkblue[3], alpha, maxColorValue = 255)
purple <- grDevices::rgb(purple[1], purple[2], purple[3], alpha, maxColorValue = 255)
pink <- grDevices::rgb(pink[1], pink[2], pink[3], alpha, maxColorValue = 255)


# https://r-charts.com/colors/
pal.cohort <- paletteer::paletteer_d("ggthemes::Tableau_10",10)
pal.cohort <- (pal.cohort[c(2,1,3)])
#pal.cohort <- rev(c("#003f5c", "#bc5090", "#ffa600"))
#pal.cohort <- rev(c("#003f5c", "#bc5090", "#ffa600"))


df.pal.study <- data.frame(
  study = c("Psomagen", "Q-1335", "Q-1580"),
  pal = c("#006400", "#E9967A", "#9400D3")
)
pal.study <- df.pal.study$pal

pal.cohort <- paletteer::paletteer_d("ggthemes::Tableau_10",10)
pal.cohort <- (pal.cohort[c(2,1,3)])


pal.cohort.n <- c("#AEC7E8", "#FFBB78", "#98DF8A")
pal.cohort.n2 <- c("#1f77b4", "#ff7f0e", "#2ca02c")
pal.cond3 <- paletteer::paletteer_d("futurevisions::atomic_orange", 3)
pal.txt_stat <- c("lightblue", "orange")



# ------- Functions imported from olink package ---# 
#This function is called by various other functions to perform checks of NPX-data. For now, it only looks for assays which have all NPX=NA, but there are other redundant tasks that could be moved here
npxCheck <- function(df){
  # # Check whether df contains NPX or QUANT
  if ('NPX' %in% colnames(df)) {
    data_type <- 'NPX'
  } else if ('Quantified_value' %in% colnames(df)) {
    data_type <- 'Quantified_value'
  } else {
    stop('Neither NPX or Quantified_value column present in the data')}
  
  #### Identify assays that have only NA:s ####
  all_nas <- df  %>%
    dplyr::group_by(OlinkID) %>%
    dplyr::summarise(n = dplyr::n(), n_na = sum(is.na(!!rlang::ensym(data_type)))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n == n_na) %>%
    dplyr::pull(OlinkID)
  
  if(length(all_nas) > 0) {
    
    warning(paste0('The assays ',
                   paste(all_nas, collapse = ', '),
                   ' have NPX=NA for all samples. They will be excluded from the analysis'),
            call. = FALSE)
    
  }
  
  return(list(all_nas = all_nas, data_type = data_type))
}

# ----------- Custom t.test ------------- #
olink_ttest_2 <-
  function (df, variable, pair_id, ...) 
{
  if (missing(df) | missing(variable)) {
    stop("The df and variable arguments need to be specified.")
  }
  df <- df %>% dplyr::filter(stringr::str_detect(OlinkID, "OID[0-9]{5}"))
  removed.sampleids <- NULL
  removed.sampleids <- unique(c(removed.sampleids, df$SampleID[is.na(df[[variable]])]))
  df <- df[!is.na(df[[variable]]), ]
  if (!is.null(removed.sampleids) & length(removed.sampleids) > 
      0) {
    message("Samples removed due to missing variable levels: ", 
            paste(removed.sampleids, collapse = ", "))
  }
  if (!missing(pair_id)) {
    missing.pair <- NULL
    missing.pair <- df$SampleID[is.na(df[[pair_id]])]
    if (!is.null(missing.pair) & length(missing.pair) > 0) {
      message("Samples removed due to missing pair ID: ", 
              paste(missing.pair, collapse = ", "))
    }
    df <- df[!is.na(df[[pair_id]]), ]
    removed.sampleids <- unique(c(removed.sampleids, missing.pair))
  }
  if (is.character(df[[variable]])) {
    df[[variable]] <- factor(df[[variable]])
    message(paste0("Variable converted from character to factor: ", 
                   variable))
  }
  else if (!is.factor(df[[variable]])) {
    stop(paste0("The grouping variable ", variable, "is neither factor nor character. Only character and factor variable types allowed."))
  }
  var_levels <- levels(df[[variable]])
  number_of_levels <- length(var_levels)
  if (!(number_of_levels == 2)) {
    stop(paste0("The number of levels in the factor needs to be 2. Your factor has ", 
                number_of_levels, " levels."))
  }
  number_of_samples_w_more_than_one_level <- df %>% dplyr::group_by(SampleID, 
                                                                    Index) %>% dplyr::summarise(n_levels = dplyr::n_distinct(!!rlang::ensym(variable), 
                                                                                                                             na.rm = TRUE)) %>% dplyr::ungroup() %>% dplyr::filter(n_levels > 
                                                                                                                                                                                     1) %>% nrow(.)
  if (number_of_samples_w_more_than_one_level > 0) {
    stop(paste0("There are ", number_of_samples_w_more_than_one_level, 
                " samples that do not have a unique level for your variable. Only one level per sample is allowed."))
  }
  npxCheck <- npxCheck(df)
  nas_in_level <- df %>% dplyr::filter(!(OlinkID %in% npxCheck$all_nas)) %>% 
    dplyr::group_by(OlinkID, !!rlang::ensym(variable)) %>% 
    dplyr::summarise(n = dplyr::n(), n_na = sum(is.na(NPX))) %>% 
    dplyr::ungroup() %>% dplyr::filter(n - n_na <= 1) %>% 
    dplyr::pull(OlinkID)
  if (length(nas_in_level) > 0) {
    warning(paste0("The assays ", paste(nas_in_level, collapse = ", "), 
                   " have too few datapoints in one level of the factor. They will not be tested."), 
            call. = FALSE)
  }
  df <- df %>% dplyr::filter(!(OlinkID %in% npxCheck$all_nas)) %>% 
    dplyr::filter(!(OlinkID %in% nas_in_level))
  if (nrow(df) == 0) {
    stop("No assays passing initial check. T-test will not be performed.")
  }
  if (!missing(pair_id)) {
    if (!pair_id %in% colnames(df)) 
      stop(paste0("Column ", pair_id, " not found."))
    if (!is_tibble(df)) {
      message("Converting data frame to tibble.")
      df <- dplyr::as_tibble(df)
    }
    ct_pairs <- df %>% dplyr::filter(!(OlinkID %in% npxCheck$all_nas)) %>% 
      dplyr::filter(!(OlinkID %in% nas_in_level)) %>% dplyr::filter(!is.na(!!rlang::ensym(variable))) %>% 
      dplyr::group_by(OlinkID, !!rlang::ensym(pair_id)) %>% 
      dplyr::summarize(n = dplyr::n())
    if (!all(ct_pairs$n <= 2)) 
      stop(paste0("Each pair identifier must identify no more than 2 unique samples. Check pairs: ", 
                  paste(unique(ct_pairs[[pair_id]][ct_pairs$n > 
                                                     2]), collapse = ", ")))
    message(paste0("Paired t-test is performed on ", var_levels[1], 
                   " - ", var_levels[2], "."))
    p.val <- df %>% dplyr::select(dplyr::all_of(c("OlinkID", 
                                                  "UniProt", "Assay", "Panel", "NPX", variable, pair_id))) %>% 
      tidyr::pivot_wider(names_from = dplyr::all_of(variable), 
                         values_from = "NPX") %>% dplyr::group_by(Assay, 
                                                                  OlinkID, UniProt, Panel) %>% dplyr::do(broom::tidy(t.test(x = .[[var_levels[1]]], 
                                                                                                                            y = .[[var_levels[2]]], paired = TRUE, ...))) %>% 
      dplyr::ungroup() %>% dplyr::mutate(Adjusted_pval = p.adjust(p.value, 
                                                                  method = "fdr")) %>% dplyr::mutate(Threshold = ifelse(Adjusted_pval < 
                                                                                                                          0.05, "Significant", "Non-significant")) %>% dplyr::arrange(p.value)
  }
  else {
    message(paste0("T-test is performed on ", var_levels[1], 
                   " - ", var_levels[2], "."))
    p.val <- df %>% dplyr::group_by(Assay, OlinkID, UniProt, 
                                    Panel) %>% dplyr::do(broom::tidy(t.test(NPX ~ !!rlang::ensym(variable),
                                                                            alternative = "greater",
                                                                            #alternative = "less",
                                                                            data = ., ...))) %>% dplyr::ungroup() %>% dplyr::mutate(Adjusted_pval = p.adjust(p.value, 
                                                                                                                                                             method = "fdr")) %>% dplyr::mutate(Threshold = ifelse(Adjusted_pval < 
                                                                                                                                                                                                                     0.05, "Significant", "Non-significant")) %>% dplyr::rename(`:=`(!!var_levels[1], 
                                                                                                                                                                                                                                                                                     estimate1)) %>% dplyr::rename(`:=`(!!var_levels[2], 
                                                                                                                                                                                                                                                                                                                        estimate2)) %>% dplyr::arrange(p.value)
  }
  return(p.val)
}


# ------- t.test with fold change  --------- #
t.test_cat <- function(dat, cat){
  
  print(paste0("t-test on ", cat))
  
  dat1 <- dat # %>%
  #dplyr::filter(!!as.name(cat) %in% c("Severe", "None_Mild") )
  
  dat2 <- 
    olink_ttest_2(df = dat1,
                variable = cat)
  
  
  # ---------- #
  # Add fold change 
  # ---------- #
  
  # Get the average expression per group
  dat1_g <- dat1 %>%
    dplyr::group_by(OlinkID, !!as.name(cat)) %>%
    dplyr::summarise(NPX_mean = mean(NPX))
  
  # Convert to wide format
  dat1_w <- dat1_g %>%
    tidyr::pivot_wider(names_from = !!as.name(cat), values_from = NPX_mean) %>%
    dplyr::rename(Mean_N = N,
                  Mean_Y= Y) %>%
    dplyr::mutate(log2FC = Mean_Y - Mean_N)
  
  
  
  dat_tmp <- dat2 %>%
    dplyr::left_join(dat1_w)
  
  
  # Simplify output
  dat3 <- dat_tmp %>%
    dplyr::select(Assay, OlinkID, UniProt, Panel, p.value, Adjusted_pval, Threshold, Mean_Y, Mean_N, log2FC)
  
  write.csv(dat3, 
            file.path(wd$outCurr, 
                      paste0("Rest_vs_", cat, ".csv")))
  
  return(dat3)
  
}


# ------- t.test with fold change and add column identifier  --------- #

t.test_cat_col <- function(dat, cat){
  
  print(paste0("t-test on ", cat))
  
  dat1 <- dat # %>%
  #dplyr::filter(!!as.name(cat) %in% c("Severe", "None_Mild") )
  
  dat2 <- 
    olink_ttest_2(df = dat1,
                variable = cat)
  
  
  # ---------- #
  # Add fold change 
  # ---------- #
  
  # Get the average expression per group
  dat1_g <- dat1 %>%
    dplyr::group_by(OlinkID, !!as.name(cat)) %>%
    dplyr::summarise(NPX_mean = mean(NPX))
  
  # Convert to wide format
  dat1_w <- dat1_g %>%
    tidyr::pivot_wider(names_from = !!as.name(cat), values_from = NPX_mean) %>%
    dplyr::rename(Mean_N = N,
                  Mean_Y= Y) %>%
    dplyr::mutate(log2FC = Mean_Y - Mean_N)
  
  
  
  dat_tmp <- dat2 %>%
    dplyr::left_join(dat1_w)
  
  
  # Simplify output
  dat3 <- dat_tmp %>%
    dplyr::select(Assay, OlinkID, UniProt, Panel, p.value, Adjusted_pval, Threshold, Mean_Y, Mean_N, log2FC) %>%
    dplyr::mutate(Cat = cat)
  
  write.csv(dat3, 
            file.path(wd$outCurr, 
                      paste0("Rest_vs_", cat, ".csv")))
  
  return(dat3)
  
}

# ------- t.test with fold change and add column identifier  --------- #
# REMOVE NA ASSAY AND REMOVE NA IN NPX

t.test_cat_col2 <- function(dat, cat){
  
  print(paste0("t-test on ", cat))
  
  dat1 <- dat  %>%
    dplyr::filter(!is.na(Assay)) 
  
  dat2 <- 
    olink_ttest_2(df = dat1,
                  variable = cat)
  
  
  # ---------- #
  # Add fold change 
  # ---------- #
  
  # Get the average expression per group
  dat1_g <- dat1 %>%
    dplyr::filter(!is.na(NPX)) %>% 
    dplyr::group_by(OlinkID, !!as.name(cat)) %>%
    dplyr::summarise(NPX_mean = mean(NPX))
  
  # Convert to wide format
  dat1_w <- dat1_g %>%
    tidyr::pivot_wider(names_from = !!as.name(cat), values_from = NPX_mean) %>%
    dplyr::rename(Mean_N = N,
                  Mean_Y= Y) %>%
    dplyr::mutate(log2FC = Mean_N - Mean_Y,
                  baseFC = 2^log2FC)
  
  
  
  dat_tmp <- dat2 %>%
    dplyr::left_join(dat1_w)
  
  
  # Simplify output
  dat3 <- dat_tmp %>%
    dplyr::select(Assay, OlinkID, UniProt, Panel, p.value, Adjusted_pval, Threshold, Mean_Y, Mean_N, log2FC, baseFC) %>%
    dplyr::mutate(Cat = cat)
  
  write.csv(dat3, 
            file.path(wd$outCurr, 
                      paste0("Rest_vs_", cat, ".csv")))
  
  return(dat3)
  
}

# Function to extract cox regression from clinical vars
cox_extract <- function(clin_var, meta_in){
  # Initialize an empty list to store results
  all_results <- list()
  
  # Loop over each protein in cd_prot
  for (eg_clin in clin_var) {
  
  # Rename var
  tmp_df <- meta_in %>% 
    dplyr::rename(var = !!as.name(eg_clin) )
    
  # Fit the Cox proportional hazards model
  model <- coxph(Surv(OS, death) ~ var, data = tmp_df)
    
  # Extract coefficients and confidence intervals
  coef_summary <- summary(model)$coefficients
  coef_summary2 <- summary(model)$conf.int
  coef_names <- row.names(coef_summary)
  coef <- coef_summary[, 1]
  lower_bound <- coef_summary2[, 3]
  upper_bound <- coef_summary2[, 4]
  p_value <- coef_summary[, 5]
  exp_coef <- coef_summary2[, 1]
  
  # Create a data frame for the current protein
  result_df <- data.frame(
    OlinkID = eg_clin,
    Coefficient = coef,
    Lower_Bound = lower_bound,
    Upper_Bound = upper_bound,
    P_Value = p_value,
    HR = exp_coef
  )
  
  row.names(result_df) <- coef_names
  
  # Append to the list of all results
  all_results[[eg_clin]] <- result_df
  }
  return(all_results)
}

# Function to extract cox regression value from proteins
analyze_proteins_bin <- function(cd_prot, npx_in, meta_in) {
  # Initialize an empty list to store results
  all_results <- list()
  
  # Loop over each protein in cd_prot
  for (eg_protein in cd_prot) {
    #print(eg_protein)
    # Extract NPX value for the current protein
    npx_exp <- npx_in %>%
      filter(OlinkID %in% eg_protein) %>%
      select(SampleID, NPX, Assay) %>%
      dplry::rename(sample_id = SampleID) %>%
      mutate(
        sample_id = as.character(sample_id),
        npx_median = median(NPX),
        npx_bin = case_when(NPX >= npx_median ~ "high", NPX < npx_median ~ "low"),
        npx_bin = factor(npx_bin, levels = c("low", "high"))
      )
    
    # Combine with metadata
    tmp_df <- meta_in %>% left_join(npx_exp)
    
    # Fit the Cox proportional hazards model
    model <- coxph(Surv(OS, death) ~ NPX, data = tmp_df)
    
    # Extract coefficients and confidence intervals
    coef_summary <- summary(model)$coefficients
    coef_summary2 <- summary(model)$conf.int
    coef_names <- row.names(coef_summary)
    coef <- coef_summary[, 1]
    lower_bound <- coef_summary2[, 3]
    upper_bound <- coef_summary2[, 4]
    p_value <- coef_summary[, 5]
    exp_coef <- coef_summary2[, 1]
    
    # Create a data frame for the current protein
    result_df <- data.frame(
      OlinkID = eg_protein,
      Coefficient = coef,
      Lower_Bound = lower_bound,
      Upper_Bound = upper_bound,
      P_Value = p_value,
      HR = exp_coef
    )
    
    row.names(result_df) <- coef_names
    
    # Append to the list of all results
    all_results[[eg_protein]] <- result_df
  }
  
  return(all_results)
}



# Function to extract cox regression value from proteins
analyze_proteins_bin <- function(cd_prot, npx_in, meta_in) {
  # Initialize an empty list to store results
  all_results <- list()
  
  # Loop over each protein in cd_prot
  for (eg_protein in cd_prot) {
    #print(eg_protein)
    # Extract NPX value for the current protein
    npx_exp <- npx_in %>%
      filter(OlinkID %in% eg_protein) %>%
      select(SampleID, NPX, Assay) %>%
      dplyr::rename(sample_id = SampleID) %>%
      mutate(
        sample_id = as.character(sample_id),
        npx_median = median(NPX),
        npx_bin = case_when(NPX >= npx_median ~ "high", NPX < npx_median ~ "low"),
        npx_bin = factor(npx_bin, levels = c("low", "high"))
      )
    
    # Combine with metadata
    tmp_df <- meta_in %>% left_join(npx_exp)
    
    # Fit the Cox proportional hazards model
    model <- coxph(Surv(OS, death) ~ NPX, data = tmp_df)
    
    # Extract coefficients and confidence intervals
    coef_summary <- summary(model)$coefficients
    coef_summary2 <- summary(model)$conf.int
    coef_names <- row.names(coef_summary)
    coef <- coef_summary[, 1]
    lower_bound <- coef_summary2[, 3]
    upper_bound <- coef_summary2[, 4]
    p_value <- coef_summary[, 5]
    exp_coef <- coef_summary2[, 1]
    
    # Create a data frame for the current protein
    result_df <- data.frame(
      OlinkID = eg_protein,
      Coefficient = coef,
      Lower_Bound = lower_bound,
      Upper_Bound = upper_bound,
      P_Value = p_value,
      HR = exp_coef
    )
    
    row.names(result_df) <- coef_names
    
    # Append to the list of all results
    all_results[[eg_protein]] <- result_df
  }
  
  return(all_results)
}



# Function to perfom and extract cox regression value from proteins
# We use the 'cluster' parameter to account for repeated non-independent measurement
analyze_proteins_cluster <- function(cd_prot, npx_in, meta_in) {
  # Initialize an empty list to store results
  all_results <- list()
  
  # Loop over each protein in cd_prot
  for (eg_protein in cd_prot) {
    #print(eg_protein)
    # Extract NPX value for the current protein
    npx_exp <- npx_in %>%
      filter(OlinkID %in% eg_protein) %>%
      select(SampleID, NPX, Assay) %>%
      rename(sample_id = SampleID) %>%
      mutate(
        sample_id = as.character(sample_id),
        npx_median = median(NPX),
        npx_bin = case_when(NPX >= npx_median ~ "high", NPX < npx_median ~ "low"),
        npx_bin = factor(npx_bin, levels = c("low", "high"))
      )
    
    # Combine with metadata
    tmp_df <- meta_in %>% left_join(npx_exp)
    
    # Fit the Cox proportional hazards model
    model <- coxph(Surv(OS, death) ~ NPX, data = tmp_df, cluster=mrn)
    
    # Extract coefficients and confidence intervals
    coef_summary <- summary(model)$coefficients
    coef_summary2 <- summary(model)$conf.int
    coef_names <- row.names(coef_summary)
    coef <- coef_summary[, 1]
    lower_bound <- coef_summary2[, 3]
    upper_bound <- coef_summary2[, 4]
    p_value <- coef_summary[, 6]
    exp_coef <- coef_summary2[, 1]
    
    # Create a data frame for the current protein
    result_df <- data.frame(
      OlinkID = eg_protein,
      Coefficient = coef,
      Lower_Bound = lower_bound,
      Upper_Bound = upper_bound,
      P_Value = p_value,
      HR = exp_coef
    )
    
    row.names(result_df) <- coef_names
    
    # Append to the list of all results
    all_results[[eg_protein]] <- result_df
  }
  
  return(all_results)
}



# # Function to extract cox regression value from proteins in a multi-variate
# analyze_proteins_multi <- function(cd_prot, npx_in, meta_in) {
#   # Initialize an empty list to store results
#   all_results <- list()
#   
#   # Loop over each protein in cd_prot
#   for (eg_protein in cd_prot) {
#     #print(eg_protein)
#     # Extract NPX value for the current protein
#     npx_exp <- npx_in %>%
#       filter(OlinkID %in% eg_protein) %>%
#       select(SampleID, NPX, Assay) %>%
#       rename(sample_id = SampleID) %>%
#       mutate(
#         sample_id = as.character(sample_id),
#         npx_median = median(NPX),
#         npx_bin = case_when(NPX >= npx_median ~ "high", NPX < npx_median ~ "low"),
#         npx_bin = factor(npx_bin, levels = c("low", "high"))
#       )
#     
#     # Combine with metadata
#     tmp_df <- meta_in %>% left_join(npx_exp)
#     
#     # Fit the Cox proportional hazards model
#     model <- coxph(Surv(OS, death) ~ NPX + psa + ldh + alk_ph, data = tmp_df)
#     
#     # Extract coefficients and confidence intervals
#     # coef_summary <- summary(model)$coefficients["npx_binhigh",]
#     # coef_summary2 <- summary(model)$conf.int["npx_binhigh",]
#     coef_summary <- summary(model)$coefficients["NPX",]
#     coef_summary2 <- summary(model)$conf.int["NPX",]
#     #coef_names <- row.names(coef_summary)
#     coef <- coef_summary[1]
#     lower_bound <- coef_summary2[3]
#     upper_bound <- coef_summary2[4]
#     p_value <- coef_summary[5]
#     exp_coef <- coef_summary2[1]
#     
#     # Create a data frame for the current protein
#     result_df <- data.frame(
#       OlinkID = eg_protein,
#       Coefficient = coef,
#       Lower_Bound = lower_bound,
#       Upper_Bound = upper_bound,
#       P_Value = p_value,
#       HR = exp_coef
#     )
#     
#     
#     # Append to the list of all results
#     all_results[[eg_protein]] <- result_df
#   }
#   
#   return(all_results)
# }



# Function to perfom and extract cox regression value from proteins
# We use the 'cluster' parameter to account for repeated non-independent measurement
# Use the median as the bin
analyze_proteins_multi_cluster_bin <- function(cd_prot, npx_in, meta_in) {
  # Initialize an empty list to store results
  all_results <- list()
  
  # Loop over each protein in cd_prot
  for (eg_protein in cd_prot) {
    #print(eg_protein)
    # Extract NPX value for the current protein
    npx_exp <- npx_in %>%
      filter(OlinkID %in% eg_protein) %>%
      select(SampleID, NPX, Assay) %>%
      dplyr::rename(sample_id = SampleID) %>%
      mutate(
        sample_id = as.character(sample_id),
        npx_median = median(NPX),
        npx_bin = case_when(NPX >= npx_median ~ "high", NPX < npx_median ~ "low"),
        npx_bin = factor(npx_bin, levels = c("low", "high"))
      )
    
    # Combine with metadata
    tmp_df <- meta_in %>% left_join(npx_exp)
    
    # Fit the Cox proportional hazards model
    model <- coxph(Surv(OS, death) ~ npx_bin + psa + ldh + alk_ph, data = tmp_df, cluster=mrn)
    
    # Extract coefficients and confidence intervals
    coef_summary <- summary(model)$coefficients["npx_binhigh",]
    coef_summary2 <- summary(model)$conf.int["npx_binhigh",]
    #coef_names <- row.names(coef_summary)
    coef <- coef_summary[1]
    lower_bound <- coef_summary2[3]
    upper_bound <- coef_summary2[4]
    p_value <- coef_summary[6]
    exp_coef <- coef_summary2[1]
    
    # Create a data frame for the current protein
    result_df <- data.frame(
      OlinkID = eg_protein,
      Coefficient = coef,
      Lower_Bound = lower_bound,
      Upper_Bound = upper_bound,
      P_Value = p_value,
      HR = exp_coef
    )
    
    
    # Append to the list of all results
    all_results[[eg_protein]] <- result_df
  }
  
  return(all_results)
}

# Function to perfom and extract cox regression value from proteins
# We use the 'cluster' parameter to account for repeated non-independent measurement
# Use the median as the bin
analyze_proteins_multi_cluster <- function(cd_prot, npx_in, meta_in) {
  # Initialize an empty list to store results
  all_results <- list()
  
  # Loop over each protein in cd_prot
  for (eg_protein in cd_prot) {
    #print(eg_protein)
    # Extract NPX value for the current protein
    npx_exp <- npx_in %>%
      filter(OlinkID %in% eg_protein) %>%
      select(SampleID, NPX, Assay) %>%
      rename(sample_id = SampleID) %>%
      mutate(
        sample_id = as.character(sample_id),
        npx_median = median(NPX),
        npx_bin = case_when(NPX >= npx_median ~ "high", NPX < npx_median ~ "low"),
        npx_bin = factor(npx_bin, levels = c("low", "high"))
      )
    
    # Combine with metadata
    tmp_df <- meta_in %>% left_join(npx_exp)
    
    # Fit the Cox proportional hazards model
    model <- coxph(Surv(OS, death) ~ NPX + psa + ldh + alk_ph, data = tmp_df, cluster=mrn)
    
    # Extract coefficients and confidence intervals
    coef_summary <- summary(model)$coefficients["NPX",]
    coef_summary2 <- summary(model)$conf.int["NPX",]
    #coef_names <- row.names(coef_summary)
    coef <- coef_summary[1]
    lower_bound <- coef_summary2[3]
    upper_bound <- coef_summary2[4]
    p_value <- coef_summary[6]
    exp_coef <- coef_summary2[1]
    
    # Create a data frame for the current protein
    result_df <- data.frame(
      OlinkID = eg_protein,
      Coefficient = coef,
      Lower_Bound = lower_bound,
      Upper_Bound = upper_bound,
      P_Value = p_value,
      HR = exp_coef
    )
    
    
    # Append to the list of all results
    all_results[[eg_protein]] <- result_df
  }
  
  return(all_results)
}




# ---- Custom wilcox test --- 
olink_wilcox_2 <- function (df, variable, pair_id, ...) 
  {
    if (missing(df) | missing(variable)) {
      stop("The df and variable arguments need to be specified.")
    }
    df <- df %>% filter(stringr::str_detect(OlinkID, "OID[0-9]{5}"))
    removed.sampleids <- NULL
    removed.sampleids <- unique(c(removed.sampleids, df$SampleID[is.na(df[[variable]])]))
    df <- df[!is.na(df[[variable]]), ]
    if (!is.null(removed.sampleids) & length(removed.sampleids) > 0) {
      message("Samples removed due to missing variable levels: ", 
              paste(removed.sampleids, collapse = ", "))
    }
    if (!missing(pair_id)) {
      missing.pair <- NULL
      missing.pair <- df$SampleID[is.na(df[[pair_id]])]
      if (!is.null(missing.pair) & length(missing.pair) > 0) {
        message("Samples removed due to missing pair ID: ", 
                paste(missing.pair, collapse = ", "))
      }
      df <- df[!is.na(df[[pair_id]]), ]
      removed.sampleids <- unique(c(removed.sampleids, missing.pair))
    }
    if (is.character(df[[variable]])) {
      df[[variable]] <- factor(df[[variable]])
      message(paste0("Variable converted from character to factor: ", 
                     variable))
    }
    else if (!is.factor(df[[variable]])) {
      stop(paste0("The grouping variable ", variable, "is neither factor nor character. Only character and factor variable types allowed."))
    }
    var_levels <- levels(df[[variable]])
    number_of_levels <- length(var_levels)
    if (!(number_of_levels == 2)) {
      stop(paste0("The number of levels in the factor needs to be 2. Your factor has ", 
                  number_of_levels, " levels."))
    }
    number_of_samples_w_more_than_one_level <- df %>% dplyr::group_by(SampleID) %>% 
      dplyr::summarise(n_levels = dplyr::n_distinct(!!rlang::ensym(variable), 
                                                    na.rm = T)) %>% dplyr::ungroup() %>% dplyr::filter(n_levels > 
                                                                                                         1) %>% nrow(.)
    if (number_of_samples_w_more_than_one_level > 0) {
      stop(paste0("There are ", number_of_samples_w_more_than_one_level, 
                  " samples that do not have a unique level for your variable. Only one level per sample is allowed."))
    }
    npxCheck <- npxCheck(df)
    data_type <- npxCheck$data_type
    nas_in_level <- df %>% dplyr::filter(!(OlinkID %in% npxCheck$all_nas)) %>% 
      dplyr::group_by(OlinkID, !!rlang::ensym(variable)) %>% 
      dplyr::summarise(n = dplyr::n(), n_na = sum(is.na(!!rlang::ensym(data_type)))) %>% 
      dplyr::ungroup() %>% dplyr::filter(n == n_na) %>% dplyr::pull(OlinkID)
    if (length(nas_in_level) > 0) {
      warning(paste0("The assays ", paste(nas_in_level, collapse = ", "), 
                     " have only NA:s in one level of the factor. They will not be tested."), 
              call. = F)
    }
    if (!missing(pair_id)) {
      if (!pair_id %in% colnames(df)) 
        stop(paste0("Column ", pair_id, " not found."))
      if (!tibble::is_tibble(df)) {
        message("Converting data frame to tibble.")
        df <- dplyr::as_tibble(df)
      }
      ct_pairs <- df %>% dplyr::filter(!(OlinkID %in% npxCheck$all_nas)) %>% 
        dplyr::filter(!(OlinkID %in% nas_in_level)) %>% dplyr::filter(!is.na(!!rlang::ensym(variable))) %>% 
        dplyr::group_by(OlinkID, !!rlang::ensym(pair_id)) %>% 
        dplyr::summarize(n = dplyr::n())
      if (!all(ct_pairs$n <= 2)) 
        stop(paste0("Each pair identifier must identify no more than 2 unique samples. Check pairs: ", 
                    paste(unique(ct_pairs[[pair_id]][ct_pairs$n > 2]), collapse = ", ")))
      message(paste0("Paired Mann-Whitney U Test is performed on ", 
                     var_levels[1], " - ", var_levels[2], "."))
      p.val <- df %>% dplyr::filter(!(OlinkID %in% npxCheck$all_nas)) %>% 
        dplyr::filter(!(OlinkID %in% nas_in_level)) %>% dplyr::select(all_of(c("OlinkID", 
                                                                               "UniProt", "Assay", "Panel", data_type, variable, pair_id))) %>% 
        tidyr::pivot_wider(names_from = all_of(variable), values_from = all_of(data_type)) %>% 
        dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>% 
        dplyr::do(tidy(stats::wilcox.test(x = .[[var_levels[1]]], y = .[[var_levels[2]]], 
                                          paired = TRUE, conf.int = TRUE, alternative = "greater", ...))) %>% 
        dplyr::ungroup() %>% dplyr::mutate(Adjusted_pval = p.adjust(p.value, method = "fdr")) %>% 
        dplyr::mutate(Threshold = ifelse(Adjusted_pval < 0.05, "Significant", "Non-significant")) %>% 
        dplyr::arrange(p.value)
    }
    else {
      message(paste0("Mann-Whitney U Test is performed on ", var_levels[1], " - ", var_levels[2], "."))
      p.val <- df %>% dplyr::filter(!(OlinkID %in% npxCheck$all_nas)) %>% 
        dplyr::filter(!(OlinkID %in% nas_in_level)) %>% dplyr::group_by(Assay, OlinkID, UniProt, Panel) %>% 
        dplyr::do(broom::tidy(stats::wilcox.test(!!rlang::ensym(data_type) ~ !!rlang::ensym(variable), 
                                                 data = ., conf.int = TRUE, alternative = "greater", ...))) %>% 
        dplyr::ungroup() %>% dplyr::mutate(Adjusted_pval = p.adjust(p.value, method = "fdr")) %>% 
        dplyr::mutate(Threshold = ifelse(Adjusted_pval < 0.05, "Significant", "Non-significant")) %>% 
        dplyr::arrange(p.value)
    }
    return(p.val)
  }


# ------- Wilcox test with fold change and add column identifier  --------- #
# REMOVE NA ASSAY AND REMOVE NA IN NPX

wilcox_cat_col2 <- function(dat, cat){
  
  print(paste0("Wilcox test on ", cat))
  
  dat1 <- dat  %>%
    dplyr::filter(!is.na(Assay)) 
  
  dat2 <- 
    olink_wilcox_2(df = dat1,
                  variable = cat)
  
  
  # ---------- #
  # Add fold change 
  # ---------- #
  
  # Get the average expression per group
  dat1_g <- dat1 %>%
    dplyr::filter(!is.na(NPX)) %>% 
    dplyr::group_by(OlinkID, !!as.name(cat)) %>%
    dplyr::summarise(NPX_mean = mean(NPX))
  
  # Convert to wide format
  dat1_w <- dat1_g %>%
    tidyr::pivot_wider(names_from = !!as.name(cat), values_from = NPX_mean) %>%
    dplyr::rename(Mean_N = N,
                  Mean_Y= Y) %>%
    dplyr::mutate(log2FC = Mean_N - Mean_Y,
                  baseFC = 2^log2FC)
  
  
  
  dat_tmp <- dat2 %>%
    dplyr::left_join(dat1_w)
  
  
  # Simplify output
  dat3 <- dat_tmp %>%
    dplyr::select(Assay, OlinkID, UniProt, Panel, p.value, Adjusted_pval, Threshold, Mean_Y, Mean_N, log2FC, baseFC) %>%
    dplyr::mutate(Cat = cat)
  
  write.csv(dat3, 
            file.path(wd$outCurr, 
                      paste0("Wilcox_test_Rest_vs_", cat, ".csv")))
  
  return(dat3)
  
}



# Function to analyze Cox regression with highest NPX per patient
analyze_proteins_highest <- function(cd_prot, npx_in, meta_in) {
  all_results <- list()
  
  for (eg_protein in cd_prot) {
    #print(eg_protein)
    # Extract NPX values for the current protein
    npx_exp <- npx_in %>%
      filter(OlinkID == eg_protein) %>%
      dplyr::rename(sample_id = SampleID) %>%
      select(sample_id, NPX) %>%
      mutate(sample_id = as.character(sample_id))
    
    # Combine with metadata
    tmp_df <- meta_in %>% left_join(npx_exp)
    
    # Select the highest NPX per patient
    filtered_data <- select_highest_expression(tmp_df)
    
    # Binarize NPX values
    filtered_data <- filtered_data %>%
      mutate(
        npx_median = median(NPX, na.rm = TRUE),
        npx_bin = case_when(NPX >= npx_median ~ "high",
                            NPX < npx_median ~ "low"),
        npx_bin = factor(npx_bin, levels = c("low", "high"))
      )
    
    # Fit Cox regression models
    model1 <- coxph(Surv(OS, death) ~ NPX, data = filtered_data)
    model2 <- coxph(Surv(OS, death) ~ npx_bin, data = filtered_data)
    
    # Extract results for model1 (continuous NPX)
    coef_summary1 <- summary(model1)$coefficients
    coef_summary2_1 <- summary(model1)$conf.int
    
    result_df1 <- data.frame(
      OlinkID = eg_protein,
      Coefficient = coef_summary1[, 1],
      Lower_Bound = coef_summary2_1[, 3],
      Upper_Bound = coef_summary2_1[, 4],
      P_Value = coef_summary1[, 5],
      HR = coef_summary2_1[, 1]
    )
    
    # Extract results for model2 (binarized NPX)
    coef_summary2 <- summary(model2)$coefficients
    coef_summary2_2 <- summary(model2)$conf.int
    
    result_df2 <- data.frame(
      OlinkID = eg_protein,
      bin_Coefficient = coef_summary2[, 1],
      bin_Lower_Bound = coef_summary2_2[, 3],
      bin_Upper_Bound = coef_summary2_2[, 4],
      bin_P_Value = coef_summary2[, 5],
      bin_HR = coef_summary2_2[, 1]
    )
    
    # Merge both results
    result_df <- left_join(result_df1, result_df2, by = "OlinkID")
    
    # Append to list
    all_results[[eg_protein]] <- result_df
  }
  
  # Combine all results into a single data frame
  final_results <- bind_rows(all_results)
  
  return(final_results)
}

# Function to analyze Cox regression with highest NPX per patient
analyze_proteins_highest_txt <- function(cd_prot, npx_in, meta_in) {
  all_results <- list()
  
  for (eg_protein in cd_prot) {
    #print(eg_protein)
    # Extract NPX values for the current protein
    npx_exp <- npx_in %>%
      filter(OlinkID == eg_protein) %>%
      dplyr::rename(sample_id = SampleID) %>%
      select(sample_id, NPX) %>%
      mutate(sample_id = as.character(sample_id))
    
    # Combine with metadata
    tmp_df <- meta_in %>% left_join(npx_exp)
    
    # Select the highest NPX per patient
    filtered_data <- select_highest_expression(tmp_df)
    
    # Binarize NPX values
    filtered_data <- filtered_data %>%
      mutate(
        npx_median = median(NPX, na.rm = TRUE),
        npx_bin = case_when(NPX >= npx_median ~ "high",
                            NPX < npx_median ~ "low"),
        npx_bin = factor(npx_bin, levels = c("low", "high"))
      )
    
    # Fit Cox regression models
    model1 <- coxph(Surv(OS, death) ~ NPX + txt_stat, data = filtered_data)
    model2 <- coxph(Surv(OS, death) ~ npx_bin + txt_stat, data = filtered_data)
    
    # Extract results for model1 (continuous NPX)
    coef_summary1 <- summary(model1)$coefficients
    coef_summary2_1 <- summary(model1)$conf.int
    
    result_df1 <- data.frame(
      OlinkID = eg_protein,
      Coefficient = coef_summary1[1, 1],
      Lower_Bound = coef_summary2_1[1, 3],
      Upper_Bound = coef_summary2_1[1, 4],
      P_Value = coef_summary1[1, 5],
      HR = coef_summary2_1[1, 1]
    )
    
    # Extract results for model2 (binarized NPX)
    coef_summary2 <- summary(model2)$coefficients
    coef_summary2_2 <- summary(model2)$conf.int
    
    result_df2 <- data.frame(
      OlinkID = eg_protein,
      bin_Coefficient = coef_summary2[1, 1],
      bin_Lower_Bound = coef_summary2_2[1, 3],
      bin_Upper_Bound = coef_summary2_2[1, 4],
      bin_P_Value = coef_summary2[1, 5],
      bin_HR = coef_summary2_2[1, 1]
    )
    
    # Merge both results
    result_df <- left_join(result_df1, result_df2, by = "OlinkID")
    
    # Append to list
    all_results[[eg_protein]] <- result_df
  }
  
  # Combine all results into a single data frame
  final_results <- bind_rows(all_results)
  
  return(final_results)
}


# Function to analyze Cox regression with both continuous and binned NPX
analyze_proteins_multi <- function(cd_prot, npx_in, meta_in) {
  all_results <- list()
  
  for (eg_protein in cd_prot) {
    # Extract NPX values for the current protein
    npx_exp <- npx_in %>%
      filter(OlinkID %in% eg_protein) %>%
      dplyr::rename(sample_id = SampleID) %>%
      select(sample_id, NPX) %>%
      mutate(sample_id = as.character(sample_id))
    
    # Combine with metadata
    tmp_df <- meta_in %>% left_join(npx_exp)
    
    # Select the highest NPX per patient
    filtered_data <- select_highest_expression(tmp_df)
    
    # Binarize NPX values
    filtered_data <- filtered_data %>%
      mutate(
        npx_median = median(NPX, na.rm = TRUE),
        npx_bin = case_when(NPX >= npx_median ~ "high",
                            NPX < npx_median ~ "low"),
        npx_bin = factor(npx_bin, levels = c("low", "high"))
      )
    
    # Fit Cox regression models with additional covariates
    model1 <- coxph(Surv(OS, death) ~ NPX + psa + ldh + alk_ph, data = filtered_data)
    model2 <- coxph(Surv(OS, death) ~ npx_bin + psa + ldh + alk_ph, data = filtered_data)
    
    # Extract results for model1 (continuous NPX)
    coef_summary1 <- summary(model1)$coefficients["NPX", ]
    coef_summary2_1 <- summary(model1)$conf.int["NPX", ]
    
    result_df1 <- data.frame(
      OlinkID = eg_protein,
      Coefficient = coef_summary1[1],
      Lower_Bound = coef_summary2_1[3],
      Upper_Bound = coef_summary2_1[4],
      P_Value = coef_summary1[5],
      HR = coef_summary2_1[1]
    )
    
    # Extract results for model2 (binarized NPX)
    coef_summary2 <- summary(model2)$coefficients["npx_binhigh", ]
    coef_summary2_2 <- summary(model2)$conf.int["npx_binhigh", ]
    
    result_df2 <- data.frame(
      OlinkID = eg_protein,
      bin_Coefficient = coef_summary2[1],
      bin_Lower_Bound = coef_summary2_2[3],
      bin_Upper_Bound = coef_summary2_2[4],
      bin_P_Value = coef_summary2[5],
      bin_HR = coef_summary2_2[1]
    )
    
    # Merge both results
    result_df <- left_join(result_df1, result_df2, by = "OlinkID")
    
    # Append to list
    all_results[[eg_protein]] <- result_df
  }
  
  # Combine all results into a single data frame
  final_results <- bind_rows(all_results)
  
  return(final_results)
}

# Selecting the higest expression from each patient
select_highest_expression <- function(df) {
  df %>%
    group_by(mrn) %>%
    slice_max(order_by = NPX, n = 1, with_ties = FALSE) %>%
    ungroup()
}

